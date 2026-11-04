# ags_steady_state_checks.jl
# =============================================================================
# Checks of ags_steady_state.jl on a grid of cases (reactivation costless or costly; vacancy
# creation elastic or nearly free; per-hire cost on or off):
#   1. Round trip: the steady state solved from the parameters returns the calibration targets
#      (f, q, u, the wage, X/Y, Λ̄_r).
#   2. Identities: vacancy and unemployment flows, job creation, the Nash wage, 𝓡.
#   3. Reactivation: the closed form of E max(Q − c, 0) against numerical integration, and the
#      factor by which Q must rise before F_r's upper bound binds.
# Output: ags_steady_state_checks.txt
# =============================================================================

cd(@__DIR__)
include("ags_steady_state.jl")

lines = String[]
say(str) = (push!(lines, str); println(str))

"""
E max(Q − c, 0) by the midpoint rule, valid for any Q. Accepted costly draws lie on [0, m],
m = min(Q, a); with c = m y^(1/ψ_r) the integrand is smooth in y:
p_r ∫_0^m (Q − c) dF = p_r (m/a)^ψ_r ∫_0^1 (Q − m y^(1/ψ_r)) dy.
"""
function Qrep_numeric(Q, cal; n = 200_000)
    @unpack p_r, ψ_r, a = cal
    m = min(Q, a)
    y = ((1:n) .- 0.5) ./ n
    cont = p_r * (m / a)^ψ_r * sum(Q .- m .* y .^ (1 / ψ_r)) / n
    return (1 - p_r) * Q + cont
end

const CASES = (
    ("GS: costless, ξ_inv 1",          (p_r = 0.0, Λr_target = 1.0)),
    ("costly 0.25, ξ_inv 1",           (;)),
    ("costly 0.25, ξ_inv 0.01",        (ξ_inv = 0.01,)),
    ("costly 0.25, ξ_inv 0.01, x_v 1", (ξ_inv = 0.01, x_v = 1.0)),
    ("costly 0.25, ψ_r 0.5",           (ψ_r = 0.5,)),
    ("costly 0.70",                    (Λr_target = 0.70,)),
    ("costly 0.95",                    (Λr_target = 0.95,)),
)

say("AGS steady-state checks (ags_steady_state.jl). Monthly; residuals in absolute terms.")
say(@sprintf("%-32s %7s %7s %7s %8s %8s %7s %9s %9s %10s", "case", "ϕ", "κ", "Q", "z̄", "a",
             "Λ_r", "bound ×", "worst", "Q^rep err"))
for (label, over) in CASES
    tg  = merge(AGS_TARGETS, over)
    cal = calibrate_ags(tg)
    st  = steady_state_ags(cal)
    @unpack δ, s, r, A, η_L, b, ϕ, κ, x_m, ξ_inv, zbar = cal
    τ  = δ + (1 - δ) * s
    Λ̄  = tg.Λr_target
    X  = st.e * st.Q / (1 + ξ_inv) + κ * st.q * st.v
    Y  = zbar * (1 - st.u)
    res = (
        f_target     = (1 - δ) * st.f - tg.f,
        q_target     = (1 - δ) * st.q - tg.q,
        u_target     = st.u - tg.sep / (tg.sep + tg.f),
        w_target     = st.w - tg.w,
        X_Y_target   = X / Y - tg.X_Y,
        Λ_r_target   = st.Λ_r - Λ̄,
        u_flow       = (1 - δ) * st.f * st.u - τ * (1 - st.u),
        v_flow       = st.v - ((1 - δ) * ((1 - st.q) * st.v + st.Λ_r * s * (1 - st.u)) + st.e),
        v_pre        = st.v_pre - (1 - δ) * ((1 - st.q) * st.v + st.Λ_r * s * (1 - st.u)),
        creation     = st.Q - x_m * st.e^ξ_inv,
        user_cost    = st.K - (st.Q - (1 - δ) * st.Q / (1 + r)),
        R_identity   = st.R_sh - (st.Q - st.K - st.H),
        jcc          = κ + st.K / st.q - (1 - δ) / (1 + r) * ((1 - ϕ) * (zbar - st.K - b)
                       - ϕ * st.θ * (st.K + st.q * κ) + (1 - s) * (κ + st.K / st.q)
                       - (1 - ϕ) * s * st.R_sh),
        wage         = st.w - (ϕ * (zbar - st.K + st.θ * (st.K + st.q * κ) - s * st.R_sh) + (1 - ϕ) * b),
    )
    worst = maximum(abs, values(res))
    qerr  = abs(Qrep_numeric(st.Q, cal) - st.Qrep) / st.Q
    say(@sprintf("%-32s %7.4f %7.4f %7.4f %8.4f %8.4f %7.4f %9.3g %9.2e %10.2e", label, ϕ, κ,
                 st.Q, zbar, cal.a, st.Λ_r, bound_factor(cal, Λ̄), worst, qerr))
    if !(worst < 1e-10)
        say("    FAILED: $(res)")
    end
    if !(qerr  < 1e-8 )
        say("    FAILED: Q^rep closed form against quadrature, relative error $qerr")
    end
    @assert worst < 1e-10 && qerr < 1e-8 "AGS steady-state checks fail for case: $label"
end

# Both branches of Λr_ags and Qrep_ags, away from the steady state: Q from below a to past it.
say("\nReactivation functions over Q/a (both branches), against quadrature:")
say(@sprintf("%-22s %7s %8s %10s %10s %10s", "case", "Q/a", "Λ_r", "Q^rep/Q", "quadrature", "rel. err"))
for (label, over) in (("p_r 0.8, ψ_r 1",   (;)),
                      ("p_r 0.8, ψ_r 0.5", (ψ_r = 0.5,)),
                      ("p_r 0.5, ψ_r 2",   (p_r = 0.5, ψ_r = 2.0, Λr_target = 0.7)))
    cal = calibrate_ags(merge(AGS_TARGETS, over))
    for ratio in (0.0625, 0.5, 1.0, 2.0, 10.0)
        Q    = ratio * cal.a
        qrep = Qrep_ags(Q, cal)
        qnum = Qrep_numeric(Q, cal)
        err  = abs(qrep - qnum) / Q
        say(@sprintf("%-22s %7.4f %8.4f %10.6f %10.6f %10.2e", label, ratio, Λr_ags(Q, cal),
                     qrep / Q, qnum / Q, err))
        if !(err < 1e-8)
            say("    FAILED: Q^rep against quadrature at Q/a = $ratio")
        end
        @assert err < 1e-8 "Qrep_ags fails against quadrature: $label, Q/a = $ratio"
    end
end

say("")
say("u_target compares with the flow-balance unemployment τ/(τ + f), f gross (JOLTS).")
say("bound ×: Q must rise by this factor before Λ_r reaches 1 (Notes/ags_equilibrium.md §5).")
say("✓ All AGS steady-state checks passed.")

open("ags_steady_state_checks.txt", "w") do io
    foreach(l -> println(io, l), lines)
end
