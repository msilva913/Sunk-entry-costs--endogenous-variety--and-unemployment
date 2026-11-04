# ags_run_solution.jl
# =============================================================================
# First-order solution of the AGS model with costly reactivation (Part 1 of D15; equations in
# Notes/ags_equilibrium.md; calibration in ags_steady_state.jl). Toolkit by Salazar-Perez and
# Seoane (Programs baseline/solution_functions.jl). Monthly.
#
# Timing within t (the baseline's, run_solution_core.jl):
#   1. z_t, δ_t, s_t are realized.
#   2. Positions vacated at the end of t−1 are reactivated (rate Λ_r,t = F_r(Q_t)) or retired;
#      new positions e_t are created.
#   3. Production; u_t workers and v_t = v_pre,t + e_t vacancies match.
#   4. Separation s_t and destruction δ_t take effect.
# States: Ẽ_t = employment entering t; Ṽ^u_t = unmatched vacancies carried over;
#   Ṽ^s_t = positions vacated by separation at the end of t−1 (reactivated at t or retired);
#   the shocks z, δ, s (multiplicative deviations around 1; δ·δ̄, s·s̄ are the rates).
# v_pre,t = Ṽ^u_t + Λ_r,t Ṽ^s_t depends on Q_t, so v_pre is a control, not a state.
#
# Sections:
#   1. The model at AGS_TARGETS: steady state, impulse responses of u, v, θ, Q, Λ_r.
#   2. Check: with nearly free vacancy creation (ξ_inv → 0), Q is constant, so job creation pins
#      θ_t as a function of the current shocks alone: θ does not load on the stocks Ẽ, Ṽ^u, Ṽ^s
#      (the analog of the draft's lem:theta; Notes/ags_equilibrium.md §3, special cases).
#   3. D15 gate 1, stage (a): σ(u) by shock and cor(u, v), quarterly averages, HP-1600, at the
#      S14 shock processes (Programs baseline/run_solution_core.jl SHOCKS_BASE).
# Output: ags_run_solution.txt, ags_irf.png
# =============================================================================

cd(@__DIR__)
const BASELINE_DIR = joinpath(@__DIR__, "..", "Programs baseline")
include(joinpath(BASELINE_DIR, "solution_functions.jl"))
include(joinpath(BASELINE_DIR, "time_series_fun.jl"))
include("ags_steady_state.jl")
using DataFrames, Statistics, Printf, Plots

lines = String[]
say(str) = (push!(lines, str); println(str))

# Same values as SHOCKS_BASE in run_solution_core.jl (S14, filter-consistent; part6b).
const SHK = (ρz = 0.9628, σz = 0.0062, ρδ = 0.5319, σδ = 0.0641, ρs = 0.8048, σs = 0.0437)
const T_IR, T_SM, LAM = 60, 60_000, 1_600.0

# =============================================================================
# Model
# =============================================================================
flag_order, flag_deviation, flag_SSsolver = 1, true, false

@syms δbar sbar r A η_L b ϕ κ x_m ξ_inv zbar a p_r ψ_r ρ_z σ_z ρ_δ σ_δ ρ_s σ_s
parameters = [δbar; sbar; r; A; η_L; b; ϕ; κ; x_m; ξ_inv; zbar; a; p_r; ψ_r;
              ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]
estimate, position, priors = [], [], (;)

"""PAR in the order of `parameters`. p_r = 0 is costless reactivation (Λ_r ≡ 1, Q^rep = Q)."""
build_PAR(cal, shk) = [cal.δ; cal.s; cal.r; cal.A; cal.η_L; cal.b; cal.ϕ; cal.κ; cal.x_m;
                       cal.ξ_inv; cal.zbar; cal.a; cal.p_r; cal.ψ_r;
                       shk.ρz; shk.σz; shk.ρδ; shk.σδ; shk.ρs; shk.σs]

@syms Et Vu Vs z δ s Etp Vup Vsp zp δp sp
@syms u θ q v v_pre e Q K Λ_r Q_rep H w
@syms up θp qp vp v_prep ep Qp Kp Λ_rp Q_repp Hp wp
x  = [Et; Vu; Vs; z; δ; s]
y  = [u; θ; q; v; v_pre; e; Q; K; Λ_r; Q_rep; H; w]
xp = [Etp; Vup; Vsp; zp; δp; sp]
yp = [up; θp; qp; vp; v_prep; ep; Qp; Kp; Λ_rp; Q_repp; Hp; wp]
variables = [x; y; xp; yp]
varnames  = vcat(Symbol.(x), Symbol.(y))
nx, ny = length(x), length(y)
nvar, ne = nx + ny, 3
eta = [0 0 0; 0 0 0; 0 0 0; σ_z 0 0; 0 σ_δ 0; 0 0 σ_s]
@syms epsilon
ex = [epsilon]

function ags_equations()
    f    = fill(Sym("x"), nvar)
    β    = 1 / (1 + r)
    S    = 1 - δ * δbar                     # survival of destruction from t to t+1
    R_sh, R_shp = Q - K - H, Qp - Kp - Hp   # 𝓡 = β(1−δ)E[Q' − Q^rep'], by f[7] and f[11]
    ψ_cr = ψ_r / (ψ_r + 1)
    # Post-reactivation stocks and matching
    f[1]  = u - (1 - Et)                                        # unemployment entering matching
    f[2]  = θ - v / u                                           # tightness
    f[3]  = q - A * θ^(-η_L)                                    # vacancy-filling rate
    f[4]  = v - (v_pre + e)                                     # E1: vacancy pool
    f[5]  = v_pre - (Vu + Λ_r * Vs)                             # inherited vacancies
    # Positions: creation, user cost, reactivation
    f[6]  = Q - x_m * e^ξ_inv                                   # creation cost Q(e)
    f[7]  = K - (Q - β * S * Qp)                                # user cost
    f[8]  = κ + K / q - β * S * (                               # E2: job creation
                (1 - ϕ) * (zbar * zp - Kp - b) - ϕ * θp * (Kp + qp * κ)
                + (1 - sp * sbar) * (κ + Kp / qp) - (1 - ϕ) * sp * sbar * R_shp)
    f[9]  = Λ_r - ((1 - p_r) + p_r * (Q / a)^ψ_r)               # Λ_r = F_r(Q); 1 if p_r = 0
    f[10] = Q_rep - Q * (Λ_r - ψ_cr * (Λ_r - 1 + p_r))          # E max(Q − c, 0); Q if p_r = 0
    f[11] = H - β * S * Q_repp                                  # H = β(1−δ)E Q^rep'
    f[12] = w - (ϕ * (zbar * z - K + θ * (K + q * κ) - s * sbar * R_sh) + (1 - ϕ) * b)  # wage
    # Laws of motion (E3, E4)
    f[13] = Etp - S * ((1 - s * sbar) * (1 - u) + A * θ^(1 - η_L) * u)
    f[14] = Vup - S * (1 - q) * v
    f[15] = Vsp - S * s * sbar * (1 - u)
    # Shocks
    f[16] = log(zp) - ρ_z * log(z)
    f[17] = log(δp) - ρ_δ * log(δ)
    f[18] = log(sp) - ρ_s * log(s)
    return f
end

model = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = ags_equations(), nf = nvar, SS = Sym.(zeros(2nvar)), PAR_SS = parameters[:],
         flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
process_model(model)

"""Log steady state in the order of `varnames` (current and next period)."""
function SS_ags(cal)
    st = steady_state_ags(cal)
    Vs = (1 - cal.δ) * cal.s * (1 - st.u)
    Vu = (1 - cal.δ) * (1 - st.q) * st.v
    blk = log.([1 - st.u, Vu, Vs, 1.0, 1.0, 1.0,
                st.u, st.θ, st.q, st.v, st.v_pre, st.e, st.Q, st.K, st.Λ_r, st.Qrep, st.H, st.w])
    return vcat(blk, blk), st
end

"""Solve at targets `tg`; errors on a steady-state residual or a failed Blanchard–Kahn count."""
function solve_ags(tg; shk = SHK)
    cal = calibrate_ags(tg)
    PAR = build_PAR(cal, shk)
    SS, st = SS_ags(cal)
    PAR_SS = Base.invokelatest(eval_PAR_SS, PAR)
    err = maximum(abs, Base.invokelatest(eval_SS_error, PAR_SS, SS))
    if !(err < 1e-10)
        error("AGS steady-state residual $err")
    end
    deriv = Base.invokelatest(eval_deriv, PAR_SS, SS)
    sol   = solve_model(model, deriv, Base.invokelatest(eval_ShockVAR, PAR))
    if !(sol.qzflag == 1)
        error("Blanchard–Kahn failed (qzflag = $(sol.qzflag))")
    end
    return (cal = cal, st = st, SS = SS, sol_mat = sol)
end

shock_col(i, sd) = (c = zeros(nx, 1); c[3 + i, 1] = sd; c)     # i = 1 z, 2 δ, 3 s
const SDS = (SHK.σz, SHK.σδ, SHK.σs)

irf(sol, i) = 100 .* DataFrame(simulate_model((; model..., ne = 1), sol.sol_mat, T_IR,
                                              shock_col(i, SDS[i]), sol.SS, true, true), varnames)

hpq(sim, v) = hp_filter(Float64.(monthly_to_quarterly(sim[!, [v]])[!, v]), LAM)
function moments(sol, cols)
    eta_m = hcat([shock_col(i, SDS[i]) for i in cols]...)
    sim = DataFrame(simulate_model((; model..., ne = length(cols)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    hu, hv = hpq(sim, :u), hpq(sim, :v)
    return (c = cor(hu, hv), su = 100std(hu), sv = 100std(hv))
end

# =============================================================================
# 1. The model at AGS_TARGETS
# =============================================================================
sol = solve_ags(AGS_TARGETS)
say("AGS model with costly reactivation (ags_run_solution.jl). Monthly; responses in % (100 × log dev.).")
say(@sprintf("Targets: Λ̄_r %.2f, p_r %.2f, ψ_r %.2f, ξ_inv %.2f, x_v %.2f, b/w %.2f, X/Y %.3f.",
             AGS_TARGETS.Λr_target, AGS_TARGETS.p_r, AGS_TARGETS.ψ_r, AGS_TARGETS.ξ_inv,
             AGS_TARGETS.x_v, AGS_TARGETS.b_ratio, AGS_TARGETS.X_Y))
say(@sprintf("Steady state: u %.4f, θ %.4f, Q %.4f, Λ_r %.3f; ϕ %.4f, κ %.4f, z̄ %.4f, a %.3f.",
             sol.st.u, sol.st.θ, sol.st.Q, sol.st.Λ_r, sol.cal.ϕ, sol.cal.κ, sol.cal.zbar, sol.cal.a))
say(@sprintf("Shocks (monthly, S14): ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f.\n", SHK...))

say("Impulse responses to one-SD shocks (month 1 = impact; peak = largest |response|):")
say(@sprintf("  %-6s %-6s %8s %8s %8s %8s", "shock", "var", "impact", "peak", "at", "month 24"))
irfs = [irf(sol, i) for i in 1:3]
for (i, name) in enumerate(("z", "δ", "s")), v in (:u, :v, :θ, :Q, :Λ_r)
    p = irfs[i][2:end, v]
    peak = argmax(abs.(p))
    say(@sprintf("  %-6s %-6s %8.3f %8.3f %8d %8.3f", name, v, p[1], p[peak], peak, p[24]))
end

plt = plot(layout = (3, 3), size = (900, 700), legend = false)
for (i, name) in enumerate(("z", "δ", "s")), (j, v) in enumerate((:u, :v, :θ))
    plot!(plt[i, j], 0:T_IR-2, irfs[i][2:end, v], title = "$name → $v", lw = 2)
    hline!(plt[i, j], [0], color = :gray, ls = :dash)
end
savefig(plt, "ags_irf.png")

# =============================================================================
# 2. Check: nearly free vacancy creation ⇒ θ depends on the current shocks only
# =============================================================================
iθ = findfirst(==(:θ), Symbol.(y))
say("\nCheck: θ's loading on the stocks Ẽ, Ṽ^u, Ṽ^s (max |gx|) relative to its loading on z, as ξ_inv → 0.")
for (label, over) in (("costless", (p_r = 0.0, Λr_target = 1.0)), ("costly 0.25", (;)))
    rel = map((1e-2, 1e-4, 1e-6)) do ξ
        g = solve_ags(merge(AGS_TARGETS, over, (ξ_inv = ξ,))).sol_mat.gx[iθ, :]
        maximum(abs, g[1:3]) / abs(g[4])
    end
    say(@sprintf("  %-12s ξ_inv 1e-2: %.1e   1e-4: %.1e   1e-6: %.1e", label, rel...))
    if !(rel[3] < 1e-2 && rel[2] < rel[1] / 10 && rel[3] < rel[2] / 10)
        error("θ's loading on the stocks does not vanish as ξ_inv → 0 ($label)")
    end
end
say("✓ θ is a function of the current shocks alone in the vacancy free entry limit.")

# =============================================================================
# 3. D15 gate 1, stage (a): σ(u) by shock, true shocks
# =============================================================================
say("\nD15 gate 1, stage (a). Simulated $(T_SM) months, quarterly averages, HP λ = $(Int(LAM)).")
say("Data, SVAR A (1992Q3–2019Q4, long run): σ(u) z 7.40, δ 4.27, s 4.48; cor(u,v) 1951–2019 −0.925.")
say(@sprintf("  %-34s %7s %7s | %6s %6s %6s %6s | %5s %5s %5s", "case", "ϕ", "cor", "σu", "z",
             "δ", "s", "%z", "%δ", "%s"))
for (label, over) in (
        ("GS: costless, ξ_inv 1",           (p_r = 0.0, Λr_target = 1.0)),
        ("GS: costless, ξ_inv 0.01",        (p_r = 0.0, Λr_target = 1.0, ξ_inv = 0.01)),
        ("costly 0.25, ξ_inv 1",            (;)),
        ("costly 0.25, ξ_inv 0.25",         (ξ_inv = 0.25,)),
        ("costly 0.25, ξ_inv 0.01",         (ξ_inv = 0.01,)),
        ("costly 0.25, ξ_inv 0.01, x_v 1",  (ξ_inv = 0.01, x_v = 1.0)),
        ("costly 0.70, ξ_inv 0.01",         (Λr_target = 0.70, ξ_inv = 0.01)),
        ("bridge: 0.25, ξ_inv 0.01, X/Y 1.67%", (ξ_inv = 0.01, X_Y = X_Y_BRIDGE)))
    r = try
        solve_ags(merge(AGS_TARGETS, over))
    catch err
        say(@sprintf("  %-34s failed: %s", label, first(sprint(showerror, err), 80)))
        continue
    end
    all3 = moments(r, 1:3)
    one  = [moments(r, [i]).su for i in 1:3]
    sh   = 100 .* one .^ 2 ./ sum(one .^ 2)
    say(@sprintf("  %-34s %7.4f %+7.3f | %6.2f %6.2f %6.2f %6.2f | %5.1f %5.1f %5.1f", label,
                 r.cal.ϕ, all3.c, all3.su, one..., sh...))
end
say("σu: all three shocks; z, δ, s: σ(u) with that shock alone; %: shares of the sum of variances.")

open("ags_run_solution.txt", "w") do io
    foreach(l -> println(io, l), lines)
end
