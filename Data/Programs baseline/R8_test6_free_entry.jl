# R8_test6_free_entry.jl
# =============================================================================
# Q1 Test 6 (Notes/Q1_R8_D1_PATHB_plan.md): the code reproduces the prop:fe closed forms
# c_δ, c_s (gs_free_entry_beveridge.jl a_note, gs_s, evaluated at the code's own steady
# state) in the AGS/GS economy where prop:fe holds:
#   no variety         ρ ≡ 1 (f[13] replaced, as in Comparison C)
#   clone replacement  f[4] (the product-line entry Euler) replaced by Ñ_{t+1} = 1, so entry
#                      replaces exactly what δ destroys
#   σ = 0, κ = 0 (x_v = 1), p_0 = 0, and vacancy free entry ξ_inv → 0.
# Clone replacement matters. With the entry Euler in place, β(λ′/λ)(1−δ_t)(ν′+d′) = ν pins the
# survival-adjusted discount factor m(1−δ_t); the vacancy's user cost K_t = Q_t − m(1−δ_t)Q_{t+1}
# uses the same factor, so with Q pinned the δ shock barely moves K and θ (the real rate falls
# with δ). The first block reports that case for comparison (σ = 1; at σ = 0 the entry Euler
# makes entry infinitely elastic and the solution is not unique).
# Output: R8_test6_free_entry.txt
# =============================================================================

include("run_solution_core.jl")
using Printf

lines = String[]
log!(s) = (push!(lines, s); println(s))
S = SHOCKS_BASE

function c_delta(c, ρ)
    η, δ, ϕ, f, q, τ, r = c.η, c.δ_e, c.ϕ, c.f, c.q, c.τ, c.r
    ω_ρ = 1 + r - ρ * (1 - τ)
    num = (1 + r) * δ / (1 - δ) + δ / (r + δ) * (ω_ρ + ρ * ϕ * (1 - δ) * f + ρ * (1 - ϕ) * (1 - δ) * q)
    return -num / (η * ω_ρ + ρ * ϕ * (1 - δ) * f)
end
c_s(c, ρ) = -ρ * (1 - c.δ_e) * c.s / (c.η * (1 + c.r - ρ * (1 - c.τ)) + ρ * c.ϕ * (1 - c.δ_e) * c.f)

function case(ξi, σv; ags)
    tg = merge(TARGETS, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0, x_v = 1.0, ξ_inv = ξi, σ = σv))
    f = gen_model_equations()
    f[13] = ρ - 1
    ags && (f[4] = Ntp - 1)
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = f, nf = nvar, SS = SS_symbolics(parameters, tg), PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
    process_model(model)
    cal = calibrate_shares(tg)
    @assert abs(cal.κ) < 1e-12 "Test 6 needs κ = 0"
    PAR = build_PAR(cal, S)
    sol = solution_interface(model, PAR, SS_numeric(PAR, tg))
    ss = sol.ss
    c = (η = cal.η_L, δ_e = ss.δ_e, ϕ = cal.ϕ, f = ss.θ * ss.q, q = ss.q,
         τ = compute_separation_rate(ss.δ_e, cal.s), r = cal.r, s = cal.s)
    ir(st) = DataFrame(simulate_model((; model..., ne = 1), sol.sol_mat, 4, shock_col(st, 1.0),
                                      sol.SS, true, true), varnames)
    θδ, θs = ir(:δ).θ[2], ir(:s).θ[2]
    cδ, cs = c_delta(c, S.ρδ), c_s(c, S.ρs)
    return (; θδ, cδ, θs, cs, gap = max(abs(θδ / cδ - 1), abs(θs / cs - 1)))
end
report(lab, r) = log!(@sprintf("  %-14s θ̂/δ̃ = %+.6f vs c_δ = %+.6f;  θ̂/s̃ = %+.6f vs c_s = %+.6f;  gap %.2e",
                               lab, r.θδ, r.cδ, r.θs, r.cs, r.gap))

log!("Test 6: impact θ̂ per unit log shock against the prop:fe closed forms (ρ ≡ 1, κ = 0, p_0 = 0)")
log!("\nWith the product-line entry Euler (σ = 1): the user-cost channel is offset by the discount factor")
for ξi in (1e-2, 1e-4)
    report(@sprintf("ξ_inv = %.0e", ξi), case(ξi, 1.0; ags = false))
end
log!("\nAGS emulation (clone replacement, σ = 0): must converge to c_δ, c_s as ξ_inv → 0")
gaps = Float64[]
for ξi in (1e-5, 1e-6, 1e-7, 1e-8)
    r = case(ξi, 0.0; ags = true)
    push!(gaps, r.gap)
    report(@sprintf("ξ_inv = %.0e", ξi), r)
end
ok6 = issorted(gaps; rev = true) && gaps[end] < 1e-3
log!(@sprintf("\nTest 6 %s (gap falls monotonically to %.2e at ξ_inv = 1e-8; threshold 1e-3)",
              ok6 ? "PASSED" : "FAILED", gaps[end]))

write("R8_test6_free_entry.txt", join(lines, "\n") * "\n")
