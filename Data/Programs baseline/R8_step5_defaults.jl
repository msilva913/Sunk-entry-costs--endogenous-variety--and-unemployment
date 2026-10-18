# R8_step5_defaults.jl
# =============================================================================
# Q1 step 5 (Notes/Q1_R8_D1_PATHB_plan.md): D1 + PATH B are the defaults in TARGETS.
# Checks that calibrate_shares(TARGETS) converges with no overrides, that δ_e/τ = 0.087,
# and that the 41-equation model's steady-state residual is below 1e-12, for costless
# reposting (the default) and the Λ̄_r grid. Reports the calibrated parameter vector.
#
# Output: R8_step5_defaults.txt
# =============================================================================

include("run_solution_core.jl")
using Printf

const SHK = (ρz = 0.902, σz = 0.0092, ρδ = 0.592, σδ = 0.0669, ρs = 0.874, σs = 0.0854)

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("Q1 step 5 — code defaults: D1 (dest_ann = $(TARGETS.dest_ann)) + PATH B " *
     "(dest_elast_target = $(TARGETS.dest_elast_target)), b/w = $(TARGETS.b_ratio), x_v = $(TARGETS.x_v)")

function ss_residual(tg)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHK)
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = gen_model_equations(), nf = nvar,
             SS = SS_symbolics(parameters, tg), PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    sol = solution_interface(model, PAR, SS_numeric(PAR, tg); ss_tol = 1e-12)
    nunstable = count(abs.(eigvals(sol.sol_mat.hx)) .>= 1 - 1e-10)
    return cal, steady_state(cal), nunstable
end

cal, st, nun = ss_residual(TARGETS)
τ = compute_separation_rate(st.δ_e, cal.s)
log!(@sprintf("\nCostless default: SS residual < 1e-12 (enforced), unstable roots %d (want 0)", nun))
log!(@sprintf("  δ_e = %.6f/month, τ = %.4f, δ_e/τ = %.4f (target 0.087)", st.δ_e, τ, st.δ_e / τ))
@assert abs(st.δ_e / τ - 0.087) < 5e-4
log!(@sprintf("  u = %.5f   v = %.5f   θ = %.4f   Q = %.5f   x_c = %.5f", st.u, st.v, st.θ, st.Q, st.x_c))
log!(@sprintf("  X^c/Y = %.4f (outcome)   X/Y = %.4f (target %.3f)   exit elasticity = %.4f (target %.1f)",
              st.X_c / st.Y, st.X / st.Y, TARGETS.X_Y, st.dest_el, TARGETS.dest_elast_target))
log!("\n  Calibrated parameters:")
for k in (:f_e, :z, :δ, :s, :b, :ϕ, :r, :σ, :ε, :A, :η_L, :κ, :ξ_inv, :x_m, :ψ, :f_m, :p_0)
    log!(@sprintf("    %-6s = %.6g", k, getfield(cal, k)))
end

for Λr in (0.95, 0.85, 0.70)
    c, s_, n = ss_residual(merge(TARGETS, (Λr_target = Λr,)))
    log!(@sprintf("\nΛ̄_r = %.2f: SS residual < 1e-12 (enforced), unstable roots %d, α = %.4g, X^r/Y = %.4f%%",
                  Λr, n, c.α, 100 * s_.X_r / s_.Y))
end

# The p_0 = 0 benchmark (AGS-style runs) must still calibrate: PATH B switches off at p_0 = 0.
c0, s0, n0 = ss_residual(merge(TARGETS, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0)))
log!(@sprintf("\np_0 = 0 benchmark: SS residual < 1e-12 (enforced), unstable roots %d, δ_e/τ = %.4f",
              n0, s0.δ_e / compute_separation_rate(s0.δ_e, c0.s)))

write("R8_step5_defaults.txt", join(lines, "\n") * "\n")
