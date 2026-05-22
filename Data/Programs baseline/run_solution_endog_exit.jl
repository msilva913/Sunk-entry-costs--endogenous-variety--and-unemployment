# run_solution_endog_exit.jl
# =============================================================================
# Mechanism comparison B: role of ENDOGENOUS EXIT
#
# Compares two calibrations that share the same EXOGENOUS destruction parameter
# δbar and the same total separation rate τ:
#
#   (1) Full baseline — endogenous exit active
#         dest_end_frac = 0.5, p_0 = 0.5
#         dest_elast_target = 5.0  ← pins x_c channel elasticity directly
#         δ_e_endog = dest_ann-implied rate (~0.653%/month)
#         δbar_endog = (1−dest_end_frac) × δ_e_endog ≈ 0.326%/month
#         τ = sep = 0.031 (target), s_endog derived
#
#   (2) Exogenous exit — endogenous exit shut off
#         dest_end_frac = 0.0, p_0 = 0.0, Xc_Y = 0.0
#         dest_ann_exog set so δ_e_exog = δbar_endog ← SAME exog. destruction
#         τ = sep = 0.031 (SAME target) → s rises to compensate for lower δ_e
#         x_c margin absent
#
# COMPARISON DESIGN (May 22, 2026):
#   Key identification principle: to isolate the x_c dynamic amplification
#   mechanism, the comparison must hold fixed everything except the presence
#   of the endogenous margin. This requires:
#
#   (a) Same δbar: the exog shock to δ fires through δ·δbar. If δbar differs
#       across models, the same proportional δ shock generates different
#       first-period δ_e impulses, confounding the comparison. Setting
#       dest_ann_exog = annual(δbar_endog) equates this direct transmission.
#
#   (b) Same τ (total separation rate): achieved by keeping sep=0.031 as a
#       target in both models. With lower δ_e_exog, s_exog rises automatically
#       inside calibrate_shares so that τ = 1−(1−δ_e)(1−s) = 0.031 in both.
#       This keeps steady-state u, v, θ approximately equal and makes pp-
#       deviation IRF comparisons valid.
#
#   The residual structural difference between models is purely the x_c margin:
#   endog has a continuation cost threshold that responds to profitability,
#   exog does not. Any IRF difference is attributable to that channel alone.
#   Note: s_exog > s_endog (separation rate compensates for absent endogenous
#   destruction) — a structural difference in the separation margin, but one
#   that is small in magnitude and directly implied by the τ-fixing condition.
#
# PREVIOUS VERSION (appendix exercise):
#   The prior design (same dest_ann, same δ_e, different δbar) is available
#   as the "exposure-effect" comparison: it shows that routing half of SS
#   destruction through the endogenous margin halves the model's exposure to
#   the δ shock. Useful as an appendix robustness check but not the clean
#   mechanism comparison for the main text.
#
# Both calibrations use:
#   b_ratio = 0.9   (same as Comparison A exercise)
#   x_v     = 0.5   (half of recruiting cost is flow, half sunk)
#
# IRFs computed for z (technology) and δ (destruction) shocks.
# Serialized output for plotting in plot_endog_exit_comparison.jl.
# =============================================================================

include("run_solution_core.jl")
using JLD2 
# ── Shared shock process parameters (monthly, from part6b) ────────────────────
ρ_z = 0.902
σ_z = 0.0092
ρ_δ = 0.592
σ_δ = 0.0669
ρ_s = 0.90
σ_s = 0.010

T_IR        = 60    # 5 years at monthly frequency
flag_IR     = true
flag_logdev = true

# =============================================================================
# Helper: build model, solve, compute IRFs for z and δ shocks
# (identical structure to run_solution_all_delta.jl)
# =============================================================================
function solve_and_irf_B(targets_variant)
    cal = calibrate_shares(targets_variant)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

    if s < -1e-8
        error("Negative s = $s: τ < δ_e implies infeasible calibration. " *
              "Check dest_ann and sep targets.")
    end
    s = max(s, 0.0)

    zbar = z; δbar = δ; sbar = s

    PAR = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]

    SS_symbolic    = SS_symbolics(parameters, targets_variant)
    SS_numeric_val = SS_numeric(PAR, targets_variant)
    f              = gen_model_equations()

    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar      = length(parameters), ns = length(estimate),
             priors    = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables,
             varnames  = varnames,
             nx = nx, ny = ny, nvar = nvar,
             e = ex, eta = eta,
             ne = ne,
             f = f,
             nf = nvar,
             SS = SS_symbolic, PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)

    sol = solution_interface(model, PAR, SS_numeric_val)
    @unpack ss, SS, sol_mat = sol

    eta_z_col = reshape([0.0, 0.0, 0.0, σ_z, 0.0, 0.0], nx, 1)
    eta_δ_col = reshape([0.0, 0.0, 0.0, 0.0, σ_δ, 0.0], nx, 1)

    model_1shock = (; model..., ne = 1)

    irf_z = simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev)
    irf_z = 100 .* DataFrame(irf_z, varnames)

    irf_δ = simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev)
    irf_δ = 100 .* DataFrame(irf_δ, varnames)

    # Append exit flow δ_e × N in log deviations
    # log dev of (δ_e × N) ≈ log_dev(δ_e) + log_dev(N)  [first-order]
    irf_z[!, :exit_flow] = irf_z.δ_e .+ irf_z.N
    irf_δ[!, :exit_flow] = irf_δ.δ_e .+ irf_δ.N

    return (ss=ss, irf_z=irf_z, irf_δ=irf_δ, cal=cal)
end

# =============================================================================
# Pre-compute δbar_endog so the exog model can match it exactly.
# δbar = (1 − dest_end_frac) × δ_e is the exogenous component of destruction
# in the endog model. Setting dest_ann_exog = annual(δbar_endog) ensures both
# models receive the same direct impulse from a proportional δ shock.
# =============================================================================
const dest_end_frac_endog = 0.5
const δ_e_endog     = 1 - (1 - TARGETS.dest_ann)^(1/12)
const δbar_endog    = (1 - dest_end_frac_endog) * δ_e_endog #exogenous component of destruction in endog model
const dest_ann_exog = 1 - (1 - δbar_endog)^12   # annual rate matching δbar_endog

println("── δbar pre-computation ─────────────────────────────────────")
println("  δ_e_endog   = $(round(δ_e_endog,    digits=6))/month " *
        "(dest_ann = $(TARGETS.dest_ann))")
println("  δbar_endog  = $(round(δbar_endog,   digits=6))/month " *
        "(exog. component, dest_end_frac = $dest_end_frac_endog)")
println("  dest_ann_exog = $(round(dest_ann_exog, digits=6)) " *
        "(annual rate for exog model; = annual(δbar_endog))")
println("  Implied s_exog > s_endog: sep=0.031 target unchanged, " *
        "lower δ_e means higher conditional separation to keep τ fixed.")

# =============================================================================
# (1) Full baseline: endogenous exit active
# =============================================================================
println("\n" * "="^60)
println("Solving FULL BASELINE (endogenous exit, b_ratio=0.9, x_v=0.5)")
println("="^60)

# dest_elast_target pins the x_c channel elasticity directly (PATH B in
# calibrate_shares). Xc_Y is an outcome, not a target.
# dest_elast = 5 → ψ ≈ 0.033, Xc_Yc ≈ 14.5% (Abraham et al. 2019 range).
targets_endog = (TARGETS..., b_ratio=0.9, x_v=0.5, dest_elast_target=5.0)

out_endog = solve_and_irf_B(targets_endog)
println("Full baseline SS: u=$(round(out_endog.ss.u,   digits=4)), " *
        "δ_e=$(round(out_endog.ss.δ_e, digits=5)), " *
        "s=$(round(out_endog.cal.s,    digits=5)), " *
        "x_c=$(round(out_endog.ss.x_c, digits=4)), " *
        "dest_el=$(round(out_endog.cal.dest_el, digits=3)), " *
        "ψ=$(round(out_endog.cal.ψ,    digits=4))")

# =============================================================================
# (2) Exogenous exit: endogenous margin shut off
#
# Design: same δbar as endog model, same τ target (sep=0.031 unchanged).
#   dest_ann_exog → δ_e_exog = δbar_endog   [same exog. destruction]
#   dest_end_frac = 0.0 → all destruction is exogenous
#   p_0 = 0.0, Xc_Y = 0.0 → no continuation cost distribution
#   sep = 0.031 (unchanged) + lower δ_e → s_exog > s_endog automatically
#
# Structural difference from endog model: x_c margin absent. All other
# steady-state labor market objects (u, v, θ) are approximately equal.
# =============================================================================
println("\n" * "="^60)
println("Solving EXOGENOUS EXIT (same δbar, same τ, b_ratio=0.9, x_v=0.5)")
println("="^60)

targets_exog = (TARGETS...,
    dest_ann      = dest_ann_exog,  # δ_e_exog = δbar_endog: same exog. destruction
    dest_end_frac = 0.0,            # all destruction is exogenous
    p_0           = 0.0,            # no continuation cost distribution
    Xc_Y          = 0.0,            # X_c = 0 when p_0 = 0
    b_ratio       = 0.9,
    x_v           = 0.5)
    # sep = 0.031 inherited from TARGETS: τ fixed, s rises automatically

out_exog = solve_and_irf_B(targets_exog)
println("Exog. exit SS: u=$(round(out_exog.ss.u,   digits=4)), " *
        "δ_e=$(round(out_exog.ss.δ_e, digits=5)), " *
        "s=$(round(out_exog.cal.s,    digits=5)), " *
        "Δs = $(round(out_exog.cal.s - out_endog.cal.s, digits=5)) " *
        "(s rises to compensate for absent endogenous exit)")

# =============================================================================
# Serialize all output for plotting
# =============================================================================
output_B = (
    endog = (ss=out_endog.ss, irf_z=out_endog.irf_z, irf_δ=out_endog.irf_δ,
             cal=out_endog.cal, targets=targets_endog),
    exog  = (ss=out_exog.ss,  irf_z=out_exog.irf_z,  irf_δ=out_exog.irf_δ,
             cal=out_exog.cal,  targets=targets_exog),
)
serialize("irf_endog_exit.jls", output_B)
println("\nSaved: irf_endog_exit.jls")
@save "model_results.jld2" output_B # for external readability
# =============================================================================
# Quick diagnostic: print key SS differences
# =============================================================================
println("\n── Steady-State Comparison ──────────────────────────────────────")
println("  [Design check: δbar should be equal; τ should be equal; δ_e should differ]")
τ_endog = compute_separation_rate(out_endog.ss.δ_e, out_endog.cal.s)
τ_exog  = compute_separation_rate(out_exog.ss.δ_e,  out_exog.cal.s)

println("  τ:       endog=$(round(τ_endog, digits=5))  " *
        "exog=$(round(τ_exog, digits=5))  ← should match sep=$(TARGETS.sep)")
println("  δ_e:     endog=$(round(out_endog.ss.δ_e, digits=5))  " *
        "exog=$(round(out_exog.ss.δ_e,  digits=5))  ← exog = δbar_endog by design")
δbar_check = (1 - dest_end_frac_endog) * out_endog.ss.δ_e
println("  δbar:    endog=$(round(δbar_check, digits=5))  " *
        "exog=$(round(out_exog.ss.δ_e, digits=5))  ← should be equal")
println("  s:       endog=$(round(out_endog.cal.s, digits=5))  " *
        "exog=$(round(out_exog.cal.s,  digits=5))  ← exog higher: compensates for lower δ_e")
println()

for field in [:u, :v, :θ, :N, :N_e, :K, :Q]
    e_val = getfield(out_endog.ss, field)
    x_val = getfield(out_exog.ss,  field)
    println("  $field:  endog=$(round(e_val, digits=5))  exog=$(round(x_val, digits=5))")
end

println("  κ:       endog=$(round(out_endog.cal.κ,       digits=5))  " *
        "exog=$(round(out_exog.cal.κ,       digits=5))")
println("  ϕ:       endog=$(round(out_endog.cal.ϕ,       digits=5))  " *
        "exog=$(round(out_exog.cal.ϕ,       digits=5))")
println("  ψ:       endog=$(round(out_endog.cal.ψ,       digits=4))  " *
        "exog=$(round(out_exog.cal.ψ,       digits=4))  (exog: placeholder, p_0=0)")
println("  dest_el: endog=$(round(out_endog.cal.dest_el, digits=3))  " *
        "exog=$(round(out_exog.cal.dest_el, digits=3))  (exog=0 by construction)")
println("\n  Calibration note:")
println("    Endog: dest_elast_target=$(targets_endog.dest_elast_target) [PATH B]")
println("    Exog:  dest_ann=$(round(dest_ann_exog,digits=6)) → δ_e=δbar_endog; " *
        "Xc_Y=0.0 [PATH A, trivial]")
