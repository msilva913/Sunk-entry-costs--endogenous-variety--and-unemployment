# run_solution_entry_elasticity.jl
# =============================================================================
# Mechanism comparison D: role of ENTRY COST ELASTICITY (xi_inv)
#
# Compares the baseline (ξ_inv = 1.0) against a near-free-entry
# counterfactual (ξ_inv = 0.1). Both are fully recalibrated to the
# same empirical targets; all parameter differences are endogenous to
# the recalibration, so IRF differences are attributable to the entry
# margin alone.
#
# ECONOMIC MECHANISM
# ------------------
# ξ_inv governs curvature of the vacancy-posting cost schedule:
#   Q_t = x_m * e_t^ξ_inv  (marginal cost of posting the e_t-th vacancy)
#
# ξ_inv = 1 (baseline): quadratic total cost; entry responds sluggishly
#   because each additional vacancy drives up marginal cost.
# ξ_inv -> 0 (near-free entry): Q approx x_m (constant marginal cost);
#   entry can respond freely to profit opportunities without congestion.
#
# For TECHNOLOGY SHOCKS, the key amplification chain is:
#   z up -> ν_f up (firm value) -> e up (cheap entry)
#        -> N up -> ρ = N^{1/(ε-1)} up -> w_int up
#        -> JCC tightens -> θ up, u falls, v rises
# With ξ_inv close to 0, e responds elastically to ν_f, amplifying the
# variety feedback. The baseline shows low tech-shock amplification;
# this comparison isolates how much the convex entry cost suppresses it.
#
# CALIBRATION IMPLICATIONS
# ------------------------
# Holding all other targets fixed (X_Y, dest_elast_target, b_ratio,
# x_v, dest_ann, sep), lower ξ_inv implies:
#   - Lower Q and K  (X_v = e/(1+ξ_inv)*Q must match X_Y; ξ_inv -> 0
#     pushes X_v toward e*Q, so Q falls to keep X_v/Y = X_Y).
#   - Lower κ        (x_v = (K/q)/(κ+K/q) fixed; K falls so κ falls).
#   - Higher ϕ       (Nash residual: lower K shrinks denominator).
#   - x_m = Q/e^ξ_inv; at ξ_inv = 0 exactly, x_m = Q (constant cost).
# All other SS quantities (u, v, θ, δ_e) are pinned by targets and
# are unchanged across the two calibrations.
#
# COMPARISON DESIGN (consistent with B and C)
# -------------------------------------------
# Shared targets: dest_end_frac=0.5, p_0=0.5, dest_elast_target=5.0,
#   b_ratio=0.9, x_v=0.5, sep=0.031. Only ξ_inv differs.
# Output serialized in irf_xi_inv.jls for plot_xi_inv_comparison.jl.
# =============================================================================

include("run_solution_core.jl")

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
# Helper: calibrate, build model, solve, compute IRFs for z and δ shocks.
# No model equation modification: ξ_inv enters only through PAR and the
# calibrated (x_m, κ, ϕ). process_model regenerates the same symbolic
# Jacobian each call (equations are ξ_inv-agnostic at the symbolic level).
# =============================================================================
function solve_and_irf_D(targets_variant)
    cal = calibrate_shares(targets_variant)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

    if s < -1e-8
        error("Negative s = $s: τ < δ_e implies infeasible calibration.")
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

    irf_z = 100 .* DataFrame(
        simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev),
        varnames)
    irf_δ = 100 .* DataFrame(
        simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev),
        varnames)

    irf_z[!, :exit_flow] = irf_z.δ_e .+ irf_z.N
    irf_δ[!, :exit_flow] = irf_δ.δ_e .+ irf_δ.N

    return (ss=ss, irf_z=irf_z, irf_δ=irf_δ, cal=cal)
end

# =============================================================================
# (1) Baseline: ξ_inv = 1.0 — quadratic entry cost
# Same targets as Comparisons B and C (endog exit, no-variety).
# =============================================================================
println("\n" * "="^60)
println("Solving BASELINE  (ξ_inv = 1.0, quadratic entry cost)")
println("="^60)

targets_base = (TARGETS...,
    dest_end_frac     = 0.5,
    p_0               = 0.5,
    dest_elast_target = 5.0,
    b_ratio           = 0.9,
    x_v               = 0.5,
    ξ_inv             = 1.0,
)

out_base = solve_and_irf_D(targets_base)
println("Baseline SS: u=$(round(out_base.ss.u, digits=4)), " *
        "δ_e=$(round(out_base.ss.δ_e, digits=5)), " *
        "ϕ=$(round(out_base.cal.ϕ, digits=4))")

# =============================================================================
# (2) Low ξ_inv: near-constant marginal entry cost
# With ξ_inv = 0.1, Q ≈ x_m (entry cost nearly flat in e),
# entry responds elastically to ν_f. ϕ rises (lower K, smaller firm surplus).
# =============================================================================
println("\n" * "="^60)
println("Solving LOW ξ_inv (ξ_inv = 0.1, near-constant entry cost)")
println("="^60)

targets_low_ξ = (TARGETS...,
    dest_end_frac     = 0.5,
    p_0               = 0.5,
    dest_elast_target = 5.0,
    b_ratio           = 0.9,
    x_v               = 0.5,
    ξ_inv             = 0.1,
)

out_low_ξ = solve_and_irf_D(targets_low_ξ)
println("Low ξ SS:    u=$(round(out_low_ξ.ss.u, digits=4)), " *
        "δ_e=$(round(out_low_ξ.ss.δ_e, digits=5)), " *
        "ϕ=$(round(out_low_ξ.cal.ϕ, digits=4))")

# =============================================================================
# Serialize for plotting
# =============================================================================
output_D = (
    base  = (ss=out_base.ss,  irf_z=out_base.irf_z,  irf_δ=out_base.irf_δ,  cal=out_base.cal,  targets=targets_base),
    low_ξ = (ss=out_low_ξ.ss, irf_z=out_low_ξ.irf_z, irf_δ=out_low_ξ.irf_δ, cal=out_low_ξ.cal, targets=targets_low_ξ),
)
serialize("irf_xi_inv.jls", output_D)
println("\nSaved: irf_xi_inv.jls")

# =============================================================================
# Steady-state comparison
# =============================================================================
println("\n── Steady-State Comparison ──────────────────────────────────")
for field in [:u, :v, :θ, :δ_e, :N, :N_e, :ν_f, :x_c, :K, :Q]
    b_val = getfield(out_base.ss,  field)
    l_val = getfield(out_low_ξ.ss, field)
    println("  $field:  baseline=$(round(b_val, digits=5))  low_ξ=$(round(l_val, digits=5))")
end
# Cal-level parameters that vary with ξ_inv (not in ss)
println("  ϕ:     baseline=$(round(out_base.cal.ϕ,     digits=4))  low_ξ=$(round(out_low_ξ.cal.ϕ,     digits=4))")
println("  κ:     baseline=$(round(out_base.cal.κ,     digits=5))  low_ξ=$(round(out_low_ξ.cal.κ,     digits=5))")
println("  x_m:   baseline=$(round(out_base.cal.x_m,   digits=4))  low_ξ=$(round(out_low_ξ.cal.x_m,   digits=4))")
println("  ξ_inv: baseline=$(out_base.cal.ξ_inv)        low_ξ=$(out_low_ξ.cal.ξ_inv)")
