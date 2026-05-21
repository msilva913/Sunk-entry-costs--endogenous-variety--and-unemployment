# run_solution_endog_exit.jl
# =============================================================================
# Mechanism comparison B: role of ENDOGENOUS EXIT
#
# Compares two calibrations:
#   (1) Full baseline — endogenous exit active
#         dest_end_frac = 0.5, p_0 = 0.5, Xc_Y = 0.10
#         δ_e = δ̄ ≈ 1.08%/month (50% endogenous, 50% exogenous)
#
#   (2) Exogenous exit — endogenous exit shut off
#         dest_end_frac = 0.0, p_0 = 0.0, Xc_Y = 0.0
#         δ_e = δ̄ (same level; only the margin changes)
#
# Both variants share the same δ_e LEVEL at SS so the comparison isolates the
# endogenous exit amplification channel (x_c margin, continuation cost
# distribution) from the pure-level effect studied in Comparison A.
#
# Both calibrations use:
#   b_ratio = 0.9   (same as Comparison A exercise)
#   x_v     = 0.5   (half of recruiting cost is flow, half sunk)
#
# IRFs computed for z (technology) and δ (destruction) shocks.
# Serialized output for plotting in plot_endog_exit_comparison.jl.
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
# (1) Full baseline: endogenous exit active
# =============================================================================
println("\n" * "="^60)
println("Solving FULL BASELINE (endogenous exit, b_ratio=0.9, x_v=0.5)")
println("="^60)

targets_endog = (TARGETS..., b_ratio=0.9, x_v=0.5)

out_endog = solve_and_irf_B(targets_endog)
println("Full baseline SS: u=$(round(out_endog.ss.u,   digits=4)), " *
        "δ_e=$(round(out_endog.ss.δ_e, digits=5)), " *
        "s=$(round(out_endog.cal.s,    digits=5)), " *
        "x_c=$(round(out_endog.ss.x_c, digits=4))")

# =============================================================================
# (2) Exogenous exit: endogenous margin shut off, same δ_e level
#
# dest_end_frac = 0.0 → δ_e = δ (purely exogenous)
# p_0 = 0.0          → no continuation cost distribution mass
# Xc_Y = 0.0         → required: cons = 0 when p_0 = 0 → X_c = 0 → Xc_Y = 0
# δ_e level unchanged (dest_ann unchanged from TARGETS) so the only difference
# from the full baseline is the removal of the endogenous margin.
# =============================================================================
println("\n" * "="^60)
println("Solving EXOGENOUS EXIT baseline (b_ratio=0.9, x_v=0.5)")
println("="^60)

targets_exog = (TARGETS..., dest_end_frac=0.0, p_0=0.0, Xc_Y=0.0,
                b_ratio=0.9, x_v=0.5)

out_exog = solve_and_irf_B(targets_exog)
println("Exog. exit SS: u=$(round(out_exog.ss.u,   digits=4)), " *
        "δ_e=$(round(out_exog.ss.δ_e, digits=5)), " *
        "s=$(round(out_exog.cal.s,    digits=5))")

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

# =============================================================================
# Quick diagnostic: print key SS differences
# =============================================================================
println("\n── Steady-State Comparison ──────────────────────────────────────")
for field in [:u, :v, :θ, :δ_e, :N, :N_e, :x_c, :K, :Q]
    e_val = getfield(out_endog.ss, field)
    x_val = getfield(out_exog.ss,  field)
    println("  $field:  endog=$(round(e_val, digits=5))  " *
            "exog=$(round(x_val, digits=5))")
end
println("  s:  endog=$(round(out_endog.cal.s, digits=5))  " *
        "exog=$(round(out_exog.cal.s,  digits=5))")
println("  κ:  endog=$(round(out_endog.cal.κ, digits=5))  " *
        "exog=$(round(out_exog.cal.κ,  digits=5))")
println("  ϕ:  endog=$(round(out_endog.cal.ϕ, digits=5))  " *
        "exog=$(round(out_exog.cal.ϕ,  digits=5))")
