# run_solution_all_delta.jl
# =============================================================================
# Mechanism comparison: role of the LEVEL of δ (destruction rate)
#
# Compares two calibrations of the baseline model:
#   (1) Baseline  — δ_e = δ̄ ≈ 1.08%/month, dest_end_frac = 0.0, p_0 = 0.0, Xc_Y = 0.0
#   (2) High-δ    — δ_e = τ ≈ 3.1%/month,   dest_end_frac = 0.0, p_0 = 0.0, Xc_Y = 0.0
#
# Both variants shut off endogenous exit (p_0=0, dest_end_frac=0) so the only
# moving part between them is the level of δ_e. Xc_Y=0 is required for internal
# consistency: p_0=0 → cons=0 → X_c=0 → Xc_Y=0 in calibrate_shares.
#
# The high-δ calibration sets δ_e = τ so that s = 0 (all separations are firm
# destruction events). This replicates the CK timing assumption — every
# separation destroys both the match and the vacancy — while holding fixed
# matching costs (κ > 0) and household preferences (σ > 0). It isolates the
# role of the destruction rate level in the vacancy and unemployment LOMs from
# other structural differences between models.
#
# Both calibrations target the same u, v, θ, f, q, w, N steady states.
# The high-δ variant recalibrates all derived parameters (z, f_e, κ, ϕ, etc.)
# to be internally consistent at the new δ_e.
#
# IRFs computed for z (technology) and δ (destruction) shocks.
# Serialized output for plotting in separate script.
# =============================================================================

include("run_solution_core.jl")

# ── Shared shock process parameters (monthly, from part6b) ────────────────────
SHK = SHOCKS_BASE   # part6b (run_solution_core.jl); never assign numbers to the symbolic shock names

T_IR        = 60    # 5 years at monthly frequency
flag_IR     = true
flag_logdev = true

# =============================================================================
# Helper: build model, solve, compute IRFs for z and δ shocks
# Returns (ss, irf_z, irf_δ) for a given targets NamedTuple
# =============================================================================
function solve_and_irf(targets_variant)
    cal = calibrate_shares(targets_variant)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

    # Check s ≥ 0 (fails if τ < δ_e, which shouldn't happen but worth catching)
    if s < -1e-8
        error("Negative s = $s: τ < δ_e implies infeasible calibration. " *
              "Check dest_ann and sep targets.")
    end
    s = max(s, 0.0)   # clamp tiny numerical negatives to zero

    zbar = z; δbar = δ; sbar = s

    PAR = build_PAR(cal, SHK)

    # Compute symbolic SS (needed by process_model to generate derivatives)
    SS_symbolic = SS_symbolics(parameters, targets_variant)
    # Compute numeric SS separately (to be passed to solution_interface)
    SS_numeric_val = SS_numeric(PAR, targets_variant)
    f   = gen_model_equations()

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

    # Pass precomputed numeric SS to skip problematic symbolic substitution
    sol = solution_interface(model, PAR, SS_numeric_val)
    @unpack ss, SS, sol_mat = sol   # do NOT unpack eta: sol.eta is Matrix{Sym} (symbolic zeros), not Float64

    # Build numeric single-shock eta columns directly from the scalar σ parameters.
    # State order: [u, N, v_pret, z, δ, s] — shocks hit rows 4 (z), 5 (δ), 6 (s).
    # Shaped as 6×1 matrices so eta[i,:] returns a 1-element Float64 vector inside
    # simulate_model, matching sim_shocks[:,1] when ne=1.
    eta_z_col = shock_col(:z, SHK.σz)
    eta_δ_col = shock_col(:δ, SHK.σδ)

    # ne=1: one shock fired at a time; sim_shocks has 1 column, matching eta_*_col
    model_1shock = (; model..., ne = 1)

    irf_z = simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev)
    irf_z = 100 .* DataFrame(irf_z, varnames)

    irf_δ = simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev)
    irf_δ = 100 .* DataFrame(irf_δ, varnames)

    # Append exit flow δ_e*N in log deviations
    # log dev of (δ_e * N) ≈ log_dev(δ_e) + log_dev(N)  [first-order]
    irf_z[!, :exit_flow] = irf_z.δ_e .+ irf_z.N
    irf_δ[!, :exit_flow] = irf_δ.δ_e .+ irf_δ.N

    return (ss=ss, irf_z=irf_z, irf_δ=irf_δ, cal=cal)
end

# =============================================================================
# (1) Baseline calibration
# =============================================================================
println("\n" * "="^60)
println("Solving BASELINE model (δ_e = δ̄ ≈ 1.08%/month)")
println("="^60)

# Remove endogenous exit and zero fixed costs (Xc_Y=0 is required when p_0=0:
# cons = (ψ/(1+ψ))*p_0 = 0 → X_c = 0 → Xc_Y = 0 for internal consistency in calibrate_shares).
targets_base = (TARGETS..., dest_end_frac=0.0, p_0=0.0, Xc_Y=0.0, b_ratio=0.9, x_v=0.5)

out_base = solve_and_irf(targets_base)
println("Baseline SS: u=$(round(out_base.ss.u, digits=4)), " *
        "δ_e=$(round(out_base.ss.δ_e, digits=5)), " *
        "s=$(round(out_base.cal.s, digits=5))")

# =============================================================================
# (2) High-δ calibration: δ_e = τ, pure exogenous destruction
#
# Key changes from baseline:
#   dest_ann  → annual rate implying monthly δ_e = sep = 0.031
#   dest_end_frac = 0.0  → all destruction exogenous, s = 0 by Stage 1 formula
#   p_0 = 0.0            → endogenous exit distribution mass = 0 (consistent)
#
# All other targets (f, q, sep, N, w, ε, η_L, r_ann, X_Y, Xc_Y) unchanged,
# so steady-state u, v, θ are matched by construction.
# =============================================================================
println("\n" * "="^60)
println("Solving HIGH-δ model (δ_e = τ ≈ 3.1%/month, pure exogenous)")
println("="^60)

sep_monthly = TARGETS.sep                        # = 0.031
# CK's δ_e = τ in the limit. Since R8 the pre-exit separated-vacancy stock
# Ṽ^s = (1−δ)·s̄·(1−u) is a log-linearized state, so s̄ must be strictly positive:
# δ_e = τ(1 − 1e-6) gives s̄ ≈ 3e-8, and the steady state differs from δ_e = τ by ~1e-6.
# Only z and δ shocks are used here, so the s state never moves.
dest_ann_high = 1 - (1 - sep_monthly * (1 - 1e-6))^12   # ≈ 0.317 annually

targets_high_δ = (TARGETS...,
    dest_ann      = dest_ann_high,   # implies δ_e = τ monthly
    dest_end_frac = 0.0,             # all destruction exogenous → s = 0
    p_0           = 0.0,             # no endogenous exit distribution mass
    Xc_Y          = 0.0,            # required: cons=0 when p_0=0 → X_c=0 → Xc_Y=0
    b_ratio       = 0.9,            # higher replacement ratio (exercise)
    x_v           = 0.5,            # same vacancy share target as baseline
)

out_high_δ = solve_and_irf(targets_high_δ)
println("High-δ SS: u=$(round(out_high_δ.ss.u, digits=4)), " *
        "δ_e=$(round(out_high_δ.ss.δ_e, digits=5)), " *
        "s=$(round(out_high_δ.cal.s, digits=5))")

# =============================================================================
# Serialize all output for plotting
# =============================================================================
output = (
    base   = (ss=out_base.ss,    irf_z=out_base.irf_z,    irf_δ=out_base.irf_δ,    cal=out_base.cal,    targets=targets_base),
    high_δ = (ss=out_high_δ.ss,  irf_z=out_high_δ.irf_z,  irf_δ=out_high_δ.irf_δ,  cal=out_high_δ.cal,  targets=targets_high_δ),
)
serialize("irf_all_delta.jls", output)
println("\nSaved: irf_all_delta.jls")

# =============================================================================
# Quick diagnostic: print key SS differences
# =============================================================================
println("\n── Steady-State Comparison ──────────────────────────────────")
for field in [:u, :v, :θ, :δ_e, :N, :N_e, :ν_f, :x_c]
    b_val = getfield(out_base.ss,    field)
    h_val = getfield(out_high_δ.ss,  field)
    println("  $field:  baseline=$(round(b_val, digits=5))  " *
            "high_δ=$(round(h_val, digits=5))")
end
println("  s:  baseline=$(round(out_base.cal.s, digits=5))  " *
        "high_δ=$(round(out_high_δ.cal.s, digits=5))")
