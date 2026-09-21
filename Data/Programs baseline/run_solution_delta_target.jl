# run_solution_delta_target.jl
# =============================================================================
# CALIBRATION DIAGNOSTIC: does the model still fit at the measured δ̄_e?
#
# This is the diagnostic called for by decision D1 (context/decisions.md) and
# task M5 (context/pending_tasks.md). It is NOT one of the paper's mechanism
# comparisons A–D — it exists to decide whether δ_e can be FIXED externally at
# the BED Deaths measurement or must be ESTIMATED as part of Θ_e.
#
# THE QUESTION
# ------------
# Measured permanent establishment exit (BLS BED, dataclass 08, employment-
# weighted, 1993–2019) is 0.812 %/qtr → dest_ann ≈ 3.2 %/yr, i.e. δ_e/τ ≈ 0.087.
# That is roughly HALF the floor of the 6–10 %/yr range that Gabrovski-Silva
# (JEDC) find delivers realistic dynamic correlations — a range the paper's own
# introduction endorses at Draft.tex:421 — and well below BGM (2012), who set
# δ = 0.025/qtr ≈ 9.6 %/yr.
#
# WHAT COMPARISONS A AND B PREDICT
# --------------------------------
# From the N law of motion, the entry coefficient is δ_e and the survival
# coefficient (1−δ_e): a LOWER δ_e means a slower-turning-over firm stock, hence
#   (a) MORE persistence  — N recovers slowly, sustaining ρ↓ → w_int↓ → θ↓,
#       which HELPS match the LP evidence (δ→u peaks at h = 17–20 quarters);
#   (b) LESS amplification — a smaller LOM multiplier on entry, which HURTS the
#       Shimer metric σ(θ)/σ(z) = 11.70.
# The two estimation blocks therefore pull in OPPOSITE directions on δ_e, which
# is good for identification but means the fixed-δ_e choice cannot be made by
# measurement alone.
#
# THE TEST
# --------
# Solve three otherwise-identical calibrations and report, for each:
#   1. steady state,
#   2. business-cycle moments (HP λ=1,600 on quarterly aggregates of a monthly
#      simulation — the Block M object),
#   3. δ-shock IRF persistence (peak horizon, peak, value at h = 20 quarters).
#
#   BED  dest_ann = 0.0320  (measured; the D1 recommendation)
#   CODE dest_ann = 0.0754  (status quo in steady_state.jl — the JF-derived
#                            figure D1 shows to be a denominator error)
#   BGM  dest_ann = 0.0963  (Bilbiie-Ghironi-Melitz 2012, δ = 0.025/qtr)
#
# DECISION RULE
# -------------
# If σ(θ)/σ(labor_prod) at BED is close to the CODE/BGM specs, δ_e can be fixed
# at the measured value and D1 is settled. If amplification collapses while the
# δ→u persistence improves, the two blocks disagree and δ_e should move into Θ_e
# with a prior anchored at the BED value as a lower bound (see D1).
#
# Everything except dest_ann is held at the Comparison B/D shared settings, so
# any difference is attributable to the destruction rate alone.
#
# Output: irf_delta_target.jls  →  plot_delta_target_comparison.jl
# =============================================================================

include("run_solution_core.jl")

# ── Shared shock process parameters (monthly, from part6b) ────────────────────
ρ_z = 0.902
σ_z = 0.0092
ρ_δ = 0.592
σ_δ = 0.0669
ρ_s = 0.8741      # part6b: AR(1) on s = (τ−δ_e)/(1−δ_e), 1992Q3–2019Q4
σ_s = 0.0854      # part6b: AR(1) residual SD. Was 0.010, a placeholder 8.5× too small

T_IR        = 63      # 63 months = 21 quarters = LP horizons h = 0..20 inclusive
T_SM        = 60_000  # months of stochastic simulation for Block M moments
flag_IR     = true
flag_logdev = true

const MOM_LAMBDA  = 1_600.0   # principles.md §18 — must match the data side
const TARGET_AMP  = 11.70     # σ(θ)/σ(z) in Table 4 (tab:smm_moments)

# =============================================================================
# Helper: calibrate, solve, and return IRFs + simulated moments.
# Mirrors solve_and_irf_D in run_solution_entry_elasticity.jl.
# =============================================================================
function build_model_once()
    f = gen_model_equations()
    # SS_symbolics needs a valid targets_variant for the symbolic expressions;
    # any spec works — process_model only uses SS to build eval_SS, which we
    # bypass via SS_precomputed.
    dummy_targets = (shared..., dest_ann = specs[1][2])
    SS_sym = SS_symbolics(parameters, dummy_targets)
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
             SS = SS_sym, PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    return model
end

function solve_and_diagnose(model, targets_variant, label)
    cal = calibrate_shares(targets_variant)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

    if s < -1e-8
        error("$label: negative s = $s — τ < δ_e implies an infeasible calibration. " *
              "dest_ann too high relative to sep.")
    end
    s = max(s, 0.0)

    zbar = z; δbar = δ; sbar = s

    PAR = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]

    SS_numeric_val = SS_numeric(PAR, targets_variant)

    sol = solution_interface(model, PAR, SS_numeric_val)
    # NB: do NOT unpack `eta` here — assigning it anywhere in this function would
    # make the name local and shadow the global `eta` used in `model` above.
    @unpack ss, SS, sol_mat = sol
    eta_full = sol.eta

    # ── IRFs: one shock at a time (state order x = [u, N, v_pret, z, δ, s]) ────
    eta_z_col = reshape([0.0, 0.0, 0.0, σ_z, 0.0, 0.0], nx, 1)
    eta_δ_col = reshape([0.0, 0.0, 0.0, 0.0, σ_δ, 0.0], nx, 1)
    model_1shock = (; model..., ne = 1)

    irf_z = 100 .* DataFrame(
        simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev),
        varnames)
    irf_δ = 100 .* DataFrame(
        simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev),
        varnames)

    # ── Block M moments: all three shocks on, monthly → quarterly → HP(1600) ──
    mom = simulated_moments(model, sol_mat, eta_full, SS, ss, η_L, sbar)

    return (ss=ss, cal=cal, irf_z=irf_z, irf_δ=irf_δ, mom=mom, label=label)
end

# =============================================================================
# Block M moment construction.
#
# This is a prototype of the m(θ) simulator that Section 5.2 requires: it
# replaces second_moments.jl. Filter symmetry with the data side is mandatory — principles.md §18.
#
# NOTE the model counterpart of "labor productivity" in tab:smm_moments is the
# VARIABLE labor_prod (= Y/(ρL), equation f[26]), NOT the technology shock z.
# The empirical series is BLS PRS85006163, output per hour.
# =============================================================================
function simulated_moments(model, sol_mat, eta_full, SS, ss, η_L, sbar)
    sim = DataFrame(
        simulate_model(model, sol_mat, T_SM, eta_full, SS, false, true),  # log deviations
        varnames)

    # Job-finding rate: f = A θ^(1-η_L)  ⇒  log-dev  f̂ = (1-η_L)·θ̂
    sim[!, :f_rate] = (1 - η_L) .* sim.θ

    # Total separation rate τ_t = δ_e,t + s_t·s̄·(1−δ_e,t)   (eq:tau_lom, f[23]).
    # Nonlinear in its components, so rebuild in LEVELS then re-log.
    # CAREFUL: the state `s` is the multiplicative SHOCK (SS = 1); the separation
    # RATE is s_t·s̄, with s̄ = sbar the calibrated parameter. ss.s is 1.0, not s̄.
    δ_e_lvl = ss.δ_e .* exp.(sim.δ_e)
    s_lvl   = sbar   .* exp.(sim.s)
    τ_lvl   = δ_e_lvl .+ s_lvl .* (1 .- δ_e_lvl)

    sim[!, :τ_rate] = log.(τ_lvl ./ mean(τ_lvl))

    mom_vars = [:u, :v, :θ, :τ_rate, :f_rate, :δ_e, :N_e, :labor_prod]
    sim_q    = monthly_to_quarterly(sim[!, mom_vars])

    sim_hp = copy(sim_q)
    for c in mom_vars
        sim_hp[!, c] = hp_filter(Float64.(sim_q[!, c]), MOM_LAMBDA)
    end

    # lags=2 is required: `moments` hard-codes two autocorrelation column labels.
    return moments(sim_hp, :labor_prod, [:u, :v, :labor_prod]; lags=2, verbose=false)
end

# =============================================================================
# IRF persistence diagnostics for the δ shock.
# Converts the monthly IRF to quarterly averages so horizons are directly
# comparable to the Bartik LP (h in quarters, peak at h = 17–20 for u).
# =============================================================================
function irf_diagnostics(irf, ss)
    to_pp(col, ss_val) = (exp.(col ./ 100) .- 1) .* ss_val .* 100

    u_pp_m = to_pp(irf.u, ss.u)
    v_pp_m = to_pp(irf.v, ss.v)
    q(x)   = [mean(x[i:i+2]) for i in 1:3:length(x)-2]   # monthly → quarterly
    u_q, v_q = q(u_pp_m), q(v_pp_m)

    u_peak_h = argmax(abs.(u_q)) - 1          # 0-indexed horizon, quarters
    v_peak_h = argmax(abs.(v_q)) - 1

    return (u_peak = u_q[u_peak_h + 1], u_peak_h = u_peak_h, u_h20 = u_q[end],
            v_peak = v_q[v_peak_h + 1], v_peak_h = v_peak_h, v_h20 = v_q[end],
            u_q = u_q, v_q = v_q)
end

# =============================================================================
# Three calibrations. Only dest_ann differs.
# =============================================================================
shared = (TARGETS...,
    dest_end_frac     = 0.5,
    p_0               = 0.5,
    dest_elast_target = 5.0,
    b_ratio           = 0.9,
    x_v               = 0.5,
    ξ_inv             = 1.0,
)

specs = [
    ("BED  (measured, D1)",  0.0320),
    ("CODE (status quo)",    0.0754),
    ("BGM  (2012)",          0.0963),
]

println("\nBuilding symbolic model (once) ...")
model = build_model_once()
println("Symbolic model ready.\n")

results = Dict{String,Any}()
for (label, dest_ann) in specs
    println("\n" * "="^68)
    println("Solving  $label   dest_ann = $dest_ann")
    println("="^68)
    results[label] = solve_and_diagnose(model, (shared..., dest_ann = dest_ann), label)
end

# =============================================================================
# Report
# =============================================================================
using Printf

sig3(x) = round(x, sigdigits=3)
fmt3(x) = x isa Number ? rpad(string(sig3(x)), 15) : rpad(string(x), 15)

function display_moments(mom)
    dm = copy(mom)
    dm.SD  = sig3.(Float64.(dm.SD)  .* 100)
    dm.RSD = sig3.(Float64.(dm.RSD))
    for c in names(dm)
        c == "Variable" && continue
        dm[!, c] = sig3.(Float64.(dm[!, c]))
    end
    show(stdout, MIME("text/plain"), dm)
    println()
end

τ_monthly = shared.sep

println("\n\n" * "="^68)
println("STEADY STATE")
println("="^68)
println(rpad("", 22), join([rpad(l, 15) for (l, _) in specs]))
for field in [:u, :v, :θ, :δ_e, :N, :N_e, :ν_f, :K, :Q]
    vals = [fmt3(getfield(results[l].ss, field)) for (l, _) in specs]
    println(rpad("  $field", 22), join(vals))
end
println(rpad("  δ_e/τ", 22),
        join([fmt3(results[l].ss.δ_e / τ_monthly) for (l, _) in specs]))
println(rpad("  ϕ (Nash resid.)", 22),
        join([fmt3(results[l].cal.ϕ) for (l, _) in specs]))

println("\n" * "="^68)
println("BLOCK M MOMENTS   (HP λ=$(Int(MOM_LAMBDA)), quarterly, T=$(T_SM) months)")
println("  SD in percent (log-dev × 100); RSD and correlations dimensionless")
println("="^68)
for (l, _) in specs
    println("\n── $l ──")
    display_moments(results[l].mom)
end

println("\n" * "="^68)
println("HEADLINE DIAGNOSTICS")
println("="^68)
println(rpad("", 34), join([rpad(l, 15) for (l, _) in specs]))

amp = Dict(l => results[l].mom[results[l].mom.Variable .== "θ", "RSD"][1] for (l, _) in specs)
println(rpad("  σ(θ)/σ(labor_prod)", 34),
        join([fmt3(amp[l]) for (l, _) in specs]),
        "   target $TARGET_AMP")

cor_uv = Dict(l => results[l].mom[results[l].mom.Variable .== "u", "Cor(x, v)"][1] for (l, _) in specs)
println(rpad("  cor(u, v)", 34),
        join([fmt3(cor_uv[l]) for (l, _) in specs]),
        "   target −0.804")

sd_u = Dict(l => Float64(results[l].mom[results[l].mom.Variable .== "u", "SD"][1]) * 100 for (l, _) in specs)
println(rpad("  σ(u) (%)", 34),
        join([fmt3(sd_u[l]) for (l, _) in specs]),
        "   target 14.96")

dg = Dict(l => irf_diagnostics(results[l].irf_δ, results[l].ss) for (l, _) in specs)
println(rpad("  δ→u peak (pp)", 34),
        join([fmt3(dg[l].u_peak) for (l, _) in specs]),
        "   LP: +1.71")
println(rpad("  δ→u peak horizon (qtrs)", 34),
        join([rpad(string(dg[l].u_peak_h), 15) for (l, _) in specs]),
        "   LP: 17–20")
println(rpad("  δ→u at h=20 (pp)", 34),
        join([fmt3(dg[l].u_h20) for (l, _) in specs]))
println(rpad("  δ→v trough (pp)", 34),
        join([fmt3(dg[l].v_peak) for (l, _) in specs]),
        "   LP: −0.58")
println(rpad("  δ→v trough horizon (qtrs)", 34),
        join([rpad(string(dg[l].v_peak_h), 15) for (l, _) in specs]),
        "   LP: 18")

println("""

INTERPRETATION GUIDE
  • If σ(θ)/σ(labor_prod) is broadly similar across specs → δ_e can be FIXED at
    the measured BED value; D1 resolves to 0.0320 and no further work is needed.
  • If amplification falls sharply in the BED arm while δ→u persistence improves
    → the two estimation blocks disagree on δ_e. Move δ_e into Θ_e with a prior
    anchored at 0.0320 as a lower bound (D1's preferred resolution), and let the
    posterior adjudicate the product-line-vs-establishment measurement gap.
  • Note the IRF scale is NOT comparable to the LP levels until D3 is settled —
    compare SHAPE and PEAK HORIZON, not magnitude, against the LP column.
""")

pack(l, d) = (ss=results[l].ss, cal=results[l].cal, irf_z=results[l].irf_z,
              irf_δ=results[l].irf_δ, mom=results[l].mom, dest_ann=d, label=l)

output_E = (bed  = pack(specs[1][1], specs[1][2]),
            code = pack(specs[2][1], specs[2][2]),
            bgm  = pack(specs[3][1], specs[3][2]))

serialize("irf_delta_target.jls", output_E)
println("Saved: irf_delta_target.jls")
