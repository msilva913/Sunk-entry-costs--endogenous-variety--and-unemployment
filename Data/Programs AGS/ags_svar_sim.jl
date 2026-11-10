# ags_svar_sim.jl
# =============================================================================
# Simulated data for the model side of SVAR A (Block B′; D13, D15). The SVAR itself runs in
# Python (Bartek analysis/q3_model_svar.py), with the same functions as the data SVAR
# (q3_data_svar.py), so the estimator is literally the same code on both sides.
#
# For each case, one simulation of the AGS model with all three shocks (T_LONG months), then:
#   - long:  quarterly averages of the whole run, HP-1600 on the whole run (the population SVAR);
#   - raw:   the same run, quarterly averages, not filtered (validation against the true IRFs);
#   - reps:  N_REP non-overlapping windows, each after BURN quarters, filtered as the data are
#            (window_filter, time_series_fun.jl): z, s, u, v HP-filtered over N_FILTER quarters
#            (the data's 1951Q1–2019Q4) and δ over the last N_KEEP only (BED starts with SVAR A's
#            window, 1992Q3–2019Q4); the last N_KEEP quarters are kept. Both lengths are read
#            from the data's cycle file. (simulate_model reseeds its generator on every call, so
#            independent samples come from one long run.)
#   - true:  the model's impulse responses to a one-SD shock in quarter 0, quarterly averages,
#            h = 0..H_Q, not filtered (for the SVAR validation step). A quarterly SVAR innovation
#            mixes shocks from the three months of the quarter, so the benchmark averages the
#            responses to a shock in month 1, 2 and 3 of quarter 0.
# Series, as in the data CSV (q3_data_svar_cycles_hp.csv): z, s, delta (measured exit δ_{t−1}),
# u, v; see ags_model.jl model_quarterly.
#
# Shock processes: Block B′'s (decisions.md D13), from part6b section 9 (shock_calibration.csv,
# columns bprime_rho_m, bprime_sigma_m): z and s re-matched on SVAR A's window, δ as in S14. The
# S14 processes (SHOCKS_BASE / SHK) fit 1951–2019 and are larger and more persistent for z and s.
#
# Output (regenerated, not tracked): ags_svar_long_<case>.csv, ags_svar_raw_<case>.csv,
# ags_svar_reps_<case>.csv, ags_svar_true_<case>.csv, and the case list ags_svar_cases.csv.
# =============================================================================

include("ags_model.jl")
using CSV

# Cases (Nov 10, 2026: exercise 1, z amplification through the per-hire cost at moderate
# inelasticity): ξ_inv ∈ {1, 0.5} × x_v ∈ {0.5, 0.25, 0.1, 0.05} (κ's share of the cost of a hire
# 50–95 %). Two keep their earlier names: "default" (ξ_inv 1, x_v 0.5, AGS_TARGETS) and
# "compromise" (ξ_inv 0.5, x_v 0.25; Notes/ags_equilibrium.md §7.6).
function case_name(ξ, xv)
    if ξ == 1.0 && xv == 0.5
        return "default"
    elseif ξ == 0.5 && xv == 0.25
        return "compromise"
    else
        return "xi$(ξ)_xv$(xv)"
    end
end
const CASES = [(case_name(ξ, xv), merge(AGS_TARGETS, (x_v = xv, ξ_inv = ξ)))
               for ξ in (1.0, 0.5) for xv in (0.5, 0.25, 0.1, 0.05)]
const DATA_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                          "q3_data_svar_cycles_hp.csv")

"""Filter span and SVAR window from the data: rows of the cycle file (1951Q1–2019Q4) and
the quarters with BED exit (SVAR A's window, 1992Q3–2019Q4)."""
function data_windows(path)
    d = CSV.read(path, DataFrame; missingstring = "")
    return nrow(d), count(!ismissing, d.delta)          # 276, 110
end
const N_FILTER, N_KEEP = data_windows(DATA_CSV)
const BURN, N_REP = 40, 200                             # quarters

const SHK_B = SHK                          # Block B′'s processes (ags_model.jl)
const SDS_B = (SHK_B.σz, SHK_B.σδ, SHK_B.σs)
println(@sprintf("Block B′ shocks (monthly): ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f", SHK_B...))
const T_LONG = 3 * N_REP * (N_FILTER + BURN)            # months
const H_Q = 12                              # quarters, as H_IRF in q3_data_svar.py
const SVAR_VARS = ["z", "s", "δ", "u", "v", "f"]   # f: the gross job-finding rate (persistence tests)

"""Rename the model columns to the data CSV's names (δ → delta, f → jf)."""
to_data_names(df) = rename(df[!, SVAR_VARS], "δ" => "delta", "f" => "jf")

manifest = DataFrame()
for (case, tg) in CASES
    r = solve_ags(tg; shk = SHK_B)
    eta_m = hcat([shock_col(i, SDS_B[i]) for i in 1:3]...)
    sim = DataFrame(simulate_model((; model..., ne = 3), r.sol_mat, T_LONG, eta_m, r.SS,
                                   false, true), varnames)
    q = model_quarterly(sim, r.cal)

    # long sample: HP on the whole run
    CSV.write("ags_svar_long_$case.csv", to_data_names(hp_frame(q)))
    CSV.write("ags_svar_raw_$case.csv", to_data_names(q))

    # replications: N_FILTER quarters after BURN, filtered as the data are, last N_KEEP kept
    reps = DataFrame()
    for k in 1:N_REP
        last_q = k * (N_FILTER + BURN)
        w = window_filter(q[last_q-N_FILTER+1:last_q, :], N_FILTER, N_KEEP; short = ["δ"])
        w = to_data_names(w)
        insertcols!(w, 1, :rep => k)
        append!(reps, w)
    end
    CSV.write("ags_svar_reps_$case.csv", reps)

    # true impulse responses, quarterly averages, not filtered: the average of the responses to a
    # shock in month 1, 2 and 3 of quarter 0 (simulate_model's impact month is row 2)
    tru = DataFrame()
    n_m = 3 * (H_Q + 1)
    for (i, name) in enumerate(("z", "delta", "s"))
        path = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, n_m + 1,
                                        shock_col(i, SDS_B[i]), r.SS, true, true), varnames)
        path = path[2:end, :]                      # month 1 = impact
        avg = nothing
        for lag in 0:2                             # shock in month 1 + lag of quarter 0
            shifted = vcat(DataFrame(zeros(lag, ncol(path)), names(path)), path[1:n_m-lag, :])
            qi = Matrix(to_data_names(model_quarterly(shifted, r.cal)))
            if lag == 0
                avg = qi ./ 3
            else
                avg .+= qi ./ 3
            end
        end
        qd = DataFrame(avg, ["z", "s", "delta", "u", "v", "jf"])
        insertcols!(qd, 1, :shock => name, :h => 0:nrow(qd)-1)
        append!(tru, qd)
    end
    CSV.write("ags_svar_true_$case.csv", tru)
    push!(manifest, (case = case, x_v = tg.x_v, xi_inv = tg.ξ_inv, sigma = tg.σ, phi = r.cal.ϕ,
                     kappa = r.cal.κ))
    println(@sprintf("%s: ϕ %.3f, %d months simulated, %d replications of %d quarters (filtered over %d).",
                     case, r.cal.ϕ, T_LONG, N_REP, N_KEEP, N_FILTER))
end
CSV.write("ags_svar_cases.csv", manifest)       # read by q3_model_svar.py
