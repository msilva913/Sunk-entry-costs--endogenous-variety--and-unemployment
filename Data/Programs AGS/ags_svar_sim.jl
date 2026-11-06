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
# ags_svar_reps_<case>.csv, ags_svar_true_<case>.csv in this folder.
# =============================================================================

include("ags_model.jl")
using CSV

const CASES = (
    ("default",    AGS_TARGETS),                                   # x_v 0.5, ξ_inv 1
    ("compromise", merge(AGS_TARGETS, (x_v = 0.25, ξ_inv = 0.5))), # Notes/ags_equilibrium.md §7.6
)
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

"""Block B′ shock processes (ρ_m, σ_m) for z, δ, s from part6b's shock_calibration.csv."""
function bprime_shocks()
    c = CSV.read(joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                          "shock_calibration.csv"), DataFrame)
    get(name, col) = c[findfirst(==(name), c.series), col]
    return (ρz = get("z", :bprime_rho_m),     σz = get("z", :bprime_sigma_m),
            ρδ = get("delta", :bprime_rho_m), σδ = get("delta", :bprime_sigma_m),
            ρs = get("s", :bprime_rho_m),     σs = get("s", :bprime_sigma_m))
end
const SHK_B = bprime_shocks()
const SDS_B = (SHK_B.σz, SHK_B.σδ, SHK_B.σs)
println(@sprintf("Block B′ shocks (monthly): ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f", SHK_B...))
const T_LONG = 3 * N_REP * (N_FILTER + BURN)            # months
const H_Q = 12                              # quarters, as H_IRF in q3_data_svar.py
const SVAR_VARS = ["z", "s", "δ", "u", "v"]

"""Rename the model columns to the data CSV's names (δ → delta)."""
to_data_names(df) = rename(df[!, SVAR_VARS], "δ" => "delta")

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
        qd = DataFrame(avg, ["z", "s", "delta", "u", "v"])
        insertcols!(qd, 1, :shock => name, :h => 0:nrow(qd)-1)
        append!(tru, qd)
    end
    CSV.write("ags_svar_true_$case.csv", tru)
    println(@sprintf("%s: ϕ %.3f, %d months simulated, %d replications of %d quarters (filtered over %d).",
                     case, r.cal.ϕ, T_LONG, N_REP, N_KEEP, N_FILTER))
end
