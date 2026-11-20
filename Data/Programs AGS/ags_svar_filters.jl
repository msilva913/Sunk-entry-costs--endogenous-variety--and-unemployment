# ags_svar_filters.jl
# =============================================================================
# Model side of SVAR A by filter (HP or Hamilton) across reactivation rates: the overidentifying
# test through the SVAR. The reactivation rate is identified by the separation responses; the
# destruction response is then the test (draft Section 9). ags_reactivation_grid.jl runs the same
# grid on true shocks; this program runs it like-for-like with the data, under both filters.
#
# Cases: Λ̄_r ∈ {1 (costless), 0.9, 0.5, 0.25, 0.1, 0.03}, as ags_reactivation_grid.jl (costly cases
# p_r = 0.99); AGS_TARGETS otherwise. Filters:
#   hp:  HP-1600, as ags_svar_sim.jl (z, s, u, v over N_FILTER quarters, δ over the last N_KEEP);
#   ham: Hamilton (h = 8, p = 4; time_series_fun.jl), applied the same way; the first 11
#        quarters of each filtered series are missing, as in the data (SVAR window 1995Q2–2019Q4).
# Shock processes: matched for each filter on SVAR A's window
# (Bartek analysis/data/results/shock_calibration_window.csv; part6e_window_processes.py).
#
# Output (regenerated, not tracked): ags_svar_{reps,long,true}_<filter>_<case>.csv and the case
# list ags_svar_filters_cases.csv, read by Bartek analysis/q3_model_svar_filters.py.
# =============================================================================

include("ags_model.jl")
using CSV

const WINDOW_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                            "shock_calibration_window.csv")
const DATA_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                          "q3_data_svar_cycles_hp.csv")

function data_windows(path)
    d = CSV.read(path, DataFrame; missingstring = "")
    return nrow(d), count(!ismissing, d.delta)          # 276, 110
end
const N_FILTER, N_KEEP = data_windows(DATA_CSV)
const BURN, N_REP = 40, 200
const T_LONG = 3 * N_REP * (N_FILTER + BURN)
const H_Q = 12
const SVAR_VARS = ["z", "s", "δ", "u", "v", "f"]
to_data_names(df) = rename(df[!, SVAR_VARS], "δ" => "delta", "f" => "jf")

const FILTERS = [("hp", x -> hp_filter(x, LAM)), ("ham", x -> hamilton_filter(x))]
const GRID = vcat([("L1", merge(AGS_TARGETS, (p_r = 0.0, Λr_target = 1.0)))],
                  [(@sprintf("L%g", Λ), merge(AGS_TARGETS, (p_r = 0.99, Λr_target = Λ)))
                   for Λ in (0.9, 0.5, 0.25, 0.1, 0.03)])

"""The filter applied to every column of a quarterly frame (the long sample)."""
filter_frame(q, f) = DataFrame([c => f(Float64.(q[!, c])) for c in names(q)])

"""True responses to each one-SD shock, quarterly averages of shocks in months 1–3 of quarter 0."""
function true_irfs(r, sds)
    tru = DataFrame()
    n_m = 3 * (H_Q + 1)
    for (i, name) in enumerate(("z", "delta", "s"))
        path = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, n_m + 1,
                                        shock_col(i, sds[i]), r.SS, true, true), varnames)
        path = path[2:end, :]
        avg = nothing
        for lag in 0:2
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
    return tru
end

function run_case(fname, f, case, tg)
    shk = shock_processes(WINDOW_CSV; prefix = fname)
    sds = (shk.σz, shk.σδ, shk.σs)
    r = solve_ags(tg; shk = shk)
    eta_m = hcat([shock_col(i, sds[i]) for i in 1:3]...)
    sim = DataFrame(simulate_model((; model..., ne = 3), r.sol_mat, T_LONG, eta_m, r.SS,
                                   false, true), varnames)
    q = model_quarterly(sim, r.cal)
    CSV.write("ags_svar_long_$(fname)_$case.csv", to_data_names(filter_frame(q, f)))
    reps = DataFrame()
    for k in 1:N_REP
        last_q = k * (N_FILTER + BURN)
        w = window_filter(q[last_q-N_FILTER+1:last_q, :], N_FILTER, N_KEEP; short = ["δ"], filter = f)
        w = to_data_names(w)
        insertcols!(w, 1, :rep => k)
        append!(reps, w)
    end
    CSV.write("ags_svar_reps_$(fname)_$case.csv", reps)
    CSV.write("ags_svar_true_$(fname)_$case.csv", true_irfs(r, sds))
    println(@sprintf("%s %s: Λ̄_r target %.2f, ϕ %.3f; shocks ρz %.4f ρδ %.4f ρs %.4f",
                     fname, case, tg.Λr_target, r.cal.ϕ, shk.ρz, shk.ρδ, shk.ρs))
    return (filter = fname, case = case, lambda_r = tg.Λr_target, phi = r.cal.ϕ)
end

function main()
    manifest = DataFrame()
    for (fname, f) in FILTERS
        for (case, tg) in GRID
            push!(manifest, run_case(fname, f, case, tg))
        end
    end
    CSV.write("ags_svar_filters_cases.csv", manifest)
end

main()
