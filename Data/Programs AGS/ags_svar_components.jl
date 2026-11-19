# ags_svar_components.jl
# =============================================================================
# Model side of SVAR A when measured δ is a persistent component plus measurement noise
# (interpretation (a); findings.md "Shock processes: one AR(1) does not fit both filters").
#
# part6d (shock_components_profile.csv) fits BED deaths as a persistent monthly AR(1) plus
# independent monthly noise, jointly to the HP and Hamilton statistics. The data pin the noise
# share w only loosely (0.1–0.4), so the model is run at each w on the profile. In each case:
#   - the economy's δ shock is the persistent component: ρδ and σδ from the profile;
#   - the measured δ fed to the SVAR is the model's measured exit plus the noise (log points,
#     monthly, before quarterly averaging), which the economy never sees;
#   - z and s keep Block B′'s processes, as in ags_svar_sim.jl, so only δ changes.
# Everything else follows ags_svar_sim.jl (default calibration, same windows and filter):
# reps (200 data-length windows), long (whole run), true (responses to the persistent δ shock,
# not filtered, no noise).
#
# Output (regenerated, not tracked): ags_svar_{reps,long,true}_<case>.csv and the case list
# ags_svar_components_cases.csv, read by Bartek analysis/q3_model_svar_components.py.
# =============================================================================

include("ags_model.jl")
using CSV, Random

const PROFILE_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                             "shock_components_profile.csv")
const DATA_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                          "q3_data_svar_cycles_hp.csv")
const NOISE_SEED = 20261119

function data_windows(path)
    d = CSV.read(path, DataFrame; missingstring = "")
    return nrow(d), count(!ismissing, d.delta)          # 276, 110
end
const N_FILTER, N_KEEP = data_windows(DATA_CSV)
const BURN, N_REP = 40, 200                             # quarters, as ags_svar_sim.jl
const T_LONG = 3 * N_REP * (N_FILTER + BURN)            # months
const H_Q = 12
const SVAR_VARS = ["z", "s", "δ", "u", "v", "f"]
to_data_names(df) = rename(df[!, SVAR_VARS], "δ" => "delta", "f" => "jf")

prof = CSV.read(PROFILE_CSV, DataFrame)
prof = prof[prof.series .== "delta", :]
tg = AGS_TARGETS                                          # the "default" case of ags_svar_sim.jl
manifest = DataFrame()
for row in eachrow(prof)
    case = @sprintf("noise_w%.1f", row.w)
    shk = merge(SHK, (ρδ = row.rho_m, σδ = row.sigma_persistent_m))
    sds = (shk.σz, shk.σδ, shk.σs)
    println(@sprintf("%s: persistent δ ρ %.4f σ %.4f (half-life %.0f months), noise σ %.4f",
                     case, shk.ρδ, shk.σδ, row.half_life_m, row.sigma_noise_m))
    r = solve_ags(tg; shk = shk)
    eta_m = hcat([shock_col(i, sds[i]) for i in 1:3]...)
    sim = DataFrame(simulate_model((; model..., ne = 3), r.sol_mat, T_LONG, eta_m, r.SS,
                                   false, true), varnames)
    noise = row.sigma_noise_m .* randn(MersenneTwister(NOISE_SEED), nrow(sim))
    q = model_quarterly(sim, r.cal; δ_noise = noise)

    CSV.write("ags_svar_long_$case.csv", to_data_names(hp_frame(q)))
    reps = DataFrame()
    for k in 1:N_REP
        last_q = k * (N_FILTER + BURN)
        w = window_filter(q[last_q-N_FILTER+1:last_q, :], N_FILTER, N_KEEP; short = ["δ"])
        w = to_data_names(w)
        insertcols!(w, 1, :rep => k)
        append!(reps, w)
    end
    CSV.write("ags_svar_reps_$case.csv", reps)

    # true responses to a one-SD persistent δ shock (index 2), quarterly averages over the
    # three months of quarter 0, as in ags_svar_sim.jl; no noise
    n_m = 3 * (H_Q + 1)
    path = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, n_m + 1,
                                    shock_col(2, sds[2]), r.SS, true, true), varnames)
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
    tru = DataFrame(avg, ["z", "s", "delta", "u", "v", "jf"])
    insertcols!(tru, 1, :shock => "delta", :h => 0:nrow(tru)-1)
    CSV.write("ags_svar_true_$case.csv", tru)
    push!(manifest, (case = case, w = row.w, rho_delta = shk.ρδ, sigma_delta = shk.σδ,
                     half_life_m = row.half_life_m, sigma_noise = row.sigma_noise_m, phi = r.cal.ϕ))
end
CSV.write("ags_svar_components_cases.csv", manifest)
