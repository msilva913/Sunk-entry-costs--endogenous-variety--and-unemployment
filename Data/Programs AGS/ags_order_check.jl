# ags_order_check.jl
# =============================================================================
# Exploratory. Does a second-order approximation change the AGS model's moments? Labor search
# models are non-linear (the matching function, the bounded reactivation rate, CRRA discounting),
# and Gabrovski and Silva report differences at higher order. TRUE-SHOCK statistics: first order
# against second order, model against model.
#
# Run twice, once per order (the toolkit's derivative functions are global, so each order needs
# its own process), from this folder:
#     AGS_ORDER=1 julia --project="../Programs baseline" ags_order_check.jl
#     AGS_ORDER=2 julia --project="../Programs baseline" ags_order_check.jl
# Each run writes ags_order_<k>.csv; when both exist, the run also writes the comparison,
# ags_order_check.txt.
#
# Second order is pruned (the toolkit's simulate_model). Its constant risk correction shifts means
# (reported as mean u, the average log deviation of unemployment from steady state, unfiltered);
# the HP filter removes it from the cyclical moments, so those differ only through the dynamics.
# Cases at σ = 1, Block B′'s shock processes (part6b section 9), quarterly averages, HP-1600:
# all shocks and each shock alone. Moments: σ(u), σ(v), ρ₁(u), cor(u, v), cor(u, z)
# (shock_moments), the skewness of u and v, and mean u.
# =============================================================================

include("ags_model.jl")
using CSV, StatsBase

const DATA_CSV  = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results", "q3_data_svar_cycles_hp.csv")
const SHK_B = SHK                          # Block B′'s processes (ags_model.jl)
const SDS_B = (SHK_B.σz, SHK_B.σδ, SHK_B.σs)
const CASES = (
    ("default (x_v 0.5, ξ_inv 1)",     AGS_TARGETS),
    ("compromise (x_v 0.25, ξ_inv 0.5)", merge(AGS_TARGETS, (x_v = 0.25, ξ_inv = 0.5))),
    ("high κ (x_v 0.1, ξ_inv 1)",      merge(AGS_TARGETS, (x_v = 0.1,))),
    ("elastic (x_v 0.25, ξ_inv 0.03)", merge(AGS_TARGETS, (x_v = 0.25, ξ_inv = 0.03))),
)
const SETS = (("all", [1, 2, 3]), ("z", [1]), ("δ", [2]), ("s", [3]))

rows = DataFrame()
for (case, tg) in CASES
    r = solve_ags(tg; shk = SHK_B)
    for (lab, cols) in SETS
        eta_m = hcat([shock_col(i, SDS_B[i]) for i in cols]...)
        sim = DataFrame(simulate_model((; model..., ne = length(cols)), r.sol_mat, T_SM, eta_m,
                                       r.SS, false, true), varnames)
        cyc = model_cycles(sim, r.cal)
        mom = shock_moments(cyc; prod = "z")
        push!(rows, (; order = flag_order, case, shocks = lab, mom...,
                     skew_u = skewness(cyc.u), skew_v = skewness(cyc.v), mean_u = 100mean(sim.u)))
    end
    println("order $flag_order: $case done")
end
CSV.write("ags_order_$(flag_order).csv", rows)

# Comparison, once both orders exist
if isfile("ags_order_1.csv") && isfile("ags_order_2.csv")
    a = CSV.read("ags_order_1.csv", DataFrame)
    b = CSV.read("ags_order_2.csv", DataFrame)
    d = CSV.read(DATA_CSV, DataFrame; missingstring = "")
    out = String[]
    note(str) = (push!(out, str); println(str))
    note("AGS: first against second order (ags_order_check.jl). TRUE shocks, σ = 1, Block B′ processes;")
    note("quarterly averages, HP-1600. Each cell: first order / second order.")
    note(@sprintf("Data skewness of the HP cycle (1951Q1–2019Q4): u %+.2f, v %+.2f.",
                  skewness(collect(skipmissing(d.u))), skewness(collect(skipmissing(d.v)))))
    for (case, _) in CASES
        note("\n$case")
        note(@sprintf("  %-6s %-15s %-15s %-15s %-15s %-15s %-15s %-15s", "shocks", "σ(u)", "ρ₁(u)",
                      "cor(u,v)", "cor(u,z)", "skew(u)", "skew(v)", "mean u, %"))
        for (lab, _) in SETS
            x = a[(a.case .== case) .& (a.shocks .== lab), :][1, :]
            y = b[(b.case .== case) .& (b.shocks .== lab), :][1, :]
            pair(f, v1, v2) = rpad(string(f(v1), " / ", f(v2)), 15)
            f2(v) = @sprintf("%.2f", v)
            f3(v) = @sprintf("%+.3f", v)
            note("  " * rpad(lab, 7) * pair(f2, x.σu, y.σu) * " " * pair(f3, x.ρu, y.ρu) * " " *
                 pair(f3, x.cuv, y.cuv) * " " * pair(f3, x.cup, y.cup) * " " *
                 pair(f3, x.skew_u, y.skew_u) * " " * pair(f3, x.skew_v, y.skew_v) * " " *
                 pair(f3, x.mean_u, y.mean_u))
        end
    end
    note("\nmean u: average log deviation of u from its steady state (×100), before filtering; zero at")
    note("first order by construction, the risk correction at second order. cor(u,z) NaN: z does not move.")
    open("ags_order_check.txt", "w") do io
        foreach(l -> println(io, l), out)
    end
end
