# ags_persistence.jl
# =============================================================================
# Exploratory, not draft-bound. TRUE-SHOCK statistics (decisions.md S16): compare models with
# each other, not with the data. AGS side of the comparison in
# Programs baseline/Q2_persistence_comparison.jl, which reads this program's CSV: does the full
# model's business-formation block add persistence or comovement relative to AGS?
#
# AGS at AGS_TARGETS (Λ̄_r 0.25, p_r 0.8, ψ_r 1, x_v 0.5, ξ_inv 1, b/w 0.9, X/Y 1.5 %), with Block
# B′'s shock processes (part6b section 9). For each shock alone: σ(u), σ(v), ρ₁(u), ρ₁(v),
# cor(u, v), cor(u, labor productivity) on quarterly averages, HP-1600 (shock_moments in
# time_series_fun.jl), and the quarter in which u's impulse response peaks (peak_quarter).
#
# Output: ags_persistence.csv (and console). Run before Q2_persistence_comparison.jl.
# =============================================================================

include("ags_model.jl")
using CSV

const SHK_B = SHK                          # Block B′'s processes (ags_model.jl)
const SDS_B = (SHK_B.σz, SHK_B.σδ, SHK_B.σs)
const T_IRF_M = 3 * 20                      # months of impulse response (20 quarters)

rows = DataFrame()
for (case, tg) in (("AGS", AGS_TARGETS),)
    r = solve_ags(tg; shk = SHK_B)
    for (i, name) in enumerate(("z", "δ", "s"))
        sim = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, T_SM, shock_col(i, SDS_B[i]),
                                       r.SS, false, true), varnames)
        mom = shock_moments(model_cycles(sim, r.cal); prod = "z")
        ir = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, T_IRF_M + 1,
                                      shock_col(i, SDS_B[i]), r.SS, true, true), varnames)
        kq, val = peak_quarter(100 .* ir.u[2:end])
        push!(rows, (; case, shock = name, mom..., u_peak_q = kq, u_peak = val))
    end
end
CSV.write("ags_persistence.csv", rows)
println(rows)
