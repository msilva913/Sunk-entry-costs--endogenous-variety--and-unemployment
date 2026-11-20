# ags_reactivation_grid.jl
# =============================================================================
# Does a low reactivation rate substitute for a propagation mechanism? TRUE-SHOCK statistics,
# model against model (decisions.md S16): AGS's responses of u and v to one-SD δ and s shocks
# as the steady-state reactivation rate Λ̄_r falls from 1 (costless, Gabrovski–Silva) toward 0
# (every separation retires its position: Coles and Kelishomi's destruction, applied to s).
# AGS_TARGETS otherwise (σ = 1, ξ_inv 1, x_v 0.5, X/Y 1.5 %), Block B′ shock processes; the costly
# cases use p_r = 0.99 so that Λ̄_r can go down to 0.03. Quarterly averages, 100 × log deviation.
# The like-for-like comparison with the data (through the SVAR) is still to be built.
#
# Output: ags_reactivation_grid.txt
# =============================================================================

include("ags_model.jl")

const NQ = 13
qa(x) = [mean(x[3k-2:3k]) for k in 1:NQ]
CASES = [("costless (Λ̄_r = 1)", merge(AGS_TARGETS, (p_r = 0.0, Λr_target = 1.0)))]
for Λ in (0.9, 0.5, 0.25, 0.1, 0.03)
    push!(CASES, ("Λ̄_r = $Λ", merge(AGS_TARGETS, (p_r = 0.99, Λr_target = Λ))))
end

lines = String[]
say(str) = (push!(lines, str); println(str))
say("AGS: responses of u and v to δ and s across reactivation rates (ags_reactivation_grid.jl).")
say("TRUE shocks (S16), one SD, quarterly averages, 100 × log deviation; q0 = the shock's quarter.")
say("Costly cases: p_r = 0.99, ψ_r = 1. Everything else at AGS_TARGETS; Block B′ shocks.")
say(@sprintf("\n%-20s %6s  %-14s %6s %6s %7s %7s   %-14s %6s %7s %7s", "case", "ϕ",
             "δ→u peak", "q4", "q8", "cum", "v q4", "s→u peak", "q4", "cum", "v q4"))
for (lab, tg) in CASES
    r = solve_ags(tg)
    cells = String[]
    for i in (2, 3)
        ir = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, 3NQ + 1, shock_col(i, SDS[i]),
                                      r.SS, true, true), varnames)[2:end, :]
        u, v = 100 .* qa(ir.u), 100 .* qa(ir.v)
        k = argmax(abs.(u))
        if i == 2
            push!(cells, @sprintf("%+.2f (q%d)     %+6.2f %+6.2f %+7.1f %+7.2f", u[k], k - 1, u[5], u[9], sum(u), v[5]))
        else
            push!(cells, @sprintf("%+.2f (q%d)     %+6.2f %+7.1f %+7.2f", u[k], k - 1, u[5], sum(u), v[5]))
        end
    end
    say(@sprintf("%-20s %6.3f  %s   %s", lab, r.cal.ϕ, cells[1], cells[2]))
end
say("\ncum: sum over quarters 0–12. v q4: v's response at quarter 4.")
open("ags_reactivation_grid.txt", "w") do io
    foreach(l -> println(io, l), lines)
end
