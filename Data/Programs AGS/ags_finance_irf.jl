# ags_finance_irf.jl
# =============================================================================
# Step 2a of the finance plan (Nov 15, 2026): TRUE-SHOCK impulse responses of AGS with the
# financial friction (ags_finance_model.jl), model against model (decisions.md S16). The early
# warning: does the premium move by the data's order of magnitude after a δ shock (EBP +4 to +8 bp
# at quarters 1–2, q3_delta_finance.txt), and does u's response to δ build to a later peak?
#
# Cases: the payout cost κ_d on a grid (0 = funds move freely, so the friction is a level wedge;
# large = retention close to fixed), at b̄ = 0.5, χ = 1; and χ = 4 at κ_d = 1. AGS_TARGETS,
# Block B′ shock processes. Quarterly averages of the response to a one-SD shock, 100 × log
# deviation (prem: annualized spread, bp; f: gross job finding), quarters 0–12.
# Reference: AGS without the friction, run separately (ags_finance_reference.csv, written by
# this program's AGS branch: julia ags_finance_irf.jl ags).
#
# Output: ags_finance_irf.txt, ags_finance_irf.csv
# =============================================================================

using CSV
const NQ = 13
if isempty(ARGS)
    MODE = "fin"
else
    MODE = ARGS[1]
end

if MODE == "ags"
    include("ags_model.jl")
    r = solve_ags(AGS_TARGETS)
    rows = DataFrame()
    for (i, shk) in enumerate(("z", "δ", "s"))
        ir = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, 3NQ + 1, shock_col(i, SDS[i]),
                                      r.SS, true, true), varnames)[2:end, :]
        uq = [mean(ir.u[3k-2:3k]) for k in 1:NQ]
        fm = log.(1 .- r.cal.δ .* exp.(ir.δ)) .- log(1 - r.cal.δ) .+ (1 - r.cal.η_L) .* ir.θ
        fq = [mean(fm[3k-2:3k]) for k in 1:NQ]
        for k in 1:NQ
            push!(rows, (; shock = shk, h = k - 1, u = 100uq[k], f = 100fq[k]))
        end
    end
    CSV.write("ags_finance_reference.csv", rows)
    println("Saved ags_finance_reference.csv")
else
    include("ags_finance_model.jl")
    CASES = [("κ_d = $k", merge(FIN_TARGETS, (κ_d = k,))) for k in (0.0, 0.1, 1.0, 10.0, 100.0)]
    push!(CASES, ("χ = 4, κ_d = 1", merge(FIN_TARGETS, (χ = 4.0,))))
    VARS = (:u, :f, :v, :e, :W, :ebp, :λ)
    rows = DataFrame()
    lines = String[]
    say(s) = (push!(lines, s); println(s))
    cal0, st0, fin0 = calibrate_fin()
    say("AGS with the financial friction: true-shock responses (ags_finance_irf.jl). AGS_TARGETS, Block B′ shocks.")
    say(@sprintf("Steady state (b̄ %.2f, χ %.1f): ψ̄ %.5f monthly = %.0f bp a year; Ψ̄^m − 1 = %.5f; ϑ̄ = %.3f (payout/funds a month);",
                 FIN_TARGETS.b_bar, FIN_TARGETS.χ, cal0.ψbar, 12e4 * cal0.ψbar, cal0.ψbar * (1 + FIN_TARGETS.χ), cal0.ϑbar))
    say(@sprintf("W/Y = %.4f, B/Y = %.4f, Π/Y = %.4f, outlays/Y = %.4f (monthly; Y = z̄(1 − u)); agency cost/Y = %.5f.",
                 fin0.W / (cal0.zbar * (1 - st0.u)), fin0.B / (cal0.zbar * (1 - st0.u)),
                 fin0.Π / (cal0.zbar * (1 - st0.u)), fin0.O / (cal0.zbar * (1 - st0.u)),
                 cal0.ψbar * fin0.B / (cal0.zbar * (1 - st0.u))))
    say("Quarterly averages of the response to a one-SD shock; 100 × log deviation (ebp: annualized spread, bp).")
    say("Each cell: peak (quarter) | value at quarters 1, 4, 8.")
    for (lab, ft) in CASES
        r = solve_fin(AGS_TARGETS, ft)
        for (i, shk) in enumerate(("z", "δ", "s"))
            ir = irf_fin(r, i, 3NQ)
            ir.f = log.(1 .- r.cal.δ .* exp.(ir.δ)) .- log(1 - r.cal.δ) .+ (1 - r.cal.η_L) .* ir.θ
            ir.ebp = 12e4 .* r.cal.ψbar .* (exp.(ir.prem) .- 1)       # bp, annualized
            for v in VARS
                if v == :ebp
                    p = quarterly(ir[!, v], NQ)
                else
                    p = 100 .* quarterly(ir[!, v], NQ)
                end
                for k in 1:NQ
                    push!(rows, (; case = lab, shock = shk, var = String(v), h = k - 1, irf = p[k]))
                end
            end
        end
    end
    CSV.write("ags_finance_irf.csv", rows)
    if isfile("ags_finance_reference.csv")
        ref = CSV.read("ags_finance_reference.csv", DataFrame)
    else
        ref = nothing
    end
    function cell(p)
        k = argmax(abs.(p))
        return @sprintf("%+.3f (q%d) | %+.3f %+.3f %+.3f", p[k], k - 1, p[2], p[5], p[9])
    end
    for shk in ("δ", "z", "s")
        say("\n$shk shock")
        if !isnothing(ref)
            g = ref[ref.shock .== shk, :]
            say("  AGS reference   u " * cell(g.u) * "    f " * cell(g.f))
        end
        say("  " * rpad("case", 16) * join([rpad(String(v), 34) for v in VARS[1:4]]))
        for (lab, _) in CASES
            cells = String[]
            for v in VARS
                p = rows[(rows.case .== lab) .& (rows.shock .== shk) .& (rows.var .== String(v)), :irf]
                push!(cells, rpad(cell(p), 34))
            end
            say("  " * rpad(lab, 16) * join(cells[1:4]))
            say("  " * rpad("", 16) * join([rpad(String(v), 34) for v in VARS[5:end]]))
            say("  " * rpad("", 16) * join(cells[5:end]))
        end
    end
    open("ags_finance_irf.txt", "w") do io
        foreach(l -> println(io, l), lines)
    end
end
