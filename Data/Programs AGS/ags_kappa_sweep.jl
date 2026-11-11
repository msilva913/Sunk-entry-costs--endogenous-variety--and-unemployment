# ags_kappa_sweep.jl
# =============================================================================
# Exploratory, not draft-bound. TRUE-SHOCK statistics (each shock simulated alone): diagnostics
# only, not comparable with the data's SVAR numbers (decisions.md S16); the like-for-like
# comparison is Bartek analysis/q3_model_svar.py.
# How much does the per-hire cost κ amplify z in the AGS model,
# with every other ingredient unchanged?
#
# κ is set through x_v, the sunk share of recruiting cost: κ = (1 − x_v)/x_v · K/q, so x_v = 1 is
# κ = 0 and a lower x_v a larger κ. X/Y stays at 1.5 %, so a lower x_v moves recruiting cost from
# the sunk position to the per hire; it does not add cost. Everything else at AGS_TARGETS
# (Λ̄_r = 0.25, p_r = 0.8, ψ_r = 1), at ξ_inv = 1 (the default) and ξ_inv = 0.01.
#
# Mechanism (Pissarides 2009): the cost of a hire, κ + K/q, has a part that rises in a boom and a
# part that does not. In DMP the rising part is c/q (q falls). Here K also jumps, because Q jumps
# and is expected to fall back (a capital loss on the position). A fixed κ rises with neither, so a
# larger κ share weakens the damping. The fundamental surplus FS stays flat across the grid because
# the calibration re-pins ϕ to keep w = 1: the effect does not run through FS here.
#
# Reported per case: κ, ϕ, κ's share of the cost of a hire, FS; σ(u) with each shock alone and
# with all three; σ(θ)/σ(z) and σ(v) under the z shock; the z-shock IRF extremes and K's impact
# response; cor(u, v) with all shocks. Quarterly averages, HP-1600, T_SM months.
# Data, SVAR A (1992Q3–2019Q4, preferred ordering, long run): σ(u) by shock = √share × σ(u) in A's
# own window (10.65): z 6.32, δ 3.64, s 3.82; three shocks together 8.24. Ad hoc (Nov 5, 2026);
# q3_data_svar.py must produce these before they are cited.
#
# Output: ags_kappa_sweep.txt, ags_kappa_irf.png.
# =============================================================================

include("ags_model.jl")

const XV_GRID = (1.0, 0.75, 0.5, 0.25, 0.1, 0.05)
const XI_GRID = (1.0, 0.01)

"""σ(u), σ(v) and σ(θ)/σ(z) under the z shock alone (quarterly, HP)."""
function z_moments(r)
    sim = DataFrame(simulate_model((; model..., ne = 1), r.sol_mat, T_SM, shock_col(1, SDS[1]),
                                   r.SS, false, true), varnames)
    hu, hv, hθ, hz = hpq(sim, :u), hpq(sim, :v), hpq(sim, :θ), hpq(sim, :z)
    return (su = 100std(hu), sv = 100std(hv), sθz = std(hθ) / std(hz))
end

out = String[]
note(str) = (push!(out, str); println(str))
note("AGS κ sweep (ags_kappa_sweep.jl). κ through x_v; X/Y 1.5 %; Λ̄_r 0.25, p_r 0.8, ψ_r 1.")
note("Data, SVAR A window (ad hoc): σ(u) z 6.32, δ 3.64, s 3.82; three shocks 8.24.")
irf_store = Dict{Tuple{Float64,Float64},DataFrame}()
for ξ in XI_GRID
    note(@sprintf("\nξ_inv = %.2f", ξ))
    note(@sprintf("  %5s %7s %6s %6s %6s | %6s %6s %6s %6s | %7s %6s | %7s %4s %7s %7s | %6s",
                  "x_v", "κ", "ϕ", "κ shr", "FS", "σu z", "σu δ", "σu s", "σu all",
                  "σθ/σz", "σv z", "u min", "at", "θ max", "K imp", "cor"))
    for xv in XV_GRID
        r = nothing
        try
            r = solve_ags(merge(AGS_TARGETS, (x_v = xv, ξ_inv = ξ)))
        catch err
            note(@sprintf("  %5.2f failed: %s", xv, first(sprint(showerror, err), 90)))
            continue
        end
        st, cal = r.st, r.cal
        τ    = cal.δ + (1 - cal.δ) * cal.s
        FS   = (cal.zbar - cal.b - (cal.r + τ) * cal.κ - cal.s * st.R_sh) / cal.zbar
        zm   = z_moments(r)
        su_d = moments(r, [2]).su
        su_s = moments(r, [3]).su
        allm = moments(r, 1:3)
        ir   = irf(r, 1)
        irf_store[(ξ, xv)] = ir
        p_u, p_θ = ir[2:end, :u], ir[2:end, :θ]
        iu   = argmin(p_u)
        note(@sprintf("  %5.2f %7.4f %6.3f %6.2f %6.3f | %6.2f %6.2f %6.2f %6.2f | %7.2f %6.2f | %7.3f %4d %7.3f %7.2f | %+6.3f",
                      xv, cal.κ, cal.ϕ, cal.κ / (cal.κ + st.K / st.q), FS, zm.su, su_d, su_s,
                      allm.su, zm.sθz, zm.sv, p_u[iu], iu - 1, maximum(p_θ), ir[2, :K], allm.c))
    end
end
note("\nκ shr = κ/(κ + K/q), κ's share of the cost of a hire (= 1 − x_v). FS = (z̄ − b − (r+τ)κ − s𝓡)/z̄.")
note("σu z, δ, s: σ(u) with that shock alone; σu all: all three. u min, θ max: z-shock IRF extremes (%),")
note("'at' = month (0 = impact); K imp = K's impact response to z (%). cor = cor(u, v), all shocks.")
open("ags_kappa_sweep.txt", "w") do io
    foreach(l -> println(io, l), out)
end

# Figure: z-shock IRFs of u, v, θ, K by x_v; ξ_inv = 1 (top row), 0.01 (bottom row).
# Colors: the dataviz reference categorical palette, in fixed order (validated; the table above
# is the text view that its contrast warning requires).
const PAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6")
fig = plot(layout = (2, 4), size = (1400, 640), background_color = "#fcfcfb",
           titlefontsize = 10, tickfontsize = 7, legendfontsize = 7, gridalpha = 0.15,
           left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
for (row, ξ) in enumerate(XI_GRID), (col, v) in enumerate((:u, :v, :θ, :K))
    for (k, xv) in enumerate(XV_GRID)
        if haskey(irf_store, (ξ, xv))
            path = irf_store[(ξ, xv)][2:end, v]
            if row == 1 && col == 1
                lab = "x_v = $xv"
            else
                lab = ""
            end
            plot!(fig[row, col], 0:length(path)-1, path, lw = 2, color = PAL[k], label = lab)
        end
    end
    hline!(fig[row, col], [0], color = "#8a8984", lw = 1, ls = :dash, label = "")
    title!(fig[row, col], @sprintf("z shock: %s  (ξ_inv = %.2f)", v, ξ))
    if row == 2
        xlabel!(fig[row, col], "months")
    end
    if col == 1
        ylabel!(fig[row, col], "%")
    end
end
savefig(fig, "ags_kappa_irf.png")
