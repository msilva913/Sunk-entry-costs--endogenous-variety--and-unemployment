# ags_run_solution.jl
# =============================================================================
# Runs the AGS model (ags_model.jl) at AGS_TARGETS. Monthly.
#
# Sections:
#   1. The model at AGS_TARGETS: steady state, impulse responses of u, v, θ, Q, Λ_r.
#   2. Check: with nearly free vacancy creation (ξ_inv → 0) and linear utility (σ = 0), Q is
#      constant and discounting is constant, so job creation pins θ_t as a function of the
#      current shocks alone: θ does not load on the stocks Ẽ, Ṽ^u, Ṽ^s (the analog of the draft's
#      lem:theta; Notes/ags_equilibrium.md §3). With σ > 0 the discount factor depends on
#      consumption, which depends on the stocks, so the result holds only at σ = 0.
#   3. D15 gate 1, stage (a): σ(u) by shock and cor(u, v), quarterly averages, HP-1600, at the
#      Block B′'s shock processes (ags_model.jl SHK).
# Output: ags_run_solution.txt, ags_irf.png
# =============================================================================

include("ags_model.jl")

lines = String[]
say(str) = (push!(lines, str); println(str))

# =============================================================================
# 1. The model at AGS_TARGETS
# =============================================================================
sol = solve_ags(AGS_TARGETS)
say("AGS model with costly reactivation (ags_run_solution.jl). Monthly; responses in % (100 × log dev.).")
say(@sprintf("Targets: Λ̄_r %.2f, p_r %.2f, ψ_r %.2f, ξ_inv %.2f, x_v %.2f, b/w %.2f, X/Y %.3f, σ %.2f.",
             AGS_TARGETS.Λr_target, AGS_TARGETS.p_r, AGS_TARGETS.ψ_r, AGS_TARGETS.ξ_inv,
             AGS_TARGETS.x_v, AGS_TARGETS.b_ratio, AGS_TARGETS.X_Y, AGS_TARGETS.σ))
say(@sprintf("Steady state: u %.4f, θ %.4f, Q %.4f, Λ_r %.3f; ϕ %.4f, κ %.4f, z̄ %.4f, a %.3f.",
             sol.st.u, sol.st.θ, sol.st.Q, sol.st.Λ_r, sol.cal.ϕ, sol.cal.κ, sol.cal.zbar, sol.cal.a))
say(@sprintf("Shocks (monthly, Block B′): ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f.\n", SHK...))

say("Impulse responses to one-SD shocks (month 1 = impact; peak = largest |response|):")
say(@sprintf("  %-6s %-6s %8s %8s %8s %8s", "shock", "var", "impact", "peak", "at", "month 24"))
irfs = [irf(sol, i) for i in 1:3]
const IRF_VARS = (:u, :v, :e, :θ, :Q, :K, :Λ_r)
for (i, name) in enumerate(("z", "δ", "s")), v in IRF_VARS
    p = irfs[i][2:end, v]
    peak = argmax(abs.(p))
    say(@sprintf("  %-6s %-6s %8.3f %8.3f %8d %8.3f", name, v, p[1], p[peak], peak, p[24]))
end

# Figure: columns = shocks, rows = the shock's own path, then the variables. One series per
# panel, so no legend; each panel has its own y-scale (responses differ by orders of magnitude).
const SHOCK_STATE = (:z, :δ, :s)
const ROWS = (:shock, IRF_VARS...)
const LINE, ZERO, SURFACE = "#2a78d6", "#8a8984", "#fcfcfb"   # dataviz reference palette
plt = plot(layout = (length(ROWS), 3), size = (1050, 1650), legend = false,
           background_color = SURFACE, titlefontsize = 10, tickfontsize = 7,
           guidefontsize = 8, gridalpha = 0.15, framestyle = :axes, left_margin = 7Plots.mm)
for (j, name) in enumerate(("z", "δ", "s")), (i, v) in enumerate(ROWS)
    if v == :shock
        col, ttl = SHOCK_STATE[j], "$name shock: $name itself"
    else
        col, ttl = v, "$name shock: $v"
    end
    path = irfs[j][2:end, col]
    plot!(plt[i, j], 0:length(path)-1, path, title = ttl, lw = 2, color = LINE)
    if i == length(ROWS)
        xlabel!(plt[i, j], "months")
    end
    if j == 1
        ylabel!(plt[i, j], "%")
    end
    hline!(plt[i, j], [0], color = ZERO, lw = 1, ls = :dash)
end
savefig(plt, "ags_irf.png")

# =============================================================================
# 2. Check: nearly free vacancy creation ⇒ θ depends on the current shocks only
# =============================================================================
iθ = findfirst(==(:θ), Symbol.(y))
say("\nCheck (σ = 0): θ's loading on the stocks Ẽ, Ṽ^u, Ṽ^s (max |gx|) relative to its loading on z, as ξ_inv → 0.")
for (label, over) in (("costless", (p_r = 0.0, Λr_target = 1.0, σ = 0.0)), ("costly 0.25", (σ = 0.0,)))
    rel = map((1e-2, 1e-4, 1e-6)) do ξ
        g = solve_ags(merge(AGS_TARGETS, over, (ξ_inv = ξ,))).sol_mat.gx[iθ, :]
        maximum(abs, g[1:3]) / abs(g[4])
    end
    say(@sprintf("  %-12s ξ_inv 1e-2: %.1e   1e-4: %.1e   1e-6: %.1e", label, rel...))
    if !(rel[3] < 1e-2 && rel[2] < rel[1] / 10 && rel[3] < rel[2] / 10)
        error("θ's loading on the stocks does not vanish as ξ_inv → 0 ($label)")
    end
end
say("✓ θ is a function of the current shocks alone in the vacancy free entry limit.")

# =============================================================================
# 3. D15 gate 1, stage (a): σ(u) by shock, true shocks
# =============================================================================
say("\nD15 gate 1, stage (a). Simulated $(T_SM) months, quarterly averages, HP λ = $(Int(LAM)).")
say("TRUE-SHOCK statistics: not comparable with the data's SVAR numbers (decisions.md S16); the")
say("like-for-like comparison is Bartek analysis/q3_model_svar.py.")
say(@sprintf("  %-34s %7s %7s | %6s %6s %6s %6s | %5s %5s %5s", "case", "ϕ", "cor", "σu", "z",
             "δ", "s", "%z", "%δ", "%s"))
for (label, over) in (
        ("GS: costless, ξ_inv 1",           (p_r = 0.0, Λr_target = 1.0)),
        ("GS: costless, ξ_inv 0.01",        (p_r = 0.0, Λr_target = 1.0, ξ_inv = 0.01)),
        ("costly 0.25, ξ_inv 1",            (;)),
        ("costly 0.25, ξ_inv 0.25",         (ξ_inv = 0.25,)),
        ("costly 0.25, ξ_inv 0.01",         (ξ_inv = 0.01,)),
        ("costly 0.25, ξ_inv 0.01, x_v 1",  (ξ_inv = 0.01, x_v = 1.0)),
        ("costly 0.70, ξ_inv 0.01",         (Λr_target = 0.70, ξ_inv = 0.01)),
        ("bridge: 0.25, ξ_inv 0.01, X/Y 1.67%", (ξ_inv = 0.01, X_Y = X_Y_BRIDGE)))
    r = try
        solve_ags(merge(AGS_TARGETS, over))
    catch err
        say(@sprintf("  %-34s failed: %s", label, first(sprint(showerror, err), 80)))
        continue
    end
    all3 = moments(r, 1:3)
    one  = [moments(r, [i]).su for i in 1:3]
    sh   = 100 .* one .^ 2 ./ sum(one .^ 2)
    say(@sprintf("  %-34s %7.4f %+7.3f | %6.2f %6.2f %6.2f %6.2f | %5.1f %5.1f %5.1f", label,
                 r.cal.ϕ, all3.c, all3.su, one..., sh...))
end
say("σu: all three shocks; z, δ, s: σ(u) with that shock alone; %: shares of the sum of variances.")

open("ags_run_solution.txt", "w") do io
    foreach(l -> println(io, l), lines)
end
