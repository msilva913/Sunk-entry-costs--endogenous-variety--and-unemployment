# plot_xi_inv_comparison.jl
# =============================================================================
# Mechanism comparison figures: role of ENTRY COST ELASTICITY (Comparison D)
#
# Loads irf_xi_inv.jls (produced by run_solution_entry_elasticity.jl) and
# produces two publication-quality figures — one per shock (z, δ) — each
# showing a 4×2 panel grid (portrait, suitable for single journal column):
#
#   Row 1 — Labor market outcomes:
#     Panel 1: u      — pp level deviation from SS
#     Panel 2: v      — pp level deviation from SS
#
#   Row 2 — Price signals:
#     Panel 3: θ      — 100 × log deviation
#     Panel 4: w_int  — 100 × log deviation
#
#   Row 3 — Extensive margin (entry and firm stock):
#     Panel 5: N      — 100 × log deviation
#     Panel 6: e      — 100 × log deviation  (entrant vacancies; direct
#                        channel through which ξ_inv operates: lower ξ_inv
#                        → flat marginal cost → large e response to K)
#
#   Row 4 — Asset values:
#     Panel 7: K      — 100 × log deviation  (per-period vacancy dividend)
#     Panel 8: Q      — 100 × log deviation  (marginal cost of entry;
#                        with low ξ_inv, Q ≈ x_m constant → small Q response
#                        even as e surges, contrasting with baseline)
#
# Lines: Baseline (solid black, ξ_inv = 1.0) vs. Low-ξ (dashed blue, ξ_inv = 0.1).
# Key payoff: Row 3 shows that low ξ_inv amplifies both N and e responses
# to z shocks; Row 4 shows Q barely moves despite large e (flat cost schedule),
# whereas in baseline Q and e move together (Q = x_m·e^ξ_inv with ξ_inv = 1).
# =============================================================================

using Serialization, DataFrames, Plots, LaTeXStrings
pgfplotsx()

import Plots: default
default(linewidth=2, grid=true, fontfamily="Helvetica", framestyle=:box)

# ── Load serialized output ────────────────────────────────────────────────────
output = deserialize("irf_xi_inv.jls")
base   = output.base
low_ξ  = output.low_ξ

# ── Helper: log-dev to pp level deviation ─────────────────────────────────────
# irf_col = 100 × log(x_t / x_bar)  (output of simulate_model × 100)
# ss_val  = x_bar (fraction, e.g. u_bar ≈ 0.06)
# Returns (x_t − x_bar) × 100  in percentage points
function to_pp(irf_col, ss_val)
    return (exp.(irf_col ./ 100) .- 1) .* ss_val .* 100
end

# ── Main plotting function ────────────────────────────────────────────────────
function plot_mechanism_D(irf_base, irf_low, ss_base, ss_low, shock_label, filename;
                          flip_sign=false)

    T   = nrow(irf_base)
    t   = 0:(T-1)
    sgn = flip_sign ? -1 : 1

    irf_b = sgn == 1 ? irf_base : DataFrame(Dict(c => sgn .* irf_base[!, c] for c in names(irf_base)))
    irf_l = sgn == 1 ? irf_low  : DataFrame(Dict(c => sgn .* irf_low[!,  c] for c in names(irf_low)))

    u_b = to_pp(irf_b.u, ss_base.u)
    u_l = to_pp(irf_l.u, ss_low.u)
    v_b = to_pp(irf_b.v, ss_base.v)
    v_l = to_pp(irf_l.v, ss_low.v)

    lbase = L"Baseline ($\xi_{\mathrm{inv}} = 1.0$)"
    llow  = L"Low-$\xi$ ($\xi_{\mathrm{inv}} = 0.1$)"
    cs    = [:black :royalblue]
    ls    = [:solid :dash]

    p = Plots.plot(
        layout        = (4, 2),
        size          = (850, 900),
        titlefontsize = 10,
        tickfontsize  = 8,
        labelfontsize = 9,
        left_margin   = 6Plots.mm,
        bottom_margin = 4Plots.mm,
        right_margin  = 2Plots.mm,
    )

    # ── Row 1: Labor market outcomes ─────────────────────────────────────────

    plot!(p[1], t, [u_b u_l],
          label   = false,
          title   = L"u",
          ylabel  = "pp deviation",
          xlabel  = "",
          color   = cs, linestyle = ls)

    plot!(p[2], t, [v_b v_l],
          label          = [lbase llow],
          legend         = :topright,
          legendfontsize = 7,
          title          = L"v",
          ylabel         = "pp deviation",
          xlabel         = "",
          color          = cs, linestyle = ls)

    # ── Row 2: Price signals ─────────────────────────────────────────────────

    plot!(p[3], t, [irf_b.θ irf_l.θ],
          label   = false,
          title   = L"\theta",
          ylabel  = "% deviation",
          xlabel  = "",
          color   = cs, linestyle = ls)

    plot!(p[4], t, [irf_b.w_int irf_l.w_int],
          label   = false,
          title   = L"w_{\mathrm{int}}",
          ylabel  = "% deviation",
          xlabel  = "",
          color   = cs, linestyle = ls)

    # ── Row 3: Extensive margin ──────────────────────────────────────────────
    # e (entrant vacancies) is the direct channel: flat cost schedule → large e.
    # N is the cumulative effect: persistent entry → higher firm mass.

    plot!(p[5], t, [irf_b.N irf_l.N],
          label   = false,
          title   = L"N",
          ylabel  = "% deviation",
          xlabel  = "",
          color   = cs, linestyle = ls)

    plot!(p[6], t, [irf_b.e irf_l.e],
          label   = false,
          title   = L"e \;(\mathrm{entrant\;vacancies})",
          ylabel  = "% deviation",
          xlabel  = "",
          color   = cs, linestyle = ls)

    # ── Row 4: Asset values ──────────────────────────────────────────────────
    # K responds similarly across models (set by JCC); Q = x_m·e^ξ_inv so
    # with ξ_inv = 0.1, Q barely moves even as e surges — revealing that
    # the flat cost schedule decouples entry volume from posting cost.

    plot!(p[7], t, [irf_b.K irf_l.K],
          label   = false,
          title   = L"K",
          ylabel  = "% deviation",
          xlabel  = "months",
          color   = cs, linestyle = ls)

    plot!(p[8], t, [irf_b.Q irf_l.Q],
          label   = false,
          title   = L"Q \;(= x_m \cdot e^{\xi_{\mathrm{inv}}})",
          ylabel  = "% deviation",
          xlabel  = "months",
          color   = cs, linestyle = ls)

    savefig(p, filename)
    println("Saved: $filename")
    return p
end

# ── Produce figures ───────────────────────────────────────────────────────────
println("\nPlotting z shock...")
p_z = plot_mechanism_D(base.irf_z, low_ξ.irf_z, base.ss, low_ξ.ss,
                       "z", "mechanism_D_z_shock.pdf"; flip_sign=true)

println("\nPlotting δ shock...")
p_δ = plot_mechanism_D(base.irf_δ, low_ξ.irf_δ, base.ss, low_ξ.ss,
                       "δ", "mechanism_D_delta_shock.pdf")

println("\nDone. Figures: mechanism_D_z_shock.pdf, mechanism_D_delta_shock.pdf")

# ── Quick numerical summary ───────────────────────────────────────────────────
println("\n── Peak responses ──────────────────────────────────────────────────")
for (label, irf_b, irf_l, ss_b, ss_l) in [
        ("z shock", base.irf_z, low_ξ.irf_z, base.ss, low_ξ.ss),
        ("δ shock", base.irf_δ, low_ξ.irf_δ, base.ss, low_ξ.ss)]
    println("  $label:")
    u_b = to_pp(irf_b.u, ss_b.u); u_l = to_pp(irf_l.u, ss_l.u)
    v_b = to_pp(irf_b.v, ss_b.v); v_l = to_pp(irf_l.v, ss_l.v)
    println("    u peak:    baseline=$(round(maximum(abs.(u_b)), digits=4)) pp  " *
            "low-ξ=$(round(maximum(abs.(u_l)), digits=4)) pp")
    println("    v trough:  baseline=$(round(minimum(v_b),       digits=4)) pp  " *
            "low-ξ=$(round(minimum(v_l),       digits=4)) pp")
    println("    θ peak:    baseline=$(round(maximum(abs.(irf_b.θ)),  digits=3))%  " *
            "low-ξ=$(round(maximum(abs.(irf_l.θ)),  digits=3))%")
    println("    N peak:    baseline=$(round(maximum(abs.(irf_b.N)),  digits=3))%  " *
            "low-ξ=$(round(maximum(abs.(irf_l.N)),  digits=3))%")
    println("    e peak:    baseline=$(round(maximum(abs.(irf_b.e)),  digits=3))%  " *
            "low-ξ=$(round(maximum(abs.(irf_l.e)),  digits=3))%")
    println("    Q peak:    baseline=$(round(maximum(abs.(irf_b.Q)),  digits=3))%  " *
            "low-ξ=$(round(maximum(abs.(irf_l.Q)),  digits=3))%")
    println("    e/Q ratio: baseline=$(round(maximum(abs.(irf_b.e))/max(maximum(abs.(irf_b.Q)),1e-8), digits=2))  " *
            "low-ξ=$(round(maximum(abs.(irf_l.e))/max(maximum(abs.(irf_l.Q)),1e-8), digits=2))  " *
            "  <- entry-cost decoupling")
end
