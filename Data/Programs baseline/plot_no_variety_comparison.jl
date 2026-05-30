# plot_no_variety_comparison.jl
# =============================================================================
# Mechanism comparison C: role of VARIETY (N → ρ externality)
#
# Loads irf_no_variety.jls (produced by run_solution_no_variety.jl) and
# produces two publication-quality figures — one per shock (z, δ) — each
# showing a 4×2 panel grid (portrait, suitable for single journal column):
#
#   Row 1 — Labor market outcomes:
#     Panel 1: u           — pp level deviation from SS
#     Panel 2: v           — pp level deviation from SS
#
#   Row 2 — Price signals:
#     Panel 3: θ           — 100 × log deviation
#     Panel 4: w_int       — 100 × log deviation
#             (key panel: w_int = ρ·z/μ in baseline, = z/μ in no-variety;
#              the gap directly shows the N→ρ→w_int amplification channel)
#
#   Row 3 — Extensive margin:
#     Panel 5: N           — 100 × log deviation
#     Panel 6: N_e         — 100 × log deviation
#
#   Row 4 — Vacancy value chain:
#     Panel 7: K           — 100 × log deviation
#     Panel 8: Q           — 100 × log deviation
#             (Q = x_m·e; the gap vs. K reveals duration amplification,
#              but the K gap itself reflects ρ→w_int→surplus)
#
# Lines: Full baseline (solid black) vs. No-variety, ρ≡1 (dotted darkorange).
#
# Color convention across comparisons:
#   Baseline (black solid) is always the richer/fuller model specification.
#   Comparison A alternative: royalblue dashed  (high-δ level)
#   Comparison B alternative: crimson dashed    (exogenous exit)
#   Comparison C alternative: darkorange dotted (no variety, ρ≡1)
#
# Note: Panel 4 (w_int) is the diagnostic payoff for Comparison C. In the
# baseline, w_int = ρ(N)·z/μ so a z↑ or N↑ pushes w_int through two channels:
# directly via z and indirectly via ρ. With ρ≡1, only the direct z channel
# operates. The gap between baseline and no-variety w_int responses isolates
# the ρ(N) amplification of labor demand.
# =============================================================================

using Serialization, DataFrames, Plots, LaTeXStrings
pgfplotsx()

import Plots: default
default(linewidth=2, grid=true, fontfamily="Helvetica", framestyle=:box)

# ── Load serialized output ────────────────────────────────────────────────────
output_C = deserialize("irf_no_variety.jls")
base     = output_C.base
novar    = output_C.novar

# ── Helper: log-dev to pp level deviation ─────────────────────────────────────
function to_pp(irf_col, ss_val)
    return (exp.(irf_col ./ 100) .- 1) .* ss_val .* 100
end

# ── Main plotting function ────────────────────────────────────────────────────
function plot_mechanism_C(irf_base, irf_novar, ss_base, ss_novar, shock_label, filename;
                          flip_sign=false)

    T = nrow(irf_base)
    t = 0:(T-1)
    sgn = flip_sign ? -1 : 1

    irf_b  = sgn == 1 ? irf_base  : DataFrame(Dict(c => sgn .* irf_base[!, c]  for c in names(irf_base)))
    irf_nv = sgn == 1 ? irf_novar : DataFrame(Dict(c => sgn .* irf_novar[!, c] for c in names(irf_novar)))

    u_b  = to_pp(irf_b.u,  ss_base.u)
    u_nv = to_pp(irf_nv.u, ss_novar.u)
    v_b  = to_pp(irf_b.v,  ss_base.v)
    v_nv = to_pp(irf_nv.v, ss_novar.v)

    lbase  = L"Baseline ($\rho_t = N_t^{1/(\varepsilon-1)}$)"
    lnovar = L"No variety ($\rho_t \equiv 1$, $\zeta=0$)"
    cs = [:black :darkorange]
    ls = [:solid :dot]

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

    plot!(p[1], t, [u_b u_nv],
          label     = false,
          title     = L"u",
          ylabel    = "pp deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    plot!(p[2], t, [v_b v_nv],
          label          = [lbase lnovar],
          legend         = :topright,
          legendfontsize = 7,
          title          = L"v",
          ylabel         = "pp deviation",
          xlabel         = "",
          color          = cs, linestyle = ls)

    # ── Row 2: Price signals ─────────────────────────────────────────────────

    plot!(p[3], t, [irf_b.θ irf_nv.θ],
          label     = false,
          title     = L"\theta",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # w_int is the key diagnostic panel: gap = N→ρ→w_int amplification
    plot!(p[4], t, [irf_b.w_int irf_nv.w_int],
          label     = false,
          title     = L"w_{\mathrm{int}} \;(= \rho z/\mu)",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 3: Extensive margin ──────────────────────────────────────────────

    plot!(p[5], t, [irf_b.N irf_nv.N],
          label     = false,
          title     = L"N",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    plot!(p[6], t, [irf_b.N_e irf_nv.N_e],
          label     = false,
          title     = L"N_e",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 4: Vacancy value chain ───────────────────────────────────────────
    # K and Q: with ρ≡1, ν_f = f_e/μ is a constant, so the entry Euler is
    # attenuated. K still reflects the JCC surplus, but the ρ→surplus channel
    # is absent. Q = x_m·e amplifies K via duration D(δ_e); the gap shows
    # how much of Comparison A's duration effect passes through variety.

    plot!(p[7], t, [irf_b.K irf_nv.K],
          label     = false,
          title     = L"K",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    plot!(p[8], t, [irf_b.Q irf_nv.Q],
          label     = false,
          title     = L"Q \;(= x_m \cdot e)",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    display(p)
    savefig(p, filename)
    println("Saved: $filename")
    return p
end

# ── Produce figures ───────────────────────────────────────────────────────────
println("\nPlotting z shock...")
p_z = plot_mechanism_C(base.irf_z, novar.irf_z, base.ss, novar.ss,
                       "z", "mechanism_C_z_shock.pdf"; flip_sign=true)

println("\nPlotting δ shock...")
p_δ = plot_mechanism_C(base.irf_δ, novar.irf_δ, base.ss, novar.ss,
                       "δ", "mechanism_C_delta_shock.pdf")

println("\nDone. Figures: mechanism_C_z_shock.pdf, mechanism_C_delta_shock.pdf")

# ── Quick numerical summary ───────────────────────────────────────────────────
println("\n── Peak responses ──────────────────────────────────────────────────")
for (label, irf_b, irf_nv, ss_b, ss_nv) in [
        ("z shock",  base.irf_z, novar.irf_z, base.ss, novar.ss),
        ("δ shock",  base.irf_δ, novar.irf_δ, base.ss, novar.ss)]
    println("  $label:")
    u_b  = to_pp(irf_b.u,  ss_b.u);  u_nv = to_pp(irf_nv.u, ss_nv.u)
    v_b  = to_pp(irf_b.v,  ss_b.v);  v_nv = to_pp(irf_nv.v, ss_nv.v)
    println("    u peak:     base=$(round(maximum(abs.(u_b)),          digits=4)) pp  " *
            "novar=$(round(maximum(abs.(u_nv)),         digits=4)) pp")
    println("    v trough:   base=$(round(minimum(v_b),                digits=4)) pp  " *
            "novar=$(round(minimum(v_nv),               digits=4)) pp")
    println("    θ peak:     base=$(round(maximum(abs.(irf_b.θ)),      digits=3))%  " *
            "novar=$(round(maximum(abs.(irf_nv.θ)),     digits=3))%")
    println("    w_int peak: base=$(round(maximum(abs.(irf_b.w_int)),  digits=3))%  " *
            "novar=$(round(maximum(abs.(irf_nv.w_int)), digits=3))%")
    println("    N trough:   base=$(round(minimum(irf_b.N),            digits=3))%  " *
            "novar=$(round(minimum(irf_nv.N),           digits=3))%")
    println("    K peak:     base=$(round(maximum(abs.(irf_b.K)),      digits=3))%  " *
            "novar=$(round(maximum(abs.(irf_nv.K)),     digits=3))%")
    println("    Q peak:     base=$(round(maximum(abs.(irf_b.Q)),      digits=3))%  " *
            "novar=$(round(maximum(abs.(irf_nv.Q)),     digits=3))%")
end

# =============================================================================
# NOTES ON IRF SCALE AND SHAPE — Comparison C (variety externality)
# Last updated: May 30, 2026
# =============================================================================
#
# ── Comparison design ────────────────────────────────────────────────────────
#
# Both models share: same τ = 0.031, same ε = 4.3 (markup μ = 1.30 unchanged),
# same b_ratio = 0.9, dest_elast_target = 5 (endogenous exit active in both).
# Only difference: baseline has ρ_t = N_t^{1/(ε−1)}, no-variety has ρ_t ≡ 1.
#
# ── z shock: mechanism ───────────────────────────────────────────────────────
#
# Baseline: z↑ → w_int = ρ·z/μ rises via two channels:
#   (a) direct: z↑ alone raises w_int proportionally
#   (b) indirect: z↑ → N rises (entry) → ρ = N^{1/(ε−1)} rises → w_int further
# No-variety: only channel (a) operates. The w_int panel gap is purely (b).
# Because surplus (w_int − b) drives the JCC and Nash wage, the missing ρ
# amplification compresses θ, u, and v responses.
#
# ── δ shock: mechanism ───────────────────────────────────────────────────────
#
# δ↑ → N falls → ρ falls (baseline) → w_int = ρ·z/μ falls → surplus falls
# → JCC shifts down → additional θ compression and u amplification.
# No-variety: N↓ does not affect ρ. Any difference in u and v is purely from
# the x_c feedback (same in both since dest_elast_target = 5) and the LOM
# coefficients. The ρ channel provides a STABILIZING amplifier for δ shocks:
# N↓ → ρ↓ → w_int↓ → surplus↓ → vacancies↓ → additional u↑.
#
# ── w_int panel interpretation ───────────────────────────────────────────────
#
# The gap in the w_int panels is the cleanest single-panel read on the
# variety mechanism. Under ρ≡1, w_int = z_t/μ moves only with z; under
# baseline, w_int = ρ_t·z_t/μ moves with both z and N. Any δ-shock gap in
# w_int in the baseline is absent in no-variety, directly quantifying how
# much of the wage amplification is driven by the variety externality.
# =============================================================================
