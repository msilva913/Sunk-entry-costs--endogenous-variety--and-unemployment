# plot_endog_exit_comparison.jl
# =============================================================================
# Mechanism comparison B: role of ENDOGENOUS EXIT
#
# Loads irf_endog_exit.jls (produced by run_solution_endog_exit.jl) and
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
#
#   Row 3 — Extensive margin:
#     Panel 5: N           — 100 × log deviation
#     Panel 6: N_e         — 100 × log deviation
#
#   Row 4 — Vacancy value and exit flow:
#     Panel 7: K           — 100 × log deviation  (per-period vacancy dividend)
#     Panel 8: exit_flow   — 100 × log deviation  (δ_e × N; shows endogenous
#                            exit amplification directly — key Comparison B panel)
#
# Lines: Full baseline / endog. exit (solid black) vs.
#        Exogenous exit (dashed crimson).
#
# Color convention across comparisons:
#   Baseline (black solid) is always the richer/fuller model specification.
#   Comparison A alternative: royalblue dashed  (high-δ level)
#   Comparison B alternative: crimson dashed    (exogenous exit)
#   Comparison C alternative: (reserved)        (no variety, future)
#
# Note: exit_flow replaces Q from Comparison A. Under endogenous exit,
# δ_e moves with x_c (the continuation cost cutoff), so exit_flow = δ_e × N
# captures both the endogenous destruction rate and the firm stock response.
# When endogenous exit is off (exog model), δ_e is fixed, so exit_flow ≡ N
# in log-dev — the panel collapses to the N response, making the difference
# between rows 3 and 4 a direct visualization of the x_c amplification.
# =============================================================================

using Serialization, DataFrames, Plots, LaTeXStrings
pgfplotsx()

import Plots: default
default(linewidth=2, grid=true, fontfamily="Helvetica", framestyle=:box)

# ── Load serialized output ────────────────────────────────────────────────────
output_B = deserialize("irf_endog_exit.jls")
endog    = output_B.endog
exog     = output_B.exog

# ── Helper: log-dev to pp level deviation ─────────────────────────────────────
# irf_col  = 100 × log(x_t / x_bar)   (output of simulate_model × 100)
# ss_val   = x_bar (fraction, e.g. u_bar ≈ 0.07)
# Returns  (x_t − x_bar) × 100  in percentage points
function to_pp(irf_col, ss_val)
    return (exp.(irf_col ./ 100) .- 1) .* ss_val .* 100
end

# ── Main plotting function ────────────────────────────────────────────────────
function plot_mechanism_B(irf_endog, irf_exog, ss_endog, ss_exog, shock_label, filename)

    T = nrow(irf_endog)
    t = 0:(T-1)

    # u and v: convert from 100×log-dev to pp level deviation
    u_e = to_pp(irf_endog.u, ss_endog.u)
    u_x = to_pp(irf_exog.u,  ss_exog.u)
    v_e = to_pp(irf_endog.v, ss_endog.v)
    v_x = to_pp(irf_exog.v,  ss_exog.v)

    lendog = L"Baseline (endog. exit)"
    lexog  = L"Exog. exit ($\bar\delta$ fixed, $\tau$ fixed)"
    cs = [:black :crimson]
    ls = [:solid :dash]

    # 4×2 grid. Do NOT set legend=false at top level — that overrides per-subplot
    # settings. Instead, pass legend=false to each plot! call individually.
    p = Plots.plot(
        layout         = (4, 2),
        size           = (700, 900),
        titlefontsize  = 10,
        tickfontsize   = 8,
        labelfontsize  = 9,
        left_margin    = 6Plots.mm,
        bottom_margin  = 4Plots.mm,
        right_margin   = 2Plots.mm,
    )

    # ── Row 1: Labor market outcomes ─────────────────────────────────────────

    # Panel 1 — u (pp dev from SS); no legend here
    plot!(p[1], t, [u_e u_x],
          label     = false,
          title     = L"u",
          ylabel    = "pp deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 2 — v: shared legend lives here (topright, both lines described)
    plot!(p[2], t, [v_e v_x],
          label          = [lendog lexog],
          legend         = :topright,
          legendfontsize = 7,
          title          = L"v",
          ylabel         = "pp deviation",
          xlabel         = "",
          color          = cs, linestyle = ls)

    # ── Row 2: Price signals ─────────────────────────────────────────────────

    # Panel 3 — θ (market tightness)
    plot!(p[3], t, [irf_endog.θ irf_exog.θ],
          label     = false,
          title     = L"\theta",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 4 — w_int (marginal revenue product = labor_prod / μ)
    plot!(p[4], t, [irf_endog.w_int irf_exog.w_int],
          label     = false,
          title     = L"w_{\mathrm{int}}",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 3: Extensive margin ──────────────────────────────────────────────

    # Panel 5 — N (firm mass / product variety stock)
    plot!(p[5], t, [irf_endog.N irf_exog.N],
          label     = false,
          title     = L"N",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 6 — N_e (new varieties created; BGM entry margin)
    plot!(p[6], t, [irf_endog.N_e irf_exog.N_e],
          label     = false,
          title     = L"N_e",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 4: Vacancy value and exit flow ───────────────────────────────────
    # Key diagnostic for Comparison B:
    # - K: shows whether the endogenous exit margin affects the vacancy value
    # - exit_flow = δ_e × N: under endogenous exit, δ_e responds to x_c, so
    #   this panel shows amplification through the destruction margin directly.
    #   Under exogenous exit, δ_e is fixed → exit_flow ≡ N in log-dev.
    #   The gap between the two lines is therefore the pure x_c contribution.

    # Panel 7 — K (per-period vacancy dividend)
    plot!(p[7], t, [irf_endog.K irf_exog.K],
          label     = false,
          title     = L"K",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    # Panel 8 — exit flow δ_e × N
    plot!(p[8], t, [irf_endog.exit_flow irf_exog.exit_flow],
          label     = false,
          title     = L"\delta_e \times N",
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
p_z = plot_mechanism_B(endog.irf_z, exog.irf_z, endog.ss, exog.ss,
                       "z", "mechanism_B_z_shock.pdf")

println("\nPlotting δ shock...")
p_δ = plot_mechanism_B(endog.irf_δ, exog.irf_δ, endog.ss, exog.ss,
                       "δ", "mechanism_B_delta_shock.pdf")

println("\nDone. Figures: mechanism_B_z_shock.pdf, mechanism_B_delta_shock.pdf")

# ── Quick numerical summary ───────────────────────────────────────────────────
println("\n── Peak responses ──────────────────────────────────────────────────")
for (label, irf_e, irf_x, ss_e, ss_x) in [
        ("z shock",  endog.irf_z, exog.irf_z, endog.ss, exog.ss),
        ("δ shock",  endog.irf_δ, exog.irf_δ, endog.ss, exog.ss)]
    println("  $label:")
    u_e = to_pp(irf_e.u, ss_e.u); u_x = to_pp(irf_x.u, ss_x.u)
    v_e = to_pp(irf_e.v, ss_e.v); v_x = to_pp(irf_x.v, ss_x.v)
    println("    u peak:        endog=$(round(maximum(abs.(u_e)),         digits=4)) pp  " *
            "exog=$(round(maximum(abs.(u_x)),         digits=4)) pp")
    println("    v trough:      endog=$(round(minimum(v_e),               digits=4)) pp  " *
            "exog=$(round(minimum(v_x),               digits=4)) pp")
    println("    θ peak:        endog=$(round(maximum(abs.(irf_e.θ)),     digits=3))%  " *
            "exog=$(round(maximum(abs.(irf_x.θ)),     digits=3))%")
    println("    N trough:      endog=$(round(minimum(irf_e.N),           digits=3))%  " *
            "exog=$(round(minimum(irf_x.N),           digits=3))%")
    println("    exit_flow pk:  endog=$(round(maximum(abs.(irf_e.exit_flow)), digits=3))%  " *
            "exog=$(round(maximum(abs.(irf_x.exit_flow)), digits=3))%")
    println("    K peak:        endog=$(round(maximum(abs.(irf_e.K)),     digits=3))%  " *
            "exog=$(round(maximum(abs.(irf_x.K)),     digits=3))%")
end

# =============================================================================
# NOTES ON IRF SCALE AND SHAPE — Comparison B (endogenous vs. exogenous exit)
# Last updated: May 22, 2026
# =============================================================================
#
# ── Comparison design (updated) ───────────────────────────────────────────────
#
# The comparison holds fixed δbar (exogenous destruction component) and τ (total
# separation rate) across both specifications. This cleanly isolates the x_c
# amplification channel:
#
#   Full baseline (endog): δ_e ≈ 0.00653/month, δbar ≈ 0.00326 (dest_end_frac=0.5)
#                          τ = 0.031, s derived; x_c margin active (dest_elast=5)
#   Exog. exit:            δ_e = δbar_endog ≈ 0.00326 (halved by design)
#                          τ = 0.031 (same), s_exog > s_endog automatically
#                          x_c margin absent (p_0=0, dest_end_frac=0)
#
# Identification logic:
#   (a) Same δbar → same direct transmission of any δ shock through δ_e.
#       A proportional shock to δ generates the same first-period δ_e impulse
#       in both models; any IRF difference is due to the x_c channel, not
#       differential shock exposure.
#   (b) Same τ → approximately same SS u, v, θ. pp-deviation IRFs are directly
#       comparable. The higher s_exog is a derived consequence of keeping τ
#       fixed with lower δ_e, not an independent structural assumption.
#
# Previous design (same dest_ann → same δ_e, different δbar) is retained as an
# appendix exercise illustrating the pure "exposure effect" of dest_end_frac.
#
# ── Calibrated parameters (b_ratio=0.9, x_v=0.5, dest_elast_target=5) ───────
#
# Parameters that should be similar (design check):
#   δbar:  0.00326 (endog) = 0.00326 (exog)   ← equal by construction
#   τ:     0.031   (endog) ≈ 0.031  (exog)    ← equal by target (sep=0.031)
#   u,v,θ: approximately equal                 ← implied by same τ, f, q
#
# Parameters that differ structurally (as designed):
#   δ_e:   0.00653 (endog) vs 0.00326 (exog)  ← exog = δbar_endog
#   s:     s_endog < s_exog                    ← exog compensates via higher s
#   ψ:     ~0.033 (endog) vs 1.5 (exog placeholder, irrelevant with p_0=0)
#   dest_el: 5.0 (endog by target) vs 0.0 (exog, no endogenous exit)
#   f_e, z: may differ due to different π_s in each model
#
# ── δ shock: mechanism ────────────────────────────────────────────────────────
#
# With δbar equated, the same proportional δ shock generates the same direct
# first-period δ_e impulse in both models. Any difference in the IRFs reflects
# the x_c dynamic amplification in the endog model:
#   - δ↑ → δ_e↑ → x_c falls (fewer firms find continuation worthwhile) →
#     endogenous exit rises → additional N destruction beyond exog impulse.
# The magnitude of this amplification depends on dest_elast = ψ×ζ/(1−ζ) = 5.
# At the current calibration ζ≈0.993, ψ≈0.033, so the density f(x_c) at the
# threshold is still low — amplification visible but moderate.
# The gap in N and u responses between endog and exog lines is the pure x_c
# contribution: endog should show larger u, deeper N trough (amplification),
# or faster recovery depending on whether x_c buffers or amplifies.
#
# ── z shock: mechanism ───────────────────────────────────────────────────────
#
# The δ shock transmission channel is irrelevant for z. Differences reflect:
#   (a) In the endog model: z↑ → π_s rises → x_c rises → F(x_c) rises →
#       endogenous exit falls → δ_e falls → N accumulates faster (positive
#       amplification through x_c). Magnitude governed by dest_elast.
#   (b) Different f_e between models (from different π_s): affects entry margin.
# The z-shock comparison is less clean because f_e differs; interpret with care.
#
# ── exit_flow (δ_e × N) panel interpretation ─────────────────────────────────
#
# δ shock: With same δbar, the direct δ_e impulse is equated. Any gap between
#   endog and exog exit_flow at h=1 is the pure x_c contribution (endogenous
#   exit responding to the shock). At dest_elast=5, this gap may be visible.
#   At h>1, persistence difference reflects slower N recovery in endog model if
#   x_c amplifies, or faster if x_c partially mean-reverts.
#
# z shock: Endog exit_flow should show more movement than exog (x_c responds to
#   profitability: z↑ → x_c↑ → F(x_c)↑ → exit rate falls → exit_flow ≈ δ_e·N
#   falls; both δ_e falling and N rising). The exog exit_flow tracks N only.
# =============================================================================
