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

    lendog = "Full baseline (endog. exit)"
    lexog  = "Exog. exit (δ_e = δ̄)"
    cs = [:black :crimson]
    ls = [:solid :dash]

    p = Plots.plot(
        layout         = (4, 2),
        size           = (700, 900),
        legend         = :topright,
        titlefontsize  = 10,
        tickfontsize   = 8,
        labelfontsize  = 9,
        legendfontsize = 7,
        left_margin    = 6Plots.mm,
        bottom_margin  = 3Plots.mm,
        right_margin   = 2Plots.mm,
    )

    # ── Row 1: Labor market outcomes ─────────────────────────────────────────

    # Panel 1 — u (pp dev from SS)
    plot!(p[1], t, [u_e u_x],
          label     = [lendog lexog],
          title     = L"u",
          ylabel    = "pp deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 2 — v (pp dev from SS)
    plot!(p[2], t, [v_e v_x],
          label     = [lendog lexog],
          title     = L"v",
          ylabel    = "pp deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 2: Price signals ─────────────────────────────────────────────────

    # Panel 3 — θ (market tightness)
    plot!(p[3], t, [irf_endog.θ irf_exog.θ],
          label     = [lendog lexog],
          title     = L"\theta",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 4 — w_int (marginal revenue product = labor_prod / μ)
    plot!(p[4], t, [irf_endog.w_int irf_exog.w_int],
          label     = [lendog lexog],
          title     = L"w_{\mathrm{int}}",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 3: Extensive margin ──────────────────────────────────────────────

    # Panel 5 — N (firm mass / product variety stock)
    plot!(p[5], t, [irf_endog.N irf_exog.N],
          label     = [lendog lexog],
          title     = L"N",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 6 — N_e (new varieties created; BGM entry margin)
    plot!(p[6], t, [irf_endog.N_e irf_exog.N_e],
          label     = [lendog lexog],
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
          label     = [lendog lexog],
          title     = L"K",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    # Panel 8 — exit flow δ_e × N
    plot!(p[8], t, [irf_endog.exit_flow irf_exog.exit_flow],
          label     = [lendog lexog],
          title     = L"\delta_e \times N",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    display(p)
    savefig(filename)
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
# Last updated: May 21, 2026
# =============================================================================
#
# ── What this comparison isolates ────────────────────────────────────────────
#
# Comparison A varied the LEVEL of δ_e (baseline δ̄ vs. high-δ = τ), holding
# the exit structure fixed (exogenous only). Comparison B holds the LEVEL of
# δ_e fixed at the same SS value and varies whether exit is endogenous or not:
#   Full baseline: δ_e = 1 − (1−δ·δ̄)·F(x_c),  dest_end_frac = 0.5
#   Exog. exit:    δ_e = δ·δ̄ (x_c irrelevant, p_0 = 0, dest_end_frac = 0)
# Both models target the same annual destruction rate dest_ann = 0.0754 (BED),
# so δ_e_SS ≈ 0.0065 in both. However, because dest_end_frac = 0.5, the
# EXOGENOUS component δbar is halved in the endog model (0.00326 vs 0.00651).
#
# ── Calibrated parameters (b_ratio=0.9, x_v=0.5) ───────────────────────────
#
# Key parameters that are similar across models:
#   ϕ:   0.659 (endog) vs 0.640 (exog)   — nearly identical; Nash wage not a confounder
#   κ:   0.058 (endog) vs 0.063 (exog)   — nearly identical; JCC not a confounder
#   z:   1.369 (endog) vs 1.375 (exog)   — same SS productivity
#
# Key parameters that differ structurally:
#   δbar: 0.00326 (endog) vs 0.00651 (exog)  — HALF as large by construction
#   ψ:    0.014  (endog) vs 1.5    (exog)     — very flat distribution in endog
#   f_m:  25.1   (endog) vs 1.0    (exog)     — scale of continuation cost distribution
#   f_e:  20.3   (endog) vs 32.7   (exog)     — lower sunk entry cost in endog
#
# ── δ shock: why endog exit shows SMALLER u, θ, N responses ─────────────────
#
# The IRFs show a counterintuitive result: the endogenous exit model has
# smaller labor market responses to the δ shock (~0.028 pp u peak at h=3)
# than the exogenous exit model (~0.055 pp at h=6). Two factors explain this:
#
# FACTOR 1 — HALVED SHOCK TRANSMISSION (dominant):
#   δ_e = 1 − (1−δ·δbar)·F(x_c). The δ shock fires through the δ·δbar term.
#   With δbar_endog = 0.00326 = δbar_exog/2, the same proportional δ shock
#   generates half the direct δ_e impulse in the endog model. This is a
#   mechanical consequence of dest_end_frac = 0.5 routing half the SS
#   destruction through x_c rather than the exogenous δ channel.
#   Result: the exog model receives twice the direct vacancy/firm destruction
#   impulse from the δ shock, leading to deeper N troughs (~0.10% vs ~0.05%)
#   and more prolonged u elevation.
#
# FACTOR 2 — NEARLY INACTIVE x_c MARGIN (ψ ≈ 0):
#   The continuation cost distribution shape parameter ψ = 0.014 in the endog
#   calibration (vs. ψ = 1.5 in the exog model, which is irrelevant since p_0=0).
#   The density at the threshold f(x_c) ∝ ψ × (x_c/f_m)^{ψ-1} → 0 when ψ → 0.
#   F(x_c) barely responds to changes in x_c with ψ this small. The "x_c buffer
#   mechanism" (endogenous exit absorbing part of the shock) is quantitatively
#   negligible at this calibration. The smaller u/N response is almost entirely
#   due to Factor 1.
#
# ROBUSTNESS: The direction (endog smaller) is mechanically robust as long as
#   dest_end_frac > 0, because that always implies δbar_endog < δbar_exog for
#   the same δ_e target. The MAGNITUDE of the difference scales with dest_end_frac.
#
# ── z shock: why endog exit shows LARGER θ, N responses ─────────────────────
#
# For the z shock, the δ shock transmission channel is irrelevant (δ doesn't
# move with z). The larger responses in the endog model reflect:
#   (a) Lower f_e (20.3 vs 32.7): same Q → more N_e → more N accumulation
#   (b) Lower δbar: the N LOM coefficient (1-δ_e) is the same at SS but the
#       exogenous drag on N is smaller (δbar = 0.00326 vs 0.00651), so N
#       accumulates faster from any given N_e impulse
#   (c) The ψ ≈ 0 x_c channel barely contributes to z amplification either
# Result: endog model N peaks at ~0.30% vs ~0.25%, θ at ~1.2% vs ~0.9%
#
# ── exit_flow (δ_e × N) panel interpretation ─────────────────────────────────
#
# δ shock: Both models show nearly identical exit_flow spikes (~6.5% at h=1).
#   Under exog exit: larger N drop × constant δ_e
#   Under endog exit: smaller N drop × similar δ_e (ψ≈0 means δ_e barely rises
#   through x_c, so total exit_flow ≈ same). The gap between models is small.
#
# z shock: Endog model exit_flow slightly larger (~0.25% peak vs ~0.20%) because
#   higher N accumulation × slightly lower δ_e (x_c rising reduces δ_e, but
#   barely with ψ≈0; the effect is mainly the N stock being larger).
#
# ── Design note for future calibrations ──────────────────────────────────────
#
# If the goal is to isolate the dynamic x_c mechanism, the calibration needs
# ψ to be larger (e.g., ψ = 1.0–1.5 as in Broer et al. 2025). At ψ ≈ 0,
# the endogenous exit margin exists on paper but is nearly dormant dynamically.
# Comparison B as currently configured primarily illustrates the EXPOSURE effect
# of dest_end_frac (half the destruction is shielded from the δ shock) rather
# than the dynamic x_c amplification story.
# =============================================================================
