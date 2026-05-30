# plot_mechanism_comparison.jl
# =============================================================================
# Mechanism comparison figures: role of the LEVEL of δ (Comparison A)
#
# Loads irf_all_delta.jls (produced by run_solution_all_delta.jl) and
# produces two publication-quality figures — one per shock (z, δ) — each
# showing a 4×2 panel grid (portrait, suitable for single journal column):
#
#   Row 1 — Labor market outcomes:
#     Panel 1: u      — pp level deviation from SS
#     Panel 2: v      — pp level deviation from SS
#
#   Row 2 — Price signals:
#     Panel 3: θ      — 100 × log deviation
#     Panel 4: w_int  — 100 × log deviation  (= labor_prod/μ; redundant with
#                        labor_prod in log-dev since μ is constant, but more
#                        directly interpretable in the JCC/Nash wage chain)
#
#   Row 3 — Extensive margin (firm stock and entry):
#     Panel 5: N      — 100 × log deviation
#     Panel 6: N_e    — 100 × log deviation
#
#   Row 4 — Asset values (vacancy creation chain):
#     Panel 7: K      — 100 × log deviation  (per-period vacancy dividend)
#     Panel 8: Q      — 100 × log deviation  (present value of K; Q = x_m·e
#                        with ξ_inv=1, so d log Q ≡ d log e)
#
# Lines: Baseline (solid black) vs. High-δ (dashed blue).
# Row 4 is the diagnostic payoff: K responds similarly across models while Q
# amplifies more under high-δ, revealing the short-duration asset mechanism
# (see notes at bottom of file).
#
# Future: add Comparison B (no-endog-exit, dash-dot red) and
#         Comparison C (no-variety, dotted orange) lines once those scripts run.
# =============================================================================

using Serialization, DataFrames, Plots, LaTeXStrings
pgfplotsx()

import Plots: default
default(linewidth=2, grid=true, fontfamily="Helvetica", framestyle=:box)

# ── Load serialized output ────────────────────────────────────────────────────
output  = deserialize("irf_all_delta.jls")
base    = output.base
high_δ  = output.high_δ

# ── Helper: log-dev to pp level deviation ─────────────────────────────────────
# irf_col  = 100 × log(x_t / x_bar)   (output of simulate_model × 100)
# ss_val   = x_bar (fraction, e.g. u_bar ≈ 0.06)
# Returns  (x_t − x_bar) × 100  in percentage points
function to_pp(irf_col, ss_val)
    return (exp.(irf_col ./ 100) .- 1) .* ss_val .* 100
end

# ── Main plotting function ────────────────────────────────────────────────────
function plot_mechanism_A(irf_base, irf_high, ss_base, ss_high, shock_label, filename;
                          flip_sign=false)

    T = nrow(irf_base)
    t = 0:(T-1)
    sgn = flip_sign ? -1 : 1

    # Scale all IRF columns by sgn (negation for negative shock convention)
    irf_b = sgn == 1 ? irf_base : DataFrame(Dict(c => sgn .* irf_base[!, c] for c in names(irf_base)))
    irf_h = sgn == 1 ? irf_high : DataFrame(Dict(c => sgn .* irf_high[!, c] for c in names(irf_high)))

    # u and v: convert from 100×log-dev to pp level deviation (then apply sign)
    u_b = to_pp(irf_b.u, ss_base.u)
    u_h = to_pp(irf_h.u, ss_high.u)
    v_b = to_pp(irf_b.v, ss_base.v)
    v_h = to_pp(irf_h.v, ss_high.v)

    lbase = L"Baseline ($\delta_e = \bar{\delta}$)"
    lhigh = L"High-$\delta$ ($\delta_e = \tau$)"
    cs = [:black :royalblue]
    ls = [:solid :dash]

    # 4×2 grid. Do NOT set legend=false at top level — that overrides per-subplot
    # settings. Instead, pass legend=false to each plot! call individually.
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

    # Panel 1 — u (pp dev from SS); legend=false here, panel 2 carries it
    plot!(p[1], t, [u_b u_h],
          label     = false,
          title     = L"u",
          ylabel    = "pp deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 2 — v: shared legend lives here (top-right, both lines described)
    plot!(p[2], t, [v_b v_h],
          label          = [lbase lhigh],
          legend         = :topright,
          legendfontsize = 7,
          title          = L"v",
          ylabel         = "pp deviation",
          xlabel         = "",
          color          = cs, linestyle = ls)

    # ── Row 2: Price signals ─────────────────────────────────────────────────

    # Panel 3 — θ (market tightness)
    plot!(p[3], t, [irf_b.θ irf_h.θ],
          label     = false,
          title     = L"\theta",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 4 — w_int (marginal revenue product of labor = labor_prod / μ;
    #            identical to labor_prod in log-dev since markup μ is constant)
    plot!(p[4], t, [irf_b.w_int irf_h.w_int],
          label     = false,
          title     = L"w_{\mathrm{int}}",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 3: Extensive margin ──────────────────────────────────────────────

    # Panel 5 — N (firm mass / product variety stock)
    plot!(p[5], t, [irf_b.N irf_h.N],
          label     = false,
          title     = L"N",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # Panel 6 — N_e (new varieties created; BGM entry margin)
    plot!(p[6], t, [irf_b.N_e irf_h.N_e],
          label     = false,
          title     = L"N_e",
          ylabel    = "% deviation",
          xlabel    = "",
          color     = cs, linestyle = ls)

    # ── Row 4: Vacancy asset values ──────────────────────────────────────────
    # Key diagnostic: K responds similarly across models; Q amplifies more
    # under high-δ because Q is a shorter-duration asset (lower Q_SS = K_SS /
    # (1 − β(1−δ_e))), so the same ΔK is a larger fraction of Q_SS.
    # Since ξ_inv = 1: Q = x_m·e ⟹ d log Q ≡ d log e (entrant vacancies).

    # Panel 7 — K (per-period vacancy dividend; drives JCC)
    plot!(p[7], t, [irf_b.K irf_h.K],
          label     = false,
          title     = L"K",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    # Panel 8 — Q (present value of K; Q = x_m·e with ξ_inv = 1)
    plot!(p[8], t, [irf_b.Q irf_h.Q],
          label     = false,
          title     = L"Q \;(= x_m \cdot e)",
          ylabel    = "% deviation",
          xlabel    = "months",
          color     = cs, linestyle = ls)

    savefig(p, filename)
    println("Saved: $filename")
    return p
end

# ── Produce figures ───────────────────────────────────────────────────────────
println("\nPlotting z shock...")
p_z = plot_mechanism_A(base.irf_z, high_δ.irf_z, base.ss, high_δ.ss,
                       "z", "mechanism_A_z_shock.pdf"; flip_sign=true)

println("\nPlotting δ shock...")
p_δ = plot_mechanism_A(base.irf_δ, high_δ.irf_δ, base.ss, high_δ.ss,
                       "δ", "mechanism_A_delta_shock.pdf")

println("\nDone. Figures: mechanism_A_z_shock.pdf, mechanism_A_delta_shock.pdf")

# ── Quick numerical summary ───────────────────────────────────────────────────
println("\n── Peak responses ──────────────────────────────────────────────────")
for (label, irf_b, irf_h, ss_b, ss_h) in [
        ("z shock",  base.irf_z, high_δ.irf_z, base.ss, high_δ.ss),
        ("δ shock",  base.irf_δ, high_δ.irf_δ, base.ss, high_δ.ss)]
    println("  $label:")
    u_b = to_pp(irf_b.u, ss_b.u); u_h = to_pp(irf_h.u, ss_h.u)
    v_b = to_pp(irf_b.v, ss_b.v); v_h = to_pp(irf_h.v, ss_h.v)
    println("    u peak:    baseline=$(round(maximum(abs.(u_b)),  digits=4)) pp  " *
            "high-δ=$(round(maximum(abs.(u_h)),  digits=4)) pp")
    println("    v trough:  baseline=$(round(minimum(v_b),        digits=4)) pp  " *
            "high-δ=$(round(minimum(v_h),        digits=4)) pp")
    println("    θ peak:    baseline=$(round(maximum(abs.(irf_b.θ)),   digits=3))%  " *
            "high-δ=$(round(maximum(abs.(irf_h.θ)),   digits=3))%")
    println("    N trough:  baseline=$(round(minimum(irf_b.N),   digits=3))%  " *
            "high-δ=$(round(minimum(irf_h.N),   digits=3))%")
    println("    K peak:    baseline=$(round(maximum(abs.(irf_b.K)),   digits=3))%  " *
            "high-δ=$(round(maximum(abs.(irf_h.K)),   digits=3))%")
    println("    Q peak:    baseline=$(round(maximum(abs.(irf_b.Q)),   digits=3))%  " *
            "high-δ=$(round(maximum(abs.(irf_h.Q)),   digits=3))%")
    println("    Q/K ratio: baseline=$(round(maximum(abs.(irf_b.Q))/max(maximum(abs.(irf_b.K)),1e-8), digits=2))  " *
            "high-δ=$(round(maximum(abs.(irf_h.Q))/max(maximum(abs.(irf_h.K)),1e-8), digits=2))  " *
            "  ← duration amplification")
end

# =============================================================================
# NOTES ON IRF SCALE AND SHAPE — Comparison A (level of δ_e)
# =============================================================================
#
# ── Why u and v responses to z are small ─────────────────────────────────────
#
# Both panels show modest u and v movements to the technology shock. This is
# not a model failure but the Shimer (2005) unemployment volatility puzzle.
# The binding constraint is the ratio b/w_int ≈ b/w (since w_int ≈ w at the
# normalized SS with w = 1). The Hagedorn-Manovskii amplification of θ to z is
# approximately w_int/(w_int − b), which with b_ratio = 0.9 and w_int ≈ 1
# gives ≈ 10×. Combined with the matching function elasticity 1/(1−η_L) = 2.5,
# the θ amplification is ≈ 20–25×, moderate relative to data.
#
# Note: the relevant outside-option ratio is b/w_int, NOT b/labor_productivity.
# Labor productivity = w_int × μ (markup-inflated), so b/labor_prod ≈ 0.69
# understates the true tightness b/w_int ≈ 0.90. Nash bargaining (ϕ ≈ 0.5–0.9)
# then means wages absorb most of the z variation, leaving little for the
# vacancy margin.
#
# ── Why κ = 0 is intentional ─────────────────────────────────────────────────
#
# The calibration sets x_v = 1 (all recruiting cost is sunk, paid at match
# formation), so the flow matching cost per period is zero by construction.
# This is the Coles-Kelishomi timing: the entire cost of a match is the sunk
# entry cost f_e, not a rental per vacancy posted. κ = 0 slightly amplifies
# θ responses relative to κ > 0 (no per-period drag on the Nash wage).
#
# ── Why e responds MORE under high-δ to z (opposite of N_e) ─────────────────
#
# e (entrant vacancies) and N_e (new product varieties) are distinct margins:
#   e   is pinned by the vacancy value: Q = x_m × e  (f[8])
#   N_e is pinned by entry-sector labor: N_e = z·z̄·L_e/f_e  (f[17])
#
# e is the CK vacancy-posting margin; N_e is the BGM variety-creation margin.
#
# For e: Q is the present value of per-period vacancy dividends K, so
#   Q_SS = K_SS / (1 − β(1−δ_e)).
# Under high-δ (δ_e = 0.031), Q_SS ≈ K_SS/0.034 — a SHORT-DURATION asset.
# Under baseline (δ_e = 0.0065), Q_SS ≈ K_SS/0.010 — a LONG-DURATION asset.
# When z rises and K increases by ΔK, the log-deviation response is:
#   d log(Q) = d log(K) × (1 − β(1−δ_e)) / (1 − β(1−δ_e)·ρ_z)
# Under high-δ: 0.034/0.129 ≈ 0.26;  baseline: 0.010/0.107 ≈ 0.09.
# The same proportional K gain represents a much larger fraction of the smaller
# Q_SS under high-δ. Short-duration assets have large proportional price swings
# from the same flow improvement — analogous to a short-term vs. long-term bond.
#
# For N_e: driven by the firm value ν_f Bellman (f[4]). The present-value
# multiplier 1/(1 − β(1−δ_e)·ρ_z) is SMALLER under high-δ (denominator 0.129
# vs. 0.107), compressing the ν_f response to z and hence the L_e and N_e
# response.
#
# ── Why N responds more under high-δ even though N_e responds less ───────────
#
# The N law of motion (log-linearized) is:
#   n̂_{t+1} = (1 − δ_e) · n̂_t + δ_e · n̂_e,t
# The coefficient on the entry flow n̂_e is δ_e itself. Under high-δ,
# δ_e = 0.031 vs. 0.0065 in baseline (4.8× larger). The long-run N multiplier
# to a persistent N_e impulse is δ_e / (1 − (1−δ_e)·ρ_z):
#   baseline: 0.0065 / 0.104 ≈ 0.063 × n̂_e
#   high-δ:   0.031  / 0.126 ≈ 0.246 × n̂_e
# So even if n̂_e is smaller under high-δ, as long as it is not more than
# ~4× smaller, the cumulative N response is larger. This is a high-turnover
# economy: a given entry flow adds proportionally more to the stock because the
# replacement rate N_e/N = δ_e is much higher.
#
# ── δ-shock IRFs vs. z-shock IRFs ────────────────────────────────────────────
#
# The δ shock operates through a structurally different channel than z. A δ↑
# destroys firm-vacancy pairs directly (through the vacancy LOM), reducing N
# and θ simultaneously, independent of the Nash wage surplus. The ρ = N^{1/(ε−1)}
# channel then amplifies: N↓ → ρ↓ → w_int = ρ·z·z̄/μ falls → entry slows →
# persistent N depression. This mechanism does not require a large firm surplus
# (1−ϕ) to generate large u and v responses, unlike z shocks. Comparison A
# isolates how the LEVEL of δ_e (not just the shock size) alters the dynamics:
# higher δ_e at SS creates more simultaneous vacancy destruction per unit shock,
# amplifying both the u rise and the N trough relative to baseline.
# =============================================================================
