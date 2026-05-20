# Inspecting the Mechanism — Design Notes
**Last updated:** May 20, 2026 (evening)

---

## Goal

The mechanism section answers the skeptic's question: **what does this model do that neither
Coles-Kelishomi (CK) nor augmented Gabrovski-Silva (AGS) can?**

The baseline model has two innovations over AGS (which itself is GS + κ > 0 + σ > 0):
1. **Variety effects** — ρ = N^{1/(ε-1)} links firm mass to the JCC via w_int = ρ·z/μ
2. **Endogenous exit** — continuation cost cutoff x_c amplifies destruction shocks

And one innovation over CK:
3. **δ < τ (level of δ)** — destruction is rare but permanent; CK implicitly sets δ_e = τ

The section shows all three channels via model counterfactuals, for **both the z shock and
the δ shock**. GS (JEDC 2025) already established the δ-level mechanism for the z shock
and (u, labor_prod) comovement — this section extends that result and adds the δ shock.

---

## Three Comparisons

### Comparison A — Level of δ  (`run_solution_all_delta.jl`) ✅ written
**Question:** What happens to IRF persistence when δ_e = τ (CK timing) vs. δ_e = δ̄ ≪ τ?

**Implementation:**
- Baseline: `dest_ann = 0.0754`, `dest_end_frac = 0.5`
- High-δ: `dest_ann` s.t. monthly δ_e = sep = 0.031, `dest_end_frac = 0.0`, `p_0 = 0.0`
- All other targets unchanged → same SS u, v, θ
- s = 0 follows automatically in Stage 1 of `calibrate_shares` (τ = δ_e → s = (τ-δ_e)/(1-δ_e) = 0)

**Mechanism:** At high δ_e, both LOMs mean-revert quickly each period — vacancy survival
is (1-δ_e) ≈ 97% vs. 99% in baseline, firm survival likewise. A δ shock is therefore
less persistent: the system rebounds fast because normal turnover is already high.
At low δ_e, destruction is rare; when it hits, N and v recover slowly because the entry
flow N_e = δ_e·N/(1-δ_e) is small. This is directly visible in the LOMs:
  - Firm LOM: N_{t+1} = (1-δ_e)(N + N_e) — high δ_e compresses N fast
  - Vacancy LOM: v_pret_{t+1} = (1-δ_e)[...] — high δ_e kills the vacancy stock each period

**Key note:** This is NOT the endogenous exit comparison. p_0 = 0 and dest_end_frac = 0
eliminates the endogenous margin entirely. The comparison is purely about the level of
exogenous destruction in the LOMs — replicating CK timing without CK's other features
(κ = 0, σ = 0), which would confound the comparison.

**Feasibility warning:** With δ_e = τ = 3.1%/month, Stage 2 profit share may bind.
If calibration fails, lower Xc_Y from 0.10 to 0.07 in `targets_high_δ`.

---

### Comparison B — Endogenous Exit  (`run_solution_no_endog_exit.jl`) ❌ to write
**Question:** What does the x_c amplification add beyond a fixed δ_e shock?

**Implementation:**
- Baseline: as above
- No-endog-exit: `(TARGETS..., p_0=0.0, dest_end_frac=0.0)` — δ_e stays at baseline
  level but is now fully exogenous (Λ = 1 always, x_c irrelevant)
- This is a re-parameterization of the baseline; same SS u, v, θ, δ_e by construction

**Mechanism:** After a δ shock, endogenous exit amplifies because falling profits push
x_c down → more marginal firms exit → δ_e rises above the exogenous shock alone. This
amplification is absent with p_0 = 0. Key visible effect: exit_flow = δ_e×N should show
a larger and more persistent spike in baseline vs. no-endog-exit, even for the same δ shock.

**Note on dest_el:** Broer et al. (IER 2025) use separation elasticity ψ as a calibration
target (their preferred ψ = 1). Your `dest_el = ψ_shape × ζ/(1-ζ)` is the analogous object.
Rather than targeting dest_el directly (infeasibility risk), report it as an outcome:
baseline dest_el vs. 0 in no-endog-exit. This provides external validation via Broer.

---

### Comparison C — Variety Effects  (`run_solution_no_variety.jl`) ❌ to write
**Question:** What does the N → ρ → w_int → JCC channel add?

**Implementation:**
- Use `run_solution_general_CES.jl` as template (already has ζ parameterization)
- No-variety: set ζ = 0 so ρ = N^0 = 1 regardless of N
- Recalibrate to same targets

**Mechanism:** After a δ shock, N falls → ρ falls → w_int = ρ·z/μ falls → JCC tightens →
vacancy creation suppressed. This is the channel that drives simultaneous u↑ and v↓.
With ζ = 0, the JCC no longer tightens through the variety channel, so vacancies should
fall less. This directly answers "what can your model do that AGS cannot?"

For the z shock, variety effects work symmetrically: z↑ → profits rise → N_e rises →
N rises → ρ rises → w_int rises further → additional vacancy expansion. Shutting ζ = 0
dampens this second-round amplification.

---

## Plot Design (settled)

**Layout:** Two figures — one for z shock, one for δ shock. Each figure is a **4×2 grid**
(7 panels, one blank or use for legend):

| Panel | Variable | Units |
|-------|----------|-------|
| 1 | u (unemployment rate) | pp levels deviation |
| 2 | v (vacancy rate) | pp levels deviation |
| 3 | θ (market tightness) | 100×log dev |
| 4 | labor_prod | 100×log dev |
| 5 | N (firm mass) | 100×log dev |
| 6 | N_e (entry flow) | 100×log dev |
| 7 | exit_flow = δ_e×N | 100×log dev |

**Why u and v in levels:** All variants share the same steady-state u and v (by calibration
design), so pp deviations are directly comparable across lines and directly comparable to
empirical IRFs in Section 4.

**Why exit_flow not δ_e alone:** δ_e×N is the exit flow — symmetric to N_e as entry flow,
directly corresponds to the BED Deaths instrument used empirically. Plots both the rate
and the mass effect simultaneously. δ_e in levels can be added as a small supplementary
panel for Comparison A (where δ_e moves a lot) if needed.

**Lines per panel (4 total):**
- Baseline (solid black)
- High-δ / Comparison A (dashed blue) — level-of-δ channel
- No-endog-exit / Comparison B (dash-dot red) — endogenous exit channel
- No-variety / Comparison C (dotted orange) — variety channel

---

## Implementation Order

1. Run `run_solution_all_delta.jl` → debug → serialize `irf_all_delta.jls`
2. Write + run `run_solution_no_endog_exit.jl` → serialize `irf_no_endog_exit.jls`
3. Write + run `run_solution_no_variety.jl` (from `run_solution_general_CES.jl`) → serialize `irf_no_variety.jls`
4. Write `plot_mechanism_comparison.jl` loading all four IRF files → produce figures

---

## Key Decisions Made (May 20, 2026)

- **dest_end_frac vs. dest_el as calibration target:** Keep dest_end_frac. Feasibility of
  dest_el cannot be verified in advance since it depends jointly on ε, r, δ_e, Xc_Y.
  Report implied dest_el as outcome for comparison with Broer et al. (2025).

- **CK comparison:** Do NOT run full CK model (κ=0, σ=0 would confound). Instead,
  Comparison A (high-δ with κ>0, σ>0) isolates the δ-level channel cleanly. Describe in
  paper as "replicating the CK timing assumption while holding fixed matching costs and
  preferences."

- **s shock:** Not included in mechanism figures. Focus is on z and δ shocks.

- **Horizon:** T_IR = 60 months (5 years). Long enough to show persistence differences
  across Comparison A variants. May extend to 80 for δ shock if convergence is slow.
