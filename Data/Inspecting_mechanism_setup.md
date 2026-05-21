# Inspecting the Mechanism — Design Notes
**Last updated:** May 21, 2026

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

### Comparison B — Endogenous Exit  (`run_solution_endog_exit.jl`) ✅ written and run
**Question:** What does the x_c amplification add beyond a fixed δ_e shock?

**Implementation:**
- Full baseline: `(TARGETS..., b_ratio=0.9, x_v=0.5)` — dest_end_frac=0.5, p_0=0.5, Xc_Y=0.10
- Exog exit: `(TARGETS..., dest_end_frac=0.0, p_0=0.0, Xc_Y=0.0, b_ratio=0.9, x_v=0.5)`
  — δ_e stays at same SS level (same dest_ann) but x_c is irrelevant
- Plot script: `plot_endog_exit_comparison.jl` → `mechanism_B_z_shock.pdf`, `mechanism_B_delta_shock.pdf`
- Calibrated parameters close: ϕ (0.659 vs 0.640), κ (0.058 vs 0.063) — comparison is clean

**IRF results (May 21, 2026):**

δ shock:
- Endog exit has SMALLER u response (peak 0.028 pp at h≈3) vs exog (peak 0.055 pp at h≈6)
- N trough: endog −0.05%, exog −0.10% (roughly half)
- exit_flow (δ_e×N): nearly identical in both models (~6.5% spike at h=1, rapid decay)
- N_e trough: endog −0.05%, exog −0.25% (larger entry depression in exog)

z shock:
- Endog exit has LARGER θ (1.2% vs 0.9% peak) and N (0.30% vs 0.25% peak)
- u and v responses similar in shape, endog slightly deeper u trough (−0.030 vs −0.022 pp)
- K and N_e responses nearly identical at impact

**Interpretation — two key findings:**

1. DOMINANT FACTOR for δ shock: The exogenous destruction parameter δbar is halved by
   construction (δbar_endog = 0.00326 vs δbar_exog = 0.00651) because dest_end_frac = 0.5
   routes half the SS destruction through x_c. Since δ_e = 1−(1−δ·δbar)·F(x_c), the δ shock
   fires through δ·δbar, giving half the direct δ_e impulse in the endog model. The smaller
   u/θ/N responses are primarily a mechanical exposure effect, not a dynamic x_c buffer story.

2. CALIBRATED ψ ≈ 0.014 in endog model: The continuation cost distribution is nearly degenerate
   (density f(x_c) → 0 as ψ → 0). F(x_c) barely responds to changes in profitability. The
   "x_c jump variable" dynamic amplification mechanism is quantitatively nearly inactive at this
   calibration. For ψ to matter, it would need to be of order 1 (Broer et al. 2025 use ψ = 1).

3. For z shock: larger endog responses come from lower f_e (20.3 vs 32.7) and lower δbar
   (less N drag in the LOM), not from the x_c channel (which is dormant with ψ ≈ 0.014).

**Design implication:** Comparison B as currently configured primarily illustrates the
EXPOSURE effect of dest_end_frac on the δ shock channel. To cleanly isolate the dynamic
x_c mechanism, either: (a) target ψ ≈ 1.0 directly in calibration, or (b) shock total δ_e
directly rather than the exogenous component δ. Current results are still valid for showing
that endogenous exit changes the propagation of δ shocks, but the mechanism is exposure-
based rather than dynamic-buffering-based.

**Note on dest_el:** Broer et al. (IER 2025) use separation elasticity ψ as a calibration
target (their preferred ψ = 1). Calibrated ψ = 0.014 here is far below their value, which
may explain the dormant x_c channel. Report implied dest_el as an outcome alongside ψ to
provide external validation context.

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

## Plot Design (settled — May 21, 2026)

Each comparison has its own pair of figures (z shock + δ shock), each a **4×2 portrait grid**
suitable for a single journal column. Variables chosen to trace the transmission chain top to bottom.

| Row | Left panel | Right panel | Units |
|-----|-----------|-------------|-------|
| 1 | u | v | pp levels deviation |
| 2 | θ | w_int | 100×log dev |
| 3 | N | N_e | 100×log dev |
| 4 | K (Comp. A) or exit_flow (Comp. B) | Q (Comp. A) or exit_flow (Comp. B) | 100×log dev |

**Why w_int not labor_prod:** labor_prod = μ × w_int; identical in log-dev since μ is constant.
w_int enters the JCC directly — the more mechanistically interpretable choice.

**Why Q in Comparison A:** K and Q side-by-side reveal the short-duration asset amplification
(same K response, much larger Q response under high-δ). Q ≡ e in log-dev with ξ_inv=1.

**Why exit_flow in Comparison B:** δ_e×N shows the total destruction flow directly. Under
exog exit, d log(exit_flow) = d log(N); under endog exit, δ_e also moves. The gap is the
x_c contribution (but note: at ψ ≈ 0.014, gap is small for δ shock).

**Color convention:**
- Baseline (solid black) — always the richer/fuller model
- Comparison A alternative: royalblue dashed   (high-δ level)
- Comparison B alternative: crimson dashed      (exogenous exit)
- Comparison C alternative: (reserved)          (no variety, future)

---

## Implementation Order

1. ✅ Run `run_solution_all_delta.jl` → `irf_all_delta.jls`
   Plot: `plot_mechanism_comparison.jl` → `mechanism_A_z_shock.pdf`, `mechanism_A_delta_shock.pdf`

2. ✅ Run `run_solution_endog_exit.jl` → `irf_endog_exit.jls`
   Plot: `plot_endog_exit_comparison.jl` → `mechanism_B_z_shock.pdf`, `mechanism_B_delta_shock.pdf`

3. ❌ Write + run `run_solution_no_variety.jl` → `irf_no_variety.jls`
   Use `run_solution_general_CES.jl` as template; set ζ=0 so ρ=1 regardless of N.
   Plot: extend `plot_mechanism_comparison.jl` or new `plot_variety_comparison.jl`

4. ❌ Consider whether to re-calibrate Comparison B targeting ψ ≈ 1.0 to activate the
   dynamic x_c mechanism (see calibration design note in Comparison B section above).

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
