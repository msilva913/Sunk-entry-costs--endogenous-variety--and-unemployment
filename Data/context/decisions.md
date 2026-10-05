# Decision Register
**Last updated:** October 5, 2026 (S12 rewritten after the Prop. 5 restatement) · **Branch:** `costly_vacancy_reposting`

One entry = one decision: a *choice*, not a task. Tasks live in
[`pending_tasks.md`](pending_tasks.md). When a decision is settled, move it to §Settled with
the date and the reason; the draft's prose depends on knowing *why*.

Status key: 🔴 blocks estimation · 🟠 blocks a paper section · 🟡 improves the paper

| # | Decision | Status |
|---|---|---|
| D2 | PATH A or PATH B as the calibration baseline | 🔴 open, recommendation PATH B |
| D3 | How β̂ and β(θ) are made comparable | 🔴 open, recommendation shape-only baseline; revisit after E9 |
| D4 | Does the LD→v IRF enter Block B? | 🟠 open, reopened by D11 |
| D10 | Empirical target for the δ→u IRF shape | 🟡 partial, waits on E9 |
| D11 | Does the reposting margin enter this paper? | 🟡 adopted in the draft; formal go/no-go after R2 |

---

## 🔴 D2. PATH A or PATH B?

`calibrate_shares()` has two mutually exclusive routes through Stages 2–3:

| | Target | ψ | `Xc_Y` |
|---|---|---|---|
| **PATH A** (`steady_state.jl` default; what §5.2 documents) | `Xc_Y = 10%` | ≈ 0.014, outcome | target |
| **PATH B** (every §5.3 runner) | `dest_elast_target` (5.0 is a placeholder) | ≈ 0.033, outcome | outcome, ≈ 14.5% |

Both ψ values are correct for their path. The question is which object carries the prior.
The mechanism runners also override `b_ratio = 0.9, x_v = 0.5` against the `TARGETS` defaults
`0.71, 1.0`, so "the baseline" currently means different things in different files.

**Recommendation: PATH B, with `dest_elast_target` estimated, and §5.2 rewritten to match.**
- The exit elasticity is externally checkable (Broer et al. IER 2025 prefer ψ near 1) and is
  an interpretable object for a prior.
- `Xc_Y = 10%` rests on one Abraham et al. number, and PATH B delivers a plausible 14.5%.
- PATH B inverts analytically where PATH A root-finds, which matters inside the sampler.
- **Entangled with D1:** `dest_ann = 0.0320` does not converge under PATH A (the Stage-4 Q
  root-find fails) but does under PATH B (δ_e/τ = 0.0873, u = 0.0703). Implement D1 by
  switching the default to PATH B.

Consequence: Θ_e = {b, x_v, ξ⁻¹, ω_δ, p_0, σ, dest_elast_target, α, (ρ_x, σ_x)×3}. Under PATH B,
p_0 is estimated because it converts the pinned ψ into `cons`, not because "only the product is
identified" (the current Stage 3 text). Write the canonical target set into
[`data_and_files.md`](data_and_files.md) once decided.

## 🔴 D3. How are β̂ and β(θ) made comparable?

β̂ is a cross-state response per 1 SD of the Bartik instrument (17.4 pp) with time fixed
effects absorbing the national component. β(θ) is an aggregate response per 1 SD of the
structural shock. The draft never says how they are reconciled.

The literature's answer (Guren, McKay, Nakamura & Steinsson, NBER Macro Annual 2021; the
"missing intercept" of Wolf 2019/2023) is to **simulate the same regression inside the
model**: cross-regional estimates discipline structural parameters but do not translate into
aggregate effects without a GE model.

Options, most to least rigorous:
1. **Model-implied relative IRF** (GMNS-faithful). This *is* task E9 if E9 is scoped to
   produce a model-implied cross-regional coefficient.
2. **Peak-normalize both sides**: Block B identifies shape and persistence, σ_δ comes from
   Block M. **Current baseline recommendation.**
3. **Free scale λ_β with a tight prior**: appendix robustness.
4. Compare β̂ to the aggregate IRF directly: **do not** (the missing-intercept error).

**Recommendation:** option 2 as baseline, option 3 in the appendix, cite GMNS in §5.2. Upgrade
to option 1 if E9 delivers it. Scope E9 with both purposes in mind.

Sources: [GMNS, NBER WP 26881](https://www.nber.org/system/files/working_papers/w26881/w26881.pdf) ·
[Wolf, "The Missing Intercept"](https://economics.mit.edu/sites/default/files/publications/missing_intercept.pdf)

## 🟠 D4. Does the LD→v IRF enter Block B?

Status quo (S6): Block B is δ-only, LD→v sits in the appendix, τ is calibrated externally.
The premise was that LD is contaminated because it mixes reposted separations with deliberate
position cuts. **Under D11 that mixture is the object the reposting margin parameterizes**, and
the persistently negative LD→v IRF (trough −0.45 pp) is the missing middle in the data. If D11
is confirmed, LD→v becomes the natural identifying moment for α and should enter Block B.
Revisit with R2.

## 🟡 D10. What is the empirical target for the δ→u IRF shape?

**E8 (Sept 22–23).** Lagged-instrument LPs (p = 4, 8, 12, 16):
- Within-state autocorrelation of B̃^δ is 0.91 at lag 1.
- The δ→u peak moves from h = 17 (baseline) to h = 10–13, with magnitude 0.95–1.86 pp.
- A non-monotone shape appears (significant h = 0–2, insignificant h = 3–8, second rise
  h = 9–13) with no theoretical precedent.
- The δ→v trough at h = 4–6 (−0.55 to −0.65 pp) is robust.

The BED four-quarter recognition rule is a confirmation lag, not a dating shift, so it cannot
explain the late peak (BLS Handbook of Methods; cite the Handbook, not Sadeghi 2008).

**Options:** (a) target only robust short-horizon features plus scale; (b) baseline-to-baseline
matching via a model-side LP (E9), which preserves the full IRF design; (c) hybrid with
horizon reweighting; (d) defend the augmented shape (high burden).

**Recommendation:** test (b) first; fall back to (a). **E9 design caveat:** the simulated
instrument must reproduce the data instrument's persistence (0.91 quarterly), not the shock's
(ρ_δ = 0.592 monthly ≈ 0.21 quarterly). That needs industry-level δ processes with their own
serial correlation and heterogeneous state exposure.

## 🟡 D11. Does the reposting margin enter this paper?

**Status (Oct 1).** Written into the draft on this branch (recruiter block, `eq:v_lom`,
`eq:agg_repost_costs`, `eq:rc`/`eq:gdp`, `def:equilibrium`, wage appendix). Props 1 and 3
audited and hold; Prop 5 restated for the costless-reposting benchmark. Not yet in the code
(R8) or in §5.2 (α has no prior). The formal go/no-go waits on the identification check R2.
Equations: [`model_equations.md`](model_equations.md) Block 2. Argument:
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).

**Why the margin is needed.**
- Model cor(u,v) is **+0.995** against **−0.804** in the data.
- Raising δ̄_e helps monotonically but only to +0.887 across the defensible range; almost all of
  δ_e's leverage lies where measurement rules it out.
- Silencing the s shock gives −0.911, so the Beveridge mechanics are sound. The s shock under
  costless reposting is the problem.
- **Inescapable:** ∂v_pre/∂s = (1−δ_e)(1−u) > 0 is an accounting identity (note §4.3a; Prop. 5
  Part 1). s cannot be silenced: CK's contribution is that separation shocks drive Beveridge
  dynamics, and the s process is data-calibrated (σ_s = 0.0854). Three shocks are required.
- The missing middle: permanent position cuts at surviving firms. GS and Shao-Silos fold it
  into δ (overstating variety destruction); costless reposting folds it into s.

**Formulation (settled Sept 28–29).**
- Each vacated position draws a reactivation cost α·χ, χ ~ F (the continuation-cost
  distribution), reactivated iff Q ≥ α·χ: Λ_r = F(Q/α) ∈ [1−p_0, 1]. One new parameter α.
- A cost, not a constant: a fixed reposting rate is not invariant to the shocks studied, and a
  homogeneous fee gives Λ_r ∈ {0,1}. Dispersion delivers an interior, procyclical rate.
- Threshold is the stock value Q (one-time keep-or-retire), not the flow value K.
- Per-position draws; the recruiter/retailer segmentation keeps the retailer block unchanged.
- **Estimate α, do not calibrate it.** No clean data object exists, and Λ_r adjudicates a
  headline moment.

**Costs.** With Λ_r < 1 both shocks can lower vacancies, so the δ–s asymmetry becomes
quantitative rather than a sign result. The variety asymmetry survives: only δ moves N.

**Channel attribution.** Baseline-vs-AGS moves three channels at once (variety, exit,
reactivation, since AGS sets p_0 = 0). One-at-a-time switches: ζ→0 (variety), ω_δ→0 (exit,
reactivation retained), α→0 (reactivation, exit retained).

**Open before go/no-go:** R2. α and δ_e act through the same channel (the share of separations
that destroy a vacancy), so cor(u,v) cannot separate them; the LD→v IRF should.

---

## Settled — do not reopen without a reason

| # | Decision | Settled | Why |
|---|---|---|---|
| D1 | **δ̄_e = BED establishment deaths, employment-weighted, 1993–2019: `dest_ann = 0.0320`, δ_e/τ = 0.087**, fixed (not estimated) for now | Sept 6, 2026 | The old 0.0754 mixed denominators (JF's 21% is a share of gross job losses, applied to τ). Employment weighting, because δ_e/τ is a share of worker flows and the δ instrument is employment-weighted BED Deaths. Establishments, not firms, for consistency with the instrument. The literature brackets δ in [3.2%, 10%]/yr (note), so describe 3.2% as a lower bound and footnote that symmetry forces one δ_e to serve both product-line and employment margins. Fixed because, without the reposting margin, an estimated δ_e would absorb the s-process misspecification; move it into Θ_e only after D11 (R6). **Not yet implemented** (D1 cascade) |
| D5 | Comparison B's χ^c channel strength is folded into D2 | Sept 5, 2026 | Estimating `dest_elast_target` makes the channel size a reported result; drop the separate ψ ≈ 1 recalibration |
| D6 | Separate LP primary, joint LP as robustness for δ | June 2026 | δ→v trough stable across both (−0.582 vs −0.535); joint zero β^LD reflects collinearity |
| D7 | Industry credit control (Route B) dropped, noted as a limitation | Recommended Sept 2026; not formally confirmed by MS | Data access is binding; the severity placebo carries the external-validity argument |
| S1 | HP λ = 1,600, not 100,000 | May 19, 2026 | λ = 100k flips the sign of cor(δ,u) and cor(δ,v). `app:filter_robustness` |
| S2 | BED **Deaths** (dataclass 08), not Closings | May 14, 2026 | Permanent exits are the model object; cut r(δ,LD) from 0.423 to 0.334 |
| S3 | Vacancy outcome = vacancy **rate** in pp, levels difference | — | Symmetry with the unemployment rate. Never log |
| S4 | LD is the primary s-type instrument; QU is the (failing) placebo; TS is appendix only | — | `principles.md` §3, §11 |
| S5 | LOO national shock rates; 2006 QCEW base-year shares; COVID cap 2019Q4 | — | `principles.md` §5, §6, §9 |
| S6 | τ calibrated externally at 3.1%/month, not identified from the LD IRF | — | ⚠️ Under review via D4 |
| S7 | Two-block posterior with degrees-of-freedom normalization | June 2026 | Each quadratic form ≈ χ²(K); stops the 42-dim IRF block swamping Block M. `eq:posterior` |
| S8 | Filter asymmetry resolved on the model side | June 2026 | m(θ) simulated then filtered; β(θ) read unfiltered off the state space (Canova-Ferroni 2011) |
| S9 | ε = 4.3 fixed externally | — | BK needs only ε > 1; Prop. 5 Part 2 needs ε > 2. [`estimation_design.md`](estimation_design.md) §Determinacy |
| S10 | **Timing: the draft is the reference model.** Convention A (matches formed at t face exit but not separation before producing); exit at Stage 2 against the current cutoff, so survival t→t+1 is (1−δ_t)F(χ^c_{t+1}) in every value function and law of motion; the aggregate state (incl. δ_t, s_t) is known at the start of t, and only the incidence of δ_t, s_t is revealed at the end of t | Oct 1, 2026 | Convention A is standard DMP, matches Gabrovski-Silva, and conserves positions; the code adopted it Sept 6. The draft's Stage-2 timing is what the propositions use (Prop. 5 Part 2's joint (χ^c, N) system). The code's two deviations (exit dating, post-exit states) are fixed in R8; spec in [`model_equations.md`](model_equations.md) "Code status" |
| S11 | **Reactivation costs X^r are paid at Stage 2 of t on positions vacated at the end of t−1, valued at the realized Q_t** | Oct 1, 2026 | Forced by S10: reactivation happens at the Stage-2 node. Paying at separation against E_t Q_{t+1} would make Λ_r condition on an expectation and lose its procyclicality |
| S12 | **Prop. 5 (restated Oct 5, wording approved by MS).** Part 1 (p_0 = 0): u and v_pre rise by the reposting inflow d up to the period-t entry response; positive comovement at h = 1 under conditions (i)–(ii). Part 2 (p_0 > 0, α → 0): to first order in F′, under `eq:exit_cond` (det J > 0) and inflow-dominates-drain. Part 3: h = 1, with the period-t entry term | Oct 5, 2026 | Under the timing convention (S10) θ_t and e_t respond to period-t shocks. At p_0 = 0, (i) holds automatically and (ii) binds at most 2.6% (0.7% at D1) (`run_prop5_s_check.jl`). Part 3 keeps h = 1 although the first-month sign is fragile at monthly frequency (entry offsets 81% at D1); vacancies fall from month 2 and in every quarterly average at D1 (`run_prop5_parts23_check.jl`). `eq:gN_cond` was dropped because it fails at the calibration |
