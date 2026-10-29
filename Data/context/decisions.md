# Decision Register
**Last updated:** October 14, 2026 (D13 adopted: SVAR Block B′ replaces Bartik; D3 moot, D10 → (a) as validation). October 12: D12 updated for the Oct 11 theory rebuild; S12 rewritten for the new `sec:mechanism` structure; D1 entry compressed). Earlier: October 9 (D12 added) · **Branch:** `costly_vacancy_reposting`

One entry = one decision: a *choice*, not a task. Tasks live in
[`pending_tasks.md`](pending_tasks.md). When a decision is settled, move it to §Settled with
the date and the reason; the draft's prose depends on knowing *why*.

Status key: 🔴 blocks estimation · 🟠 blocks a paper section · 🟡 improves the paper

| # | Decision | Status |
|---|---|---|
| D2 | PATH A or PATH B as the calibration baseline | ✅ **settled Oct 18, 2026: PATH B** (code default since Q1 step 5) |
| D3 | How β̂ and β(θ) are made comparable | ⚪ moot under D13 (Oct 14): indirect inference uses the same SVAR on data and simulations |
| D4 | Does the LD→v IRF enter Block B? | 🟠 open, reopened by D11 |
| D10 | Empirical target for the δ→u IRF shape | 🟡 resolved to (a) as validation under D13 (Oct 14); E9 no longer needed |
| D11 | Does the reposting margin enter this paper? | 🟡 in the draft and the code (Q1, Oct 19); go/no-go after Q2 + R2/Q4 |
| D12 | Venue and the theory-first pivot (JET?) | 🟠 open (MS), recommendation: no-go for JET; theory-forward paper for RED / JME / AEJ: Macro |
| D13 | Conditional moments without Bartik: Block B′ = aggregate SVAR responses | 🔴 adopted as the working design Oct 14 (MS approved the plan); not yet built |
| D14 | Is the s innovation orthogonalized against z and δ? | 🟡 low stakes after Q2 (Oct 20): z and δ explain 1.1 % of the s innovation, so (a) and (b) coincide; take (a). New question: the dropped cross-lag δ_{t−1} → s_t (t = 2.8) |

---

## ✅ D14. Is the s innovation orthogonalized against z and δ? — **settled Oct 28 (MS): no; independent AR(1)s, σ_s not purged** (see the Oct 28 update below)

**Raised Oct 19, 2026 (Q2 scoping).** `part6b_shock_persistence_cyclical.py` orthogonalizes
the δ innovation against z (Cholesky; σ_δ falls) but fits s as a **univariate AR(1)**, so
σ_s = 0.0854 is a reduced-form residual SD that still contains whatever coincides with the z
and δ innovations. The model treats all three innovations as uncorrelated.

**Why the z → δ argument does not carry over.** For δ, stripping the z-correlated part avoids
double counting: the model generates z → χ^c → δ_e itself. The model has **no** endogenous s.
Stripping the z- and δ-correlated part of s therefore discards real countercyclical separation
variation that the model cannot reproduce. The script's line "the same argument applies to s"
is wrong on this point.

**Options.**
- (a) **Orthogonalize, z → δ → s.** Smaller σ_s, so a smaller wrong-sign Beveridge contribution
  from s, and consistency with the B′ identification (D13). Cost: the model undershoots σ(τ)
  and the cyclicality of separations, which must be reported.
- (b) **Correlated innovations** (as CK for (p, δ)). Faithful to the data, but it puts cross-shock
  correlation back into a paper whose argument separates δ from s.
- (c) **Endogenous separations** responding to z. A model change; out of scope.

**Also relevant.** Noise in BED Deaths enters s = (τ − δ)/(1 − δ) with weight ≈ δ/τ ≈ 0.09, which
builds a spurious *negative* s–δ innovation correlation. Only the joint system shows its size.

**Oct 20 update (Q2, `part6c_shock_var_trivariate.txt`).**
- The premise is rejected. The s innovation correlates −0.105 with z and +0.041 with δ, and
  orthogonalizing removes 1.1 % of its variance. (a) and (b) differ in the third decimal of
  σ_s, and no BED-noise correlation appears.
- Take (a), for consistency with B′.
- The live issue is a different one: lagged δ predicts s (+0.362, t = 2.81). The model lacks a
  δ → s channel, and the diagonal convention drops it. Whether to model it (an exogenous
  δ_{t−1} term in the s process, or an endogenous separation margin) is open, and belongs with
  the Beveridge-gap question in D11.

**Oct 28 update (MS; `shock_calibration.txt` §5).**
- **Choice:** the model takes three independent AR(1)s, and s is **not** purged. The model has
  no endogenous s, so purging would discard variation it cannot regenerate. This replaces the
  Oct 20 "take (a)". ρ_s now comes from the univariate AR(1), not the VAR diagonal.
- **Rule fixed in advance:** revisit if the z innovation explains ≥ 10 % of the s innovation
  and |t| ≥ 2.
- **Result on the S14 window (1951–2019), after both τ/jf fixes: the rule passes, narrowly.**
  R² = 9.2 % (threshold 10 %), corr −0.30, t = −4.38. 1951–84 R² 13.3 % (loading −2.13);
  1985–2019 R² 2.4 % (t = −1.44). Before the jf/τ timing fix it read 19.1 %. By the rule the
  baseline is option (i), independent s with σ_s not purged. The correlation is significant, so
  report it, and (ii) remains a natural robustness check.
- **Robustness options** (the baseline is (i)):
  - (i) keep independent s and report the missing z–s comovement;
  - (ii) load s on ε^z, i.e. s_t = ρ_s s_{t−1} + b·ε^z_t + ε^s_t;
  - (iii) estimate the s process on 1985–2019 only, where independence holds;
  - (iv) endogenous separations, which is a model change.

  Caveat for (ii): it raises z's model share of var(u) through an exogenous channel. It must
  be justified by the data, not by the fit (z is under-amplified 3–6×).

**Original recommendation (Oct 19):** (a) for the baseline, (b) as robustness, with the σ(τ)
shortfall reported.
Decide after the Q2 numbers (share of s innovation variance explained by z and δ).

---

## 🔴 D13. What replaces the Bartik IRF block?

**Adopted Oct 14, 2026 as the working design (MS approved the plan).** The Bartik LP paths are
unstable beyond short horizons (E8). They therefore leave the objective and become
**out-of-sample validation**: δ→u at h = 0–2, the δ→v trough at h = 4–6, the severity
placebo, and the LD→v IRF.

**Block B′: aggregate SVAR conditional responses.**
- A VAR in (z, δ, s, u, v), with the shocks ordered first, as in the shock calibration.
- Targets: the δ-shock and s-shock responses of u and v, h = 0–12 quarters.
- Model counterpart: the same VAR estimated on simulated data (indirect inference), so the
  same estimator runs on both sides. This settles D3 by construction.
- Ω from a VAR bootstrap.

Why Block B′ is still needed: unconditional moments cannot separate δ from s propagation. The
s-shock response of v identifies α; the δ-shock response identifies ξ jointly with δ_e.

**Consequences:**
- E9, Ω_β (`part5_wcrb.py`) and the β(θ) extractor in LP units drop off the critical path.
- D10 resolves to (a) as validation.
- D4's LD→v becomes a validation moment.
- The shock processes come from a trivariate (z, δ, s) VAR (Q2), not a bivariate VAR plus an
  AR(1) for s.

The estimator stays limited-information (Bayesian SMM / indirect inference) over Θ_e, with
the shock processes fixed in a first step. Full-information estimation is deferred: it needs
measurement error to avoid stochastic singularity.

The full sequence is tasks Q1–Q6 in [`pending_tasks.md`](pending_tasks.md). D12's
simple-model-first question is decided at the Q5 gate, against a threshold agreed before Q5
runs.

---

## 🟠 D12. Venue, and whether to pivot to theory-first for JET

**Raised Oct 9, 2026 (MS).** The proposal: lead with theory (`prop:bgm_nest`, `prop:equilibria`,
`prop:fe`, `prop:exit_s`), keep the SVAR-calibrated "Inspecting the mechanism" section,
calibrate or estimate the rest, and move the Bartik/LP work to a separate empirical paper.
Full assessment: [`session_handout_20261009.md`](session_handout_20261009.md) §2.

**Recommendation: no-go for JET; adopt the theory-forward structure for a field macro journal.**
- **The theory reads as applied theory.**
  - Four propositions are nesting or isomorphism results; `prop:curves` is now a remark.
  - `prop:fe` is clean but log-linear, at the ξ → ∞ corner, with variety, exit and reposting off.
  - The finite-ξ result is a numerical remark (`rem:finite_xi`), not a theorem.
- **Oct 12 update (after the Oct 11 rebuild).** The repairs strengthen the field-journal case
  without changing the category. `prop:equilibria` now holds under stated matching regularity
  (`ass:matching` + (C)), with the SMCD class as a nesting example; `prop:double_limit` is
  corrected and motivates clone replacement; `sec:mechanism` leads with `lem:theta` +
  `prop:fe`, followed by `lem:accounting`, `rem:finite_xi`, `prop:exit_s`. A new implication:
  under costless reposting, s shocks give positive comovement at every persistence, so the
  reposting margin (D11) is necessary, not optional, for the Beveridge sign.
- **Every dynamic theorem switches the business formation block off**, so the CK departure in
  the theory rests on GS's δ/s split plus the BED measurement.
- **Raising the odds would take months of high-risk theory.** That means finite-ξ results or
  variety inside the dynamic theorem. R17 would help but does not change the category.
- **A JET version contradicts the stated objective** (JME / AEJ: Macro / JPE, `../CLAUDE.md`).
  After a rejection, rebuilding the quantitative part would cost months.
- **The empirical paper is not viable alone.** The QU placebo fails, LD is not separately
  identified, and the δ→u peak is unstable (E8).

**Recommended structure:**
- lead with `prop:fe` (R15 option (b); **in the draft since Oct 11**);
- keep the SVAR-calibrated mechanism (the SVAR is bivariate z, δ; s from a univariate AR(1));
- cut Bartik to the robust short-horizon facts (δ→u at h = 0–2, the δ→v trough at h = 4–6, the
  severity placebo) as motivating evidence;
- defer full-information estimation.

Venue ranking: RED, JME, AEJ: Macro, then JEDC as the safe option. TE is not a fallback.

**Consequences if adopted:**
- D10 → option (a).
- D3 and D4 become largely moot, because IRF matching leaves.
- Ω_β and E9 drop off the critical path.
- R8, D1, PATH B and D11 stay on it, because the quantitative Beveridge sign must be fixed for
  any venue.

---

## ✅ D2. PATH A or PATH B? — settled Oct 18, 2026: PATH B

**Implemented in Q1 step 5** (branch `r8_timing_reposting`). `TARGETS` in `steady_state.jl` now
holds `dest_ann = 0.0320` (D1), `dest_elast_target = 5.0` (placeholder, estimated later),
`b_ratio = 0.9`, `x_v = 0.5` (the mechanism-runner values, so "the baseline" means one thing in
every file) and `Λr_target = nothing` (costless reposting). PATH A remains reachable by dropping
`dest_elast_target`. At the defaults (`R8_step5_defaults.txt`): δ_e/τ = 0.0873, ψ = 0.0136,
X^c/Y = 13.5 % (outcome), SS residual < 1e-12. PATH B switches itself off at p_0 = 0. **Owed:**
§5.2 still documents PATH A with b/w = 0.71, x_v = 1.0 (P4).

The record below is the pre-decision analysis (ψ values there were at `dest_ann = 0.0754`).

`calibrate_shares()` has two mutually exclusive routes through Stages 2–3:

| | Target | ψ | `Xc_Y` |
|---|---|---|---|
| **PATH A** (old default; what §5.2 documents) | `Xc_Y = 10%` | ≈ 0.014, outcome | target |
| **PATH B** (every §5.3 runner; default since Q1 step 5) | `dest_elast_target` (5.0 is a placeholder) | ≈ 0.033, outcome | outcome, ≈ 14.5% |

Both ψ values are correct for their path. The question is which object carries the prior.
Before Q1 step 5 the mechanism runners also overrode `b_ratio = 0.9, x_v = 0.5` against the
`TARGETS` defaults `0.71, 1.0`.

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

**Oct 9, 2026 — linked to D12.** If D12's recommended structure is adopted (no IRF matching in
this paper), take option (a). The robust short-horizon facts become motivating evidence, not
estimation targets, and E9 is no longer needed.

## 🟡 D11. Does the reposting margin enter this paper?

**Oct 19 record (Q1 closed).** Implementation is settled; the go/no-go is not.
- **Settled (implementation):** reposting enters through a target Λ̄_r (`Λr_target`) with α
  backed out; shared F with no separate reactivation shape (MS, Oct 14); threshold Q;
  per-position draws; the switch `rp` gives the costless benchmark (`Λr_target = nothing`, the
  default). Code: `run_solution_core.jl` f[24], f[39]–f[41], 𝓡 = Q − K − H in f[3] and f[7],
  X^r in f[16]. Tests 3–5 verify the steady state (14 identities) and position conservation
  with retirement (`R8_test_steady_state.txt`, `R8_test_reposting.txt`).
- **Θ_e:** Λ̄_r replaces α ([`parameters.md`](parameters.md)).
- **Open (go/no-go):** whether a Λ̄_r inside (1 − p_0, 1) delivers cor(u,v) ≈ −0.8 once the
  shock processes come from the Q2 VAR, and whether Λ̄_r is identified separately from σ_s and
  p_0 (R2/Q4). At the part6b shocks the grid reaches only +0.461 at Λ̄_r = 0.70.
- **Related:** R19 (pricing replacement offsets `prop:fe`'s user-cost channel) is a separate
  theory decision and does not reopen the formulation.

**Oct 18 status.** In the code on `r8_timing_reposting` (Q1 steps 3–4, Oct 14): controls Λ_r,
Q_rep, H; Λ̄_r is targeted and α backed out. Shared F (MS, Oct 14): reactivation is nearly
inelastic (d ln Λ_r/d ln Q ≈ 0.006), so the margin acts through the level Λ̄_r, not its
cyclicality. First simulated cor(u,v) at the defaults (`R8_step7_baseline.txt`, part6b shocks):
+0.995 costless, +0.968 (0.95), +0.812 (0.85), +0.461 (0.70), against −0.804. Direction right,
size not yet; the go/no-go waits on Q2 (shock processes) and R2/Q4 (identification).

**Status (Oct 1).** Written into the draft on this branch (recruiter block, `eq:v_lom`,
`eq:agg_repost_costs`, `eq:rc`/`eq:gdp`, `def:equilibrium`, wage appendix). Props 1 and 3
audited and hold; `lem:accounting` and `prop:exit_s` cover the costless-reposting benchmark. Not yet in the code
(R8) or in §5.2 (α has no prior). The formal go/no-go waits on the identification check R2.
Equations: [`model_equations.md`](model_equations.md) Block 2. Argument:
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).

**Why the margin is needed.**
- Model cor(u,v) is **+0.995** against **−0.804** in the data.
- Raising δ̄_e helps monotonically but only to +0.887 across the defensible range; almost all of
  δ_e's leverage lies where measurement rules it out.
- Silencing the s shock gives −0.911, so the Beveridge mechanics are sound. The s shock under
  costless reposting is the problem.
- **Inescapable:** ∂v_pre/∂s = (1−δ_e)(1−u) > 0 is an accounting identity (note §4.3a;
  `lem:accounting` Part 1). At free entry `prop:fe` part 3 makes it a theorem: the s response is
  positively sloped at every persistence (D1: (1−δ̄)f̄[η_L+ϕ(1−δ̄)f̄/τ̄] = 4.05 ≥ 1). s cannot be silenced: CK's contribution is that separation shocks drive Beveridge
  dynamics, and the s process is data-calibrated (σ_s = 0.0854). Three shocks are required.
- The missing middle: permanent position cuts at surviving firms. GS and Shao-Silos fold it
  into δ (overstating variety destruction); costless reposting folds it into s.

**Formulation (settled Sept 28–29).**
- Each vacated position draws a reactivation cost α·χ, χ ~ F (the continuation-cost
  distribution), reactivated iff Q ≥ α·χ: Λ_r = F(Q/α) ∈ [1−p_0, 1]. One new parameter α.
- A cost, not a constant: a fixed reposting rate is not invariant to the shocks studied, and a
  homogeneous fee gives Λ_r ∈ {0,1}. Dispersion delivers an interior rate that moves with Q
  (⚠ Oct 14: not procyclical in general. Q falls after z and s shocks, rises after δ shocks at
  finite ξ, and is fixed under vacancy free entry; see `model_equations.md` Block 2).
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
| D1 | **δ̄_e = BED establishment deaths, employment-weighted, 1993–2019: `dest_ann = 0.0320`, δ_e/τ = 0.087**, fixed (not estimated) for now | Sept 6, 2026 | The old 0.0754 mixed denominators (JF's 21% is a share of gross job losses, applied to τ). Employment weighting, because δ_e/τ is a share of worker flows and the δ instrument is employment-weighted BED Deaths. Establishments, not firms, for consistency with the instrument. The literature brackets δ in [3.2%, 10%]/yr (note), so describe 3.2% as a lower bound and footnote that symmetry forces one δ_e to serve both product-line and employment margins. Fixed because, without the reposting margin, an estimated δ_e would absorb the s-process misspecification; move it into Θ_e only after D11 (R6). **Not yet implemented** (D1 cascade, critical path 0b). A lower δ_e helps the δ-shock Beveridge sign through the stock elasticity ξ·ē/v̄ (finite ξ, `rem:finite_xi`) and the user cost of a durable vacancy (free entry, `prop:fe`: slope −0.24 at D1 vs +0.17 at δ_e = τ, ρ_δ = 0.592). (The Oct 5 threshold-T argument was withdrawn Oct 6; see [`../Notes/predetermined_vacancies_and_the_CK_special_case.md`](../Notes/predetermined_vacancies_and_the_CK_special_case.md).) |
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
| S9 | ε = 4.3 fixed externally | — | BK holds at ε = 4.3 with a wide margin: a bounded equilibrium exists for ε > 2.40 at the defaults (`bk_epsilon_threshold.jl`, Oct 18; the earlier "ε > 1" from `Notes/Baseline_Blanchard_Kahn.md` holds for the LR-BGM block but not the full model, which fails below 1.40 even at p_0 = 0). No proposition needs ε > 2 any more: `eq:gN_cond` was dropped Oct 5, and `prop:exit_s`'s `eq:exit_cond` needs no sign on g_N. [`estimation_design.md`](estimation_design.md) §Determinacy |
| S10 | **Timing: the draft is the reference model.** Convention A (matches formed at t face exit but not separation before producing); exit at Stage 2 against the current cutoff, so survival t→t+1 is (1−δ_t)F(χ^c_{t+1}) in every value function and law of motion; the aggregate state (incl. δ_t, s_t) is known at the start of t, and only the incidence of δ_t, s_t is revealed at the end of t | Oct 1, 2026 | Convention A is standard DMP, matches Gabrovski-Silva, and conserves positions; the code adopted it Sept 6. The draft's Stage-2 timing is what the propositions use (`prop:exit_s`'s joint (χ^c, N) system). The code's two deviations (exit dating, post-exit states) are fixed in R8; spec in [`model_equations.md`](model_equations.md) "Code status" |
| S11 | **Reactivation costs X^r are paid at Stage 2 of t on positions vacated at the end of t−1, valued at the realized Q_t** | Oct 1, 2026 | Forced by S10: reactivation happens at the Stage-2 node. Paying at separation against E_t Q_{t+1} would make Λ_r condition on an expectation and lose its procyclicality |
| S13 | **Sample endpoint 2019Q4 for every series, applied *before* filtering** (shock processes, Block M moments, Block B′ VAR) | Oct 20, 2026 (MS) | Rules out pandemic effects. HP filtering through 2020+ leaks the 2020 spike into the 2017–19 cycles, and the old moment targets (σ(u) 14.96, cor(u,v) −0.804) included 2020–24; the pre-pandemic values are 12.48 and −0.925 (1951–2019). The **start** of the window is still open (findings.md "Sample windows") |
| S14 | **Sample window, option B: 1951Q1–2019Q4 for z, s and the labor-market moments (u, v, θ, f, s, z); 1992Q3–2019Q4 only for δ and δ-related moments** | Oct 20, 2026 (MS) | The literature standard (Shimer 2005, Hagedorn–Manovskii 2008, Coles–Kelishomi 2018) and the comparison referees will make. 1992–2019 alone sits entirely in the post-1985 era, when labor productivity is acyclical with u (cor +0.13 vs −0.37 before 1985), so z looks weak by construction. **Spec:** (i) s before 1992Q3 built from τ with δ held at a constant δ̄, (τ − δ̄)/(1 − δ̄); **amended Oct 28 (MS): δ̄ = mean of the first 10 BED quarters (1992Q3–1994Q4), which removes the step at the splice**; the measured identity from 1992Q3 (cyclically s ≈ τ: correlation 0.99, δ/τ ≈ 0.09). (ii) Shock processes, **amended Oct 28 (MS)**: z an AR(1) on 1951–2019; δ from the (z, δ) VAR on 1992–2019, keeping ρ_δ and σ_δ ⊥ z (the part6b double-counting argument); s a univariate AR(1) on 1951–2019, σ_s not purged (D14). Off-diagonals excluded. All in `part6b_shock_persistence_cyclical.py`. (iii) Block M and B′ moments on the same windows. (iv) Robustness: 1985–2019. (v) End cut at 2019Q4 before filtering (S13). **Implemented Oct 28** (`shock_calibration.txt`) |
| S15 | **Data SVAR identification (Block B′): z first; preferred ordering z, s, δ, then u, v; z, δ, s reported for completeness; u and v residual innovations combined (their relative order is arbitrary)** | Oct 29, 2026 (MS) | z is the most plausibly exogenous series, and the model generates z → δ_e within the period, so z-driven exit counts toward z. s is purely exogenous in the model, while measured exit responds to the cycle within the period, so δ goes after s. The δ–s order moves only the δ–s split (`q3_data_svar.py`) |
| S12 | **Structure of `sec:mechanism` (R15 option (b), MS, Oct 11).** (1) `lem:theta` + `prop:fe`: free entry (ξ = ∞), closed-form δ–s contrast, rigorous statement; informal discussion with `tab:fe_irf`. (2) `lem:accounting`: old Prop. 5 Parts 1 (s) and 3 (δ), h = 1 accounting at any entry technology, p_0 = 0. (3) `rem:finite_xi`: the finite-ξ entry cushion, numerical. (4) `prop:exit_s`: old Part 2 (p_0 > 0, α → 0), to first order in F′, under `eq:exit_cond` and inflow-dominates-drain | Oct 5 (statements), Oct 11 (structure) | `prop:fe` is sharper where both apply and carries the δ–s asymmetry as a theorem; the accounting lemma covers finite ξ; exit is a separate mechanism. Under S10 θ_t and e_t respond to period-t shocks. At p_0 = 0 condition (ii) binds at most 2.6% (0.7% at D1) (`run_prop5_s_check.jl`); the δ part's first-month sign is fragile at monthly frequency (entry offsets 81% at D1) but vacancies fall from month 2 (`run_prop5_parts23_check.jl`). The old Part 3 sufficient condition (R12c) was wrong and became `rem:finite_xi` |
