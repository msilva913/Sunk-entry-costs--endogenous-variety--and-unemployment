# Decision Register
**Last updated:** September 22, 2026 (D10 updated with E8 results) · **Branch:** `Organize_Project_State_Estimation`

Every open decision that must be settled before estimation, plus the ones already settled
that agents keep re-litigating. One entry = one decision. When you settle one, move it to
§Settled with the date and the reason — do not delete it, because the draft's prose depends
on knowing *why*.

Status key: 🔴 blocks estimation · 🟠 blocks a paper section · 🟡 improves the paper · ⚪ optional

---

## ✅ D1. Which δ̄_e? — SETTLED September 6, 2026

**Decision: BED establishment deaths, employment-weighted, 1993–2019.**
`dest_ann = 0.0320` (3.2 %/yr, 0.271 %/month), **δ_e/τ = 0.087.**
Window chosen by MS for the longer sample; 2001–2019 would give 0.077.

### Why this number

1. **The status quo was an arithmetic error.** `dest_ann = 0.0754` came from
   Jaimovich-Floetotto's 21% × τ. Their 21% is a share of *gross job losses*; τ is the
   *total separation rate*. Mixing the denominators inflates δ̄_e by roughly 1.8×. Applied to
   its own base, JF's share gives ~1.3 %/qtr, close to the BED measurement. Settled
   independently of any question about model fit.

2. **Employment weighting, not counts.** δ_e/τ is by construction a share of separations, so
   both must be worker flows. Decisively, the δ Bartik instrument *is* employment-weighted
   BED Deaths, so calibrating to anything else would have the empirical section and the
   calibration describing different objects.

3. **Establishments, not firms.** BDS firm deaths (2.26 %/yr employment-weighted,
   δ_e/τ = 0.061) are *lower* than BED establishment deaths, since a firm death requires
   every establishment to exit. Consistency with the instrument decides it. Note this
   reverses an earlier argument in this file that BED was a lower bound because product
   lines die inside surviving firms — the firm-level number is lower, not higher.

### The weighting tension, and how the paper handles it

Under DS-CES symmetry one δ_e serves as both the fraction of *product lines* destroyed
(driving ρ(N)) and the fraction of *employment* destroyed by exit (driving τ). In the data
these differ by 3.7×: firm deaths are 2.26% of employment but 8.33% of firms, because dying
firms are small. The model has no size heterogeneity and cannot match both.

**MS decision: footnote the choice.** State that symmetry forces one δ_e to serve both
margins, that we match the employment-weighted flow, and that this understates variety
destruction. Do not bury it.

### Measurement menu, for the record

| measure | %/yr | δ_e/τ |
|---|---|---|
| BDS firm deaths, employment-weighted | 2.26 | 0.061 |
| BED estab deaths, emp-weighted, 2001–19 | 2.84 | 0.077 |
| **BED estab deaths, emp-weighted, 1993–19 — CHOSEN** | **3.21** | **0.087** |
| BDS job destruction from estab deaths | 4.50 | 0.124 |
| BDS estabs at dying firms / all estabs | 6.57 | 0.182 |
| Code's old JF-derived value | 7.6 | 0.210 |
| BDS firm deaths, count-weighted | 8.33 | 0.233 |
| BGM 2012 | 9.6 | 0.269 |
| BDS establishment exits, count-weighted | 9.93 | 0.280 |
| Coles-Kelishomi (δ_e = τ) | 32.3 | 1.000 |

### Why fixed rather than estimated

The register previously recommended moving δ_e into Θ_e. Rejected, on the corrected
diagnostic (`findings.md` §D1/M5 re-run):

- δ_e moves amplification by a factor of ~2.9 across the candidate range
  (σ(θ)/σ(labor_prod) 1.02 → 2.93) but cannot reach the 11.70 target even at the BGM value,
  falling 4× short.
- The δ→u peak sits at h = 1 quarter for *every* δ_e. The hump-shape gap (**M8**) is
  invariant to it.
- The binding Block M problem is elsewhere: with ρ_s, σ_s at placeholder values the
  unconditional cor(u,v) is **+0.76** against −0.80 in the data. Silencing the s shock gives
  −0.911. Every correlation in the moment table is contaminated until the s process is
  calibrated.

Freeing a cleanly measured parameter while known misspecification sits upstream would let
δ_e absorb blame for the s process and for M8. Fix it at the measurement, resolve those,
then revisit whether it needs to be free.

### Consequences to carry out

- `steady_state.jl` `TARGETS.dest_ann`: 0.0754 → 0.0320, and correct the source comment.
- Re-run Comparisons A–D, regenerate the eight figures, re-run `mechanism_stats.jl`, and
  update §5.3 against its output.
- Draft §5.2 external block and `tab:calib_targets` row 1.
- The "δ̄_e/τ̄ ≈ 0.21" claims at `Draft.tex:1481` and in `app:comparison_D` become ≈ 0.087.
  Prop. 5 Part 3's condition is "δ̄_e/τ̄ small", so the lower value *strengthens* it.
- **P3b**: the intro endorses Gabrovski-Silva's 6–10 %/yr range at `Draft.tex:421`. At 3.2%
  we are at half the floor. GS needed that range in a model without endogenous exit, the
  variety externality, or finitely elastic entry. Whether those compensate is a result the
  estimation delivers, not an assumption — write it as an open question.
- **N15**: the chosen number must be emitted by a program. `observables.py` currently computes
  δ over a COVID-contaminated window (`final = 2025-10-01`). Add a pre-COVID calibration
  target computation there before the number is treated as sourced.

## 🔴 D2. Which calibration path is the estimation baseline — PATH A or PATH B?

**The conflict.** `calibrate_shares()` has two mutually exclusive routes through Stages 2–3:

- **PATH A** (`steady_state.jl` default; the one the draft documents). Target `Xc_Y = 10%`
  → root-find X^c/Y^c → π_s → `cons` → ψ. Here ψ and p_0 are not separately identified from
  steady-state targets — only the product `cons = (ψ/(1+ψ))·p_0` is pinned — which is
  exactly the draft's stated reason for putting **p_0 in Θ_e**.
- **PATH B** (used by every mechanism comparison: `run_solution_endog_exit.jl`,
  `run_solution_entry_elasticity.jl`). Target `dest_elast_target = 5.0` → ψ analytically
  (≈0.033) → `cons` → π_s. Here **`Xc_Y` stops being a target and becomes an outcome**
  (≈14.5%, not 10%), and ψ is pinned rather than estimated.

Draft §5.2 Stages 2–3 and `tab:calib_targets` describe **PATH A**. Every figure in §5.3 was
produced with **PATH B**. The two imply different contents for Θ_e.

**Also note** the mechanism runs override the defaults with `b_ratio=0.9, x_v=0.5`, whereas
`TARGETS` defaults to `b_ratio=0.71, x_v=1.0`. So "the baseline" currently means three
different things depending on which file you open. Whatever is decided here, write the
single canonical target set into [`data_and_files.md`](data_and_files.md).

### Clarified September 5, 2026 — yes, and 5.0 is a placeholder

**Yes: `dest_elast_target` should be estimated, and 5.0 is a placeholder chosen for the
mechanism comparisons, not an estimate.** It appears nowhere in `const TARGETS`; it is
supplied explicitly by `run_solution_endog_exit.jl:164` and `run_solution_entry_elasticity.jl`
so that those comparisons hold the χ^c channel elasticity fixed at a known value while other
things vary.

This means **D2 is not "estimate ψ or don't" — it is a reparameterization question: which
object carries the prior?** In both paths ψ ends up in Θ_d, determined by the targets plus a
draw. What differs is the free object drawn alongside p_0:

| | free object | ψ at the calibrated point | `Xc_Y` |
|---|---|---|---|
| PATH A | `Xc_Y` (fixed at 10%) | ≈ 0.014 (outcome) | target |
| PATH B | `dest_elast_target` | ≈ 0.033 (outcome) | outcome, ≈ 14.5% |

⚠️ This also resolves the apparent ψ conflict in the old notes: **0.014 and 0.033 are both
correct, for PATH A and PATH B respectively.** Neither is an error.

**Recommendation: PATH B, with `dest_elast_target` estimated rather than fixed at 5.0, and
§5.2 rewritten to match.** Three reasons: the exit elasticity is externally checkable
(Broer et al. IER 2025 prefer an implied ψ near 1, so reporting it invites exactly that
comparison, and it is a far more interpretable object to put a prior on than a fixed-cost
share); `Xc_Y = 10%` rests on a single Abraham et al. number and comes out at a still-plausible
14.5% as an outcome anyway; and PATH B replaces a root-find with an analytic inversion in the
inner loop, which matters when the calibration is re-solved at every draw.

**Consequence:** Θ_e becomes {b, x_v, ξ⁻¹, ω_δ, p_0, σ, dest_elast_target, (ρ_x, σ_x)×3} —
13 parameters, one more than the draft's current list. The draft's justification for
estimating p_0 also needs restating: under PATH B, p_0 is what converts the pinned ψ into
`cons`, so it remains estimated but for a different reason than the "only the product is
identified" argument in the current Stage 3 text.

---

## 🔴 D3. How are β̂ and β(θ) made comparable?

**The gap.** `eq:ql_irf` compares β̂ to β(θ), but they are not the same object:

- **β̂** — cross-state response, in pp, per one standard deviation *of the Bartik
  instrument* (σ̂ = 17.4 pp), estimated **with time fixed effects**, so the national
  component of the shock is absorbed.
- **β(θ)** — aggregate response, per one standard deviation of the *structural* δ shock,
  with no cross-sectional dimension at all.

Nothing in the draft states how these are reconciled. Left implicit, the IRF block
identifies only the *shape* of the response, not σ_δ — and possibly not even that, if
time-FE absorption bites differently across horizons.

### What the literature does (reviewed September 5, 2026)

This is a known and well-developed problem, and the literature's answer is **not** "rescale
the cross-sectional estimate." Two strands:

**1. Guren, McKay, Nakamura & Steinsson, "What Do We Learn from Cross-Regional Empirical
Estimates in Macroeconomics?" (NBER Macro Annual 2021)** — the methodological reference.
Their prescription: **simulate the same regression inside the model.** Write down the
multi-region model, run the identical experiment, construct the *model-implied cross-regional
coefficient*, and match that to β̂ — rather than comparing β̂ to the model's aggregate IRF.
The reason is that with time fixed effects the cross-regional estimate differences out
national general-equilibrium effects, so it is a *relative* object; regional responses "can
only be translated into aggregate responses through the lens of a fully specified general
equilibrium model." Cross-regional estimates are excellent for *discriminating between
models and disciplining structural parameters*, and poor for reading off aggregate effects.

**2. The "missing intercept" literature** — Wolf (2019, 2023), and Adão-Arkolakis-Esposito —
names precisely the wedge: aggregate GE effects and the aggregate policy response are
absorbed by the time fixed effects, so scaling up a cross-sectional estimate implicitly
assumes the intercept away. Wolf's demand-equivalence approach recovers the aggregate by
combining cross-sectional with time-series evidence.

**What this means here.** The good news is that the paper only wants Block B to *discipline
structural parameters*, which is exactly what GMNS say cross-regional estimates are good for.
The bad news is that the by-the-book solution — build a multi-region version of the model and
compute the model-implied Bartik LP coefficient — is a serious undertaking and probably out
of scope. So the realistic options, in descending order of rigour:

1. **Model-implied relative IRF (GMNS-faithful, expensive).** Construct a two-region or
   continuum-of-regions version, run the shift-share experiment, form the model analogue of
   β̂ with the aggregate component differenced out. Most defensible; likely a separate paper's
   worth of work.
2. **Peak-normalize both sides** and state that Block B identifies the *shape and persistence*
   of the δ transmission, not its scale, with σ_δ identified from Block M. This is the honest,
   cheap version of the GMNS point: shape is the part of the relative IRF least contaminated
   by the missing intercept. **Recommended baseline.**
3. **Free scale parameter λ_β with a tight prior**, estimated as a nuisance. Retains level
   information; the prior is doing real work and a referee will say so. **Recommended as the
   appendix robustness**, since it also reveals how much the shape-only restriction costs.
4. Compare β̂ directly to the aggregate model IRF. **Do not do this** — it is precisely the
   missing-intercept error, and it is the current implicit specification.

**Recommendation: option 2 as the baseline, option 3 in the appendix, and cite GMNS explicitly
in §5.2** when stating that the cross-regional design identifies relative propagation. Framing
it as a deliberate, literature-grounded choice converts an obvious referee objection into a
paragraph that shows command of the method.

Sources: [GMNS, NBER WP 26881](https://www.nber.org/system/files/working_papers/w26881/w26881.pdf) ·
[Wolf, "The Missing Intercept"](https://economics.mit.edu/sites/default/files/publications/missing_intercept.pdf)

---

## 🟠 D4. Does the s/LD IRF enter the quasi-likelihood block?

Block B is currently δ→u and δ→v only. LD is not separately identified (r(δ,LD)=0.334; the
joint LP drives β^LD to zero at all horizons), and `findings.md` records the decision to move
LD→vacancy to the appendix and calibrate τ externally at 3.1%/month.

**Consequence of the status quo:** ρ_s and σ_s lean entirely on Block M — the very moments
§5.2 argues cannot separately discipline δ vs. s propagation. That argument cuts against the
paper if it is not addressed head-on.

**Recommendation: keep Block B δ-only, and say plainly in §5.2 that s-shock propagation is
identified from second moments plus the externally calibrated τ, and that this is a
limitation.** Feeding a contaminated LD IRF into the likelihood would be worse than omitting
it.

---

## 🟠 D5. Comparison B — is the χ^c channel adequately exercised?

⚠️ **Rewritten September 5, 2026.** The old framing here described the *superseded*
Comparison B design. `run_solution_endog_exit.jl` was redesigned on May 22, 2026 to equate
δ̄ *and* τ across specifications, so Comparison B no longer demonstrates an "exposure effect"; it
isolates the χ^c margin, and its finding is that the endogenous-exit specification **recovers faster**
because a higher δ_e means a higher replacement rate. The draft §5.3 reports the new design.
The old "same dest_ann, different δ̄" version survives only in the runner's header as a
possible appendix robustness check.

**The live question** is whether the χ^c channel is quantitatively interesting at the
calibrated exit elasticity. Broer et al. (IER 2025) prefer an implied ψ near 1; PATH B at
`dest_elast_target = 5.0` gives ψ ≈ 0.033. Under [D2](decisions.md) the elasticity becomes an
*estimated* parameter, which largely dissolves this decision — the data will say how large the
channel is, and the answer is reportable rather than assumed.

**Recommendation: fold this into D2.** Estimate `dest_elast_target`, report the implied ψ and
exit elasticity in the posterior table next to the Broer et al. value, and drop the separate
ψ≈1.0 recalibration exercise unless a referee asks.

---

## 🟡 D6. Primary LP specification — separate or joint?

**Standing recommendation (unchanged since June 2026): separate LP as primary, joint as a
robustness check for δ.** The δ vacancy IRF is qualitatively stable across both (trough
−0.582 separate vs −0.535 joint, overlapping bands). The joint spec's zero β^LD reflects
collinearity, not confirmed reposting — `principles.md` §14 already says this.

Treat as settled unless the industry-credit control (D7) changes the correlation.

---

## 🟡 D7. Industry credit control (Route B) — implement or drop?

Would test the leading hypothesis for the residual r(δ,LD)=0.334 and could validate the
separate LPs. Requires BofA/ICE OAS by sector (Bloomberg) or Compustat leverage; data access
is the binding constraint.

**Recommendation: drop for this paper and note it as a limitation.** The severity placebo
already carries the external-validity argument, and the QU placebo failure is not fixable
this way either (`principles.md` §11).

---

## 🟡 D8. v3 residualization (+ monetary-policy sensitivity φ_j × ΔMP_t)

`part2b_residualize_shocks_v3.py` is written but has never been run. The draft already
promises it as an appendix robustness.

**Recommendation: run it once and report it, or delete the promise from the draft.** Cheap
either way; right now the draft writes a cheque the pipeline has not cashed.

---

## ⚪ D9. Free-entry steady state in `steady_state_checks.jl`

The solver lands on the θ=1.80 branch instead of ≈0.51; the diagnosed fix is a bracketed
solve over `(log(0.1), log(0.6))`. Purpose is expositional only — showing ξ_inv has no handle
on z-shock amplification — and Comparison D already makes the substantive point.

**Recommendation: fix if it takes an hour, otherwise drop.** Not on the critical path.

---

## 🟡 D10. What is the correct empirical target for the δ→u IRF shape?

**Original question (Sept 21).** Is the Bartik LP peak horizon (h=17–20) a propagation
fact or a persistence artifact?

**E8 result (Sept 22–23).** The lagged-instrument test was run (p=4, 8, 12, 16). Key findings:

1. **Instrument persistence confirmed.** Within-state autocorrelation of B̃^δ (after time
   FEs) has median ρ=0.91 at lag 1, 0.79 at lag 8. 90–100% of states significant at all lags.
2. **δ→u peak horizon is not stable.** Shifts from h=17 (baseline) to h=10 (p=4–12) to
   h=13 (p=16). Peak magnitude ranges 0.95–1.86 pp — the p=16 coefficient *exceeds* the
   baseline (1.86 vs 1.71). This drift is symptomatic of extracting an innovation signal
   from a near-unit-root process: the identifying variation shrinks with each added lag.
3. **Non-monotone shape emerges:** significant at h=0–2, insignificant at h=3–8, second
   rise at h=9–13. Stable across all four lag orders, but unprecedented in theory.
4. **δ→v early trough at h=4–6** is robust (−0.55 to −0.65 pp, p<0.01 at all lag orders).
5. **Asymmetry:** unemployment peak location unstable, vacancy early trough rock-solid.

**The core tension.** The baseline LP is not robust to persistence controls — the peak
horizon moves substantially. But the augmented LP produces a shape no known model
generates and whose peak location and magnitude drift with the lag order. Neither spec
gives a clean IRF target for estimation.

**Remaining concern: BED Deaths as a lagging indicator.** An establishment "death" in
BED data requires zero employment for two consecutive quarters. The official death
therefore lags the economic process. This measurement lag could contribute to both the
baseline and augmented IRF shapes.

**Options going forward:**
- **(a)** Target only robust short-horizon features (h=0–2 for u, h=4–6 for v) plus scale.
  Treat the peak horizon as uninformative. Conservative but defensible.
- **(b)** **Baseline-to-baseline matching.** Run the same LP specification on model-simulated
  panel data. If the model's shock persistence (ρ_δ=0.592) produces comparable Wold
  contamination to the data's instrument persistence (ρ=0.91), the comparison is
  apples-to-apples and the full baseline IRF is a valid target — both sides are "wrong"
  in the same way. **This preserves the original Bayesian IRF matching design.** Must be
  tested: if the model LP peaks at h=1–2 even in the Wold representation, the asymmetry
  invalidates this approach.
- **(c)** Hybrid: match full baseline path but verify with model-side LP that the Wold
  contamination is comparable; horizon-reweight Block B if not.
- **(d)** Present the augmented LP and defend the non-monotone shape — requires finding
  a theoretical or econometric explanation for the dip at h=3–8. High burden of proof.

**Recommendation (tentative):** Test option (b) first — it has the highest payoff (preserves
the full IRF matching design) and requires only one model simulation with the LP run on
top. If it fails (model LP peaks at h=1–2 even in Wold form), fall back to option (a).

**Status:** 🟡 partially resolved. E8 ran (p=4, 8, 12, 16) and identified the problem.
The estimation design can proceed under option (a) immediately, or wait for the option (b)
test to determine if the full IRF path is usable.

---

## Settled — do not re-open without a reason

| # | Decision | Settled | Why |
|---|---|---|---|
| S1 | HP λ=1,600, not λ=100,000 | May 19, 2026 | λ=100k flips the sign of cor(δ,u) and cor(δ,v); it is the outlier at 7/18 cells against HP-1600 and Hamilton. Over-smoothing artifact. `app:filter_robustness` |
| S2 | BED **Deaths** (dataclass 08), not Closings (06) | May 14, 2026 | Deaths = permanent exits, which is the model object. Dropping temporary shutdowns cut r(δ,LD) from 0.423 to 0.334 |
| S3 | Vacancy outcome = vacancy **rate** in pp, levels difference | — | Symmetry with the unemployment rate. Never log. `principles.md` §2 |
| S4 | LD is the primary s-type instrument; QU is the (failing) placebo; TS is appendix only | — | `principles.md` §3, §11 |
| S5 | LOO national shock rates; 2006 QCEW base-year shares; COVID cap 2019Q4 | — | `principles.md` §5, §6, §9 |
| S6 | τ calibrated externally at 3.1%/month, not identified from the LD IRF | — | JOLTS LD mixes reposting with deliberate reduction and partial closure; unfixable without firm-level data. `findings.md` §LD Vacancy Interpretation |
| S7 | Two-block posterior with degrees-of-freedom normalization | June 2026 | Stops the 42-dim IRF block swamping Block M; each quadratic form is ≈χ²(K) with mean K. `eq:posterior` |
| S8 | Filter asymmetry resolved on the model side, not the data side | June 2026 | m(θ) simulated-then-filtered; β(θ) read unfiltered off the state space. Canova-Ferroni (2011) |
| S9 | ε = 4.3 fixed externally | — | Determinacy is not binding: BK needs only ε>1, Prop. 5 Part 2 needs ε>2, both satisfied. See [`estimation_design.md`](estimation_design.md) §Determinacy |
