# Decision Register
**Last updated:** September 23, 2026 · **Branch:** `instruments_LP_coefficient`

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

### ⚠️ Reframed September 23, 2026 — the measurements bracket δ_e, they do not identify it

The literature review in
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md)
changes how this entry should be read. No measurement identifies the model's δ:

- **BED deaths (3.2 %/yr)** is a directly observed **lower bound**. It misses product lines
  and positions destroyed inside surviving firms.
- **Permanent layoffs (GS, ~10 %/yr)** is an **upper bound**. Recall rates measure
  worker-job attachment, not position survival, so a position refilled by a different worker
  counts as permanent.
- **Product destruction (BGM, 8.8–10 %/yr)** is an upper bound for the exit interpretation:
  roughly two-thirds of it occurs at continuing firms. Netting that out gives ≈2.7 %/yr from
  exiting firms, close to our 3.2 %. **BGM and this paper do not disagree about measurement.**
- **Shao-Silos (0.9 %/month)** identifies nothing about destruction; it is a residual from an
  unemployment-rate target.

So the literature supplies an interval of roughly **[3.2 %, 10 %]/yr**, not a point. δ_e is
properly *prior information*, not a calibrated constant.

**This does not reopen D1 now, and the sequencing matters.** With no reposting margin in the
model, an estimated δ_e would be pulled upward by the estimator trying to rescue cor(u,v),
landing on a value that reflects the missing margin rather than the destruction rate. That is
exactly the absorption problem the section below describes. So:

- **Now:** hold δ_e at the observed lower bound, and describe it in the draft *as* a lower
  bound rather than as the truth.
- **Once [D11](#-d11-does-the-reposting-margin-enter-this-paper) is resolved and the margin
  exists:** move δ_e into Θ_e with a prior supported on [3.2 %, 10 %], mass near the lower
  end. The posterior then reports the size of the middle margin, which is a result worth
  having.

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

⏳ **None of these are done as of September 23, 2026.** `steady_state.jl:613` still holds
`dest_ann = 0.0754` behind a warning comment, so every mechanism figure and every §5.3 number
is still at δ_e/τ = 0.210. The work is tracked as the "Step 0 cascade" in
[`pending_tasks.md`](pending_tasks.md).

- `steady_state.jl:613` `TARGETS.dest_ann`: 0.0754 → 0.0320, and correct the source comment.
- Re-run Comparisons A–D, regenerate the eight figures, re-run `mechanism_stats.jl`, and
  update §5.3 against its output.
- Draft §5.2 external block and `tab:calib_targets` row 1.
- The "δ̄_e/τ̄ ≈ 0.21" claims become ≈ 0.087. There are **four** of them, verified
  September 23, 2026: `Draft.tex:1481` (`sec:mechanism`), `:2787` (§5.3), `:3646` and
  `:3921` (both in `app:comparison_D`). Prop. 5 Part 3's condition is "δ̄_e/τ̄ small", so the
  lower value *strengthens* it.
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

### ⚠️ Overlap with D10 — noted September 23, 2026

**Option 1 above and D10 option (b), below, are the same exercise.** D10 option (b) — run the identical LP specification on
model-simulated panel data — *is* the GMNS "simulate the regression inside the model"
prescription. Task **E9** in [`pending_tasks.md`](pending_tasks.md) builds it.

Consequence: if E9 is built properly, option 1 stops being "probably out of scope" and
becomes available, and D3's recommendation should be revisited and likely upgraded from
option 2 (shape-only) to option 1. Conversely, if E9 is scoped narrowly — just to answer the
Wold-contamination question — it will not deliver a D3-grade model-implied cross-regional
coefficient. **Scope E9 with both purposes in mind.**

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

### ⚠️ [D11](#-d11-does-the-reposting-margin-enter-this-paper) cuts against this — September 23, 2026

The premise above is that LD is *contaminated*, because it mixes match dissolution with
reposting against deliberate reduction and partial closure. Under the three-margin reading
that mixture is not contamination. It is exactly the object a reposting margin
parameterizes, and the persistently negative LD→vacancy IRF is the middle margin showing up
in the data.

**If D11 resolves in favor of the reposting margin, this decision must be revisited**: the
LD→vacancy IRF becomes the natural identifying moment for the reposting rate and its
elasticity, and should enter Block B rather than the appendix. If D11 resolves against, the
recommendation above stands unchanged.

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

**~~Remaining concern: BED Deaths as a lagging indicator.~~ RETIRED September 23, 2026.**
The earlier entry said a BED death "requires zero employment for two consecutive quarters"
and inferred that the official death lags the economic process, contributing to the IRF
shape. **Both halves are wrong.**

1. **The rule is four quarters, not two.** BLS defines a death as an establishment with zero
   or no employment reported in the third month of *four consecutive quarters following the
   last quarter with positive employment*. Writing q₀ for the last quarter with positive
   employment: the establishment is a **closing** in q₀+1, and a **death** only if q₀+1
   through q₀+4 are all zero. That is the closing quarter plus three further confirmed
   quarters.
2. **The death is dated to the quarter of closure, not to the quarter of confirmation.**
   The four-quarter rule is a *confirmation requirement*, not a timing shift. A death
   occurring in 2008Q4 is recorded in 2008Q4; it simply could not be known until 2009Q4.

**Consequence for D10:** the recognition rule produces a publication lag of roughly a year
and vintage revisions. It does **not** displace the δ series in time and therefore **cannot
explain the δ→u peak at h = 17–20.** This candidate explanation is removed. E8's
instrument-persistence finding remains the live one.

Sources: [BLS BED concepts](https://www.bls.gov/opub/hom/bdm/concepts.htm); Sadeghi,
"The births and deaths of business establishments in the United States," *Monthly Labor
Review*, December 2008. ⚠️ Cite the Handbook, not the MLR article, if the rule is ever
stated in the paper: the MLR article describes the *proposed* 2008 methodology (zero in a
quarter plus the next four), which is not the rule the published series uses.

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

### ⚠️ Design caveat for E9 — added September 23, 2026

Option (b) compares two *different* kinds of persistence, and the test is only informative if
the simulation reproduces the right one.

- The data-side persistence is in the **instrument**: ρ(B̃^δ) = 0.91 at quarterly lag 1. It
  comes from slow-moving industry composition times serially correlated *industry-level*
  death rates.
- The model-side persistence is in the **shock**: ρ_δ = 0.592 monthly, i.e. ≈ **0.21
  quarterly** — a factor of four lower.

If E9 simulates 50 independent draws of the aggregate model and builds a shift-share
instrument on top, the simulated instrument will inherit ρ ≈ 0.21, the model LP will show no
buildup, and option (b) will be rejected for a reason that has nothing to do with the model's
propagation. To be informative the simulated panel needs **industry-level δ processes with
their own serial correlation and heterogeneous state exposure**, so that the simulated
instrument's autocorrelation is comparable to 0.91. Match the instrument's persistence, not
the shock's. See also the D3 overlap note above.

**Status:** 🟡 partially resolved. E8 ran (p=4, 8, 12, 16) and identified the problem.
The estimation design can proceed under option (a) immediately, or wait for the option (b)
test to determine if the full IRF path is usable.

---

## 🔴 D11. Does the reposting margin enter this paper?

**Opened September 23, 2026.** Full argument:
[`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).

**The problem.** The model's unconditional Beveridge correlation is **+0.995** against
**−0.804** in the data. Two findings locate the cause:

1. **Raising δ̄_e helps monotonically but cannot get there.** cor(u,v) improves from +0.995 to
   **+0.887** as δ̄_e/τ̄ rises from 0.087 to 0.271. The direction is right and the mechanism is
   real. But that is 0.108 of a gap of roughly 1.8, across the entire defensible range.
   The relationship is also strongly convex: at δ_e = τ the identity forces s̄ = 0 and
   cor(u,v) must equal the s-silenced −0.911, so the average slope to that endpoint is about
   −2.1 per unit of δ̄_e/τ̄ against a local slope near −0.59 in the measured region.
   **Almost all of the leverage δ_e offers lies in a region the measurement rules out.**
2. **Silencing the s shock fixes it completely** (−0.911). The Beveridge mechanics are sound.
   What breaks them is the s shock under *costless reposting*.
3. **δ_e and λ̄ are complements working through one channel** — the share of separations that
   destroy a vacancy. So cor(u,v) alone cannot separate them, which extends the
   λ̄-versus-σ_s identification concern below to a three-way problem and argues for holding
   δ_e fixed while λ̄ is estimated rather than freeing both.

**Why no antecedent has this problem.** BGM, Shao-Silos, and GS are all driven by a
productivity shock as the single source of exogenous volatility; they hold δ and s fixed. CK
does shock separations but has one margin with δ_e ≡ τ, so every separation destroys a
vacancy and the Beveridge curve is mechanical.

**Why we cannot avoid it.** CK's contribution is that *separation shocks* drive Beveridge
dynamics. Decomposing that one shocked margin into δ and s obliges us to shock **both**.
Shocking δ while holding s fixed would assert the match-dissolution component is acyclical,
which is false (σ_s = 0.0854, cor(cycle_s, cycle_τ) = 0.997) and is an arbitrary restriction
on exactly the margin CK identify as the driver. Three shocks is required, not chosen.

**The missing middle.** Three margins exist in the data; the literature and this paper fold
the middle one into opposite neighbors.

| Margin | Destroys product line? | Destroys position? | Rate |
|---|---|---|---|
| Establishment death | Yes | Yes | 3.2 %/yr |
| **Permanent position cut at a surviving firm** | **No** | **Yes** | *the missing middle* |
| Match separation, reposted | No | No | remainder |

GS and Shao-Silos fold the middle into δ, overstating variety destruction. We fold it into s,
understating vacancy destruction. Our own LD→vacancy IRF (persistently negative, trough
−0.45 pp) is that margin showing up in the data.

**The fix, and why it must be a cost not a constant.** A reduced-form reposting fraction λ
will not survive review. It is the only mechanical margin in a model where entry, exit, and
vacancy creation all optimize, and — decisively — it is **not invariant to the shocks being
studied**. If reposting is a decision, firms repost less when position values fall, i.e. in
recessions, which is exactly when the Beveridge curve shifts. A constant λ holds it fixed
there.

A *homogeneous* reposting fee also fails: all positions being identical, the decision is
all-or-nothing, λ ∈ {0,1}. **An interior reposting rate requires dispersion in the reposting
cost**, exactly as the firm-level exit margin requires F to deliver an interior χ^c. So:
draw c from a distribution G, repost iff Q_t ≥ c, giving λ_t = G(Q_t), procyclical, with
countercyclical position destruction. λ = 1 is nested.

**Identification.** The LD→vacancy IRF identifies both the level λ̄ (sign and magnitude) and
the elasticity (shape across horizons). This gives the LD arm of Block B a structural job for
the first time and cuts against [D4](#-d4-does-the-sld-irf-enter-the-quasi-likelihood-block)
/ S6, which currently retire that IRF.

**Recommendation: estimate, do not calibrate.** No clean data object exists — recall rates
measure worker attachment not position survival, BED contraction nets within the quarter and
is frequency-dependent, JOLTS reports a stock. λ̄ adjudicates a headline moment, so fixing it
assumes the answer. And the elasticity cannot be calibrated at all; it exists only inside the
model. Note this is the mirror of the D1 logic: δ_e is fixed *because* it is cleanly measured,
while λ *is* the misspecification and is the right thing to free.

**Costs, honestly.** Prop. 5 Part 1 becomes conditional on λ̄ above a threshold λ̄*, and the
proof must be redone (`lem:vpre` and the Part 2 Jacobian both touch the channel). More
seriously, with λ̄ < 1 both shocks lower vacancies, so the headline δ–s asymmetry becomes
quantitative rather than a sign result. Mitigations: the *variety* asymmetry is untouched,
since only δ moves N and so only δ produces a persistent outward shift; and λ̄* becomes a
reportable, testable threshold.

**Open sub-question:** verify λ̄ and σ_s are separately identified on simulated data before
committing. Both move σ(v) and cor(u,v). If not, λ̄ must be fixed on a reported grid.

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
