# Decision Register
**Last updated:** September 5, 2026 · **Branch:** `Organize_Project_State_Estimation`

Every open decision that must be settled before estimation, plus the ones already settled
that agents keep re-litigating. One entry = one decision. When you settle one, move it to
§Settled with the date and the reason — do not delete it, because the draft's prose depends
on knowing *why*.

Status key: 🔴 blocks estimation · 🟠 blocks a paper section · 🟡 improves the paper · ⚪ optional

---

## 🔴 D1. Which δ̄_e? BED Deaths (2.83%/yr) or Jaimovich-Floetotto (7.54%/yr)

### Canonical derivation — read this first

```
1. SOURCE   BLS BED, dataclass 08 "deaths" = establishments with zero employment
            in the third month of 4 consecutive quarters (i.e. no re-entry within
            a year). Total private. Published RATE series:
            BDS0000000000000000110008RQ5
              = "percentage employment lost from establishment deaths", quarterly.
            Already employment-weighted, with employment in the denominator —
            the same denominator as tau, so the two are directly comparable.

2. WINDOW   Sample mean over 2001Q1-2019Q4 (matches the LP sample; pre-COVID per
            principles.md 6)                                   ->  0.718 %/qtr

3. FREQ     Quarterly -> monthly. Flows are additive within the quarter
            (principles.md 9), so divide by 3                  ->  0.2393 %/mo

4. ANNUAL   Compound: 1-(1-0.002393)^12                        ->  2.83 %/yr
            This is `dest_ann`. steady_state.jl Stage 1 inverts it exactly:
            delta_e = 1-(1-dest_ann)^(1/12) = 0.2393 %/mo  [round-trips]

5. RATIO    delta_e / tau = 0.2393 / 3.10                      ->  0.077
            tau = 3.10 %/mo (TARGETS.sep, Shimer/JOLTS total separation rate)
```

**Target: `dest_ann` = 0.0283, δ_e/τ = 0.077.** On the longer 1993–2019 window:
0.812 %/qtr → 0.2707 %/mo → `dest_ann` = 0.0320, δ_e/τ = 0.087. Prefer the ratio as the
primary target and derive `dest_ann` from it (see recommendation below).

Three judgment calls are embedded, all flagged below: the **window** (2001–2019 vs 1993–2019,
worth ~0.01 in the ratio), **deaths not closings** (settled — closings include temporary
shutdowns), and **establishment vs firm** (open — BED is establishment-level, BDS is firm-level).

---

### Where the target sits in the literature

| | %/qtr | %/yr | δ_e/τ |
|---|---|---|---|
| Coles-Kelishomi (δ_e = τ) | 9.30 | 32.3 | 1.000 |
| BGM 2012 (δ = 0.025/qtr) | 2.50 | 9.6 | 0.269 |
| Gabrovski-Silva JEDC preferred range | 1.54–2.60 | **6–10** | 0.165–0.280 |
| Code now (`dest_ann = 0.0754`) | 1.95 | 7.6 | 0.210 |
| **Our target (BED deaths, 1993–2019)** | **0.81** | **3.2** | **0.087** |

⚠️ **This creates a new inconsistency with the paper's own introduction.** `Draft.tex:421`
cites Gabrovski-Silva approvingly: annual destruction rates "in the range of 6–10%" generate
dynamic correlations that closely match the data. Our target is **3.2%, about half the floor
of that range** — not marginally lower but outside it. The intro currently endorses a range
the calibration would abandon. This must be addressed explicitly, not quietly.

**Why the model may nonetheless fit at 3.2% where GS needed 6–10%.** The offsetting forces
are exactly the channels this model adds and GS lacks: the endogenous-exit margin ω_δ (a
profit-sensitive destruction rate that supplies cyclical exit without requiring a high
*exogenous* rate), the variety externality ρ(N) → w^int → JCC, and finitely elastic entry ξ.
This is a genuinely attractive story — the richer model matches the data at an empirically
defensible δ where a leaner one needed an inflated one — but **it is a conjecture until the
estimation runs.** Do not write it into the draft as established.

**The countervailing measurement argument — take seriously.** The model's δ destroys *product
lines*; BED measures *establishment* deaths. A surviving multi-product firm that discontinues
a product line is δ in the model and invisible in BED. So the BED deaths rate is arguably a
**lower bound** on δ, not a point estimate, and something between 3.2% and 10% could be
defensible on measurement grounds alone. (This cuts against the D1 correction and partly
rehabilitates a higher number — though not via the JF route, which remains an arithmetic
error regardless.)

**Preferred resolution: move δ_e into Θ_e and estimate it**, with a prior anchored at the BED
deaths rate as a lower bound and enough mass to reach the BGM/GS range. Three advantages: it
converts the tension into a result rather than an assumption; the posterior for δ_e becomes a
reportable finding directly comparable to CK, BGM, and GS; and it lets the data adjudicate the
product-line-vs-establishment gap that no measurement can settle. Fixing δ_e externally at
either 3.2% or 7.5% asserts an answer the paper is well placed to estimate.

**Diagnostic — RUN September 5, 2026, and it does not settle D1.**
`run_solution_delta_target.jl` compares dest_ann ∈ {0.0320, 0.0754, 0.0963}. Full results in
[`findings.md`](findings.md) §D1/M5. Summary:

- ✅ The mechanism ordering is confirmed: σ(θ)/σ(labor_prod) = 0.48 / 0.85 / 1.11, monotone in
  δ_e, exactly as Comparisons A and B predict.
- ❌ **All three specifications miss the amplification target (11.70) by an order of magnitude, and all
  three put the δ→u peak at h = 1 quarter against the LP's h = 17–20.** δ_e is not the binding
  constraint, so the decision rule above does not discriminate.

**D1 therefore stays open, but its priority drops.** Three prior problems have to be fixed
before any δ_e choice can be evaluated: the observable mapping for labor productivity (model
`labor_prod` is 3× too volatile and 0.9975-correlated with N^e — it is tracking entry, not
technology), the fact that Θ_e are hand-set rather than estimated, and the hump-shape gap in
the δ→u response. See `findings.md` for the ranked list.

The **measurement** argument for the BED number is unaffected by any of this — 0.0754 remains
an arithmetic error regardless of how the model fits. If a number must be written into
`steady_state.jl` today, use 0.0320 (or 0.0320's 1993–2019 sibling) and revisit after the
observable mapping is fixed.

---

**The conflict.** Three places, two numbers:

| Where | Value | Implied δ_e/τ |
|---|---|---|
| Draft §5.1 (`sec:ck_comparison`) | δ̄ ≈ 1.079%/qtr = **4.25%/yr** (BED Deaths, emp-weighted) | 0.116 |
| Draft §5.2 external block + `tab:calib_targets` | δ_e^ann = **7.54%/yr** (Jaimovich-Floetotto 21% × τ=3.1%/mo) | 0.210 |
| `steady_state.jl:566` | `dest_ann = 0.0754`, commented `[BED Deaths, emp-weighted]` | 0.210 |

The code's *value* is Jaimovich; the code's *comment* is BED. Every mechanism figure
(Comparisons A–D) and the "δ̄_e/τ̄ ≈ 0.21" claim at `Draft.tex:1481` were produced at 0.210.

**Why it matters.** §5.1 states that setting δ from BED Deaths rather than from the
aggregate separation rate is one of the paper's two departures from Coles-Kelishomi. As
things stand, the paper claims that departure and does not implement it. Separately, the
entry-cushion coefficient in Proposition 5 Part 3 is δ̄_e/(r+δ̄_e), so the choice scales the
paper's central quantitative claim.

### Verified September 5, 2026 — the two numbers are not rival estimates; 0.0754 is a denominator error

**What the δ series actually is.** `part2_shock_rates.py` and `observables.py` both use BLS
BED **dataclass 08 (Deaths)**, series `BDS0000000000000000110008LQ5`: employment at
establishments with zero employment in the third month of **four consecutive quarters**
following their last quarter with positive employment. So yes — the instrument's δ counts
only establishments that do **not** re-enter within at least a year. Deaths ⊂ Closings;
closings include temporary and seasonal shutdowns that reopen. `observables.py` forms
`delta = deaths_employment_t / payems_t`, so it is employment-weighted with total employment
in the denominator — the same denominator as τ. **δ_e and τ are already commensurate.**

**What Jaimovich & Floetotto's 21% is.** JF Table 2, column 2 is the share of *gross job
losses* from closing establishments: `sum(L_C)/sum(L)`, BED elem0006 over elem0004. The
denominator is gross job losses — net employment declines at contracting and closing
establishments — **not** total separations. Total separations are far larger because they
include quits and replacement churn that produce no net job destruction. Note also that JF
use **closings**, not deaths.

**The arithmetic.** Published BED rates, total private, % of employment per quarter, fetched
from the BLS API (v1, series `BDS...11{0008,0006,0004}RQ5`):

| window | deaths | closings | gross job losses | δ_e/τ if deaths | δ_e/τ if closings | closings as % of gross losses |
|---|---|---|---|---|---|---|
| 1993–2019 | 0.812 | 1.315 | 6.681 | **0.087** | 0.141 | 19.7 % |
| 2001–2019 (LP sample) | 0.718 | 1.205 | 6.362 | **0.077** | 0.130 | 18.9 % |
| 2010–2019 | 0.603 | 1.073 | 5.815 | **0.065** | 0.115 | 18.5 % |

The closings share of gross losses (19.7 % over 1993–2019) reproduces both JF's ~21 % and
`part11`'s 19.5 %, which confirms the identification of their statistic. **The deaths share of
gross losses is only ~12 %.**

Two conclusions:

- **21 % × τ is not a valid construction of δ̄_e.** It applies a share whose denominator is
  gross job losses to the total separation rate, inflating δ̄_e to 1.95 %/qtr — more than
  double the measured deaths rate (0.72–0.81) and still well above the closings rate
  (1.21–1.32). Nothing in BED supports 1.95 %/qtr as an exit rate.
- **Even applied correctly, JF's share targets *closings*, not deaths.** Closings include
  temporary and seasonal shutdowns that reopen within a year — which are not permanent
  product-line destruction, and not what the Bartik instrument measures.

> ⚠️ **Correction, September 5, 2026.** An earlier version of this entry claimed
> "19.5 % × 5.38 = 1.05 %/qtr → 4.13 %/yr" as independent confirmation of the BED figure.
> That was wrong twice over: `BDS...110005RQ5` is *losses at contracting establishments*, not
> total gross losses (the repo labels it correctly in `BED_data_construct.py:16`); and the
> calculation is circular in any case, since JF's share was computed from that same base, so
> share × base merely recovers the closings rate. The 4.13 % figure has no standing — use the
> direct measurements in the table above.

**Why JF's choice was reasonable for them.** Closings, not deaths, is the natural object for
a paper about the cyclicality of the operating-establishment margin and markup variation, and
deaths require a four-quarter confirmation lag, so in 2008 the deaths series was both shorter
and published with substantial delay. Their choice is sound for their question; it just does
not transfer to a δ_e/τ ratio.

**Recommendation: target δ_e/τ ≈ 0.08 directly, from BED Deaths on the 2001–2019 LP sample.**
That is δ_e = 0.718 %/qtr = 0.239 %/month, δ_e^ann = **2.84 %**. Use 1993–2019 (δ_e/τ = 0.087,
δ_e^ann = 3.21 %) if the longer window is preferred for a steady-state target; the choice
between them is second-order, and both are far from 0.21.

Rationale beyond the arithmetic: deaths is the measurement the empirical section is already
built on (the δ Bartik instrument *is* dataclass 08), and a lower δ̄_e/τ̄ *strengthens* Prop. 5
Part 3, whose sufficient condition is "δ̄_e/τ̄ small". The §5.2 paragraph should cite the BED
Deaths rate directly and drop the JF routing entirely.

**Prefer targeting the ratio, not `dest_ann`.** δ_e/τ is the economically meaningful object —
it is the share of separations attributable to permanent firm exit, and it is what
Proposition 5 Part 3 and the Coles-Kelishomi comparison both turn on (CK set δ_e = τ, i.e.
δ_e/τ = 1). Making it the primary target and deriving `dest_ann` from it also makes the
departure from CK legible at a glance.

⚠️ **Why the draft's 1.079 %/qtr / 11.6 % is too high — resolved.** `observables.py` sets
`final = '2025-10-01'`, so its δ series runs to the end of the BED cache (2021Q4) and the mean
**includes the 2020 COVID quarters**. It also divides by PAYEMS (total nonfarm, including
government) rather than private employment, which biases the *level* the other way. A
steady-state target should be computed pre-COVID, consistent with the LP's 2019Q4 cap
(`principles.md` §6). Recompute on 2001Q1–2019Q4 before the number enters `tab:calib_targets`.

⚠️ **Open sub-question: establishments or firms?** BED measures *establishment* deaths, so a
multi-establishment firm closing one location counts as exit. The model's object is a firm /
product line. The superseded 0.940 %/qtr figure came from BDS, which is firm-level. If the
firm-level rate is materially lower, δ_e/τ falls further — again in the direction that helps
Prop. 5. `bds2022.csv` and `BDS_Extension.md` are in the repo; worth a check before finalizing.

The cost of switching is re-running every mechanism figure.

**If chosen, propagate to:** `steady_state.jl:566`; the §5.2 external-block paragraph;
`tab:calib_targets` row 1; the ≈0.21 claims at `Draft.tex:1481` and in `app:comparison_D`;
Comparisons A–D figures; the calibration table in [`data_and_files.md`](data_and_files.md).

---

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
