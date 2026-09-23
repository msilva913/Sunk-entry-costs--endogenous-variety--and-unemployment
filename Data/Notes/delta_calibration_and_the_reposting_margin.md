# Calibrating δ: product lines, positions, and the reposting margin
**Written:** September 23, 2026 · **Branch:** `instruments_LP_coefficient`
**Status:** analysis and recommendation. Nothing here is settled policy.

**Purpose.** δ̄_e = 3.2 %/yr looks low against every antecedent, and the model's unconditional
Beveridge correlation has the wrong sign. These two facts are usually treated as one problem.
They are not. This note establishes what δ measures, argues the case for our choice on the
merits, compares the literature honestly, shows why raising δ̄_e cannot fix the Beveridge
curve, and locates the fix in a reposting margin that should be modeled as a cost.

---

## 1. What δ is being asked to do

A single parameter δ_e carries two jobs.

1. **Variety.** The fraction of product lines withdrawn each period, driving ρ(N) and the
   N → ρ → w^int → θ chain.
2. **Labor.** The fraction of employment destroyed by exit, entering τ through the identity
   τ = δ_e + (1−δ_e)s, and destroying pre-committed vacancies through `lem:vpre`.

Under DS-CES symmetry these coincide. In the data they do not, for two independent reasons.
Dying establishments are smaller than average, so a count-weighted rate exceeds an
employment-weighted one by roughly 3.7× (BDS firm deaths: 8.33 % of firms, 2.26 % of
employment). And product lines die inside surviving firms, which no establishment-exit
measure observes.

Both wedges push the same way. The variety margin wants a **larger** number than the labor
margin, so any single-δ_e calibration must choose which job to serve.

### 1.1 The case for serving the labor margin

This is a choice made on the merits, not an inherited convention. Five arguments, roughly in
order of force.

**(a) The central proposition depends on the employment-weighted rate, not the count.**
Proposition 5 Part 3 turns on `lem:vpre`: a δ shock destroys pre-committed vacancies at
exiting firms. The vacancies destroyed are *those firms' vacancies*. Dying establishments are
small and hold correspondingly few positions. Feeding the model a count-weighted δ_e would
overstate vacancy destruction by the same factor of roughly 3.7, and vacancy destruction is
precisely the object the paper's headline asymmetry rests on. Getting this margin right is
not a side concern; it is the mechanism.

**(b) A count-based δ_e would contradict the paper's own measurement.** The empirical section
measures employment lost to establishment death directly, from BED, employment-weighted. If
the calibration set δ_e near 10 %/yr, the model would attribute roughly three times as much
job destruction to firm exit as the paper's own data section reports. Sections 4 and 5 would
describe different economies.

**(c) The targeted moments are labor moments.** Block M is built from {u, v, s, f, δ, N^e, z}
and Block B from the δ→u and δ→v impulse responses. Every moment the estimator sees is a
labor market object, and the decomposition of τ is first-order for all of them. The variety
channel enters these moments only indirectly, through ρ(N).

**(d) The instrument is employment-weighted.** The Bartik shock is built from
employment-weighted BED deaths, so the IRF being matched is the response to an
employment-weighted destruction shock. If the middle margin discussed in §5 transmits
differently from outright death, which is plausible since the firm survives and may repost
later, then the LP identifies the death component specifically. The calibrated δ_e should
match what the LP identifies.

**(e) It is the conservative error.** Understating variety destruction understates
amplification, which makes the paper's quantitative claims a lower bound. Overstating job and
vacancy destruction would inflate the δ–s asymmetry, which is the result the paper is
selling. Given a choice between two biases, we take the one that runs against us.

The cost is real and is now stated in §5.2 of the draft: our results are a lower bound on the
amplification contributed by the variety channel.

---

## 2. What we currently do

`dest_ann = 0.0320`: BED establishment deaths, employment-weighted, 1993–2019.

| Quantity | Value |
|---|---|
| δ̄_e | 3.2 %/yr, 0.271 %/month, 0.812 %/qtr |
| τ̄ | 3.1 %/month, 9.30 %/qtr |
| **δ̄_e/τ̄** | **0.087** |
| Deaths as a share of BED gross job losses | 12.0 % (closings: 19.5 %) |

The measure is unambiguous. An establishment with zero third-month employment for four
consecutive quarters is gone, and the job losses attributed to it are permanent. It is
observed rather than inferred. That is its virtue, and §6 argues it is also its limitation.

---

## 3. How the literature calibrates the same object

| Paper | δ (annual) | δ/τ | Basis | Shocks |
|---|---|---|---|---|
| Coles & Kelishomi (2018) | 32.3 % | **1.00** | δ_e ≡ τ by construction; one margin | separation |
| Bilbiie, Ghironi & Melitz (2012) | 10 % | n/a | Product destruction incl. within-firm churn | z only |
| Shao & Silos (2013) | ~10.8 % | **0.26** | Residual from Σ = 0.035 after choosing s to hit u ≈ 5.4 % | z only |
| Gabrovski & Silva (JEDC) | 10 % | **0.26** | JOLTS layoffs and discharges net of recalls | z only |
| **This paper** | **3.2 %** | **0.087** | BED establishment deaths, employment-weighted | z, δ, s |

### Coles & Kelishomi

Sets δ_e = τ: every separation destroys the position. This delivers the Beveridge curve
mechanically, because nothing is ever reposted. It is the upper bound, it is not defensible
as measurement, and relaxing it is why this paper exists. But note that CK is also the only
antecedent that shocks the separation margin at all, which matters in §4.

### Bilbiie, Ghironi & Melitz

BGM set δ = 0.025 quarterly, stating it as a 10 % annual rate holding both as a share of
products and as market share, and validate against a minimum product destruction rate of
8.8 % measured as a market share. Two observations.

First, the 10 % is **not** a count-weighted artifact, contrary to the natural guess. BGM claim
it simultaneously as a product share, a market share, and (footnote 12) a job destruction
rate. Symmetry makes the three identical in their model.

Second, and decisively, their target covers product destruction **by existing and exiting
firms**. Over five years it accounts for 44 % of output, of which 30.4 points occur at firms
that continue operating. Netting that out leaves roughly 2.7 %/yr from exiting firms, close to
our 3.2 %.

**BGM and this paper do not disagree about measurement.** Their δ is larger because a single
parameter in a model without a separation margin must carry churn inside surviving firms.
This is the cleanest defense of our number, and it is now in §5.2.

### Gabrovski & Silva

GS reason that δ is the loss of a business opportunity, so it should be calibrated from job
loss due to firm exit *and product obsolescence*. Lacking a direct series they take JOLTS
layoffs and discharges (1.4 %/month), apply recall rates from Fujita-Moscarini and Lam-Qiu
(30–40 % recalled), and obtain 9.84–11.5 %/yr, settling on 10 %.

**The concept is right and the criticism of our measure lands.** GS state explicitly that
there may be δ-type destruction at firms neither exiting nor closing production lines.
Establishment deaths miss all of it.

**The implementation overshoots, for a reason GS do not address.** Layoffs net of recalls
measures *worker–job attachment*, not *position survival*. A worker never recalled may be
replaced by a new hire into the same position, which is a separation with reposting and an s
event in the model. GS themselves note that firing for cause belongs with the separation
shock, yet the recall adjustment does not remove it. Their 60–70 % permanent share is an
upper bound on positions destroyed, not an estimate of it. Their own reported BDS figure of
about 14 %/yr for separations at exiting *and shrinking* establishments, which they call
likely an over-estimate, points the same way.

Our 3.2 % is too narrow for the concept. Their 10 % is too broad for the measurement.

### Shao & Silos

Share our accounting exactly: Σ = τ + (1−τ)s with τ the no-repost rate. They set Σ = 0.035
from Shimer, choose s so steady-state unemployment matches the free-entry benchmark (≈5.4 %),
and take τ = 0.009 as the residual. The split is disciplined by the unemployment rate, not by
destruction data. They also assume one job per firm, so exit destroys exactly one position by
construction.

Their τ/Σ ≈ 0.26 coincides with GS's, but this is not corroboration. One number is a residual
from an unemployment target; the other an inference from recall rates.

### Jaimovich & Floetotto

Not a calibration antecedent. JF have no unemployment margin and no separation shock. Their
contribution here is measurement: the 21.17 % closings share of gross job losses over
1992Q3–2005Q2. They could not have used deaths, first published May 2009. Applying their
gross-job-loss share to a *separation* rate is what produced the erroneous
`dest_ann = 0.0754`.

---

## 4. Why we shock both margins, and why that creates the tension

### 4.1 The consistency requirement

CK's contribution is that **separation shocks** drive Beveridge curve dynamics. This paper
accepts that and rejects CK's conflation of firm destruction with match dissolution. But
decomposing one shocked margin into two parts carries an obligation: **both parts must be
shocked.**

Shocking δ while holding s fixed would assert that the match-dissolution component of
separations is acyclical. That is false in the data — the calibrated s process has
σ_s = 0.0854 and ρ_s = 0.874, and cor(cycle_s, cycle_τ) = 0.997 — and it would be an arbitrary
restriction on exactly the margin CK treat as the driver. The decomposition would not be a
decomposition. It would be a relabeling with an unmotivated acyclicality assumption on one
piece.

Shocking s alone discards the δ mechanism and the paper with it.

**So three shocks is not a discretionary enrichment. It is what taking CK seriously while
refusing CK's conflation requires.**

### 4.2 The tension this creates

With δ̄_e/τ̄ = 0.087, roughly 91 % of separations are match separations. Under costless
reposting (Prop. 5 Part 1) an s shock raises u and v together, so the dominant separation
margin pushes the Beveridge correlation the wrong way:

| | BED (0.0320) | code (0.0754) | BGM (0.0963) | **s silenced** | data |
|---|---|---|---|---|---|
| cor(u, v) | **+0.995** | +0.947 | +0.887 | **−0.911** | **−0.804** |

**Raising δ̄_e does not fix this.** At BGM's 9.6 %/yr — three times our value, above both GS
and Shao-Silos — cor(u,v) is still +0.887. Recalibrating δ upward buys almost nothing on the
moment it is supposed to buy.

**Silencing the s shock fixes it completely**, at −0.911. The Beveridge mechanics are sound.
What breaks them is the s shock under costless reposting.

### 4.3 Why no antecedent faces this

**None of the single-margin papers shocks separations at all.** BGM, Shao-Silos, and GS are
driven by a productivity shock as the single source of exogenous volatility. GS hold δ and s
fixed and match the Beveridge curve with z alone, exactly as our model does with s silenced.
CK does shock separations, but has only one margin and sets δ_e = τ, so every separation
destroys a vacancy and the Beveridge curve is mechanical.

The tension is therefore not a defect of our calibration. It is the first substantive thing
the decomposition reveals: **once you admit that most separations are not firm exits, and you
shock that margin as consistency requires, the costless-reposting assumption inherited from
the literature becomes untenable.** That is a finding, and the paper should present it as one.

---

## 5. The missing middle

The literature and this paper make opposite errors about the same margin.

| Margin | Destroys product line? | Destroys position/vacancy? | Rate |
|---|---|---|---|
| Establishment death | **Yes** | **Yes** | 3.2 %/yr |
| **Permanent position cut at a surviving firm** | **No** | **Yes** | *the missing middle* |
| Match separation, reposted | No | No | remainder |

GS and Shao-Silos fold the middle into δ, overstating variety destruction: a surviving firm
that cuts a position does not withdraw its product. We fold it into s, understating vacancy
destruction, which is why cor(u,v) is positive.

**The middle margin is visible in our own data.** The LD→vacancy IRF is persistently negative,
trough −0.45 pp at h = 20, and `findings.md` explains it as JOLTS LD mixing match dissolution
with reposting, deliberate workforce reduction, and partial establishment closure. The last
two are the middle margin. [D4](../context/decisions.md) and S6 currently treat this as
contamination and retire the IRF. Under the three-margin reading it is not contamination. It
is the margin showing up in the data.

---

## 6. Can calibration pin δ down, or should it shape priors?

This deserves a direct answer, because §3 is usually read as a menu to choose from.

**No measurement in §3 identifies the model's δ. They bracket it.**

- **BED deaths (3.2 %/yr)** cleanly measures a well-defined quantity: employment lost at
  establishments that permanently die. But that quantity is a **lower bound** on δ, because
  it misses product lines and positions destroyed inside surviving firms.
- **Permanent layoffs (GS, ~10 %/yr)** is an **upper bound**, because recall rates cannot
  distinguish a destroyed position from a position refilled by a different worker.
- **Product destruction (BGM, 8.8–10 %/yr)** is an upper bound for the exit interpretation,
  since roughly two-thirds of it occurs at continuing firms.
- **Shao-Silos (0.9 %/month)** identifies nothing about destruction; it is a residual from an
  unemployment target.

So the literature supplies an interval of roughly **[3.2 %, 10 %] per year**, with the lower
endpoint directly observed and the upper endpoint inferred. It does not supply a point.

**The honest implication is that δ_e is prior information, not a calibrated constant.** The
measurement gives support and shape for a prior, and the data should choose within it. Fixing
δ_e at an endpoint asserts precision the measurement does not have.

**But this should be sequenced, not done now.** [D1](../context/decisions.md) fixed δ_e partly
because freeing a parameter while known misspecification sits upstream lets it absorb the
blame. That argument is correct and currently binding: with no reposting margin, an estimated
δ_e would be pulled upward by the estimator trying to rescue cor(u,v), and would land at a
value reflecting the missing margin rather than the destruction rate. The recommendation is
therefore:

- **Now, with no reposting margin:** hold δ_e fixed at the directly observed lower bound.
  Report it as a lower bound rather than as the truth.
- **Once the reposting margin is in the model:** move δ_e into Θ_e with a prior supported on
  roughly [3.2 %, 10 %], mass concentrated near the lower end, since only that endpoint is
  observed. The posterior then reports what the data say about the middle margin's size, which
  is a result worth having.

This reframes §3 from a choice among rival calibrations into the construction of a prior, and
it is the more defensible position to take in print.

---

## 7. The reposting margin, modeled as a cost

### 7.1 Why a reduced-form λ will not survive review

The natural first pass is an exogenous constant λ ∈ [0,1], the fraction of match separations
after which the firm reposts, so the v law of motion's inflow becomes λ·s_t·(1−u_t). λ = 1
nests the current model.

A skeptical referee will say this should be written as a cost. **The referee is right, and
the objection is substantive rather than stylistic.**

- **It is inconsistent with the model's own logic.** Entry, exit, and vacancy creation are all
  optimizing margins. Reposting would be the only mechanical one, and it sits in the middle of
  the mechanism the paper is about.
- **It is not invariant to the shocks being studied.** This is the decisive objection. If
  reposting is a decision, firms repost less when the value of a position falls, which is
  in recessions. A constant λ holds the reposting rate fixed exactly in the state where the
  Beveridge curve shifts. The parameter would not be structural in the Lucas sense, and the
  paper would be estimating a reduced form while claiming a mechanism.
- **It gives up the most interesting implication.** Endogenous reposting makes position
  destruction countercyclical for free, which is the empirically relevant pattern and a
  genuine prediction rather than an assumption.

### 7.2 A homogeneous fee does not work either

Worth stating explicitly, because it is the referee's likely first suggestion. Suppose every
repost costs the same fee c_r. Then the firm reposts if and only if Q_t ≥ c_r, and since all
positions are identical the decision is all-or-nothing: λ_t ∈ {0, 1}. A homogeneous fee does
not generate partial reposting; it simply scales the value of a vacancy.

**An interior reposting rate requires dispersion in the reposting cost.** This is exactly why
the firm-level exit margin needs the distribution F to deliver an interior cutoff χ^c.

### 7.3 The formulation

On separation, the firm draws a reposting cost c from a distribution G, and reposts if the
value of a vacancy covers it:

```
repost  ⟺  Q_t ≥ c ,        λ_t = G(Q_t)
```

The reposting rate is now endogenous and **procyclical**, so position destruction is
countercyclical. The steady-state rate is λ̄ = G(Q̄), and the elasticity of λ with respect to
Q is the new structural content.

Three properties recommend this.

1. **It reuses machinery the model already has.** The endogenous exit margin draws
   continuation costs from F with a mass point p_0 and a Pareto tail, and compares them to a
   cutoff. G is the same construct one level down, at the position rather than the firm. The
   same estimation treatment, the same class of distribution, and the same interpretation of
   the mass point carry over.
2. **It restores an internal symmetry.** The model currently lets firms decide whether to
   continue a product line but forces them to refill every position. That asymmetry is hard
   to defend on its own terms.
3. **λ = 1 is nested**, at G degenerate below Q̄, so the cost of the current assumption is
   directly reportable.

### 7.4 Identification

Two objects: the level λ̄, and the elasticity of λ_t with respect to Q_t.

1. **The LD→vacancy IRF identifies both.** In the model an s shock raises v when λ̄ is high
   and lowers it when λ̄ is low, with a threshold λ̄* where the sign flips. The empirical
   response is negative and persistent, placing λ̄ below the threshold; its magnitude
   identifies the level, and its *shape* over horizons identifies the elasticity, because
   endogenous reposting makes the response deepen as Q falls. This gives the LD arm of Block B
   a structural job for the first time and repays [D4](../context/decisions.md)/S6.
2. **cor(u,v) in Block M** moves monotonically in λ̄. Powerful, but it is the moment being
   explained, so it should not be the sole source.
3. **BED gross job flows** bound λ̄ loosely. Gross job losses split into losses at closing
   establishments and at contracting establishments, and contraction is net position
   destruction at surviving firms. Treat as a bound, not a target, for the reasons in §7.5.

### 7.5 Calibrate or estimate?

**Estimate, with an informative prior.** The first reason is dispositive.

**1. There is no clean data object.** Reposting is not measured, and every proxy measures
something adjacent with a known bias. Recall rates measure worker attachment, not position
survival. BED contraction nets hires against separations within an establishment-quarter, so
it misses positions destroyed and recreated inside the quarter, and its value is
frequency-dependent while our model is monthly. JOLTS reports the vacancy stock, so
repostings are not separable from new positions.

**2. It adjudicates a headline result.** λ̄ determines whether the model reproduces the
Beveridge curve. Fixing it externally would assume the answer to the question being asked.
Note this is the mirror of the D1 logic and points the other way: δ_e is fixed *because* it is
cleanly measured and freeing it would let it absorb misspecification. λ *is* the
misspecification, so it is the right thing to free.

**3. The elasticity cannot be calibrated at all.** Even granting an external value for λ̄,
nothing in the data gives the responsiveness of reposting to the value of a position. That
object exists only inside the model and must come from the likelihood.

**Prior.** On λ̄, a Beta centered well below 1 with wide support. The recall-rate literature
is the only external anchor and it is an upper bound on non-reposting, so it should shape the
prior loosely rather than pin it. On the elasticity, a diffuse prior; the sign is the
restriction worth imposing.

**Identification caveat to check before committing.** λ̄ and σ_s both move σ(v) and cor(u,v)
and may be weakly separated by unconditional moments alone. The LD IRF is what should separate
them. Verify on simulated data before the sampler runs. If they are not separately identified,
λ̄ must be fixed on a reported grid, and the paper should say so plainly.

---

## 8. Recommendation

**1. Keep δ̄_e fixed at 3.2 %/yr for now, and describe it as a lower bound.** The GS critique
of the concept lands, but their measure is an upper bound and ours is the only directly
observed quantity. More importantly, raising δ̄_e does not deliver the moment it is supposed
to. Reopening this now would cost a full regeneration and buy +0.887 instead of +0.995.

**2. Reframe §3 of the paper's calibration discussion around bracketing, not choosing.** The
literature supplies an interval of roughly [3.2 %, 10 %], not a point. Say so, and say that
we take the observed lower endpoint deliberately because it is conservative for our result.

**3. Plan to move δ_e into Θ_e once the reposting margin exists**, with a prior on that
interval. The posterior on δ_e then reports the size of the middle margin, which is a result.

**4. Model reposting as a cost with dispersion, not as a constant λ.** Estimate the level and
the elasticity. The constant-λ version should appear, if at all, as the special case that
shows what the assumption was costing.

**5. Put the CK-consistency argument in the paper.** That three shocks are required rather
than chosen is the answer to the obvious referee question about why the antecedents do not
have this problem, and it converts an apparent weakness into the paper's first finding.

---

## 9. Costs and risks

**Proposition 5 Part 1 becomes conditional.** It currently proves an s shock raises u and v
together under costless reposting. With endogenous reposting it holds only above a threshold
λ̄*. This is arguably a stronger proposition, since it states *when* the Beveridge-inconsistent
response occurs and makes the threshold testable, but the proof must be redone. `lem:vpre` and
the Part 2 Jacobian both touch the reposting channel.

**The δ–s asymmetry becomes quantitative rather than qualitative.** With λ̄ < 1 both shocks
lower vacancies, and the distinction is degree and persistence rather than sign. This is the
most serious cost, because the sign asymmetry is the headline.

Two mitigations, both honest. The *variety* asymmetry is untouched: only δ moves N, so only δ
produces a persistent outward Beveridge shift, while a non-reposted s separation destroys a
position but leaves the product line and the entry margin intact. And the threshold λ̄* is
reportable: the paper can state where the sign flips and where the data place λ̄ relative to
it, which is more informative than asserting the sign.

**Scope.** This touches the v law of motion, Prop. 5, Θ_e, and the estimation design. It is
not a calibration tweak, and it should be settled before Block B is built, because it changes
what Block B targets.

---

## 10. Open items this note creates

| # | Item |
|---|---|
| 1 | Decide whether the reposting margin enters this paper. Gates Block B alongside [D10](../context/decisions.md) |
| 2 | If yes: open a decision entry, and revisit [D4](../context/decisions.md)/S6, which retire the LD→vacancy IRF the margin needs for identification |
| 3 | Verify λ̄ and σ_s are separately identified on simulated data |
| 4 | Add the GS, Shao-Silos, and CK-consistency discussion to §5.2 |
| 5 | Derive λ̄* and redo Prop. 5 Part 1 |
| 6 | Decide whether δ_e moves into Θ_e in this paper or is deferred |

## Sources

All in `Key papers/`, with a greppable text mirror under `Key papers/markdown/`
(see [`pipeline.md`](../context/pipeline.md) Part 3).

- Bilbiie, Ghironi & Melitz (2012), *JPE* 120(2) — calibration section and the introduction's
  product-destruction accounting
- Gabrovski & Silva (JEDC) — destruction-rate calibration; single-shock design
- Shao & Silos (2013), *EER* 63, 243–255, DOI 10.1016/j.euroecorev.2013.07.009 — §2.3
- Jaimovich & Floetotto (2008), *JME* 55(7), 1238–1252, DOI 10.1016/j.jmoneco.2008.08.008
- Bernard, Redding & Schott (2010), *AER* 100(1), 70–97, DOI 10.1257/aer.100.1.70 — cited here
  via BGM's reporting; the paper itself is not in `Key papers/`
- Coles & Kelishomi (2018) — δ_e ≡ τ, and the separation-shock design
- Model-side numbers: [`findings.md`](../context/findings.md) §D1/M5 with calibrated s
