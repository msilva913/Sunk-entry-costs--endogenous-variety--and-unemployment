# Calibrating δ: product lines, positions, and the reposting margin
**Written:** September 23, 2026 · **Revised:** September 28, 2026 · **Branch:** `costly_vacancy_reposting`
**Status:** analysis and recommendation. Nothing here is settled policy.
**Sept 28 revision:** Added §4.3a verifying that partial reposting is inescapable: the
reposting inflow is an accounting identity (Prop. 5 Part 1, Channel 1), so no parameter escapes
it, and every other lever is closed. Rewrote §7: the reposting cost is drawn from F, the
continuation cost distribution, scaled by α, as an independent **per-position** draw (each
vacated position draws its own χ ~ F). By the law of large numbers every surviving firm reposts
the same fraction Λ_r = F(Q/α), which preserves the single firm size the DS-CES aggregation
needs — a one-draw-per-firm version breaks it. The reposting rate is written **Λ_r** (not λ,
which is marginal utility). Supersedes both the separate-distribution G framing and the
same-draw α·χ framing. Added §7.3b settling the reactivation threshold: it is the stock value
the stock value Q (one-time reactivation, permanent retirement of non-reactivated slots), not the flow
value K (a recurring active-search-maintenance model — considered and rejected). Aligned with
[`../context/model_equations.md`](../context/model_equations.md) Block 2. See §4.3a, §7.3, §7.3a,
§7.3b.

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

**Raising δ̄_e helps, monotonically, but nowhere near enough.** cor(u,v) improves steadily,
from +0.995 to +0.887 as δ̄_e/τ̄ goes from 0.087 to 0.271. The direction is right and the
mechanism is real: a larger δ_e means more separations destroy vacancies. But the improvement
over the entire empirically defensible range is 0.108, against a gap of about 1.8 to the
data. Raising δ̄_e is a partial remedy, not a solution, and even the most aggressive defensible
value leaves the correlation strongly positive.

**The relationship is also strongly convex, and the measured range sits in its flat part.**
This follows from the model's own accounting. At δ_e = τ the identity s = (τ−δ_e)/(1−δ_e)
gives s̄ = 0: the s shock has nothing to scale, and cor(u,v) must equal the s-silenced value
of −0.911. So the full path runs from +0.995 at δ̄_e/τ̄ = 0.087 to −0.911 at δ̄_e/τ̄ = 1. The
average slope over that path is about −2.1 per unit of δ̄_e/τ̄, while the local slope in the
measured region is about −0.59, roughly 28 % of it. **Almost all of the Beveridge improvement
available from δ_e is concentrated in a region the measurement rules out.** That is why
δ_e cannot be the instrument of the fix, even though it pushes the right way.

**Silencing the s shock fixes it completely**, at −0.911. The Beveridge mechanics are sound.
What breaks them is the s shock under costless reposting.

**δ_e and Λ̄_r are complements, not substitutes, and this creates an identification concern.**
Both work through the same channel: the share of separations that destroy a vacancy rather
than recycle it. Raising δ_e moves separations into the destroying margin; lowering Λ̄_r makes
the non-destroying margin destroy. Because they act on the same moment through the same
channel, cor(u,v) alone will not separate them. This extends the Λ̄_r-versus-σ_s identification
concern in §7.5 to a three-way problem, and it is a further argument for holding δ_e fixed
while Λ̄_r is estimated, rather than freeing both at once.

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

### 4.3a Why no parameter escapes the s-reposting identity

The previous subsection asserts that costless reposting is untenable. This subsection verifies
it is inescapable, not merely inconvenient. The claim is that at reasonable δ, with the s
process the data delivers, no choice of parameters fits the Beveridge curve while separations
are costlessly reposted. Some separations must fail to repost.

**The reposting inflow is an accounting identity, not a parameter-dependent result.** Under
costless reposting, the response of pre-committed vacancies to an s shock is (Prop. 5 Part 1,
Channel 1):

```
∂v_pre,t+1 / ∂s_t = (1 − δ_{e,t+1})(1 − u_t) > 0
```

Every separation that is costlessly reposted returns next period as a vacancy. Proposition 5
Part 1 proves u↑ and v↑ together for **any** ρ_s ∈ [0,1), and the accompanying footnote shows
Channel 1 survives even under **free entry**, where the entry channel vanishes. So there is no
region of (ξ_inv, ε, b/w, entry elasticity, …) in which the s shock lowers v while separations
repost. The positive s contribution to cor(u,v) is structural, not calibrational.

**Every other lever is closed:**

- **Silence or shrink s.** Ruled out. σ_s = 0.0854 is data-pinned and s is strongly cyclical,
  cor(cycle_s, cycle_τ) = 0.997. And s carries most of the business-cycle volatility, so
  shrinking it to fix the sign destroys σ(u) and σ(v). This is the core tension: s must be
  large for volatility and is the same thing breaking the sign. Moving σ_s cannot resolve it.
- **Raise δ_e.** Ruled out quantitatively (§4.2). The whole defensible range buys 0.108 of a
  gap of about 1.8, and Part 3 already relies on δ_e/τ small to get the *right* sign. δ works;
  it is far too small to outweigh s.
- **Amplify through b/w.** Scales all shocks, including s's positive contribution. It does not
  touch the sign.
- **Entry.** Channel 2 can only add to ∂v/∂s. Removing it, at free entry, still leaves
  Channel 1 positive.

The only quantity that multiplies Channel 1 is the reposting rate. Attenuating it, by letting
some separations fail to repost, is the sole lever that reaches the culprit. **That is partial
reposting by definition.**

**What is and is not inescapable.** Two things must be separated:

1. **Λ_r < 1 in some form is inescapable** for the Beveridge sign, given reasonable δ and a
   data-calibrated s. This is what the identity forces.
2. **The endogenous per-position machinery is not required to fit cor(u,v) alone.** A constant
   Λ_r < 1 would hit that single moment. The endogenous version earns its extra structure
   elsewhere: matching the conditional LD→v IRF, surviving the Lucas critique (§7.1), and
   delivering the countercyclical position-destruction prediction.

**The genuine alternatives are larger changes, not smaller.** A fourth shock tuned to negative
comovement is a bigger conceptual change with weaker motivation, and it papers over the s
puzzle rather than resolving it. A wage-block overhaul to make z dominate is untested and may
not flip the sign, since real rigidity amplifies s's contribution too. Reinterpreting measured
s as partly position destruction is partial reposting under another name. None is smaller or
better-motivated.

**The required intervention is partial, not a sledgehammer.** With s silenced, z and δ alone
give cor(u,v) = −0.911, slightly *more* negative than the data's −0.804. So the data wants s
to contribute a small *positive* amount, not zero. Only a moderate interior Λ_r is needed to
attenuate the reposting inflow to that level, which is exactly what an interior α delivers.

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

### 7.1 Why a reduced-form Λ_r will not survive review

The natural first pass is an exogenous constant Λ_r ∈ [0,1], the fraction of match separations
after which the firm reposts, so the v law of motion's inflow becomes Λ_r·s_t·(1−u_t). Λ_r = 1
nests the current model.

A skeptical referee will say this should be written as a cost. **The referee is right, and
the objection is substantive rather than stylistic.**

- **It is inconsistent with the model's own logic.** Entry, exit, and vacancy creation are all
  optimizing margins. Reposting would be the only mechanical one, and it sits in the middle of
  the mechanism the paper is about.
- **It is not invariant to the shocks being studied.** This is the decisive objection. If
  reposting is a decision, firms repost less when the value of a position falls, which is
  in recessions. A constant Λ_r holds the reposting rate fixed exactly in the state where the
  Beveridge curve shifts. The parameter would not be structural in the Lucas sense, and the
  paper would be estimating a reduced form while claiming a mechanism.
- **It gives up the most interesting implication.** Endogenous reposting makes position
  destruction countercyclical for free, which is the empirically relevant pattern and a
  genuine prediction rather than an assumption.

### 7.2 A homogeneous fee does not work either

Worth stating explicitly, because it is the referee's likely first suggestion. Suppose every
repost costs the same fee c_r. Then the firm reposts if and only if Q_t ≥ c_r, and since all
positions are identical the decision is all-or-nothing: Λ_r ∈ {0, 1}. A homogeneous fee does
not generate partial reposting; it simply scales the value of a vacancy.

**An interior reposting rate requires dispersion.** A homogeneous fee compared against a
homogeneous value gives a corner. For a smooth interior Λ_r, either the reposting cost or the
value it is compared against must vary across positions. This is the same reason the
firm-level exit margin needs the distribution F to deliver an interior cutoff χ^c. The next
subsection shows which source of dispersion is available here.

### 7.3 The formulation: per-position reposting cost drawn from F

*This supersedes two earlier drafts. The first drew the reposting cost from a separate
distribution G tied to the entry cost; the second used one draw per firm, α·χ, which breaks
firm symmetry (below). The resolution draws the cost from F at the **position** level.*

**Notation.** The reposting rate is written **Λ_r**, capital, parallel to the survival
probability Λ = F(χ^c). It is not λ: lowercase λ is the marginal utility of consumption, which
appears in every Euler equation as the SDF ratio λ'/λ. Λ̄_r denotes the steady-state reposting
rate. See [`../context/model_equations.md`](../context/model_equations.md).

**Value dispersion is not available as the source.** The continuation cost χ is an iid
per-period draw, following the Broer (2025) formulation the exit margin already uses. Every
firm that continues this period faces the same distribution next period, so all continuing
firms value an additional filled position identically. A homogeneous reposting fee then gives
bang-bang reposting, Λ_r ∈ {0,1}, not a smooth rate. There is no forward-looking firm-value
dispersion for the reposting decision to exploit. The dispersion must come from the cost side.

**The cost should not be tied to G.** The sunk cost x drawn from G pays for office space and
HR infrastructure. That asset is durable and position-specific. It is not un-built when a
match dissolves. So reposting a vacated slot does not re-incur x and does not scale with the
original draw. G is the cost of creating capacity. Reposting uses capacity that already
exists. The two are different objects, and tying reposting to G has no economic content.

**The cost is drawn from F, per position, scaled by α.** On separation, each vacated position
independently draws a reactivation cost from F, scaled by α. This is a one-time decision:
reactivate the slot for search (obtain a durable vacancy) or retire it permanently. Reactivate
iff the value of the vacancy obtained covers the cost:

```
cost of reactivating a vacated position  =  α · χ ,   χ ~ F,   α > 0
reactivate  ⟺  Q_t ≥ α · χ  ⟺  χ ≤ Q_t/α
```

The threshold is the value of the vacancy obtained, the **gross/stock** value Q, not the flow
value K. The R7 derivation (`model_equations.md` f[3], against the draft's recruiter block
`eq:value_recruiter_matched`) confirms this exactly: on separation the recruiter obtains an
unfilled vacancy worth Q_{t+1}, which the Bellman already discounts, so reactivation compares α·χ
to Q directly. (An earlier draft wrote the threshold as Q − K; that was an over-refinement — the
exact object is Q, and K ≈ 0.6% of Q made the difference negligible anyway.) Whether the
threshold is the stock value Q or the flow value K is the substance of the
one-time-versus-recurring choice settled in §7.3b. This adds one parameter, the scale α, and no
new distribution.

**Per-position draws, not one per firm — this preserves the single firm size.** Broer's iid χ
keeps all firms identical at production because exit is *terminal*: survivors reset each period
and the draw leaves no trace on anyone who stays. Reposting instead lands on a *continuing*
state, the firm's stock of positions. One draw per firm would make each firm repost all or none
of its separations, so firms would carry different stocks forward, firm-size heterogeneity would
accumulate, and the DS-CES symmetric aggregation (ρ = N^(1/(ε−1))) would break. Drawing per
position fixes this: with each firm holding many positions, the law of large numbers makes every
surviving firm repost the identical fraction

```
Λ_r = F(Q_t/α) = min[ (1 − p_0) + p_0·(Q_t/(α·f_m))^ψ , 1 ]
```

so all firms shrink together and stay symmetric replicas of the aggregate. One firm size
survives, exactly as Broer intended, now at the position level.

**The vacancy law of motion** becomes the baseline with Λ_r inserted on the separation inflow:

```
v_pre,t+1 = (1 − δ_{e,t})·[ (1 − q_t)·v_t + Λ_r·s_t·sbar·(1 − u_t) ]
```

Only δ_e appears, because reposting is a sub-event of survival: a position reposts only if its
firm continues, probability 1 − δ_e, so (1 − δ_e) factors out and Λ_r sits inside. The
non-reposted fraction 1 − Λ_r of separations is the missing middle of §5 — position destruction
at surviving firms — now spread uniformly across firms rather than concentrated. See
[`../context/model_equations.md`](../context/model_equations.md) Block 2 for the full equation
set.

**Properties.** Λ_r is procyclical: Q rises in booms, so Q/α rises, more positions repost, and
position destruction is countercyclical. Bounds: Λ_r ∈ [1 − p_0, 1]. The atom 1 − p_0 always
reposts (zero-cost positions), so the most that can ever be destroyed is fraction p_0 of
separations — **p_0 caps the margin's strength.** Full costless reposting is nested at
α ≤ Q̄/f_m, where Q/α ≥ f_m ⟹ F = 1 ⟹ Λ_r = 1; the margin bites (Λ_r < 1) only for α > Q̄/f_m.

**Identification.** α is what the cor(u,v) and LD→v moments identify, and the estimate reports
how much reposting friction the Beveridge curve demands. The durable-asset logic says reposting
should be cheap; the Beveridge curve needs it costly enough to bite (α > Q̄/f_m). The posterior
on α tells us whether these are reconcilable.

### 7.3a Why the reposting and continuation costs share F

Drawing both costs from the same distribution needs an economic rationale, not just a plea to
parsimony. The rationale is a two-technology split.

- **G is the distribution of costs to _create_ productive capacity.** Office space, capital,
  HR infrastructure. One-time, sunk, durable.
- **F is the distribution of idiosyncratic frictions in _operating and maintaining_ capacity
  that already exists.** Recurring, drawn anew each period.

Continuation and reposting are both operate-existing-capacity decisions. Continuation asks
whether it is worth keeping the product line running given this period's friction draw.
Reposting asks whether it is worth restoring a vacated position to operation given a draw from
the same law. Entry is a create-capacity decision, so it draws from G. Reposting is not, so it
draws from F.

The friction F represents is the firm's period-by-period idiosyncratic cost of devoting
organizational capacity to keeping units in operation. Managerial attention, compliance and
certification, internal coordination, access to congested service and hiring markets. In a
disrupted firm-period all of this is expensive. In a slack period it is free. The mass point
1 − p_0 has one interpretation on both margins: periods with no binding disruption, where
continuing and re-staffing are both costless. The power-law tail is the severity of
disruption when it bites.

The paper already commits to F for the exit margin, following Broer (2025). Once idiosyncratic
fixed costs of adjustment are drawn from a common law each period, the consistent default is
that the same law governs the position margin. Positing a second, unrelated cost distribution
for reposting is the move that needs a positive justification. There is no evidence that the
severity law of re-staffing frictions differs in shape from that of continuation frictions.

**Same shape, different scale, independent draws.** The two margins differ in scale, not shape.
Maintaining a product line is a larger object than restoring one position, and α captures that
wedge. The family and the tail index ψ are shared; the level is not. The draws are independent:
exit is one firm-level draw χ compared to χ^c; reposting is many position-level draws from the
same family compared to Q/α. Independence across positions is what makes reposting uniform
across firms and preserves the single firm size. Position-level idiosyncratic reactivation
frictions are, if anything, the more natural reading of a per-slot friction than one firm-wide
draw.

This is a discipline, not a free assumption. Sharing F's shape forces the exit rate and the
aggregate position-destruction rate to comove. Both rise with idiosyncratic-friction severity.
Both respond to aggregate conditions, but through their own value objects, χ^c for exit and Q/α
for reposting. The data can reject a common ψ.

### 7.3b One-time reactivation versus recurring maintenance: which reposting story?

The reactivation threshold is the value of the vacancy obtained. There are two economically
distinct ways to model what reposting *is*, and they imply different thresholds — the stock
value Q or the flow value K. Each has a plausible story.

**Story A — one-time reactivation (threshold ≈ Q).** When a worker leaves, the firm makes a
discrete, lumpy decision about that specific slot: reactivate it for the search market, or retire
it. The physical capacity is durable and already built (office, equipment, HR infrastructure from
the sunk cost x ~ G), but bringing the slot back to hiring-ready status is a one-time job:
rewrite and re-approve the requisition, repost the ad, re-brief recruiters, reconfigure the
workstation. That cost is α·χ, idiosyncratic because it depends on the firm's organizational
state at the moment of separation. A reactivated slot becomes a durable vacancy that searches
until it fills, at no further reposting cost. A retired slot is gone: the headcount line is cut,
and to have it back later the firm must create a fresh position at the full cost Q. Reactivate
iff Q ≥ α·χ. Empirical counterpart: requisition and headcount decisions — firms freeze
or eliminate reqs after departures in bad times.

**Story B — recurring active-search maintenance (threshold ≈ K).** The slot persists costlessly
once built, but being *actively in the hiring market* is a recurring burden: keeping the req live
and refreshed, paying recruiters to screen and interview this month, holding the budget line
open, ongoing compliance. Each period the firm pays α·χ (fresh draw) to keep the position
actively recruiting, or lets it go idle — dormant, not matching this period, but not destroyed.
Idling forfeits only this period's matching shot, worth K ≈ (r+δ_e)·Q, while the slot survives.
So the firm searches iff K ≥ α·χ. Empirical counterpart: firms throttle recruiting *intensity* —
reqs go "on hold" in bad months and reactivate later — without the position being eliminated.

**Verdict: Story A.** Both are endogenous and procyclical, so both clear the §7.1 Lucas-critique
bar. Story A wins on four dimensions that matter for this paper:

1. **Concept.** The mechanism is the missing middle of §5 — *permanent* position destruction at
   surviving firms. That is Story A. Story B is temporary withdrawal (a recruiting pause), which
   quietly redefines the mechanism away from destruction.
2. **Persistence.** The reposting channel's main dynamic prize is persistence via stock
   depletion: destroyed positions must be rebuilt at the full cost Q, slowing recovery and
   bearing on the δ→u peak-horizon gap (M8). A delivers it; B undercuts it, because idle
   positions reactivate cheaply and the withdrawal is transitory.
3. **Parsimony and fit.** A stays inside the model's "one-time posting cost Q + per-match cost κ"
   structure and adds no state. B introduces a per-period flow vacancy cost the model
   deliberately lacks, *and* a dormant-position state, and risks double-counting against κ.
4. **Empirics.** The persistent LD→v trough (negative through h≈20) that identifies the margin
   (§7.4) matches permanent destruction (A); temporary withdrawal (B) predicts a shallower,
   faster-reverting trough.

B has the more vivid micro-story — recruiting intensity is genuinely recurring — so it is worth a
sentence as an alternative interpretation, but its temporariness, extra machinery, and weaker
persistence make it the wrong choice here. The model uses Story A: threshold Q (the stock value;
R7 confirms it is exactly Q, not Q − K), one-time reactivation at separation, permanent
retirement of non-reactivated positions.

### 7.4 Identification

The margin adds one free parameter, the scale α. The level is Λ̄_r = F(Q̄/α). The
responsiveness of Λ_r to Q_t is not a separate parameter. It is the slope of F at Q̄/α, and
F's shape (p_0, ψ, f_m) is already identified by the exit margin. So the reposting margin asks
the data for one number, α, plus a consistency check that F's shape fits both margins.

1. **The LD→vacancy IRF identifies α.** In the model an s shock raises v when Λ̄_r is high and
   lowers it when Λ̄_r is low, with a threshold Λ_r* where the sign flips. The empirical response
   is negative and persistent, placing Λ̄_r below the threshold. Its magnitude pins the level,
   hence α, and its *shape* over horizons disciplines whether the shared F is consistent,
   because reposting deepens the response as Q falls. This gives the LD arm of Block B a
   structural job for the first time and repays [D4](../context/decisions.md)/S6.
2. **cor(u,v) in Block M** moves monotonically in Λ̄_r. Powerful, but it is the moment being
   explained, so it should not be the sole source.
3. **BED gross job flows** bound Λ̄_r loosely. Gross job losses split into losses at closing
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

**2. It adjudicates a headline result.** Λ̄_r determines whether the model reproduces the
Beveridge curve. Fixing it externally would assume the answer to the question being asked.
Note this is the mirror of the D1 logic and points the other way: δ_e is fixed *because* it is
cleanly measured and freeing it would let it absorb misspecification. Reposting *is* the
misspecification, so it is the right thing to free.

**3. The responsiveness is inherited from F, not free.** Under the per-position formulation the
responsiveness of reposting to the value of a position is the slope of F at Q̄/α, and F's shape
is already pinned by the exit margin. So the only new object to estimate is the scale α. It
must be estimated for the same two reasons: no clean data object measures it, and it
adjudicates the Beveridge result.

**Prior.** On α, a prior with support in (Q̄/f_m, ∞), the range that makes the margin bite
(Λ_r < 1); the plausible region also keeps reposting cheaper than continuation (α < 1). The
recall-rate literature bounds non-reposting from above and should shape the prior loosely rather
than pin it. F's shape parameters keep their existing priors from the exit margin. The reposting
margin adds no new shape prior.

**Identification caveat to check before committing.** Two checks. First, α and σ_s both move
σ(v) and cor(u,v) and may be weakly separated by unconditional moments alone; the LD IRF is
what should separate them. Second, α and F's scale f_m both affect the position-destruction
rate, and they must be separated by the fact that exit compares χ to χ^c while reposting
compares it to Q, so the two cutoffs move differently over the cycle. Verify both on simulated
data before the sampler runs. If α is not separately identified, fix it on a reported grid and
say so plainly.

---

## 8. Recommendation

**1. Keep δ̄_e fixed at 3.2 %/yr for now, and describe it as a lower bound.** The GS critique
of the concept lands, but their measure is an upper bound and ours is the only directly
observed quantity. Raising δ̄_e does push cor(u,v) the right way, and that should be said
rather than denied — but the whole defensible range buys 0.108 of a gap of about 1.8, and
reopening this now would cost a full regeneration of every mechanism figure to move from
+0.995 to +0.887. The return does not justify it, and the estimation is the right place to
revisit the level once the reposting margin exists.

**2. Reframe §3 of the paper's calibration discussion around bracketing, not choosing.** The
literature supplies an interval of roughly [3.2 %, 10 %], not a point. Say so, and say that
we take the observed lower endpoint deliberately because it is conservative for our result.

**3. Plan to move δ_e into Θ_e once the reposting margin exists**, with a prior on that
interval. The posterior on δ_e then reports the size of the middle margin, which is a result.

**4. Model reposting as a cost drawn from F at scale α, not as a constant Λ_r.** Each vacated
position draws a reposting cost α·χ, χ ~ F, and reposts if the vacancy value Q covers it, giving
Λ_r = F(Q/α). This reuses the exit margin's distribution, adds one parameter, preserves the
single firm size (per-position draws), and generates the missing middle of §5 as a by-product.
Estimate α. The constant-Λ_r version should appear, if at all, as the special case that shows
what the assumption was costing.

**5. Put the CK-consistency argument in the paper.** That three shocks are required rather
than chosen is the answer to the obvious referee question about why the antecedents do not
have this problem, and it converts an apparent weakness into the paper's first finding.

---

## 9. Costs and risks

**Proposition 5 Part 1 becomes conditional.** It currently proves an s shock raises u and v
together under costless reposting. With endogenous reposting it holds only above a threshold
Λ_r*. This is arguably a stronger proposition, since it states *when* the Beveridge-inconsistent
response occurs and makes the threshold testable, but the proof must be redone. `lem:vpre` and
the Part 2 Jacobian both touch the reposting channel.

**The δ–s asymmetry becomes quantitative rather than qualitative.** With Λ̄_r < 1 both shocks
lower vacancies, and the distinction is degree and persistence rather than sign. This is the
most serious cost, because the sign asymmetry is the headline.

Two mitigations, both honest. The *variety* asymmetry is untouched: only δ moves N, so only δ
produces a persistent outward Beveridge shift, while a non-reposted s separation destroys a
position but leaves the product line and the entry margin intact. And the threshold Λ_r* is
reportable: the paper can state where the sign flips and where the data place Λ̄_r relative to
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
| 3 | Verify on simulated data that α is separately identified from both σ_s and F's scale f_m (see §7.5) |
| 4 | Add the GS, Shao-Silos, and CK-consistency discussion to §5.2 |
| 5 | Derive Λ_r* and redo Prop. 5 Part 1 |
| 6 | Decide whether δ_e moves into Θ_e in this paper or is deferred |
| 7 | Confirm the per-position formulation against the firm Bellman: firm-level χ governs continuation (cutoff χ^c); position-level draws from F govern reposting (cutoff Q/α); check the calibrated point sits in the biting range α > Q̄/f_m |
| 8 | Re-derive the job creation condition with the reposting option (see [`../context/model_equations.md`](../context/model_equations.md) Block 2, f[3]) |

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
- Broer (2025), *The Unemployment Risk Channel in Business Cycle Fluctuations* — source of the
  iid per-period fixed-cost formulation used for the exit margin, and reused here for reposting.
  Citation unverified; full author list, venue, and DOI to be added
- Model-side numbers: [`findings.md`](../context/findings.md) §D1/M5 with calibrated s
