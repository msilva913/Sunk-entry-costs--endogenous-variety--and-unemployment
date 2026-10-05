# Can entry outrun destruction? The predetermined vacancy share and Coles–Kelishomi as a special case

**Written:** October 5, 2026
**Numbers from:** `../Programs baseline/predetermined_vacancy_share.jl` →
`predetermined_vacancy_share.txt` (self-contained; mirrors the closed-form steady-state
algebra of `steady_state.jl`, with every identity checked by assertion)
**Draft objects used:** `eq:v_lom`, `eq:v_ll`, `eq:e_ss`, `app:comparison_D`, `app:loglin` ¶7

---

## 1. The question

A positive δ shock does two things to vacancies at once. It destroys the vacancies carried
into the period, and it changes the incentive to post new ones. For the Beveridge curve what
matters is not whether entry responds but whether it responds **enough** to cover the loss.

Coles and Kelishomi need finitely elastic vacancy creation to stop vacancies rising with
unemployment after a destruction shock. We claim not to. The question this note settles is
*why*, since the obvious argument cuts both ways: a higher destruction rate raises the entry
incentive but also destroys more vacancies, so it is not obvious which side wins.

The answer is that the destruction term is proportional to the **predetermined** share of the
vacancy stock, not to the whole stock, and that share collapses as δ_e rises. At CK's
calibration there is almost no inherited stock left for a destruction shock to remove.

---

## 2. The accounting

The vacancy law of motion (`eq:v_lom`) is

```
v_t = (1 − δ_e,t)[(1 − q(θ_{t−1})) v_{t−1} + Λ_r,t s_{t−1}(1 − u_{t−1})] + e_t
```

Log-linearized (`eq:v_ll`, in **level** deviations x̃ = x − x̄):

```
ṽ_t = (1 − δ̄_e)[(1 − q̄) ṽ_{t−1} + s̄(1 − ū) s̃_{t−1}] − v̄^pre · δ̃_e,t + ẽ_t
```

with **v̄^pre ≡ v̄ − ē**, the steady-state predetermined stock.

On impact of a pure δ shock, with ṽ_{t−1} = 0 and s̃ = 0, everything but two terms drops:

```
ṽ_t = ẽ_t − v̄^pre · δ̃_e,t
```

So **vacancies rise if and only if new postings exceed the destroyed predetermined stock**:

```
ẽ_t > v̄^pre · δ̃_e,t
```

This is the whole question, and it is exact. Note what the destruction term is proportional
to: `v̄^pre`, not `v̄`. A destruction shock can only remove vacancies that were already
standing at the start of the period. Vacancies posted within the period are not exposed.

---

## 3. The threshold, in three equivalent forms

Converting to log deviations (ê = ẽ/ē, δ̂_e = δ̃_e/δ̄_e) gives the threshold entry elasticity

```
ê_t / δ̂_e,t  >  T  ≡  v̄^pre / ē  =  v̄/ē − 1
```

Using the steady-state entry relation ē = δ̄_e(v̄ + L̄) (`eq:e_ss`, `steady_state.jl:217`):

```
T = [v̄/(v̄ + L̄)] / δ̄_e − 1                                    (form 2)
```

and since the f, q and sep targets pin both v̄/(v̄+L̄) and τ̄, writing
κ_v ≡ [v̄/(v̄+L̄)]/τ̄:

```
T = κ_v / (δ̄_e/τ̄) − 1 ,        κ_v = 1.2034                   (form 3)
```

**The threshold is inversely proportional to the destruction rate, and given the labor-market
targets it is a function of δ̄_e/τ̄ alone.** The same algebra gives the two statistics
directly:

```
ē/v̄      = (δ̄_e/τ̄) / κ_v          entrant share of the vacancy stock
v̄^pre/v̄  = 1 − (δ̄_e/τ̄) / κ_v      predetermined share
```

So the predetermined share and δ̄_e/τ̄ are the same object up to the constant κ_v. There is
only one free statistic here, not two.

### Why the obvious scaling arguments fail

Two tempting arguments both cancel and must be set aside.

- *"A δ shock destroys δ̄_e × v̄ vacancies, so the hole scales with δ̄_e; but entry also scales
  with δ̄_e, so they cancel."* The first clause is wrong. The hole is v̄^pre × δ̃_e, and
  v̄^pre itself depends on δ̄_e. That dependence is the entire result.
- *"Higher δ̄_e raises unemployment more, which slackens the market more, which draws in more
  entry."* True, but it also raises the replacement flow and the destroyed flow. The
  incentive channel does not break the tie; the composition channel does.

---

## 4. The numbers

Targets held fixed across all rows: f = 0.41, q = 0.80, sep = 0.031 (monthly). θ, u, L and v
are **invariant** to δ̄_e by construction, because s adjusts to absorb the change — this is
the design of Comparison A ("all other targets unchanged → same SS u, v, θ").

```
θ = 0.5125   u = 0.070295   L = 0.929705   v = 0.036026
v̄/(v̄+L̄) = 0.037304        κ_v = 1.203369
```

| Calibration | δ̄_e (monthly) | δ̄_e/τ̄ | ē/v̄ | **v̄^pre/v̄** | **T** |
|---|---|---|---|---|---|
| **D1 settled** (BED Deaths, 3.20 %/yr) | 0.002707 | 0.0873 | 0.073 | **0.927** | **12.78** |
| Code baseline (superseded, 7.54 %/yr) | 0.006512 | 0.2100 | 0.175 | 0.825 | 4.73 |
| BGM / GS / Shao–Silos (10 %/yr) | 0.008742 | 0.2820 | 0.234 | 0.766 | 3.27 |
| **Coles–Kelishomi** (δ_e = τ) | 0.031000 | 1.0000 | 0.831 | **0.169** | **0.203** |

**The threshold is 62.9 times higher at the D1 calibration than at CK's.** At CK's, 83 % of
the standing vacancy stock is new postings and only 17 % is inherited, so a destruction shock
has almost nothing to remove and a 0.2 % entry response suffices to turn the vacancy
response positive. At D1's value 93 % of the stock is inherited and entry would have to rise
by nearly 13 % per 1 % destruction shock.

---

## 5. Coles–Kelishomi as a special case

CK set δ_e = τ: every separation destroys the position. Two things follow.

**It is the maximum admissible value.** s ≥ 0 requires δ_e ≤ τ, so δ_e = τ is the boundary of
the parameter space. By form 3, T is decreasing in δ̄_e/τ̄, so **CK sit exactly at the point
that minimizes the entry threshold over the entire admissible range.** There is no
calibration of this model in which entry has an easier time outrunning destruction.

**It is why they need a friction.** With T = 0.203, almost any positive entry response makes
vacancies rise with unemployment. Making vacancy creation finitely elastic is what holds the
entry response below the threshold. So the assumption is not doing work about the labor
market in general; it is repairing a consequence of having one separation margin.

Conversely, our result does not come from the business formation block, from variety, or from
anything on the preference side. It comes from measuring δ_e and finding it an order of
magnitude below τ, which raises the threshold by roughly the same factor. **Gabrovski–Silva
already have the δ/s split, so an augmented GS model gets this too** — the contribution here
is the measurement, not the structure. (The structure earns its place on other grounds; see
`role_of_ingredients.md` and the intro's necessity argument.)

---

## 6. What this does and does not establish

**Establishes.** The bar entry must clear, exactly, as a function of one calibrated
statistic. Why that bar moves by two orders of magnitude between CK's calibration and ours.
Why CK's calibration is the most favorable possible case for positive comovement.

**Does not establish.** Whether the model's entry response actually clears the bar. That
needs ê itself, which comes from ê = ξQ̂ together with Q̂ from the asset-pricing relation and
the JCC (`app:loglin` ¶9, ¶9a). An indicative check: at the baseline ξ_inv = 1.0, so ξ = 1
and ê = Q̂. Clearing T = 12.78 would need Q̂ ≈ 12.8 % for a 1 % δ shock, which is far outside
anything the IRFs show; at CK's T = 0.203 it needs only Q̂ ≈ 0.2 %. That is consistent with
the Comparison D figures, where vacancies fall at every ξ we run. ⚠️ Indicative only — Q̂ is
not computed in this note.

**Relation to the draft's existing argument.** `app:comparison_D` argues through v = θu and
the free-entry limit, which is a different decomposition and reaches the same conclusion. It
also contains a non sequitur: it says that at low δ̄_e "the proportional u rise is small ...
*and the duration shortening effect on K is correspondingly weak*. The θ decline therefore
dominates." A weaker effect on K means θ falls *less*, which makes vacancies *more* likely to
rise. The conclusion is right but that sentence argues against it. The accounting in this
note is cleaner because the destruction term appears explicitly instead of being hidden
inside θ.

**A notation bug to fix.** `eq:u_ll` writes the destruction impulse as δ̃_e,t (level
deviation, coefficient ūf̄); `eq:v_ll` writes δ̂_e,t with coefficient v̄^pre. Since
∂/∂δ_e of (1−δ_e)·stock is −stock, the coefficient v̄^pre is correct for a **level**
deviation, so `eq:v_ll`'s hat should be a tilde. As printed, a reader who takes the hat as a
log deviation overstates the destruction term by a factor 1/δ̄_e ≈ 370.

---

## 7. For the draft

1. **The intro paragraph on CK** should run on the predetermined share, not on duration or on
   the relative size of the u and θ responses. The sentence to make is that a destruction
   shock can only remove vacancies that were already standing, and that at δ_e = τ almost
   none are.
2. **`app:loglin` ¶7** is the natural home for the threshold. It already introduces
   v̄^pre = v̄ − ē and says the destruction impulse "directly destroys the predetermined
   stock" without drawing the consequence.
3. **Fix the `app:comparison_D` non sequitur** and the `eq:v_ll` hat.
4. **Numbers are program-sourced** (N15 satisfied) via `predetermined_vacancy_share.jl`. The
   two worth quoting are v̄^pre/v̄ and T at D1 against CK.
5. ⚠️ The code still runs `dest_ann = 0.0754`, so the middle row is what the current figures
   embody. T = 4.73 there against 12.78 at D1: **the D1 cascade strengthens this argument**,
   and §5.3 understates it today.
