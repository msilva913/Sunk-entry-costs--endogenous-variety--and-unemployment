# The δ shock at free entry (ξ → ∞): a step-by-step derivation in the GS economy

**Written:** October 6, 2026. Not for the draft yet.
**Numbers from:** `../Programs baseline/gs_free_entry_beveridge.jl` →
`gs_free_entry_beveridge.txt` (last block, "Exposition tables"). Every closed form below is
checked there by assertion against a simulated IRF of the same linear system.
**Setting:** our model under the conditions of `prop:independence` (σ = 0, ζ = 0, p₀ = 0).
The labor-market block then closes on its own, and it is the Gabrovski–Silva economy with a
time-varying destruction rate δ_t.

**One simplification for exposition: κ = 0.** All recruiting costs are sunk (our x_v = 1;
GS's γ = 0; CK's c = 0). §7 says what κ > 0 changes. Numbers barely move.

---

## 1. The nonlinear system

### Timing (the draft's convention, S10)

- δ_t is realized at the start of t. Every period-t decision conditions on it.
- The destruction it causes happens at the end of t.
- So u_t and the inherited vacancies are fixed at t. Tightness θ_t and entry e_t respond
  at t. The destroyed jobs and positions show up at t+1.

### Matching

f(θ) = Aθ^(1−η) and q(θ) = Aθ^(−η), with θ = v/u.

### Equations

```
(N1) Free entry         Q_t = x_m
(N2) Vacancy rental     K_t = Q_t − β(1 − δ_t) Q_{t+1}
(N3) Job creation       K_t / q(θ_t) = β(1 − δ_t) E_t S_{t+1}
        with  S_{t+1} = (1−ϕ)(y − b − K_{t+1}) + (1 − s − ϕ f(θ_{t+1})) K_{t+1}/q(θ_{t+1})
(N4) Unemployment       u_{t+1} = u_t − (1 − δ_t) f(θ_t) u_t + [δ_t + s(1 − δ_t)](1 − u_t)
(N5) Vacancies          v_{t+1} = (1 − δ_t)[ (1 − q(θ_t)) v_t + s(1 − u_t) ] + e_{t+1}
(N6) Tightness          θ_t = v_t / u_t
```

### What each equation says

- **(N1) Free entry.** The cost of creating a vacancy is drawn from G on [0, x_m]. With
  infinitely elastic creation (ξ → ∞), the value of a vacancy is pinned at x_m whenever some
  entry occurs.
- **(N2) Rental.** K_t is the flow value of holding a vacancy for one period (eq:Kdef). The
  firm holds an asset worth Q_t. Next period the asset is worth Q_{t+1}, but only if the
  business opportunity survives the end-of-period destruction draw, with probability
  1 − δ_t.
- **(N3) Job creation** (eq:jcc_eq).
  - **Left side: the expected cost of a hire.** The rental K_t times the expected time to
    fill, 1/q(θ_t).
  - **Right side: the discounted value of a match.** It counts only if the job survives the
    destruction draw at the end of t, since a match formed at t first produces at t+1.
  - **S: the firm's surplus from a match**, after substituting the Nash wage. It has two
    parts:
    - (1−ϕ)(y − b − K): the firm's share of output net of the worker's outside option and
      of the forgone rental. A matched firm no longer holds an open vacancy.
    - (1 − s − ϕf)·K/q: the hiring cost the firm saves next period because the match
      continues. The term ϕf·K/q is the part of that saving the worker captures through the
      wage, via her outside option of finding another job.
- **(N4) Unemployment** (eq:u_lom). Next period's unemployed come from two groups:
  - today's unemployed who were not hired, or were hired at a business that was destroyed
    before producing;
  - the employed who lose their job, either to destruction δ_t or to separation s at a
    surviving business.
- **(N5) Vacancies** (eq:v_lom). Next period's vacancies come from three sources:
  - today's unfilled vacancies that survive destruction;
  - positions vacated by separations at surviving businesses and reposted at no cost;
  - new entry e_{t+1}.
- **(N6)** defines tightness.

### How free entry orders the system

- (N1)–(N3) contain no stock. Given the path of δ, they determine the path of θ by
  themselves.
- (N4) then gives u.
- (N6) gives v = θu.
- (N5) becomes the equation that determines entry e. Entry is whatever it takes to deliver
  v = θu, given the inherited vacancies.

**With finite ξ, the causation runs the other way.** e = G(Q) is set by Q. (N5) determines
v, (N6) gives θ, and (N3) gives Q. This reversal is why the inherited-vacancy accounting
decides nothing at ξ = ∞.

---

## 2. Three steady-state facts used below

Write r ≡ 1/β − 1 for the discount rate, and B ≡ β(1 − δ̄).

1. **Rental.** From (N1)–(N2), K̄ = x_m[1 − β(1−δ̄)] = x_m(r + δ̄)/(1 + r). This is a user
   cost: the asset price times the discount rate plus the depreciation rate.
2. **Job creation.** K̄/q̄ = B·S̄.
3. **Unemployment.** From (N4), ū = τ̄/(τ̄ + (1−δ̄)f̄), with τ̄ ≡ δ̄ + s̄(1−δ̄). Equivalently,
   (1 − ū)/ū = (1−δ̄)f̄/τ̄.

---

## 3. Log-linearization, equation by equation

Hats are log deviations from the steady state. For the destruction rate we use the level
deviation d_t ≡ δ_t − δ̄, since δ_t is a small rate.

### 3.1 Rental

```
(L1)   K̂_t = d_t / (r + δ̄)
```

**Derivation.** With Q fixed, dK_t = βx_m·d_t, and βx_m/K̄ = 1/(r + δ̄).

**Interpretation.** A one-point rise in the destruction rate raises the vacancy's rental by
1/(r + δ̄) in proportion. The longer-lived the vacancy (small r + δ̄), the larger the
proportional jump.

| | D1 | CK-like (δ = τ) |
|---|---|---|
| 1/(r + δ̄) | **167** | **29** |

Among the weights that matter for u and v at free entry, this is the only one that differs
much across the two calibrations. (The vacancy-equation weights in §3.4 differ too, but at
free entry they only determine entry.)

### 3.2 Job creation

**Derivation.**

- Left side: ln(K/q) moves with K̂_t + ηθ̂_t, since q ∝ θ^(−η).
- Right side: ln β(1−δ_t) moves with −d_t/(1−δ̄).
- For S, differentiate term by term:
  - dK = K̄·K̂;
  - df = f̄(1−η)θ̂;
  - d(K/q) = (K̄/q̄)(K̂ + ηθ̂).
- Divide by S̄ = K̄/(q̄B) (fact 2).

The result is

```
(L2)   K̂_t + ηθ̂_t  =  − d_t/(1 − δ̄)
                       + B·E_t{ [ (1 − s̄ − ϕf̄) − (1−ϕ)q̄ ] K̂_{t+1}
                              + [ (1 − s̄ − ϕf̄)η − ϕf̄(1−η) ] θ̂_{t+1} }
```

**Interpretation of each weight.**

- **η on θ̂_t (left).** Tighter markets lengthen the expected time to fill, so a hire
  costs more.
- **1/(1−δ̄) on d_t.** A hire made at t pays off only if the job survives destruction at
  the end of t. This is about 1.
- **(1 − s̄ − ϕf̄) on K̂_{t+1}.** A higher rental next period makes it more valuable to
  already have a worker, because the firm saves next period's hiring cost. It is weighted by
  the probability the match continues (1 − s̄), net of the share the worker captures (ϕf̄).
- **−(1−ϕ)q̄ on K̂_{t+1}.** A matched firm forgoes the vacancy's rental, so a higher
  rental lowers the firm's surplus.
- **(1 − s̄ − ϕf̄)η − ϕf̄(1−η) on θ̂_{t+1}.**
  - Tighter markets next period raise the hiring cost a continuing match saves (+).
  - They also raise the worker's outside option and hence the wage (−).

After (L1) is substituted, **(L2) involves only θ̂ and d.** Neither u nor v appears.

### 3.3 Unemployment

**Derivation.** Differentiate (N4) with respect to u_t, δ_t and ln θ_t:

- ∂u_{t+1}/∂u_t = 1 − (1−δ̄)f̄ − τ̄;
- ∂u_{t+1}/∂δ_t = f̄ū + (1−s̄)(1−ū);
- ∂u_{t+1}/∂ln θ_t = −(1−δ̄)f̄(1−η)ū.

Divide by ū. Fact 3 gives [f̄ū + (1−s̄)(1−ū)]/ū = f̄[τ̄ + (1−s̄)(1−δ̄)]/τ̄ = f̄/τ̄, because
(1−s̄)(1−δ̄) = 1 − τ̄. Define

```
λ ≡ 1 − (1−δ̄)f̄ − τ̄          ω ≡ (1−δ̄) f̄ (1−η)
```

Then

```
(L3)   û_{t+1} = λ û_t  +  (f̄/τ̄) d_t  −  ω θ̂_t
```

**Interpretation of each weight.**

- **λ: persistence.** An extra unemployed worker today is still unemployed next month
  unless hired, with probability 1 − (1−δ̄)f̄. The τ̄ adjustment reflects that one fewer
  employed worker means one fewer separation.
- **f̄/τ̄: the destruction weight.** A one-point rise in δ_t adds two groups to unemployment:
  - the employed who lose their job, (1−s̄)(1−ū);
  - the new matches lost before producing, f̄ū.

  Relative to the stock ū, that is f̄/τ̄. **The weight is large because unemployment is a
  small stock, ū ≈ τ̄/(τ̄ + f̄), so a given inflow is a large proportional change.** It
  depends on τ̄ and f̄, not on how τ̄ splits into δ and s.
- **ω: the hiring channel.** Lower tightness at t means fewer hires at t, and those workers
  are still unemployed at t+1.

| | D1 | CK-like |
|---|---|---|
| λ | 0.559 | 0.559 |
| f̄/τ̄ | 13.26 | 13.65 |
| ω | 0.164 | 0.164 |

**The unemployment equation is essentially the same in both economies**, because f̄, q̄
and τ̄ are held fixed.

### 3.4 Vacancies

**Derivation.** Differentiate (N5). The only non-obvious term is the fill rate:
dq = −ηq̄θ̂, so −(1−δ̄)v̄·dq = (1−δ̄)ηq̄v̄θ̂. Divide by v̄, use ū/v̄ = 1/θ̄, and write
X̄ ≡ (1−q̄)v̄ + s̄(1−ū) for the stock carried into the next period before destruction:

```
(L4)   v̂_{t+1} = − (X̄/v̄) d_t                     destruction of the inherited stock
                 + (1−δ̄)(1−q̄) v̂_t                 unfilled vacancies that survive
                 + (1−δ̄)ηq̄ θ̂_t                     slower filling leaves more unfilled
                 − (1−δ̄)(s̄/θ̄) û_t                 fewer employed, fewer reposts
                 + (ē/v̄) ê_{t+1}                    entry
```

**Interpretation of each weight.**

- **X̄/v̄: the share of next period's stock that comes from this period's inherited
  vacancies.** Destruction removes these one for one. It is large at D1 because most
  vacancies are reposted positions.
- **(1−δ̄)(1−q̄)** is the fraction of vacancies that are neither filled nor destroyed.
- **(1−δ̄)ηq̄.** Tighter markets lower the fill rate, so more vacancies carry over. After a
  δ shock θ̂ < 0, so this term is negative: vacancies fill faster.
- **(1−δ̄)s̄/θ̄.** Reposts come from the employed. More unemployment means fewer of them.
- **ē/v̄: entry's share of the stock.** It is a thin slice at D1.

| | D1 | CK-like |
|---|---|---|
| X̄/v̄ | **0.930** | **0.174** |
| (1−δ̄)(1−q̄) | 0.197 | 0.169 |
| (1−δ̄)ηq̄ | 0.480 | 0.480 |
| (1−δ̄)s̄/θ̄ | 0.055 | 0 |
| ē/v̄ | **0.073** | **0.831** |

**At free entry, (L4) does not determine v. It determines ê**, given v̂ from (L5).

### 3.5 Tightness

```
(L5)   v̂_t = θ̂_t + û_t
```

---

## 4. iid destruction shock (d_0 > 0, d_t = 0 for t ≥ 1)

### Step 1: tightness

- For t ≥ 1 there is no shock, so θ̂_t = 0 satisfies (L1)–(L2). (L2) is forward-looking and
  has no state variable, so this is the unique bounded solution.
- At t = 0, E_0θ̂_1 = 0 and K̂_1 = 0, so (L2) reduces to K̂_0 + ηθ̂_0 = −d_0/(1−δ̄). With
  (L1):

```
θ̂_0 = −(1/η) [ 1/(r+δ̄) + 1/(1−δ̄) ] d_0
```

**Why tightness falls.** The value of a hire is barely affected (the survival term,
1/(1−δ̄) ≈ 1). The rental rises by 1/(r+δ̄), which is 167 at D1. For the expected cost of a
hire, K/q, to stay equal to the value of a hire, the fill rate q must rise. So θ must fall.
The rental term dominates.

### Step 2: date 0

- u_0 is predetermined, so û_0 = 0.
- By (L5), v̂_0 = θ̂_0 < 0.
- The inherited stock is fixed, so entry absorbs the cut: ê_0 = (v̄/ē)θ̂_0.

### Step 3: date 1

- θ̂_1 = 0, because the rental is back to normal.
- By (L3), û_1 = (f̄/τ̄)d_0 + ω|θ̂_0| > 0. Two sources: the jobs destroyed at the end of
  date 0, and the hires that did not happen at date 0.
- By (L5), **v̂_1 = û_1 > 0**.
- By (L4), entry rises to deliver it. The inherited stock is lower, mostly because
  vacancies were cut and filled faster at date 0, not because of destruction.

### Step 4: date 2 on

θ̂_t = 0, û_{t+1} = λû_t and v̂_t = û_t. The economy returns along the ray v = θ̄u.

### Numbers (per 1 % rise in δ_0, in %)

| t | D1: θ̂ | û | v̂ | ê | CK-like: θ̂ | û | v̂ | ê |
|---|---|---|---|---|---|---|---|---|
| 0 | −0.76 | 0 | **−0.76** | −10.5 | −1.56 | 0 | **−1.56** | −1.9 |
| 1 | 0 | 0.16 | **+0.16** | +9.3 | 0 | 0.68 | **+0.68** | +2.0 |
| 2 | 0 | 0.09 | +0.09 | +0.9 | 0 | 0.38 | +0.38 | +0.3 |

At D1, of û_1 = 0.16, destruction contributes (f̄/τ̄)d_0 = 0.04 and the date-0 hiring
collapse ω|θ̂_0| = 0.12.

### Conclusion for iid shocks

- Vacancies fall only at date 0, when unemployment has not yet moved.
- From date 1 on, u and v move together along the ray v = θ̄u.
- **The Beveridge response is positively sloped after impact, at every δ̄**, including D1.
  In the (u, v) plane the economy drops straight down, jumps up and to the right, and
  slides back along the ray.

---

## 5. Persistent destruction shock (d_t = ρ^t d_0, 0 < ρ < 1)

### Step 1: tightness

Guess θ̂_t = a·d_t. Insert it, (L1) and E_t d_{t+1} = ρd_t into (L2), and solve for a:

```
        1/(1−δ̄)  +  [ 1 − Bρ( (1−s̄−ϕf̄) − (1−ϕ)q̄ ) ] / (r+δ̄)
a = − ─────────────────────────────────────────────────────────
              η  −  Bρ[ (1−s̄−ϕf̄)η − ϕf̄(1−η) ]
```

With ρ = 0 this is the iid case. Persistence adds two forces, both through the rental
staying high next period:

- It raises the value of already having a worker, which saves next period's higher hiring
  cost. This damps the fall in θ.
- It erodes next period's match surplus. This amplifies the fall.

The denominator also accounts for tightness being low next period too.

| | D1 | CK-like |
|---|---|---|
| a at ρ = 0.592 | **−292** | **−47** |
| a, iid | −280 | −50 |

Tightness at each date is θ̂_t = a·ρ^t·d_0. It recovers at rate ρ.

### Step 2: unemployment

From (L3) with û_0 = 0:

```
û_{t+1} = λ û_t + (f̄/τ̄ + ω|a|) ρ^t d_0
```

Each period adds an inflow from destruction (f̄/τ̄) and from lost hiring (ω|a|). Both decay
at rate ρ, and unemployment carries each one forward at rate λ.

### Step 3: vacancies, horizon by horizon

By (L5), v̂_t = −|a|ρ^t d_0 + û_t.

- **t = 0:** v̂_0 = −|a|d_0 < 0, with û_0 = 0.
- **t = 1:** v̂_1 = [f̄/τ̄ − |a|(ρ − ω)] d_0. So **vacancies fall at h = 1 iff
  |a|(ρ − ω) > f̄/τ̄.**
  - D1: 292 × 0.428 = 125 > 13.3.
  - CK-like: 47 × 0.428 = 20 > 13.7.
  - Both hold. h = 1 does not separate the two economies.
- **Later horizons.** Divide v̂_t < 0 by ρ^t:

  |a|ρ > (f̄/τ̄ + ω|a|) · Σ_{i=0}^{t−1} (λ/ρ)^i.

  The right side rises with t. So once vacancies turn positive, they stay positive.
- **Vacancies stay below steady state forever iff ρ > λ and |a|(ρ − λ − ω) ≥ f̄/τ̄.** This
  needs ρ > λ + ω = 0.72. At ρ = 0.592 vacancies eventually turn positive in both
  economies, but by very different amounts.

**Numbers (per 1 % rise in δ_0, ρ = 0.592, in %):**

| t | D1: θ̂ | û | **v̂** | CK-like: θ̂ | û | **v̂** |
|---|---|---|---|---|---|---|
| 0 | −0.79 | 0 | **−0.79** | −1.46 | 0 | **−1.46** |
| 1 | −0.47 | 0.17 | **−0.30** | −0.87 | 0.66 | **−0.20** |
| 2 | −0.28 | 0.19 | **−0.09** | −0.51 | 0.76 | **+0.25** |
| 3 | −0.16 | 0.16 | **+0.00** | −0.30 | 0.66 | **+0.36** |
| 4 | −0.10 | 0.13 | **+0.03** | −0.18 | 0.51 | **+0.33** |

- **At D1**, vacancies are below steady state while unemployment is rising (months 0–2).
  The later positive part is about 0.03 %, against a trough of −0.79 %.
- **At δ = τ**, vacancies are well above steady state for most of the period in which
  unemployment is high.

### Step 4: one number for the whole path

Regress ṽ on ũ along the impulse response (levels). Using the solutions above and two
geometric sums, the slope is

```
slope = θ̄ [ 1 − |a| ρ(1−λ²) / ( (f̄/τ̄ + ω|a|)(1 + ρλ) ) ]
```

It is negative iff |a|[ρ(1−λ²) − ω(1+ρλ)] > (f̄/τ̄)(1+ρλ).

| | iid | ρ = 0.592 | ρ needed for slope < 0 |
|---|---|---|---|
| D1 | +0.51 | **−0.24** | ρ > 0.37 |
| CK-like | +0.51 | **+0.17** | none in (0, 1) |

---

## 6. Where δ̄ relative to τ̄ enters

Walk back through the weights:

- **The unemployment equation (L3) is the same in both economies.** λ, f̄/τ̄ and ω depend
  on f̄ and τ̄, which are held fixed, not on the δ/s split.
- **The non-rental terms of the job creation condition (L2)** also differ little.
- **What differs is the rental weight 1/(r + δ̄)** in (L1): 167 at D1, 29 at δ = τ. It
  scales the tightness response a, and nothing else.

So the comparison is between two sensitivities to the same one-point rise in δ:

- **Tightness** falls by roughly 1/(r + δ̄), times a factor that is about the same in both
  economies.
- **Unemployment** rises by f̄/τ̄, because unemployment is a small stock, of order τ̄.

Vacancies, v = θu, fall when the first beats the second. That is the comparison of r + δ̄
with τ̄. Two things follow:

- At D1, r + δ̄ is about a fifth of τ̄. The rental of a long-lived vacancy is far more
  sensitive to destruction than unemployment is.
- With δ = τ (CK, s = 0), r + δ̄ exceeds τ̄. The vacancy dies with the match, so its rental
  can never be the more sensitive of the two.

Holding everything else fixed, the program reports the path condition in the form
(r + δ̄)/τ̄ < Γ_P at ρ = 0.592:

| | (r + δ̄)/τ̄ | Γ_P |
|---|---|---|
| D1 | 0.19 | 0.60 |
| CK-like | 1.11 | 0.52 |

Γ_P barely moves with δ̄, so the side of the boundary each economy lands on is decided by
(r + δ̄)/τ̄.

---

## 7. Remarks

1. **κ > 0** (a per-match cost, our x_v < 1). The rental terms in (L2) are scaled by
   x_v = (K̄/q̄)/(κ + K̄/q̄). At x_v = 0.5 the numbers barely change: the path slope at
   ρ = 0.592 is −0.23 at D1 and +0.17 at δ = τ.
2. **Why the earlier steady-state version was misleading.** A permanent shock (ρ → 1) hides
   the fact that the iid response is positively sloped after impact. The dynamic version
   shows that persistence is essential.
3. **For later (not in the draft yet).**
   - `prop:ds_asymmetry` Part 3's sufficient condition ("δ̄_e/τ̄ small") fails in the iid
     free-entry case of §4, where entry refills the stock at t+1.
   - The additional ingredients are σ > 0, ζ > 0 (variety lowers S), p₀ > 0, and the e ≥ 0
     bound.
