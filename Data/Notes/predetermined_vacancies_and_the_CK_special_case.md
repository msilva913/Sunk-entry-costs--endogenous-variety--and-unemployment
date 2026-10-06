# A downward-sloping Beveridge response to destruction at ξ → ∞

> **Start with [`beveridge_free_entry_GS_proposition.md`](beveridge_free_entry_GS_proposition.md).**
> It states and proves the result in the minimal Gabrovski–Silva economy, with the
> sufficient condition r + δ < τ (r the discount rate). This note covers the richer nested
> economy, dynamics and robustness.

**Written:** October 5, 2026. **Rewritten:** October 6, 2026. The rewrite supersedes the
Oct 5 argument, whose threshold had a units error (Appendix B).
**Numbers from:** `../Programs baseline/beveridge_free_entry.jl` →
`beveridge_free_entry.txt`. The program is self-contained. It mirrors `calibrate_shares`
PATH B with the mechanism-runner targets (b_ratio = 0.9, x_v = 0.5, dest_elast_target = 5).
It asserts the steady-state JCC, checks the dynamic tightness coefficient at ρ = 1 against
a direct steady-state derivative, and checks the closed-form sign condition against the
simulated IRF.
**Draft objects:** `eq:delta_e_lom`, `eq:u_lom`, `eq:v_lom`, `eq:Kdef`, `eq:jcc_eq`,
`eq:entry_eq`, `app:loglin` ¶9b–9c, `app:comparison_D`.

---

## 0. The question and the answer

**Question.** We know a negative u–v response to a δ shock is possible when vacancy
creation is finitely elastic (CK). Is it possible at ξ → ∞?

**Answer: yes, in a nested version of our model. It requires three things.**

1. **Vacancies are durable assets.** The JCC then carries a user-cost term that a
   flow-cost DMP model lacks.
2. **Destruction is a small share of separations.** The path slope is negative iff
   δ̄_e/τ̄ < 0.42 at the calibrated persistence. At D1, δ̄_e/τ̄ = 0.087. It is never negative
   at δ_e = τ.
3. **The shock is not too transitory.** At D1 the condition is ρ_δ > 0.37 monthly (0.05
   quarterly), against a calibrated ρ_δ = 0.592.

At D1, the free-entry path slope is −0.23. A flow-cost DMP model with the same steady state
gives +0.46.

---

## 1. Timing

`eq:delta_e_lom` gives δ_e,t = 1 − (1 − δ_{t−1})Λ_t. An innovation to δ at date 0 raises
δ_e,1.

- **Date 0.** Nothing is destroyed. The JCC and `eq:Kdef` price 1 − δ_e,1, so θ_0 moves.
  u_0 and the inherited vacancy stock do not.
- **Date 1.** Matches and positions are destroyed at δ_e,1.

**Units.** d_t ≡ δ̃_e,t+1 is the level deviation of the destruction rate known at t, with
d_t = ρ^t d_0. A one-standard-deviation innovation is about 3.3 % of δ̄_e.

---

## 2. Why the vacancy accounting identity cannot answer this

`eq:v_lom` gives v_t = (1 − δ_e,t)X_t + e_t, with X_t the inherited stock plus reposts. At
ξ → ∞, e_t is whatever closes the identity once θ_t and u_t are known. The identity holds
by construction and signs nothing. The predetermined share therefore does not decide the
sign at free entry.

**The question reduces to: does θ̂ + û < 0?**

---

## 3. The free-entry equilibrium in five lines

```
(i)   Q_t = x_m                                     free entry pins the vacancy value
(ii)  K_t = x_m [1 − β(1 − δ_e,t+1)]                eq:Kdef  ⇒  K̂_t = d_t /(r + δ̄_e)
(iii) κ + K_t/q(θ_t) = β(1 − δ_e,t+1) S_{t+1}        JCC sets θ_t given K_t
(iv)  u_t from eq:u_lom
(v)   v_t = θ_t u_t ,   e_t = v_t − (1 − δ_e,t) X_t
```

Line (ii) is the crux. A vacancy is a durable asset with price x_m and rental
(r+δ_e)x_m/(1+r). A rise d in the destruction rate raises the rental by d/(r+δ̄_e) in
proportion. In a flow-cost DMP model, K is a constant and line (ii) is absent.

**The nested economy.** It keeps the draft's timing, laws of motion and JCC. It shuts off
risk aversion (m = β), the variety channel (w_int fixed), endogenous exit (Λ ≡ 1) and
reposting costs (Λ_r = 1). Each δ̄_e is calibrated as in PATH B. The free-entry economy
shares that steady state (x_m = Q̄).

---

## 4. Tightness

Guess θ̂_t = a·d_t. Write B ≡ β(1 − δ̄_e), C ≡ κ + K/q, χ_K ≡ (K̄/q̄)/C̄ = x_v and
S ≡ (1−ϕ)(w_int − K − b) + (1 − s − ϕf)C. Log-linearizing (iii) gives

```
a = − (F1 + F2 + F3) / Den

F1  = 1/(1 − δ̄_e)                                     the new match may not survive
F2  = χ_K/(r + δ̄_e) · [1 − Bρ(1 − s̄ − ϕf̄)]          the slot's rental rises, net of the
                                                      higher hiring cost a match saves later
F3  = Bρ(1 − ϕ)(K̄/C̄) / (r + δ̄_e)                     the matched recruiter bears the higher
                                                      rental, eroding the match surplus
Den = χ_K η + Bρ[ϕ f̄(1 − η) − (1 − s̄ − ϕf̄) χ_K η]
```

| Calibration | δ̄_e/τ̄ | F1 | F2 | F3 | **a** | a, flow-cost DMP |
|---|---|---|---|---|---|---|
| **D1** (3.20 %/yr) | 0.087 | 1.00 | 47.5 | 16.6 | **−287** | −4.4 |
| Code (7.54 %/yr) | 0.210 | 1.01 | 30.2 | 8.0 | −162 | −4.2 |
| BGM/GS (10 %/yr) | 0.282 | 1.01 | 24.8 | 6.1 | −130 | −4.1 |
| **CK-like** (δ_e = τ) | 1.000 | 1.03 | 8.8 | 1.8 | **−46** | −4.1 |

**Reading.**

- F1 is the only force in flow-cost DMP. It is about 1, so tightness barely responds.
- F2 + F3 are 64 at D1 and 11 at δ_e = τ. Both carry 1/(r+δ̄_e): **a long-lived asset's
  rental is more sensitive to its depreciation rate.**
- The user cost dominates F1 whenever χ_K ≫ r + δ̄_e. That is about 0.006 monthly at D1, so
  even a small sunk share of recruiting costs is enough (§7).

---

## 5. The destruction-date sign condition

Unemployment on date 1, from `eq:u_lom`:

```
ũ_1 = [ū f̄ + (1 − ū)(1 − s̄)] d_0  −  ū(1 − δ̄_e) f̄ (1 − η) θ̂_0
```

The first term is the separation effect. The second is the date-0 fall in job finding.
With θ̂_1 = aρ d_0 and v̂_1 = θ̂_1 + û_1,

```
v̂_1 / d_0  =  a · [ρ − (1 − δ̄_e) f̄ (1 − η)]  +  c_u/ū ,      c_u ≡ ū f̄ + (1 − ū)(1 − s̄)
```

Line by line:

1. **c_u/ū ≈ 13 is the separation push on u, and through v = θu on v.** Per unit of d it
   is set by τ̄ and f̄. It hardly depends on the δ/s split.
2. **a·ρ is tightness on date 1.** It is still depressed because the shock persists.
3. **−a(1−δ̄_e)f̄(1−η) is the extra unemployment from date-0 job finding.** It raises u
   and so offsets part of the tightness fall in v.
4. If ρ > (1−δ̄_e)f̄(1−η) ≈ 0.16, **vacancies fall on date 1 iff |a| > a\* ≡
   (c_u/ū)/[ρ − (1−δ̄_e)f̄(1−η)]**.

| Calibration | \|a\| | a\* | \|a\|/a\* | v̂_1/d_0 | **path slope** |
|---|---|---|---|---|---|
| **D1** | 287 | 31 | 9.3 | −110 | **−0.23** |
| Code | 162 | 31 | 5.2 | −56 | −0.12 |
| BGM/GS | 130 | 31 | 4.2 | −42 | −0.08 |
| **CK-like** | 46 | 32 | 1.45 | −6 | **+0.17** |

The path slope is Σṽũ/Σũ² over 48 months.

**On the destruction date, vacancies fall in every row, even at δ_e = τ.** The difference
across rows is the margin. At D1 the condition holds nine times over; at δ_e = τ it barely
holds.

---

## 6. From the destruction date to the path

After date 1, tightness decays at ρ_δ while unemployment decays at roughly
1 − (1−δ̄_e)f̄ − τ̄ ≈ 0.56. If the shock is short-lived, u outlasts the tightness fall. v = θu
then turns positive while u is still high, and the path slope can be positive even when
v̂_1 < 0.

- **D1:** vacancies are below steady state in months 0–2 and above from month 3. The path
  slope is −0.23.
- **CK-like:** the late positive phase dominates, giving +0.17.

**The path slope is negative iff ρ_δ exceeds:**

| Calibration | ρ_δ threshold (monthly) | quarterly |
|---|---|---|
| **D1** | **0.369** | 0.050 |
| Code | 0.444 | 0.088 |
| BGM/GS | 0.491 | 0.118 |
| CK-like | none in [0.05, 0.999] | — |

The calibrated ρ_δ = 0.592 (0.207 quarterly) is already low because BED Deaths are noisy,
and the D1 threshold sits well below it. So the result does not rest on the persistence
estimate.

---

## 7. How a lower δ_e (relative to τ) helps, and robustness

**The mechanism in one sentence.** Holding τ̄ and f̄ fixed, the unemployment push per unit
of destruction is fixed, while the tightness pull scales with 1/(r+δ̄_e). So the relevant
statistic is **(r + δ̄_e)/τ̄**: 0.19 at D1, 0.32, 0.39, and 1.11 at δ_e = τ.

- At ρ_δ = 0.592, the path slope is negative iff **δ̄_e < 0.0131 monthly (14.7 %/yr,
  δ̄_e/τ̄ = 0.42)**.
- **The gain saturates.** At D1, δ̄_e = 0.0027 is already below r ≈ 0.0033 monthly, so
  lowering δ̄_e further adds little.
- **δ_e = τ is the worst case.** It never gives a negative path slope at free entry, for
  any ρ_δ we tried.

**Robustness.** One input at a time, from the mechanism targets:

| varied | value | \|a\| D1 | slope D1 | slope CK-like |
|---|---|---|---|---|
| x_v (sunk share) | 0.05 | 165 | −0.13 | +0.20 |
| | 0.25 | 269 | −0.22 | +0.18 |
| | 0.5 | 287 | −0.23 | +0.17 |
| | 1.0 | 292 | −0.24 | +0.17 |
| b_ratio | 0.71 | 231 | −0.20 | +0.20 |
| | 0.95 | 346 | −0.26 | +0.14 |
| ρ_δ | 0.3 | 284 | **+0.08** | +0.31 |
| | 0.9 | 292 | −0.49 | +0.08 |
| | 0.99 | 294 | −0.57 | +0.06 |

- **x_v and b_ratio** move the size, not the sign. This holds down to a 5 % sunk share, as
  §4 predicts.
- **ρ_δ is the input that can flip the sign at D1**, and only below 0.37.

---

## 8. Contrast with flow-cost DMP and with CK

**Flow-cost DMP** (same steady state, K fixed) gives a path slope of +0.46 to +0.47 at every
δ̄_e. Tightness barely moves (a ≈ −4), so u and v rise together. This is the problem CK
describe for free entry. Their comparison case ("H/M with JD") is this model, not their
durable-vacancy model at ξ = ∞. In their words: "With free entry, large job separation
shocks generate a counterfactual positive correlation" (Coles and Kelishomi 2018, p. 120,
DOI 10.1257/mac.20150040).

**CK's own structure** (c = 0, all creation costs sunk, every separation destroys the job)
is our durable vacancy at δ_e = τ. As far as the converted text shows, they do not report
it at ξ → ∞ under separation shocks. In our nested economy that case has a positive path
slope at every persistence. The vacancy-value channel is weakest there because
r + δ_e is largest. That is consistent with CK needing ξ < ∞, though their calibration
(large surplus, Hosios) differs from ours.

**So the departure from CK is this.** The structure is the same: durable vacancies with
sunk creation costs. With δ_e measured from establishment deaths, free entry no longer
breaks the Beveridge curve, because the user cost of a long-lived vacancy reacts strongly
to destruction. This comes from GS's δ/s split plus measurement, not from the business
formation block.

---

## 9. What is not established

1. **The full model at ξ_inv → 0.** The nested economy drops variety, endogenous exit, the
   SDF and reposting costs.
   - Variety lowers the surplus after a δ shock (`app:loglin` ¶9b), which should push θ
     down further.
   - Endogenous exit scales the δ_e response, but scales u and θ alike.
   - The SDF and reposting effects are unsigned.

   `calibrate_shares_free_entry` and `steady_state_free_entry` exist. The dynamic run does
   not.
2. **Comparison D's ξ_inv = 0.1 run** shows vacancies rising for two months, at 7.54 %/yr
   with code timing. Within the period, the vacancy stock's elasticity to Q there is
   ξ·ē/v̄ = 1.75 (Appendix A), so ξ = 10 may not be close enough to free entry. The code also
   dates exit one period early (R8). Run the free-entry case after R8.
3. **e ≥ 0.** At date 0, free entry wants v_0 to fall with the inherited stock fixed, so
   entry absorbs the whole cut. At D1, entry would turn negative beyond a shock of 9.3 % of
   δ̄_e (58 % at δ_e = τ). A 1-SD shock (3.3 %) fits; a 3-SD shock does not. Past that point
   Q < x_m and the economy behaves as if ξ were finite. The linearization ignores it.

---

## Appendix A. Finite ξ: the stock elasticity (secondary)

At finite ξ, ẽ = ē·ξ·Q̂, so the vacancy stock's within-period elasticity to Q is ξ·ē/v̄,
with ē/v̄ = (δ̄_e/τ̄)/κ_v and κ_v = 1.2034.

| | ξ = 0.265 | ξ = 1 | ξ = 10 | ξ for unit elasticity |
|---|---|---|---|---|
| D1 | 0.019 | 0.073 | 0.726 | 13.8 |
| Code | 0.046 | 0.175 | 1.746 | 5.7 |
| CK-like | 0.220 | 0.831 | 8.310 | 1.2 |

At D1 and ξ = 1, the stock is less elastic than CK's estimated ξ = 0.265 makes it at
δ_e = τ. The reason is that reposting, not entry, supplies most vacancies. With costly
reposting (this branch), Λ_r responds to Q and adds a second term. That term is not
computed here.

---

## Appendix B. Revision log, and the units error

- **The units error.** The Oct 5 note set ê/δ̂_e > v̄^pre/ē. The right-hand side is the bar
  per unit *level* of δ_e. Per 1 % it is X̄/(v̄+L̄): 0.035 at D1 and 0.0065 at δ_e = τ.
  "62.9 times higher", "13 % per 1 % shock" and "Q̂ ≈ 12.8 %" are withdrawn.
- **Kept:** destruction hits only the inherited stock; ē/v̄ = (δ̄_e/τ̄)/κ_v.
- **Draft items affected.**
  - `app:loglin` ¶9c compares per-log loadings. The unemployment loading "vanishing in
    δ̄_e/τ̄" is a units artifact. Rewrite around §5 and §7.
  - In `app:comparison_D`, "ξ → ∞ ⇒ ê = ξQ̂ → 0" is wrong: ξQ̂ is ∞·0, and entry is the
    residual. Its free-entry reasoning via the inherited share should be replaced by §4–§5.
  - `eq:v_ll` omits −(1−δ̄_e)[v̄ q̃_{t−1} + s̄ ũ_{t−1}].
- **Context files flagged as superseded:** `CLAUDE.md`, `findings.md`, `decisions.md` D1,
  `README.md`.
- `predetermined_vacancy_share.jl` output label corrected; numbers unchanged.
