# Model Equations — Baseline and Partial-Reposting Extension
**Created:** September 28, 2026 · **Timing audit:** October 1, 2026 · **Branch:** `costly_vacancy_reposting`

This file lists every model equation with an economic interpretation, **dated explicitly** in
the timing of the draft (`Draft/Draft.tex`, `fig:Timing`). The draft is the reference model.
`Programs baseline/run_solution_core.jl` implements the same 33 equations but **deviates from
the draft's timing in two places** (see "Code status" below). Each equation carries a
**Code** line saying whether the current code matches. This file is the specification the code
restructuring in task R8 should implement.

- **Block 1** is the model with costless full reposting (Λ_r = 1).
- **Block 2** gives the equations that change under partial reposting, the extension argued for
  in [`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).

**Notation.** Subscripts are dates. E_t conditions on information at the *decision stage* of
period t (Stage 2 onward, below). Symbol map draft ↔ code:

| Draft | Code | Note |
|---|---|---|
| δ_t | `δ*δbar` | code shocks are multiplicative, SS = 1 |
| s_t | `s*sbar` | |
| χ^c_t, χ_m | `x_c`, `f_m` | χ_m is the bound of F, **not** `x_m` (bound of G) |
| Λ_t = F(χ^c_t) | `Λ` | survival probability |
| v_{pre,t} | `v_pret` | |
| m_{t+1} | `β*λp/λ` | |
| ν^f_t | `ν_f` | |

Composite parameters: β = 1/(1+r), ξ = 1/ξ_inv, μ = ε/(ε−1), ψ_c = ψ/(ψ+1).

---

# BLOCK 1 — Baseline model (costless full reposting, Λ_r = 1)

## Timing within period t (monthly; draft `fig:Timing`)

1. **Stage 1 — aggregate state.** z_t, δ_t, s_t are realized and known; every date-t decision
   conditions on them (the paper's convention throughout). δ_t and s_t are *rates* that take
   effect at Stage 5.
2. **Stage 2 — exit, reactivation, entry.** Each retailer draws χ ~ F and exits if χ > χ^c_t.
   Survival into t is (1−δ_{t−1})·F(χ^c_t) = 1 − δ_{e,t}. Recruiters at surviving product lines
   decide which positions vacated at the end of t−1 to reactivate (Block 2). Vacancy entry
   e_t = G(Q_t); new product lines N^e_t form. **u_t, v_{pre,t}, N_t are post-exit stocks.**
3. **Stage 3 — bargaining and production** with N_t product lines and L_t = 1 − u_t workers.
4. **Stage 4 — matching** out of v_t = v_{pre,t} + e_t and u_t.
5. **Stage 5 — end of t.** Separation at rate s_t and exogenous destruction at rate δ_t take
   effect; only *which* matches and product lines they hit is revealed here. s_t hits only
   matches that produced in t (Convention A): matches formed at Stage 4 are exposed to δ_t and
   to the t+1 exit draw, but first face s at the end of t+1.

**Two dating rules follow.**

- **1 − δ_{e,t+1} = (1 − δ_t)·F(χ^c_{t+1}).** The exogenous factor is dated t because δ_t is
  the period-t rate, known at t and applied at the end of t. The endogenous factor is dated t+1 because the continuation draw
  is compared with χ^c_{t+1}, which depends on period-t+1 profits and firm value
  (`eq:cutoff_eq`). Both act before production in t+1, so δ_e carries the date of the period
  in which survivors produce. Survival from t to t+1 is the same object in every value function
  and every law of motion.
- **Stocks dated t+1 are not predetermined at t.** N_{t+1}, u_{t+1}, v_{pre,t+1} depend on
  χ^c_{t+1}, which depends on N_{t+1} (the fixed point in the Prop. 5 Part 2 proof). The
  predetermined objects are the **pre-exit** stocks: B_{t+1} ≡ N_t + N^e_t for product lines,
  and (u_t, v_t) for the labor market, from which pre-exit employment (1−s_t)(1−u_t) + f_t u_t
  and pre-exit vacancies (1−q_t)v_t + s_t(1−u_t) are built.

**Exit and reactivation share Stage 2, but not an agent or a draw.** Exit is the retailer's
decision on its product line against χ^c_t. Reactivation is the recruiter's decision on a
vacated position against Q_t. Recruiters hold vacancies, so reactivation is theirs. Retailer
exit gates reactivation, because a position serving a withdrawn product line is destroyed
regardless; that is the (1 − δ_{e,t}) factor in f[22]. Since Q_t is set at the same node
(e_t = G(Q_t)), reactivation conditions on *realized* Q_t, which is what makes Λ_{r,t}
procyclical. The reposting option accrues to the recruiter's J, not the retailer's ν^f, so the
dependence stays triangular and f[1]/f[4]/f[18] are unchanged in Block 2.

## Code status — the two timing deviations (October 1, 2026)

The code is GS-style end-of-period timing in its laws of motion but draft timing in its value
equations. Fixing this is part of R8.

1. **Exit dating.** Code f[2] builds `δ_e = 1 − (1−δ·δbar)·Λ` with the *current* cutoff Λ_t and
   applies it to the t→t+1 transition in f[22]–f[24]. In draft terms that is (1−δ_t)·F(χ^c_t)
   where the draft has (1−δ_t)·F(χ^c_{t+1}): the exogenous factor is dated right, the
   endogenous factor one period early. The value equations f[3]–f[5] use
   `SDF_surv = (1−δ·δbar)·Λp`, i.e. F(χ^c_{t+1}), as the draft does. So the code prices
   survival at one date and depletes stocks at another. Code production also uses all N_t,
   including firms whose period-t draw failed the cutoff.
2. **Predetermined stocks.** Code states are the post-exit stocks [u, N, v_pret]. Under draft
   timing they must be pre-exit stocks, because the t+1 stocks depend on χ^c_{t+1}.

**Not a deviation: shock observation.** Code δ and s are period-t states, so θ_t, e_t, x_c,t
respond to them. That is the paper's convention (aggregate state known at the start of t;
corrected Oct 1, after an earlier version of this file called it a deviation).

**Restructuring spec (R8).** States are the stocks after the end-of-(t−1) incidence of δ_{t−1}
and s_{t−1} but before the period-t exit draw:
Ñ_t = (1−δ_{t−1})(N_{t−1} + N^e_{t−1}),
Ẽ_t = (1−δ_{t−1})[(1−s_{t−1})(1−u_{t−1}) + f_{t−1}u_{t−1}],
Ṽ^u_t = (1−δ_{t−1})(1−q_{t−1})v_{t−1},  Ṽ^s_t = (1−δ_{t−1})s_{t−1}(1−u_{t−1}),
plus z_t, δ_t, s_t. Then at t: N_t = Λ_tÑ_t, u_t = 1 − Λ_tẼ_t,
v_{pre,t} = Λ_t(Ṽ^u_t + Λ_{r,t}Ṽ^s_t), X^r_t = Λ_tṼ^s_t·α·M_t, all same-period controls solved
jointly with χ^c_t. The lagged rates are absorbed into the states, so no lag needs tracking.
Every §5.3 number must be regenerated afterwards; bundle with the D1 cascade.

## Shocks

| Shock | Meaning | Process | Realizes |
|---|---|---|---|
| z | technology | log z_{t+1} = ρ_z log z_t + ε^z_{t+1} | Stage 1 of t+1 |
| δ | exogenous product-line destruction | log δ_t = ρ_δ log δ_{t−1} + ε^δ_t | Stage 1 of t; applied at Stage 5 |
| s | match separation | log s_t = ρ_s log s_{t−1} + ε^s_t | Stage 1 of t; applied at Stage 5 |

Orthogonal innovations (Cholesky, z first; `part6b`). The z→δ_e link runs through χ^c, not
through shock covariance. Monthly calibration from `part6b`: ρ_z = 0.902, σ_z = 0.0092,
ρ_δ = 0.592, σ_δ = 0.0669, ρ_s = 0.8741, σ_s = 0.0854. **Code:** processes match; the
observation date differs (deviation 2).

## Continuation-cost distribution

```
F(χ) = (1 − p_0) + p_0·(χ/χ_m)^ψ   for χ ∈ [0, χ_m];   F(χ) = 1 for χ > χ_m
```

Mass 1 − p_0 at zero cost; a power-law tail with shape ψ up to χ_m. Λ_t ≡ F(χ^c_t).

## Equations f[1]–f[33]

### Exit — f[1]–f[2]

**f[1] Exit threshold** (`eq:cutoff_eq`)
```
χ^c_t = Y^c_t·(μ−1)/(μ·N_t) + ν^f_t
```
A product line continues iff its draw is below current gross profit per line plus its
post-dividend value. The draw is paid at Stage 2 of t, so the comparison uses period-t profit.
**Code:** matches in form (f[1]), but in the code N_t, Y^c_t do not respond to χ^c_t (deviation 1).

**f[2] Total exit rate** (`eq:delta_e_lom`)
```
δ_{e,t} = 1 − (1 − δ_{t−1})·Λ_t,     Λ_t = F(χ^c_t)
```
Exogenous destruction realized at the end of t−1, plus endogenous exit at Stage 2 of t.
**Code:** `δ_e = 1 − (1 − δ·δbar)·Λ` pairs the code's period-t δ with Λ_t and uses it for the
t→t+1 transition (deviation 1).

### Asset-pricing equations — f[3]–f[5]

Survival factor for t→t+1: **1 − δ_{e,t+1} = (1 − δ_t)·Λ_{t+1}**, inside E_t because both δ_t
and χ^c_{t+1} are unknown at the decision stage of t.

**f[3] Job creation condition** (`eq:jcc`, `eq:jcc_eq`)
```
κ + K_t/q_t = E_t[ m_{t+1}·(1−δ_t)Λ_{t+1}·(J_{t+1} − Q_{t+1}) ]
J_t − Q_t  = (1−ϕ)(w^int_t − K_t − b) − ϕ·θ_t(K_t + q_t κ) + (1 − s_t)(κ + K_t/q_t)
```
The average hiring cost equals the discounted, survival-weighted surplus of a match that
produces from t+1. The surplus is the recruiter's share of current net revenue plus the
re-hiring cost saved if the match survives separation. In Block 1 the surplus can be
substituted into the first line without nesting. **Code:** matches; f[3] is the substituted
form, with `SDF_surv = (1−δ·δbar)·Λp`.

**f[4] Business-formation Euler** (`eq:N_euler_eq`)
```
ν^f_t = E_t[ m_{t+1}·(1−δ_t)Λ_{t+1}·(ν^f_{t+1} + d^f_{t+1}) ]
```
**Code:** matches.

**f[5] Flow value of a vacancy** (`eq:K`)
```
K_t = Q_t − E_t[ m_{t+1}·(1−δ_t)Λ_{t+1}·Q_{t+1} ]
```
**Code:** matches.

### Wages and vacancy creation — f[6]–f[8]

**f[6]** `w^int_t = ρ(N_t)·z_t/μ` (`eq:recruiter_compensation`). **Code:** matches.

**f[7] Nash wage** (`eq:wage_eq`; derivation in the wage appendix)
```
w_t = ϕ·(w^int_t − K_t + f(θ_t)(κ + K_t/q_t)) + (1−ϕ)·b
```
f(θ)(κ + K/q) = θ(K + qκ). Under Convention A the household's continuation weight is
(1−δ_{e,t+1})(1 − s_t − f_t), matching the recruiter's, so no (1−s_t) multiplies the f term.
**Code:** matches.

**f[8] Vacancy creation** (`eq:Q`) `Q_t = x_m·e_t^(1/ξ)`, equivalently e_t = G(Q_t).
**Code:** matches. This is the position-level *creation* cost from G; reposting does not reuse
it (Block 2).

### Labor market — f[9]–f[12]

**f[9]** `θ_t = v_t/u_t` · **f[10]** `q_t = A·θ_t^(−η_L)`, `f_t = A·θ_t^(1−η_L)` ·
**f[11]** `L_t = 1 − u_t` · **f[12]** `L_t = L^c_t + L^e_t`. All period-t, post-exit.
**Code:** match in form.

### Goods market — f[13]–f[15]

**f[13]** `ρ_t = N_t^(1/(ε−1))` · **f[14]** `λ_t = C_t^(−σ)` · **f[15]** `Y^c_t = ρ_t·z_t·L^c_t`.
**Code:** match in form; N_t is the post-exit stock in the draft (deviation 1).

### Resource constraint and identities — f[16]–f[21]

**f[16] Resource constraint** (`eq:rc`, `eq:agg_fixed_costs`, `eq:X_total`)
```
Y^c_t = C_t + X_t + X^c_t
X_t   = e_t·Q_t·ξ/(ξ+1) + κ·q_t·v_t
X^c_t = N_t·∫_0^{χ^c_t} χ dF = N_t·ψ_c·χ^c_t·p_0·(χ^c_t/χ_m)^ψ
```
Continuation costs are paid at Stage 2 of t by the N_t survivors. ξ/(ξ+1) = 1/(1+ξ_inv).
**Code:** `X_c = N·p_0·ψ_c·x_c` drops the (χ^c/χ_m)^ψ factor. Relative error p_0/(Λ−1+p_0),
0.27% at Λ ≈ 0.9986; the draft now uses the exact form. Upgrade in R8.

**f[17]** `N^e_t = z_t·L^e_t/f_e` · **f[18]** `ν^f_t = ρ_t·f_e/μ` (`eq:free_entry`) ·
**f[19]** `Y_t = C_t + ν^f_t·N^e_t` (`eq:gdp`) · **f[20]** `N_t·d^f_t = Y^c_t(μ−1)/μ − X^c_t` ·
**f[21]** `v_t = v_{pre,t} + e_t`. **Code:** match.

### Laws of motion — f[22]–f[24] (Convention A)

**f[22] Pre-entry vacancies** (`eq:v_lom`)
```
v_{pre,t+1} = (1 − δ_{e,t+1})·[ (1 − q_t)·v_t + s_t·(1 − u_t) ]
```
Unmatched vacancies plus positions vacated by separations of matches that produced in t, both
surviving into t+1. Under full reposting every such position returns as a vacancy. The
reposting base is the same stock that separations are drawn from in f[23].

**f[23] Unemployment** (`eq:u_lom`)
```
u_{t+1} = [1 − (1 − δ_{e,t+1})·f_t]·u_t + τ_{t+1}·(1 − u_t),   τ_{t+1} = δ_{e,t+1} + s_t(1 − δ_{e,t+1})
```
Equivalently L_{t+1} = (1 − δ_{e,t+1})[(1 − s_t)L_t + f_t u_t]. New matches face exit but not
separation before they produce. With f[22], positions are conserved up to exit:
v_{pre,t+1} + L_{t+1} = (1 − δ_{e,t+1})(v_t + L_t).

**f[24] Product lines** (`eq:firm_law_motion`)
```
N_{t+1} = (1 − δ_{e,t+1})·(N_t + N^e_t)
```
Incumbents and entrants face exogenous destruction δ_t and the t+1 exit draw.

**Code (f[22]–f[24]):** Convention A is implemented (since Sept 6, 2026), but with δ_e built
from Λ_t rather than Λ_{t+1}, and with u, v_pret, N as predetermined states (deviations 1 and 3).

### Labor share, observables, shocks — f[25]–f[33]

**f[25]** `ls_t = w_t L_t/Y_t` · **f[26]** `labor_prod_t = Y_t/(ρ_t L_t)` ·
**f[27]–f[30]** `C/ρ`, `Y/ρ`, `Y^c/ρ`, `w/ρ` (price-deflated) · **f[31]–f[33]** shock
processes (table above). **Code:** match in form.

---

# BLOCK 2 — Partial-reposting extension

## The margin

On separation, each vacated position independently draws a reactivation cost α·χ, with χ ~ F,
the **same** distribution as the continuation cost. The decision is one-time: reactivate the
slot (pay α·χ, obtain a vacancy worth Q) or retire it permanently (to have it back later, the
recruiter must create a fresh position at cost drawn from G). Reactivate iff Q ≥ α·χ.

- **Threshold is Q, not Q − K.** The recruiter obtains an unfilled vacancy worth Q, and the
  Bellman already discounts it. A flow-value K threshold would pin Λ_r at its floor 1 − p_0;
  the recurring-maintenance story that would justify it was rejected (note §7.3b).
- **Per-position draws.** Under the recruiter/retailer segmentation, positions belong to
  recruiters and retailers rent labor competitively, so only the aggregate vacancy stock
  enters matching and the law of large numbers delivers Λ_r. Per-position draws are the
  natural reading (each position is a separate asset). The earlier argument that per-firm
  draws would break DS-CES symmetry attributed positions to retailers; it is not needed.
- **One new parameter α**, no new distribution (note §7.3a: F prices *operating* capacity,
  G prices *creating* it).
- **Notation.** Λ_r parallels the survival probability Λ. Not λ, which is marginal utility.
  Draft symbols: Q^rep for the option value; shortfall and mean-cost integral written inline.
  Do not reintroduce Ψ, D, or M in the draft (clashes with ψ, D^int, 𝒟, m). This file uses
  M_t ≡ ∫_0^{Q_t/α} χ dF as shorthand only.

## New equation

**f[new] Reposting rate** (draft: `eq:Lambda_r`, in the recruiter block; cited in `def:equilibrium`)
```
Λ_{r,t} = F(Q_t/α) = min[ (1 − p_0) + p_0·(Q_t/(α·χ_m))^ψ , 1 ]
```
Positions vacated at the end of t−1 are reactivated at Stage 2 of t against the realized Q_t.
Procyclical, so position destruction is countercyclical. Λ_r ∈ [1 − p_0, 1]: the zero-cost
atom always reactivates, so p_0 caps the margin. Saturates at 1 for α ≤ Q̄/χ_m. Track it as a
control so its impulse response can be reported (R9).

**Two routes to the costless benchmark.** α → 0 shrinks the cost scale. p_0 → 0 collapses F to
a point mass at zero, so every draw is free: Λ_r = 1, Q^rep = Q, X^r = 0. The p_0 route is
why `prop:independence` needs no fourth condition. Note Λ_r = 1 with α > 0 and p_0 > 0 is
*quantity* saturation, not the costless model: all positions repost but still pay α·χ.

## Modified equations

**f[22] Pre-entry vacancies — MODIFIED** (`eq:v_lom`)
```
v_{pre,t+1} = (1 − δ_{e,t+1})·[ (1 − q_t)·v_t + Λ_{r,t+1}·s_t·(1 − u_t) ]
```
Λ_{r,t+1} and (1 − δ_{e,t+1}) share a date: both are Stage-2 decisions of t+1, by different
agents against different thresholds. Factoring (1 − δ_{e,t+1}) out front is legitimate because
reactivation is a sub-event of survival. The non-reactivated fraction 1 − Λ_{r,t+1} is the
missing middle: position destruction at surviving product lines, so position conservation now
leaks by exactly (1 − δ_{e,t+1})(1 − Λ_{r,t+1})s_t(1 − u_t). Scope choice: only the separation
inflow is attenuated, not the unmatched carry-forward (1 − q_t)v_t.
**Code:** not implemented (R8).

**f[3] Recruiter block — MODIFIED** (R7; `eq:value_recruiter_matched`, `eq:Qrep`,
`eq:value_recruiter_surplus`, `eq:jcc`, `eq:surplus_wage`)

Only the matched recruiter's separation branch changes:
```
J_t = w^int_t − w_t + E_t[ m_{t+1}(1−δ_t)Λ_{t+1}·( s_t·Q^rep_{t+1} + (1 − s_t)·J_{t+1} ) ]
Q^rep_t = E[max(Q_t − αχ, 0)] = Q_t·Λ_{r,t} − α·M_t
```
Recruiter surplus, with the expected discounted shortfall
𝓡_t ≡ E_t[ m_{t+1}(1−δ_t)Λ_{t+1}·(Q_{t+1} − Q^rep_{t+1}) ]:
```
J_t − Q_t = w^int_t − w_t − K_t + (1 − s_t)(κ + K_t/q_t) − s_t·𝓡_t
κ + K_t/q_t = E_t[ m_{t+1}(1−δ_t)Λ_{t+1}·(J_{t+1} − Q_{t+1}) ]
```
The shortfall Q − Q^rep = Q(1 − Λ_r) + α·M is the value lost relative to free reposting:
forgone vacancies plus reactivation costs paid. **Nesting:** α → 0 or p_0 → 0.

**Representation (revised Oct 1).** Substitute the wage-substituted surplus at t+1 into the job
creation condition and **track 𝓡_t as a variable** with its own one-period equation. J_t and w_t
then drop out of the system (they can be recovered from the surplus and wage equations). The
nested expectation arises only if 𝓡_{t+1} is itself expanded, which composes two one-period
lookaheads into a t+2 term. Tracking either J − Q or 𝓡 avoids it with one forward variable;
tracking 𝓡 is preferred because the job creation condition then keeps its Block 1 form plus one
term, which is also the smallest change to code f[3]. This is how `def:equilibrium` is written.

**Deriving M.** With threshold x = Q/α, only the power-law part integrates:
```
M = ∫_0^x χ dF = (p_0ψ/χ_m^ψ)∫_0^x χ^ψ dχ = ψ_c·x·p_0·(x/χ_m)^ψ = ψ_c·x·(Λ_r − 1 + p_0)
α·M = ψ_c·Q·(Λ_r − 1 + p_0)        (interior, Q/α ≤ χ_m; α cancels except inside Λ_r)
```

| Object | Value | Per what |
|---|---|---|
| α·M | ψ_c·Q·(Λ_r − 1 + p_0) | per **vacated** position (what X^r needs) |
| α·M/Λ_r | ψ_c·Q·(Λ_r − 1 + p_0)/Λ_r | per **reactivated** position |
| ψ_c·Q | ψ_c·Q | per position that actually **pays** |

α·M already contains Λ_r (partial expectation), so Λ_r must not be applied again in X^r. As
Λ_r → 1 − p_0, α·M → 0. At Q = α·χ_m the interior and saturated forms (α·ψ_c·χ_m·p_0) agree.

**f[7] Nash wage — MODIFIED**
```
w_t = ϕ·(w^int_t − K_t + f(θ_t)(κ + K_t/q_t) − s_t·𝓡_t) + (1−ϕ)·b
```
The household surplus has no reposting term (a separated worker is unemployed either way), so
the shortfall enters bargaining only through the recruiter's surplus, and the worker bears
share ϕ of it. It works like a firing cost, except that separation is exogenous. Derived in
the draft's wage appendix.

**Wage-substituted surplus** (`eq:surplus_wage`; substituted at t+1 into the job creation condition)
```
J_t − Q_t = (1−ϕ)(w^int_t − K_t − b) − ϕ·θ_t(K_t + q_t κ) + (1 − s_t)(κ + K_t/q_t) − (1−ϕ)·s_t·𝓡_t
```
⚠️ The coefficient on 𝓡_t is **(1−ϕ)**, not 1. In the draft (Oct 1): 𝓡_t is `eq:repost_shortfall`,
and `eq:wage_eq`, `eq:surplus_wage`, the substituted `eq:jcc_eq`, the wage appendix (`app:wage`), and the
steady-state appendix (`eq:jcc_wage_ss`, wage, 𝓡, e, X^r) all carry it.

**f[16] Resource constraint — MODIFIED** (`eq:agg_repost_costs`, `eq:rc`, `eq:gdp`)
```
Y^c_t = C_t + X_t + X^c_t + X^r_t
X^r_t = (1 − δ_{e,t})·s_{t−1}(1 − u_{t−1})·α·M_t = (1 − δ_{e,t})·s_{t−1}(1 − u_{t−1})·ψ_c·Q_t·(Λ_{r,t} − 1 + p_0)
```
Costs are paid at Stage 2 of t, when reactivation happens against the realized Q_t, on the
positions vacated at the end of t−1 at surviving product lines. The count is the reposting
inflow of f[22] *without* Λ_r, because α·M already embeds it. **This timing is forced by the
Stage-2 reactivation node.** Paying at separation against E_t Q_{t+1} would make reactivation
condition on an expectation and lose the procyclicality the margin rests on. Under the R8
state vector, s_{t−1}(1 − u_{t−1}) is a function of states, so no extra lag is needed.
Use the exact integral: the X^c-style shortcut (p_0 for Λ_r − 1 + p_0) overstates X^r by 25%
at Λ_r = 0.9 and 2.5× at 0.7, and is wrong at the floor, where exact X^r = 0.

## Unchanged in form

- **f[1], f[4], f[18], f[20] (retailer block).** The firm Bellman (`eq:firm_bellman`) has only
  retail profit, the draw χ, and the continuation value, with no Q, K, J, or reposting object.
  Recruiter profits (w^int − w)L are a separate flow to the household. Reposting is downstream
  of exit, so these move only through general equilibrium.
- **f[2] δ_e.** Reposting destroys positions at surviving product lines, not product lines.
- **f[5] K.** The unfilled-vacancy Bellman does not involve reactivation.
- **f[23] unemployment.** Workers separate at τ regardless of reactivation.
- **f[24] product lines.** The variety margin is untouched: only δ moves N, which keeps the
  core δ–s asymmetry.
- **f[6], f[8]–f[15], f[17], f[19], f[21], f[25]–f[33].**

## Open items

1. **R8:** implement the restructured state vector (Code status above) together with Block 2.
   Upgrade X^c to the exact integral at the same time.
2. **R2:** verify on simulated data that α is separately identified from σ_s and from χ_m.
3. **R3 / Prop. 5 under partial reposting:** Parts 1–2 now assume costless reposting
   (p_0 = 0 or α → 0). A result for α > 0 (e.g. a threshold Λ_r* at which the s→v sign
   flips) is not derived.
4. **Draft:** extend the log-linear appendix (`app:loglin`), which
   is stated for costless reposting, if it is used for the partial-reposting model.
