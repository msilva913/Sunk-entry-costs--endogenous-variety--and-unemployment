# Model Equations — Baseline and Partial-Reposting Extension
**Created:** September 28, 2026 · **Branch:** `costly_vacancy_reposting`
**Source of truth for Block 1:** `Programs baseline/run_solution_core.jl` (33-equation system,
audited May 20, 2026; SS-timing fixes September 5, 2026).

This file lists every model equation with an economic interpretation. It has two blocks.

- **Block 1** is the current model, with *costless full reposting* (Λ_r = 1). It transcribes
  `run_solution_core.jl` equation for equation.
- **Block 2** rewrites only the equations that change under *partial reposting*, the extension
  argued for in [`../Notes/delta_calibration_and_the_reposting_margin.md`](../Notes/delta_calibration_and_the_reposting_margin.md).
  Equations whose modified algebra requires a formal re-derivation are marked **⚠ to be
  derived**; their economic content is stated but the exact terms are not asserted, per N1/N15.

Notation throughout: unprimed = period *t*; primed = period *t+1* (expected). Shocks
z, δ, s are multiplicative deviations from SS = 1. Composite parameters: β = 1/(1+r),
ξ = 1/ξ_inv (entry elasticity), μ = ε/(ε−1) (markup).

---

# BLOCK 1 — Baseline model (full costless reposting, Λ_r = 1)

## Timing (monthly, Coles-Kelishomi stage structure)

⚠️ **Corrected Sept 30, 2026 against the draft's own timing figure (`fig:Timing`).** The
earlier version of this list placed endogenous firm exit at Stage 5. The draft puts it at
**Stage 2**, early in the period, and puts only the *exogenous* destruction shock δ_t and
match separation s_t at Stage 5. Getting this wrong would put the exit node in the wrong
place in the code (R8).

1. **Stage 1** — aggregate state realized: new realizations of (z_t, δ_t, s_t).
2. **Stage 2** — **endogenous exit and reactivation.** Each *retailer* draws χ ~ F and exits
   if χ > x_c,t, giving δ_e,t = 1 − (1 − δ_{t−1})·F(x_c,t) — note the **lagged** δ. Entry:
   e_t = G(Q_t), N_e,t forms. Simultaneously, **recruiters** whose product line survived decide
   which positions vacated at the end of t−1 to reactivate, drawing χ ~ F *per position* and
   reactivating iff α·χ ≤ Q_t, giving Λ_r,t = F(Q_t/α). u_t is post-exit, pre-match.
3. **Stage 3** — bargaining and production: wages, output, profits.
4. **Stage 4** — matching: m_t matches formed, using total v_t = v_pret + e_t.
5. **Stage 5** — end of t: match separation s_t and exogenous product destruction δ_t;
   u_{t+1}, v_{t+1}, N_{t+1} updated. The δ_t realized here enters δ_e,t+1 at the next
   Stage 2.

**Why exit and reactivation share Stage 2.** Both compare a draw from F to a threshold, and
both need the period-t state. Q_t is determined at this node (it is the same node that sets
e_t = G(Q_t)), so reactivation can condition on *realized* Q_t rather than on E_{t−1}Q_t.
This is what makes Λ_r,t a function of current conditions and therefore procyclical — the
property the whole margin rests on. It also confirms the `Λ_r,t` dating in f[22] and its
pairing with (1 − δ_e,t).

**Two things simultaneity does NOT mean.**

- **Not the same agent, and not the same draw.** Exit is the *retailer's* decision on its
  product line; reactivation is the *recruiter's* decision on a vacated position, since
  recruiters are the ones holding vacancies. Two distinct draws against two distinct
  thresholds: the retailer compares its χ to x_c,t, the recruiter compares α·χ to Q_t.
  Retailer exit still gates reactivation, because a position serving a withdrawn product line
  is destroyed regardless — that is the (1 − δ_e,t) factor in f[22].
- **Not simultaneous determination.** Within Stage 2 exit resolves first, and only recruiters
  attached to surviving product lines face reactivation decisions. The reposting option
  accrues to the recruiter's J, not the retailer's ν_f, so the dependence stays triangular and
  f[1]/f[4]/f[18] remain unchanged (R7). The segmentation, not the timing, is what delivers
  that.

⚠️ **Open: the per-position-draw justification in Block 2 needs restating.** It argues that a
single draw per *firm* would make firms carry different position stocks forward and break the
DS-CES symmetric aggregation. That argument attributes positions to retailers. Under the
segmentation positions belong to *recruiters*, and retailers rent labor services competitively,
so the distribution of vacated positions across recruiters does not reach retailers at all:
only the aggregate vacancy stock enters matching, and the law of large numbers delivers Λ_r
either way. **If so, the symmetry concern dissolves and per-position draws are simply the
natural reading — each position is a separate asset with its own reactivation cost — rather
than a requirement for aggregation.** Worth confirming before the argument is relied on in
print; the draft's §Environment now states the aggregation point in the segmentation form.

## Shocks (ne = 3)

| Shock | Meaning | Process | Innovation SD |
|---|---|---|---|
| z | technology | AR(1) in logs, ρ_z | σ_z |
| δ | permanent firm / product-line destruction | AR(1) in logs, ρ_δ | σ_δ |
| s | idiosyncratic worker separation | AR(1) in logs, ρ_s | σ_s |

Innovations are orthogonal by construction (Cholesky, z ordered first, `part6b`). The z→δ_e
link runs through the equilibrium (x_c responds to z), not through shock covariance.

## Variables

**States** x = [u, N, v_pret, z, δ, s] (predetermined, known at the start of *t*):

| Symbol | Meaning |
|---|---|
| u | unemployment entering *t* (pre-matching) |
| N | mass of incumbent firms/varieties entering *t* |
| v_pret | surviving vacancies from *t−1* (pre-entry); total v_t = v_pret + e_t |
| z, δ, s | exogenous shock levels (SS = 1) |

**Controls** y (27): θ, q, L, v, e, K, Q, ρ, N_e, ν_f, d_f, w_int, w, L_e, L_c, Y_c, C, λ,
Y, x_c, δ_e, labor_prod, C_R, Y_R, Y_cR, w_R, ls.

| Symbol | Meaning |
|---|---|
| θ, q | market tightness v/u; vacancy-filling rate A·θ^(−η_L) |
| L, L_c, L_e | total employment; production workers; recruiters |
| v, e | total vacancies; new entrant vacancies |
| K, Q | net value of a vacancy; sunk cost of posting a vacancy |
| ρ | DS-CES relative price N^(1/(ε−1)) |
| N_e | new firms (entrants) |
| ν_f, d_f | firm value; firm dividend |
| w_int, w | marginal revenue product (recruiter compensation base); Nash wage |
| Y_c, Y, C | retail output; GDP; consumption |
| λ | marginal utility of consumption C^(−σ) |
| x_c | continuation-cost cutoff (firms with cost > x_c exit) |
| δ_e | endogenous destruction rate 1 − (1−δ·δbar)·F(x_c) |
| labor_prod, C_R, Y_R, Y_cR, w_R | data-consistent (price-deflated) observables |
| ls | labor share w·L/Y |

**Continuation-cost distribution** (heterogeneous exit friction), the model's one genuine
cost distribution besides the entry-cost convexity:

```
F(x) = (1 − p_0) + p_0·(x/f_m)^ψ        for x ∈ [0, f_m]
F(x) = 1                                 for x > f_m
```

Mass 1 − p_0 at zero cost (these firms never endogenously exit); a power-law tail with shape
ψ up to f_m. Λ ≡ F(x_c) is the survival probability. ψ_c ≡ ψ/(ψ+1).

Survival factor for *t*→*t+1*: `SDF_surv = (1 − δ·δbar)·Λ'`, using the *predetermined* current
δ and the *next-period* survival Λ' = F(x_c').

## Equations f[1]–f[33]

### Endogenous exit — f[1]–f[2]

**f[1] Exit threshold** (`eq:cutoff_eq` / `eq:cutoff_free_entry`)
```
x_c = Y_c·(μ−1)/(μ·N) + ν_f
```
The cost cutoff below which a firm continues equals current per-firm gross profit R^f =
Y_c·(μ−1)/(μN) plus the firm's continuation value ν_f. A firm draws a continuation cost each
period and continues iff the draw is below this threshold. The cutoff rises with profitability
and with the value of a product line, so in good times marginal firms that would otherwise
exit stay.

**f[2] Endogenous destruction rate** (`eq:delta_e_lom`)
```
δ_e = 1 − (1 − δ·δbar)·Λ,     Λ = (1 − p_0) + p_0·(x_c/f_m)^ψ
```
The fraction of firms/jobs destroyed this period. A firm is destroyed either by the exogenous
shock (rate δ·δbar) or by drawing a continuation cost above the cutoff (prob 1 − Λ). δ is
predetermined; Λ responds to the cutoff, which is where technology and demand feed the exit
margin.

### Euler / asset-pricing equations — f[3]–f[5]

**f[3] Job creation condition** (`eq:jcc_eq`)
```
κ + K/q = β·(λ'/λ)·SDF_surv·[ (1−ϕ)(w_int'−K'−b) − ϕ·θ'(K'+q'κ) + (1−s'·sbar)(κ+K'/q') ]
```
The marginal cost of filling a job, κ + K/q (flow matching cost plus the vacancy value spread
over the filling probability), equals its discounted expected return: the firm's share of match
surplus, minus the worker's outside-option term, plus the re-hiring cost saved when the match
survives separation (prob 1 − s'). Survival-weighted by SDF_surv and discounted by the
stochastic discount factor β·λ'/λ. **This is the equation the partial-reposting extension
changes most (Block 2).**

**f[4] Business-formation Euler** (`eq:firm_value_char`, BGM)
```
ν_f = β·(λ'/λ)·SDF_surv·(ν_f' + d_f')
```
The value of a firm equals the discounted, survival-weighted sum of next period's dividend
d_f' and continuation value ν_f'. With free entry (f[18]) this makes the firm mass N
forward-looking: entry occurs until firm value equals the sunk entry cost.

**f[5] Vacancy (capital) value** (`eq:Kdef`)
```
K = Q − β·(λ'/λ)·SDF_surv·Q'
```
The net value of a posted vacancy is its sunk posting cost Q minus the discounted,
survival-weighted resale/continuation value of an unfilled vacancy next period. Q is what the
firm sinks; K is what remains after accounting for survival.

### Wage block — f[6]–f[8]

**f[6] Marginal revenue product** (`eq:recruiter_compensation`)
```
w_int = ρ·z·zbar/μ
```
The interior value of a worker's output, DS-CES: relative price ρ times productivity z·zbar,
divided by the markup μ.

**f[7] Nash bargaining wage**
```
w = ϕ·(w_int − K + θ(K + q·κ)) + (1−ϕ)·b
```
The wage splits match surplus: worker share ϕ of the joint value (output net of the firm's
capital, plus the tightness-weighted hiring-cost term) plus (1−ϕ) of the outside option b.

**f[8] Vacancy creation** (`eq:entry_eq`)
```
Q = x_m·e^(ξ_inv)
```
The marginal cost of posting the e_t-th vacancy. Convex in the entry flow e (exponent ξ_inv),
scaled by x_m. This convexity is the reduced form of a *distribution* of sunk vacancy-creation
costs; only its shape parameter ξ_inv enters. This is the position-level *creation* cost, the
sibling the reposting margin does **not** reuse (see Block 2 and note §7.3).

### Labor market — f[9]–f[12]

**f[9] Tightness** `θ = v/u` — vacancies per unemployed worker, using total v.

**f[10] Vacancy-filling rate** `q = A·θ^(−η_L)` — Cobb-Douglas matching; a tighter market
fills vacancies more slowly.

**f[11] Employment** `L = 1 − u` — labor-market clearing.

**f[12] Labor split** `L = L_c + L_e` — employment divides into production workers L_c and
recruiters L_e.

### Goods market — f[13]–f[15]

**f[13] Relative price** `ρ = N^(1/(ε−1))` — DS-CES: more varieties N raise the relative price
of each, the variety channel that drives amplification.

**f[14] Household Euler** `λ = C^(−σ)` — marginal utility of consumption.

**f[15] Retail production** `Y_c = ρ·z·zbar·L_c` — retail output from production labor at
productivity z·zbar and relative price ρ.

### Resource constraint, identities, laws of motion — f[16]–f[24]

**f[16] Resource constraint** (`eq:rc`)
```
Y_c = C + X + X_c
X   = e/(1+ξ_inv)·Q + κ·q·v          (sunk posting costs + matching costs)
X_c = N·p_0·ψ_c·x_c                    (aggregate continuation costs)
```
Retail output funds consumption, vacancy investment, and continuation costs. X splits into the
integral of the marginal posting-cost schedule, ∫₀^e x_m·u^ξ_inv du = e·Q/(1+ξ_inv), and the
fixed matching cost κ per match (q·v matches, Pissarides 2009). X_c is the total continuation
cost paid by incumbents.

**f[17] New entrants** `N_e = z·zbar·L_e/f_e` — entrant firms produced by recruiter labor L_e
at sunk entry cost f_e in labor units.

**f[18] Free entry** (`eq:free_entry`) `ν_f = ρ·f_e/μ` — entry until firm value equals the sunk
entry cost expressed in consumption units.

**f[19] GDP** (`eq:gdp`) `Y = C + ν_f·N_e` — output is consumption plus the value of newly
created firms.

**f[20] Income identity** `N·d_f = Y_c·(μ−1)/μ − X_c` — aggregate dividends equal retail
profits net of continuation costs.

**f[21] Total vacancies** (`eq:v_lom` split) `v = v_pret + e` — surviving pre-entry vacancies
plus new entrant vacancies; both participate in matching this period.

**f[22] Pre-entry vacancy law of motion** (`eq:v_lom`) — **modified in Block 2**
```
v_pret' = (1 − δ_e)·[ (1 − q)·v + s·sbar·(1 − u) ]
```
Vacancies carried into *t+1* are surviving firms' (factor 1 − δ_e) unmatched vacancies
(1 − q)v plus **reposted separations** s·sbar·(1 − u). Under full reposting *every* separation
returns as a vacancy. This is the term the Beveridge-curve problem lives in: the reposting
inflow (1 − δ_e)·s·sbar·(1 − u) > 0 means an s shock raises v (Prop. 5 Part 1, Channel 1).

**f[23] Unemployment law of motion** (`eq:u_lom`, Convention A pre-matching)
```
τ_t = δ_e + s·sbar·(1 − δ_e)
u'  = [1 − (1 − δ_e)·f(θ)]·u + τ_t·(1 − u),     f(θ) = A·θ^(1−η_L)
```
Unemployment next period is job-losers who did not find work plus newly separated workers.
Matches formed at *t* reach *t+1* only if the firm survives, hence (1 − δ_e) on the job-finding
outflow. τ_t applies to the inherited stock (1 − u) only, so a match cannot dissolve before it
produces. This base is what f[22] reposts on, so v_pre' + L' = (1 − δ_e)(v + L): positions are
conserved exactly.

**f[24] Firm law of motion** (`eq:N_lom_eq`) `N' = (1 − δ_e)·(N + N_e)` — incumbents plus
entrants, both facing current-period destruction δ_e. This is the *variety* margin.

### Labor share — f[25]

**f[25]** `ls = w·L/Y`.

### Data-consistent observables — f[26]–f[30]

**f[26]** `labor_prod = Y/(ρ·L)` — real output per worker.
**f[27]–f[30]** `C_R = C/ρ`, `Y_R = Y/ρ`, `Y_cR = Y_c/ρ`, `w_R = w/ρ` — price-deflated
aggregates, so model series match empirically deflated data.

### Exogenous shock processes — f[31]–f[33]

**f[31]** `log z' = ρ_z·log z` — technology.
**f[32]** `log δ' = ρ_δ·log δ` — structural exit (⊥ z).
**f[33]** `log s' = ρ_s·log s` — worker separation.

---

# BLOCK 2 — Partial-reposting extension

**What changes and why.** Under full reposting every match separation returns as a vacancy, so
the s shock mechanically raises v (Block 1, f[22]; Prop. 5 Part 1). To fit the Beveridge curve
at reasonable δ, some separations must fail to repost (note §4.3a: this is inescapable). On
separation, each vacated position independently draws a reposting cost from the **same**
continuation-cost distribution F, scaled by α. This is a **one-time reactivation decision**:
reactivate the slot for search (pay α·χ, obtain a durable vacancy) or retire it permanently
(headcount cut; to have it back later, create a fresh position at full cost Q). Reactivate iff
the value of the vacancy obtained covers the cost:

```
cost of reactivating a vacated position = α·χ,   χ ~ F,   α > 0
reactivate  ⟺  Q ≥ α·χ  ⟺  χ ≤ Q/α
```

**Threshold is Q** (the value of an unfilled vacancy). Confirmed by the R7 derivation (f[3]
below) against the draft's recruiter block: on separation the recruiter obtains an unfilled
vacancy worth Q_{t+1} (`eq:value_recruiter_matched`), which the Bellman already discounts, so
reactivation compares α·χ to Q directly. *(An earlier draft of this file used Q − K; that was an
over-refinement — the exact object is Q. K ≈ 0.6% of Q, so the numerical difference was
negligible.)* Q is the only threshold that yields a non-degenerate cyclical Λ_r; a flow-value K
threshold would pin Λ_r at its floor 1 − p_0. The alternative that would justify a K threshold —
reposting as a *recurring per-period active-search maintenance* decision — was considered and
rejected: it models temporary withdrawal rather than the permanent position destruction the
missing middle requires, needs a new dormant-position state, and delivers less persistence. See
note §7.3b for the two stories and the verdict.

**Per-position draws, not one draw per firm.** The reposting cost is drawn independently for
each of a firm's many positions, not once per firm. This is deliberate and preserves the
model's single firm size. Broer's iid continuation cost keeps all firms identical at production
because *exit is terminal* — survivors reset each period, leaving no trace. Reposting instead
lands on a *continuing* state (the firm's position stock), so a single per-firm draw would make
firms carry different stocks forward, and firm-size heterogeneity would accumulate and break the
DS-CES symmetric aggregation (`ρ = N^(1/(ε−1))`). With per-position draws, the law of large
numbers makes **every** surviving firm repost the identical fraction Λ_r = F(Q/α), so all firms
shrink together and stay symmetric replicas of the aggregate — one firm size, exactly as Broer
intended, now at the position level. (An earlier draft used one draw per firm with cost α·χ and
a three-region "missing middle" sorted by χ; that version breaks firm symmetry and is
superseded. The missing middle survives here, spread uniformly across firms.)

**New parameter:** α (reposting-cost scale), estimated. **Distributions unchanged:** F keeps
(p_0, ψ, f_m); reposting adds no new distribution. Rationale for drawing from F is in note
§7.3a.

**Notation.** The reposting rate is written **Λ_r**, capital, parallel to the survival
probability Λ = F(x_c) — both are F(·) evaluated at a threshold. It is **not** λ: lowercase λ is
already the marginal utility of consumption (f[14]), which appears in every Euler equation as
the SDF ratio λ'/λ. Capital Λ and lowercase λ coexist in the current code; Λ_r extends that.

## New equation

**f[new] Reposting rate** — give Λ_r its own defining equation (a tracked control, like δ_e and
x_c, so its impulse response can be reported; or substitute inline as Λ is)
```
Λ_r = F(Q/α) = min[ (1 − p_0) + p_0·(Q/(α·f_m))^ψ , 1 ]
```
The fraction of a surviving firm's vacated positions that are reactivated (threshold Q/α, derived
in f[3]). Procyclical: when Q rises in booms, Q/α rises, Λ_r rises, so position destruction is
countercyclical. Bounds: Λ_r ∈ [1 − p_0, 1]. The atom 1 − p_0 always reactivates (zero-cost
positions), so the most that can ever be destroyed is fraction p_0 of separations — **p_0 caps
the reposting margin's strength.** Λ_r *saturates* at 1 (all vacated positions repost) for
α ≤ Q̄/f_m; the margin bites (Λ_r < 1) only for α > Q̄/f_m. Note this is quantity saturation, not
the costless Block-1 benchmark: at Λ_r = 1 with α > 0 and p_0 > 0 all positions repost but
still pay α·χ.

**Two routes to the costless benchmark, not one.** An earlier version of this line said the
costless model is recovered "only as α → 0." That is incomplete. **p_0 → 0 also recovers it**,
and by a different mechanism: F collapses to a point mass at zero, so every reactivation cost
draw is α·χ = 0. Then Λ_r = 1, α·M = 0, Q^rep = Q, and X_r = 0. The α → 0 route shrinks the
cost scale; the p_0 → 0 route removes the cost mass.

The p_0 route is the one that matters for **Proposition 1** (`prop:independence`). Its third
condition is p_0 = 0, which removes endogenous exit — and because reposting draws from the
*same* F, it removes the reposting margin at the same time. So Proposition 1 holds with its
three original conditions and needs no fourth. Had reactivation costs come from a separate
distribution, p_0 = 0 would not have shut reposting down and the proposition would have
required an extra condition. This is a second payoff of reusing F, beyond the
one-fewer-parameter argument in note §7.4.

## Modified equations

**f[22] Pre-entry vacancy law of motion — MODIFIED**
```
v_pret' = (1 − δ_e)·[ (1 − q)·v + Λ_r'·s·sbar·(1 − u) ],     Λ_r' = F(Q'/α)
```

⚠️ **Λ_r is PRIMED here — corrected Sept 30, 2026.** An earlier version wrote `Λ_r` unprimed,
which dates the reposting rate one period before the vacancy it values. The position separates
at end of *t* (hence `s`, `(1−u)` unprimed), and the vacancy obtained by reactivating it is a
*t+1* vacancy worth `Q'`. So the threshold is `Q'/α` and the rate is `Λ_r' = F(Q'/α)`. This
matches the recruiter Bellman in f[3] below, where a period-*t* separation is evaluated with
`Λ_r,t+1`, and the draft's `eq:Qrep`. In draft notation (with `v_t` on the left) the same rate
is written `Λ_{r,t}` — the rate always carries the period of the vacancy it values, which is
the LHS stock's period in both notations.

⚠️ **Open: the δ_e dating convention differs between draft and code.** The draft writes
`v_t = (1−δ_{e,t})[…]`, dating survival at the LHS period; this code equation writes
`(1 − δ_e)` unprimed with `v_pret'` on the left, dating it at the RHS period. Under the draft's
convention `Λ_{r,t}` and `(1−δ_{e,t})` pair at the same date; under the code's they do not
(`Λ_r'` vs `δ_e`). One of the two conventions is wrong, or they differ deliberately. Resolve
before R8 — see [`../Notes/LOM_timing_consistency.md`](../Notes/LOM_timing_consistency.md),
which settled the u-LOM convention and is the right place to record this.

Otherwise identical to the baseline f[22] with the reposting rate inserted on the separation
inflow.
Under full reposting Λ_r = 1 and it collapses to Block 1. The non-reposted fraction (1 − Λ_r) of
separations is the missing middle: those positions are destroyed. This is the single change that
lets an s shock lower v, delivering the Beveridge curve. *Scope choice:* attenuation is applied
to the separation inflow only, not to the unmatched-vacancy carry-forward (1 − q)v; a firm's
already-posted unmatched vacancies are sunk this period. Revisit if the maintenance decision
should also cover them.

**Why the LOM needs only δ_e for destruction.** The "survive and repost" event is a sub-event
of survival: a position is reposted only if its firm continues, and continuing has probability
exactly 1 − δ_e. So the reposting flow always sits inside the survival flow, and one can pull
(1 − δ_e) out front, leaving the reposting rate Λ_r inside the bracket. Writing the weight as
(1 − δ·δbar)·F(Q/α) would spell out "survive exogenous exit **and** clear the reposting cutoff"
as one joint probability; factoring uses (1 − δ_e) = (1 − δ·δbar)·Λ to separate the two. The
factoring is legitimate precisely because reposting ⊂ survival — if an *exiting* firm could
repost, the weight would not nest inside 1 − δ_e and the two destruction events would have to be
written separately.

**f[3] Job creation condition — MODIFIED (derived R7, Sept 29)**

Derivation grounds in the draft's recruiter block (`eq:value_recruiter_unmatched`,
`eq:value_recruiter_matched`, `eq:jcc`). The **only** value function that changes is the matched
recruiter's, and only its *separation branch*. In the draft, on separation (prob s_t) the
recruiter obtains an unfilled vacancy worth Q_{t+1} — reposted for free. Under partial reposting
that becomes the **reactivation option value**

```
Q^rep_{t+1} = E_χ[ max(Q_{t+1} − α·χ, 0) ] = Q_{t+1}·Λ_r,t+1 − α·M_{t+1}
Λ_r = F(Q/α),   M = ∫_0^{Q/α} χ dF(χ);   interior (Λ_r<1): α·M = ψ_c·Q·(Λ_r−1+p_0), ψ_c = ψ/(ψ+1)
```

**Notation (aligned with `Draft.tex` `eq:Qrep`).** The draft writes this option value as **Q^rep**,
renamed from an earlier `Ψ` that clashed visually with the shape index `ψ`. In the *draft* the mean
integral is written inline (not given a symbol — a capital `M` would collide with the matching
function `m`/the SDF), and the shortfall is written inline as `Q − Q^rep`, *not* `D` (plain `D`
collides with the intermediary dividend `D^{int}` and the duration multiplier `\mathcal{D}`). This
doc uses `M ≡ ∫_0^{Q/α} χ dF` as shorthand only. The closed form
`α·M = ψ_c·Q·(Λ_r−1+p_0)` holds only in the **interior** regime `Q/α ≤ χ_m` (`Λ_r < 1`); in the
saturated regime it is `α·∫_0^{χ_m} χ dF`.

### Deriving M, and what its factors mean

Worth writing out once, because three different "per-position cost" objects live here and
conflating them has already caused one error (see f[16]). With threshold `x = Q/α`, the atom
at zero contributes nothing, so only the power-law part integrates:

```
dF = p_0*ψ*χ^(ψ-1)/χ_m^ψ dχ        on (0, χ_m]
M  = ∫_0^x χ dF = (p_0*ψ/χ_m^ψ) * ∫_0^x χ^ψ dχ
   = p_0 * ψ/(ψ+1) * x^(ψ+1)/χ_m^ψ
   = ψ_c * x * p_0*(x/χ_m)^ψ
   = ψ_c * x * (Λ_r - 1 + p_0)          since Λ_r = F(x) = (1-p_0) + p_0*(x/χ_m)^ψ
α*M = ψ_c * Q * (Λ_r - 1 + p_0)     using x = Q/α
```

**α cancels.** It scales the cost and the threshold in opposite directions, so α enters
α*M only through Λ_r.

**The two factors.**

| Factor | Meaning |
|---|---|
| `ψ_c*Q = α*ψ_c*x` | mean cost among positions that pay a **positive** cost. Independent of p_0 |
| `(Λ_r - 1 + p_0)` | probability a vacated position is a *paying* reposter: draws from the power law **and** clears the threshold |

So `α*M` = (probability of paying) × (mean payment) = expected cost **per vacated position**.

**Three distinct objects — do not substitute one for another:**

| Object | Value | Per what |
|---|---|---|
| `α*M` | `ψ_c*Q*(Λ_r-1+p_0)` | per **vacated** position — the one f[16] needs |
| `α*M/Λ_r` | `ψ_c*Q*(Λ_r-1+p_0)/Λ_r` | per **reactivated** position (includes the free atom) |
| `ψ_c*Q` | `ψ_c*Q` | per position that actually **pays** |

**Limits.** As `Λ_r → 1 - p_0` the threshold goes to zero, no power-law draw clears, and
`α*M → 0`: only free positions repost, so no costs are paid. At `Λ_r = 1`,
`α*M = ψ_c*Q*p_0`. At the saturation boundary `Q = α*χ_m` the interior form and the
saturated form `α*ψ_c*χ_m*p_0` agree, so α*M is continuous across regimes.


**Threshold is Q, not Q − K** (corrects an earlier draft of this file). The recruiter obtains an
unfilled vacancy worth Q_{t+1}; the draft's Bellman already discounts it by m(1−δ)F, so the
reactivation compares α·χ to Q directly. Reactivate iff **Q ≥ α·χ ⟺ χ ≤ Q/α**.

**Matched recruiter value** (modified `eq:value_recruiter_matched`):
```
J_t = w_int − w + β·(λ'/λ)·(1 − δ)·F(x_c')·[ s·Q'^rep + (1 − s)·J' ]
```

**Recruiter surplus** (`eq:value_recruiter_surplus` gains one term):
```
J_t − Q_t = w_int − w − K + (1 − s)(κ + K/q)  −  s·β(λ'/λ)(1−δ)F(x_c')·(Q' − Q'^rep)
```
The shortfall `Q' − Q'^rep = Q'(1 − Λ_r') + α·M'` is the per-separated-position value lost relative
to free full reposting: `α·M'` = reactivation costs paid by reposters, plus `Q'(1−Λ_r')` = vacancy
value forgone by non-reposters.

**Modified JCC.** Substituting the surplus into the unchanged average-hiring-cost relation
`κ + K/q = E m(1−δ)F(J' − Q')`:
```
κ + K/q = β(λ'/λ)(1−δ)F(x_c')·[ w_int'−w'−K' + (1−s')(κ+K'/q') − s'·β(λ''/λ')(1−δ')F(x_c'')·(Q''−Q''^rep) ]
```
The new term `− s'·(discounted (Q''−Q''^rep))` is the extra penalty: a job created now may separate
at t+1 (prob s'), and that separation's shortfall `Q−Q^rep` is realized on the vacancy obtained at
t+2, so it enters with a **nested** expectation. **Implementation note:** rather than substitute,
keep the matched value J (or the surplus S = J − Q) as a tracked jump variable with its own
Bellman (the modified `eq:value_recruiter_matched` above); this keeps the system first-order
Markov and avoids the nested E. **Nesting:** the *costless* current model is recovered as **α → 0**
(then Q^rep → Q and the shortfall → 0); note this differs from the f[22] LOM, whose *quantity* of
reposted vacancies saturates at Λ_r = 1. At Λ_r = 1 with α > 0 all vacancies repost but costs are
still paid (Q^rep < Q), so only α → 0 gives the fully costless benchmark.

**f[5] Vacancy (capital) value — UNCHANGED in form.** The unfilled-vacancy Bellman and K =
Q − β(λ'/λ)(1−δ)F(x_c')Q' do not involve the reposting decision (reposting acts on the *matched*
value's separation branch). K inherits new equilibrium values through Q and the surplus, but its
equation is untouched. Confirmed by the derivation.

**f[16] Resource constraint — MODIFIED (add reactivation costs actually paid)**
```
Y_c = C + X + X_c + X_r
X_r = (1 − δ_e)·s_{−1}·sbar·(1 − u_{−1}) · α·M,    α·M = ψ_c·Q·(Λ_r − 1 + p_0)
```

⚠️ **Timing follows from the f[22] correction above — settled Sept 30, 2026.** The cost is paid
when the position is reactivated, and reactivation happens at the period of the vacancy
obtained. So a period-*t* resource constraint carries costs on positions vacated at end of
*t−1*, valued at the *t* threshold `Q_t`: hence the **lagged** separation flow with the
**current** `α·M`. This pairs term by term with `eq:v_lom`'s reposting inflow, which in draft
notation is `(1−δ_{e,t})·Λ_{r,t}·s_{t−1}(1−u_{t−1})` — the same count, with `α·E[χ|repost]`
per position instead of one vacancy per position. In draft notation:
`X_{r,t} = (1−δ_{e,t})·s_{t−1}(1−u_{t−1})·α∫_0^{Q_t/α} χ dF`.

**R8 implementation note.** The draft's LOM already carries lagged separation flows, so this is
natural there. The code's f[22] uses *unprimed* `s·sbar·(1−u)` for `v_pret'`, so a t-dated X_r
needs a lagged flow the code does not currently track. Two options: (a) add the lag, matching
the draft; or (b) re-date to *pay at separation* — the cost is sunk at end of *t* against the
anticipated `Q'`, giving `X_r = (1 − δ_e)·s·sbar·(1 − u)·α·M'` with no new state. Option (b) is
cheaper and arguably the more natural reading of a sunk reactivation decision, but it is a
different model, not a notational variant. **Pick one and record it here.**
**α·M is expected reactivation cost per VACATED position, not per reposter.** `M` is a *partial*
(unconditional) expectation, so `α·M = α·Λ_r·E[χ | χ ≤ Q/α]` — it already contains Λ_r. That is
what makes it the right multiplicand for the count of vacated positions, and it is why **Λ_r
must not appear again** in X_r. Full derivation and the three distinct per-position objects:
§"Deriving M, and what its factors mean" under f[3] above. The conditional mean among reposters
is `α·M/Λ_r`; the two coincide only at Λ_r = 1, which is why conflating them is easy.

⚠️ **Do NOT apply the draft's X_c shortcut to X_r — corrected Sept 30, 2026.** The draft writes
`X^c = N·p_0·ψ_c·χ^c` (`eq:agg_fixed_costs`), dropping the `(χ^c/χ_m)^ψ` factor the exact
integral carries. Substituting `p_0` for `(Λ − 1 + p_0)` costs a relative error of
`p_0/(Λ − 1 + p_0)`:

- **For X_c it is negligible.** From f[2], Λ = (1−δ_e)/(1−δ·δbar) = 0.99729/0.998645 ≈ **0.9986**
  at `dest_ann = 0.0320`, so with p_0 = 0.5 the error is **0.27 %**.
- **For X_r it is not.** Λ_r < 1 by construction — that is the entire point of the margin. At
  Λ_r = 0.9 the shortcut overstates X_r by 25 %; at Λ_r = 0.7, by **2.5×**.
- **At the floor Λ_r = 1 − p_0 it is qualitatively wrong**: exact X_r = 0 (only free positions
  repost) while the shortcut gives ψ_c·Q·p_0 > 0.

The approximation is good for X_c *because* survival is near 1, and bad for X_r *because*
reposting is not. Use the exact integral for X_r, as written above. There is no reason to
approximate either: the exact forms are `α·M = ψ_c·Q·(Λ_r − 1 + p_0)` and
`X_c = N·ψ_c·x_c·p_0·(x_c/f_m)^ψ`, no harder to write than the shortcuts. Upgrading X_c costs
0.27 % and removes the inconsistency.

**f[1] / f[4] / f[18] Exit threshold, firm value, free entry — UNCHANGED in form (confirmed Sept 29).**
Checked against the draft's retailer block. The firm Bellman (`eq:firm_bellman`) contains only
retail operating profit R^f = (μ−1)/μ·ρy (`eq:retail_profits`), the continuation-cost draw χ, and
the discounted continuation value — no Q, K, J, or reposting object. Hence the cutoff
x_c = R^f + ν_f (f[1]), the dividend d_f = Π^f/N with Π^f = Y_c/ε − X_c, and the
business-formation Euler (f[4]) carry no reposting term. Recruiter profits Π^int = (w_int − w)L
are a *separate* flow to the household (`Y = wL + Π^f + Π^int`), not part of ν_f or d_f. The
recruiter/retailer segmentation makes the dependence triangular: the recruiter's Q, J, and the
reposting option depend on retailer survival (1−δ)F(x_c), but the retailer's x_c, ν_f, d_f do
**not** depend on reposting. Reposting is downstream of the exit decision. So the reposting
option lives entirely in the recruiter's J → f[3]; f[1]/f[4]/f[18] move only through general
equilibrium (aggregate Y_c, w_int, N → R^f), with no new term.

## Equations that do NOT change, and why

- **f[2] δ_e (destruction rate):** firm exit is unaffected. Reposting destroys *positions* at
  *surviving* firms, not firms/varieties.
- **f[23] Unemployment LOM:** workers separate into unemployment at τ_t regardless of whether
  the firm reposts the vacated position. Worker flows are unchanged; only the vacancy return
  changes (f[22]). Note the position-conservation identity v_pre' + L' = (1 − δ_e)(v + L) now
  breaks by exactly the non-reposted mass — that leakage *is* the missing middle.
- **f[24] Firm LOM (N):** the variety margin is untouched. This is the core asymmetry the paper
  keeps: only δ moves N, so only δ produces a persistent outward Beveridge shift; a non-reposted
  s separation destroys a position but leaves the product line and the entry margin intact.
- **f[6]–f[21], f[25]–f[33]:** wage block, goods market, entry, observables, and shock processes
  are unchanged in form (they inherit new equilibrium values through Q, K, and the LOMs).

## Open derivation items (see note §9 and §10)

1. Re-derive f[3] (JCC) and confirm f[5], f[1], f[4] with the reposting option; restate
   `lem:vpre` and the Prop. 5 Part 2 Jacobian, which touch the reposting channel.
2. ✅ **Resolved Sept 30, 2026.** X_r uses the **exact** integral, not the X_c shortcut — see
   f[16]. Remaining sub-item: choose the X_r timing convention (lagged flow vs pay-at-separation)
   and record it. The δ_e dating convention between draft and code also needs settling (f[22]).
3. Derive the threshold Λ_r* where the s→v sign flips; Prop. 5 Part 1 becomes conditional on Λ_r.
4. Verify on simulated data that α is separately identified from σ_s and from f_m (note §7.5).
