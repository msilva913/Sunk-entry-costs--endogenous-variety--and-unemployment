# The AGS model with costly reactivation: equilibrium conditions
**Written:** November 4, 2026 · **For:** Part 1 of the two-part plan (decision D15)
**Status:** a derivation for discussion. Code: [`../Programs AGS/`](../Programs%20AGS/)
(`ags_steady_state.jl`, `ags_steady_state_checks.jl`, `ags_run_solution.jl`). Corrected Nov 4:
the state vector (§3) and the vacancy free entry special case. The AGS restrictions are imposed on the equations in
[`../context/model_equations.md`](../context/model_equations.md) (Block 1 plus the Block 2
reposting margin), and reactivation is given its own cost distribution. Check the notation
against the draft's `prop:ags` before any of this enters the paper.

---

## 1. The model in one paragraph

Start from DMP. Add three things.
1. **Positions are durable assets** (Coles and Kelishomi; Gabrovski and Silva). Creating a position costs Q, and an unfilled position carries over. The flow cost of holding it is a user cost K, not an exogenous c.
2. **Destruction δ removes jobs *and* vacancies.** Product lines are destroyed at rate δ, and with them every position attached, filled or vacant. Separation s ends a match but leaves the position.
3. **Reactivation is costly.** A position vacated by a separation is reactivated only if its reactivation cost is below the value of a vacancy. Otherwise it is retired. With costless reactivation, the model is Gabrovski–Silva with aggregate shocks, and with s ≡ 0 it is Coles–Kelishomi.

Households have CRRA utility (σ = 1 by default, as in the full model) and own the firms, so jobs and positions are valued with the household's discount factor m_{t+1} = β(C_{t+1}/C_t)^{−σ} (§2, "The household and its savings decision"). σ = 0 is linear utility and constant discounting, the Gabrovski–Silva benchmark. Destroyed product lines are replaced by clones at no cost ("clone replacement"), so there is no business-formation block. The shocks are productivity z, destruction δ and separation s.

## 2. Environment

**Timing within a period (monthly):**
1. z_t, δ_t and s_t are realized and known.
2. Positions vacated at the end of t − 1 are reactivated or retired. New positions e_t are created.
3. Production; u_t unemployed workers and v_t vacancies match.
4. At the end of the period, separations s_t and destruction δ_t take effect.

**Matching:** a constant-returns matching function. With θ = v/u, the job-filling rate is q(θ) = Aθ^{−η} and the job-finding rate is f(θ) = θq(θ).

**Positions:**
- **Creation:** creating e_t new positions costs Q(e_t) = x_m e_t^{1/ξ} at the margin. Heterogeneous creation costs give an upward-sloping supply of positions, with elasticity ξ. Vacancy free entry is the limit ξ → ∞, where Q = x_m.
- **Hiring:** there is a per-hire cost κ ≥ 0, which is zero in Gabrovski–Silva.

**Reactivation:** a position vacated by a separation at the end of t − 1 draws a reactivation cost c at t, i.i.d. across positions and periods, from F_r. It is reactivated iff c ≤ Q_t (paying c and keeping a position worth Q_t beats retiring it), and retired otherwise.

*The distribution F_r.* A mixture of an atom at zero and a bounded power law:
- with probability 1 − p_r, c = 0 (the position is reactivated free);
- with probability p_r, c is drawn from a power law on [0, a]: Pr(c ≤ x | costly) = (x/a)^{ψ_r} for 0 ≤ x ≤ a, with density ψ_r x^{ψ_r − 1}/a^{ψ_r}.

So the CDF is F_r(x) = (1 − p_r) + p_r (x/a)^{ψ_r} for 0 ≤ x ≤ a, and 1 for x > a. The three parameters:
- **p_r ∈ [0, 1]**, the share of positions with a positive cost; p_r = 0 is costless reactivation;
- **a > 0**, the largest possible cost, which sets the scale; calibrated from Λ̄_r;
- **ψ_r > 0**, the shape. ψ_r = 1 is uniform on [0, a]; ψ_r < 1 puts more mass near zero; ψ_r > 1 puts more mass near a. ψ_r is also the elasticity of Λ_r − (1 − p_r) with respect to Q, so it governs how strongly reactivation responds to Q.

A useful property: conditional on being costly and below x ≤ a, the draw averages E[c | costly, c ≤ x] = ψ_{c,r}·x, with ψ_{c,r} = ψ_r/(1 + ψ_r) (one half when uniform).

*The reactivation rate* is the probability that c ≤ Q:

  Λ_r(Q) = F_r(Q) = (1 − p_r) + p_r (Q/a)^{ψ_r}  if Q ≤ a,
  Λ_r(Q) = 1                                       if Q > a  (every costly draw is accepted).

*The expected net value of a vacated position* is what it is worth, Q, if reactivated, less the cost paid:

  Q^rep(Q) = E max{Q − c, 0} = Q·Λ_r(Q) − E[c · 1{c ≤ Q}].

The atom costs nothing, so only accepted costly draws pay. Their probability is p_r·min{Q/a, 1}^{ψ_r}, and their average cost is ψ_{c,r}·min{Q, a}. Hence:
- **Q ≤ a (interior):** the accepted costly draws are a share Λ_r − (1 − p_r) and cost ψ_{c,r}Q on average, so

  Q^rep(Q) = Q·[ Λ_r(Q) − ψ_{c,r}(Λ_r(Q) − 1 + p_r) ];

- **Q > a (bound):** every position is reactivated (Λ_r = 1), and the costly ones cost ψ_{c,r}·a on average, so

  Q^rep(Q) = Q − p_r ψ_{c,r} a.

The two branches meet at Q = a, where both give a(1 − p_r ψ_{c,r}). The interior branch, rearranged by type of position:

  Q^rep/Q = (1 − p_r) · 1 + (Λ_r − 1 + p_r)(1 − ψ_{c,r}) + (1 − Λ_r) · 0:

free positions keep all of Q, costly positions that are reactivated keep Q net of their average cost, and retired positions keep nothing.

*Checks.*
- p_r = 0 gives Λ_r = 1 and Q^rep = Q.
- If no costly draw is accepted (Λ_r = 1 − p_r), Q^rep = (1 − p_r)Q.
- At the working point (Λ̄_r = 0.25, p_r = 0.8, ψ_r = 1, so Q̄/a = 1/16), Q^rep/Q = 0.20 + 0.05 × 0.5 = 0.225: a vacated position keeps 22.5 % of its value on average.

`ags_steady_state_checks.jl` verifies both branches against numerical integration for Q/a from 1/16 to 10. The first-order solution uses the interior branch, which is exact as long as Q stays below a (§5).

**Wages:** Nash bargaining, with worker share ϕ and flow value of unemployment b.

**The household and its savings decision** (added Nov 6, 2026; `ags_model.jl` f[19]–f[21]).
- **Members and income.** A representative household of measure one. Each member supplies one unit of labor inelastically. Employed members earn w_t; unemployed members have the non-market flow b, which sets the outside option in bargaining but is not a market resource (as in the full model). The household pools income and consumes C_t, with utility E_0 Σ_t β^t C_t^{1−σ}/(1 − σ).
- **Assets.** (i) Shares in a mutual fund that owns every firm, in fixed supply normalized to 1, with ex-dividend price V_t and dividends D_t; (ii) a one-period riskless bond in zero net supply with gross return R^f_t.
- **Budget:** C_t + V_t S_{t+1} + B_{t+1} = w_t(1 − u_t) + (V_t + D_t) S_t + R^f_{t−1} B_t.
- **First-order conditions,** with m_{t+1} ≡ β(C_{t+1}/C_t)^{−σ}:
  - shares: V_t = E_t[ m_{t+1}(V_{t+1} + D_{t+1}) ];
  - bond: 1 = R^f_t E_t[ m_{t+1} ].
- **Market clearing:** S = 1 and B = 0, so C_t = w_t(1 − u_t) + D_t. With dividends D_t = (z̄z_t − w_t)(1 − u_t) − X_t − X^r_t, this is the resource constraint
  C_t = z̄z_t(1 − u_t) − X_t − X^r_t,
  where X_t = e_t Q_t/(1 + ξ_inv) + κ q_t v_t is creating positions and hiring, and X^r_t = Ṽ^s_t(Q_t Λ_{r,t} − Q^rep_t) is reactivation costs.

**How to read the savings decision.**
1. **The household cannot change the economy's saving by trading.** Shares are in fixed supply and bonds net to zero. Aggregate saving happens inside firms, as investment in durable positions: creation (e_t at marginal cost Q_t), the per-hire cost κ and reactivation costs. These positions are the economy's capital.
2. **The household's Euler equations set the price of saving instead.** The discount factor m_{t+1}, equivalently the riskless rate R^f_t, adjusts until the household is willing to consume exactly what firms do not invest.
3. **Firms act for their shareholders,** valuing matches (E2) and positions (K, 𝓡) with m_{t+1}. So the household's willingness to shift consumption over time disciplines firms' investment.
4. **After a transitory shock that makes current resources scarce relative to the future** (output falls, or investment must rise), R^f rises and m falls. That lowers the value of a new match and of a vacant position, so firms cut creation and hiring. This frees resources and smooths consumption, at the cost of a deeper fall in vacancies and job finding. The more investment can be postponed (elastic creation), the smaller the rise in R^f.
5. **σ = 0:** the household is indifferent to the timing of consumption, R^f = 1/β at all dates, and firms' investment does not feed back on discounting (Gabrovski–Silva).
6. **The share price is determined but feeds back on nothing here:** with clone replacement there is no entry decision. In the full model the same equation, with free entry of product lines, becomes the entry Euler, which pins m_{t+1}(1 − δ_t): priced entry is a flexible investment margin that holds the discount factor in place (R19).

## 3. Equilibrium

Three objects collect the dynamics, all valued with the household's discount factor m_{t+1} = β(C_{t+1}/C_t)^{−σ}:
- the **user cost** of a vacancy, K_t ≡ Q_t − E_t[ m_{t+1}(1 − δ_t) Q_{t+1} ];
- the **hiring cost** per match, Γ_t ≡ κ + K_t / q(θ_t);
- the **reactivation shortfall**, 𝓡_t ≡ E_t[ m_{t+1}(1 − δ_t)( Q_{t+1} − Q^rep(Q_{t+1}) ) ], the expected value lost because reactivation isn't free.

(The first version of this note wrote C_t for the hiring cost; it is Γ_t now that C_t is consumption.)

**Definition (equilibrium).** Given u_0, Ṽ^u_0, Ṽ^s_0 and the shock processes for (z_t, δ_t, s_t), an equilibrium is a stochastic sequence {θ_t, e_t, C_t, u_{t+1}, Ṽ^u_{t+1}, Ṽ^s_{t+1}}, with Q_t = Q(e_t), such that for all t:

**(E1) Vacancies:**
  θ_t u_t = v_{pre,t} + e_t,   v_{pre,t} = Ṽ^u_t + Λ_r(Q_t) Ṽ^s_t.

**(E2) Job creation:**
  Γ_t = E_t{ m_{t+1}(1 − δ_t)[ (1 − ϕ)(z_{t+1} − K_{t+1} − b) − ϕ θ_{t+1} K_{t+1} − ϕ f(θ_{t+1}) κ + (1 − s_{t+1}) Γ_{t+1} − (1 − ϕ) s_{t+1} 𝓡_{t+1} ] }.

**(E3) Unemployment:**
  u_{t+1} = u_t − (1 − δ_t) f(θ_t) u_t + τ_t (1 − u_t),   τ_t = δ_t + (1 − δ_t) s_t.

**(E4) Carried-over positions:**
  Ṽ^u_{t+1} = (1 − δ_t)(1 − q(θ_t)) v_t   (unmatched vacancies),
  Ṽ^s_{t+1} = (1 − δ_t) s_t (1 − u_t)      (positions vacated by separation, reactivated at t+1 or retired).

**(E5) Resources:** C_t = z̄ z_t(1 − u_t) − X_t − X^r_t, with X_t and X^r_t as in §2.

K_t, Γ_t and 𝓡_t are the functions of (θ, e, C, δ) defined above. That's three state variables (u, Ṽ^u, Ṽ^s) and three jump variables (θ, e, C), against one state (u) and one jump (θ) in DMP. With σ = 0, C drops out of E2 and the user cost, and the system is E1–E4 with constant β. The share price V_t and the riskless rate R^f_t follow from the Euler equations in §2 and feed back on nothing.

*Why v_pre is not a state.* Inherited vacancies v_{pre,t} = Ṽ^u_t + Λ_r(Q_t)Ṽ^s_t depend on Q_t, a jump variable, through the reactivation rate. The two predetermined stocks must therefore be carried separately. (The first version of this note wrote E4 for v_pre with Λ_r(Q_{t+1}) and called v_pre a state, which is inconsistent.) With costless reactivation (Λ_r ≡ 1), v_pre = Ṽ^u + Ṽ^s is predetermined and the state collapses to (u, v_pre).

**Recovered afterwards** (they don't enter E1–E4):
- the wage, w_t = ϕ[ z_t − K_t + θ_t(K_t + q(θ_t)κ) − s_t 𝓡_t ] + (1 − ϕ)b;
- the match surplus;
- the share price V_t and the riskless rate R^f_t (§2).

### How it relates to DMP, term by term

| | DMP | This model |
|---|---|---|
| Vacancy cost | flow cost c | user cost K_t = Q_t − β(1−δ_t)E_t Q_{t+1}: depreciation through δ, plus capital gains if e moves |
| Discounting | constant β (risk-neutral) | m_{t+1} = β(C_{t+1}/C_t)^{−σ}, the owners' discount factor; β when σ = 0 |
| Survival in discounting | 1 − s | 1 − δ discounts the match; s enters the continuation (1 − s)Γ |
| Separation | costless: the firm re-posts | the shortfall s·𝓡 enters the surplus like a firing cost, borne in part (1 − ϕ) by the recruiter |
| Vacancy stock | jump variable | states: carried-over vacancies Ṽ^u and vacated positions Ṽ^s, plus creation e |
| Destruction | not separate from s | removes jobs *and* the vacancies attached (E3 and E4 both carry 1 − δ) |

### Special cases

- **Vacancy free entry (ξ → ∞), with σ = 0:** Q = x_m is constant, so Λ_r and Q^rep are constants, K_t = x_m[1 − β(1 − δ_t)] depends only on δ_t, and 𝓡_t = β(1 − δ_t)(x_m − Q^rep) is proportional to survival. E2 then pins θ_t as a function of the current shocks alone, the analog of `lem:theta` with the 𝓡 term: θ does not load on the stocks (u, Ṽ^u, Ṽ^s). `ags_run_solution.jl` section 2 checks this numerically. In this limit, reactivation is a constant attenuation, not a cyclical margin. s shocks still move θ, through two terms in E2: the expected continuation (1 − s_{t+1})C_{t+1}, as in DMP, and the shortfall (1 − ϕ)s_{t+1}𝓡_{t+1}. (The first version said only the shortfall term; that was wrong.) **With σ > 0 the result fails:** m_{t+1} depends on consumption, which depends on the stocks, so θ loads on them too (`ags_run_solution.jl` section 2 therefore checks it at σ = 0).
- **Linear utility (σ = 0):** m = β, R^f = 1/β; the household block drops out of E2 and the user cost. This is the version Gabrovski and Silva, Shimer and Hagedorn–Manovskii work with.
- **Costless reactivation (p_r = 0, so Λ̄_r = 1):** the formulas give Λ_r ≡ 1 and Q^rep = Q, so 𝓡 = 0, and the scale a has no effect. The limit a → 0 with p_r > 0 is the same economy. (Λ̄_r = 1 with p_r > 0 is not: every position is reactivated, but the costly draws still pay, so 𝓡 > 0.) With κ = 0 this is Gabrovski–Silva with aggregate shocks, and additionally with s ≡ 0 it is Coles–Kelishomi (`rem:nesting`).

### Steady state

**Notation.** f, q and θ are rates conditional on survival (JOLTS gross rates divided by 1 − δ), τ = δ + (1 − δ)s, and β = 1/(1 + r). **In steady state m = β whatever σ,** so the steady state and the calibration below do not depend on σ; σ affects only the dynamics. The household side adds C = z̄(1 − u) − X − X^r, V = D/r and R^f = 1 + r. The suffix _Q denotes a steady-state ratio per unit of Q: K_Q = K/Q, κ_Q = κ/Q, R_Q = 𝓡/Q. These ratios are code shorthand (`calibrate_ags`), not model objects: κ is the per-hire cost (a parameter), K the user cost and 𝓡 the reactivation shortfall (both endogenous).

#### Forward: from parameters to the steady state (`steady_state_ags`)

Given (δ, s, r, A, η, b, ϕ, κ, x_m, ξ_inv, z̄, a, p_r, ψ_r), the steady state reduces to one equation in θ, with a scalar equation for Q nested inside.

**Step 1: flows (the Beveridge curve).** For a given θ:

  q = Aθ^{−η},  f = θq,  u = τ / (τ + (1 − δ)f),  v = θu.

**Step 2: positions.** Creation replaces the positions destroyed with product lines and the vacated positions that are retired:

  e(Q) = δ(v + 1 − u) + (1 − δ)(1 − Λ_r(Q)) s (1 − u),  Q = x_m e(Q)^{ξ_inv}.

Q − x_m e(Q)^{ξ_inv} is increasing in Q (a higher Q raises reactivation, which lowers the creation needed), so Q is unique. With costless reactivation, e is explicit and so is Q; under vacancy free entry (ξ_inv = 0), Q = x_m.

**Step 3: prices of a position.**

  K = Q (r + δ)/(1 + r),  𝓡 = (1 − δ)(Q − Q^rep(Q))/(1 + r).

K is interest plus the risk that destruction removes the position. 𝓡 is the expected value lost at separation, through retirement or reactivation costs.

**Step 4: job creation.** E2 in steady state, with the wage substituted:

  (1 − ϕ)(z̄ − K − b − s𝓡) − ϕθ(K + qκ)  =  (r + τ)/(1 − δ) · (κ + K/q).

The left side is the recruiter's share of the flow surplus. The right side is the cost of a hire, κ + K/q, annuitized over the match's expected life: from E2, (κ + K/q)[1 + r − (1 − δ)(1 − s)] = (1 − δ) × flow, and 1 + r − (1 − δ)(1 − s) = r + τ. Solve for θ; the wage w then follows from the Nash equation.

**Uniqueness.** As θ rises, the left side falls and the right side rises, so the two cross once:
- θ(K + qκ) = θK + fκ rises;
- K rises, because v rises, so e rises, so Q rises;
- 𝓡 rises with Q, because d(Q − Q^rep)/dQ = 1 − Λ_r ≥ 0;
- K/q rises, because q falls and K rises.

This is DMP's job-creation logic, except that the vacancy cost K rises with θ, through creation, instead of being a constant c.

#### Backward: from targets to parameters (`calibrate_ags`)

Calibration runs the same chain in reverse, and every unknown then enters linearly:

| Step | Forward (given parameters) | Backward (given targets) |
|---|---|---|
| Flows | θ unknown → q, f, u, v | f, q from JOLTS → θ, u, v and A directly |
| Positions | solve Q = x_m e(Q)^{ξ_inv} | Λ̄_r is a target, so e is explicit; Q from X/Y; then x_m = Q/e^{ξ_inv} and a from Λ̄_r |
| Job creation | solve for θ | solve for z̄ |
| Wage | gives w | w = 1 is a target, so it gives ϕ |

With the _Q ratios, K_Q = (r + δ)/(1 + r), κ_Q = (1 − x_v)/x_v · K_Q/q and R_Q = (1 − δ)·gap/(1 + r), where gap = (1 − Λ̄_r) + ψ_r/(1 + ψ_r)·(Λ̄_r − 1 + p_r). Job creation then reads

  z̄ − w = zw_Q · Q,  zw_Q = K_Q + s·R_Q + (r + τ)/(1 − δ) · (κ_Q + K_Q/q):

flow profit per worker covers the user cost, the expected reactivation loss and the annuitized hiring cost. Recruiting costs are X = X_Q · Q with X_Q = e/(1 + ξ_inv) + κ_Q·q·v. Since Y = z̄(1 − u) is also linear in Q, the target X/Y pins Q in closed form:

  Q = X_Y · w · (1 − u) / [X_Q − X_Y · zw_Q · (1 − u)].

**Where the cost of a match comes from** (approximate, worked by hand from the calibrated values in `ags_steady_state_checks.txt`; not produced by a program):

| | user cost K_Q | reactivation s·R_Q | annuitized hiring | zw_Q | Q |
|---|---|---|---|---|---|
| Costless, ξ_inv 1 | 92 % | 0 | 8 % | 0.0065 | 9.74 |
| Λ̄_r 0.25, ξ_inv 1 | 21 % | 77 % | 2 % | 0.0283 | 1.27 |

With Λ̄_r = 0.25, the expected reactivation loss is about three quarters of the cost of a match. A fixed X/Y then forces Q down, the flow surplus z̄ − w shrinks, and Nash bargaining hits w = 1 only with ϕ near 1 (see "Implied bargaining power" in §6). The small steady-state share of the hiring cost understates κ's role over the cycle: κ is fixed while K moves with Q, so a larger κ share makes the hiring cost less procyclical (Pissarides 2009). Setting κ = 0 (x_v = 1) lowers σ(u) from z from 5.76 to 3.88 in `ags_run_solution.txt`.

#### Comparison with the baseline (`steady_state.jl`)

| | AGS | Baseline |
|---|---|---|
| Calibration | closed form | Stage 4 root-finds Q: output includes variety and fixed costs, so it is not linear in Q (PATH A adds two more root-finds) |
| Forward solve | one equation in θ, with a scalar equation for Q nested inside | a multi-dimensional system (θ with N and the exit cutoff x_c, through variety and survival), plus the Λ_r ↔ Q iteration in `repost_state` |

**Which one is used when.** In estimation, each draw goes through the calibration: the estimated parameters (ξ_inv, x_v, b/w, Λ̄_r, p_r, ψ_r and the shocks) plus the fixed targets back out the rest, so the long-run targets hold exactly at every draw. In AGS that step is closed form. The forward solve serves two purposes:
- checks: `ags_steady_state_checks.jl` solves forward from the calibrated parameters and recovers the targets;
- counterfactuals that hold parameters fixed (for example, a higher δ̄ with x_m and z̄ unchanged), which calibration cannot do, since it would re-target and absorb the change.

## 4. Observables (for the SVAR)

| Model object | Data counterpart |
|---|---|
| z_t (measured productivity Y/L = z exactly, since there is no variety) | output per worker |
| δ_{e,t} = δ_{t−1} (measured exit is destruction, lagged one period) | BED Deaths rate |
| s_t | s from τ and δ (part6b) |
| u_t, v_t | unemployment and vacancy rates |

On the model side, the SVAR's shock block is therefore observed almost exactly, so the z-first, s, δ ordering recovers the model's shocks up to δ's dating.

## 5. Is reactivation away from its bound?

Λ_r < 1 iff Q_t < a. Two facts settle it.

**In steady state it is interior by construction.** a is calibrated so that Λ_r(Q̄) equals a target Λ̄_r < 1, giving

  Q̄/a = [ (Λ̄_r − 1 + p_r)/p_r ]^{1/ψ_r}.

**Over the cycle, it binds only if Q rises by the factor**

  Q_t/Q̄ > [ p_r / (Λ̄_r − 1 + p_r) ]^{1/ψ_r}.

| Λ̄_r | p_r | ψ_r | Q must rise by | Binds? |
|---|---|---|---|---|
| 0.25 | 0.8 | 1 | 16× | never in practice |
| 0.25 | 0.8 | 0.5 | 256× | never |
| 0.70 | 0.8 | 1 | 1.6× (+60 %) | only in very large booms |
| 0.95 | 0.8 | 1 | 1.07× (+7 %) | **often:** an occasionally binding constraint |

- **At the low reactivation rates the data favor** (v falls after s shocks, so Λ̄_r is low), the bound is far away, and first-order perturbation is safe.
- **It becomes a real concern only if estimation pushes Λ̄_r toward 1.** In that case:
  - report the simulated share of periods with Q_t > a as a diagnostic;
  - if that share isn't negligible, treat the bound as occasionally binding, with a piecewise-linear solution in the style of OccBin (citation unverified).
- **Near vacancy free entry, Q barely moves anyway** (Q_t/Q̄ ≈ (e_t/ē)^{1/ξ}), so the bound is even less of a concern.

**The lower end is not a constraint:** Λ_r ≥ 1 − p_r always (the zero-cost atom), and the power law is smooth for Q > 0.

**Why decoupling from the exit distribution matters here.** Before Nov 4, the baseline code drew reactivation costs from the exit distribution, with shape ψ, about 0.0085 at the working point. Hitting Λ̄_r = 0.25 then requires Q̄/a = 0.0625^{118}. That is why α runs to 10^93–10^178 in `Q2_p0_floor.txt`. It also makes reactivation almost unresponsive to Q: the elasticity of Λ_r − 1 + p_r with respect to Q is ψ_r. With its own ψ_r of order 1, reactivation responds meaningfully, for instance with an elasticity of Λ_r of 0.2 at Λ̄_r = 0.25, p_r = 0.8 and ψ_r = 1. The scale a stays reasonable.

## 6. Parameters

| Group | Parameters | Pinned by |
|---|---|---|
| Preferences, matching | β, A, η, ϕ, b | standard targets (r, f̄, η from the literature, b/w) |
| Positions | x_m, ξ, κ | Q̄ from the steady state; ξ estimated (or ξ → ∞); κ = 0 in GS |
| Reactivation | p_r, ψ_r, a | a from Λ̄_r; p_r and Λ̄_r estimated; ψ_r fixed (for example at 1) or estimated |
| Shocks | δ̄, s̄; (ρ, σ) for z, δ, s | BED and CPS means; part6b's matched processes |

### Calibration (`ags_steady_state.jl`, `calibrate_ags`)

The baseline's `calibrate_shares`, adapted. The long-run targets hold exactly at every parameter draw.
- **Stage 1, conversions:** as in the baseline. f and q are divided by 1 − δ, and e replaces destroyed and retired positions.
- **δ = δ_e:** with no endogenous exit, all of BED Deaths is destruction (MS, Nov 4).
- **Stages 2–3 dropped:** no exit, no variety (ε, N, p_0, ω_δ, the exit elasticity, X^c/Y).
- **Stage 4 in closed form:** z̄ and recruiting costs X are both linear in Q, so the target X/Y pins Q without root-finding.
  - Part 1 targets **X/Y = 1.5 %**, the data target.
  - The Part 1 → Part 2 bridge targets **X/Y = 1.67 %**, the baseline's X/(w_int·L) = 1.5 % × 1.1165 (MS, Nov 4). AGS output is w_int·L, so this is the same recruiting cost relative to match revenue.
- **Stage 5:** ϕ from the Nash wage, x_m = Q/e^{ξ_inv}, a from Λ̄_r.
- **Reactivation costs** are excluded from X, as in the baseline.

**Implied bargaining power.** At Λ̄_r = 0.25, ϕ = 0.92 (0.96 at ξ_inv = 0.01), against 0.61 with costless reactivation. The baseline working point has the same feature (ϕ = 0.914): a low Λ̄_r raises e, so a fixed X/Y forces Q down, leaving little surplus for Nash bargaining to split. Whether to accept ϕ near 1, or target ϕ and let X/Y be an outcome, is open.

## 7. Transmission: when does a shock that raises u also raise v?

> ⚠️ **Nov 6, 2026: the tables in this section are at σ = 0** (linear utility, the default before Nov 6). The default is now σ = 1; §8 gives what changes.
>
> ⚠️ **Nov 6, 2026 (decisions.md S16):** the model numbers in this section are **true-shock** statistics (one shock simulated at a time). The data numbers they are set against (σ(u) by shock 6.32, 3.64, 3.82) come through the SVAR, which does not recover the model's shocks exactly (δ about 55 %). The like-for-like comparison runs the same SVAR on simulated data: `../Bartek analysis/data/results/q3_model_svar.txt` §5. The mechanisms below stand; the model–data gaps quoted here are to be restated from that output.

*Added Nov 5, 2026. Evidence from `../Programs AGS/ags_xi_sweep.txt` (ξ_inv grid at x_v = 0.5 and 0.25, everything else at `AGS_TARGETS`) and `ags_kappa_sweep.txt`. The account of Coles and Kelishomi below is a hypothesis about their setting, not a summary of their paper (citation unverified).*

### 7.1 The accounting

From v = θu,

  d log v = d log θ + d log u.

u is a slow-moving stock; θ is a jump variable. A shock that raises u produces a Beveridge-curve movement (v falls) only if it lowers θ by more, in percent, than it raises u, over the horizon the quarterly HP cycle picks up. Otherwise v rises with u. Two forces lower θ:
- **incentives (E2):** a higher user cost K, a lower match value (more expected separation, through (1 − s_{t+1})), or a larger reactivation shortfall s𝓡;
- **the stock of positions:** with inelastic creation, destroyed or retired positions are replaced slowly, so v falls directly.

### 7.2 Separation shocks (s)

- **Incentives are weak.** Near vacancy free entry, Q and K are almost constant, and s lowers θ only through (1 − s_{t+1}) and s𝓡.
- **The stock effect needs inelastic creation.** At Λ̄_r = 0.25, three quarters of vacated positions are retired. With elastic creation they are replaced at once; with inelastic creation replacement needs Q to rise, which is slow and costly.
- **Unemployment rises strongly** (σ(u) from s is 4.6–6.4).

So with elastic creation d log u dominates and v rises with u: this is Shimer's (2005) separation-shock problem in DMP, where free entry pins θ. Inelastic creation restores the Beveridge relation.

| ξ_inv (x_v = 0.5) | 4 | 1 | 0.5 | 0.25 | 0.1 | 0.01 |
|---|---|---|---|---|---|---|
| cor(u, v), s alone | −0.65 | −0.42 | −0.17 | **+0.19** | +0.64 | +0.98 |

The sign flips at ξ_inv ≈ 0.3 (x_v = 0.5) and ≈ 0.45 (x_v = 0.25). s is the largest shock, so this drags cor(u, v) with all shocks to −0.34 to −0.41 at low ξ_inv; the best value in the grid is −0.66 at ξ_inv = 4 (data −0.925).

### 7.3 Destruction shocks (δ)

δ has a channel that s lacks: it enters the user cost directly. In the vacancy free entry limit, K_t = x_m[1 − β(1 − δ_t)], so a higher δ_t raises the cost of holding every position and lowers θ through E2 even with perfectly elastic creation (the user-cost channel of `prop:fe` in the draft). Destruction also removes vacancies directly, which matters with inelastic creation.

| ξ_inv (x_v = 0.5) | 4 | 1 | 0.5 | 0.25 | 0.1 | 0.01 |
|---|---|---|---|---|---|---|
| cor(u, v), δ alone | −0.85 | −0.81 | −0.79 | −0.79 | −0.89 | −0.49 |
| v trough after a one-SD δ shock, % (month) | −0.32 (7) | −0.25 (7) | −0.21 (6) | −0.17 (6) | −0.12 (5) | −0.76 (0) |

- **The Beveridge relation holds at every ξ_inv.** Inelasticity changes its timing: with inelastic creation v falls gradually to a trough at month 6–7; with elastic creation it drops on impact and recovers fast.
- **It weakens as creation becomes elastic,** because destroyed vacancies are replaced at once, and because the δ shock is transitory (ρ_δ = 0.53 a month): θ recovers within months while u stays high for a year or more.
- **σ(u) from δ stays at 0.52–0.79 throughout.** δ's size problem (§4 of the empirical note's comparison) is separate from ξ.

### 7.4 One destruction shock (Coles and Kelishomi) against two

With s ≡ 0, every separation destroys the position, and the model's single destruction shock does the work of both δ and s. *Hypothesis:* if that shock is calibrated to the total separation rate, with elastic creation and transitory shocks, a shock moves u a lot while θ's fall is brief, which is the mechanism of §7.2. A single shock can then produce a counterfactual (positive) Beveridge relation. Separating δ (destroys jobs and positions, raises the user cost) from s (ends matches, leaves the positions) is what lets δ keep the Beveridge relation and confines the problem to s at elastic creation.

### 7.5 Link to R19

In the full model, the product-line entry Euler pins m(1 − δ), which cancels δ's effect on the user cost (the R19 offset). If so, priced entry removes exactly the channel of §7.3: δ would then behave like s, and could produce positive comovement at elastic creation. That would be a further reason for clone replacement in Part 1, and a sense in which R19 does more than damp volatility. *Not yet checked.*

### 7.6 The trade-off in ξ

- **z amplification wants elastic creation.** At x_v = 0.25, σ(u) from z reaches the data's 6.32 near ξ_inv ≈ 0.07.
- **Three moments want inelastic creation (ξ_inv ≳ 0.5):** a mild cor(u, z) (data −0.22, 1951–2019; Gabrovski and Silva's point), the s Beveridge relation (§7.2), and the δ Beveridge relation's gradual timing (§7.3).
- **The best compromise in the grid** (x_v = 0.25, ξ_inv = 0.5) gives σ(u) from z 3.81 (data 6.32), cor(u, z) −0.27 (−0.22), cor(u, v | s) −0.04, cor(u, v) −0.47 (−0.925).

### 7.7 Checks not yet run

1. AGS with a CK-style destruction shock (δ̄ at the total separation rate, s̄ ≈ 0), at ξ_inv = 0.01 and 1: does cor(u, v | δ) turn positive (§7.4)?
2. The full model at the working point, priced entry against clone replacement: does R19 weaken or flip δ's Beveridge relation (§7.5)?
3. A two-way grid x_v ∈ {0.25, 0.1, 0.05} × ξ_inv ∈ {1, 0.5, 0.25}: does a high per-hire share at moderate inelasticity satisfy §7.6?

## 8. Endogenous discounting: σ = 1 against σ = 0

*Added Nov 6, 2026. Sources: `../Programs AGS/ags_run_solution.txt`, `ags_xi_sweep.txt`, `ags_persistence.csv` (true shocks, model against model); `../Bartek analysis/data/results/q3_model_svar.txt` (through the SVAR, model against data). Block B′'s shock processes in the SVAR comparison.*

### 8.1 What σ = 1 changes in AGS

| | σ = 0 | σ = 1 |
|---|---|---|
| σ(u) from z, through the SVAR: default / compromise (data 6.29) | 1.19 / 2.09 | **1.55 / 2.59** |
| σ(u) from δ, through the SVAR (data 3.63) | 0.79 / 0.77 | 0.81 / 0.80 |
| σ(u) from s, through the SVAR (data 3.80) | 4.44 / 4.13 | 4.50 / 4.21 |
| cor(u, v \| s), through the SVAR (data −0.96) | −0.40 / −0.03 | −0.45 / −0.17 |
| distance J, one-SD responses: default / compromise | 362 / 351 | **325 / 329** |
| cor(u, v \| s), true shocks, ξ_inv = 0.01, x_v = 0.5 | +0.98 | **+0.42** |
| cor(u, v \| δ), true shocks, ξ_inv = 0.01, x_v = 0.5 | −0.49 | **−0.94** |

1. **z gains 25–30 %.** The channel is the discount factor's (§2), but its direction has not yet been traced through the responses of C, m and R^f; that check is pending.
2. **Endogenous discounting disciplines elastic creation.** At σ = 0 and elastic creation, retired and destroyed positions are re-created at once, which produced the anti-Beveridge movement after s shocks (+0.98) and weakened δ's (−0.49). At σ = 1 that re-creation must be paid for out of consumption, R^f rises and slows it, so cor(u, v | s) falls to +0.42 and cor(u, v | δ) tightens to −0.94. The sign of cor(u, v | s) now flips only below ξ_inv ≈ 0.2 (σ = 0: ≈ 0.3).
3. **δ's size and persistence barely change** (σ(u) 0.57 true-shock at ξ_inv = 1; ρ₁(u) ≈ 0.75; u peaks in quarter 1). Neither does any shock's persistence.

### 8.2 Why the full model with clone replacement gains more for δ

In the full model without exit and with clone replacement, σ = 1 doubles δ's effect (σ(u) 1.10 against 0.58 at σ = 0.03); in AGS it does not (0.57 at both). The difference is **forced, costly replacement**. In the full model, destroyed product lines must be re-created at once and entry uses labor, so a δ shock makes current resources scarce whatever firms do, consumption dips, and R^f rises. In AGS clones are free, so a δ shock only removes jobs and positions, and creation can be postponed. The earlier reading that entry labor is negligible holds only at σ ≈ 0: **entry labor matters through its interaction with σ.** Under priced entry, entry is a flexible investment margin and the entry Euler pins m(1 − δ), so σ barely matters (0.54 against 0.56; R19).

### 8.3 Implications

- AGS at σ = 1 is the right Part 1 benchmark for the bridge: it shares the household and the discount factor with Part 2.
- The Part 1 → Part 2 comparison has three steps, not two: AGS (free clones); the full model with clone replacement (costly, forced replacement: adds the consumption–discount channel for δ); the full model with priced entry (flexible replacement: removes it, R19).
- Persistence after δ and s is still missing at σ = 1 (ρ₁(u) 0.73–0.80 through the SVAR, against 0.94–0.96 in the data). Endogenous discounting adds size and comovement, not delay.
