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

Utility is linear (σ = 0), so the discount factor is a constant β. Destroyed product lines are replaced by clones at no cost ("clone replacement"), so there is no business-formation block. The shocks are productivity z, destruction δ and separation s.

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

**Reactivation:** each vacated position draws a cost c from F_r.
- F_r has a mass 1 − p_r at zero, and the rest follows a power law with shape ψ_r on [0, a].
- The position is reactivated iff c ≤ Q_t. So the reactivation probability is

  Λ_r(Q) = min{ (1 − p_r) + p_r (Q/a)^{ψ_r}, 1 },

  and the expected net value of a vacated position is

  Q^rep(Q) = E max{Q − c, 0} = Q·[Λ_r(Q) − ψ_{c,r}(Λ_r(Q) − 1 + p_r)],   ψ_{c,r} = ψ_r/(1 + ψ_r).

**Wages:** Nash bargaining, with worker share ϕ and flow value of unemployment b.

## 3. Equilibrium

Three objects collect the dynamics:
- the **user cost** of a vacancy, K_t ≡ Q_t − β(1 − δ_t) E_t Q_{t+1};
- the **hiring cost** per match, C_t ≡ κ + K_t / q(θ_t);
- the **reactivation shortfall**, 𝓡_t ≡ β(1 − δ_t) E_t[ Q_{t+1} − Q^rep(Q_{t+1}) ], the expected value lost because reactivation isn't free.

**Definition (equilibrium).** Given u_0, Ṽ^u_0, Ṽ^s_0 and the shock processes for (z_t, δ_t, s_t), an equilibrium is a stochastic sequence {θ_t, e_t, u_{t+1}, Ṽ^u_{t+1}, Ṽ^s_{t+1}}, with Q_t = Q(e_t), such that for all t:

**(E1) Vacancies:**
  θ_t u_t = v_{pre,t} + e_t,   v_{pre,t} = Ṽ^u_t + Λ_r(Q_t) Ṽ^s_t.

**(E2) Job creation:**
  C_t = β(1 − δ_t) E_t[ (1 − ϕ)(z_{t+1} − K_{t+1} − b) − ϕ θ_{t+1} K_{t+1} − ϕ f(θ_{t+1}) κ + (1 − s_{t+1}) C_{t+1} − (1 − ϕ) s_{t+1} 𝓡_{t+1} ].

**(E3) Unemployment:**
  u_{t+1} = u_t − (1 − δ_t) f(θ_t) u_t + τ_t (1 − u_t),   τ_t = δ_t + (1 − δ_t) s_t.

**(E4) Carried-over positions:**
  Ṽ^u_{t+1} = (1 − δ_t)(1 − q(θ_t)) v_t   (unmatched vacancies),
  Ṽ^s_{t+1} = (1 − δ_t) s_t (1 − u_t)      (positions vacated by separation, reactivated at t+1 or retired).

K_t, C_t and 𝓡_t are the functions of (θ, e, δ) defined above. That's three state variables (u, Ṽ^u, Ṽ^s) and two jump variables (θ, e), against one state (u) and one jump (θ) in DMP.

*Why v_pre is not a state.* Inherited vacancies v_{pre,t} = Ṽ^u_t + Λ_r(Q_t)Ṽ^s_t depend on Q_t, a jump variable, through the reactivation rate. The two predetermined stocks must therefore be carried separately. (The first version of this note wrote E4 for v_pre with Λ_r(Q_{t+1}) and called v_pre a state, which is inconsistent.) With costless reactivation (Λ_r ≡ 1), v_pre = Ṽ^u + Ṽ^s is predetermined and the state collapses to (u, v_pre).

**Recovered afterwards** (they don't enter E1–E4):
- the wage, w_t = ϕ[ z_t − K_t + θ_t(K_t + q(θ_t)κ) − s_t 𝓡_t ] + (1 − ϕ)b;
- the match surplus;
- consumption.

### How it relates to DMP, term by term

| | DMP | This model |
|---|---|---|
| Vacancy cost | flow cost c | user cost K_t = Q_t − β(1−δ_t)E_t Q_{t+1}: depreciation through δ, plus capital gains if e moves |
| Survival in discounting | 1 − s | 1 − δ discounts the match; s enters the continuation (1 − s)C |
| Separation | costless: the firm re-posts | the shortfall s·𝓡 enters the surplus like a firing cost, borne in part (1 − ϕ) by the recruiter |
| Vacancy stock | jump variable | states: carried-over vacancies Ṽ^u and vacated positions Ṽ^s, plus creation e |
| Destruction | not separate from s | removes jobs *and* the vacancies attached (E3 and E4 both carry 1 − δ) |

### Special cases

- **Vacancy free entry (ξ → ∞):** Q = x_m is constant, so Λ_r and Q^rep are constants, K_t = x_m[1 − β(1 − δ_t)] depends only on δ_t, and 𝓡_t = β(1 − δ_t)(x_m − Q^rep) is proportional to survival. E2 then pins θ_t as a function of the current shocks alone, the analog of `lem:theta` with the 𝓡 term: θ does not load on the stocks (u, Ṽ^u, Ṽ^s). `ags_run_solution.jl` section 2 checks this numerically. In this limit, reactivation is a constant attenuation, not a cyclical margin. s shocks still move θ, through two terms in E2: the expected continuation (1 − s_{t+1})C_{t+1}, as in DMP, and the shortfall (1 − ϕ)s_{t+1}𝓡_{t+1}. (The first version said only the shortfall term; that was wrong.)
- **Costless reactivation (p_r = 0, so Λ̄_r = 1):** the formulas give Λ_r ≡ 1 and Q^rep = Q, so 𝓡 = 0, and the scale a has no effect. The limit a → 0 with p_r > 0 is the same economy. (Λ̄_r = 1 with p_r > 0 is not: every position is reactivated, but the costly draws still pay, so 𝓡 > 0.) With κ = 0 this is Gabrovski–Silva with aggregate shocks, and additionally with s ≡ 0 it is Coles–Kelishomi (`rem:nesting`).

### Steady state

**Notation.** f, q and θ are rates conditional on survival (JOLTS gross rates divided by 1 − δ), τ = δ + (1 − δ)s, and β = 1/(1 + r). The suffix _Q denotes a steady-state ratio per unit of Q: K_Q = K/Q, κ_Q = κ/Q, R_Q = 𝓡/Q. These ratios are code shorthand (`calibrate_ags`), not model objects: κ is the per-hire cost (a parameter), K the user cost and 𝓡 the reactivation shortfall (both endogenous).

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
