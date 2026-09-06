# Laws of Motion: Timing Conventions and Position Conservation
**Author:** Mario Silva (with Claude) · **Date:** September 5, 2026
**Purpose:** Resolve an inconsistency between the unemployment and vacancy laws of motion
discovered when the steady-state guard was added to `solution_interface`. Establishes which
equation is the odd one out, and what it implies for Gabrovski-Silva (JEDC).

---

## 1. Timing within a period

From `run_solution_core.jl` (monthly, following the CK stage structure):

| Stage | Event |
|---|---|
| 1 | Shocks realize: (z_t, δ_t, s_t) |
| 2 | Bargaining and **production** with the inherited workforce L_t = 1 − u_t |
| 3 | Vacancy investment: entrants e_t post |
| 4 | **Matching**: m_t = q(θ_t)v_t = f(θ_t)u_t matches form, where v_t = v^pre_t + e_t |
| 5 | **Separation**: firms die at δ_{e,t}; surviving matches separate at s_t |

Two facts follow and drive everything below.

- **u_t is pre-matching unemployment** — the stock entering t, before stage 4. Hence
  L_t = 1 − u_t is the workforce that *produces* at stage 2.
- **v_t is the stock at the time of matching** — after entry, before matching. Matches are
  formed out of v_t, and (1 − q(θ_t))v_t remain unfilled.

Define the **post-matching pool** of filled jobs, i.e. the set of matches in existence when
stage 5 arrives:

$$E_t \;\equiv\; \underbrace{(1-u_t)}_{\text{inherited}} \;+\; \underbrace{f(\theta_t)u_t}_{\text{just matched at stage 4}}$$

Everything hinges on one question: **at stage 5, is the newly matched flow f(θ_t)u_t exposed
to the separation rate s_t, or only to firm death δ_{e,t}?**

---

## 2. The position-conservation identity

Count *positions* — filled plus vacant — rather than workers. At stage 5:

- A firm dies (prob. δ_{e,t}): **all** its positions vanish, filled and vacant alike. This is
  the paper's core mechanism and is why (1−δ_{e,t}) multiplies the whole bracket in the
  vacancy LOM.
- A match separates at a **surviving** firm (rate s_t): the position is *not* destroyed. It
  converts from filled to vacant and is reposted at no additional sunk cost.

So positions are conserved up to firm death, and any correct pair of laws of motion must
satisfy

$$\boxed{\;v^{pre}_{t+1} \;+\; L_{t+1} \;=\; (1-\delta_{e,t})\bigl(v_t + L_t\bigr)\;}\tag{PC}$$

This is the diagnostic. A leak in (PC) means positions are disappearing without a firm dying —
a worker leaves employment and no vacancy comes back.

---

## 3. Convention A — pre-matching (Gabrovski-Silva, old Dynare, `compute_unemployment`)

New matches face **firm death only** in their first period; they cannot separate at rate s
until they have produced.

Unemployment (matches formed at t survive to t+1 only if the firm does, hence the
$(1-\delta_{e,t})$ on the outflow; separations hit the inherited stock $1-u_t$ only):

$$u_{t+1} = \bigl[1-(1-\delta_{e,t})f(\theta_t)\bigr]u_t \;+\; \tau_t\,(1-u_t)$$

Vacancies:

$$v^{pre}_{t+1} = (1-\delta_{e,t})\bigl[(1-q(\theta_t))v_t + s_t(1-u_t)\bigr]$$

The separation base and the reposting base are **the same object**, $1-u_t$. Employment:

$$L_{t+1} = (1-\tau_t)(1-u_t) + (1-\delta_{e,t})f(\theta_t)u_t$$

**Check (PC).** Using τ = 1 − (1−δ_e)(1−s) so that 1 − τ = (1−δ_e)(1−s), and q v = f u:

$$
\begin{aligned}
v^{pre}_{t+1} + L_{t+1}
&= (1-\delta_e)(1-q)v \;+\; (1-\delta_e)s(1-u) \;+\; (1-\delta_e)(1-s)(1-u) \;+\; (1-\delta_e)fu \\[2pt]
&= (1-\delta_e)\bigl[(1-q)v + (1-u) + qv\bigr] \\[2pt]
&= (1-\delta_e)\,(v + L) \quad\checkmark
\end{aligned}
$$

**Conserved exactly.** Steady state:

$$(1-\delta_e)fu = \tau(1-u) \quad\Longrightarrow\quad u = \frac{\tau}{\tau + (1-\delta_e)f}$$

which is exactly the original `compute_unemployment`, draft `eq:u_ss`, and the closed form
`e = δ_e(v + 1 − u)` (`eq:e_ss`). **The whole block is mutually consistent.**

---

## 4. Convention B — post-matching (current draft `eq:u_lom` and model `f[23]`)

New matches are exposed to the full separation rate τ in their first period. The draft states
this explicitly: τ_t applies to last period's employed "*including those who matched in t−1*."

$$u_{t+1} = (1-f(\theta_t))u_t \;+\; \tau_t\bigl[(1-u_t) + f(\theta_t)u_t\bigr] = (1-f)u_t + \tau_t E_t$$

Employment is then simply $L_{t+1} = (1-\tau_t)E_t$. For (PC) to hold, the reposting base must
be the *same* pool the separations were drawn from:

$$v^{pre}_{t+1} = (1-\delta_{e,t})\bigl[(1-q(\theta_t))v_t + s_t\,E_t\bigr]$$

**Check (PC).**

$$
\begin{aligned}
v^{pre}_{t+1} + L_{t+1}
&= (1-\delta_e)(1-q)v \;+\; (1-\delta_e)sE \;+\; (1-\delta_e)(1-s)E \\[2pt]
&= (1-\delta_e)\bigl[(1-q)v + E\bigr]
\end{aligned}
$$

and since $E = L + qv$, this is $(1-\delta_e)(v+L)$ ✓. Steady state:

$$fu = \tau E \quad\Longrightarrow\quad u = \frac{\tau}{\tau + (1-\tau)f}$$

Note $qv = fu = \tau E$, from which `eq:e_ss` follows exactly:

$$e = \delta_e v + (1-\delta_e)qv - (1-\delta_e)sE = \delta_e v + (1-\delta_e)(\tau-s)E = \delta_e\bigl[v + (1-\tau)E\bigr] = \delta_e(v + 1 - u)$$

using τ − s = δ_e(1−s). So **`eq:e_ss` holds under both conventions** — it is not diagnostic.

---

## 5. The actual defect: the two conventions are mixed

The current model uses **Convention B for unemployment** (`f[23]`, `eq:u_lom`) and
**Convention A for vacancies** (`f[22]`, `eq:v_lom`), plus Convention A's steady-state
formula (`compute_unemployment`, `eq:u_ss`). That mixture violates (PC):

$$
v^{pre}_{t+1} + L_{t+1} = (1-\delta_e)(v + L) \;-\; \underbrace{(1-\delta_e)\,s\,q\,v}_{\text{leak}}
$$

**Numerically** (baseline calibration, `dest_ann` = 0.0754):

| | value |
|---|---|
| $(1-\delta_e)(v+L)$ | 0.95864439 |
| Convention B v LOM: $v^{pre}+L'$ | 0.95864439 ✓ |
| Mixed (current): $v^{pre}+L'$ | 0.95791733 ✗ |
| leak $=(1-\delta_e)s\,q\,v$ | 7.2707e-4 |

**Interpretation of the leak.** It is exactly the workers matched in period t who separate at
the end of period t. The u LOM counts their separation — they reappear as unemployed in t+1 —
but the v LOM does not repost their job. The position evaporates.

The mixture is also what produced the 8.18e-4 residual on `f[23]` once the free-entry bug in
`SS_numeric` was fixed and the guard was added.

---

## 6. Question 1 — Was the vacancy LOM inconsistent all along?

**No.** The vacancy LOM has been unchanged and correct within its own convention throughout:

| Source | u LOM | v LOM | (PC) |
|---|---|---|---|
| `GS.jl` f[8]–f[9] (Gabrovski-Silva JEDC) | Convention A | Convention A | ✅ exact |
| `model.mod:184,187` (old Dynare) | Convention A | Convention A | ✅ exact |
| Current `f[22]`, `f[23]` | **Convention B** | Convention A | ❌ leaks |

Verified numerically: GS.jl gives leak = 0.0 identically, and its steady-state
u = 0.070294 — precisely the value the original `compute_unemployment` returns.

So the inconsistency was **introduced when the unemployment LOM was rewritten to Convention B**
without propagating the change to the vacancy LOM or to the steady-state formulas. The vacancy
LOM is the *unchanged* equation; it only looks wrong when read against the new u LOM.

⚠️ This corrects an earlier conclusion in this investigation that "the typo is fundamentally in
the vacancy law of motion." That was wrong. The v LOM is the survivor of a consistent older
system; the u LOM is what moved.

---

## 7. Question 2 — Implications for Gabrovski-Silva (JEDC)

**None. GS requires no correction.** It uses Convention A consistently across both laws of
motion, and satisfies (PC) exactly. Its steady-state unemployment formula
$u = \tau/(\tau + (1-\delta)f)$ is the correct solution of *its* law of motion, not an
approximation to Convention B's.

Two consequences worth recording for the current paper:

1. **The correction is internal to this paper**, not an erratum on GS. Any text implying the
   GS formulation was flawed would be wrong.
2. **The two papers use different timing conventions**, and the current draft nests GS-adjacent
   benchmarks (`prop:bgm_nest`, the CK comparison in §5.1). If Convention B is retained, the
   difference should be stated where the benchmarks are introduced — under Convention B a match
   can dissolve before ever producing, under Convention A it cannot. It is a substantive
   modeling assumption about whether very-short-tenure separations exist, not a normalization.

---

## 8. Two coherent ways to fix

**Option A — revert to Convention A (minimal).** Restore the GS-form u LOM in `f[23]` and
`eq:u_lom`. Then `compute_unemployment`, `eq:u_ss`, `eq:v_lom`, `eq:e_ss` and all four
`e = δ_e(v+1−u)` sites are correct exactly as originally written, and the September 5 change to
`compute_unemployment` should be reverted. Cost: contradicts the draft's explicit sentence that
τ applies to employed workers "including those who matched in t−1," which would have to go.

**Option B — complete the move to Convention B.** Keep `f[23]`/`eq:u_lom` and fix the three
places still on Convention A:
- `eq:v_lom` (`Draft.tex:633`): `s_{t-1}(1-u_{t-1})` → `s_{t-1}[(1-u_{t-1}) + f(\theta_{t-1})u_{t-1}]`
- `f[22]` (`run_solution_core.jl`): `s*sbar*(1 - u)` → `s*sbar*((1 - u) + f_t*u)`
- `steady_state_checks.jl:95`: same substitution in the assert

plus `eq:u_ss` (`Draft.tex:3273`) and `compute_unemployment` — **already done** — and the
derivation at `Draft.tex:3289`, whose substitution $(1-\delta_e)qv = \tau(1-u)$ becomes
$\tau(1-u)/(1-s)$ under Convention B (the endpoint `eq:e_ss` is unaffected).

**Steady-state impact of the choice:** u = 7.029% (A) vs 7.194% (B), a 2.3% relative
difference, with corresponding small shifts in v, θ and every pp normalization in the IRF
figures.

**Effect on the propositions:** under Convention B the reposting flow is larger by f·u (≈3.2%
of the base), so Proposition 5 Part 1 Channel 1 (∂v_pre/∂s_t > 0) and `lem:vpre`'s bracket both
gain magnitude. Neither changes sign, and no stated condition changes. Convention B mildly
strengthens both results; Convention A leaves them exactly as published.

---

## 9. Status

- ✅ `compute_unemployment`, `L_fun`, three inline copies, `SS_symbolics` — moved to Convention B
- ✅ `e_fun` — restored to `δ_e(v + 1 − u)` (valid under **both** conventions; an intermediate
  "exact LOM" rewrite was wrong and is reverted)
- ❌ `f[22]`, `eq:v_lom`, `steady_state_checks.jl:95` — still Convention A; **decision pending**
- Current state: `steady_state_checks.jl` passes every assert except line 95, and N = 1.0 exactly.
