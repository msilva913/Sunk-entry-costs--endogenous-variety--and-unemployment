# Q1 plan: R8 + D1 + PATH B in one regeneration pass
**Drafted:** October 14, 2026 · **Branch:** `costly_vacancy_reposting` · **Status:** approved Oct 14 (MS), not started

> **Choices confirmed by MS (Oct 14):**
> 1. Reposting is calibrated and estimated through the target **Λ̄_r**, with α backed out.
>    Λ̄_r ∈ [1−p_0, 1) must be strictly interior for perturbation, so the costless benchmark is
>    a switch that sets Λ_r ≡ 1, not Λ̄_r = 1.
> 2. **Lagged δ is added as a state**, so δ_e stays a model variable.
> 3. **§5.3 regeneration is deferred to Q5.**
> 4. Comparison B is re-specified in code in Q1 and run in Q5.

**Goal.** Bring the code in line with the draft (`Draft.tex`, the reference model, S10) so
every later quantitative step (Q2–Q6, `context/pending_tasks.md`) runs on the model the paper
states. Three changes ship together because each invalidates every §5 number:
- **R8** — draft timing plus the partial-reposting margin;
- **D1** — `dest_ann = 0.0320`;
- **PATH B** — the default calibration path, which D1 needs in order to converge.

**Spec.** `context/model_equations.md` ("Code status", Block 1 f[1]–f[33], Block 2).
**Files.** `Programs baseline/run_solution_core.jl` (model), `Programs baseline/steady_state.jl`
(calibration and SS), `Programs baseline/steady_state_checks.jl` (certification), and the
runners that build `TARGETS`.

---

## What changes, and why

### A. Timing (R8 (i)) — `run_solution_core.jl`
Today the code depletes stocks with Λ_t but prices survival with Λ_{t+1}. Its states are
post-exit stocks, which cannot be predetermined, because N_{t+1}, u_{t+1} and v_{pre,t+1}
depend on χ^c_{t+1}.

- **New states** (pre-exit, after the end-of-(t−1) incidence):
  - Ñ_t = (1−δ_{t−1})(N_{t−1}+N^e_{t−1})
  - Ẽ_t = (1−δ_{t−1})[(1−s_{t−1})(1−u_{t−1}) + f_{t−1}u_{t−1}]
  - Ṽ^u_t = (1−δ_{t−1})(1−q_{t−1})v_{t−1}
  - Ṽ^s_t = (1−δ_{t−1})s_{t−1}(1−u_{t−1})
  - plus z_t, δ_t, s_t, and **δ_{t−1}** (see D-choice 2).
- **N_t, u_t, v_{pre,t} become same-period controls**, solved jointly with χ^c_t:
  - N_t = Λ_tÑ_t
  - u_t = 1 − Λ_tẼ_t
  - v_{pre,t} = Λ_t(Ṽ^u_t + Λ_{r,t}Ṽ^s_t)
- **Laws of motion** (four new f's) replace f[22]–f[24], with δ_t the current state:
  - Ñ' = (1−δ)(N+N^e)
  - Ẽ' = (1−δ)[(1−s)(1−u)+f u]
  - Ṽ^u' = (1−δ)(1−q)v
  - Ṽ^s' = (1−δ)s(1−u)
- **Value equations f[3]–f[5] are unchanged.** They already use (1−δ_t)Λ_{t+1}.
- **f[2]:** δ_{e,t} = 1 − (1−δ_{t−1})Λ_t, using the lagged-shock state.

### B. Reposting margin (R8 (ii)) — `run_solution_core.jl`
- **New controls:**
  - Λ_{r,t} = min{(1−p_0) + p_0(Q_t/(αχ_m))^ψ, 1}; the interior branch suffices for perturbation;
  - 𝓡_t = E_t[m_{t+1}(1−δ_t)Λ_{t+1}(Q_{t+1} − Q^rep_{t+1})], with
    Q − Q^rep = Q[(1−Λ_r) + ψ_c(Λ_r − 1 + p_0)] on the interior branch.
- **f[3] (JCC)** gains −(1−ϕ)s'𝓡' inside the expectation (`eq:jcc_eq`).
- **f[7] (wage)** gains −ϕ s 𝓡 (`eq:wage_eq`). Note the coefficients: (1−ϕ) in the
  wage-substituted surplus, ϕ in the wage.
- **f[16] (resource constraint)** gains
  X^r_t = Λ_tṼ^s_t·ψ_c·Q_t(Λ_{r,t} − 1 + p_0) (`eq:agg_repost_costs`). Do **not** multiply by
  Λ_r again: the partial expectation already embeds it.
- **New parameter α**, appended to the `parameters` vector.

### C. Calibration and steady state — `steady_state.jl`
- **Steady state.** The SS itself is unaffected by the timing change (Λ_t = Λ_{t+1}), but it
  changes with reposting:
  - entry e includes replacement of retired positions (`eq:e_ss`/`eq:e_theta`): +(1−δ_e)(1−Λ_r)s(1−u);
  - the JCC and the wage carry −s𝓡;
  - the resource constraint carries X^r.
  - Transcribe from the draft's steady-state appendix (`app:steady_state`, reposting added
    Oct 1).
- **D-choice 1 (recommended): calibrate to a reposting-rate target Λ̄_r, recover α.** With Λ̄_r
  given, e and 𝓡 are known functions of Q, so the Stage 4 root-find on Q keeps its form. Then
  α = Q̄ / (χ_m·((Λ̄_r − 1 + p_0)/p_0)^{1/ψ}) in Stage 5, once χ_m (`f_m`) is known. Λ̄_r ∈
  [1−p_0, 1] is interpretable and has a data analog (share of separations reposted). It also
  becomes the estimated parameter in Q6 in place of α.
- **PATH B as the default (settles D2).** `TARGETS` gains `dest_elast_target = 5.0`, and takes
  the §5.3 values `b_ratio = 0.9`, `x_v = 0.5`, plus `dest_ann = 0.0320`.
  - PATH A stays available behind the existing switch.
  - **Placeholder** `Λr_target`: 1.0 (costless) for the regression tests, then a reported grid
    (principle 21: label it a placeholder; the value is estimated in Q6).
- **`steady_state(para)`** (the round-trip solver used by the checks) gains Λ_r. It is
  implicit through Q = e^{ξ_inv}x_m and e(Λ_r), so add it as a third unknown, or solve an
  inner 1-D fixed point.
- **`SS_numeric`** maps the new states and controls:
  - Ñ = N/Λ, Ẽ = (1−u)/Λ, Ṽ^u = (1−δ)(1−q)v, Ṽ^s = (1−δ)s(1−u);
  - Λ_r, 𝓡, and δ_lag = 1.
- **Retire or fence `SS_symbolics`.** It is latent and already flagged as defective; it must
  not silently diverge further.

### D. Design choices for you to confirm
1. **Reparametrize reposting by the target Λ̄_r** (recommended above) rather than α.
2. **δ_e as a model variable** by adding the lagged shock δ_{t−1} as a state (one extra state).
   Block M targets the δ series, and `run_solution_delta_target.jl` reads `:δ_e`.
   - The alternative is to post-process δ_e from simulated paths.
   - τ_t stays post-processed as now, or add s_{t−1} likewise.
3. **Defer the §5.3 prose regeneration to Q5.** Q1 delivers verified code and baseline IRFs.
   Comparisons A–D are re-run in Q5 at the Q2 shock processes and over the sensitivity grid,
   so regenerating §5.3 now would mean doing it twice. §5.3 numbers are already stale and
   flagged.
4. **Re-specify Comparison B in code now (R8 (v)), run it in Q5.** With α present, p_0 = 0
   removes exit *and* reactivation. The exit-only switch is p_0 > 0 with ω_δ = 0.

---

## Order of work and tests
Each step ends with a passing test before the next starts.

1. **Baseline snapshot.** Run the current code at p_0 = 0 and at the current default. Save the
   IRFs (`irf_baseline_pre_R8.jls`) for regression.
2. **Timing restructure only** (A), with reposting off (α absent, Λ_r ≡ 1).
   - **Test 1 (exact):** at **p_0 = 0**, Λ ≡ 1 and both timings coincide, so the new IRFs must
     equal the snapshot to ~1e-10.
   - **Test 2:** at p_0 > 0, SS residual < 1e-12, BK root count correct, and the position
     identity v_{pre,t+1} + L_{t+1} = (1−δ_{e,t+1})(v_t + L_t).
3. **Steady state with reposting** (C).
   - **Test 3:** `steady_state_checks.jl` extended to the new JCC, wage, e and X^r, at 1e-12.
     Round-trip calibrate → solve recovers the targets, including Λ̄_r.
4. **Dynamics with reposting** (B).
   - **Test 4:** Λ̄_r = 1 (α → 0) and p_0 = 0 reproduce step 2 exactly.
   - **Test 5:** positions are conserved up to exit and the retired mass
     (1−δ_{e,t+1})(1−Λ_{r,t+1})s_t(1−u_t).
   - **Test 6:** in the κ = σ = 0, ξ → ∞, p_0 = 0 limit the code reproduces the `prop:fe`
     coefficients c_δ, c_s from `gs_free_entry_beveridge.jl`. Use a small ξ_inv, as
     `gs_finite_xi_entry_offset.jl` does.
5. **D1 + PATH B defaults** (C). Confirm convergence at `dest_ann = 0.0320`. Report the
   calibrated parameter vector and δ_e/τ = 0.087.
6. **Update the dependents.**
   - Runners that build `TARGETS` (`run_solution_*.jl`, `run_prop5_*.jl`, `run_xi_*`): check
     each still runs. Variables are accessed by name, so most need only the target changes.
   - Comparison B per D-choice 4.
7. **Baseline output.** One run at the new defaults: IRFs to z, δ, s, plus cor(u,v) over a
   Λ̄_r grid as a first look (not draft-bound until Q5).

## Afterwards (context, same commit as the code)
- `model_equations.md`: "Code status" → resolved; record the new state vector.
- `pipeline.md`: the new `TARGETS`; `SS_symbolics` status.
- `decisions.md`: D2 settled (PATH B); D11 implemented; Λ̄_r parametrization.
- `parameters.md`: α → Λ̄_r in Θ_e; the new defaults.
- `pending_tasks.md`: Q1 / R8 done; Q5 inherits the §5.3 regeneration and Comparison B runs.

## Estimated effort
Two to three sessions. Steps 1–2 take one session (the restructure is mechanical but touches
the state vector everywhere). Steps 3–4 take one. Steps 5–7 take a half to one.

## Risks
- **BK at the new timing.** The fixed point N_t = F(χ^c_t)Ñ_t with χ^c_t depending on N_t
  adds same-period simultaneity. `eq:exit_cond` says it is well behaved when the exit margin is
  inelastic enough. Check the root count at the PATH B calibration, where the left side was
  0.002–0.005.
- **The Λ_r saturation kink** (Λ_r = 1 for Q ≥ αχ_m). Perturbation needs the steady state
  strictly interior: Λ̄_r < 1. The Λ̄_r = 1 regression must be run with α → 0 directly, not at
  the kink.
- **Symbolic Jacobian size** grows by about 6 variables. This is acceptable, but recompilation
  time rises.
