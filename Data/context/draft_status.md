# Draft Status — `Draft/Draft.tex`
**Last updated:** October 5, 2026 (Prop. 5 restated; earlier Oct 1 work: wage dependencies, `def:equilibrium`, timing prose)

**Build:** compiles with 0 errors and 0 missing graphics. The 4 undefined references and 1
duplicate label (P2 in [`pending_tasks.md`](pending_tasks.md)) predate the reposting work.
Compile from inside `Data/Draft/` (the graphicspath leads with relative entries), two passes,
with `-synctex=1`. If an editor auto-compiles on save, the `.aux` can be truncated
mid-build; delete it and rebuild.

**Notation.** χ^c_t (not x^c_t); Λ_t ≡ F(χ^c_t); ς for the survival quantile. Reposting:
Λ_{r,t} = F(Q_t/α), option value Q^rep_t (`eq:Qrep`), shortfall written inline as Q − Q^rep;
the symbols Ψ, D, M are retired (clashes with ψ, D^int, 𝒟, m). 𝓡_t is the expected
discounted shortfall in the wage appendix. Label prefixes `eq:`/`tab:`/`fig:`/`sec:`/`app:`.

**Timing.** The draft is the reference model ([S10](decisions.md)); the code deviates in two
places until R8 ([`model_equations.md`](model_equations.md) "Code status"). Until R8 and the D1
cascade, every §5.3 number comes from the code's timing.

---

## Section status

| Section | Status | Outstanding |
|---|---|---|
| 1 Introduction (incl. §1.1 evidence) | ✅ | Prose says reposting makes vacancies "self-correct" (~L154, ~L308); qualify for partial reposting (R11) |
| 2 Environment | ✅ | Convention A `eq:u_lom`; `eq:v_lom` with Λ_{r,t}. Timing prose (Oct 1, late) states the paper's convention: the aggregate state incl. s_t, δ_t is known at the start of t; only incidence is revealed at the end of t. `fig:Timing` unchanged (MS: consistent with the convention); it notes that s hits only incumbent matches |
| 3 Equilibrium | ✅ | Wage dependencies fixed Oct 1: 𝓡_t defined at `eq:repost_shortfall`; `eq:wage_eq` carries −s_t𝓡_t; `eq:surplus_wage` carries −(1−ϕ)s_t𝓡_t. `eq:Lambda_r` defines the reactivation rate. `def:equilibrium` audited for completeness Oct 1: 4 equations for (N, C, χ^c, Y^c), 6 for (θ, K, Q, e, v, u) with `eq:theta_eq` added. The job creation condition `eq:jcc_eq` substitutes the surplus at t+1, so J_t and w_t drop out (`eq:surplus_eq` removed); 𝓡_t is carried with its own equation, so there is no nested expectation. Auxiliary list cites `eq:Lambda_r`, the matching rates, and 𝓡_t. Prop. 1 proof updated to match |
| 3.x Limiting cases | ✅ | — |
| 3.x Shock transmission (`sec:mechanism`) | ✅ | Prop. 5 and its discussion rewritten Oct 5; the δ̄_e/τ̄ ≈ 0.21 claim here was removed |
| 3.x Steady state | ✅ | — |
| 4.1–4.2 Instruments, LP specification | ✅ | — |
| 4.3 Results | ⚠️ | δ–LD paragraph; joint-LP sentence; malformed `\ref` at the end (P1) |
| 5.1 Relation to Coles and Kelishomi | ✅ | δ̄ figures in the D1 cascade |
| 5.2 Calibration and estimation | ⚠️ contested | Documents PATH A, figures use PATH B ([D2](decisions.md)); β̂ ↔ β(θ) unstated ([D3](decisions.md)); α absent; `app:weighting_robustness` unwritten |
| 5.3 Inspecting the mechanism | ✅ | ⏳ Regeneration owed (R8 + D1 + PATH B). Numbers are checked against `mechanism_stats.txt`, which is the authority. The "Choice of shocks" paragraph now says the comparisons assume costless reposting |
| 5.4 Posterior estimates (`sec:posterior`, L3310) | ❌ empty stub | Blocked on the estimation build |
| 6 Conclusion (L3314) | ❌ not started | Label is `ref:conclusion`; references use `sec:conclusion` |

## Appendix status

| Appendix | Status |
|---|---|
| A Relation to literature | ✅ |
| B Data sources | ✅ |
| C Instrument diagnostics | ✅ |
| D Steady state | ✅ — reposting added Oct 1 (shortfall in the job creation condition and wage, Λ_r, Q^rep, X^r, entry `eq:e_ss` with the non-reactivated term); P5: X^c still uses the p_0ψ_cχ^c shortcut |
| E Derivations — profit share, resource constraint, labor share, **wage determination** (rederived Oct 1: Convention A employment law with δ_{e,t+1}, the reposting shortfall 𝓡_t, the wage-substituted surplus with (1−ϕ)), log-linear propagation (`eq:u_ll` rederived for Convention A) | ✅ |
| F Additional comparisons (incl. Comparison D) | ✅ |
| G Limiting-case proofs (`prop:independence`, `prop:double_limit`, `prop:ags`, `prop:bgm_nest`) | ✅ |
| H Proof of `prop:ds_asymmetry` (`app:proof_ds`) | ✅ restated Oct 1 |
| I Steady-state equilibria (`prop:equilibria`) | ✅ |
| `app:weighting_robustness`, `app:robustness` | ❌ referenced, do not exist |

---

## Propositions

| Proposition | Status |
|---|---|
| `prop:independence` (1) | ✅ Holds with its three conditions (σ = 0, ζ = 0, p_0 = 0). p_0 = 0 collapses F to a point mass, so Λ_r = 1, Q^rep = Q, X^r = 0: reposting vanishes with exit, so no fourth condition is needed. J_t is in the stated labor block |
| `prop:double_limit` (2) | ✅ Unaffected (assumes p_0 = 0) |
| `prop:ags` (3) | ✅ Holds; states the clone-replacement motivation and that p_0 = 0 is assumed (no-exit also follows from f_e ≥ χ_m, which would leave reactivation on). AGS differs from the baseline on three channels |
| `prop:bgm_nest` | ⚠️ Not audited for reposting/timing (R12). At p_0 = 0 the conservation leak is zero |
| `prop:equilibria`, `prop:curves` | ⚠️ Not audited (R12) |
| `prop:ds_asymmetry` (5) | ✅ Restated Oct 5 under the timing convention (R13); Part 2 numerics wait on R8 |

### Proposition 5 (`prop:ds_asymmetry`, proof in `app:proof_ds`) — restated Oct 5, 2026

**Framing.** Part 1 is a benchmark: with p_0 = 0 there is neither endogenous exit nor costly
reactivation, and s shocks produce positive u–v comovement. Negative comovement must come from
δ shocks with a large enough role (Part 3) or from partial reposting. All three parts use the
paper's timing convention: the aggregate state is known at the start of t, so θ_t and e_t
respond on impact.

- **Part 1** (p_0 = 0): u_t unchanged; v_t moves only through e_t. At h = 1,
  Δu = d − (1−δ)(1−η_L)q·Δe_t and Δv_pre = d + (1−δ)[1−(1−η_L)q]·Δe_t, with
  d = (1−δ_t)(1−u_t). Positive comovement under (i) and (ii) in the statement. Proof in four
  steps; Step 4 gives two exact cases (free entry; θ̄ = 1), in both of which Δe_t = 0.
  Discussion: (i) holds automatically and (ii) binds at most 2.6% at the p_0 = 0 benchmark
  (0.7% at D1); Shimer (2005) for the quantitative nature of the s-shock comovement.
- **Part 2** (p_0 > 0, α → 0): to first order in F′(χ^c), under `eq:exit_cond` — det J > 0
  of the 3×3 system (χ^c, N, u) at t+1, which adds the exit → unemployment → profit loop the
  old 2×2 system omitted — and inflow-dominates-drain, where the drain includes the period-t
  exit term. `eq:gN_cond` (g_N < 0) was dropped: it fails at the calibration (g_N > 0); the
  left-hand side of `eq:exit_cond` is 0.002–0.005 there. ρ_s extension by continuity.
- **Part 3** (δ shock): negative comovement at h = 1 provided surviving period-t entry plus
  t+1 entry do not offset the destruction (`lem:vpre`, now a decomposition). Sufficient
  condition δ̄_e/τ̄ small. Discussion: at monthly frequency the first-month sign is fragile
  (entry offsets 81% at D1, more than 100% at 0.21), but vacancies fall from month 2 and in
  every quarterly average at D1.
- **Remark** after the proposition: g(N,u); at standard calibrations g_N > 0 (option value
  dominates); `eq:exit_cond` needs no sign for g_N.

**Prop. 7** in old notes is not referenced in the draft; treat as withdrawn.

## Companion theory notes (uncited)

`Notes/Baseline_Blanchard_Kahn.md` and `Notes/AGS_Blanchard_Kahn.md`: BK holds under ε > 1 and
endogenous exit makes it easier to satisfy, so ε > 2 is a Prop. 5 requirement, not a
determinacy one (P3).
