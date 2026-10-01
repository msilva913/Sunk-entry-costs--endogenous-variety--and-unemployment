# Draft Status — `Draft/Draft.tex`
**Last updated:** October 1, 2026 (cleanup; timing, wage appendix, Prop. 5 restatement)

**Build:** compiles with 0 errors and 0 missing graphics. The 4 undefined references and 1
duplicate label (P2 in [`pending_tasks.md`](pending_tasks.md)) predate the reposting work.
Compile from inside `Data/Draft/` (the graphicspath leads with relative entries), two passes,
with `-synctex=1`. If an editor auto-compiles on save, the `.aux` can be truncated
mid-build; delete it and rebuild.

**Notation.** χ^c_t (not x^c_t); Λ_t ≡ F(χ^c_t); ς for the survival quantile. Reposting:
Λ_{r,t} = F(Q_t/α), option value Q^rep_t (`eq:Qrep`), shortfall written inline as Q − Q^rep;
the symbols Ψ, D, M are retired (clashes with ψ, D^int, 𝒟, m). 𝓡_t is the expected
discounted shortfall in the wage appendix. Label prefixes `eq:`/`tab:`/`fig:`/`sec:`/`app:`.

**Timing.** The draft is the reference model ([S10](decisions.md)); the code deviates in three
places until R8 ([`model_equations.md`](model_equations.md) "Code status"). Until R8 and the D1
cascade, every §5.3 number comes from the code's timing.

---

## Section status

| Section | Status | Outstanding |
|---|---|---|
| 1 Introduction (incl. §1.1 evidence) | ✅ | Prose says reposting makes vacancies "self-correct" (~L154, ~L308); qualify for partial reposting (R11) |
| 2 Environment | ✅ | — Convention A `eq:u_lom` and `eq:v_lom` with Λ_{r,t}; timing text says which stocks each shock reaches; `fig:Timing` notes that s hits only incumbent matches |
| 3 Equilibrium | ✅ | R11: (1−ϕ) on the shortfall in `eq:surplus_wage`/`eq:surplus_eq`; −s_t𝓡_t in `eq:wage_eq`; prose at the job creation condition and the wage-substituted surplus; Λ_r has no labelled equation |
| 3.x Limiting cases | ✅ | — |
| 3.x Shock transmission (`sec:mechanism`) | ✅ | δ̄_e/τ̄ ≈ 0.21 at L1716 → 0.087 in the D1 cascade |
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
| D Steady state | ✅ — P5: `X^c` lines at L3647/L3721 use the shortcut (L3647 also lacks p_0) |
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
| `prop:ds_asymmetry` (5) | ✅ Restated Oct 1 for costless reposting; see below |

### Proposition 5 (`prop:ds_asymmetry`, proof in `app:proof_ds`)

**Framing.** Part 1 is a benchmark: with p_0 = 0 there is neither endogenous exit nor costly
reactivation, and s shocks produce positive u–v comovement. Negative comovement must come from
δ shocks with a large enough role (Part 3) or from partial reposting. The paragraph after the
proposition says Parts 1–2 assume costless reposting and that partial reposting gives s shocks a
destructive component. The §5.3 "Choice of shocks" paragraph is qualified the same way.

- **Part 1** (p_0 = 0): an s_t shock raises u_{t+1} and v_{pre,t+1} by the same
  d = (1−δ_t)(1−u_t), so vacancies rise at h = 1 provided −∂e_{t+1}/∂s_t < d, for any ρ_s.
  Proof in four steps: (1) ∂u = d; (2) ∂v_pre = d, for any ρ_s and entry technology;
  (3) ∂v = d + ∂e = θ̄d + ū∂θ, so the condition is necessary and sufficient; (4) two exact cases
  under ρ_s = 0 and the Prop. 1 conditions (free entry: θ pinned, ∂v = θ̄d; θ̄ = 1: equal
  deviations propagate with tightness unchanged). The discussion cites Shimer (2005) for
  separation shocks giving positive u–v correlation only quantitatively, and reports the
  numerical check: entry offsets at most 1.5% of the inflow at the calibrated p_0 = 0 benchmark
  (`run_prop5_s_check.jl`, [`findings.md`](findings.md)).
- **Part 2** (p_0 > 0, α → 0, DS-CES, `eq:gN_cond`): the result extends if the reposting inflow
  (1−δ_t)Λ_{t+1}(1−u_t) weakly dominates the exit drain. ∂χ^c_{t+1}/∂s_t < 0 via the 2×2
  Jacobian of the joint (χ^c, N) system (`eq:gN_cond` ⟺ g_N < 0 ⟺ det J > 1). The proof
  closes with Part 1's Step 3, with ∂v_pre in place of d. Open: R3.
- **Part 3** (δ shock): u↑ and v↓ if entry does not fully offset pre-committed vacancy
  destruction (`lem:vpre`); sufficient condition δ̄_e/τ̄ small, since the entry cushion carries
  δ̄_e/(r+δ̄_e). The calibrated value 0.087 (D1) strengthens it; the draft still says 0.21 in
  four places (D1 cascade). Under Convention A, ∂u/∂δ_e = (1−s̄)(1−u) + f·u.
- **Remark** after the proposition: g(N,u) ≡ d^f(N,u) + f_eρ(N)/μ; `eq:gN_cond` is the
  condition for profit dilution to dominate the option-value rise.

**Prop. 7** in old notes is not referenced in the draft; treat as withdrawn.

## Companion theory notes (uncited)

`Notes/Baseline_Blanchard_Kahn.md` and `Notes/AGS_Blanchard_Kahn.md`: BK holds under ε > 1 and
endogenous exit makes it easier to satisfy, so ε > 2 is a Prop. 5 requirement, not a
determinacy one (P3).
