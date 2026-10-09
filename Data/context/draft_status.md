# Draft Status — `Draft/Draft.tex`
**Last updated:** October 9, 2026 (`prop:equilibria` proof gap → R16; Prop. 5 vs `prop:fe` division of labor; D12). Oct 7: R15 candidate restructured; Oct 6: Prop. 5 Part 3 qualified in ξ; Oct 5: Prop. 5 restated; Oct 1: wage dependencies, `def:equilibrium`, timing prose)

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
| 1 Introduction (incl. §1.1 evidence) | ✅ | Prose says reposting makes vacancies "self-correct" (~L154, ~L308); qualify for partial reposting (R11). **Oct 5:** endogenous exit no longer claimed as a recession amplifier (contradicted by §5.3 — Flag 1), reworded to the δ_e decomposition plus the host distribution F for reactivation draws (`ee5a2a5`); new paragraph giving the measurement argument for why the business formation block is needed given GS (`bcee2cb`); Fact 1 gains the entry margin, the C7–C10 correlations and the one-for-one clause (`1e94613`). **Still owed:** a paragraph on the ξ/CK contrast built on the predetermined share (text drafted in session, not spliced), and the ξ sentence (gated on R14) |
| 2 Environment | ✅ | Convention A `eq:u_lom`; `eq:v_lom` with Λ_{r,t}. Timing prose (Oct 1, late) states the paper's convention: the aggregate state incl. s_t, δ_t is known at the start of t; only incidence is revealed at the end of t. `fig:Timing` unchanged (MS: consistent with the convention); it notes that s hits only incumbent matches |
| 3 Equilibrium | ✅ | Wage dependencies fixed Oct 1: 𝓡_t defined at `eq:repost_shortfall`; `eq:wage_eq` carries −s_t𝓡_t; `eq:surplus_wage` carries −(1−ϕ)s_t𝓡_t. Oct 5: prose notes s_t𝓡_t enters alongside K_t (an extra flow cost of a filled position); recruiter profits Π^int now net of X_t and X^r_t, which fixes the `eq:gdp` income identity (it was off by X even before reposting); `eq:phi_calib` gains −s𝓡. `eq:Lambda_r` defines the reactivation rate. `def:equilibrium` audited for completeness Oct 1: 4 equations for (N, C, χ^c, Y^c), 6 for (θ, K, Q, e, v, u) with `eq:theta_eq` added. The job creation condition `eq:jcc_eq` substitutes the surplus at t+1, so J_t and w_t drop out (`eq:surplus_eq` removed); 𝓡_t is carried with its own equation, so there is no nested expectation. Auxiliary list cites `eq:Lambda_r`, the matching rates, and 𝓡_t. Prop. 1 proof updated to match |
| 3.x Limiting cases | ✅ | **Oct 5:** post-Prop-3 passage rewritten (`2fe27da`). Four channels over AGS, not three — **profit dilution** added (Y^c/N in `eq:euler_bf`, distinct from variety, survives ρ ≡ 1). Second tier added: AGS does not price replacement, and AGS collapses measured and structural destruction (Λ = 1 ⇒ δ_e,t = δ_{t−1}). Ablation paragraph now names Comparison C only, glosses ω_δ at first use, and states that exit and reactivation must be switched separately rather than jointly through p_0. Prop. 2 discussion gains a knife-edge caveat. ⚠ **Commits §5.3 to three runs that do not exist**: novar-vs-AGS, α → 0, and a re-specified Comparison B (R8 clause v) |
| 3.x Shock transmission (`sec:mechanism`) | ✅ | Prop. 5 and its discussion rewritten Oct 5; the δ̄_e/τ̄ ≈ 0.21 claim here was removed |
| 3.x Steady state | ✅ | Oct 5: `eq:e_theta` gains the entry that replaces non-reactivated positions (implicit in θ when α > 0); `eq:curve_jcc` gains −s𝓡(θ); `eq:destruction_curve` corrected for the atom (was valid only for p_0 = 1) |
| 4.1–4.2 Instruments, LP specification | ✅ | — |
| 4.3 Results | ⚠️ | δ–LD paragraph; joint-LP sentence; malformed `\ref` at the end (P1) |
| 5.1 Relation to Coles and Kelishomi | ✅ | δ̄ figures in the D1 cascade |
| 5.2 Calibration and estimation | ⚠️ contested | Documents PATH A, figures use PATH B ([D2](decisions.md)); β̂ ↔ β(θ) unstated ([D3](decisions.md)); α absent; `app:weighting_robustness` unwritten |
| 5.3 Inspecting the mechanism | ✅ | ⏳ Regeneration owed (R8 + D1 + PATH B). Numbers are checked against `mechanism_stats.txt`, which is the authority. The "Choice of shocks" paragraph now says the comparisons assume costless reposting |
| 5.4 Posterior estimates (`sec:posterior`, L3310) | ❌ empty stub | Blocked on the estimation build |
| 6 Conclusion (L3314) | ❌ not started | Label is `ref:conclusion`; references use `sec:conclusion` |

## Appendix status

**`app:loglin` — rewritten Oct 5, 2026 (`9582da1`, `10f3b60`).** ¶9 previously signed the entry response by reading one equation twice: `eq:K_decomp` and `eq:Q_ll_delta` are the same asset-pricing identity (φ_D = D·δ̄_e/(r+δ̄_e) identically), so it cannot sign K̂ and Q̂ independently. New ¶9a derives the dividend from the JCC (`eq:K_jcc_ll`): survival (−), surplus (−), and **congestion relief** (+), the last being the actual entry cushion. ¶9b records that at free entry Q is pinned so the identity *does* sign K (`eq:K_freeentry`, `eq:theta_freeentry`). ¶9c derives why δ̄_e/τ̄ decides the sign of v̂ (`eq:tau_loading`). ê^δ > 0 is now stated as a numerical result from Comparison D rather than derived. `eq:v_ll`'s destruction impulse corrected from a hat to a tilde (level deviation, matching its coefficient v̄^pre). **Owed (R12b):** a figure or `mechanism_stats.jl` number for ê^δ > 0, and the threshold T from the predetermined-vacancies note belongs in ¶7

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
| I Steady-state equilibria (`prop:equilibria`) | ❌ **Proof gap found Oct 9 (R16).** Step 4 does not deliver "generically exactly two" (zero is even; existence is not shown; a tangency is codimension one). Step 5 uses the degenerate ε → ∞ limit with f_e fixed |
| `app:weighting_robustness`, `app:robustness` | ❌ referenced, do not exist |

---

## Propositions

| Proposition | Status |
|---|---|
| `prop:independence` (1) | ✅ Holds with its three conditions (σ = 0, ζ = 0, p_0 = 0). p_0 = 0 collapses F to a point mass, so Λ_r = 1, Q^rep = Q, X^r = 0: reposting vanishes with exit, so no fourth condition is needed. J_t is in the stated labor block |
| `prop:double_limit` (2) | ✅ Unaffected (assumes p_0 = 0) |
| `prop:ags` (3) | ✅ Holds; states the clone-replacement motivation and that p_0 = 0 is assumed (no-exit also follows from f_e ≥ χ_m, which would leave reactivation on). AGS differs from the baseline on three channels |
| `prop:bgm_nest` | ⚠️ Not audited for reposting/timing (R12). At p_0 = 0 the conservation leak is zero |
| `prop:equilibria` | ❌ Proof does not establish the statement (R16; see appendix I above). Not audited for reposting (R12) |
| `prop:curves` | ⚠️ Not audited (R12). Remark-level: the draft says it "follows directly from the expressions above" |
| `prop:ds_asymmetry` (5) | ✅ Restated Oct 5 under the timing convention (R13); Part 2 numerics wait on R8 |

### Candidate: Beveridge response at free entry (R15) — **not in the draft**

Drafted in [`../Notes/beveridge_free_entry_GS.tex`](../Notes/beveridge_free_entry_GS.tex), the authoritative source. Stated in the AGS economy of `prop:ags` with σ = κ = 0, which by `rem:nesting` is Gabrovski–Silva, at ξ → ∞, with log shocks δ_t = δ̄·exp(δ̃_t) and s_t = s̄·exp(s̃_t), one at a time. **Structure as of Oct 7 (s-shock part added):** `lem:theta` (θ̂_t = a·δ̃_t + b·s̃_t, common denominator D, b = −ρ(1−δ̄)s̄/D; determinacy iff ϕ(1−δ̄)f̄ < η_L(2+r−τ̄), automatic when ϕ ≤ η_L), `rem:recursive` (block recursivity, BK root count, why finite ξ breaks it), `prop:fe` for x ∈ {δ, s} with R_x ≡ ℓ_x/|c_x| (ℓ_δ = f̄δ̄/τ̄, ℓ_s = (1−δ̄)f̄(1−δ̄/τ̄)): part 1, h = 1 iff ρ > ω + R_x; part 2, path slope θ̄[1 − ρ(1−λ²)/((1+ρλ)(ω+R_x))], negative iff ρ(1−λ²)/(1+ρλ) > ω + R_x; part 3, if (1−δ̄)f̄[η_L + ϕ(1−δ̄)f̄/τ̄] ≥ 1 an s shock gives v̂_t > 0 for all t ≥ 1 and a positive slope at every ρ; iid ⇒ slope θ̄ for both. `cor:allh` (below steady state at all horizons iff ρ − λ ≥ ω + R_x). All iff and all in closed form. The δ–s asymmetry is R_δ = 0.045 vs R_s = 9.56 at D1, ρ = 0.592. **Blocked on a design decision, not on work:** see R15 in [`pending_tasks.md`](pending_tasks.md). Option (b), the self-contained δ–s theorem, is now written in the note; the choice between (a) and (b) is MS's. Do not splice until that is settled, since option (b) changes what Prop. 5 is for.

**Division of labor with Prop. 5 (Oct 8, 2026).**
- `prop:fe` is sharper where both apply. It gives iff conditions in parameters and covers the
  whole path. It also shows Part 3's "δ̄_e/τ̄ small" fails at ξ → ∞ for an iid shock.
- Prop. 5 alone covers finite ξ (the calibrated case, where v_pre is a state and a small
  entrant weight ē/v̄ = δ̄(v̄+1−ū)/v̄ is what matters), endogenous exit (Part 2), and variety
  (Part 1 needs only p_0 = 0).
- The two bracket ξ.
- **Recommendation: option (b).** `prop:fe` leads. Part 2 stays as a standalone exit result.
  Parts 1 and 3 become a short general-ξ statement, with R12c fixed.

**Theory-bar assessment (Oct 9, [D12](decisions.md)).** At a theory journal `prop:fe` would read
as applied theory: log-linear, at the ξ → ∞ corner, with variety, exit and reposting off. A
global version is R17.

### Proposition 5 (`prop:ds_asymmetry`, proof in `app:proof_ds`) — restated Oct 5, 2026

**Part 3 qualified Oct 6, 2026.** The statement now reads "a sufficient condition, at a finite entry elasticity ξ, is that δ̄_e/τ̄ is small", and the discussion gains a paragraph explaining the restriction: as ξ → ∞ the entry cushion is no longer held back by the convexity of the posting cost and, for an iid shock, it overturns the destruction at every horizon after impact whatever δ̄_e/τ̄ (the iid statement of `prop:fe` in the R15 note). The draft's own numerics were never general in ξ — `run_prop5_parts23_check.jl` runs ξ_inv ∈ {0.5, 1, 2} only — but the text did not say so. The phrase "the free-entry case is treated separately below" is a placeholder; it becomes a cross-reference when R15 lands.

⚠ **Proof bug, logged as R12c (Oct 5).** The Part 3 proof argues that the entry response vanishes because "the duration-shortening contribution to K̂^δ carries coefficient δ̄_e/(r+δ̄_e)" — the same circular reading of `eq:K_decomp` that was removed from `app:loglin`. The conclusion looks right and is corroborated twice over (the destruction term is bounded away from zero; the threshold derivation in [`../Notes/predetermined_vacancies_and_the_CK_special_case.md`](../Notes/predetermined_vacancies_and_the_CK_special_case.md) reaches the same place from the vacancy law of motion), but the justification is unsound. Rewrite via `eq:K_jcc_ll`. Care: the proof differentiates w.r.t. the **level** of δ_t, where ∂θ/∂δ_t is bounded rather than vanishing, so the loading shortcut does not transfer

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
