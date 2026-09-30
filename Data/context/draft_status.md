# Draft Status — `Draft/Draft.tex`
**Last updated:** September 30, 2026 (reposting-margin revision landed; Props 1 and 3 audited;
graphics path fixed. Previous update September 23, 2026; last full section-by-section pass
September 5, 2026)

**Reposting notation (settled Sept 28–30, 2026).** Reactivation rate `Λ_{r,t} = F(Q_t/α)`;
option value `Q_t^{rep}` (`eq:Qrep`). The shortfall is written **inline** as `Q − Q^{rep}` and
the mean-cost integral inline — the symbols **`Ψ`, `D`, `M` are retired** and must not be
reintroduced. Full equation set: [`model_equations.md`](model_equations.md) Block 2.

**Build state:** compiles with **0 errors, 0 missing graphics**. The 4 undefined references and
1 multiply-defined label below are pre-existing and predate the reposting work. Compile from
inside `Data/Draft/` — the graphicspath now leads with relative entries.

Notation: χ_t^c (not x_t^c), Λ_t ≡ F(χ_t^c), ς for the survival quantile (renamed from ζ,
which clashed with the variety taste parameter). Label prefixes `eq:`/`tab:`/`fig:`/`sec:`/
`app:` throughout. The June 7, 2026 commit renamed F_χ → F everywhere.

---

## Section status

| Section | Status | Outstanding |
|---|---|---|
| 1 Introduction (incl. §1.1 motivating reduced-form evidence) | ✅ complete | — |
| 2 Environment | ✅ complete | Reposting added Sept 30: `eq:v_lom` carries Λ_{r,t}; timing prose and `fig:Timing` describe the shared Stage-2 exit/reactivation node |
| 3 Equilibrium (household, recruiters, wages, retailers, entry/exit, aggregation, definition) | ✅ complete | Reposting added Sept 30: `eq:Qrep`, `eq:surplus_wage`, `eq:jcc` restructured to reference the surplus (single one-period E), `eq:agg_repost_costs` with derivation, `eq:rc`/`eq:gdp` carry X_t^r, `def:equilibrium` threaded (J_t, Λ_{r,t}, Q_t^{rep}, X_t^r). ⚠️ Λ_r still has **no labelled equation** — inline in the auxiliary list |
| 3.x Limiting cases | ✅ complete | — |
| 3.x Shock transmission and the δ–s asymmetry (`sec:mechanism`) | ✅ complete | δ_e/τ ≈ 0.21 claim at `Draft.tex:1481` → 0.087 in the [D1](decisions.md) cascade |
| 3.x Steady state | ✅ complete | — |
| 4.1 Instrument construction | ✅ complete | — |
| 4.2 LP specification | ✅ complete | — |
| 4.3 Results | ⚠️ mostly drafted | δ–LD correlation paragraph; an in-text joint-LP sentence; the section ends with a bare, malformed `\ref{app:diagnostics}` |
| 5.1 Relationship to Coles and Kelishomi | ✅ complete (compressed) | δ̄ figure updates in the [D1](decisions.md) cascade |
| 5.2 Calibration and estimation | ⚠️ design complete, contested | Documents PATH A while the figures use PATH B ([D2](decisions.md)); β̂ ↔ β(θ) mapping unstated ([D3](decisions.md)); `app:weighting_robustness` promised but unwritten |
| 5.3 Inspecting the mechanism | ✅ complete (Comparisons A–D) | ⏳ **Regeneration owed.** Every figure and all ~24 numbers were produced at δ_e/τ = 0.210; [D1](decisions.md) is settled at 0.087. Re-run the four runners and `mechanism_stats.jl`, then diff before editing. [D2](decisions.md) may force a second pass |
| 5.4 Posterior estimates (`sec:posterior`) | ❌ **empty stub** — `Draft.tex:3012–3015`, three comment lines (posterior-prior plots, estimates table, moments at the posterior mean) | Blocked on the whole estimation build |
| 6 Conclusion | ❌ not started — `Draft.tex:3016` is a bare `\section{Conclusion}` followed by `\newpage\appendix` | Also mislabelled `ref:conclusion`; see LaTeX problems below |

## Appendix status

| Appendix | Status |
|---|---|
| A Relation to literature (`tab:ingredients_comparison`) | ✅ |
| B Data sources (`app:data_sources`) | ✅ |
| C Instrument diagnostics (`app:diagnostics`) — Deaths/Closings by supersector, residual δ–LD correlation, LD LPs, joint LP, QU placebo, filter robustness | ✅ |
| D Steady state (`app:steady_state`) | ✅ |
| E Derivations — profit share/ψ, resource-constraint curve, labor share, wage determination, log-linearized propagation | ✅ |
| F Impulse responses: additional comparisons (`app:additional_comparisons`, incl. Comparison D) | ✅ |
| G Limiting cases — proofs of `prop:independence`, `prop:double_limit`, `prop:ags`, `prop:bgm_nest` | ✅ |
| H Proof of `prop:ds_asymmetry` (`app:proof_ds`) | ✅ |
| I Steady-state equilibria — proof of `prop:equilibria` + comparative statics | ✅ |
| `app:weighting_robustness` | ❌ promised in §5.2, does not exist |
| `app:robustness` | ❌ referenced, does not exist |

## Known LaTeX problems (`Draft.log`, re-verified September 23, 2026)

- Undefined: `sec:conclusion`, `app:robustness`, `app:weighting_robustness`, `eq:labor_C_N`
- Multiply defined: `eq:profit_share`
- Malformed `\ref` at the end of §4.3 and in the §5.2 footnote (`documented in\ref{...}`)

**`sec:conclusion` has a diagnosed one-line cause.** `Draft.tex:3016` reads
`\section{Conclusion}\label{ref:conclusion}` — the label is `ref:conclusion` while every
reference points at `sec:conclusion`. Renaming the label fixes it. The section body is still
empty, so the fix is cosmetic until §6 is written.

---

## Propositions — inventory

All proofs complete except where noted. Numbering follows the draft.

### Reposting audit — September 30, 2026

| Proposition | Status under partial reposting |
|---|---|
| `prop:independence` (1) | ✅ **Holds, conditions unchanged.** `p_0 = 0` collapses F to a point mass, so Λ_r = 1, Q^rep = Q, X_r = 0: reposting vanishes *with* exit. No fourth condition needed — a payoff of drawing reactivation costs from the same F. Necessity unaffected: Λ_r and Q^rep depend only on Q_t, and X_r flows labor→formation like X_t. **`J_t` added to the stated labor block**, matched by `eq:surplus_eq` |
| `prop:double_limit` (2) | ✅ Unaffected — same argument; it assumes `p_0 = 0` |
| `prop:ags` (3) | ✅ Holds, and **gained the clone-replacement motivation**, which was never in any draft version (checked `git log -S`); only `Notes/AGS_Blanchard_Kahn.md` named the device. Proof now states `p_0 = 0` is *assumed*: no-exit follows from either `p_0 = 0` or `f_e ≥ χ_m` (since `χ_t^c = f_e`), and only the first kills reactivation. Under the second, `eq:v_ags`/`eq:rc_limit`/`eq:jcc_ags` would need Λ_r, −X_r, and the shortfall |
| `prop:bgm_nest` | ⚠️ **Not yet audited.** Sets `p_0 = 0` so it likely survives, but it leans on position conservation, and `model_equations.md` f[23] notes conservation "now breaks by exactly the non-reposted mass." At `p_0 = 0` that leakage is zero. **Verify rather than assume** |
| `prop:equilibria`, `prop:curves` | ⚠️ Not yet audited |
| `prop:ds_asymmetry` (5) | 🔴 **Part 1 breaks** — see the inventory entry below and task R3 |

**AGS now differs from the baseline on three channels, not two** (variety, endogenous exit,
*and* reactivation, since clone replacement sets `p_0 = 0`). The draft states this at the AGS
paragraph together with the one-at-a-time switches: `ζ→0` variety, `ω_δ→0` exit with
reactivation retained, **`α→0`** reactivation with exit retained.

**`prop:independence`, `prop:double_limit`, `prop:ags`** — limiting cases; proofs in the
limiting-cases appendix.

**`prop:bgm_nest` (Nesting LR-BGM).** Three-block proof:
1. Setup: p_0 = 0 → F(0) = 0 → X_t^c = 0, Λ_t = 1, δ_{e,t} = δ_{t−1}
2. Firm LOM (`eq:N_lom_eq`): substituting Λ = 1 plus the gross-output identity
   Y^c = z(1−u)ρ/μ makes the parenthesized term N^e_{t−1}, the BGM entry flow
3. Business-formation Euler (`eq:N_euler_eq`): Λ_{t+1} = 1, X^c = 0 → μ-ratio form with
   ν^f = f_e ρ(N), i.e. BGM equity pricing
4. Resource constraint (`eq:gdp`, **not** `eq:rc`): under p_0 = 0, Y^c = C + X_t; netting
   X_t gives Y_t = C_t + ν^f N^e, BGM eq. (3)
5. Closing remark: the isomorphism holds for any realization of {1−u_t, X_t}, not that the
   equilibrium distributions coincide

Key point: the isomorphism is at the **GDP** level (`eq:gdp`), not gross output (`eq:rc`).
X_t nets out as an intermediate input on both sides.

**`prop:equilibria` (DS-CES steady-state equilibria).** Restated for DS-CES with a
five-step proof in its own appendix section.

**`prop:ds_asymmetry` = Proposition 5 (δ–s asymmetry).** The paper's centerpiece; three
parts, plus a Remark and `lem:vpre`, proved in `app:proof_ds`.

- **Part 1** (s shock, exogenous exit, p_0 = 0): for any ρ_s ∈ [0,1), u↑ and v↑ together at
  h = 1 — positive u–v comovement.
  🔴 **Broken by partial reposting, not merely made conditional (Sept 30, 2026).** The proof in
  `app:proof_ds` argues that at ρ_s = 0 the JCC's right-hand side at t+1 "is invariant to `s_t`
  at leading order." `eq:surplus_wage` now carries
  `−s_t E_t m(1−δ)Λ(Q_{t+1}−Q_{t+1}^{rep})`, which depends on `s_t` directly, so that step
  fails. A dangling `eq:jcc_wage` citation there was repointed so the draft compiles, but the
  **argument is not fixed**. Locate it by searching `Draft.tex` for `invariant to $s_t$`.
  Part 1's conclusion should become conditional on a threshold Λ_r*. This is task **R3**, and it
  is more urgent than the September 30 session-A handout assumed.
- **Part 2** (s shock, endogenous exit, p_0 > 0, DS-CES): if `eq:gN_cond` holds —
  f_e N̄ < z̄(1−ū)(ε−2)/(ε−1), which requires **ε > 2** and small entry costs — and the
  reposting inflow weakly dominates the endogenous-exit drain on pre-committed vacancies,
  then positive comovement for ρ_s ∈ [0, ρ̄_s).
- **Part 3** (δ shock): u↑ and v↓ — negative comovement — if entry does not fully offset
  pre-committed vacancy destruction. Sufficient condition: δ̄_e/τ̄ small, so the entry-cushion
  coefficient δ̄_e/(r+δ̄_e) → 0. ⚠️ The quantitative claim is still stated at δ_e/τ ≈ 0.21 in
  four places (`Draft.tex:1481`, `:2787`, `:3646`, `:3921`). [D1](decisions.md) settled it at
  **0.087**, which *strengthens* the result. The retracted 0.116 has no standing.

Proof structure:
- `lem:vpre` (pre-committed vacancy destruction):
  ∂v_pre,t+1/∂δ_t = −Λ_{t+1}[(1−q(θ_t))v_t + s_t(1−u_t)] < 0; the bracket is predetermined
  with respect to δ_t.
- Part 1: two channels. Channel 1 (reposting) ∂v_pre/∂s_t = (1−δ_{e,t+1})(1−u_t) > 0.
  Channel 2 (entry), sub-step (a) pre-entry tightness θ^pre falls because u rises more than
  v_pre; sub-step (b) K_{t+1} rises because q(θ_{t+1}) > q̄, so Q_{t+1} > Q̄ and entry rises.
  Combined, ∂v_{t+1}/∂s_t > 0. A footnote notes that under free entry Channel 2 vanishes but
  Channel 1 survives.
- Part 2: signs ∂χ^c_{t+1}/∂s_t via the 2×2 Jacobian of the joint (χ^c, N) system. g_u < 0
  (higher u → lower profits → lower cutoff); `eq:gN_cond` ⟺ g_N < 0 ⟺ det J > 1; hence
  ∂χ^c/∂s_t < 0 and ∂Λ/∂s_t < 0. The drain rises in ρ_s, so dominance at ρ_s = 0 suffices.
- Part 3: uses `lem:vpre` for |∂v_pre/∂δ_t| > 0; the entry cushion from the K decomposition
  carries coefficient δ̄_e/(r+δ̄_e) → 0 as δ̄_e/τ̄ → 0.

Remark following the proposition defines g(N,u) ≡ d^f(N,u) + f_e ρ(N)/μ as the
continuation-profit threshold. g_N < 0 balances the dilution effect
(d^f ∝ N^{(2−ε)/(ε−1)}, net exponent negative for ε > 2) against the option-value rise
(f_e ρ(N)/μ ∝ N^{1/(ε−1)}); `eq:gN_cond` is exactly the condition for dilution to dominate.

**Prop 7** — listed as pending in the old notes, but nothing in the draft references it.
Treat as withdrawn unless a use appears.

---

## Companion theory notes (not yet cited in the draft)

`Notes/Baseline_Blanchard_Kahn.md` (June 6, 2026) and `Notes/AGS_Blanchard_Kahn.md`
(June 4, 2026) derive Blanchard-Kahn conditions for the full baseline and the AGS σ = 0 case.
Headline: BK holds under **ε > 1** — so ε > 2 is a Proposition 5 requirement, not a
determinacy requirement — and endogenous exit only makes BK easier to satisfy. Worth a
footnote or a short appendix subsection; see [`estimation_design.md`](estimation_design.md)
§Determinacy for why it also matters to the sampler.
