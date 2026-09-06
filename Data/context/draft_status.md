# Draft Status — `Draft/Draft.tex`
**Last updated:** September 5, 2026 (verified section by section against the `.tex` and
`Draft.log`; supersedes the "as of May 19" status that used to live in `findings.md`)

Notation: χ_t^c (not x_t^c), Λ_t ≡ F(χ_t^c), ς for the survival quantile (renamed from ζ,
which clashed with the variety taste parameter). Label prefixes `eq:`/`tab:`/`fig:`/`sec:`/
`app:` throughout. The June 7, 2026 commit renamed F_χ → F everywhere.

---

## Section status

| Section | Status | Outstanding |
|---|---|---|
| 1 Introduction (incl. §1.1 motivating reduced-form evidence) | ✅ complete | — |
| 2 Environment | ✅ complete | — |
| 3 Equilibrium (household, recruiters, wages, retailers, entry/exit, aggregation, definition) | ✅ complete | — |
| 3.x Limiting cases | ✅ complete | — |
| 3.x Shock transmission and the δ–s asymmetry (`sec:mechanism`) | ✅ complete | δ_e/τ ≈ 0.21 claim at `Draft.tex:1481` depends on [D1](decisions.md) |
| 3.x Steady state | ✅ complete | — |
| 4.1 Instrument construction | ✅ complete | — |
| 4.2 LP specification | ✅ complete | — |
| 4.3 Results | ⚠️ mostly drafted | δ–LD correlation paragraph; an in-text joint-LP sentence; the section ends with a bare, malformed `\ref{app:diagnostics}` |
| 5.1 Relationship to Coles and Kelishomi | ✅ complete (compressed) | δ̄ figure depends on [D1](decisions.md) |
| 5.2 Calibration and estimation | ⚠️ design complete, contested | Documents PATH A while the figures use PATH B ([D2](decisions.md)); β̂ ↔ β(θ) mapping unstated ([D3](decisions.md)); `app:weighting_robustness` promised but unwritten |
| 5.3 Inspecting the mechanism | ✅ complete (Comparisons A–D) | Figures must be regenerated if [D1](decisions.md) or [D2](decisions.md) change |
| 5.4 Posterior estimates (`sec:posterior`) | ❌ **empty stub** — three TODO comments | Blocked on the whole estimation build |
| 6 Conclusion | ❌ not started | `sec:conclusion` is referenced but undefined |

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

## Known LaTeX problems (`Draft.log`)

- Undefined: `sec:conclusion`, `app:robustness`, `app:weighting_robustness`, `eq:labor_C_N`
- Multiply defined: `eq:profit_share`
- Malformed `\ref` at the end of §4.3 and in the §5.2 footnote (`documented in\ref{...}`)

---

## Propositions — inventory

All proofs complete except where noted. Numbering follows the draft.

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
- **Part 2** (s shock, endogenous exit, p_0 > 0, DS-CES): if `eq:gN_cond` holds —
  f_e N̄ < z̄(1−ū)(ε−2)/(ε−1), which requires **ε > 2** and small entry costs — and the
  reposting inflow weakly dominates the endogenous-exit drain on pre-committed vacancies,
  then positive comovement for ρ_s ∈ [0, ρ̄_s).
- **Part 3** (δ shock): u↑ and v↓ — negative comovement — if entry does not fully offset
  pre-committed vacancy destruction. Sufficient condition: δ̄_e/τ̄ small, so the entry-cushion
  coefficient δ̄_e/(r+δ̄_e) → 0. ⚠️ The quantitative claim is stated at δ_e/τ ≈ 0.21; under
  [D1](decisions.md) it would become 0.116, which *strengthens* the result.

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
