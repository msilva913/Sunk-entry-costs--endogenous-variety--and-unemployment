# Session Handout — October 11–12, 2026
**Branch:** `costly_vacancy_reposting` · **Last pushed commit:** `e81caf3` (handout and context; before it `523cedd` and MS's rebuild `8882a42`)
**Previous handout:** [`session_handout_20261009.md`](session_handout_20261009.md) (strategy, D12)

These two sessions rebuilt the theory in §3 and Appendices G–J, made the continuation-cost
aggregate X^c exact in draft and code, and tightened the prose of Sections 1–3. No full-model
runs. **No §5 number has changed** (the X^c fix moves them slightly; the R8 + D1 regeneration
covers that).

---

## 1. What was done

### Theory (R15 option (b), R16, R12 audit, R12c)
- **`sec:mechanism` rebuilt** around free entry. Order: the free-entry paragraph with the user
  cost K_t = x_m(r+δ_t)/(1+r) (`eq:K_user_cost`) → `lem:theta` (signs only; closed forms in
  `app:proof_fe`, `eq:c_delta`) → `prop:fe` (rigorous; no Γ in the main text) → `tab:fe_irf`
  (δ at ρ = 0.592, s at ρ = 0.874) → `lem:accounting` (old Prop. 5 Parts 1 and 3, any ξ) →
  `rem:finite_xi` (numerical; replaces the wrong R12c sufficient condition) → `prop:exit_s`
  (old Part 2).
  - Conditions: |c_x|(ρ_x − ω_u) > ℓ_x at h = 1; |c_x|[ρ_x(1−λ_u²)/(1+ρ_xλ_u) − ω_u] > ℓ_x along
    the path.
  - The discussion is informal and leads with the two-channel contrast: tightness responds 0.79
    vs 0.04, unemployment 0.036 vs 0.37 per percent of the δ and s shocks.
- **Limiting cases restructured (MS's sequencing).** The opening question (keep the destruction
  flow, abstract from entry and exit) → `prop:independence` → clone replacement (`prop:ags`,
  `rem:nesting`) → `prop:double_limit`, which justifies clone replacement as a limit. The single
  limit ε → ∞ fails because business formation requires profits.
- **Audit corrections.**
  - `prop:double_limit` (i) **was wrong**. It is now N̄ = z̄L̄(1−δ̄)/(f̄_e(r+δ̄)) plus the rule
    `eq:N_rule_limit`, and the stability clause is dropped.
  - `prop:ags`: the BK claim is corrected.
  - `prop:bgm_nest`: the proof is rewritten.
  - `prop:curves` is now a remark.
  - `prop:independence`: necessity is made exact.
  - `rem:nesting`: CK and GS set the flow vacancy cost to zero.
  - `app:loglin`: re-dated to S10, and the `eq:tau_ll` coefficient corrected.
  - Full record: [`draft_status.md`](draft_status.md) "Propositions".
- **`prop:equilibria` (R16)** now holds under `ass:matching`:
  - CRS;
  - (M1) f ≤ min{1, θ};
  - (M2) ε_f ∈ (0,1) and nonincreasing;
  - plus condition (C), ln Φ convex in ln θ.

  Proof: ln N_JC − ln N_RC is convex, so there are at most two crossings. **Oct 12:**
  - Φ (break-even recruiter compensation) is defined in the main text (`eq:Phi`).
  - (C) is interpreted: the elasticity of Φ does not fall as the market tightens.

  The SMCD footnote cites Silva (2026), SSRN 6956640, as a working paper (**do not mention the
  submission**). The low equilibrium is genuine under bounded matching (q ≈ 1, u ≈ 1/3).

### X^c correction (P5)
- **Rule:** use the partial expectation × the pre-draw count, or the conditional mean × the
  post-draw count, never the partial expectation × N_t.
- **Result:** X^c = (1−δ_{t−1})(N_{t−1}+N^e_{t−1})∫χdF = N_t·E(χ | χ ≤ χ^c) = Nψ_cωχ^c in steady
  state, with ω = 1−(1−p_0)(1−δ)/(1−δ_e).
- **Draft:** `eq:agg_fixed_costs`, `eq:exit_threshold_curve`, the steady-state appendix (N, L^e
  and L^c were also wrong), and `eq:profit_share_ss`.
- **Code:** `payer_share` in `steady_state.jl`; f[16] and the SS block in
  `run_solution_core.jl`; three mirrors.
- **Verified:**
  - SS residual 5e-14 (PATH B, D1) and 4e-15 (PATH A).
  - The Bellman cutoff matches to 1e-15.
  - PATH A is unchanged; under PATH B, x_c rises 0.07%.
- **Error origin:** see the memory note (prior mixture weights in the code; Sept 30 pairing in
  the draft).

### Prose
- Intro accuracy:
  - The LP numbers match §4 (+0.22 pp on impact, ≈ 1.71 pp after four to five years, vacancy
    rate −0.58 pp).
  - The instrument is built from deaths.
  - The roadmap is fixed.
- §2: reactivation is stated once (MS's paragraph).
- §3: the conditional mean is corrected and X^r compressed.
- The DMP paragraph is kept (MS).

### New programs (all in `Programs baseline/`, outputs `.txt` alongside)
`gs_finite_xi_entry_offset.jl`, `double_limit_check.jl`, `equilibria_convexity_check.jl` (needs
the SMCD repo as a sibling checkout). Inventory: [`pipeline.md`](pipeline.md).

### Commits (Oct 11–12)
`8a42ca6` R16 · `0642487` R12 audit · `1628f30` context · `e1a2ab0` rebuild · `2c5b78e` limiting
cases · `f5cce4d` free-entry discussion · `1e63cd7` Sections 1–3 pass · `3c58bb7` X^c draft ·
`50ce448` X^c code · `523cedd` `eq:Phi` / Assumption 1. All pushed to origin. This handout and
the context cleanup were committed as `e81caf3` (Oct 12).

---

## 2. Implications to carry forward
- **The reposting margin is necessary, not optional.** Under costless reposting, s shocks give
  positive u–v comovement at every persistence (`prop:fe` part 3 at free entry,
  `lem:accounting` at finite ξ). Model cor(u,v) = +0.995 cannot be fixed without D11/R8.
- **The δ result is non-monotone in ξ** (`rem:finite_xi`): δ̄*(ξ) = 3.8% / 2.7% / 1.6% per year at
  ξ = 1 / 2 / 10 (ρ_δ = 0.592).
- **Paper prospects (D12).** The repairs make the theory sound for a field journal (RED / JME /
  AEJ: Macro) but do not move it into JET's category. Every dynamic theorem still switches
  variety off.

---

## 3. Open items (owner → where tracked)

| Item | Owner | Where |
|---|---|---|
| D12 venue/structure; D2 PATH A vs B | **MS** | `decisions.md` |
| D10(a) adoption → intro LP sentence to h = 0–2 / h = 4–6 (P7) | MS, then edit | `pending_tasks.md` |
| Add `silva_2026`, `den_haan_2000` to `references.bib` (P6) | **MS** | `pending_tasks.md` |
| **R8 + D1 + PATH B regeneration**, then D11 go/no-go via R2 | next model task | critical path 0b |
| Three δ̄_e/τ̄ ≈ 0.21 claims (L3272, L4132, L4459) | with the regeneration | 0b |
| `prop:exit_s` and `lem:accounting` numerics after R8 | with the regeneration | 0b |
| R18: remark on selection and stability of the low steady state | theory | R18 |
| R17: nonlinear `prop:fe` (handle e ≥ 0) | optional theory | R17 |
| R12b: cite for ê^δ > 0 in `app:loglin` | after R8 | R12b |
| R11: X^r/Y in §5.2's Y^Gross/Y identity | with P4 | R11 |
| `role_of_ingredients.md` §2 rewrite | before R14 | R14 |
| P1 (§4.3), P3/P3b, `app:robustness`, `app:weighting_robustness` | writing | P-tasks |
| M2 free-entry bracket in `steady_state_checks.jl` (fails at `HEAD` too) | low | M2 |

---

## 4. Environment notes for another device
- **Julia** runs as `julia --startup-file=no` from `Data/`. Packages load. The SMCD repo must be
  cloned as a sibling (`GitHub/Soft-Min-Cobb-Douglas-Matching-Function`) for
  `equilibria_convexity_check.jl`.
- **LaTeX:**
  - Build from `Data/Draft/` with three `pdflatex -synctex=1` passes. **Never `bibtex`.**
  - Close TeXstudio before scripted edits.
  - If the `.aux` gets truncated, delete it and rebuild.
  - Expect two undefined references (`app:robustness`, `app:weighting_robustness`) and no
    duplicates.
- **Files:**
  - `Draft.tex` contains 124 NUL bytes (pre-existing), so use `grep -a`.
  - Several `.jl` and all context `.md` files use CRLF.
  - Script edits should normalize line endings on read and restore them on write. Use `raw""`
    strings for LaTeX in Julia.

## 5. Reading order
1. This handout.
2. `../CLAUDE.md` "Priority tasks right now".
3. [`pending_tasks.md`](pending_tasks.md) critical path.
4. [`decisions.md`](decisions.md) D12 and D2.
5. [`draft_status.md`](draft_status.md) before touching `Draft.tex`.
6. [`model_equations.md`](model_equations.md) "Code status" before starting R8.
