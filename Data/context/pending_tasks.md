# Pending Tasks
**Last updated:** October 1, 2026 (cleanup) · **Branch:** `costly_vacancy_reposting`

Actionable work only. Decisions live in [`decisions.md`](decisions.md), draft section status in
[`draft_status.md`](draft_status.md), results in [`findings.md`](findings.md), the model
specification in [`model_equations.md`](model_equations.md).

The draft is theory-complete and empirics-complete. Section 5 (estimation) is the remaining
critical path, and the code must be brought in line with the draft (R8) before any number is
regenerated.

---

## Critical path

**0a. E9 — model-side LP test.** ⛔ Gates the Block B design ([D10](decisions.md)) and may
settle [D3](decisions.md). Run the baseline LP specification on a model-simulated 50-state
panel. If the model LP also builds monotonically, the full IRF path is a usable target;
if it peaks at h = 1–2, target short horizons. ⚠️ The simulated instrument must reproduce the
data instrument's persistence (0.91 quarterly), not the shock's (≈ 0.21 quarterly): use
industry-level δ processes with their own serial correlation and heterogeneous state
exposure. Templates: `simulate_model` (`solution_functions.jl`), `simulated_moments` in
`run_solution_delta_target.jl`.

**0b. One regeneration pass: R8 + D1 + PATH B.** ⛔ Three changes each invalidate every §5.3
number; do them together so the mechanism runs are regenerated once.
- **R8** (code restructuring and reposting margin), below.
- **D1:** `dest_ann = 0.0320` (still 0.0754 at `steady_state.jl:613`). Converges only under
  PATH B, so implement by switching the default path ([D2](decisions.md)).
- Then: re-run Comparisons A–D, regenerate the eight `mechanism_*.pdf`, re-run
  `mechanism_stats.jl`, diff `mechanism_stats.txt`, update §5.3. Also §5.2's external block,
  `tab:calib_targets` row 1, and the four δ̄_e/τ̄ ≈ 0.21 claims (now ≈ 0.087): `Draft.tex`
  L1716 (`sec:mechanism`), L3085 (§5.3), L3945 (`app:loglin`), L4223 (`app:comparison_D`).
  Re-run `run_prop5_s_check.jl` and update the 1.5% bound in the draft if it moves.

**1. Settle [D2](decisions.md) and [D3](decisions.md).** D3 after E9. Write the canonical
target set into `data_and_files.md`.

**2. `part5_wcrb.py` — wild cluster bootstrap → full Ω_β.** `part5_lp.py` saves only scalar
per-horizon SEs; `eq:ql_irf` needs the full covariance across horizons and outcomes
(dimension set by D10/E9). Persist the matrix (`omega_beta.npy`).

**3. `moments_bootstrap.py` — block bootstrap → Ω_m** (block length 8 quarters).

**4. Refresh `second_moments.jl` for Block M:** λ = 1,600 and all seven series
{u, v, s, f, δ, N^e, z}, in the row order of `tab:smm_moments`. Depends on M7.

**5. `model_irf.jl` — β(θ) extractor** in the LP's units and normalization (per D3).

**6. Objective, priors (including α), sampler.** `run_solution_core.jl` still has
`estimate = []`, `priors = (;)`. Evaluate the log-posterior once at the calibrated point and
confirm both quadratic forms are O(K) before sampling.

**7. Write §5.4 `sec:posterior`, `app:weighting_robustness`, §6 Conclusion**, plus the AGS
counterfactual.

---

## Reposting margin and timing

| # | Task | Status |
|---|---|---|
| R2 | **Verify α, σ_s, δ_e are separately identified** on simulated data. α and δ_e move cor(u,v) through the same channel; α and χ_m both move position destruction. The LD→v IRF should separate them (exit compares χ to χ^c, reactivation to Q). Otherwise fix α on a reported grid. Gates the D11 go/no-go | ❌ needs R8 |
| R3 | **Prop. 5 for partial reposting.** Parts 1–2 now hold for costless reposting only (Oct 1 restatement). Open: a result for α > 0, e.g. a threshold Λ_r* at which the s→v sign flips, and whether Part 2's claim that the drain rises in ρ_s survives | ❌ |
| R4 | **Revisit [D4](decisions.md)/S6**: LD→v as the identifying moment for α | after R2 |
| R5 | **Add to §5.2:** the GS and Shao-Silos comparison, the δ_e bracketing argument, and why three shocks are required (CK consistency) | ⚠️ do regardless |
| R6 | **Decide whether δ_e moves into Θ_e** with a prior on [3.2%, 10%]/yr. Only once the margin exists, or δ_e absorbs its blame | after R8 |
| R8 | **Restructure the code to the draft's timing and add the reposting margin.** Spec: [`model_equations.md`](model_equations.md) "Code status" and Block 2. (i) States x_t = [u_{t−1}, v_{t−1}, B_t, z_t, δ_{t−1}, s_{t−1}]; N_t, u_t, v_{pre,t} become same-period controls solved with χ^c_t; survival (1−δ_t)F(χ^c_{t+1}) in both value equations and laws of motion. (ii) Add α, Λ_{r,t} as a tracked control, modified f[22], X^r in f[16], the recruiter surplus as a tracked jump (no nested E), and the modified wage (−ϕ·s_t·𝓡_t; surplus coefficient 1−ϕ). (iii) Exact X^c integral. (iv) Regression tests: α → 0 and p_0 = 0 reproduce Block 1; positions conserved up to exit and the non-reactivated mass | ❌ **next model task** |
| R9 | **Simulate the dynamic payoffs:** amplification (procyclical Λ_r, reposting option in the job creation condition) and persistence (stock depletion; does reposting move the δ→u peak?). Use one-at-a-time switches, not baseline-vs-AGS: ζ→0 (variety), ω_δ→0 (exit, reactivation retained), α→0 (reactivation, exit retained). Implies a fourth mechanism comparison, an α→0 run | after R8 |
| R11 | **Finish the draft's reposting edits:** (1−ϕ) coefficient on the shortfall in `eq:surplus_wage` and `eq:surplus_eq`; −s_t𝓡_t inside the φ(·) bracket of `eq:wage_eq` (appendix already done); prose at the JCC (says the surplus is "substituted"; the RHS description is the old substituted form) and at the wage-substituted surplus (says "into `eq:jcc`"); give Λ_r a labelled equation; intro prose that says reposting makes vacancies "self-correct" (~L154, ~L308) | ⚠️ MS writes |
| R12 | **Audit the remaining propositions for reposting and timing:** `prop:bgm_nest` (leans on position conservation; leak is zero at p_0 = 0, but verify), `prop:equilibria`, `prop:curves` | ❌ |

## Empirical code tasks

| # | Task | Status |
|---|---|---|
| E1 | Re-run `part7b_sloos.py` against the current instrument (output from 2026-04-09, pre-BED-Deaths) | ⚠️ stale |
| E2 | Re-run the GFC diagnostic (`lp_irf_delta_gfc.csv` from 2026-05-12, pre-Deaths) | ⚠️ stale |
| E3 | `part7d_ld_gfc.py`: GFC_{t+h} outcome dummy, to test whether the GFC drives the monotone LD unemployment IRF | ❌ |
| E4 | Pre-GFC sample option (`max_qt="2007Q4"`) in `run_lp()` | ❌ |
| E5 | Refresh the BED cache to 2024Q4 (`refresh_bed_cache.py`, run locally) | ❌ |
| E6 | Run `part2b_residualize_shocks_v3.py` once, or delete the appendix promise | ❌ |
| E7 | `build_report_html.py` (needs pandoc) | optional |
| E9 | Model-side LP test — critical path item 0a | ❌ **high priority** |

## Paper writing tasks

| # | Task | Status |
|---|---|---|
| P1 | §4.3: δ–LD correlation paragraph, in-text joint-LP sentence, fix the malformed `\ref` at the section end | ⚠️ |
| P2 | Broken references: undefined `sec:conclusion` (label is `ref:conclusion` at L3314), `app:robustness`, `app:weighting_robustness`, `eq:labor_C_N`; duplicate `eq:profit_share` | ⚠️ |
| P3 | Cite the Blanchard-Kahn notes (`Notes/Baseline_Blanchard_Kahn.md`): BK under ε > 1, strengthened by endogenous exit | ❌ |
| P3b | GS's 6–10%/yr range vs δ_e = 3.2%: both sentences (L443, L468) are commented out. Check that no active claim remains; if the comparison is restored, explain why the added channels permit a lower δ_e | ⚠️ check |
| P4 | §5.2 rewrite to the D2 path, with α in the estimated block, `tab:calib_targets`, and priors; fix the superseded 7.54% δ_e row | after D2 |
| P5 | Steady-state appendix: `X^c = Nψ_cχ^c` (L3647) is missing p_0 and the (χ^c/χ_m)^ψ factor; L3721 uses the shortcut. Align with the exact `eq:agg_fixed_costs` | ⚠️ |

## Model tasks

| # | Task | Status |
|---|---|---|
| M2 | Free-entry steady state in `steady_state_checks.jl` lands on the θ = 1.80 branch; bracketed solve over `(log(0.1), log(0.6))`. Expositional only | low priority |
| M4 | `eval_SS` toolkit bug (`for ip in npar` iterates once), so callers need not pass a precomputed SS. Matters for the sampler | ❌ |
| M7 | 🔴 **Observable mapping for labor productivity.** Model `labor_prod = Y/(ρL)` has SD 0.0387 vs 0.0128 in data and correlates 0.9975 with N^e. Every RSD in Block M has it in the denominator. Candidates: Y_c/L_c or a different deflation; then an observable-mapping table in §5.2 | ❌ blocks Block M |
| M8 | Hump-shape gap: the model peaks at h = 1–2. Resolved through E9/D10, not a new mechanism; R9 checks whether reposting persistence matters | waits on E9 |
| M10 | ⚠️ Latent: `SS_symbolics` uses the gross Lerner share (μ−1) for f_e, i.e. the p_0 = 0 model. Inert because every runner passes `SS_precomputed`. Fix before reviving the symbolic SS path; never use it as a linearization point | dormant |

---

## Completed (details in git history and [`findings.md`](findings.md))

- **Sept 5–6:** non-steady-state linearization fixed (free-entry mismatch in `SS_numeric`;
  `solution_interface` now errors if the SS residual exceeds tolerance); f[16] posting-cost
  bug; `hp_filter` replaced by a true HP filter (model second moments before Sept 5 are
  invalid); Convention A restored in the code; `mechanism_stats.jl` (N15 for §5.3).
- **Sept 21–23:** s process calibrated from part6b (ρ_s = 0.874, σ_s = 0.0854); D1/M5
  diagnostic; E8 instrument-persistence test.
- **Sept 28–30:** reposting formulation (D11), job creation condition with reposting (R7), draft
  chunks A–C (R10), Props 1 and 3 audited.
- **Oct 1:** wage appendix rederived (Convention A employment law, reposting term);
  `eq:u_lom` and its restatements moved to Convention A; timing text and `fig:Timing`;
  `def:equilibrium` initial conditions; Prop. 5 Part 1 restated and reproved, with
  `run_prop5_s_check.jl`; `model_equations.md` re-dated to the draft's timing.
