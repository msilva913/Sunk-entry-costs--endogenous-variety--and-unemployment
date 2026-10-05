# Pending Tasks
**Last updated:** October 5, 2026 (R13 done; Part 2 numerics wait on R8) · **Branch:** `costly_vacancy_reposting`

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
  `tab:calib_targets` row 1, and the **three remaining** δ̄_e/τ̄ ≈ 0.21 claims (now ≈ 0.087).
  Locate them by searching `Draft.tex` for `0.21` rather than by line number — as of Oct 5,
  2026 they sit at L3156 (§5.3), L4020 (`app:loglin`) and L4298 (`app:comparison_D`), but these
  drift with every edit — the figures previously recorded here were 15 lines stale. The fourth, in `sec:mechanism`, was resolved Oct 5 in the Prop. 5
  rewrite: it now leads with the measured 0.087 and uses 0.21 only as an explicit comparison.
  Re-run `run_prop5_s_check.jl` and `run_prop5_parts23_check.jl` after the regeneration; the
  bounds currently in the draft and [`findings.md`](findings.md) are **2.6 % / 0.7 %** for
  condition (ii) and **0.81** for the Part 3 entry offset (c3_v). Update if they move.

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
| R13 | ✅ **Prop. 5 restated for start-of-period shock observation — done Oct 5** (wording approved by MS). Part 1: accounting with period-t entry, conditions (i)–(ii). Part 2: `eq:gN_cond` replaced by `eq:exit_cond` (det J > 0 of the 3×3 exit–variety–unemployment system; `eq:gN_cond` failed at the calibration, g_N > 0), period-t exit term added to the drain, B_t no longer predetermined, monotonicity-in-ρ_s claim replaced by continuity, first order in F′. Part 3: period-t entry term; `lem:vpre` restated as a decomposition. Discussion: conditions (i)–(ii) very easily satisfied at p_0 = 0 (≤ 2.6%, 0.7% at D1); Part 3 h = 1 sign fragile at monthly frequency but vacancies fall from month 2 and in every quarterly average at D1. **Remaining:** redo the Part 2 numerics after R8 (`run_prop5_parts23_check.jl` is indicative for p_0 > 0 because of the code's exit dating) | ⏳ numerics after R8 |
| R14 | **Ingredient ablation** (see [`../Notes/role_of_ingredients.md`](../Notes/role_of_ingredients.md)): after R8, one-at-a-time ζ→0, ω_δ→0, α→0, ξ→∞, κ=0, σ=0 on one calibration; report Block M moments and δ IRFs. **Two cells added Oct 5 (note §2):** (a) **ξ→∞ under an *s* shock**, not only δ — Comparison D tested δ only, and CK's result concerns a position-destroying shock with no variety channel, which is our s margin; (b) **ξ→∞ jointly with α > 0**, since under free entry Q is pinned and entry may replace retired positions, undoing the reposting margin. If (b) bites, the ingredients are not separable and that pair must be reported together rather than one-at-a-time. **The markup is excluded** from the list: μ = 1 removes the profits sustaining a finite N, so it eliminates the object under study rather than switching off a channel, and under DS-CES ζ = μ−1 so ε cannot move one without the other. Also: decide whether to add the exit-on-productivity VAR moment (â_21 ≈ −1.05) to Block M so `dest_elast_target` is disciplined; fix the elasticity print in `run_prop5_parts23_check.jl` (divide by σ_z, not 100σ_z) | after R8 |
| R2 | **Verify α, σ_s, δ_e are separately identified** on simulated data. α and δ_e move cor(u,v) through the same channel; α and χ_m both move position destruction. The LD→v IRF should separate them (exit compares χ to χ^c, reactivation to Q). Otherwise fix α on a reported grid. Gates the D11 go/no-go | ❌ needs R8 |
| R3 | **Prop. 5 for partial reposting.** Parts 1–2 now hold for costless reposting only (Oct 1 restatement). Open: a result for α > 0, e.g. a threshold Λ_r* at which the s→v sign flips, and whether Part 2's claim that the drain rises in ρ_s survives | ❌ |
| R4 | **Revisit [D4](decisions.md)/S6**: LD→v as the identifying moment for α | after R2 |
| R5 | **Add to §5.2:** the GS and Shao-Silos comparison, the δ_e bracketing argument, and why three shocks are required (CK consistency) | ⚠️ do regardless |
| R6 | **Decide whether δ_e moves into Θ_e** with a prior on [3.2%, 10%]/yr. Only once the margin exists, or δ_e absorbs its blame | after R8 |
| R8 | **Restructure the code to the draft's timing and add the reposting margin.** Spec: [`model_equations.md`](model_equations.md) "Code status" and Block 2. (i) States are the pre-exit stocks after the end-of-(t−1) incidence (Ñ_t, Ẽ_t, Ṽ^u_t, Ṽ^s_t) plus z_t, δ_t, s_t (known at t); N_t, u_t, v_{pre,t} become same-period controls solved with χ^c_t; survival (1−δ_t)F(χ^c_{t+1}) in both value equations and laws of motion. (ii) Add α, Λ_{r,t} as a tracked control, modified f[22], X^r in f[16], 𝓡_t as a tracked variable with its own equation, and f[3] gaining −(1−ϕ)·s'·𝓡' (as in `eq:jcc_eq`; no J needed, no nested E); the wage equation f[7] gains −ϕ·s·𝓡. (iii) Exact X^c integral. (iv) Regression tests: α → 0 and p_0 = 0 reproduce Block 1; positions conserved up to exit and the non-reactivated mass. (v) **Re-specify Comparison B** (`run_solution_endog_exit.jl`, described in `Draft.tex` as p_0 = 0 and ω_δ = 0): under Block 1 the two restrictions are one experiment because there is no reactivation margin to lose, but once α exists p_0 = 0 removes exit **and** reactivation (settled Oct 1, see R9). The exit-only switch is p_0 > 0 with ω_δ = 0; B as currently coded becomes the both-margins-off run. Flagged Oct 5 while revising the post-Proposition-3 passage, which now names B as the exit-only comparison | ❌ **next model task** |
| R9 | **Simulate the dynamic payoffs:** amplification (procyclical Λ_r, reposting option in the job creation condition) and persistence (stock depletion; does reposting move the δ→u peak?). Use one-at-a-time switches, not baseline-vs-AGS: ζ→0 (variety), ω_δ→0 (exit, reactivation retained), α→0 (reactivation, exit retained). Implies a fourth mechanism comparison, an α→0 run | after R8 |
| R11 | **Finish the draft's reposting edits.** Wage, surplus, job creation prose, and steady-state appendix done Oct 1. `eq:Lambda_r` added and `def:equilibrium` completed Oct 1 (θ_t = v_t/u_t, matching rates, 𝓡_t; J_t eliminated by substituting the surplus at t+1). Remaining: intro prose that says reposting makes vacancies "self-correct" (~L154, ~L308); §5.2's Y^Gross/Y = 1 + X/Y + X^c/Y needs X^r/Y (with P4) | ⚠️ |
| R12c | **The proof of Prop 5 Part 3 uses the circular K-decomposition step.** `app:proof_ds` (Draft.tex around L5164) argues that the entry response vanishes because "the duration-shortening contribution to K-hat^delta carries coefficient deltabar_e/(r+deltabar_e)", which is the same reading of `eq:K_decomp` holding Q-hat fixed that R12b removed from `app:loglin`. **The conclusion looks right and is corroborated independently** -- the destruction term converges to (1-qbar)vbar + sbar(1-ubar) > 0, bounded away from zero, while every channel through which a delta shock moves the vacancy value vanishes with deltabar_e, and the threshold derivation in [`../Notes/predetermined_vacancies_and_the_CK_special_case.md`](../Notes/predetermined_vacancies_and_the_CK_special_case.md) reaches the same place from the vacancy law of motion. But the justification as written is not sound, and it is inside a proof rather than prose, so **rewrite it via `eq:K_jcc_ll`** (survival, surplus, congestion) rather than via the identity. Care needed: the proof differentiates with respect to the LEVEL of delta_t, so the "every loading is proportional to deltabar_e" shortcut does not transfer directly -- d(theta)/d(delta_t) is bounded, not vanishing. Not attempted Oct 5 because it needs Part 3's proof worked through properly | after R8 |
| R12b | **Circular derivation in `app:loglin` par.9 -- FIXED Oct 5, but the numerical claim still needs a cite.** `eq:K_decomp` and `eq:Q_ll_delta` were the same asset-pricing identity (phi_D = D * deltabar_e/(r+deltabar_e) identically), yet the draft signed K-hat from delta-hat_e holding Q-hat fixed and then signed Q-hat from that K-hat -- one equation used twice. The entry cushion is **congestion relief** (theta down => q up => higher dividend), not duration shortening; duration signs K only at free entry, where Q is pinned by the entry technology. New par.9a/9b/9c add `eq:K_jcc_ll`, `eq:K_freeentry`, `eq:theta_freeentry`, `eq:tau_loading`. **Still owed:** e-hat^delta > 0 is now stated as a numerical result from Comparison D rather than derived, so point it at a specific figure or a `mechanism_stats.jl` number. Also compute the deltabar_e/taubar vs deltabar_e/(r+deltabar_e) loading comparison in a program so the magnitudes can enter the text under N15 | after R8 |
| R12 | **Audit the remaining propositions for reposting and timing:** `prop:bgm_nest` (leans on position conservation; leak is zero at p_0 = 0, but verify. ⚠️ Its proof says the bracket in `eq:N_lom_eq` "equals N^e_{t−1}"; it equals N_{t−1} + N^e_{t−1}, and the N^e formula and the "G(Q_t) at Q_t = f_e/ν" phrase look garbled), `prop:equilibria`, `prop:curves` | ❌ |

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
- **Oct 5:** wage/𝓡 dependency audit: `eq:gdp` income identity (Π^int net of X, X^r), `eq:curve_jcc`, `eq:e_theta`, `eq:phi_calib`, `eq:destruction_curve` (p_0 atom).
- **Oct 5:** Prop. 5 restated under the timing convention (R13): statement, proof, `lem:vpre`,
  discussion; `run_prop5_s_check.jl` extended; new `run_prop5_parts23_check.jl`.
- **Oct 1 (late):** wage dependencies carried through the main text and steady-state
  appendix; `eq:repost_shortfall`, `eq:Lambda_r`, `eq:theta_eq` added; `def:equilibrium`
  completed with J_t eliminated; §2 timing prose corrected to start-of-period shock
  observation.
- **Oct 1:** wage appendix rederived (Convention A employment law, reposting term);
  `eq:u_lom` and its restatements moved to Convention A; timing text and `fig:Timing`;
  `def:equilibrium` initial conditions; Prop. 5 Part 1 restated and reproved, with
  `run_prop5_s_check.jl`; `model_equations.md` re-dated to the draft's timing.
