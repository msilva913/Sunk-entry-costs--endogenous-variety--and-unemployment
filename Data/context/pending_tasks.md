# Pending Tasks
**Last updated:** October 19, 2026 (Q1 closed and merged; Q2 scoped, D14 opened). Before: October 18, 2026 (Q1 steps 5–7 done, close-out pending; R19 added: pricing replacement offsets the user-cost channel; BK ε-threshold located; `rem:dilution_variety` + App. K). Earlier: October 14, 2026 (quantitative plan Q0–Q6 adopted, D13; E9 superseded). Oct 13: P8 theory-section cleanup logged; R3 flagged as the main theory gap. Oct 12: cleanup after the Oct 11–12 theory sessions: R15, R16, R12, R12c, P5 done and compressed into §Completed). · **Branch:** `costly_vacancy_reposting`

Actionable work only. Decisions live in [`decisions.md`](decisions.md), draft section status in
[`draft_status.md`](draft_status.md), results in [`findings.md`](findings.md), the model
specification in [`model_equations.md`](model_equations.md).

The theory sections (§§1–3, Appendices G–J) were audited and corrected Oct 11–12. What remains:
bring the code in line with the draft (R8), regenerate every §5 number once, and write §5.4 and §6.

---

## Critical path — the quantitative plan adopted Oct 14, 2026 ([D13](decisions.md))

Bartik IRFs leave the objective and become validation. Block B′ is aggregate SVAR responses
matched by indirect inference. Full plan: `~/.claude/plans/wiggly-gliding-pelican.md`
(Oct 14). **Items 0a (E9), 2 and 5 below are superseded**; 0b is Q1.

| # | Task | Status |
|---|---|---|
| Q0 | **Presentation, in parallel.** (1) Notation audit of §§2–3 and the appendices; propose renames, apply after MS approves; add a notation table. Known clashes: e_t vs e_{zt}; ρ(N) vs ρ_x; ω vs ω_u vs ω_δ; Φ(θ) vs Φ_θ; η_L vs translog η; s in the SMCD footnote. (2) A "model at a glance" paragraph opening §2. **Done Oct 14:** Tier 1 renames plus the state index s → i (`25c2327`); "model at a glance" paragraph and notation table, Appendix K `app:notation` (`816b73b`). Left as is by choice: m(u,v) vs m_{t+1}, the λ family (HP, marginal utility, λ_u), the ε elasticity family | ✅ |
| Q1 | **R8 + D1 + PATH B, one regeneration pass.** Plan: [`../Notes/Q1_R8_D1_PATHB_plan.md`](../Notes/Q1_R8_D1_PATHB_plan.md). **Steps 1–7 done on `r8_timing_reposting` (Oct 14 and Oct 18); Tests 1–6 pass.** Step 5: `TARGETS` = D1 (`dest_ann = 0.0320`) + PATH B (`dest_elast_target = 5.0`), b/w = 0.9, x_v = 0.5, `Λr_target = nothing`; δ_e/τ = 0.0873 (`R8_step5_defaults.txt`); D2 settled. Step 6: `SHOCKS_BASE`, `shock_col`, qzflag guard; ten runners use `build_PAR` and named shock columns and all run; Comparison A's CK spec uses δ_e = τ(1−1e-6); legacy sensitivity runners marked unmaintained. Step 7: `R8_step7_baseline.txt` (cor(u,v) grid) and `R8_test6_free_entry.txt` (Test 6 passes under AGS emulation; gap 1.9e-4; found R19). **Closed Oct 19:** branches reconciled (`853b053`: cleaned core, `sigma_determinacy.jl`, `R8_test_timing` on `build_PAR`, JLD2 restored; all 8 R8 programs reproduce, all 10 runners pass); `model_equations.md` "Code status" resolved with the f[k] map; `parameters.md` Λ̄_r replaces α in Θ_e; `pipeline.md` conventions; D11 record; merged into `costly_vacancy_reposting` | ✅ |
| Q2 | **Trivariate (z, δ, s) VAR** for the shock processes. Extend `part6b_shock_persistence_cyclical.py`; z → δ → s; HP-1600, 1992Q3–2019Q4; off-diagonals excluded from the exogenous processes as for â₂₁. Replaces the separate AR(1) for s (corr(η^δ, η^LD) ≈ 0.49). Output `.txt` for §5.2 (N15). **More urgent after Q1 step 7 (Oct 18):** with the part6b AR(1) s process, cor(u,v) reaches only +0.46 at Λ̄_r = 0.70; the s shock sets the sign, so the Λ̄_r required cannot be read before the VAR. **Scoped Oct 19 (see [D14](decisions.md)).** Gaps in the part6b approach: (1) σ_s is a univariate residual SD, never orthogonalized against z and δ; (2) the double-counting argument used for z → δ does not apply to s, because the model has no endogenous s, so orthogonalizing is a choice (D14), not a correction; (3) BED noise enters s with weight ≈ δ/τ, which builds a spurious negative s–δ correlation; (4) ρ_s absorbs any lagged z, δ → s feedback, which is never reported; (5) the shock processes must come from the same identification as Block B′. **Deliverables:** full A₁ with t-stats, residual correlation matrix, σ_s before and after orthogonalization, the share of s innovation variance explained by z and δ, monthly (ρ, σ); then rerun the Step-7 cor(u,v) grid and report σ(τ) at the new processes. Output `.txt` for §5.2 and an updated `SHOCKS_BASE`. Environment: `Data/environment.yml` (Oct 20). **Done Oct 20:** `part6c_shock_var_trivariate.py` → `shock_var_trivariate.txt/.csv`; `Q2_shock_processes_check.jl/.txt`. Results (findings.md "Q2"): the s innovation is nearly orthogonal (1.1 % explained), so concern (1) is rejected; ρ_s/σ_s fall to 0.8387/0.0773 because lagged δ predicts s (t = 2.81). **cor(u,v) is unchanged** (+0.483 at Λ̄_r = 0.70, +0.114 at 0.55; floor 0.5): the Beveridge gap is not a shock-process artifact. **Still owed:** (i) adopt (a) in `SHOCKS_BASE` (with `SHK` in `sigma_determinacy.jl`) after MS approval, and rerun the R8 programs whose outputs depend on it; (ii) §5.2 text. **Oct 20 (MS): S13 (end 2019Q4, cut before filtering; applied to part6b/6c) and S14 (window option B).** **Next:** build the S14 processes (a (z, s) VAR on 1951–2019 with s from τ before 1992; the δ equation on 1992–2019), then refresh `SHOCKS_BASE` and rerun the `Q2_*.jl` programs. Block M targets move to the same windows (pre-pandemic cor(u,v) ≈ −0.93)  **Oct 28: S14 processes built in `part6b_shock_persistence_cyclical.py` (`shock_calibration.txt`): monthly z (0.9121, 0.0128), δ (0.5824, 0.0663), s (0.8418, 0.0812). The D14 independence rule triggers on 1951–2019 (z explains 15.7 % of the s innovation, all from before 1985); D14 reopened, MS to choose. Then refresh `SHOCKS_BASE` and rerun the `Q2_*.jl` programs (needs MS approval).**  **Oct 28: τ fixed** (1994 CPS correction now ×1.1 after 1994, `observables.py`; `raw_data.pkl` s and jf only, option (b)); part6b rerun: s (0.8364, 0.0751) monthly; D14 rule still triggers (19.1 %). **Owed:** MS decides D14 and the jf timing bug (findings.md); recompute τ̄ and the τ/jf moments; align `part8` (1.16); move the FRED key out of `observables.py:24` and rotate it. | ⏳ D14, jf timing, then SHOCKS_BASE |
| Q3 | **Targets.** Block M: refresh `second_moments.jl` (λ = 1,600, all seven series) plus lead/lag cross-correlations of δ_e with u and v; Ω_m via `moments_bootstrap.py`. Block B′: a VAR in (z, δ, s, u, v), δ- and s-shock responses of u and v at h = 0–12, Ω by VAR bootstrap; model side = the same VAR on simulated data. Validation set: Bartik δ→u h = 0–2, δ→v h = 4–6, severity placebo, LD→v | ❌ after Q1, Q2 |
| Q4 | **Identification check** (extends R2). Jacobian of [M; B′] w.r.t. Θ_e = {b, ξ⁻¹, x_v, ω_δ, p_0, α, (σ)} at the calibrated point: rank, condition number, loadings; Monte Carlo recovery. Weakly identified parameters are fixed on a reported grid | ❌ after Q3 |
| Q5 | **Hypotheses H1–H5, then calibrated mechanism comparisons over a sensitivity grid.** H1 δ → negative u–v comovement and persistent u. H2 α > 0 necessary and sufficient for cor(u,v) ≈ −0.8. H3 business formation adds δ→u persistence beyond AGS. H4 endogenous exit is minor. H5 reactivation amplifies z and s, damps δ. Grid: ξ⁻¹ ∈ {0.5, 1, 2}, b/w ∈ {0.7, 0.8, 0.9}, ω_δ, and α bracketing cor(u,v). Comparisons A–D plus α → 0; `mechanism_stats.jl`. **Gate for D12 / simple-model-first:** agree a threshold first (e.g. the C-vs-AGS gap moves the δ→u half-life by < 25 % or cor(u,v) by < 0.1 ⇒ lead with the simple model) | ❌ after Q1, Q2 |
| Q6 | **Estimation.** Bayesian limited-information SMM / indirect inference over Θ_e, shocks fixed from Q2 (two-step), `eq:posterior` with B′ in place of B; RWMH. Full-information estimation deferred | ❌ after Q4 |

**0a. E9 — model-side LP test.** ⚪ *Superseded Oct 14 by D13.* ⛔ Gates the Block B design ([D10](decisions.md)) and may
settle [D3](decisions.md). Run the baseline LP specification on a model-simulated 50-state
panel. If the model LP also builds monotonically, the full IRF path is a usable target; if it
peaks at h = 1–2, target short horizons. ⚠️ The simulated instrument must reproduce the data
instrument's persistence (0.91 quarterly), not the shock's (≈ 0.21 quarterly): use
industry-level δ processes with their own serial correlation and heterogeneous state exposure.
Templates: `simulate_model` (`solution_functions.jl`), `simulated_moments` in
`run_solution_delta_target.jl`.

**0b. One regeneration pass: R8 + D1 + PATH B.** ⛔ Each invalidates every §5.3 number; do them
together so the mechanism runs are regenerated once.
- **R8** (code restructuring and reposting margin), below.
- **D1:** `dest_ann = 0.0320` (still 0.0754 at `steady_state.jl:623`). Converges only under
  PATH B, so implement by switching the default path ([D2](decisions.md)).
- Then: re-run Comparisons A–D, regenerate the eight `mechanism_*.pdf`, re-run
  `mechanism_stats.jl`, diff `mechanism_stats.txt`, update §5.3. Also §5.2's external block,
  `tab:calib_targets` row 1, and the **three remaining** δ̄_e/τ̄ ≈ 0.21 claims (now ≈ 0.087).
  Find them by searching `Draft.tex` for `0.21` (use `grep -a`: the file has NUL bytes). As of
  Oct 12 they sit at L3272 (§5.3), L4132 (`app:loglin`) and L4459 (`app:comparison_D`). The
  hits at L1674 and L1784 (`sec:mechanism`) are correct as written.
- Re-run `run_prop5_s_check.jl` and `run_prop5_parts23_check.jl`. The bounds in the draft
  (`lem:accounting` discussion) are **2.6 % / 0.7 %** for condition (ii) and **0.81** for the
  δ-part entry offset. Update if they move. Also redo the `prop:exit_s` numerics (the current
  program is indicative for p_0 > 0 because of the code's exit dating).

**1. Settle [D2](decisions.md) and [D3](decisions.md).** D3 after E9. Write the canonical
target set into `data_and_files.md`.

**2. `part5_wcrb.py` — wild cluster bootstrap → full Ω_β.** `part5_lp.py` saves only scalar
per-horizon SEs; `eq:ql_irf` needs the full covariance across horizons and outcomes (dimension
set by D10/E9). Persist the matrix (`omega_beta.npy`).

**3. `moments_bootstrap.py` — block bootstrap → Ω_m** (block length 8 quarters).

**4. Refresh `second_moments.jl` for Block M:** λ = 1,600 and all seven series
{u, v, s, f, δ, N^e, z}, in the row order of `tab:smm_moments`. Depends on M7.

**5. `model_irf.jl` — β(θ) extractor** in the LP's units and normalization (per D3).

**6. Objective, priors (including α), sampler.** `run_solution_core.jl` still has
`estimate = []`, `priors = (;)`. Evaluate the log-posterior once at the calibrated point and
confirm both quadratic forms are O(K) before sampling.

**7. Write §5.4 `sec:posterior`, `app:weighting_robustness`, `app:robustness`, §6
Conclusion** (`ref:conclusion` exists as an empty stub), plus the AGS counterfactual.

---

## Reposting margin and theory

| # | Task | Status |
|---|---|---|
| R8 | **Restructure the code to the draft's timing and add the reposting margin.** Spec: [`model_equations.md`](model_equations.md) "Code status" and Block 2. (i) States are the pre-exit stocks after the end-of-(t−1) incidence (Ñ_t, Ẽ_t, Ṽ^u_t, Ṽ^s_t) plus z_t, δ_t, s_t (known at t); N_t, u_t, v_{pre,t} become same-period controls solved with χ^c_t; survival (1−δ_t)F(χ^c_{t+1}) in both value equations and laws of motion. (ii) Add α, Λ_{r,t} as a tracked control, modified f[22], X^r in f[16], 𝓡_t as a tracked variable with its own equation, and f[3] gaining −(1−ϕ)·s'·𝓡' (as in `eq:jcc_eq`); the wage equation f[7] gains −ϕ·s·𝓡. (iii) ~~Exact X^c~~ done Oct 11 (P5). (iv) Regression tests: α → 0 and p_0 = 0 reproduce Block 1; positions conserved up to exit and the non-reactivated mass. (v) **Re-specify Comparison B** (`run_solution_endog_exit.jl`): once α exists, p_0 = 0 removes exit **and** reactivation. The exit-only switch is p_0 > 0 with ω_δ = 0; B as currently coded becomes the both-margins-off run. The post-Prop-3 passage already names B as the exit-only comparison. ⚠ §5.3 prose is committed to three runs that do not exist yet (novar-vs-AGS, α → 0, re-specified B); check whether any script implements `prop:ags` as stated | ❌ **next model task** |
| R2 | **Verify α, σ_s, δ_e are separately identified** on simulated data. α and δ_e move cor(u,v) through the same channel; α and χ_m both move position destruction. The LD→v IRF should separate them. Otherwise fix α on a reported grid. Gates the D11 go/no-go | ❌ needs R8 |
| R3 | **Theory for partial reposting.** `lem:accounting` and `prop:exit_s` assume costless reposting. Open: a result for α > 0, e.g. a threshold Λ_r* at which the s→v sign flips. Note (Oct 11): under costless reposting the theory implies s shocks give positive comovement at every persistence, so D11/R8 is *necessary* to fix the +0.995 cor(u,v). **Oct 13 review: this is now the main gap in the theory.** The draft names partial reposting as one of the two sources of negative comovement, yet no result covers it. Suggested route: infinitely elastic vacancy creation (ξ → ∞, σ = κ = 0). There the labor market block stays block recursive, so derive a threshold Λ_r* (or α*) at which the s-shock path slope changes sign, in the style of `prop:fe` | ❌ |
| R4 | **Revisit [D4](decisions.md)/S6**: LD→v as the identifying moment for α | after R2 |
| R5 | **Add to §5.2:** the GS and Shao-Silos comparison, the δ_e bracketing argument, and why three shocks are required (CK consistency) | ⚠️ do regardless |
| R6 | **Decide whether δ_e moves into Θ_e** with a prior on [3.2%, 10%]/yr. Only once the margin exists | after R8 |
| R9 | **Simulate the dynamic payoffs:** amplification (Λ_r moves with Q, whose cyclicality is shock-specific: report z, s and δ separately, since reposting may amplify z and s shocks and damp δ shocks; reposting option in the job creation condition) and persistence (stock depletion; does reposting move the δ→u peak?). One-at-a-time switches: ζ→0 (variety), ω_δ→0 (exit, reactivation retained), α→0 (reactivation, exit retained). Implies a fourth mechanism comparison, an α→0 run | after R8 |
| R11 | **Last reposting edit:** §5.2's Y^Gross/Y = 1 + X/Y + X^c/Y (`eq:calib_stage2` text, ~L2955) needs X^r/Y (with P4). The intro "self-correct" prose was fixed Oct 11–12; the L350 footnote ("partially self-corrects") is accurate | ⚠️ |
| R12b | **`app:loglin` ¶9 numerical claim needs a cite.** ê^δ > 0 is stated as a numerical result from Comparison D; point it at a figure or a `mechanism_stats.jl` number. Also compute the δ̄_e/τ̄ vs δ̄_e/(r+δ̄_e) loading comparison in a program (N15) | after R8 |
| R14 | **Ingredient ablation** (see [`../Notes/role_of_ingredients.md`](../Notes/role_of_ingredients.md)). ⚠ **Rewrite that note's §2 first:** variety is *not* what makes the δ sign survive vacancy free entry (`prop:fe` shows the user-cost channel does). After R8: one-at-a-time ζ→0, ω_δ→0, α→0, ξ→∞, κ=0, σ=0 on one calibration; report Block M moments and δ IRFs. Cells: ξ→∞ under an *s* shock; ξ→∞ jointly with α > 0 (under vacancy free entry, entry may replace retired positions). The markup is excluded (μ = 1 removes the profits that sustain N). Also: decide on the exit-on-productivity VAR moment (â_21 ≈ −1.05) for Block M; fix the elasticity print in `run_prop5_parts23_check.jl` (divide by σ_z, not 100σ_z) | after R8 |
| R17 | **Test a nonlinear `prop:fe` (hypothesis).** At ξ → ∞, Q = x_m is pinned and the job creation condition has no endogenous states, so with Markov (δ, s) the exact policy should be θ_t = Θ(δ_t, s_t) globally. Steps: existence and uniqueness of Θ (contraction); monotonicity in δ and s; signs of u, v at h = 1 and along the path. Check numerically against `gs_free_entry_beveridge.jl` first. **Constraint:** vacancy free entry binds only while e_t > 0. At D1 ê_0 ≈ −10.9 per 1% destruction shock (κ = 0), so entry hits zero for shocks above ≈ 9% of δ̄ (9.3% at x_v = 0.5, `beveridge_free_entry.txt`), about 1.4 σ_δ; a global Θ must allow Q_t < x_m there. 4–8 weeks, medium risk | ❌ optional |
| R18 | **Remark on selection and stability of the low steady state** of `prop:equilibria` (genuine under bounded matching: q ≈ 1, u ≈ 1/3). Which one the calibration selects and why; local stability of each | ❌ |
| R19 | **Pricing replacement offsets the vacancy user cost (found Oct 18, Test 6).** With product-line entry, the entry Euler β(λ′/λ)(1−δ_t)(ν′+d′) = ν pins the survival-adjusted discount factor m(1−δ_t). The vacancy user cost K_t = Q_t − m(1−δ_t)Q_{t+1} uses the same factor, so with Q pinned a δ shock barely moves K or θ: the real rate falls with δ. θ̂/δ̃ ≈ −0.003 (ρ ≡ 1, σ = 1, ξ_inv = 1e-4) against c_δ = −0.72 under clone replacement; the baseline δ→θ response is ≈ −0.06 per unit at month 12. `prop:fe` stays correct as a statement about AGS, but its mechanism is nearly switched off once business formation is priced. **Steps:** (a) **double-limit check first (cheap):** does the offset survive ε → ∞ with f_e = f̄_e/ε (`prop:double_limit`), where μ → 1 and the entry Euler degenerates? If it vanishes, AGS is still the right limit and the gap is about the order of limits; if not, the AGS–baseline gap is real; (b) a short remark after `prop:fe`, the dynamic counterpart of "AGS does not price replacement"; (c) revisit how `sec:mechanism` and D12 frame `prop:fe`, and which channels carry the δ sign in the baseline (finite-ξ predetermined stock, variety). Evidence: `R8_test6_free_entry.txt`, first block | ❌ **next theory task** |

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
| P2 | Two undefined references remain in `Draft.log`: `app:robustness`, `app:weighting_robustness` (both appendices unwritten; critical path item 7). No duplicate labels | ⚠️ |
| P3 | Cite the Blanchard-Kahn notes (`Notes/Baseline_Blanchard_Kahn.md`): BK under ε > 1, strengthened by endogenous exit | ❌ |
| P3b | GS's 6–10%/yr range vs δ_e = 3.2%: both sentences are commented out. Check that no active claim remains; if restored, explain why the added channels permit a lower δ_e | ⚠️ check |
| P4 | §5.2 rewrite to the D2 path, with α in the estimated block, `tab:calib_targets`, and priors; fix the superseded 7.54% δ_e row | after D2 |
| P6 | **`references.bib` (MS):** add `silva_2026` (SSRN 6956640) and `den_haan_2000` (doi 10.1257/aer.90.3.482). Both were added by hand to `Draft.bbl`; the next bibtex run drops them otherwise | ⚠️ owed |
| P7 | If D10(a) is adopted, switch the intro LP sentence (+0.22 pp on impact, ≈ 1.71 pp after four to five years) to the robust h = 0–2 and h = 4–6 features | after D10 |
| P8 | **Theory-section cleanup (Oct 13 review).** All in `Draft.tex` §3 and its appendices; plan at `~/.claude/plans/wiggly-gliding-pelican.md`.<br>• **B1** New `ass:vacancy_fe` (infinitely elastic vacancy creation) collecting Q = x_m, e_t > 0, 0 < δ̄ ≤ τ̄, (1−δ̄)f̄ + τ̄ < 1. `prop:fe`'s proof uses λ_u ≥ 0, which its hypotheses did not contain.<br>• **B2** `rem:finite_xi`: say the switch to "every δ̄" lies between ξ = 10 and 100; call 1 %/yr to τ̄ the calibratable range; label the 85 % (GS economy, ξ = 1) and 81 % (full model, p_0 = 0) offsets.<br>• **B3** Rename the resource-constraint bound N̄ to N^max (it clashed with Prop. 3's N̄). "N_RC → N̄" is false; the right statement is N_RC < N^max.<br>• **B4** SMCD footnote: "that equilibrium" → "the low equilibrium".<br>• **B5** Steady-state figure captions explain mechanisms (N6); trim them.<br>• **B6** CK pinpoint cite (p. 125, "analogous to the free entry case").<br>• **B7** `app:comparison_D`: "duration shortening" → user cost.<br>• **B8** `rem:nesting` consumption identity omits x_m^{−ξ}.<br>• **B9** Limiting cases: "two respects" lists three; the comparison-design paragraph duplicates §5.3 and becomes a pointer.<br>• **B10** `prop:bgm_nest` → remark `rem:bgm_nest` (MS: it does connect the model to BGM, but the section is overloaded).<br>• **B11** `sec:mechanism` opener: step count; "slope negatively".<br>• **B12** Labels `prop:curves` → `rem:curves`, `ref:conclusion` → `sec:conclusion`.<br>• **B13** Terminology: never bare "free entry" for vacancies; say "vacancy free entry" or "infinitely elastic vacancy creation" (MS, Oct 13).<br>**Done Oct 13**: built with 0 errors and only the two known undefined references. New numbering: Assumption 1 `ass:vacancy_fe`, Assumption 2 `ass:matching`; Props 1–6 (fe = 4, exit_s = 5, equilibria = 6); Remarks 1–5 (nesting, bgm_nest, dmp, finite_xi, curves). Also fixed: §5.2 cited `prop:bgm_nest` for the CK nesting. Three new issues are logged in `draft_status.md` "Remaining known issues": the garbled §5.3 Comparison A sentence, the Comparison B description mismatch, and the claim at ~L902 that infinitely elastic vacancy creation nests DMP. The last was resolved the same day by the new `rem:dmp`: isomorphism with DMP when every position bears the flow user cost K_t, except where e_t ≥ 0 binds (MS's conjecture, verified algebraically) | ✅ |

## Model tasks

| # | Task | Status |
|---|---|---|
| M2 | Free-entry steady state in `steady_state_checks.jl` fails to bracket (also fails at `HEAD` before the P5 fix); bracketed solve over `(log(0.1), log(0.6))`. Expositional only | low priority |
| M4 | `eval_SS` toolkit bug (`for ip in npar` iterates once), so callers need not pass a precomputed SS. Matters for the sampler | ❌ |
| M7 | 🔴 **Observable mapping for labor productivity.** Model `labor_prod = Y/(ρL)` has SD 0.0387 vs 0.0128 in data and correlates 0.9975 with N^e. Every RSD in Block M has it in the denominator. Candidates: Y_c/L_c or a different deflation; then an observable-mapping table in §5.2 | ❌ blocks Block M |
| M8 | Hump-shape gap: the model peaks at h = 1–2. Resolved through E9/D10, not a new mechanism; R9 checks whether reposting persistence matters | waits on E9 |
| M10 | ⚠️ Latent: `SS_symbolics` uses the gross Lerner share (μ−1) for f_e, i.e. the p_0 = 0 model. Inert because every runner passes `SS_precomputed`. Fix before reviving the symbolic SS path | dormant |

---

## Completed (details in git history, [`findings.md`](findings.md), [`draft_status.md`](draft_status.md))

- **Sept 5–6:** non-steady-state linearization fixed; f[16] posting-cost bug; true HP filter
  (model second moments before Sept 5 are invalid); Convention A restored; `mechanism_stats.jl`.
- **Sept 21–23:** s process calibrated from part6b (ρ_s = 0.874, σ_s = 0.0854); D1/M5
  diagnostic; E8 instrument-persistence test.
- **Sept 28–Oct 1:** reposting formulation (D11), job creation with reposting (R7), draft
  chunks A–C (R10), wage appendix and `def:equilibrium` completed, `model_equations.md` re-dated.
- **Oct 5:** Prop. 5 restated under the timing convention (R13); `app:loglin` ¶9 circularity
  fixed; wage/𝓡 dependency audit.
- **Oct 6–7:** free-entry δ–s result (`Notes/beveridge_free_entry_GS.tex`, now superseded by the draft).
- **Oct 11 (R15, R16, R12, R12c):** `lem:theta` + `prop:fe` lead `sec:mechanism`; old Prop. 5
  split into `lem:accounting`, `rem:finite_xi`, `prop:exit_s`. `prop:equilibria` restated under
  `ass:matching` + (C). Every proposition in §3 and Appendices G–J audited and corrected
  (`prop:double_limit` (i) was wrong). Programs: `gs_finite_xi_entry_offset.jl`,
  `double_limit_check.jl`, `equilibria_convexity_check.jl`.
- **Oct 11 (P5):** X^c = N·E(χ | χ ≤ χ^c) = Nψ_cωχ^c, ω = 1−(1−p_0)(1−δ)/(1−δ_e), in the draft
  (`eq:agg_fixed_costs`, `eq:exit_threshold_curve`, steady-state appendix) and the code
  (`payer_share` in `steady_state.jl`; `run_solution_core.jl` f[16] and SS; three mirrors).
  SS residual 5e-14 (PATH B, D1); PATH A unchanged; PATH B x_c +0.07%.
- **Oct 11–12 (prose):** limiting cases led by clone replacement; free-entry discussion recast
  around the δ–s contrast with `tab:fe_irf`; Sections 1–3 efficiency and accuracy pass;
  `ass:matching` tightened with Φ defined in the main text (`eq:Phi`).
- **Oct 18:** Q1 steps 5–7 (D1 + PATH B defaults, D2 settled; runners converted; reactivation
  grid; Test 6). BK threshold located (ε* = 2.40; `bk_epsilon_threshold.jl`) with the loops
  behind it; `rem:dilution_variety` (Remark 3) + Appendix K (`app:determinacy`, `tab:eps_threshold`).
  `costly_vacancy_reposting` merged into `r8_timing_reposting` (`1a25965`).
