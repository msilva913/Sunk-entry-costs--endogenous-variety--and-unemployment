# Pending Tasks
**Last updated:** September 5, 2026 · **Branch:** `Organize_Project_State_Estimation`

Actionable work only. **Decisions** (things to choose, not do) live in
[`decisions.md`](decisions.md); **draft section status** lives in
[`draft_status.md`](draft_status.md); **results** live in [`findings.md`](findings.md).

The paper is theory-complete and empirics-complete. Section 5 is the entire remaining
critical path.

---

## 🛑 STOP — the perturbation solutions are linearized around a non-steady state

**Discovered September 5, 2026. This precedes everything else on the critical path.**

`solution_interface` prints `Max SS residual` and then proceeds unconditionally — it never
asserts the residual is small (`run_solution_core.jl:85`). The residual is **not small**:

| runner / case | Max SS residual |
|---|---|
| `run_solution_entry_elasticity.jl` baseline (Comparison D, ξ_inv=1.0) | **14.468** |
| `run_solution_entry_elasticity.jl` low-ξ (Comparison D, ξ_inv=0.1) | **14.178** |
| `run_solution_delta_target.jl` BED / CODE / BGM | **22.074 / 14.468 / 12.209** |

**Root cause — it is the free-entry condition `f[18]`.** Verified to all printed digits:
residual = ν_f(SS) − ρ·f_e(cal)/μ.

`SS_numeric` (`run_solution_core.jl:367`) recomputes the steady state from scratch instead of
calling `steady_state(cal)`. It derives its own entry cost
`f_e_s = zbar_s·(μ−1)·L_s/N_s·(1−δ_e)/(δ_e·μ + r)` and sets `ν_f_s = ρ·f_e_s/μ` — but the
`PAR` vector handed to the solver carries `cal.f_e` from `calibrate_shares`. The two disagree
by roughly 2.4× (CODE case: implied f_e = 32.5 in `SS_numeric` vs `cal.f_e` = 13.66), so free
entry fails by exactly their difference.

`steady_state.jl` itself is **fine** — `steady_state_checks.jl` asserts `f[16]`, `f[18]`,
`f[19]`, `f[20]`, the JCC, the Nash wage and both LOMs at 1e-12, and all pass under
`dest_ann` ∈ {0.0320, 0.0754, 0.0963}. The bug is only in the second, redundant SS
construction used as the linearization point.

**Scope — every runner that passes `SS_numeric` to `solution_interface`:**
`run_solution_all_delta.jl` (A), `run_solution_endog_exit.jl` (B), `run_solution_no_variety.jl`
(C), `run_solution_entry_elasticity.jl` (D), `run_solution_delta_target.jl`. **All four
mechanism comparisons in §5.3 are affected**, as are the D1/M5 diagnostic numbers in
`findings.md` (now retracted — see that file).

**Fix, in order:**
1. **Guard first**: ` SS_max < 1e-8` in `solution_interface` so this can never pass
   silently again. Two lines, and it will immediately fail every mechanism runner — which is
   the correct behavior.
2. **Replace `SS_numeric`** with a thin mapper that calls `steady_state(cal)` and orders its
   output into the SS vector. That construction is already verified. Deleting the duplicate
   removes the whole class of bug.
3. Re-run Comparisons A–D and regenerate `mechanism_*.pdf`.
4. Only then revisit D1/M5.

**LOM timing inconsistency — DECISION PENDING.** The u LOM (`f[23]`, `eq:u_lom`) uses a post-matching separation base while the v LOM (`f[22]`, `eq:v_lom`) uses a pre-matching reposting base. The mixture leaks job positions at rate (1-δ_e)·s·q·v. Gabrovski-Silva and the old Dynare model use the pre-matching convention consistently and conserve positions exactly, so the v LOM is the unchanged equation and the u LOM is what moved. Full analysis, both fix options, and the implications for GS: [`../Notes/LOM_timing_consistency.md`](../Notes/LOM_timing_consistency.md).

**Related bug — FIXED Sept 5, 2026.** `run_solution_core.jl` f[16] computed the vacancy
creation cost as `X = e·ξ_inv/(1+ξ_inv)·Q + κqv`. Correct form is `e/(1+ξ_inv)·Q + κqv`
(draft eq. 29, p. 20): sunk posting costs are the integral of the marginal schedule
Q = x_m·e^ξ_inv, giving ∫₀^e x_m u^ξ_inv du = e·Q/(1+ξ_inv); the second term is fixed
matching costs. Provenance confirmed by MS: it came from the original ξ/(ξ+1) transcribed
symbol-for-symbol with ξ→ξ_inv rather than inverted (ξ = 1/ξ_inv ⇒ ξ/(ξ+1) = 1/(1+ξ_inv)).
`steady_state.jl:362,742` and `SS_numeric` (line 439) always had it right, so this was an
inconsistency between f[16] and every other use. The two expressions coincide at ξ_inv = 1,
so only **Comparison D's ξ_inv = 0.1 case** was affected — by 10× in the sunk-cost term.

---

## Critical path to a complete draft

Steps 1–3 are gates; nothing downstream is worth doing until they clear.
**But fix the steady-state bug above first — none of these numbers mean anything until then.**

**1. Settle [D1](decisions.md) (which δ̄_e), [D2](decisions.md) (PATH A or B), and
[D3](decisions.md) (β̂ ↔ β(θ) scaling).** ⛔
Then write the winning target set into `data_and_files.md` as canonical, and regenerate any
mechanism figure affected by D1/D2.

**2. `part5_wcrb.py` — wild cluster bootstrap → full Ω_β.** ⛔
`part5_lp.py` runs each horizon as a separate regression and saves only a scalar clustered
`se`. `eq:ql_irf` needs the **42×42** covariance across h = 0..20 and across both outcomes;
LP coefficients are strongly correlated across horizons by construction. Persist the matrix
(`omega_beta.npy`), not bootstrap SEs. n = 50 clusters is at the lower bound of asymptotic
reliability, so this doubles as the paper's inference robustness.

**3. `moments_bootstrap.py` — block bootstrap → Ω_m.** ⛔
Block length 8 quarters. Nothing currently computes the covariance of the empirical moments.

**4. Refresh `second_moments.jl` for Block M.** ⛔
Currently λ = 100,000 and only 8 series. Needs λ = 1,600 and all seven series
{u, v, s, f, δ, N^e, z}, returning m(θ) in exactly the row order of `tab:smm_moments`.
Validate against the empirical table at the calibrated point.

**5. `model_irf.jl` — β(θ) extractor.** ⛔
Map the state-space solution to a model IRF in the LP's units and normalization (per D3):
1-SD δ shock, u and v in pp, quarterly, h = 0..20, level difference vs. t−1.

**6. Objective + priors, then the sampler.** ⛔
`run_solution_core.jl:112` still reads `estimate = []`; `priors = (;)`. Evaluate the
log-posterior kernel once at the calibrated point and confirm both quadratic forms are
O(K) before sampling.

**7. Write §5.4 `sec:posterior`, `app:weighting_robustness`, and §6 Conclusion.** ⛔
Plus the AGS two-model counterfactual.

---

## Empirical code tasks (secondary)

| # | Task | Status |
|---|---|---|
| E1 | **Re-run `part7b_sloos.py`** against the current instrument — existing output is from 2026-04-09 with `instr_sd`=11.16 (pre-BED-Deaths) | ⚠️ stale output |
| E2 | **Re-run the GFC diagnostic** — `lp_irf_delta_gfc.csv` is from 2026-05-12, also pre-Deaths | ⚠️ stale output |
| E3 | **`part7d_ld_gfc.py`** — GFC_{t+h} outcome-quarter dummy, to test whether the GFC drives the monotone LD unemployment IRF (inconsistent with ρ_LD ≈ 0.3) | ❌ never written |
| E4 | **Pre-GFC sample restriction** — add a `max_qt="2007Q4"` option to `run_lp()` in `part5_lp.py`; same question as E3 from the other side | ❌ |
| E5 | **Refresh the BED cache to 2024Q4** — run `refresh_bed_cache.py` locally (BLS rate limits block it in a sandbox); the `ext_2024` sample in part11 currently truncates at 2021Q4, and `raw_data.pkl` δ ends 2021Q4 | ❌ |
| E6 | **Run `part2b_residualize_shocks_v3.py` once**, or delete the appendix promise from the draft — see [D8](decisions.md) | ❌ |
| E7 | **`build_report_html.py`** — needs `conda install -c conda-forge pandoc` locally | ❌ optional |

## Paper writing tasks (secondary)

| # | Task | Status |
|---|---|---|
| P1 | **§4.3** — add the δ–LD correlation paragraph and an in-text joint-LP sentence; fix the malformed `\ref` that ends the section | ⚠️ |
| P2 | **Fix broken cross-references** — undefined `sec:conclusion`, `app:robustness`, `app:weighting_robustness`, `eq:labor_C_N`; multiply-defined `eq:profit_share` | ⚠️ |
| P3b | **Reconcile `Draft.tex:421` with the new δ_e**: the intro endorses Gabrovski-Silva's 6–10%/yr range while the calibration would use 3.2%. Rewrite to explain why the added channels (ω_δ, variety, ξ) permit a lower empirically-grounded δ — or report an estimated δ_e instead | ⚠️ **new, created by [D1](decisions.md)** |
| P3 | **Cite the Blanchard-Kahn notes** — `Notes/Baseline_Blanchard_Kahn.md` establishes BK under ε > 1 and that endogenous exit strengthens BK. Currently uncited; worth a footnote or a short appendix subsection | ❌ |
| P4 | **§5.2 rewrite** to match whichever calibration path D2 selects | blocked on D2 |

## Model tasks (secondary)

| # | Task | Status |
|---|---|---|
| M1 | ~~Comparison B extension at ψ ≈ 1.0~~ — folded into [D2](decisions.md): estimate `dest_elast_target` and report the implied ψ, rather than recalibrating | ✅ superseded |
| M2 | **Free-entry steady state** in `steady_state_checks.jl` — solver lands on the θ=1.80 branch; fix is a bracketed solve over `(log(0.1), log(0.6))`. See [D9](decisions.md) | ❌ low priority |
| M3 | ~~Verify the calibrated ψ~~ ✅ **resolved Sept 5, 2026**: 0.014 is the PATH A outcome, 0.033 the PATH B outcome. Only 0.033 applies to the current Comparison B | ✅ |
| M5 | **Pre-estimation diagnostic** — `run_solution_delta_target.jl` ✅ **written and run Sept 5, 2026**. Compares dest_ann ∈ {0.0320 BED, 0.0754 code, 0.0963 BGM} on σ(θ)/σ(labor_prod) vs 11.70 and δ→u persistence vs the LP peak at h=17–20. Serializes `irf_delta_target.jls`. See [D1](decisions.md) | ✅ **run Sept 5, 2026** — see findings.md §D1/M5. Does NOT settle D1: all three specifications miss amplification by ~10× and peak at h=1 vs the LP h=17–20 |
| M7 | 🔴 **Fix the observable mapping for labor productivity.** Model `labor_prod = Y/(ρL)` has SD 0.0387 vs 0.0128 in data and correlates 0.9975 with N^e — it tracks the ν_f·N^e entry term in Y, not technology. Every RSD in Block M has this in the denominator. Candidates: Y_c/L_c, or a differently deflated series. Then write an explicit observable-mapping table into §5.2 | ❌ **blocks Block M** |
| M8 | 🔴 **Chase the hump-shape gap**: model δ→u peaks at h=1 quarter, LP at h=17–20. If Θ_e cannot close this, Block B and Block M will fight. Establish before building the sampler | ❌ |
| M6 | **`hp_filter` was not the HP filter** — fixed in `time_series_fun.jl` Sept 5, 2026 (was a first-difference/Whittaker smoother: ~5× the correct cycle SD, corr 0.26 with true HP at λ=1600). Now pentadiagonal, matches `statsmodels.hpfilter` to 1e-12. **Any model-side second moment computed before this date is invalid**, incl. anything from `second_moments.jl`/`second_moments_GS.jl` | ✅ fixed |
| M9 | 🔴 **Reproducibility remediation (principles.md N15).** The Sept 5, 2026 rewrite of §5.3 introduced ~24 numbers (firm-stock and shock half-lives, f_e levels, IRF peaks/troughs, exit-flow decomposition, w_int gap) computed in ad-hoc Julia sessions, not emitted by any program. ✅ **Closed Sept 6, 2026.** `mechanism_stats.jl` emits 66 quantities to `mechanism_stats.txt` and `mechanism_stats.tex`, all reproducing the values in §5.3. **Decision (MS):** the prose keeps the literals rather than `\input`-ing the macros. The numbers are checked against `mechanism_stats.txt` during draft updates, and that file is the authority. Re-run the script after any mechanism runner and diff it before touching §5.3 | ✅ |
| M4 | **Consider fixing the `eval_SS` toolkit bug** (`for ip in npar` iterates once) so callers need not pass a pre-computed SS — matters once the sampler calls the solver thousands of times | ❌ |

---

## Completed

Comparisons A–D are complete; their findings are in [`findings.md`](findings.md) and their
files in [`pipeline.md`](pipeline.md). The full completed-work log through June 7, 2026 —
`run_solution_core.jl` audit, the λ=1,600 switch, the BED Deaths switch, Propositions 1–6
and their proofs, the §4.2/§4.3 and §5.1/§5.2 drafting, `app:filter_robustness` — is in the
git history (`git log --since=2026-04-01`) and reflected in the status tables of
[`draft_status.md`](draft_status.md). It is no longer duplicated here.
