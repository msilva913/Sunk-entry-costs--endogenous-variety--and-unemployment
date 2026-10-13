# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** October 12, 2026 · **Author:** Mario Silva
**Current branch:** `costly_vacancy_reposting` (branched from `instruments_LP_coefficient`; the most up-to-date version of the project)
**Project paused June 7, 2026 – September 5, 2026.** Status re-verified against the repo on
resumption; see [`context/pending_tasks.md`](context/pending_tasks.md) §Critical path.
**Publication objective:** applied theory, quantitative macro paper publishable in a top field journal such as the *Journal of Monetary Economics* or *AEJ: Macroeconomics*, and potentially the *Journal of Political Economy*.

**Priority:** core findings must mark a clear departure from Coles and Kelishomi (2018).

> ### ⛔ Read before doing anything
> [`context/principles.md`](context/principles.md) closes with **15 non-negotiable rules
> (N1–N15)** covering conduct, writing voice, code style, LaTeX conventions, citation
> integrity, and reproducibility. They override convenience and any instruction that
> conflicts with them. Four that bite constantly:
>
> - **N1/N2** — treat every request as a hypothesis; audit model-wide before structural
>   edits, and push back with a diagnosis rather than complying silently.
> - **N4/N5** — applied-macroeconomist voice, **American English**, short sentences,
>   sparing dashes.
> - **N11–N14** — never fabricate a citation; DOIs required; say
>   "Citation unverified; source missing" rather than guess.
> - **N15** — every number in the draft must come out of a program in the replication
>   workflow. Never paste a figure computed in an ad-hoc session.
---

## Project in one paragraph

Empirical component of a DSGE paper extending the Gabrovski-Silva framework (Bilbiie, Ghironi & Melitz 2012 + Coles & Kelishomi 2018). Three aggregate shocks: **δ** (permanent firm/product-line destruction), **s** (match separation), **z** (technology). Core prediction: δ destroys both product lines and vacancies simultaneously → persistent Beveridge curve dynamics; s shocks self-correct via costless reposting. Bartik shift-share instruments built from BLS BED Deaths (δ) and JOLTS layoffs+discharges (LD) identify shock-specific IRFs via panel local projections across 50 US states.

---

## Working directory

```
Data/Bartek analysis/
  data/cache/        — parquet cache files
  data/instruments/  — Bartik instrument CSVs and parquets
  data/results/      — LP output CSVs and all plots
```

---

## Context files (read these for details)

Start at [`context/README.md`](context/README.md) — it indexes the folder and says which
file owns what.

| File | Contents |
|------|----------|
| [`context/README.md`](context/README.md) | **Index and reading order.** Start here |
| [`context/decisions.md`](context/decisions.md) | **Decision register**: D1–D12 (D1 settled, D2/D3/D11/D12 open, D10 partial), S1–S12 settled, each with a recommendation and a reason |
| [`context/pending_tasks.md`](context/pending_tasks.md) | Actionable work, ordered; the critical path to a complete draft |
| [`context/estimation_design.md`](context/estimation_design.md) | How conditional (IRF) and unconditional moments combine, determinacy, what is missing, build order |
| [`context/findings.md`](context/findings.md) | All empirical and model-mechanism results with numbers |
| [`context/draft_status.md`](context/draft_status.md) | Section-by-section draft state; proposition inventory and proof structures |
| [`context/data_and_files.md`](context/data_and_files.md) | Data sources, FRED series IDs, file index, empirical targets, LP spec |
| [`context/parameters.md`](context/parameters.md) | All structural parameters: Θ_ext / Θ_d / Θ_e, economic meaning, identifying moments |
| [`context/model_equations.md`](context/model_equations.md) | All 33 model equations with interpretation (Block 1 = current; Block 2 = reposting-margin modifications, [D11](context/decisions.md)) |
| [`context/nesting_map.md`](context/nesting_map.md) | How the baseline nests or approximates DMP, CK, GS and BGM: restrictions, equivalence type, labels, proofs, caveats. **Read before editing the limiting cases or citing a nesting result** |
| [`context/pipeline.md`](context/pipeline.md) | Both pipelines (Python empirical + Julia model): scripts, run order, code conventions |
| [`context/principles.md`](context/principles.md) | 21 standing rules every agent must follow |

---

## Other key files:
See Key Papers for most important/closely related papers.

## Quick-reference: most important facts

**Baseline instrument:** δ Bartik (BED Deaths), v2 residualized, σ̂ = 17.4 pp, sample 2001Q1–2019Q4
**Headline IRFs:** δ→u peak +1.71 pp (h=17); δ→v trough −0.58 pp (h=18); sig from h=0/h=2 respectively.
⚠️ The **peak horizon is not robust** to instrument persistence (E8, Sept 22–23): with lagged
instruments the δ→u peak drifts to h=10–13 and 0.95–1.86 pp. Stable features are δ→u at h=0–2
and the δ→v trough at h=4–6. See [D10](context/decisions.md)
**Headline validation:** severity placebo (u^nat × B^δ) insignificant at ALL h=0..20 for vacancy outcome
**Instrument correlation:** r(δ,LD) = 0.334 after v2 + Deaths switch; industry-level financial conditions are leading explanation
**QU placebo fails:** corr(β_unemp, β_vac) = −0.971 across horizons — demand contamination signature
**Vacancy convention:** outcome is vacancy RATE = V/LF × 100 in pp, levels difference. Never log. Symmetric to unemployment rate.
**Calibration:** τ̄ = 9.30%/qtr (3.10%/month). **δ̄_e is settled** ([D1](context/decisions.md), September 6, 2026): `dest_ann = 0.0320` (BED Deaths, employment-weighted, 1993–2019: 0.812%/qtr, 3.21%/yr), giving **δ_e/τ = 0.087**. ⏳ **Not yet implemented** — `Programs baseline/steady_state.jl:623` still holds the superseded `0.0754` (δ_e/τ = 0.210), which applied a gross-job-loss share to the separation rate. The previously recorded 11.6% is COVID-contaminated and has no standing.
**Filter:** HP λ=1,600 on log-levels, applied identically to data and simulations. Not λ=100,000 — that flips the sign of cor(δ,u) and cor(δ,v).
**Paper draft:** `Draft/Draft.tex` (renamed from `Draft_11May2026.tex` on June 5, 2026). Theory and empirics complete through §5.3 (theory audited and corrected Oct 11–12); §5.4 is an empty stub and §6 is unstarted. Build: three `pdflatex` passes, **never `bibtex`**. Section-by-section status and the proposition inventory: [`context/draft_status.md`](context/draft_status.md).
**Model code:** 33-equation system in `Programs baseline/run_solution_core.jl`, audited May 20, 2026; X^c exact since Oct 11 (`payer_share` in `steady_state.jl`). Conventions, equation order, and `steady_state.jl` architecture: [`context/pipeline.md`](context/pipeline.md) Part 2.
**No estimation code exists yet.** `run_solution_core.jl:131` still reads `estimate = []`; `:133` reads `priors = (;)`.

---

## Priority tasks right now

**Resume:** [`context/session_handout_20261012.md`](context/session_handout_20261012.md)
(Oct 11–12: theory rebuild and X^c fix). Previous: `context/session_handout_20261009.md`
(strategy, D12). Older handouts are in `context/archive/`.

**Oct 13, 2026 — review and cleanup (P8, done).** Theory-section presentation fixes. New
`ass:vacancy_fe`; `prop:bgm_nest` is now `rem:bgm_nest`; new `rem:dmp`. Under vacancy free
entry with e_t > 0, the model is isomorphic to DMP in which every position bears the flow user
cost K_t, so v_pre matters only if e_t ≥ 0 binds. Propositions renumbered (fe = 4, exit_s = 5,
equilibria = 6). ⚠ **Terminology rule (MS):** never write bare "free entry" for
vacancies. Say "vacancy free entry" or "infinitely elastic vacancy creation"; bare "free entry"
means product lines (`eq:free_entry`). The main remaining theory gap is **R3** (no result for
partial reposting). Details: `context/draft_status.md` "Propositions".

The theory sections are audited and corrected (Oct 11–12). **Section 5 is the remaining
critical path**, and the code must first be brought in line with the draft. Decisions that gate
the rest (all MS's; [`context/decisions.md`](context/decisions.md)):

1. **[D12] Venue and structure.** Recommendation: no JET pivot; a theory-forward paper for
   RED / JME / AEJ: Macro that keeps the SVAR-calibrated mechanism and cuts Bartik to the
   robust short-horizon facts. If adopted, E9, Ω_β and D3 leave the critical path and D10 takes
   option (a). **The publication objective above is unchanged.**
2. **[D2] PATH A or PATH B.** The draft documents PATH A; every §5.3 figure uses PATH B. D1
   converges only under PATH B. Recommendation PATH B.
3. **[D10]/[D3]** (only if IRF matching stays): E9, the model-side LP test, decides both.

**Next model task: one regeneration pass R8 + D1 + PATH B** (code restructuring to the draft's
timing plus the reposting margin; `dest_ann = 0.0320`). Spec in
[`context/model_equations.md`](context/model_equations.md) "Code status" and Block 2. Then
re-run Comparisons A–D, `mechanism_stats.jl`, the `run_prop5_*` checks, and update §5.3, §5.2,
`tab:calib_targets` row 1 and the three remaining `δ̄_e/τ̄ ≈ 0.21` claims (`grep -a 0.21`).
This is also the D11 go/no-go: under costless reposting s shocks give positive u–v comovement at
every persistence (`prop:fe` part 3, `lem:accounting`), so the reposting margin is *necessary*
to fix model cor(u,v) = +0.995.

**Theory structure in the draft (Oct 12).** §3 limiting cases: `prop:independence`, clone
replacement (`prop:ags`, `rem:nesting`), `prop:double_limit` (corrected), `prop:bgm_nest`.
`sec:mechanism`: `lem:theta` + `prop:fe` (free entry, δ–s theorem, `tab:fe_irf`),
`lem:accounting`, `rem:finite_xi`, `prop:exit_s`. Steady state: `prop:curves` (a remark),
`ass:matching` (with `eq:Phi`) + `prop:equilibria`. Inventory and audit record:
[`context/draft_status.md`](context/draft_status.md) "Propositions".

**Standing warnings.**
- ⚠ **Build with `pdflatex` only** (three passes, `-synctex=1`, from `Data/Draft/`). Running
  `bibtex` destroys the tracked `Draft.bbl`, because `references.bib` is not in `Data/Draft/`.
  **Owed (MS):** add `silva_2026` (SSRN 6956640) and `den_haan_2000` to `references.bib`; both
  exist only in `Draft.bbl`.
- ⚠ Close TeXstudio before scripted edits to `Draft.tex` (EPERM on rename otherwise).
  `Draft.tex` contains NUL bytes (pre-existing); use `grep -a`.
- ⚠ The SMCD paper is cited as a working paper only. Do not mention where it is submitted.
- ⚠ **`Notes/role_of_ingredients.md` §2 is wrong** (variety is not what makes the δ sign survive
  free entry). Rewrite before designing R14.
- ⚠ §5.3 prose is committed to three runs that do not exist: novar-vs-AGS, α → 0, and a
  re-specified Comparison B (R8 clause v).
- **Code timing (Oct 1).** The draft is the reference model ([S10](context/decisions.md)). The
  code dates exit one period early and uses post-exit stocks as states; fixed in R8.

After the regeneration, build in order: Ω_β (`part5_wcrb.py`, full matrix; dimension set by
D10/E9) → Ω_m (block bootstrap) → refreshed m(θ) at λ = 1,600 over all seven series → β(θ)
extractor → objective and priors → sampler → §5.4, `app:weighting_robustness`,
`app:robustness`, §6. Full ordered list: [`context/pending_tasks.md`](context/pending_tasks.md).
