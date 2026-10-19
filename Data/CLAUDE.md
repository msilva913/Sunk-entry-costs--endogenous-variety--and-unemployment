# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** October 18, 2026 · **Author:** Mario Silva
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
**Calibration:** τ̄ = 9.30%/qtr (3.10%/month). **δ̄_e is settled** ([D1](context/decisions.md)): `dest_ann = 0.0320` (BED Deaths, employment-weighted, 1993–2019), **δ_e/τ = 0.087**. **Implemented Oct 18 on `r8_timing_reposting`** (Q1 step 5): `TARGETS` holds D1 + PATH B (`dest_elast_target = 5.0`), b/w = 0.9, x_v = 0.5, `Λr_target = nothing`. D2 settled (PATH B).
**Filter:** HP λ=1,600 on log-levels, applied identically to data and simulations. Not λ=100,000 — that flips the sign of cor(δ,u) and cor(δ,v).
**Paper draft:** `Draft/Draft.tex` (renamed from `Draft_11May2026.tex` on June 5, 2026). Theory and empirics complete through §5.3 (theory audited and corrected Oct 11–12); §5.4 is an empty stub and §6 is unstarted. Build: three `pdflatex` passes, **never `bibtex`**. Section-by-section status and the proposition inventory: [`context/draft_status.md`](context/draft_status.md).
**Model code:** `Programs baseline/run_solution_core.jl`, restructured by R8 (Oct 14–18, branch `r8_timing_reposting`): draft timing, states [Ñ, Ẽ, Ṽ^u, Ṽ^s, z, δ, δl, s], partial reposting (Λ_r, Q_rep, H), 41 equations; X^c exact. Build PAR with `build_PAR(cal, SHK)`, shocks with `shock_col` and `SHOCKS_BASE`. Conventions and the Q1 program list: [`context/pipeline.md`](context/pipeline.md) Part 2.
**No estimation code exists yet.** `run_solution_core.jl:131` still reads `estimate = []`; `:133` reads `priors = (;)`.

---

## Priority tasks right now

> **Oct 14, 2026 — the critical path is now Q0–Q6** ([`context/pending_tasks.md`](context/pending_tasks.md),
> [D13](context/decisions.md)). Bartik IRFs leave the objective and become validation.
> Block B′ = aggregate SVAR responses by indirect inference. Shock processes come from a
> trivariate (z, δ, s) VAR. **E9, D3 and Ω_β are superseded**, so the D10/D3/E9 list below is
> history. The order is: Q0 notation and presentation (in parallel) → Q1 = R8 + D1 + PATH B →
> Q2 VAR → Q3 targets → Q4 identification → Q5 hypotheses and mechanism comparisons (the D12
> gate) → Q6 estimation.

**Resume:** [`context/session_handout_20261018.md`](context/session_handout_20261018.md)
(Oct 18: Q1 steps 5–7 done on **`r8_timing_reposting`**, with `costly_vacancy_reposting` merged
in; Test 6 found R19; BK ε-threshold; Remark 3 + App. K). Previous: `session_handout_20261014.md`
(Q0, D13, Q1 steps 1–4), `20261012` (theory rebuild), `20261009` (strategy, D12). Older handouts
are in `context/archive/`.

**Next, in order** (details in [`context/pending_tasks.md`](context/pending_tasks.md)):
1. **Close Q1:** `model_equations.md` Code status, `parameters.md` (Λ̄_r for α), `pipeline.md`
   conventions, D11 record; merge `r8_timing_reposting` → `costly_vacancy_reposting`.
2. **R19 (theory, MS):** pricing replacement offsets `prop:fe`'s user-cost channel in the
   baseline. Double-limit check first, then a remark after `prop:fe`, then the framing.
3. **Q2:** trivariate (z, δ, s) VAR. More urgent: with the part6b s process cor(u,v) reaches only
   +0.46 at Λ̄_r = 0.70 (data −0.80).

**Theory structure in the draft (Oct 18).** §3 limiting cases: `prop:independence`, clone
replacement (`prop:ags`, `rem:nesting`), `prop:double_limit`, `rem:bgm_nest`,
`rem:dilution_variety` (Remark 3; argument in App. K `app:determinacy`). `sec:mechanism`:
`lem:theta` + `prop:fe` (Prop. 4, `tab:fe_irf`), `lem:accounting`, `rem:finite_xi`,
`prop:exit_s` (Prop. 5). Steady state: `rem:curves`, `ass:matching` + `prop:equilibria`
(Prop. 6). Inventory: [`context/draft_status.md`](context/draft_status.md) "Propositions".
⚠ **Terminology (MS):** never write bare "free entry" for vacancies ("vacancy free entry" or
"infinitely elastic vacancy creation"); bare "free entry" means product lines.

**Standing warnings.**
- ⚠ **Build with `pdflatex` only** (three passes, `-synctex=1`, from `Data/Draft/`). Running
  `bibtex` destroys the tracked `Draft.bbl`. **Owed (MS):** add `silva_2026` (SSRN 6956640) and
  `den_haan_2000` to `references.bib`.
- ⚠ Close TeXstudio before scripted edits to `Draft.tex`. `Draft.tex` contains NUL bytes; use
  `grep -a`. Edit `Draft.tex` in the current checkout (both branches now carry the same draft).
- ⚠ Shell gotchas: the Bash tool collapses backslashes in heredocs and inline `-e` strings; write
  Julia/LaTeX edit scripts with the Write tool and `raw""` strings. `@sprintf` needs one literal
  format string. Never assign numbers to the symbolic model names (ρ_z, s, b, …), including inside
  REPL loops (Oct 18: `s = r.irf_s` broke `gen_model_equations` in VS Code).
- ⚠ No GitHub credentials in the agent shell: MS pushes from their own Git client.
- ⚠ The SMCD paper is cited as a working paper only. Do not mention where it is submitted.
- ⚠ **`Notes/role_of_ingredients.md` §2 is wrong**; rewrite before designing R14.
- ⚠ §5.3 prose is committed to runs that do not exist yet (novar-vs-AGS, α → 0, re-specified
  Comparison B); they belong to Q5.
