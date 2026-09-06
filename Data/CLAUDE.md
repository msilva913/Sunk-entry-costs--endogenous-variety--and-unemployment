# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** September 5, 2026 · **Author:** Mario Silva
**Project paused June 7, 2026 – September 5, 2026.** Status re-verified against the repo on
resumption; see [`context/pending_tasks.md`](context/pending_tasks.md) §Blockers.
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
| [`context/decisions.md`](context/decisions.md) | **Decision register**: D1–D9 open, S1–S9 settled, each with a recommendation and a reason |
| [`context/pending_tasks.md`](context/pending_tasks.md) | Actionable work, ordered; the critical path to a complete draft |
| [`context/estimation_design.md`](context/estimation_design.md) | How conditional (IRF) and unconditional moments combine, determinacy, what is missing, build order |
| [`context/findings.md`](context/findings.md) | All empirical and model-mechanism results with numbers |
| [`context/draft_status.md`](context/draft_status.md) | Section-by-section draft state; proposition inventory and proof structures |
| [`context/data_and_files.md`](context/data_and_files.md) | Data sources, FRED series IDs, file index, empirical targets, LP spec |
| [`context/pipeline.md`](context/pipeline.md) | Both pipelines (Python empirical + Julia model): scripts, run order, code conventions |
| [`context/principles.md`](context/principles.md) | 21 standing rules every agent must follow |

---

## Other key files:
See Key Papers for most important/closely related papers.

## Quick-reference: most important facts

**Baseline instrument:** δ Bartik (BED Deaths), v2 residualized, σ̂ = 17.4 pp, sample 2001Q1–2019Q4
**Headline IRFs:** δ→u peak +1.71 pp (h=17); δ→v trough −0.58 pp (h=18); sig from h=0/h=2 respectively
**Headline validation:** severity placebo (u^nat × B^δ) insignificant at ALL h=0..20 for vacancy outcome
**Instrument correlation:** r(δ,LD) = 0.334 after v2 + Deaths switch; industry-level financial conditions are leading explanation
**QU placebo fails:** corr(β_unemp, β_vac) = −0.971 across horizons — demand contamination signature
**Vacancy convention:** outcome is vacancy RATE = V/LF × 100 in pp, levels difference. Never log. Symmetric to unemployment rate.
**Calibration:** τ̄ = 9.30%/qtr (3.10%/month). **δ_e/τ ≈ 0.08** is the recommended primary target (BED Deaths, 2001–2019: 0.718%/qtr, 2.84%/yr). The code's `dest_ann = 0.0754` (δ_e/τ = 0.210) applies a gross-job-loss share to the separation rate and should be retired; the previously recorded 11.6% is COVID-contaminated. See [D1](context/decisions.md).
**Filter:** HP λ=1,600 on log-levels, applied identically to data and simulations. Not λ=100,000 — that flips the sign of cor(δ,u) and cor(δ,v).
**Paper draft:** `Draft/Draft.tex` (renamed from `Draft_11May2026.tex` on June 5, 2026). Theory and empirics complete through §5.3; §5.4 is an empty stub and §6 is unstarted. Section-by-section status and the proposition inventory: [`context/draft_status.md`](context/draft_status.md).
**Model code:** 33-equation system in `Programs baseline/run_solution_core.jl`, audited May 20, 2026. Conventions, equation order, and `steady_state.jl` architecture: [`context/pipeline.md`](context/pipeline.md) Part 2.
**No estimation code exists yet.** `run_solution_core.jl:112` still reads `estimate = []`.

---

## Priority tasks right now

The paper is theory-complete and empirics-complete. **Section 5 is the whole remaining
critical path.** Three decisions gate it — settle these before writing any estimation code,
because each one changes what the code has to do:

1. **[D1] Which δ̄_e?** Recommended: target **δ_e/τ ≈ 0.08** directly (BED Deaths, 2001–2019:
   0.718%/qtr, 2.84%/yr). The code's `dest_ann = 0.0754` (δ_e/τ = 0.210) comes from applying
   a *gross-job-loss* share to the *separation* rate, and the draft's 11.6% is COVID-
   contaminated. Every mechanism figure and Prop. 5 Part 3's δ̄_e/τ̄ claim rests on this.
2. **[D2] PATH A or PATH B calibration?** The draft documents PATH A (`Xc_Y` target);
   every §5.3 figure was produced with PATH B (`dest_elast_target`). This changes the
   contents of Θ_e, so it must precede the priors.
3. **[D3] How do β̂ and β(θ) become comparable?** β̂ is a cross-state response per 1 SD of the
   Bartik instrument with time FEs absorbed; β(θ) is an aggregate response per 1 SD of the
   structural shock. The draft never says how. Determines whether Block B identifies σ_δ or
   only the shape of the response.

Then build, in order: Ω_β (`part5_wcrb.py`, full 42×42 — not per-horizon SEs) → Ω_m (block
bootstrap) → refreshed m(θ) at λ=1,600 over all seven series → β(θ) extractor → objective and
priors → sampler → §5.4, `app:weighting_robustness`, §6.

Lower priority: δ–LD paragraph in §4.3; four undefined `\ref`s and one duplicate label in
`Draft.log`; re-run the stale SLOOS and GFC LPs against the current 17.4 pp instrument; cite
the Blanchard-Kahn notes.

See [`context/decisions.md`](context/decisions.md) for the full decision register with
recommendations, and [`context/pending_tasks.md`](context/pending_tasks.md) for the task list.
