# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** May 14, 2026 · **Author:** Mario Silva
**Code style:** clean, succinct, interpretable by empirical macroeconomist

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

| File | Contents |
|------|----------|
| [`context/pipeline.md`](context/pipeline.md) | Full script table, run order, technical conventions |
| [`context/findings.md`](context/findings.md) | All empirical results, IRF numbers, interpretation, paper draft status |
| [`context/data_and_files.md`](context/data_and_files.md) | Data sources, FRED series IDs, key file index, calibration targets, LP spec reference |
| [`context/principles.md`](context/principles.md) | 16 principles every agent must follow |
| [`context/pending_tasks.md`](context/pending_tasks.md) | Prioritized code + paper writing tasks |

---

## Quick-reference: most important facts

**Baseline instrument:** δ Bartik, v2 residualized, σ̂ = 14.0 pp, sample 2001Q1–2019Q4
**Headline IRFs:** δ→u peak +1.50 pp (h=20); δ→v trough −0.44 pp (h=14); both highly significant
**Headline validation:** severity placebo (u^nat × B^δ) insignificant at ALL h=0..20 for vacancy outcome
**Instrument correlation:** r(δ,LD) = 0.423 survives v2; industry-level financial conditions are leading explanation
**QU placebo fails:** corr(β_unemp, β_vac) = −0.971 across horizons — demand contamination signature
**Vacancy convention:** outcome is vacancy RATE = V/LF × 100 in pp, levels difference. Never log. Symmetric to unemployment rate.
**Calibration:** δ̄ = 0.940%/qtr, τ̄ = 9.34%/qtr, δ/τ = 10.1%
**Paper draft:** `Draft/Draft_11May2026.tex` — Intro + Sections 2–4.2 complete; prop:bgm_nest + proof complete; 4.3 partially drafted

---

## Priority tasks right now

1. Draft Section 4.3 Results (δ→u, δ→v IRF narrative; severity placebo paragraph; joint LP robustness)
2. Draft Section 3.10 (δ-vs-s asymmetry mechanism proposition)
3. Wild cluster bootstrap (`part5_wcrb.py`) — n=50 clusters at lower reliability bound
4. GFC diagnostic for LD unemployment (`part7d_ld_gfc.py`)

See [`context/pending_tasks.md`](context/pending_tasks.md) for full list.
