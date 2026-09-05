# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** September 5, 2026 · **Author:** Mario Silva
**Project paused June 7, 2026 – September 5, 2026.** Status re-verified against the repo on
resumption; see [`context/pending_tasks.md`](context/pending_tasks.md) §Blockers.
**Publication objective** applied theory, quantitative macro paper publishable in top field journal like Journal of Monetary Economics/Americal Economic Journal:Macroeconomics and potentially in Journal of Political Economy.

**Code style:** clean, succinct, interpretable by empirical macroeconomist. 

**Draft writing style: **clean, succinct, applied macroeconomist/theorist. Avoid writing Appendix before appendix reference because that causes "Appendix to be printed twice.

**Priority**: core findings must mark clear departure with respect to Coles and Kelishomi (2018)
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
| [`context/pending_tasks.md`](context/pending_tasks.md) | Prioritized code + paper writing tasks; **blockers at top** |
| [`context/estimation_design.md`](context/estimation_design.md) | How conditional (IRF) and unconditional moments combine; what is missing to run the estimation |

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
**Calibration:** δ̄ = 1.079%/qtr (BED Deaths) = 4.25%/yr, τ̄ = 9.34%/qtr, δ/τ = 11.6% — **CONTRADICTED by the code and by §5.2**, which use `dest_ann = 0.0754` (Jaimovich-Siu 21% × τ), i.e. δ_e/τ = 0.210. Resolve before estimating (blocker B1).
**Paper draft:** `Draft/Draft_11May2026.tex` — Intro + Sections 2–4.2 complete; prop:bgm_nest + proof complete; 4.3 partially drafted. D5 extended: K decomposition eq:K_decomp–eq:e_ll_delta. New appendix app:additional_comparisons with Comparison D (ξ role, CK connection). **prop:ds_asymmetry complete** (3 parts: s/exog, s/endog, δ): Part 2 condition eq:gN_cond ($f_e\bar N < \bar z(1-\bar u)(\varepsilon-2)/(\varepsilon-1)$, requires ε>2 + small entry costs); Remark defines g(N,u), explains dilution-vs-variety tradeoff. Proof: Channel 2 formal (sub-steps a/b via N_{t+1} predetermination); Part 3 uses K-decomp coefficient $\bar\delta_e/(r+\bar\delta_e)\to 0$ as δ/τ→0.
**Model code:** `Programs baseline/run_solution_core.jl` audited May 20: r=discount rate, ρ=relative price (consistent with steady_state.jl); κ=matching cost paid per match (κ·q·v in RC); u LOM uses total v_t=v_pret+e_t; BFE at f[4] after JCC. Comparison D: `run_solution_entry_elasticity.jl` / `plot_xi_inv_comparison.jl`; ξ_inv=1.0 vs 0.1; serializes `irf_xi_inv.jls`.

---

## Priority tasks right now

The paper is theory-complete and empirics-complete; **Section 5 is the whole remaining
critical path**, and it is blocked on three things, in this order:

1. **Fix the δ_e calibration contradiction** (BED Deaths 4.25%/yr vs. `dest_ann = 0.0754`).
   Everything downstream — mechanism IRFs, Prop. 5 Part 3's δ̄_e/τ̄ claim, Comparison D — is
   conditioned on this number.
2. **`part5_wcrb.py`** — wild cluster bootstrap producing the **full 42×42 Ω_β** across
   horizons and outcomes, not just per-horizon SEs. `eq:ql_irf` cannot be evaluated without it.
3. **Write the estimation code.** None exists: `run_solution_core.jl:112` still has
   `estimate = []`. Need Ω_m (block bootstrap), an m(θ) simulator at λ=1,600 covering
   {u,v,s,f,δ,N^e,z}, a β(θ) extractor matching the LP object, and a sampler.
   See [`context/estimation_design.md`](context/estimation_design.md), including the
   unresolved scale/units mapping between β̂ and β(θ).

Lower priority / opportunistic: δ–LD paragraph in §4.3; `app:weighting_robustness`; four
undefined `\ref`s and one duplicate label in `Draft.log`; re-run stale SLOOS and GFC LPs
against the current 17.4 pp instrument; §6 Conclusion.

See [`context/pending_tasks.md`](context/pending_tasks.md) for full list.
