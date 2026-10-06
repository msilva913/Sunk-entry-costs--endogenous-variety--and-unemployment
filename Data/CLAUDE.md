# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)
**Last updated:** October 5, 2026 · **Author:** Mario Silva
**Current branch:** `costly_vacancy_reposting` (branched from `instruments_LP_coefficient`; the most up-to-date version of the project)
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
| [`context/decisions.md`](context/decisions.md) | **Decision register**: D1–D11 (D1 settled, D2/D3/D11 open, D10 partial), S1–S9 settled, each with a recommendation and a reason |
| [`context/pending_tasks.md`](context/pending_tasks.md) | Actionable work, ordered; the critical path to a complete draft |
| [`context/estimation_design.md`](context/estimation_design.md) | How conditional (IRF) and unconditional moments combine, determinacy, what is missing, build order |
| [`context/findings.md`](context/findings.md) | All empirical and model-mechanism results with numbers |
| [`context/draft_status.md`](context/draft_status.md) | Section-by-section draft state; proposition inventory and proof structures |
| [`context/data_and_files.md`](context/data_and_files.md) | Data sources, FRED series IDs, file index, empirical targets, LP spec |
| [`context/parameters.md`](context/parameters.md) | All structural parameters: Θ_ext / Θ_d / Θ_e, economic meaning, identifying moments |
| [`context/model_equations.md`](context/model_equations.md) | All 33 model equations with interpretation (Block 1 = current; Block 2 = reposting-margin modifications, [D11](context/decisions.md)) |
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
**Calibration:** τ̄ = 9.30%/qtr (3.10%/month). **δ̄_e is settled** ([D1](context/decisions.md), September 6, 2026): `dest_ann = 0.0320` (BED Deaths, employment-weighted, 1993–2019: 0.812%/qtr, 3.21%/yr), giving **δ_e/τ = 0.087**. ⏳ **Not yet implemented** — `Programs baseline/steady_state.jl:613` still holds the superseded `0.0754` (δ_e/τ = 0.210), which applied a gross-job-loss share to the separation rate. The previously recorded 11.6% is COVID-contaminated and has no standing.
**Filter:** HP λ=1,600 on log-levels, applied identically to data and simulations. Not λ=100,000 — that flips the sign of cor(δ,u) and cor(δ,v).
**Paper draft:** `Draft/Draft.tex` (renamed from `Draft_11May2026.tex` on June 5, 2026). Theory and empirics complete through §5.3; §5.4 is an empty stub and §6 is unstarted. Section-by-section status and the proposition inventory: [`context/draft_status.md`](context/draft_status.md).
**Model code:** 33-equation system in `Programs baseline/run_solution_core.jl`, audited May 20, 2026. Conventions, equation order, and `steady_state.jl` architecture: [`context/pipeline.md`](context/pipeline.md) Part 2.
**No estimation code exists yet.** `run_solution_core.jl:131` still reads `estimate = []`; `:133` reads `priors = (;)`.

---

## Priority tasks right now

The paper is theory-complete and empirics-complete. **Section 5 is the whole remaining
critical path.** D1 is settled (but not implemented). Three open decisions gate the rest —
settle these before writing any estimation code, because each one changes what the code has
to do:

1. **[D10] What is the Block B IRF target?** 🟡 The single most urgent item. E8 showed the
   baseline LP's peak horizon is an artifact of near-unit-root instrument persistence
   (ρ = 0.91 at lag 1), and the augmented LP's peak is itself unstable across lag orders.
   The next step is **E9**: run the identical LP on a model-simulated 50-state panel. If the
   model LP also builds monotonically, baseline-to-baseline matching preserves the full
   42-dimensional IRF target; if it peaks at h=1–2, fall back to short-horizon targeting.
2. **[D2] PATH A or PATH B calibration?** The draft documents PATH A (`Xc_Y` target);
   every §5.3 figure was produced with PATH B (`dest_elast_target`). This changes the
   contents of Θ_e, so it must precede the priors.
3. **[D3] How do β̂ and β(θ) become comparable?** β̂ is a cross-state response per 1 SD of the
   Bartik instrument with time FEs absorbed; β(θ) is an aggregate response per 1 SD of the
   structural shock. The draft never says how. Determines whether Block B identifies σ_δ or
   only the shape of the response. **Note the overlap with D10:** E9 *is* the
   Guren-McKay-Nakamura-Steinsson "simulate the regression inside the model" prescription,
   which is D3's most defensible option. Building E9 well upgrades D3.

**Resume:** [`context/session_handout_20261001.md`](context/session_handout_20261001.md) (R13 is **done**; read its Oct 5 update).

**October 5, 2026 session — theoretical positioning and a derivation bug.** No model
runs; draft prose, one new program, one new note.

- **Intro positioning.** Endogenous exit is no longer sold as a recession amplifier
  (contradicted by §5.3, Flag 1). A new paragraph gives the measurement argument for
  the business formation block: δ is identified from establishment deaths, so lines must
  hold many positions, which needs profits, which under DS-CES *is* the taste for
  variety (ζ = μ − 1). Fact 1 now discusses the entry margin it always claimed to.
- **Post-Prop-3 passage** rewritten: four channels over AGS, not three (**profit
  dilution** added), plus the point that AGS does not *price* replacement.
- **`app:loglin` ¶9 was circular** — `eq:K_decomp` and `eq:Q_ll_delta` are one
  asset-pricing identity, used twice to sign both K̂ and Q̂. The entry cushion is
  **congestion relief**, not duration shortening. New ¶9a/9b/9c.
- **New result.** The threshold entry elasticity for vacancies to rise after a δ shock
  is T = κ_v/(δ̄_e/τ̄) − 1, so **CK's δ_e = τ minimizes it over the admissible range**
  (T = 0.203 against 12.78 at D1). Their need for finitely elastic vacancy creation
  repairs a consequence of having one separation margin. See
  [`Notes/predetermined_vacancies_and_the_CK_special_case.md`](Notes/predetermined_vacancies_and_the_CK_special_case.md).
  ⚠ This is **GS's δ/s split plus our measurement**, not the business formation block.
  ⛔ **Superseded Oct 6, 2026.** T = 12.78/0.203 is the bar per unit *level* of δ_e,
  not per 1 %. Per 1 % it is 0.035 against 0.0065, tiny in both cases. The note is
  rewritten around two other channels: the stock elasticity ξ·ē/v̄ (finite ξ), and the
  user cost of a durable vacancy (free entry, where the nested model's u–v slope is
  −0.23 at D1 and +0.17 at δ_e = τ). CK's free-entry case is flow-cost DMP. Program:
  `Programs baseline/beveridge_free_entry.jl`.
- **New tasks.** R8 clause (v) (re-specify Comparison B), **R12b** (ê^δ > 0 needs a
  program cite), **R12c** (the Prop. 5 Part 3 proof uses the same circular step — not
  fixed, needs Part 3 worked through).
- ⚠ **`Notes/role_of_ingredients.md` §2 is wrong** and flagged inline: variety is *not*
  what makes the δ sign survive free entry. Flag 3 and the ξ/variety substitutes
  conclusion rest on it. **Rewrite before designing R14.**
- ⚠ **Build with `pdflatex` only.** Running `bibtex` destroys `Draft.bbl`, which is
  tracked, because `references.bib` is not in `Data/Draft/`.
- ⚠ §5.3 is now committed in prose to three runs that do not exist: novar-vs-AGS,
  α → 0, and a re-specified Comparison B. Check whether any script implements
  `prop:ags` as stated before R8 starts.

**Code timing (Oct 1, 2026).** The draft is the reference model ([S10](context/decisions.md)).
The code dates exit one period early in its laws of motion and uses post-exit
stocks as states. Fix in task R8 per [`context/model_equations.md`](context/model_equations.md)
"Code status", in the same regeneration pass as the D1 cascade and PATH B.

**Also owed: the D1 cascade.** D1 is settled but unimplemented. Set `dest_ann = 0.0320`,
re-run Comparisons A–D, regenerate the eight `mechanism_*.pdf`, re-run `mechanism_stats.jl`,
and update §5.3, §5.2, `tab:calib_targets` row 1, the **three remaining** `δ̄_e/τ̄ ≈ 0.21`
claims in `Draft.tex` (§5.3, `app:loglin`, `app:comparison_D` — find them with a search for
`0.21`, since line numbers drift; the fourth, in `sec:mechanism`, was resolved Oct 5), and
**P3b**.

Then build, in order: Ω_β (`part5_wcrb.py`, a full matrix — not per-horizon SEs; **its
dimension is set by D10/E9**, 42×42 only if the full IRF path survives) → Ω_m (block
bootstrap) → refreshed m(θ) at λ=1,600 over all seven series → β(θ) extractor → objective and
priors → sampler → §5.4, `app:weighting_robustness`, §6.

Lower priority: δ–LD paragraph in §4.3; four undefined `\ref`s and one duplicate label in
`Draft.log`; re-run the stale SLOOS and GFC LPs against the current 17.4 pp instrument; cite
the Blanchard-Kahn notes.

See [`context/decisions.md`](context/decisions.md) for the full decision register with
recommendations, and [`context/pending_tasks.md`](context/pending_tasks.md) for the task list.
