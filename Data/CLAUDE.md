# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)

**Last updated:** March 27, 2026
**Author:** Mario Silva
**Purpose:** Persistent project context for fresh Cowork sessions. Paste this into any new session to restore full project state.

---

## 1. Project Overview

Empirical component of a DSGE paper extending the **Gabrovski-Silva** framework (building on Bilbiie, Ghironi & Melitz 2012). The model adds full business formation and exit to a search-and-matching labor market with finitely elastic vacancy creation (Coles & Kelishomi 2018). Three aggregate shocks: **δ** (firm/product line destruction), **s** (match separation), **z** (technology).

**Model's core prediction (the Beveridge asymmetry):**
- δ shock → destroys both product lines and vacancies simultaneously → persistent Beveridge curve dynamics
- s shock → destroys only matches; surviving firms repost vacancies at zero sunk cost → vacancy stock partially self-corrects

**Role of LP exercise:** Provide external empirical discipline on the δ vs. s transmission asymmetry using Bartik (shift-share) instruments. Panel local projections (LPs) use state-level variation to identify shock-specific IRFs.

---

## 2. Working Directory

```
Data/Bartek analysis/
```

All scripts run from this directory. Key subdirectories:
```
data/cache/          — parquet cache files (BLS raw data, pre-processed)
data/instruments/    — Bartik instrument CSVs and parquets, shares, outcomes
data/results/        — LP output CSVs and all IRF plots
```

---

## 3. Pipeline — Full Status

### Scripts (run in order):

| Script | Function | Status | Key Output |
|--------|----------|--------|-----------|
| `part1_shares.py` | Base-year (2006) state×supersector employment shares ω_{s,j} from QCEW | ✅ Complete | `data/instruments/shares_base2006.parquet` |
| `part2_shock_rates.py` | BED establishment closing rates → LOO national shock rates g^δ | ✅ Complete | `data/instruments/shock_rates_1992Q3_2024Q2.parquet` |
| `part2_shock_rates_s.py` | JOLTS separations (total, LD, quits) → LOO national shock rates | ✅ Complete | `shock_rates_s_*.parquet`, `shock_rates_ld_*.parquet`, `shock_rates_qu_*.parquet` |
| `part2b_shock_comovement.py` | Regress g^δ and g^s on Δlog(productivity) → residuals ν^δ, ν^s | ✅ Complete | `shock_rates_delta_resid.parquet`, `shock_rates_ld_resid.parquet`, `shock_rates_qu_resid.parquet` |
| `part3_instrument.py` | Bartik aggregation for δ | ✅ Complete | `data/instruments/delta_instrument_base2006.csv` |
| `part3_instrument_s.py` | Bartik aggregation for total separations (s) | ✅ Complete | `data/instruments/s_instrument_base2006.csv` |
| `part3_resid_instruments.py` | Bartik aggregation for residualized ν^δ, ν^{LD}, ν^{quits} | ✅ Complete | `delta_instrument_resid_base2006.csv`, `ld_instrument_resid_base2006.csv`, `qu_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | Fetches LAUS unemployment + JOLTS state vacancies → quarterly panel | ✅ Complete | `data/instruments/laus_quarterly.parquet` (includes `vacancies` column) |
| `part5_lp.py` | All panel LP regressions: unemployment + vacancy outcomes, δ/s/LD/quits | ✅ Complete | IRF CSVs + plots in `data/results/` |
| `part6_shock_persistence.py` | AR(1) estimation of shock persistence | ✅ Complete | `data/results/shock_persistence.csv` |
| `part7_nfci.py` | Chicago Fed NFCI → quarterly, headline + risk subindex | ✅ Complete | `data/cache/nfci_quarterly.parquet` |

### Notable additions in March 2026 session:
- **Priority 1 completed:** `part4_outcomes.py` now fetches JOLTS state-level job openings (SA, stock → quarterly average). Series ID format verified: `JTS000000{ST}0000000JOL` (21 chars). Data available from Dec 2000 for all 50 states.
- **Priority 2 completed:** JOLTS decomposed into layoffs+discharges (LD) and quits (quits) via `part2_shock_rates_s.py` and `part3_resid_instruments.py`. Separate Bartik instruments built and residualized for each.
- **Vacancy LP outcome** added to `part5_lp.py` via `outcome="vacancy"` parameter. Uses log-change specification; lagged log vacancies replace lagged unemployment rate as control.

---

## 4. Key Empirical Findings

### 4.1 Shock Persistence (part6)
- **ρ_δ = 0.617** (half-life ≈ 1.44 quarters)
- **ρ_s = 0.751** (half-life ≈ 2.42 quarters)
- Both shocks essentially dissipate by h=10 quarters

### 4.2 Unemployment IRFs (baseline, part5)

**δ shock (baseline, full sample 1992Q3–2019Q4):**
- Significant positive unemployment response, plateaus around h=5–9
- Broadly consistent with model predictions and ρ_δ half-life
- Post-2001 subsample produces similar shape (sample comparability confirmed)

**s shock (baseline, sample 2001Q1–2019Q4):**
- **Anomaly:** IRF rises monotonically through h=16 even though shock dissipates by h=8–10
- Mathematically inconsistent with linear transmission from a shock with ρ_s=0.751
- Anomaly survives NFCI interaction — financial amplification does not explain it

### 4.3 NFCI Interaction Results (part5)

**δ shock × NFCI:**
- β_h (at average financial conditions) is **insignificant at all horizons beyond h=0**
- δ_h (interaction coefficient) is **highly significant (p<0.001) throughout**
- **Interpretation:** The entire δ unemployment effect operates through tight financial conditions. At average conditions, δ shocks have no detectable unemployment effect.
- This is a **puzzle** relative to the model, which has no financial friction.

**s shock × NFCI:**
- β_h grows steadily and becomes significant from h=7 onward even at average financial conditions
- δ_h also significant but smaller
- Monotonically rising s IRF **survives** the interaction — not explained by financial amplification

### 4.4 Residualized Instruments (part2b + part3)

- After purging Δlog(productivity) comovement, δ and s IRFs are on similar orders of magnitude
- Productivity component was inflating raw δ estimates
- Cross-instrument correlation: raw δ vs. s ≈ 0.54; reduced but remains meaningfully positive after residualization
- LD (layoffs+discharges) and δ more correlated than quits and δ, as expected from cyclicality

### 4.5 JOLTS Decomposition — LD vs. Quits (Priority 2 — NEW)

Three unemployment LP specifications now available:
1. **Total separations (s, baseline):** Monotonically rising IRF anomaly (as above)
2. **Layoffs+discharges (LD):** Primary specification; closest to model's s_t. IRF shape closer to economic prior.
3. **Quits (placebo):** Model predicts different profile. Use as identification check.

Cross-instrument correlations and productivity loadings now available from `part2b` — LD and δ share countercyclical variation; quits load positively on productivity (as expected).

### 4.6 Vacancy IRFs (Priority 1 — NEW, KEY RESULTS)

**δ shock → vacancies (residualized instrument, `lp_irf_delta_vacancy.csv`):**
- Vacancies fall **persistently and significantly** from h=1 onward
- β grows in magnitude: h=1: −0.111, h=4: −0.260, h=8: −0.314, h=12: −0.423, h=16: −0.422 log points per 1-SD shock
- All statistically significant (p<0.01 from h=2 onward; p<0.05 at h=1)
- ✅ **Confirms model prediction:** δ destroys vacancies persistently

**LD shock → vacancies (residualized, `lp_irf_ld_vacancy.csv`):**
- Vacancies also fall, significantly (h=0: −0.095, h=5: −0.253, h=6: −0.295)
- Effect weakens after h=8 and becomes less significant at longer horizons
- ⚠️ **Mixed finding:** Model predicts flat or rising vacancies after LD (reposting channel). Empirically vacancies fall, though less persistently than after δ.
- Possible interpretation: reposting is partial, not full. Or LD instrument still captures some product-line destruction.

**Quits shock → vacancies (`lp_irf_qu_vacancy.csv`):**
- Coefficients small and mostly insignificant
- ✅ **Consistent with placebo role:** quits do not strongly move vacancies

**Total separations (s) → vacancies (`lp_irf_ts_vacancy.csv`):**
- Negative and significant through h=3–8; similar to LD pattern
- Confirms contamination from quits reduces the signal

**Key plot:** `data/results/lp_irf_vacancy_delta_ld.png` — side-by-side δ vs. LD vacancy IRFs (the core asymmetry test)

### 4.7 Central Tension (unchanged from prior sessions)

Model predicts: s shocks partially self-correct via reposting; δ generates persistent Beveridge curve dynamics. Empirically:
1. δ effects on unemployment concentrate almost entirely in tight financial conditions (not in model)
2. s (LD) generates monotonically rising unemployment IRF even at average financial conditions (inconsistent with shock persistence)
3. Vacancies fall after both δ and LD shocks — reposting channel appears empirically weaker than model implies

These tensions need direct engagement in the paper text.

---

## 5. Remaining Tasks (Updated)

### Remaining from original Priority List:

**PRIORITY 3 — SLOOS C&I Index as Alternative Interaction Variable** ❌ Not yet done
- Fetch Senior Loan Officer Opinion Survey C&I net tightening index from FRED
- Re-run δ NFCI interaction replacing NFCI_risk with SLOOS C&I
- Tests whether NFCI result is genuinely credit-supply vs. recession severity
- Implementation: `part7b_sloos.py` or extend `part7_nfci.py`

**PRIORITY 4a — Recession-severity placebo for interaction term** ❌ Not yet done
- Run three parallel interaction variants: (1) B × NFCI, (2) B × Δu_nat, (3) both simultaneously
- Tests if NFCI is proxying recession depth rather than financial amplification
- Implementation: extend `part5_lp.py` with new interaction columns

**PRIORITY 4b — Wild cluster bootstrap p-values** ❌ Not yet done
- Current state-clustered SEs with n=50 clusters are at lower bound of reliability
- Supplement existing results with wild cluster bootstrap
- Python: `wildboottest` or `linearmodels` bootstrap option

**PRIORITY 4c — Pre-GFC sample restriction for s** ❌ Not yet done
- Restrict s sample to pre-GFC quarters to test whether monotonically rising IRF is GFC-driven
- Implementation: add `max_qt="2007Q4"` option to `run_lp()` in `part5_lp.py`

### New task arising from vacancy results:

**PRIORITY NEW — Interpret and document LD vacancy finding**
- LD shock lowers vacancies (not flat/rising as model predicts)
- Consider: (a) partial reposting interpretation, (b) LD instrument still correlated with δ shock, (c) financial frictions mute reposting
- Draft paper text explaining the asymmetry between δ and LD vacancy IRFs (δ: persistent decline; LD: moderate, fading decline)
- The δ–LD vacancy asymmetry IS visible — the fall is larger and more persistent for δ

---

## 6. LP Specification Reference

### Unemployment LP (baseline)
```
y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h · B^(k)_{s,t}
                         + γ₁ · u_{s,t-1} + γ₂ · log(LF_{s,t-1}) + ε
```

### Vacancy LP
```
log(v_{s,t+h}) - log(v_{s,t-1}) = α_s + α_t + β_h · B^(k)_{s,t}
                                   + γ₁ · log(v_{s,t-1}) + γ₂ · log(LF_{s,t-1}) + ε
```

### NFCI interaction variant
```
... + β_h · B_{s,t} + δ_h · B_{s,t} × NFCI_risk_dm_t + ...
```
NFCI_risk_dm_t demeaned within estimation sample so β_h = IRF at average financial conditions.

**Key conventions:**
- Base year for shares: 2006
- LOO (leave-one-out) national shock rates to avoid mechanical correlation
- Instruments rescaled ×1000 for correct quarterly rate units
- IRF coefficients standardized by cross-sectional SD of instrument (1-SD shock interpretation)
- Outcome window capped at **2019Q4** (COVID exclusion)
- LP horizons: h = 0, 1, ..., 16 (or 20) quarters
- State FE + time FE; SE clustered at state level (n=50)
- 50 states (DC handled case-by-case depending on data availability)
- State FIPS as zero-padded 2-digit strings throughout

---

## 7. Data Sources

| Source | Variable | Coverage | Notes |
|--------|----------|----------|-------|
| BLS QCEW (annual bulk CSVs) | Employment by state × industry | 2006 (base year) | agglvl=54, 12 supersectors |
| BLS BED | Establishment closings (δ instrument) | 1992Q3–2024Q2 | LOO national rates |
| BLS JOLTS (national) | Total separations, LD, quits (s instrument) | 2001Q1–2023Q1 | LOO national rates |
| BLS JOLTS (state, model-imputed) | Job openings stock (vacancy outcome) | Dec 2000–present, 50 states + DC | SA series; quarterly = average of 3 months |
| BLS LAUS | State unemployment rate + labor force | 1990M1–present, 50 states | Monthly → quarterly |
| Chicago Fed NFCI | Financial conditions, risk subindex | 1990Q1+ | Weekly → quarterly; demeaned in sample |
| BLS OPHNFB | Nonfarm business productivity | Quarterly | Used in part2b productivity residualization |

---

## 8. File Index (Key Files)

### Instruments
- `data/instruments/delta_instrument_resid_base2006.csv` — primary δ instrument (residualized)
- `data/instruments/ld_instrument_resid_base2006.csv` — primary LD instrument (residualized)
- `data/instruments/qu_instrument_resid_base2006.csv` — quits instrument (residualized, placebo)
- `data/instruments/s_instrument_base2006.csv` — total separations instrument (raw)
- `data/instruments/laus_quarterly.parquet` — state outcomes: unemployment + labor force + vacancies

### LP Results
- `data/results/lp_irf_delta.csv` — δ → unemployment (full sample)
- `data/results/lp_irf_delta_post2001.csv` — δ → unemployment (post-2001)
- `data/results/lp_irf_delta_resid.csv` — δ residualized → unemployment
- `data/results/lp_irf_delta_nfci.csv` — δ × NFCI interaction
- `data/results/lp_irf_s.csv` — total s → unemployment
- `data/results/lp_irf_s_resid.csv` — s residualized → unemployment
- `data/results/lp_irf_s_nfci.csv` — s × NFCI interaction
- `data/results/lp_irf_ld_resid.csv` — LD residualized → unemployment
- `data/results/lp_irf_qu_resid.csv` — quits residualized → unemployment
- `data/results/lp_irf_delta_vacancy.csv` — δ → log vacancies ⭐
- `data/results/lp_irf_ld_vacancy.csv` — LD → log vacancies ⭐
- `data/results/lp_irf_qu_vacancy.csv` — quits → log vacancies (placebo) ⭐
- `data/results/lp_irf_ts_vacancy.csv` — total s → log vacancies

### Key Plots
- `data/results/lp_irf_vacancy_delta_ld.png` — core asymmetry test (δ vs. LD vacancy IRFs) ⭐
- `data/results/lp_irf_vacancy_decomp_overlay.png` — all shock types vacancy IRFs overlaid
- `data/results/lp_irf_combined.png` — unemployment IRFs combined
- `data/results/lp_irf_delta_ld_qu_comparison.png` — unemployment IRFs for LD decomposition
- `data/results/lp_irf_beveridge_asymmetry.png` — Beveridge asymmetry summary
- `data/results/shock_persistence.csv` — ρ_δ=0.617, ρ_s=0.751

### Model files (outside working directory)
- `Draft_21Feb2026.tex` — model draft
- `steady_state_refactored.jl` — Julia steady state code

---

## 9. Principles for Any Agent Working on This Project

1. **δ vs. s asymmetry is the core.** Every empirical choice should be evaluated against testing this asymmetry.
2. **Vacancy LP is now available.** Use it. The δ vacancy decline is the strongest confirmation of the model's mechanism; the LD partial decline is a puzzle to address.
3. **LD is the primary s-type instrument** (not total separations). Quits are placebo.
4. **No financial frictions in the model.** Do not introduce them beyond the NFCI interaction term.
5. **LOO always.** National shock rates must leave out the state being instrumented.
6. **COVID cap: 2019Q4** for all LP outcomes. Instruments can extend further.
7. **Single-responsibility pipeline.** Data fetching, instrument construction, and regression are in separate scripts.
8. **State FIPS as zero-padded 2-digit strings** throughout.
9. **Vacancy aggregation = average** (stock measure), not sum. Separation aggregation = sum (flow measure).
10. **The δ financial state-dependence and the s monotonic rise are puzzles**, not confirmations. Flag them in the paper text.
