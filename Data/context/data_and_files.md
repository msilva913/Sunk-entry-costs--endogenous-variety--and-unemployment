# Data Sources and File Index
**Last updated:** September 5, 2026

## Data Sources

| Source | Variable | Coverage | Notes |
|--------|----------|----------|-------|
| BLS QCEW (annual bulk CSVs) | Employment by state×industry | 2006 (base year) | agglvl=54, 12 supersectors |
| BLS BED dataclass=06 | Establishment closings | 1992Q3–2024Q2 | Includes temporary shutdowns |
| BLS BED dataclass=08 | Establishment deaths | 1992Q3–2021Q4 | Permanent exits; Deaths ⊆ Closings |
| BLS JOLTS (national) | Total separations, LD, quits | 2001Q1–2023Q1 | LOO national rates |
| BLS JOLTS (state) | Job openings stock (vacancy outcome) | Dec 2000–present, 50 states | SA; quarterly = avg of 3 months |
| BLS LAUS | State unemployment rate + labor force | 1990M1–present | Monthly → quarterly |
| BLS OPHNFB | Nonfarm business productivity | Quarterly | v2 residualization control |
| BEA via FRED (RVAM, RVAC, …) | Real value-added by supersector | 2005Q1–present | Extended to 1992Q1 via Chow-Lin |
| Chicago Fed NFCI | Financial conditions, risk subindex | 1990Q1+ | Weekly → quarterly; demeaned in sample |
| Barnichon (2010) Composite HWI | Pre-JOLTS vacancy proxy | 1951–2007 | Spliced with JOLTS at 2001Q1 |

**BEA API is blocked** (both sandbox and local). Use FRED exclusively. FRED series IDs:

| BLS Code | FRED Series | Sector |
|----------|------------|--------|
| 10 | RVAM | Mining |
| 20 | RVAC | Construction |
| 30 | RVAMA | Manufacturing |
| 41 | RVAW | Wholesale |
| 42 | RVAR | Retail |
| 43 | RVAT + RVAU | Transport+Utilities |
| 50 | RVAI | Information |
| 55 | RVAFI + RVARL | Finance+RE |
| 60 | RVAPBS | Prof & Business Services |
| 65 | RVAES + RVAHC | Edu+Health |
| 70 | RVAER + RVAAF | Entertainment+Food |
| 80 | RVAOSEG | Other Services |

## LP Specification Reference

Generic estimating equation (both outcomes):
```
y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h · B̃^δ_{s,t} + γ₁ · y_{s,t-1} + γ₂ · log(L_{s,t-1}) + ε
```
- `y = u_{s,t}`: unemployment rate (LAUS, pp) — eq:lp in draft
- `y = v_{s,t} = 100 × V/L`: vacancy rate (JOLTS/LAUS, pp) — **NOT log vacancy**
- Both outcomes in percentage points; β_h is pp change per 1-SD shock
- β_h scaled by `instr_sd = 17.4 pp` (SD of B̃^δ; the 14.0 value was the pre-May-14
  closings×π instrument and is stale wherever it still appears)
- Sample: 2001Q1–2019Q4; state+time FE; SE clustered at state (n=50)
- h = 0, 1, …, 20 quarters

NFCI interaction variant:
```
… + β_h · B_{s,t} + δ_h · B_{s,t} × NFCI_risk_dm_t + …
```
(NFCI_risk_dm demeaned in sample → β_h = IRF at average financial conditions)

## Key Instruments (data/instruments/)

| File | Description |
|------|-------------|
| `delta_instrument_resid_base2006.csv` | Primary δ (v2 residualized) |
| `ld_instrument_resid_base2006.csv` | Primary LD (v2 residualized) |
| `qu_instrument_resid_base2006.csv` | QU placebo (v2 residualized) |
| `s_instrument_base2006.csv` | Total separations (raw, appendix only) |
| `laus_quarterly.parquet` | Outcomes: u_rate, vac_rate, labor_force, 50 states |

## Key LP Results (data/results/)

| File | Description |
|------|-------------|
| `lp_irf_delta_resid.csv` | **BASELINE**: δ → unemployment (2001Q1+, v2) |
| `lp_irf_delta_vacancy.csv` | **BASELINE**: δ → vacancy rate (2001Q1+, v2) |
| `lp_irf_delta_uv_baseline.png` | **PAPER FIGURE** fig:delta_uv_irf — side-by-side u+v |
| `lp_irf_delta_resid_1997.csv` | δ → unemployment (1997Q1+, sample extension) |
| `lp_irf_ld_resid.csv` | LD → unemployment |
| `lp_irf_ld_vacancy.csv` | LD → vacancy rate |
| `lp_irf_qu_resid.csv` | QU → unemployment (placebo — FAILS) |
| `lp_irf_qu_vacancy.csv` | QU → vacancy rate (placebo — FAILS) |
| `lp_joint_delta_ld_unemp.csv` | Joint LP: δ+LD → unemployment |
| `lp_joint_delta_ld_vacancy.csv` | Joint LP: δ+LD → vacancy rate |
| `placebo_sev_delta_unemp.csv` | Severity placebo: δ → u |
| `placebo_sev_delta_vac.csv` | Severity placebo: δ → v (cleanest result) |
| `shock_persistence.csv` | AR(1): ρ_δ=0.617/0.647, ρ_LD=0.489/0.305 |
| `var_calibration.csv` | VAR(1): ρ_δ=0.600, ρ_LD=0.510, β(endex)=0.228 |
| `var_calibration_zd.csv` | **MODEL CALIB**: bivariate VAR(1) on HP-1600 log(z),log(δ); Cholesky z-first; ρ_z^Q=0.735→ρ_z^m=0.902, ρ_δ^Q=0.207→ρ_δ^m=0.592, σ_z^m=0.0092, σ_δ^m=0.0669 |
| `recession_scatter_primary.png` | Cross-recession fig (paper) |
| `state_scatter_primary.png` | Cross-state fig (paper) — fig:state_scatter |
| `jf_table2_ext_2019_latex.tex` | JF Table 2 extension (preferred sample) |

## SMM Moment Construction

**Baseline filter: HP with λ=1,600** (log-levels; applied identically to data and model simulations).
- Switched from λ=100,000 (Shimer/CK convention) on May 19, 2026
- Reason: λ=100,000 reverses sign of cor(δ,u) and cor(δ,v) — artifact of over-smoothing, not data feature
- Robustness documented in `app:filter_robustness` (Table B.6): HP-1600 and Hamilton agree; HP-100k is outlier at 7/18 cells
- Script: `observables_moments.py` using `filter_transform(..., lamb=1600)` from `time_series_functions.py`
- Output pickle: `raw_data.pkl` (delta ends 2021Q4 — stale BED cache; refresh when re-estimating)

**New helper script:** `labor_market_dyn_corr_by_filter.py` — produces `Bartek analysis/data/results/dyn_corr_delta_u_v_by_filter.csv` and plots

## Calibration Targets

| Parameter | Value | Source |
|-----------|-------|--------|
| τ̄ | 9.30%/qtr (3.10%/month) | Shimer (2012) total separation rate; `TARGETS.sep = 0.031` |
| δ̄_e (deaths, 2001–2019) | **0.718%/qtr (0.239%/month; 2.84%/yr)** | BLS BED published rate `BDS…110008RQ5`, total private |
| δ̄_e (deaths, 1993–2019) | 0.812%/qtr (3.21%/yr) | same series, longer window |
| **δ_e/τ** | **0.077 (2001–2019) / 0.087 (1993–2019)** | **Recommended primary target — see [D1](decisions.md)** |
| BED closings rate | 1.205%/qtr (2001–2019) | `…110006RQ5`. Includes temporary shutdowns — *not* the model object |
| BED gross job losses | 6.362%/qtr (2001–2019) | `…110004RQ5`. Denominator of JF Table 2 col 2 |
| closings / gross losses | 18.9% (2001–2019), 19.7% (1993–2019) | Reproduces JF's ~21% and `part11`'s 19.5% |
| deaths / gross losses | ~12% | The deaths analogue of JF's statistic |

> ⚠️ The previously recorded δ̄ = 1.079%/qtr (δ/τ = 11.6%) is **too high**: `observables.py`
> runs to `final = '2025-10-01'`, so its mean includes the 2020 COVID quarters, and it divides
> by PAYEMS (total nonfarm) rather than private employment. Recompute pre-COVID before the
> number enters `tab:calib_targets`. The 4.13%/yr figure that briefly appeared here was an
> arithmetic error and has no standing — see the correction note in [D1](decisions.md).

| ρ_δ | 0.617 (raw) / 0.647 (resid) | part6, 2001Q1+ window, Bartik agg |
| ρ_LD | 0.489 (raw) / 0.305 (resid) | part6, 2001Q1+ window |
| corr(η^δ, η^LD) | +0.397 (raw) / +0.488 (resid) | part6 |
| ρ_δ (VAR industry) | 0.600 | part2d, industry-level panel |
| β (LD←δ) | +0.228 | part2d companion matrix |
| **ρ_z^m (model)** | **0.902** | part6b; VAR(1) HP-1600 log(z), Cholesky, 1992Q3–2019Q4 |
| **ρ_δ^m (model)** | **0.592** | part6b; VAR(1) HP-1600 log(δ) ⊥ z, Cholesky z-first |
| **σ_z^m (model)** | **0.0092** | part6b; unconditional-variance-matched from quarterly |
| **σ_δ^m (model)** | **0.0669** | part6b; structural δ shock ⊥ z innovation |
| β(δ←z) quarterly | −1.050 (SE=0.830) | part6b; endogenous exit channel — NOT in shock process |
| corr(u^z, u^δ) | −0.284 | part6b; pre-Cholesky; absorbed by model equilibrium |

> ⚠️ **δ̄ is contested — see [D1](decisions.md).** The table above reports the BED Deaths
> measurement (4.25%/yr, δ_e/τ = 0.116). `steady_state.jl:566` instead sets
> `dest_ann = 0.0754` — labelled "[BED Deaths, emp-weighted]" but actually the Jaimovich-Floetotto
> 21% × τ figure (0.651%/month, 1.94%/qtr, δ_e/τ = 0.210) — and every mechanism figure was
> produced at that value. Draft §5.1 cites the former, §5.2 the latter. Resolve before
> estimating.

> ⚠️ **There is no single canonical target set — see [D2](decisions.md).**
> `TARGETS` in `steady_state.jl` defaults to PATH A with `b_ratio = 0.71, x_v = 1.0`;
> the mechanism comparisons override to PATH B (`dest_elast_target = 5.0`) with
> `b_ratio = 0.9, x_v = 0.5`; the draft documents PATH A. Once D2 is settled, write the
> winning target set here and treat it as canonical.

## Model code

Moved to [`pipeline.md`](pipeline.md) Part 2 (model conventions, equation order,
`steady_state.jl` architecture, PATH A/B), so that all "how the code works" material lives
in one file. This file is data sources, file locations, and empirical targets only.
