# Data Sources and File Index
**Last updated:** May 12, 2026

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
- β_h scaled by `instr_sd = 14.0 pp` (cross-sectional SD of B̃^δ)
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
| δ̄ | 0.940%/qtr (0.313%/month) | BDS employment-weighted exit rate |
| τ̄ | 9.34%/qtr (3.11%/month) | Shimer (2012) total separation rate |
| δ/τ | 10.1% | Ratio above |
| ρ_δ | 0.617 (raw) / 0.647 (resid) | part6, 2001Q1+ window |
| ρ_LD | 0.489 (raw) / 0.305 (resid) | part6, 2001Q1+ window |
| corr(η^δ, η^LD) | +0.397 (raw) / +0.488 (resid) | part6 |
| ρ_δ (VAR) | 0.600 | part2d, industry-level panel |
| β (LD←δ) | +0.228 | part2d companion matrix |
