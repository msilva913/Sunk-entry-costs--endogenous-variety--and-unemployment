# Pipeline — Full Status
**Last updated:** May 12, 2026

Working directory: `Data/Bartek analysis/`
Subdirectories: `data/cache/` · `data/instruments/` · `data/results/`

## Scripts (run in order)

| Script | Function | Status | Key Output |
|--------|----------|--------|------------|
| `part1_shares.py` | 2006 QCEW state×supersector shares ω_{s,j} | ✅ | `shares_base2006.parquet` |
| `part2_shock_rates.py` | BED deaths → LOO national g^δ | ✅ | `shock_rates_1992Q3_2024Q2.parquet` |
| `part2_shock_rates_s.py` | JOLTS separations → LOO g^LD, g^QU, g^TS | ✅ | `shock_rates_{ld,qu,s}_*.parquet` |
| `part2b_residualize_shocks.py` | v2 residualization: regress on Δlog p_t + Δlog VA_{j,t-1} | ✅ | `shock_rates_{delta,ld,qu}_resid.parquet` |
| `part2c_granger_lp.py` | Panel LP Granger causality δ vs LD | ✅ | `granger_lp_*.png` |
| `part2d_var_calibration.py` | Bivariate panel VAR(1) on (ν^δ, ν^LD) | ✅ | `var_calibration.csv`, `var_irf_chol.png` |
| `part3_instrument.py` | Bartik aggregation for δ (raw) | ✅ | `delta_instrument_base2006.csv` |
| `part3_instrument_s.py` | Bartik aggregation for TS (raw) | ✅ | `s_instrument_base2006.csv` |
| `part3_resid_instruments.py` | Bartik aggregation for ν^δ, ν^LD, ν^QU | ✅ | `{delta,ld,qu}_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | LAUS unemployment + JOLTS state vacancies | ✅ | `laus_quarterly.parquet` |
| `part5_lp.py` | All panel LP regressions, both outcomes | ✅ | IRF CSVs + plots in `data/results/` |
| `part5_joint_lp.py` | Joint LP: δ and LD in same regression | ✅ | `lp_joint_delta_ld_{unemp,vacancy}.{csv,png}` |
| `part6_shock_persistence.py` | AR(1) shock persistence | ✅ | `shock_persistence.csv` |
| `part7_nfci.py` | Chicago Fed NFCI → quarterly | ✅ | `nfci_quarterly.parquet` |
| `part7c_recession_placebo.py` | Dual-severity placebo (u^nat level + change) | ✅ | `placebo_sev_*.csv`, `placebo_sev_*.png` |
| `part9_recession_scatter.py` | Cross-recession motivating scatter | ✅ | `recession_scatter_{primary,appendix}.png` |
| `part10_state_scatter.py` | Cross-state motivating scatter (50 states) | ✅ | `state_scatter_{primary,appendix}.png`, `state_scatter.csv` |
| `part10b_state_scatter_sensitivity.py` | X-axis sensitivity: full vs recession-quarter | ✅ | `state_scatter_sensitivity.png` |
| `part11_jf_table2.py` | JF (2008) Table 2 replication + extension | ✅ | `jf_table2_*.{parquet,tex,txt}` |
| `construct_delta_instrument.py` | BED Deaths vs Closings diagnostic plot | ✅ | `bed_closings_vs_deaths.png` |

## Run Order for Fresh Start
```
part1 → part2 → part2b → part3_resid → part4 → part5 → part6
```
For v2 residualization (requires FRED key):
```powershell
$env:FRED_API_KEY = "<key>"
python run_v2_locally.py   # orchestrates part2b → part3_resid → part5
```
After first run, `data/cache/bea_va_quarterly_12ind.parquet` is cached; key not needed again.

## Key Technical Notes
- **BED Deaths (dataclass=08)**: establishments absent 4+ consecutive quarters — permanent exits. Used as δ numerator. Deaths ⊆ Closings by construction. Deaths/Closings wtd mean = 0.90 (2001–2019).
- **LOO**: national shock rates exclude the state being instrumented (denominator only — state-level deaths by supersector not publicly available).
- **Base year**: 2006 QCEW shares (agglvl=54, private sector).
- **v2 residualization**: `log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t-1} + ν^k_{j,t}`. VA from FRED (RVAM, RVAC, etc.), extended back to 1992Q1 via Chow-Lin + Denton adjustment.
- **State FIPS**: zero-padded 2-digit strings throughout.
- **Vacancy aggregation**: average of monthly stocks (not sum). Separation aggregation: sum (flow).
- **COVID cap**: 2019Q4 for all LP outcomes.
