# Data Sources and File Index
**Last updated:** September 23, 2026

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
| `shock_calibration.csv` / `.txt` | **MODEL CALIB** (S14, Oct 28): one row per series for z, δ, τ, s, with each row's sample window. z, τ, s univariate AR(1) on 1951–2019 (s spliced before 1992Q3); δ from the (z, δ) VAR on 1992Q3–2019Q4, σ_δ ⊥ z. Monthly: ρ_z 0.9121, σ_z 0.0128; ρ_δ 0.5824, σ_δ 0.0663; ρ_s 0.8054, σ_s 0.0711. τ is a diagnostic, not a shock. The `.txt` is the full report (method, splice checks, D14 test) |
| `q3_data_svar.txt` / `_fevd.csv` / `_irf.csv` / `_irf.png` | **Q3 Block B′ data side** (Oct 29): SVARs A/B/C, variance shares of u and v with ordering bounds, IRFs with bootstrap bands. See findings.md "Q3: data SVAR" |
| `cycles_panel.png` | HP-1600 log cycles for all four series, each on its own window (z, τ, s 1951–2019; δ 1992Q3–2019Q4); dotted line at the splice |
| `separation_decomposition.png` | τ, δ, s in levels and cycles, 1951–2019, with the δ̄ stand-in before the 1992Q3 splice; τ = δ + (1−δ)s. cor(cycle_s, cycle_τ) = 0.995 on the BED window |
| `recession_scatter_primary.png` | Cross-recession fig (paper) |
| `state_scatter_primary.png` | Cross-state fig (paper) — fig:state_scatter |
| `jf_table2_ext_2019_latex.tex` | JF Table 2 extension (preferred sample). LaTeX tables unchanged by the Oct 5 edit |
| `jf_table2_extension.txt` | Text report. **Oct 5, 2026:** gains C7–C10, the entry-comovement correlations (openings vs deaths, closings, net job change; and deaths vs net job change), under both filters. Correlations are **text-output only** — deliberately not added to the LaTeX tables, which are already ten columns wide |

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
| τ̄ | 9.30%/qtr (3.10%/month) | Shimer (2012) total separation rate; `TARGETS.sep = 0.031`. ⚠ **To recompute (Oct 28):** τ was rebuilt (CPS 1994 factor and timing fixes); means are now 0.0356 (1951–93) and 0.0258 (1992Q3–2019Q4) |
| δ̄_e (deaths, 2001–2019) | **0.718%/qtr (0.239%/month; 2.84%/yr)** | BLS BED published rate `BDS…110008RQ5`, total private |
| δ̄_e (deaths, 1993–2019) | 0.812%/qtr (3.21%/yr) | same series, longer window |
| **δ_e/τ** | **0.087 (1993–2019) — CANONICAL** | **Settled by [D1](decisions.md), Sept 6, 2026.** `dest_ann = 0.0320`. The 2001–2019 window gives 0.077; MS chose the longer window |
| BED closings rate | 1.205%/qtr (2001–2019) | `…110006RQ5`. Includes temporary shutdowns — *not* the model object |
| BED gross job losses | 6.362%/qtr (2001–2019) | `…110004RQ5`. Denominator of JF Table 2 col 2 |
| closings / gross losses | 18.9% (2001–2019), 19.7% (1993–2019) | Reproduces JF's ~21% and `part11`'s 19.5% |
| deaths / gross losses | ~12% | The deaths analogue of JF's statistic |

> ⚠️ **Retracted measurements, kept so they are not resurrected.** δ̄ = 1.079%/qtr
> (δ/τ = 11.6%) is **too high**: `observables.py` runs to `final = '2025-10-01'`, so its mean
> includes the 2020 COVID quarters, and it divides by PAYEMS (total nonfarm) rather than
> private employment. The 4.13%/yr figure that briefly appeared here was an arithmetic error
> and has no standing. See the correction note in [D1](decisions.md).
>
> ⏳ **Open N15 obligation.** The canonical 0.087 must be emitted by a program before it
> enters `tab:calib_targets`. Add a pre-COVID calibration-target computation to
> `observables.py`; today the number is sourced only to the BLS published series, not to the
> replication workflow.

| ρ_δ | 0.617 (raw) / 0.647 (resid) | part6, 2001Q1+ window, Bartik agg |
| ρ_LD | 0.489 (raw) / 0.305 (resid) | part6, 2001Q1+ window |
| corr(η^δ, η^LD) | +0.397 (raw) / +0.488 (resid) | part6 |
| ρ_δ (VAR industry) | 0.600 | part2d, industry-level panel |
| β (LD←δ) | +0.228 | part2d companion matrix |
| **ρ_z^m (model)** | **0.9121** | part6b S14 (Oct 28); AR(1) on HP-1600 log z, 1951Q1–2019Q4 |
| **σ_z^m (model)** | **0.0128** | part6b S14; AR(1) residual SD, unconditional-variance-matched to monthly |
| **ρ_δ^m (model)** | **0.5824** | part6b S14; own lag of the (z, δ) VAR(1), 1992Q3–2019Q4 |
| **σ_δ^m (model)** | **0.0663** | part6b S14; VAR residual SD ⊥ the z innovation (Cholesky z-first; 4.3 % smaller than reduced form) |
| **ρ_s^m (model)** | **0.8054** | part6b S14; univariate AR(1) on spliced s, 1951Q1–2019Q4 |
| **σ_s^m (model)** | **0.0711** | part6b S14; AR(1) residual SD, not purged (D14) |
| s̄, τ̄, δ̄_e (BED window) | 0.02356 / 0.02576 / 0.00225 per month | part6b S14, 1992Q3–2019Q4 means (pre-1992 s levels are imputed: never use them for targets) |
| β(δ←L.z) quarterly | −1.192 (t = −1.36) | part6b S14 VAR; endogenous exit channel, NOT in the shock process |
| corr(u^z, u^δ) | −0.290 | part6b S14; before the Cholesky step |

> **`SHOCKS_BASE` holds these S14 values since Oct 29** (before: ρ_z 0.902, σ_z 0.0092;
> ρ_δ 0.592, σ_δ 0.0669; ρ_s 0.8741, σ_s 0.0854). The Oct 20 interim
> numbers (S13 rerun, part6c trivariate) are superseded. Method:
> [`session_handout_20261028.md`](session_handout_20261028.md) §2.
>
> **LP series:** PRS85006163 is nonfarm-business **output per worker** (Shimer/HM/CK
> standard), not per hour.
>
> **Block M targets:** the σ(u) 0.1496 / cor(u,v) −0.804 values cover 1951–2024, including the
> pandemic. Under S13/S14 they become 1951–2019: σ(u) 0.1248, σ(v) 0.1382, cor(u,v) −0.925.
> The full table has not been regenerated yet.

> ✅ **δ̄ is settled — [D1](decisions.md), September 6, 2026.** BED establishment deaths,
> employment-weighted, 1993–2019: `dest_ann = 0.0320` (3.21%/yr), **δ_e/τ = 0.087**.
>
> ⏳ **But not yet implemented.** `steady_state.jl:623` still sets `dest_ann = 0.0754` — the
> Jaimovich-Floetotto 21% × τ figure (δ_e/τ = 0.210), which mixes a gross-job-loss share with
> a separation rate — behind a comment saying so. Every mechanism figure and every §5.3
> number is still at 0.210. The switch is the "Step 0 cascade" in
> [`pending_tasks.md`](pending_tasks.md); it is bundled with the ρ_s/σ_s calibration so one
> regeneration covers both.

> ⚠️ **There is no single canonical target set — see [D2](decisions.md).**
> `TARGETS` in `steady_state.jl` defaults to PATH A with `b_ratio = 0.71, x_v = 1.0`;
> the mechanism comparisons override to PATH B (`dest_elast_target = 5.0`) with
> `b_ratio = 0.9, x_v = 0.5`; the draft documents PATH A. Once D2 is settled, write the
> winning target set here and treat it as canonical.

## Model code

Moved to [`pipeline.md`](pipeline.md) Part 2 (model conventions, equation order,
`steady_state.jl` architecture, PATH A/B), so that all "how the code works" material lives
in one file. This file is data sources, file locations, and empirical targets only.
