# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)

**Last updated:** April 14, 2026
**Author:** Mario Silva
**Purpose:** Persistent project context for fresh Cowork sessions. Paste this into any new session to restore full project state.
**Code structure:** clean, succinct code easily interpretable by empirical macroeconomist.

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
| `part2b_shock_comovement.py` | v3 residualization: regress g^k on controls → residuals ν^k | ✅ Complete | `shock_rates_delta_resid.parquet`, `shock_rates_ld_resid.parquet`, `shock_rates_qu_resid.parquet` |
| `part2b_extend_va.py` | Constrained Chow-Lin backcast of sector VA pre-2005; annual benchmarks from BEA CSV | ✅ Complete | `data/cache/bea_va_quarterly_12ind_extended.parquet` |
| `part3_instrument.py` | Bartik aggregation for δ (raw, kept for archival) | ✅ Complete | `data/instruments/delta_instrument_base2006.csv` |
| `part3_resid_instruments.py` | Bartik aggregation for residualized ν^δ, ν^{LD}, ν^{QU} | ✅ Complete | `delta_instrument_resid_base2006.csv`, `ld_instrument_resid_base2006.csv`, `qu_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | Fetches LAUS unemployment + JOLTS state vacancies → quarterly panel | ✅ Complete | `data/instruments/laus_quarterly.parquet` (includes `vacancies` column) |
| `part5_lp.py` | All panel LP regressions: residualized instruments only; GFC interaction diagnostic added | ✅ Complete | IRF CSVs + plots in `data/results/` |
| `part6_shock_persistence.py` | AR(1) estimation of shock persistence | ✅ Complete | `data/results/shock_persistence.csv` |
| `part7_nfci.py` | Chicago Fed NFCI → quarterly, headline + risk subindex | ✅ Complete | `data/cache/nfci_quarterly.parquet` |
| `part2c_granger_lp.py` | Panel LP Granger causality: δ vs LD (raw + residualized) | ✅ Complete | `data/results/granger_lp_*.png` |
| `part2d_var_calibration.py` | Bivariate panel VAR(L) on (ν^δ, ν^{LD}): companion matrix, Σ, Cholesky IRFs | ✅ Complete | `data/results/var_calibration.csv`, `var_irf_chol.png` |
| `part5_joint_lp.py` | Joint LP with δ and LD simultaneously (v2-residualized); robustness check for δ | ✅ Complete | `lp_joint_delta_ld_{unemp,vacancy}.csv` + `.png` |
| `part7b_sloos.py` | SLOOS C&I net tightening × δ interaction LP; unemployment and vacancy outcomes | ✅ Complete | `lp_irf_delta_sloos_lm_{unemp,vacancy}.csv` + comparison plots |

### Notable additions in April 14, 2026 session:
- **Constrained Chow-Lin backcast** (`part2b_extend_va.py`) fully rewritten with annual benchmark constraints from BEA annual industry VA file. See section 15 for full methodology. Annual benchmark CSV: `bea_va_annual_by_supersector.csv` (BEA Table 1.3.6, April 9 2026, coverage 1997–2025). Extended quarterly cache: `data/cache/bea_va_quarterly_12ind_extended.parquet`.
- **FRED series corrections** in `BLS_TO_FRED` mapping: `RVAESHS` replaces `RVAES` for BLS65 (Education+Health+Social Assistance full aggregate; `RVAES` was education subsector only, ~16% of total). `RVAAER+RVAAF` replaces `RVAER+RVAAF` for BLS70 (`RVAER` does not exist on FRED).
- **`part5_lp.py` major refactor**: (1) all raw instrument LPs removed — raw δ, raw s-TS, NFCI interaction on raw instruments all dropped; residualized instruments are now the only specification; (2) GFC interaction diagnostic added — two new regressors B_{s,t}×GFC_t (shock-quarter dummy) and GFC_{t+h} (outcome-quarter dummy) — see section 15.4 for rationale.
- **QU placebo failure structural interpretation**: QU vacancy effect is large (−0.049 at h=15, larger than δ peak), onset delayed to h=8 — consistent with industry demand contamination via GFC timing, not structural reposting channel. See section 15.3.
- **Annual benchmark file**: `bea_va_annual_by_supersector.csv` lives in `Data/Bartek analysis/` (same directory as scripts). Required by `part2b_extend_va.py`. Do not move.

### Notable additions in April 9, 2026 session:
- **`part7b_sloos.py`** written and run. Key finding: SLOOS interaction near-zero for unemployment h=0–8, weakly positive h=12–16; vacancy SLOOS interaction uniformly insignificant. NFCI shows sign-switching (recession-severity confound). SLOOS confirms δ effect operates at average credit conditions — consistent with frictionless model. See section 14.1.
- **v3 residualization** adds φ_j × ΔMP_t. r(δ,LD) unchanged (0.423→0.4243). Monetary policy differential sensitivity does not explain residual correlation. See section 14.2.
- **Wu-Xia shadow FFR** extended with FEDFUNDS post-2022M2. Cached at `data/cache/wuxia_quarterly.parquet`. Delete cache to rebuild.
- **Census CBP establishment size** (year 2000): φ_j range 1.990–3.839 log pts. Cached at `data/cache/cbp_estab_size.parquet`.

### Notable additions in April 8, 2026 session:
- **`part5_joint_lp.py`**: joint LP with δ and LD. LD vacancy collapses to zero once δ controlled. δ vacancy stable (diff < 0.13 pp). β^δ shrinks in unemployment equation → δ absorbing shared variation, not clean structural separation. See section 13.

---

## 4. Key Empirical Findings

### 4.1 Shock Persistence
- ρ_δ = 0.479 (residualized; raw 0.617), ρ_LD = 0.288 (residualized; raw 0.489)
- corr(η^δ, η^LD) = +0.444 (residualized) — model zero-correlation assumption violated

### 4.2 Unemployment IRFs (residualized, post-2001)
- **δ**: significant positive unemployment response, peak ~+2.8 pp at h=12, significant h=0–20
- **LD**: peak ~+1.7 pp at h=17, significant throughout, but monotonically rising past shock half-life (anomalous — likely GFC-dominated)
- **QU (placebo)**: peak +1.3 pp at h=15, significant from h=4 — placebo fails

### 4.3 Vacancy IRFs (residualized, post-2001) ⭐ Key results
- **δ**: persistent significant decline h=1–20, peak −0.032 pp at h=12; significant p<0.001 throughout
- **LD**: significant decline h=0–12, peak −0.016 pp at h=5–6; fades h=13–14, returns h=18–20 (partial reposting pattern but contaminated)
- **QU (placebo)**: FAILS — near-zero h=0–6, then large and rising h=8–20, peak −0.049 pp at h=15 (LARGER than δ); onset delay is structural clue (see section 15.3)

### 4.4 SLOOS and Credit Supply Tests
- SLOOS δ_h: near-zero h=0–8, weakly positive h=12–16 (small, delayed credit channel)
- SLOOS vacancy δ_h: uniformly insignificant — vacancy destruction not credit-mediated
- Conclusion: δ effect operates at average credit conditions; consistent with frictionless model

### 4.5 Instrument Correlations (v3 residualization)
- r(δ,LD) = 0.424 — survives v1→v2→v3 residualization unchanged
- r(δ,QU) = 0.048 — QU orthogonalized successfully
- r(LD,QU) = −0.076

---

## 5. Remaining Tasks (Updated April 14)

**PRIORITY 1 — Run GFC interaction LPs and interpret results**
- `part5_lp.py` now produces GFC-controlled IRFs for all three instruments × two outcomes
- Key diagnostic question: does QU vacancy IRF go flat under GFC controls?
  - If YES → placebo failure is entirely GFC-driven, paper can state this cleanly
  - If NO → contamination is structural, pre-dates crisis; also an informative finding
- Compare β_h(QU, vacancy) baseline vs GFC-controlled; plot in `lp_irf_gfc_qu_comparison.png`

**PRIORITY 2 — Wild cluster bootstrap (Priority 4b)**
- n=50 clusters at lower bound of reliability for cluster-robust SEs
- Use `wildboottest` Python package

**PRIORITY 3 — Recession-severity placebo (Priority 4a)**
- Run (1) B×NFCI, (2) B×Δu_nat, (3) both simultaneously
- Formally confirms NFCI proxies recession depth; partially superseded by SLOOS but still useful

**PRIORITY 4 — Paper text**
- Draft empirical section: v1→v2→v3 as systematic confound elimination → r(δ,LD) is structural
- SLOOS robustness confirms frictionless transmission; GFC diagnostic characterises QU contamination
- δ vacancy finding is the headline structural result

---

## 6. LP Specification Reference

### Baseline unemployment LP
```
y_{s,t+h} - y_{s,t-1} = α_s + α_t
                         + β_h · B^(k)_{s,t}
                         + γ₁ · u_{s,t-1} + γ₂ · log(LF_{s,t-1}) + ε_{s,t,h}
```

**Economic rationale:** β_h is the cumulative impulse response of state unemployment to a 1-SD Bartik shock at horizon h quarters. α_s absorbs permanent state-level differences (right-to-work laws, industry composition). α_t absorbs aggregate time series variation common to all states (business cycle, national policy). State-clustered SEs account for within-state serial correlation in outcomes.

### Baseline vacancy LP
```
vac_rate_{s,t+h} - vac_rate_{s,t-1} = α_s + α_t
                                       + β_h · B^(k)_{s,t}
                                       + γ₁ · vac_rate_{s,t-1} + γ₂ · log(LF_{s,t-1}) + ε
```

vac_rate = vacancies × 1000 / labor_force × 100 (pp units). β_h < 0 means shocks destroy vacancy slots.

### GFC interaction LP (diagnostic for QU placebo failure)
```
y_{s,t+h} - y_{s,t-1} = α_s + α_t
                         + β_h · B^(k)_{s,t}
                         + δ_h · B^(k)_{s,t} × GFC_t        [shock-quarter dummy]
                         + φ_h · GFC_{t+h}                   [outcome-quarter dummy]
                         + γ₁ · y_{s,t-1} + γ₂ · log(LF_{s,t-1}) + ε
```

GFC_t = 1 if shock quarter ∈ {2008Q3–2009Q4}; GFC_{t+h} = 1 if outcome quarter ∈ {2008Q3–2009Q4}.

**Economic rationale for two separate dummies:**
- `B × GFC_t` absorbs differential amplification of shocks occurring *during* the GFC (e.g. exceptionally large layoff shocks in 2009 driving unusually large IRFs)
- `GFC_{t+h}` absorbs the aggregate vacancy/unemployment collapse *when outcomes land at the GFC trough* — the critical mechanism for the QU placebo failure. A shock in 2007Q1 at h=8 has its outcome measured in 2009Q1; states exposed to high-quit industries (leisure, retail) experience large vacancy declines there regardless of the original shock value. The time FE α_t absorbs this at the shock date t but not at the outcome date t+h. Without GFC_{t+h}, the QU instrument's cross-state exposure to consumer-facing industries spuriously predicts vacancy declines at the trough.

β_h in the GFC-controlled spec identifies the IRF in non-GFC shock AND non-GFC outcome quarters — isolating the structural transmission channel from crisis-period dynamics.

### SLOOS interaction LP (credit supply robustness)
```
... + β_h · B_{s,t} + δ_h · B_{s,t} × SLOOS^{dm}_t + ...
```
SLOOS^{dm}_t demeaned within estimation sample. β_h = IRF at average credit conditions. Implemented in `part7b_sloos.py`.

**Key conventions:**
- Base year for employment shares: 2006
- LOO (leave-one-out) national shock rates to avoid mechanical correlation
- Instruments rescaled ×100 for pp units (not ×1000 as in early versions — verify)
- IRF coefficients standardized by cross-sectional SD of instrument (1-SD shock interpretation)
- Outcome window capped at **2019Q4** (COVID exclusion)
- LP horizons: h = 0, 1, ..., 20 quarters
- State FE + time FE; SE clustered at state level (n=50)
- State FIPS as zero-padded 2-digit strings throughout

---

## 7. Data Sources

| Source | Variable | Coverage | Notes |
|--------|----------|----------|-------|
| BLS QCEW (annual bulk CSVs) | Employment by state × industry | 2006 (base year) | agglvl=54, 12 supersectors |
| BLS BED | Establishment closings (δ instrument) | 1992Q3–2024Q2 | LOO national rates |
| BLS JOLTS (national) | Total separations, LD, quits | 2001Q1–2023Q1 | LOO national rates |
| BLS JOLTS (state) | Job openings stock (vacancy outcome) | Dec 2000–present, 50 states | SA series; quarterly = avg of 3 months |
| BLS LAUS | State unemployment rate + labor force | 1990M1–present, 50 states | Monthly → quarterly |
| Chicago Fed NFCI | Financial conditions, risk subindex | 1990Q1+ | Weekly → quarterly; used in part7_nfci.py only |
| BLS OPHNFB | Nonfarm business productivity | Quarterly | Used in part2b productivity residualization |
| FRED RVA series | Quarterly industry real VA by supersector | 2005Q1+ | `RVAM`, `RVAC`, `RVAMA`, `RVAW`, `RVAR`, `RVAT+RVAU`, `RVAI`, `RVAFI+RVARL`, `RVAPBS`, **`RVAESHS`** (not RVAES), **`RVAAER+RVAAF`** (not RVAER) |
| BEA Table 1.3.6 (Excel) | Annual industry real VA | **1997–2025** | `bea_va_annual_by_supersector.csv`; used as annual benchmark constraints in Chow-Lin backcast |
| Atlanta Fed / FRED | Wu-Xia shadow FFR + FEDFUNDS extension | 1960Q1–present | `data/cache/wuxia_quarterly.parquet` |
| Census CBP API | Establishment size by supersector | Year 2000 | `data/cache/cbp_estab_size.parquet`; φ_j for v3 residualization |

---

## 8. File Index (Key Files)

### Annual benchmark (NEW)
- `bea_va_annual_by_supersector.csv` — BEA annual real VA by BLS supersector, 1997–2025, billions of 2017 chained dollars. One column per supersector (BLS10–BLS80). Used by `part2b_extend_va.py` as annual constraint targets for Chow-Lin disaggregation. **File must remain in `Data/Bartek analysis/`**.

### Instruments
- `data/instruments/delta_instrument_resid_base2006.csv` — primary δ instrument (residualized v3)
- `data/instruments/ld_instrument_resid_base2006.csv` — primary LD instrument (residualized v3)
- `data/instruments/qu_instrument_resid_base2006.csv` — quits instrument (residualized v3, placebo)
- `data/instruments/laus_quarterly.parquet` — state outcomes: unemployment + labor force + vacancies

### Caches
- `data/cache/bea_va_quarterly_12ind.parquet` — FRED quarterly VA, 2005Q1+
- `data/cache/bea_va_quarterly_12ind_extended.parquet` — Chow-Lin extended quarterly VA, 1992Q1+; **part2b_shock_comovement.py uses this automatically when present**
- `data/cache/wuxia_quarterly.parquet` — Wu-Xia shadow FFR + FEDFUNDS splice, quarterly; delete to rebuild
- `data/cache/cbp_estab_size.parquet` — log avg establishment size φ_j, 12 supersectors, time-invariant

### LP Results (baseline)
- `data/results/lp_irf_delta_resid.csv` — δ → unemployment (residualized)
- `data/results/lp_irf_ld_resid.csv` — LD → unemployment (residualized)
- `data/results/lp_irf_qu_resid.csv` — QU → unemployment (residualized, placebo)
- `data/results/lp_irf_delta_vacancy.csv` — δ → vacancy rate ⭐
- `data/results/lp_irf_ld_vacancy.csv` — LD → vacancy rate ⭐
- `data/results/lp_irf_qu_vacancy.csv` — QU → vacancy rate (placebo fails) ⭐

### LP Results (GFC-interaction diagnostic, NEW)
- `data/results/lp_irf_delta_gfc.csv` — δ → unemployment (+GFC controls)
- `data/results/lp_irf_ld_gfc.csv` — LD → unemployment (+GFC controls)
- `data/results/lp_irf_qu_gfc.csv` — QU → unemployment (+GFC controls)
- `data/results/lp_irf_delta_vac_gfc.csv` — δ → vacancy (+GFC controls)
- `data/results/lp_irf_ld_vac_gfc.csv` — LD → vacancy (+GFC controls)
- `data/results/lp_irf_qu_vac_gfc.csv` — QU → vacancy (+GFC controls) ⭐ key diagnostic
- `data/results/lp_irf_gfc_qu_comparison.png` — QU vacancy baseline vs GFC-controlled comparison ⭐

### Key Plots
- `data/results/lp_irf_vacancy_decomp_overlay.png` — δ vs LD vs QU vacancy IRFs overlaid
- `data/results/lp_irf_beveridge_asymmetry.png` — δ/LD unemployment and vacancy side-by-side
- `data/results/lp_irf_beveridge_path.png` — (u_h, v_h) trajectory in UV space
- `data/results/lp_irf_delta_ld_qu_overlay.png` — unemployment IRF comparison
- `data/results/va_chowlin_backcast_validation.png` — Chow-Lin sector-by-sector validation plot

---

## 9. Principles for Any Agent Working on This Project

1. **δ vs. s asymmetry is the core.** Every empirical choice should be evaluated against testing this asymmetry.
2. **Residualized instruments only in part5.** Raw instrument LPs have been removed. Do not re-add them.
3. **LD is the primary s-type instrument.** QU is the placebo — it currently fails the placebo test (large delayed vacancy effect). This is a key identification concern, not a resolved finding.
4. **No financial frictions in the model.** SLOOS/NFCI analyses live in `part7b_sloos.py`. Do not add financial interactions to `part5_lp.py`.
5. **LOO always.** National shock rates must leave out the state being instrumented.
6. **COVID cap: 2019Q4** for all LP outcomes.
7. **Single-responsibility pipeline.** Data fetching, instrument construction, and regression are in separate scripts.
8. **FRED series for BLS65 is `RVAESHS`** (Education+Health+Social Assistance aggregate). Not `RVAES` (education alone, ~16% of total). This was corrected April 14.
9. **FRED series for BLS70 is `RVAAER+RVAAF`**. `RVAER` does not exist on FRED. Corrected April 14.
10. **Annual VA benchmarks** (`bea_va_annual_by_supersector.csv`) cover **1997–2025** only. Chow-Lin constrained years: 1997–2004. Pre-1997 backcast (1992–1996) is unconstrained.
11. **Chow-Lin must include intercept** in the GDP regression. Omitting the intercept produces negative centered-R² for sectors where mean growth ≠ mean GDP growth (Mining, Finance, Information). See section 15.1.
12. **r(δ,LD) = 0.424 is structural**, not a contamination artifact. All aggregate confound hypotheses tested across v1→v2→v3 have been rejected. The shared variation reflects common industry-specific fundamentals driving both margins simultaneously. See sections 13.3 and 14.3.
13. **Joint LP is a robustness check for δ, not a structural LD result.** δ vacancy stable across joint/separate. LD zero in joint spec reflects identification failure, not confirmed reposting.
14. **QU placebo failure**: large delayed vacancy effect (h=8–20) likely driven by GFC timing. GFC interaction diagnostic in `part5_lp.py` tests this directly. Results pending.
15. **SLOOS series**: use `DRTSCILM` (C&I tightening standards). NOT `DRSDCILM` (that is loan demand — wrong concept entirely).

---

## 10–14. Prior Session Findings

*[Sections 10–14 unchanged from April 9 version. See prior CLAUDE.md for full text of sections 10 (demand contamination diagnosis), 11 (QU placebo failure), 12 (Granger/VAR calibration), 13 (joint LP), 14 (SLOOS + v3 residualization).]*

---

## 15. New Findings and Discussion — April 14, 2026

### 15.1 Constrained Chow-Lin Backcast of Sector VA Pre-2005

**Problem:** The FRED quarterly RVA series (RVAM, RVAC, etc.) all start 2005Q1, limiting the v2/v3 residualization sample to 2005Q2+ and the LP sample to 2005Q3+. The BED δ instrument covers 1992Q3+, leaving a 12-year window unused.

**Solution:** GDP-anchored constrained Chow-Lin temporal disaggregation (Chow & Lin 1971, JASA), implemented in `part2b_extend_va.py`. Extends quarterly VA back to 1992Q1 for all 12 BLS supersectors. Combined with JOLTS vacancy data starting 2001Q1, the effective LP sample for residualized instruments extends from ~2005Q3 to **2001Q1** — adding four years of pre-GFC identifying variation including the 2001 recession.

**Methodology — Step 1: GDP elasticity regression (post-2005 sample)**

For each BLS supersector j, estimate by OLS using post-2005 quarterly data:

```
Δlog(VA_{j,t}) = α_j + β_j · Δlog(GDP_t) + ε_{j,t}
```

Key: the intercept α_j is **required**. Sectors have sector-specific mean growth that differs from GDP (e.g., Information has secular tech-driven growth above GDP; Mining has secular decline). Omitting the intercept forces the regression through the origin and produces negative centered-R² for these sectors — a mathematical impossibility for OLS with an intercept. The intercept captures trend growth; β_j captures the cyclical GDP elasticity.

R² diagnostic by sector (low R² → β suppressed in backcast):

| BLS | Sector | α | β_OLS | β_used | R² | Ann. constrained |
|-----|--------|---|-------|--------|----|-----------------|
| 10 | Mining | +0.009 | +0.136 | **0.000** | 0.001 | 8 yrs |
| 20 | Construction | −0.005 | +0.848 | +0.848 | 0.287 | 8 yrs |
| 30 | Manufacturing | −0.005 | +1.540 | +1.540 | 0.759 | 8 yrs |
| 41 | Wholesale | −0.004 | +1.334 | +1.334 | 0.580 | 8 yrs |
| 42 | Retail | −0.000 | +1.100 | +1.100 | 0.424 | 8 yrs |
| 43 | Transport+Utils | −0.002 | +1.453 | +1.453 | 0.602 | 8 yrs |
| 50 | Information | +0.012 | +0.678 | +0.678 | 0.248 | 8 yrs |
| 55 | Finance+RE | +0.004 | +0.230 | **0.000** | 0.050 | 8 yrs |
| 60 | Prof+Business | +0.004 | +0.930 | +0.930 | 0.633 | 8 yrs |
| 65 | Education+Health | −0.000 | +1.455 | +1.455 | 0.769 | 8 yrs |
| 70 | Leisure+Hosp | −0.024 | +5.442 | +5.442 | 0.855 | 8 yrs |
| 80 | Other Services | −0.012 | +2.045 | +2.045 | 0.808 | 8 yrs |

R² threshold for β suppression: 0.15. Mining and Finance+RE suppressed to intercept-only backcast — GDP explains near-zero variation in these sectors. Annual constraints do all the meaningful work for them.

**Methodology — Step 2: Unconstrained recursive backcast**

For each pre-2005 quarter, recursively back-propagate from the 2005Q1 anchor level:
```
log(VA_{j,t}) = log(VA_{j,t+1}) - [α_j + β_j · Δlog(GDP_{t→t+1})]
```

**Methodology — Step 3: Annual benchmark constraint (proportional Denton)**

BEA annual VA = mean of 4 quarterly SAAR values. For each year A in 1997–2004 with annual benchmark Y_{j,A}:

```
scale_A = Y_{j,A} / mean(VA^unconstrained_{j,Q1:Q4 of A})
VA^constrained_{j,q} = VA^unconstrained_{j,q} · scale_A   for q ∈ {Q1..Q4 of A}
```

Guarantees mean(VA^constrained) = Y_{j,A} exactly. Years 1992–1996 (no annual benchmark) use unconstrained backcast only. Annual benchmarks sourced from `bea_va_annual_by_supersector.csv` (coverage 1997–2025, see section 8).

**Annual benchmark data — `bea_va_annual_by_supersector.csv`:**
- Source: BEA GDP-by-Industry Table 1.3.6, April 9 2026 release
- Units: billions of chained 2017 dollars (matches FRED quarterly RVA series)
- Coverage: 1997–2025, 12 BLS supersectors
- Aggregation notes:
  - BLS43 = BEA Transport+Warehousing (line 40) + Utilities (line 10) summed
  - BLS55 = BEA Finance+Insurance (line 55) + Real Estate+Rental (line 60) summed
  - BLS65 = BEA line 74 aggregate (matches `RVAESHS` exactly, verified post-2005)
  - BLS70 = BEA line 81 aggregate (matches `RVAAER+RVAAF` exactly, verified post-2005)
- Cross-verification: 11/12 sectors match FRED quarterly mean to ≤0.02% post-2005; all 12 correct

**Important caveat:** The Chow-Lin extension does NOT expand the LP outcome sample window — that remains constrained to 2001Q1 by JOLTS vacancy data. What it expands is the residualization sample for the control variable Δlog VA_{j,t-1}, allowing better residualization of the shock rates for quarters 2001Q1–2004Q4. The impact on IRF results is empirically small (differences at 4th decimal place), confirming the pre-2005 residualization change is a methodological improvement rather than a material change to findings.

---

### 15.2 Updated FRED Series Corrections (BLS_TO_FRED)

Two errors in the BLS_TO_FRED mapping corrected in both `part2b_shock_comovement.py` and `part2b_extend_va.py`:

**BLS65 — Education, Health, and Social Assistance:**
- WRONG (previously): `RVAES` = "Educational Services" subsector alone (BEA line 75), ~16% of full supersector
- CORRECT: `RVAESHS` = "Educational Services, Health Care, and Social Assistance" full aggregate (BEA line 74)
- Impact: the residualization regression for BLS65 was using a series covering only 16% of employment. Now using full aggregate matching BEA annual line 74 exactly (R=1.0000 post-2005).

**BLS70 — Arts, Entertainment, Accommodation, and Food:**
- WRONG (previously): `RVAER` — does not exist on FRED
- CORRECT: `RVAAER` + `RVAAF` summed = Arts+Entertainment + Accommodation+Food (matches BEA line 81 exactly)
- Impact: BLS70 was silently failing or using fallback data in prior sessions. Now correctly specified.

---

### 15.3 QU Placebo Failure — Structural Interpretation

**The finding:** QU vacancy IRF is near-zero h=0–6, then rises monotonically to −0.049 pp at h=15 (larger in magnitude than δ peak of −0.032 pp). Partial-F statistics reach 17–18 at h=14–15. This is not a weak or marginal effect — QU is the strongest instrument for vacancies at long horizons, which is the opposite of what the model predicts.

**Why the delayed onset is the key diagnostic clue:**
The QU instrument loads heavily on high-quit industries: Leisure+Hospitality (β=5.44, R²=0.86), Retail, and Other Services. These industries had the highest quit rates *and* the deepest GFC contractions. The delayed onset of the QU vacancy effect (~7–8 quarters) corresponds exactly to the 2006–2007 shock date cohort whose LP horizon h=8 outcomes land in 2008–2009 — the GFC trough. States with high employment shares in Leisure+Hospitality received large QU Bartik values in 2006–2007 (tight labor market, high quits) and then experienced disproportionately large vacancy declines 8 quarters later as the GFC hit consumer-facing industries hardest.

**Why residualization did not fix it:**
Residualization removes aggregate and industry-specific demand trends from national shock rates (part2b). But the contamination here operates at the *LP identification level*: it is the cross-state covariance between industry employment shares and *subsequent* local economic conditions that generates the spurious correlation. This is a Roth-Sant'Anna (2023) concern about Bartik instruments — the identifying variation comes from employment shares, and if those shares are correlated with future conditions, no upstream residualization can address it.

**Why δ vacancy survives the same concern:**
δ (establishment exits) has an independent structural interpretation — firm exits mechanically destroy vacancy slots. This accounting identity holds regardless of whether demand or credit conditions drive the exit. The δ vacancy effect therefore has causal content even if the instrument is partly demand-contaminated. QU has no analogous independent mechanism: quits do not mechanically destroy vacancies; they trigger reposting. Its vacancy effect must come entirely from the structural transmission channel, so demand contamination fully invalidates it.

**The GFC diagnostic:** `part5_lp.py` now runs GFC-controlled vacancy LPs (section 15.4). If QU vacancy β_h goes flat under GFC controls, the placebo failure is confirmed as GFC-timing contamination. Results pending from next run.

---

### 15.4 GFC Interaction Diagnostic in part5_lp.py

**Motivation:** The QU vacancy placebo failure has a delayed onset (h=8+) that coincides with GFC timing for pre-2008 shock cohorts. To test whether this drives the result, `part5_lp.py` now includes GFC-controlled LP variants alongside baseline.

**Implementation:** `attach_gfc_interaction()` adds `instr_x_gfc` and `gfc_shock` columns to each panel before estimation. Inside `run_lp_horizon()` at each h, `gfc_outcome` is constructed as a dummy for the *outcome quarter* t+h falling in 2008Q3–2009Q4. The two GFC regressors are:

```
B_{s,t} × GFC_t       — differential effect for shocks occurring during GFC
GFC_{t+h}             — aggregate outcome-quarter level shift at GFC trough
```

The outcome-quarter dummy (`GFC_{t+h}`) is the critical addition. Time FEs absorb aggregate variation at the *shock date* t, not the outcome date t+h. For observations with shock in 2006–2007 and h=8–10, the outcome lands in 2008–2009 when all consumer-facing industries contracted. `GFC_{t+h}` absorbs this aggregate collapse, allowing β_h to be identified from non-GFC outcome quarters.

**Output files** (in addition to baseline):
- `lp_irf_{delta,ld,qu}_gfc.csv` — GFC-controlled unemployment IRFs
- `lp_irf_{delta,ld,qu}_vac_gfc.csv` — GFC-controlled vacancy IRFs
- `lp_irf_gfc_qu_comparison.png` — overlay of QU vacancy baseline vs GFC-controlled ⭐

**Removed from part5_lp.py (April 14):**
- All raw instrument LPs (`irf_delta`, `irf_s`, `irf_delta_post`)
- All NFCI interaction LPs from part5 (these live in `part7b_sloos.py`)
- `load_instruments()`, `load_nfci()`, `attach_nfci_interaction()` functions
- `DELTA_INSTR_FILE`, `S_INSTR_FILE`, `NFCI_FILE` path constants
- Dead code: NFCI was being attached to residualized panels but never used in any `run_lp()` call

**Rationale for removing raw instruments:** Raw instrument LPs served as initial diagnostics. The residualized instruments are the paper's identification strategy. Keeping raw LPs adds runtime, output clutter, and dependency on `part7_nfci.py` for main results. NFCI analysis lives in `part7b_sloos.py` where it belongs.

---

### 15.5 Priority Order for Next Session

1. **Run `part5_lp.py` and examine GFC-controlled QU vacancy IRF** — the key pending diagnostic. Compare `lp_irf_qu_vacancy.csv` vs `lp_irf_qu_vac_gfc.csv`. Does QU vacancy go flat? If yes, characterise as GFC contamination; if no, characterise as structural demand contamination pre-dating the GFC.

2. **Wild cluster bootstrap (Priority 4b)** — n=50 clusters at lower bound. Use `wildboottest` Python package to supplement cluster-robust SEs.

3. **Recession-severity placebo (Priority 4a)** — run B×NFCI and B×Δu_nat interactions to formally confirm NFCI proxies recession depth.

4. **Paper text** — empirical section: v1→v2→v3 as systematic confound elimination; δ vacancy is headline result; QU failure characterized via GFC diagnostic; SLOOS confirms frictionless transmission.

---

### 15.6 Updated Data Caches and Files

| File | Source | Coverage | Notes |
|------|--------|----------|-------|
| `bea_va_annual_by_supersector.csv` | BEA Table 1.3.6, April 9 2026 | 1997–2025 | Annual benchmark CSV for Chow-Lin; keep in `Data/Bartek analysis/` |
| `bea_va_annual_supersectors.xlsx` | Same BEA source | 1997–2025 | Formatted Excel with verification sheet; for reference only |
| `data/cache/bea_va_quarterly_12ind_extended.parquet` | part2b_extend_va.py | 1992Q1–2025Q4 | Constrained Chow-Lin backcast; auto-used by part2b_shock_comovement.py when present |
| `data/results/lp_irf_{delta,ld,qu}_gfc.csv` | part5_lp.py (new) | h=0..20 | GFC-controlled unemployment IRFs |
| `data/results/lp_irf_{delta,ld,qu}_vac_gfc.csv` | part5_lp.py (new) | h=0..20 | GFC-controlled vacancy IRFs |
| `data/results/lp_irf_gfc_qu_comparison.png` | part5_lp.py (new) | — | QU vacancy: baseline vs GFC-controlled ⭐ |
| `data/results/va_chowlin_backcast_validation.png` | part2b_extend_va.py | — | Per-sector backcast quality check |
