# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)

**Last updated:** April 1, 2026
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

### Notable additions in April 1, 2026 session:
- **part2b rewritten** to use FRED as data source for industry VA (BEA direct API confirmed blocked; FRED RVAM/RVAC/... series confirmed as only available source, starting 2005Q1).
- **v1 residualization** (fallback, currently active) runs successfully: δ on `(dlog_p,)`, LD and QU on `(dlog_p, log_theta_lag)`.
- **v2 enriched residualization** code is complete in part2b; awaits local run with `FRED_API_KEY` set.
- **`run_v2_locally.py`** written to orchestrate part2b → part3 → part5 in sequence with FRED fetch. Run from `Data/Bartek analysis/` with `$env:FRED_API_KEY = "<key>"` set.
- **CRITICAL NEW FINDING:** QU placebo fails completely. All three instruments produce nearly identical IRFs (see section 11).

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

### 4.4 Residualized Instruments (part2b + part3) — UPDATED April 1

**v1 residualization (currently active):**
- r(δ, LD) = **0.389** (down from raw ~0.54)
- r(δ, QU) = **−0.004**
- r(LD, QU) = **−0.492**
- instr_sd: δ = 13.34, LD = 17.02, QU = 8.995

Note: r(δ,LD) = 0.389 after v1 residualization; target is to drive this further down with v2 (adding lagged industry VA growth).

### 4.5 JOLTS Decomposition — LD vs. Quits (Updated April 1 — PLACEBO FAILS)

**⚠️ CRITICAL: QU placebo fails completely.** All three instruments produce nearly identical unemployment and vacancy IRFs (see section 11 for details). This is a smoking gun for demand contamination, not structural δ/LD/QU differences.

Three unemployment LP specifications (v1 residualized):
1. **δ (residualized):** Peak β ≈ +1.57 pp at h=19, significant throughout h=0–20
2. **LD (residualized):** Peak β ≈ +1.68 pp at h=17, significant throughout h=0–20
3. **QU (residualized, placebo):** Peak β ≈ +1.32 pp at h=15, **significant from h=4 onward** — PLACEBO FAILURE

LD is primary s-type instrument. QU should be placebo (flat/zero IRF) but is not — see section 11.

### 4.6 Vacancy IRFs (v1 residualized) — UPDATED April 1

**δ shock → vacancies (`lp_irf_delta_vacancy.csv`, instr_sd=13.34):**
- Vacancies fall persistently: h=1: −0.121, h=4: −0.253, h=8: −0.305, h=12: −0.409, h=14: −0.426 (peak)
- Significant from h=1 onward (p<0.01 from h=2)

**LD shock → vacancies (`lp_irf_ld_vacancy.csv`, instr_sd=17.02):**
- Peak β = −0.445 at h=20; significant through h=0–11, fades h=12–15, returns h=16–20
- Pattern similar to δ, not the flat/rising profile the model predicts

**QU shock → vacancies (`lp_irf_qu_vacancy.csv`, instr_sd=8.995, PLACEBO — FAILS):**
- Peak β = −0.440 at h=15, **significant from h=7–8 onward, rising monotonically** through h=15
- **QU placebo fails for vacancies too** — near-identical magnitude to δ and LD

**Smoking gun:** corr(QU unemployment betas h=0..20, QU vacancy betas h=0..20) = −0.971. Same for δ: −0.948. Both instruments move unemployment up and vacancies down with nearly perfect mechanical proportionality — characteristic of demand contamination, not structural identification.

**Key plot:** `data/results/lp_irf_vacancy_delta_ld.png` — core asymmetry test (δ vs. LD vacancy IRFs)

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
2. **Vacancy LP is now available.** Use it. But as of April 1, QU placebo fails — vacancy IRFs not yet separately identified across instruments.
3. **LD is the primary s-type instrument** (not total separations). QU is placebo — but QU currently fails the placebo test (see section 11).
4. **No financial frictions in the model.** Do not introduce them beyond the NFCI interaction term.
5. **LOO always.** National shock rates must leave out the state being instrumented.
6. **COVID cap: 2019Q4** for all LP outcomes. Instruments can extend further.
7. **Single-responsibility pipeline.** Data fetching, instrument construction, and regression are in separate scripts.
8. **State FIPS as zero-padded 2-digit strings** throughout.
9. **Vacancy aggregation = average** (stock measure), not sum. Separation aggregation = sum (flow measure).
10. **The δ financial state-dependence and the s monotonic rise are puzzles**, not confirmations. Flag them in the paper text.
11. **QU placebo failure is now the leading diagnostic.** Until QU vacancy IRF is flat/insignificant, instrument identification is not established. v2 enriched residualization is the proposed fix.
12. **Data source for industry VA is FRED** (`FRED_API_KEY` env var). BEA direct API is blocked. Coverage starts 2005Q1 — v2 LP sample will be shorter than v1.

---

## 10. New Insights from Session — March 31, 2026

### 10.1 Demand Contamination Diagnosis (CRITICAL)

The most important new finding is that both δ and LD instruments are likely contaminated by **industry-specific demand shocks** that survive productivity residualization. This is now the leading explanation for why both instruments show declining vacancies and similar unemployment IRFs, contrary to the model's predicted asymmetry.

**The mechanism:**
When a demand contraction hits a specific industry nationally (e.g., manufacturing), three things happen simultaneously: establishment closings rise (raising g^δ), layoffs at surviving firms rise (raising g^LD), and firms in that industry reduce vacancy posting. States with high exposure to that industry receive large Bartik values for both instruments in the same quarter when their vacancies are falling for demand reasons.

**Two clarifications on this channel:**
1. Contamination does NOT require the shock to be state-concentrated — it only needs to be industry-specific. The industry dimension alone suffices because the Bartik construction weights states by industry composition.
2. The contamination applies symmetrically to both δ and LD because they share identical employment weights and both underlying series respond to industry demand conditions. Aggregate productivity residualization does not remove industry-specific demand because it only controls for the aggregate nonfarm component.

This explains why δ and LD remain correlated at r=0.44 after productivity residualization — the shared industry-specific demand component survives.

**Implication:** The similar IRFs for δ and LD likely reflect both instruments predominantly identifying the effect of industry-specific demand contractions on state labor markets, which directly suppress vacancy creation and raise unemployment regardless of the structural δ vs. LD distinction.

---

### 10.2 Decision: Drop TS as Primary Instrument

**Total separations (TS) should be removed from primary LP analysis** and retained only as:
- A robustness/comparability specification (to connect to earlier results)
- An internal consistency check (TS ≈ LD + QU + other in the pipeline)
- A bridge to pre-2001 data if needed later

**Reason:** TS = LD + QU + other separations. Since LD and QU are nearly orthogonal after residualization (r = −0.035) and have structurally distinct implications for the vacancy channel, running TS averages two mechanically opposite signals. The LD and QU instruments separately dominate TS.

**Primary instrument set going forward:**
- δ (residualized) — product-line destruction
- LD (residualized) — layoffs and discharges at continuing establishments
- QU (residualized) — quits (placebo; no reposting prediction in model)

---

### 10.3 Enriched Residualization Specification (NEXT TASK)

The current residualization in part2b uses only aggregate productivity growth:
```
log g^k_{j,t} = α_j + γ_k Δlog p_t + ν^k_{j,t}
```

**The new specification adds lagged industry value-added growth and lagged market tightness:**
```
log g^k_{j,t} = α_j + γ_k Δlog p_t + λ_k Δlog VA_{j,t-1} + μ_k log θ_{t-1} + ν^k_{j,t}
```

where:
- `Δlog VA_{j,t-1}` = real quarterly value-added growth in BEA industry j, **lagged one quarter** (predetermined w.r.t. current shock)
- `log θ_{t-1}` = national market tightness (v/u), **lagged one quarter**
- Both regressors lagged to ensure predetermination — contemporaneous values are endogenous to the shock

**Why each control:**
- `Δlog p_t`: aggregate macro cycle (common to all industries)
- `Δlog VA_{j,t-1}`: industry-specific demand conditions — the KEY new addition to purge demand contamination
- `log θ_{t-1}`: aggregate labor market conditions affecting exit decisions and layoff propensity; especially important for QU (quits are primarily driven by tightness)

**Why NOT use per capita GDP:** Industry GDP has no natural population denominator; total real chained-dollar value-added growth is the right measure.

---

### 10.4 Industry VA Data Source: FRED (BEA Direct API Blocked)

**BEA direct API (`apps.bea.gov`) is blocked** — confirmed both in the Cowork sandbox AND on the user's local machine. Switched to FRED as data source.

**FRED series IDs for BLS supersectors (BLS_TO_FRED mapping in part2b):**
| BLS Code | FRED Series IDs | Description |
|---|---|---|
| 10 | RVAM | Mining |
| 20 | RVAC | Construction |
| 30 | RVAMA | Manufacturing |
| 41 | RVAW | Wholesale trade |
| 42 | RVAR | Retail trade |
| 43 | RVAT, RVAU | Transport+Warehousing + Utilities (sum) |
| 50 | RVAI | Information |
| 55 | RVAFI, RVARL | Finance+Insurance + Real estate (sum) |
| 60 | RVAPBS | Professional & business services |
| 65 | RVAES, RVAHC | Education + Health & social assistance (sum) |
| 70 | RVAER, RVAAF | Entertainment/recreation + Accommodation/food (sum) |
| 80 | RVAOSEG | Other services |

**Critical limitation:** All FRED BEA Real VA series start **2005Q1** — no quarterly industry VA available before 2005 from any accessible source. This limits v2 residualization sample to 2005Q2 onward (one lag needed).

**v2 sample implication:** v2 residuals cover 2005Q2–2023Q1 (JOLTS end). LP sample will be 2005Q2–2019Q4 for v2 specifications. Shorter than v1 (2001Q3–2019Q4) but still covers GFC and post-GFC period.

---

### 10.5 Current Pipeline Status (April 1, 2026)

**What is complete and working:**
- `part2b_shock_comovement.py` — rewritten to use FRED fetch (`_fetch_fred_va()`). v1 fallback (no FRED key) and v2 path (with FRED key) both coded. File uses `FRED_API_KEY` env var (not BEA key).
- `part3_resid_instruments.py` — runs on existing v1 parquets, produces valid instrument CSVs.
- `part5_lp.py` — runs successfully, produces all 6 key IRF CSVs.
- `run_v2_locally.py` — orchestrator script written; runs part2b → part3 → part5 with FRED fetch.

**What requires a local run to complete:**
- v2 enriched residualization: requires `FRED_API_KEY` env var set locally.
- Run from `Data/Bartek analysis/` in PowerShell:
  ```
  $env:FRED_API_KEY = "<your_key>"
  python run_v2_locally.py
  ```
- Free FRED API key: https://fred.stlouisfed.org/docs/api/api_key.html
- After first run, `data/cache/bea_va_quarterly_12ind.parquet` is cached; FRED key not needed again.

**Key tests once v2 runs:**
1. Does r(δ,LD) drop below 0.39? (demand contamination confirmation)
2. Is QU vacancy IRF now flat/insignificant? (instrument non-identification resolved)
3. Do δ and LD vacancy IRFs diverge? (δ down, LD flat/up = model prediction confirmed)

---

### 10.6 Interpretation Framework for Paper

The emerging interpretation:

1. **Endogenous exit feedback is NOT the main explanation for similar LD IRFs.** Endogenous exit requires surviving firms to have depressed continuation values — but LD firms have already survived. The timing is wrong for endogenous exit to dominate reposting.

2. **Demand contamination IS the leading explanation.** Both instruments pick up industry-specific demand contractions. The enriched residualization (adding lagged industry VA growth) is designed to test and address this.

3. **The δ–LD vacancy asymmetry IS visible but partial.** δ: persistent significant decline; LD: moderate decline, fading after h=8. This partial asymmetry is consistent with reposting being real but weaker than the model implies, possibly because firms in recession-adjacent industries do not fully repost even when they survive.

4. **If enriched residualization sharpens the δ–LD vacancy asymmetry** (δ down, LD flat/up), this is the paper's key empirical finding: the δ/s decomposition matters, and the demand channel was obscuring it.

5. **If enriched residualization does NOT sharpen the asymmetry**, the reposting channel is genuinely weak quantitatively at empirically relevant parameter values — also a publishable finding that disciplines the model's calibration.

---

### 10.7 Vacancy IRF Interpretation (Updated)

The declining vacancy IRF after LD (despite model predicting flat/rising) was initially attributed to:
- Instrument collinearity (r=0.44 with δ) — now seen as INSUFFICIENT at this correlation level
- Endogenous exit feedback — now seen as IMPLAUSIBLE for LD timing reasons

**The new leading explanation:** demand contamination in the LD instrument causes states hit by industry demand shocks to show both high LD and falling vacancies, even though the structural reposting channel would predict rising vacancies. The enriched residualization test will resolve this.

---

## 11. New Findings from Session — April 1, 2026

### 11.1 Instrument Non-Identification (CRITICAL NEW FINDING)

**The three instruments (δ, LD, QU) produce nearly identical IRFs in both unemployment and vacancy dimensions.** This is the most important empirical finding of the April 1 session.

**Unemployment IRFs (v1 residualized, peak β per 1-SD shock):**
- δ: peak +1.57 pp at h=19 (significant throughout)
- LD: peak +1.68 pp at h=17 (significant throughout)
- QU: peak +1.32 pp at h=15 (**significant from h=4 — PLACEBO FAILS**)

**Vacancy IRFs (v1 residualized, peak β per 1-SD shock):**
- δ: peak −0.426 pp at h=14 (significant h=1–20)
- LD: peak −0.445 pp at h=20 (significant, fades then returns)
- QU: peak −0.440 pp at h=15 (**significant from h=7–8 — PLACEBO FAILS**)

The QU instrument was intended as a structural placebo — quits are voluntary separations that should trigger immediate reposting. Both the theoretical prior and the model predict a flat or positive vacancy response to QU. Empirically QU produces a vacancy IRF indistinguishable from δ. This is not consistent with QU identifying a distinct structural margin.

### 11.2 Smoking Gun for Demand Contamination

**Within-instrument correlation between unemployment and vacancy IRF coefficients across horizons:**
- QU: corr(β_unemp_h, β_vac_h) for h=0..20 = **−0.971**
- δ: corr(β_unemp_h, β_vac_h) for h=0..20 = **−0.948**

When unemployment goes up by x pp at horizon h, vacancies go down by approximately proportional amount — with near-perfect mechanical consistency. This is the signature of instruments that primarily identify **demand-side labor market contractions** (which simultaneously raise unemployment and lower vacancies through the Beveridge curve) rather than distinct structural shocks with different vacancy implications.

This pattern is inconsistent with structural identification of δ vs. s vs. quits. It is consistent with all three instruments being correlated with industry-specific demand conditions in the overlapping employment-share weighted states.

### 11.3 Updated Instrument Correlations (v1 Spec)

| Pair | Correlation |
|------|-------------|
| r(δ, LD) | 0.389 |
| r(δ, QU) | −0.004 |
| r(LD, QU) | −0.492 |

The near-zero r(δ,QU) and strongly negative r(LD,QU) suggest partial separation by construction, but the QU IRF non-flatness shows correlation ≠ identification.

### 11.4 Revised Interpretation of Prior "Asymmetry" Finding

The earlier finding that "δ vacancy IRF is more persistent than LD" (section 4.6, 10.7) likely reflects a difference in **demand shock exposure** rather than structural reposting behavior:
- δ instrument loads more heavily on industries prone to long-lived demand contractions (durable manufacturing, construction)
- LD instrument has broader sectoral coverage, diluting the persistent demand effects

The v2 enriched residualization (lagged industry VA growth) is specifically designed to remove this industry-specific demand component and test whether any instrument-specific structural signal survives.

### 11.5 Data Source Note (part2b)

The `part2b_shock_comovement.py` script now uses:
- `FRED_API_KEY` (not `BEA_API_KEY`) for industry VA fetch
- Hardcoded `BLS_TO_FRED` dictionary with series IDs listed in section 10.4
- `fredapi` Python package (`pip install fredapi`)
- Cache path: `data/cache/bea_va_quarterly_12ind.parquet`
- v1 fallback activates automatically if FRED fetch fails or key is absent

### 11.6 Priority Order for Next Session

1. **Run v2 locally** (see section 10.5) — this is the critical pending step
2. If v2 sharpens asymmetry: draft paper text on demand contamination and cleaned identification
3. If v2 does NOT sharpen asymmetry: consider pre-GFC sample restriction (Priority 4c) and wild bootstrap (Priority 4b) as next diagnostic steps
4. SLOOS interaction (Priority 3) and recession-severity placebo (Priority 4a) remain pending regardless of v2 outcome
