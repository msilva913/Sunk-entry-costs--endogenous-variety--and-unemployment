# CLAUDE.md — Empirical LP Project (Firm Entry/Exit DSGE)

**Last updated:** May 8, 2026 (session 3)
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
| `part2b_shock_comovement.py` | Regress g^δ and g^s on Δlog(productivity) → residuals ν^δ, ν^s | ✅ Complete | `shock_rates_delta_resid.parquet`, `shock_rates_ld_resid.parquet`, `shock_rates_qu_resid.parquet` |
| `part3_instrument.py` | Bartik aggregation for δ | ✅ Complete | `data/instruments/delta_instrument_base2006.csv` |
| `part3_instrument_s.py` | Bartik aggregation for total separations (s) | ✅ Complete | `data/instruments/s_instrument_base2006.csv` |
| `part3_resid_instruments.py` | Bartik aggregation for residualized ν^δ, ν^{LD}, ν^{quits} | ✅ Complete | `delta_instrument_resid_base2006.csv`, `ld_instrument_resid_base2006.csv`, `qu_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | Fetches LAUS unemployment + JOLTS state vacancies → quarterly panel | ✅ Complete | `data/instruments/laus_quarterly.parquet` (includes `vacancies` column) |
| `part5_lp.py` | All panel LP regressions: unemployment + vacancy outcomes, δ/s/LD/quits | ✅ Complete | IRF CSVs + plots in `data/results/` |
| `part6_shock_persistence.py` | AR(1) estimation of shock persistence | ✅ Complete | `data/results/shock_persistence.csv` |
| `part7_nfci.py` | Chicago Fed NFCI → quarterly, headline + risk subindex | ✅ Complete | `data/cache/nfci_quarterly.parquet` |
| `part2c_granger_lp.py` | Panel LP Granger causality: δ vs LD (raw + residualized) | ✅ Complete | `data/results/granger_lp_*.png` |
| `part2d_var_calibration.py` | Bivariate panel VAR(L) on (ν^δ, ν^{LD}): companion matrix, Σ, Cholesky IRFs | ✅ Complete | `data/results/var_calibration.csv`, `var_irf_chol.png` |
| `part5_joint_lp.py` | Joint LP with δ and LD in same regression (v2-residualized); pre-standardizes instruments; outputs comparison tables and joint vs separate plots | ✅ Complete | `lp_joint_delta_ld_{unemp,vacancy}.csv` + `.png` |
| `part9_recession_scatter.py` | Cross-recession motivating scatter: Bartik cum. exit rate vs. labor market outcomes across 3 NBER recessions | ✅ Complete | `recession_scatter_primary.png`, `recession_scatter_appendix.png`, `recession_scatter_latex.tex` |
| `part10_state_scatter.py` | Cross-state motivating scatter: avg. exit rate vs. CUG/CVS (primary) and recovery speed (appendix) for 50 states | ✅ Complete | `state_scatter_primary.png`, `state_scatter_appendix.png`, `state_scatter.csv`, `state_scatter_latex.tex` |
| `part10b_state_scatter_sensitivity.py` | Sensitivity check: full-sample vs. recession-quarter X-axis for state recovery speed panels | ✅ Complete | `state_scatter_sensitivity.png`, `state_scatter_sensitivity_latex.tex` |
| `part11_jf_table2.py` | BED-based replication and extension of Jaimovich-Floetotto (2008) Table 2: entry/exit share of gross flows and cyclical comovement, 12 Bartik supersectors + total private, 3 sample periods | ✅ Complete | `jf_table2_results.parquet`, `jf_table2_extension.txt`, `jf_table2_{sname}_latex.tex` |

### Notable additions in March 2026 session:
- **Priority 1 completed:** `part4_outcomes.py` now fetches JOLTS state-level job openings (SA, stock → quarterly average). Series ID format verified: `JTS000000{ST}0000000JOL` (21 chars). Data available from Dec 2000 for all 50 states.
- **Priority 2 completed:** JOLTS decomposed into layoffs+discharges (LD) and quits (quits) via `part2_shock_rates_s.py` and `part3_resid_instruments.py`. Separate Bartik instruments built and residualized for each.
- **Vacancy LP outcome** added to `part5_lp.py` via `outcome="vacancy"` parameter. Uses log-change specification; lagged log vacancies replace lagged unemployment rate as control.

### Notable additions in April 8, 2026 session:
- **`part5_joint_lp.py`** written and run: joint LP with δ and LD instruments in the same regression at each horizon h, on v2-residualized instruments. Instruments pre-standardized to 1-SD units before regression so coefficients are directly comparable to separate LPs. Outputs: `lp_joint_delta_ld_unemp.csv`, `lp_joint_delta_ld_vacancy.csv`, comparison plots `lp_joint_delta_ld_{unemp,vacancy}.png` (red = joint, blue = separate).
- **Key finding — joint LP:** LD vacancy coefficient collapses to near-zero and insignificant at all horizons once δ is controlled. δ vacancy coefficient is stable across joint vs separate specs. Replicates colleague's result.
- **Key diagnostic:** δ unemployment coefficient shrinks noticeably in joint vs separate (by ~0.5–1.3 pp across horizons), indicating δ is absorbing shared variation rather than the joint LP cleanly separating two structural effects. The zero LD vacancy result is not straightforwardly interpretable as structural confirmation of reposting — see section 13.
- **Post-v2 instrument correlation confirmed:** r(δ, LD) = 0.423 after v2 enriched residualization — essentially unchanged from v1 (0.389). V2 did not reduce the shared variation, which implies the residual correlation is not driven by industry demand cycles or aggregate productivity.
- **Best path forward discussion completed** — see section 13 for full treatment.

### Notable additions in April 7, 2026 session:
- **`part2c_granger_lp.py`** completed: panel LP Granger causality (δ vs LD, raw + residualized specs). Key finding: δ → LD F-stats = 4.9–8.9 throughout h=0–12; LD → δ fades to 1.1–1.8 by h=8–12 (see section 12).
- **`part2d_var_calibration.py`** completed: bivariate panel VAR(1) on (ν^δ, ν^{LD}) recovering companion matrix, innovation covariance, Cholesky IRFs. Key calibration targets (see section 12): ρ_δ=0.600, ρ_LD=0.510, β(endex)=0.228, corr(u^δ,u^{LD})=0.229.

### Notable additions in April 1, 2026 session:
- **part2b rewritten** to use FRED as data source for industry VA (BEA direct API confirmed blocked; FRED RVAM/RVAC/... series confirmed as only available source, starting 2005Q1).
- **v1 residualization** (fallback, currently active) runs successfully: δ on `(dlog_p,)`, LD and QU on `(dlog_p, log_theta_lag)`.
- **v2 enriched residualization** code is complete in part2b; awaits local run with `FRED_API_KEY` set.
- **`run_v2_locally.py`** written to orchestrate part2b → part3 → part5 in sequence with FRED fetch. Run from `Data/Bartek analysis/` with `$env:FRED_API_KEY = "<key>"` set.
- **CRITICAL NEW FINDING:** QU placebo fails completely. All three instruments produce nearly identical IRFs (see section 11).

---

## 4. Key Empirical Findings

### 4.1 Shock Persistence (part6)

**Raw rates (contaminated by demand/productivity):**
- ρ_δ = 0.617 (half-life 1.44 qtrs, N=117, 1992Q4–2021Q4)
- ρ_LD = 0.489 (half-life 0.97 qtrs, N=88, 2001Q2–2023Q1)

**Residualized rates (preferred calibration targets):**
- ρ_δ = 0.479 (half-life 0.94 qtrs, N=65, 2005Q4–2021Q4)
- ρ_LD = 0.288 (half-life 0.56 qtrs, N=65, 2005Q4–2021Q4)
- Bias from demand/productivity contamination: Δρ_δ = +0.138, Δρ_LD = +0.202

**Innovation correlation (both specs):**
- corr(η^δ, η^{LD}) = +0.397 (raw), +0.444 (residualized)
- Model zero-correlation assumption is violated in both specs
- Motivates allowing corr(ε^δ, ε^s) ≠ 0 akin to Coles-Kelishomi (2018)

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

### 4.4 Residualized Instruments (part2b + part3) — UPDATED April 8

**v1 residualization:**
- r(δ, LD) = **0.389** (down from raw ~0.54)
- r(δ, QU) = **−0.004**
- r(LD, QU) = **−0.492**
- instr_sd: δ = 13.34, LD = 17.02, QU = 8.995

**v2 enriched residualization (run confirmed April 8):**
- r(δ, LD) = **0.423** — essentially unchanged from v1. Adding lagged industry VA growth and lagged tightness did not reduce the shared variation.
- instr_sd: δ = 11.16 pp (×100 units), LD = 18.67 pp
- **Critical implication:** The residual correlation between δ and LD is not driven by industry demand cycles or aggregate productivity. It reflects something that survives both controls — most likely industry-level financial conditions or Schumpeterian reallocation episodes (see section 13).

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

### New tasks arising from April 8 session:

**PRIORITY NEW-A — Decide on primary LP specification (joint vs separate)**
- Joint LP with δ and LD simultaneously is now implemented in `part5_joint_lp.py`
- Separate LP on v2-residualized instruments remains the alternative
- Decision hinges on whether residual r(δ,LD)=0.423 reflects (a) structural endogenous exit feedback (endogenous to δ shock, so separate LP is correct) or (b) a common confound (joint LP addresses OVB)
- See section 13 for the full argument. Current recommendation: present joint LP for δ as robustness check; acknowledge LD is not separately identified in joint spec.

**PRIORITY NEW-B — Industry-level financial conditions as residualization control**
- The leading hypothesis for why r(δ,LD) survives v2 is industry-specific credit supply shocks (consistent with NFCI state-dependence finding)
- Potential fix: add industry-level credit spread or default rate series to part2b residualization
- Sources to explore: BofA/ICE option-adjusted spreads by sector (Bloomberg), or Compustat-based leverage ratios by industry
- This is the cleanest empirical fix and would directly address the instrument non-identification problem
- Implementation: extend `part2b_shock_comovement.py` with industry credit variable; re-run `part3_resid_instruments.py` and `part5_lp.py`

**PRIORITY 3 — SLOOS C&I Index as Alternative Interaction Variable** ❌ Not yet done
- Fetch Senior Loan Officer Opinion Survey C&I net tightening index from FRED
- Re-run δ NFCI interaction replacing NFCI_risk with SLOOS C&I
- Tests whether NFCI result is genuinely credit-supply vs. recession severity
- Implementation: `part7b_sloos.py` or extend `part7_nfci.py`

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

### Joint LP Results (part5_joint_lp.py) ⭐ NEW April 8
- `data/results/lp_joint_delta_ld_unemp.csv` — joint LP: δ and LD coefficients → unemployment, h=0..20
- `data/results/lp_joint_delta_ld_vacancy.csv` — joint LP: δ and LD coefficients → vacancy rate, h=0..20
- `data/results/lp_joint_delta_ld_unemp.png` — joint vs separate comparison plot, unemployment (red=joint, blue=separate)
- `data/results/lp_joint_delta_ld_vacancy.png` — joint vs separate comparison plot, vacancy (red=joint, blue=separate)

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
11. **QU placebo failure is now the leading diagnostic.** Until QU vacancy IRF is flat/insignificant, instrument identification is not established. v2 enriched residualization did not fix this — industry credit conditions are the next candidate control.
12. **Data source for industry VA is FRED** (`FRED_API_KEY` env var). BEA direct API is blocked. Coverage starts 2005Q1 — v2 LP sample will be shorter than v1.
13. **r(δ,LD) = 0.423 survives v2 residualization.** The shared variation is not industry demand or aggregate productivity. Industry-specific financial conditions (credit supply) are the leading candidate. See section 13.3.
14. **Joint LP (part5_joint_lp.py) is a robustness check for δ, not a structural result for LD.** The δ vacancy IRF is stable across joint vs. separate specs. The LD zero in the joint spec reflects identification limits, not confirmed reposting. See section 13.2.

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

---

## 12. New Findings from Session — April 7, 2026

### 12.1 Panel LP Granger Causality (part2c)

Tests whether δ and LD are causally ordered or jointly endogenous, conditional on productivity variation. Two-way within panel (industry + time FEs). HC1 robust SEs. Joint Wald F-test on 4 cross-lags at each horizon h=0..12. Sample: residualized spec 2005Q3–2021Q4.

**Key results (preferred residualized spec, F-statistics):**

| Horizon | δ → LD | LD → δ |
|---------|--------|--------|
| h=0  | 4.89 | 4.68 |
| h=2  | 7.40 | 5.96 |
| h=4  | 8.08 | 3.43 |
| h=8  | 8.89 | 1.82 |
| h=12 | 3.42 | 1.07 |

10% critical value ≈ 2.0–2.5 for joint F(4) test.

**Interpretation:** δ → LD is persistent and significant at ALL horizons through h=12. LD → δ is initially significant (h=0–4) but fades below critical value by h=8–12 and collapses to F=1.07 at h=12. This supports δ as the more primitive shock: establishment exits drive future layoffs (endogenous exit mechanism), but lagged layoffs do not predict future exit rates once shocks have had time to dissipate.

### 12.2 Bivariate Panel VAR Calibration (part2d)

**Specification:** Bivariate panel VAR(L) on (ν^δ_{j,t}, ν^{LD}_{j,t}) with industry + time FEs, two-way within transformation. Cholesky identification with δ ordered first (justified by part2c). Sample: 2005Q3–2021Q4, 12 industries, N=780 obs at L=1.

**Lag selection:** AIC selects L=4; BIC selects L=2. Disagreement → fall back to L=1 for parsimony and interpretability. (AIC/BIC gap is small: all log|Σ| values in range −8.62 to −8.73.)

**Companion matrix A_1 at selected VAR(1):**

|            | lag ν^δ | lag ν^{LD} |
|------------|---------|----------|
| ν^δ eq     | +0.600  | +0.132   |
| ν^{LD} eq  | +0.228  | +0.510   |

**Calibration parameters:**
- ρ_δ = **0.600** (model default from part6 aggregate: 0.617 — close match)
- ρ_LD = **0.510** (model default from part6 aggregate: 0.751 — significantly lower in industry-level within-estimator)
- β (LD←δ lag-1, endogenous exit) = **+0.228**
- α (δ←LD lag-1, reverse) = **+0.132**
- β/α ratio = 1.7× (δ → LD cross-persistence exceeds reverse; VAR evidence directionally consistent with part2c)

**Innovation covariance:**
- corr(u^δ, u^{LD}) = **+0.229** — model zero-correlation assumption is violated
- SD(u^δ) = 0.091, SD(u^{LD}) = 0.150

**Cholesky (δ first):**
- SD(ε^δ structural) = 0.091
- Loading of LD on δ shock (contemporaneous): 0.034
- SD(ε^{LD} orthogonal) = 0.146

**Structural IRFs on shock series:** Both shocks are transitory. ε^δ own-response starts at 0.091, decays to 0.006 by h=8 (half-life ≈ 2 quarters). ε^{LD} own-response starts at 0.146, decays similarly. Cross-response of LD to ε^δ is non-negligible through h=4 (0.034→0.026→0.013), consistent with endogenous exit mechanism.

**Note on ρ_LD discrepancy (0.510 vs 0.751):**
The industry-level panel estimator uses within-industry variation after two-way demeaning, which removes industry and aggregate time effects. The aggregate part6 AR(1) reflects national time-series persistence, which includes aggregate LD cycles. The within-estimator likely understates true persistence because aggregate macro cycles (which generate the most autocorrelation) are absorbed by time FEs. The part6 aggregate estimates remain the preferred calibration target for the model's AR(1) assumption; the panel VAR provides the cross-dynamics β, α, and the innovation correlation.

### 12.3 Priority Order for Next Session

1. **Run v2 locally** (see section 10.5) — still the critical pending step if not yet done
2. **SLOOS interaction** (Priority 3): `part7b_sloos.py`
3. **Recession-severity placebo** (Priority 4a): extend `part5_lp.py`
4. **Wild cluster bootstrap** (Priority 4b): supplement LP SEs
5. **Pre-GFC sample restriction** (Priority 4c): add `max_qt` option to `run_lp()`

---

## 13. New Findings and Discussion — April 8, 2026

### 13.1 Joint LP Implementation and Results (part5_joint_lp.py)

**What was implemented:** `part5_joint_lp.py` runs δ and LD Bartik instruments simultaneously in a single OLS regression at each horizon h, on v2-residualized instruments. Both instruments are pre-standardized to unit SD before the regression so that β^δ and β^LD are directly on the 1-SD scale (matching part5_lp.py's post-multiply scaling). Outputs include comparison tables (joint vs. separate) and plots (red = joint, blue = separate) for both unemployment and vacancy outcomes.

**Key numerical results — unemployment (joint):**
- β^δ: significant throughout h=2–20, peak +2.73 pp at h=13 (vs. +3.29 pp separate — shrinks by ~0.5–1.3 pp)
- β^LD: effectively zero and insignificant at all horizons from h=5 onward; peak +0.75 pp at h=17, p=0.22

**Key numerical results — vacancy (joint):**
- β^δ: significant h=0–8, peaks at −0.53 pp (h=5); fades h=9–13; recovers h=14–19. Pattern robust to controlling for LD.
- β^LD: near-zero and insignificant throughout (h=0: −0.070, p=0.33; h=1: +0.033, p=0.74)

**Replication of colleague's result confirmed.** The LD vacancy effect collapses to zero once δ is controlled. The specific numbers match the pattern reported by the colleague.

### 13.2 Structural Interpretability of the Joint LP

The key diagnostic is that **β^δ shrinks noticeably in the joint spec** (by ~0.5–1.3 pp for unemployment, smaller for vacancy). Under a clean structural separation, β^δ_joint should equal β^δ_separate or be larger (since OVB from positive r(δ,LD) biases β^δ_separate upward if β^LD_h > 0). The shrinkage indicates δ is absorbing the shared variation between the two instruments, not that the joint LP has cleanly resolved two orthogonal structural effects.

**The zero LD coefficient in the joint spec has two competing interpretations:**
1. **(Optimistic — colleague's reading):** The separate LD LP was contaminated by δ variation. Once controlled, the clean LD effect is zero — confirming the model's reposting prediction that layoffs at surviving firms trigger immediate vacancy reposting.
2. **(Pessimistic — preferred reading given diagnostic):** δ absorbs all vacancy-relevant variation because it is the dominant instrument (higher first-stage relevance in the vacancy equation). LD gets residual near-zero signal — not because LD structurally has no effect, but because the joint LP cannot separately identify it given the residual correlation r(δ,LD) = 0.423.

**The δ coefficient stability in the vacancy equation** (joint ≈ separate, differences < 0.13 pp) is the key supporting fact for the optimistic reading for *vacancy*, but less convincing for *unemployment* where δ shrinks more. The honest conclusion is that the joint LP confirms the δ vacancy result is robust, but does not establish that LD truly has zero vacancy effect.

### 13.3 Why r(δ, LD) = 0.423 Survives v2 Enriched Residualization

V2 added lagged industry VA growth and lagged market tightness to the residualization. The instrument correlation is essentially unchanged from v1 (0.389 → 0.423). This is a strongly informative null result: the shared variation between δ and LD is **not** explained by:
- Aggregate productivity cycles (controlled in both v1 and v2)
- Industry-specific demand cycles via value-added growth (controlled in v2)
- Aggregate labor market conditions via tightness (controlled in v2)

**What can drive the residual correlation? Candidates ranked by plausibility:**

1. **Industry-level financial conditions (most likely).** When credit tightens industry-specifically (e.g., construction/real estate post-2007), firms simultaneously exit (raising δ) and lay off workers at surviving firms (raising LD). This is not a product demand shock — it operates through the industry's funding channel. Neither lagged VA growth nor aggregate tightness captures industry-specific credit conditions. Consistent with the NFCI state-dependence finding: the δ effect concentrating in tight financial conditions suggests a credit channel that is industry-specific rather than aggregate.

2. **Schumpeterian reallocation waves.** Rapid technological displacement in a specific industry (e.g., retail/e-commerce, energy) simultaneously raises exit rates and layoff rates at surviving firms that are downsizing to the new equilibrium. This is not a business cycle phenomenon and would survive all current controls.

3. **Firm-specific productivity shocks (less likely as primary driver).** Firm-level heterogeneity aggregates out at the industry×state level that the Bartik construction uses. Firm-specific shocks would only survive aggregation if systematically correlated within industries — making them effectively industry-specific, which collapses to candidates 1 or 2.

4. **Measurement error correlation.** BED (UI records) and JOLTS (establishment survey) are separate surveys but both measure establishment-level outcomes. Shared cyclical mismeasurement in the same direction during industry downturns could generate spurious residual correlation. Harder to rule out but also harder to address.

### 13.4 Alternative Residualization Approach Considered (and Rejected)

**Proposed:** residualize on industry demand only (not productivity), then let productivity drive both δ and LD shocks with its own AR(1) process. Allows a structural decomposition into z-driven and idiosyncratic components.

**Why rejected as primary approach:**
- Aggregate TFP (OPHNFB) is a Solow residual containing cyclical mismeasurement — projecting it out of δ and LD just moves the contamination problem upstream
- Identifying the AR(1) factor process jointly with the loadings φ^δ, φ^LD from two series is barely identified without additional restrictions
- In practice this collapses to v2 residualization: the cleanest implementation is to control for productivity directly in the LP second stage rather than in the instrument construction
- Does not address the r(δ,LD) = 0.423 residual correlation since productivity is not the driver

### 13.5 Best Path Forward

**Two defensible routes:**

**Route A — Joint LP as primary for δ; acknowledge LD is not separately identified.**
- Present `part5_joint_lp.py` results as the preferred δ specification
- δ vacancy IRF is robust to controlling for LD (stable coefficient), establishing it is not an artifact of LD contamination
- Explicitly state that the LD coefficient in the joint spec is not interpretable as the structural LD effect due to identification limitations
- Paper's empirical contribution: industry-level firm destruction shocks (orthogonal to contemporaneous layoffs at continuing firms) generate persistent Beveridge curve dynamics
- Available now; no additional data work required

**Route B — Add industry-level financial conditions control to residualization (preferred if feasible).**
- The leading hypothesis for why r(δ,LD) survives v2 is industry-specific credit supply shocks
- Adding industry credit spread or default rate to part2b residualization would further reduce r(δ,LD) and potentially validate separate LPs
- Possible sources: BofA/ICE option-adjusted spreads by sector (Bloomberg/FRED), Compustat-based leverage ratios by industry
- This is the cleanest empirical fix and would directly address the instrument non-identification problem
- Also connects naturally to the NFCI interaction finding — if industry credit is the confound, controlling for it should simultaneously sharpen the δ–LD asymmetry and explain the financial state-dependence
- Implementation: extend `part2b_shock_comovement.py` with industry credit variable; re-run part3 and part5

**The joint LP and v2 residualization are complementary, not alternatives.** Route B + joint LP is the strongest possible specification. Route A alone is defensible now.

### 13.6 Priority Order for Next Session

1. **Route A immediately actionable:** decide whether to present joint LP as primary or robustness; draft paper text on δ vacancy finding using joint LP results
2. **Route B investigation:** search for accessible industry-level credit spread data (FRED first; Bloomberg if available)
3. **SLOOS interaction** (Priority 3): `part7b_sloos.py` — tests aggregate credit supply vs. recession severity for the NFCI finding; directly relevant to Route B hypothesis
4. **Recession-severity placebo** (Priority 4a): extend `part5_lp.py`
5. **Wild cluster bootstrap** (Priority 4b): supplement LP SEs with wildboottest
6. **Pre-GFC sample restriction** (Priority 4c): add `max_qt` option to `run_lp()`

---

## 14. New Findings and Completed Work — April 17, 2026

### 14.1 Pipeline Fixes Completed This Session

**part2b_residualize_shocks.py — critical panel split fix:**
- Root cause of N=65 clipping: TS (total separations, JOLTS) was being inner-joined with δ inside `_load_panel_delta()`. Fix: separated into two independent panel builders:
  - `_load_panel_delta(min_quarter="1997Q1")` — BED only, no JOLTS; covers 1997Q1+
  - `_load_panel_ldqu()` — TS, LD, QU together; covers 2001Q1+
- TS moved from δ panel to LD/QU panel (it is JOLTS-based, starts 2001Q1)
- Shared controls (productivity, VA lag, tightness lag) loaded once via helper functions
- All matplotlib calls wrapped in `_plot_diagnostics()` with `plt.ioff()` and `_ =` assignments to suppress REPL repr output

**part3_resid_instruments.py:**
- State name NaN bug fixed: `FIPS2D_TO_STATE` keys are zero-padded strings; `state_fips` must be converted to zfill(2) string BEFORE int conversion for the map lookup
- Added quarter range to post-aggregation print; WARNING if δ instrument still starts 2001Q1 (stale parquet signal)
- Correlation matrix now explicitly labels the 2001Q1+ overlap window

**part5_lp.py:**
- Added `delta_resid_1997` panel (1997Q1+ extended) for sample extension comparison
- `irf_delta_1997` computed; `lp_irf_delta_sample_extension.png` overlays 1997Q1+ vs 2001Q1+
- Fixed: duplicate `# ---` separator, redundant import removed, section [10] vac_rate guard added

**part6_shock_persistence.py:**
- Fixed path bug (hyphens → underscores in repo name)
- Replaced `_build_residualized_panel()` (which re-did the join and re-clipped to 2005Q4) with `_load_resid_series()` that reads part2b output parquets directly
- Added 2001Q1+ restriction to both residualized series for common window (pre-2001 expansion years inflate AR(1) mechanically)
- Updated N values: ν^δ residualized N≈87 (was 65), ν^LD N≈87

**Updated persistence estimates (2001Q1+ common window):**
- ρ_δ residualized = 0.647 (unexpectedly above raw 0.617 — see Section 14.2)
- ρ_LD residualized = 0.305
- Innovation correlation residualized = +0.488

### 14.2 Sample Extension Findings (δ Unemployment IRF)

Comparing 1997Q1+ vs 2001Q1+ samples for δ unemployment:
- 2001Q1+: β_h ≈ 0.107 at h=20, N=2800 per horizon (GFC-dominated identification)
- 1997Q1+: β_h ≈ 0.071 at h=20, N=3550 per horizon (attenuation from pre-2001 expansion years)
- Interpretation: pre-2001 expansion (low u, tight labor markets) shows weaker δ→u transmission — workers displaced by firm destruction find jobs faster when u^nat is low. Both results highly significant. The 2001Q1+ estimate is not biased upward by composition — it reflects genuine state-dependence in the transmission channel.

### 14.3 ρ_δ Residualized > ρ_δ Raw: Interpretation

After the panel split and 2001Q1+ restriction, ρ_δ residualized (0.647) exceeds the raw rate (0.617). This is counterintuitive but has a coherent explanation: residualization removes industry-demand variation from the δ rate. If lagged industry VA growth is *less* persistent than the structural firm destruction process, projecting it out increases the residual's persistence. The remaining structural variation in ν^δ (idiosyncratic firm destruction orthogonal to demand cycles) is genuinely more autocorrelated than the raw rate. Both estimates bracket a narrow range (0.617–0.647), so calibration is not sensitive to the choice.

### 14.4 Recession-Severity Placebo Results (part7c) ⭐ Key New Findings

**Script:** `part7c_recession_placebo.py` — estimates augmented LP:
$$y_{s,t+h} - y_{s,t-1} = \alpha_s + \alpha_t + \beta_h B^k_{s,t} + \delta^{lev}_h B^k_{s,t} \times \bar{u}^{nat}_t + \delta^{chg}_h B^k_{s,t} \times \Delta u^{nat}_t + \phi^{lev}_h \bar{u}^{nat}_t + \phi^{chg}_h \Delta u^{nat}_t + \text{controls} + \varepsilon$$

Two severity measures: $\bar{u}^{nat}_t$ = demeaned national unemployment level (primary: regime state); $\Delta u^{nat}_t$ = demeaned change (secondary: momentum). Chosen over NFCI because there is no model in which structural firm destruction effects should depend on national unemployment level — any significant interaction is unambiguously a confound.

**Key results by instrument-outcome pair:**

| Instrument | Outcome | β_h survives? | δ^lev verdict | δ^chg verdict | Conclusion |
|---|---|---|---|---|---|
| δ | unemp | ✅ Yes | Insignificant | Anti-confound at h≥17 | ✅ Clean |
| δ | vacancy | ✅ Yes | **Insignificant ALL h** | **Insignificant ALL h** | ✅ Cleanest result in paper |
| LD | unemp | ✅ Yes | Anti-confound h=2–12 | Momentum confound h≥8 | ✅ Survives |
| LD | vacancy | ⚠️ Marginal | Insignificant | Momentum confound h=0–1, 11–15 | ⚠️ Partial confound |
| QU | unemp | ❌ Delayed onset remains | Mostly insig. | Some negative | ❌ Structural |
| QU | vacancy | ❌ Persists | Level confound h=6–15 | Mostly insig. | ❌ Structural |

**The δ→vacancy finding is the paper's headline:** both severity interactions uniformly insignificant across all 21 horizons. This is the most demanding available placebo and it passes decisively.

**Anti-confound pattern (δ and LD):** δ^chg_h turns *negative* and significant at long horizons for both δ and LD unemployment — meaning high-instrument states suffer *less* additional unemployment when national deterioration is rapid. This is the opposite of a severity confound. For LD, the level interaction is also negative (high-LD states hit less hard when u^nat is elevated). These patterns are consistent with structural channel + state-dependent amplification, not instrument failure.

**QU structural contamination confirmed:** neither severity measure resolves the delayed-onset anomaly or long-horizon vacancy decline. QU contamination is structural (industry-level quit dynamics carry own persistent vacancy implications) and cannot be fixed with macro time-series controls.

**Outputs:**
- CSVs: `data/results/placebo_sev_{delta,ld,qu}_{unemp,vac}.csv`
- Plots: `data/results/placebo_sev_{delta,ld,qu}_{unemp,vac}_overlay.png`, `data/results/placebo_sev_interaction_summary.png`

### 14.5 LD Vacancy Result: Strategic Decision and Economic Interpretation

**Decision:** LD→vacancy results are moved to the appendix. Paper's primary empirical contribution is the δ channel. Aggregate separation rate calibrated as moment (τ = 3.1%) rather than identified from IRF.

**Why LD→vacancy is persistently negative despite residualization:**

JOLTS LD mixes three structurally distinct events:
- **(a) Exogenous match dissolution** — firm survives, reposts vacancy at zero sunk cost. Model's s shock. Predicts no or brief vacancy decline.
- **(b) Deliberate workforce reduction** — firm lays off to desired lower headcount. No reposting. Survives residualization because it is firm-level idiosyncratic, orthogonal to industry VA averages.
- **(c) Partial establishment closure** — plant shuts without full firm exit. Workers laid off, vacancies withdrawn. Concentrated in Manufacturing, Construction. Survives residualization for same reason.

**Why v2 residualization cannot fully fix (b) and (c):** Residualization controls remove the *industry-level average* component of g^LD driven by demand cycles. Cases (b) and (c) operate at the firm level — a corporate restructuring at one large manufacturer elevates the industry LD rate without being predicted by industry-average VA growth (other firms are fine). After residualization, ν^LD still contains these events, and because they involve simultaneous layoffs + vacancy withdrawal, the vacancy decline persists.

The one-quarter lag compounds this: contemporaneous demand signals that trigger preemptive layoffs are invisible to the predetermined lagged control.

**Resolution requires firm-level data** identifying whether each separation was followed by a vacancy posting at the same establishment in the subsequent quarter — not available at required frequency and geography.

### 14.6 Pending Work

Priority order:
1. **Wild cluster bootstrap** (part5_wcrb.py) — n=50 clusters at lower bound; asymptotic SEs may be too narrow at long horizons
2. **GFC diagnostic for LD unemployment** (part7d_ld_gfc.py) — LD unemployment rises monotonically through h=20; apply GFC_{t+h} outcome-quarter control
3. **build_report_html.py** — run locally to generate self-contained HTML report with embedded figures (requires pandoc; `conda install -c conda-forge pandoc`)

### 14.7 Key Files Updated This Session

| File | Change |
|------|--------|
| `part2b_residualize_shocks.py` | Panel split fix; TS moved to LD/QU panel; plt.ioff() |
| `part3_resid_instruments.py` | FIPS state name fix; quarter range diagnostics |
| `part5_lp.py` | 1997Q1+ sample extension; bug fixes |
| `part6_shock_persistence.py` | Path fix; reads part2b parquets directly; 2001Q1+ restriction |
| `part7c_recession_placebo.py` | New script: dual-severity placebo (level + change) |
| `build_report_html.py` | New utility: pandoc-based Markdown→HTML with base64 images |
| `report.md` | Sections 5.6, 5.7 added (severity placebo + LD interpretation); 7.1–7.4 updated |

---

## 15. New Findings and Completed Work — May 5, 2026

### 15.1 Motivating Scatter: Cross-Recession (part9_recession_scatter.py)

**Purpose:** Descriptive figure relating each recession's cumulative permanent-adjusted establishment exit rate (X-axis) to labor market outcomes (Y-axes) across the three NBER recessions (2001, 2007–09, 2020). Motivates the δ channel before the LP analysis.

**Key technical decisions:**
- X-axis: cumulative Bartik delta during recession quarters (raw, non-residualized — goal is illustrative variation, not causal identification). Units corrected: BED closings are in thousands of workers, so multiplied by 1000 before dividing by QCEW employment level.
- V/LF aggregation: recomputed directly from LAUS state panel (sum vac_persons / sum LF × 100) bypassing the tightness cache, which had a mean-vs-sum bug giving log_theta ≈ 3.74 instead of the correct V/LF ≈ 3–4%.
- HWI splice: Barnichon Composite HWI spliced with JOLTS V/LF at 2001Q1 (scale factor ≈ 0.998) for pre-JOLTS vacancy coverage. Path: `../CompositeHWI.csv`.
- Trough depth: arithmetic PP deviations (u_trough − pre_u), not log-deviations.
- Recovery metrics: (a) recession-end anchor — clock starts at NBER recession end; (b) trough anchor — clock starts at series nadir; (c) cumulative gap — area under positive deviation curve from recession START, H=20 quarters, clipped at zero.

**Figure split:**
- **Primary** (`recession_scatter_primary.png`): 2×2 — recession-end recovery speed (row 0) and cumulative gaps CUG/CVS (row 1). No suptitle; captions in LaTeX only.
- **Appendix** (`recession_scatter_appendix.png`): 2×2 — trough depth (row 0) and trough-anchored recovery speed (row 1).
- LaTeX snippets in `data/results/recession_scatter_latex.tex`.

**Key finding:** Unemployment cumulative gap scales monotonically with cumulative exit rate across the three recessions. Vacancy pattern is less clean — 2001 vacancy formation was unusually slow relative to its exit rate (acknowledged in paper text without overdoing it; likely reflects JOLTS series immaturity and structural vacancy posting differences pre-GFC).

---

### 15.2 Motivating Scatter: Cross-State (part10_state_scatter.py)

**Purpose:** 50-state scatter linking structural exit rate exposure (full-sample average Bartik delta) to recession labor market outcomes averaged across three recessions. Complements the cross-recession figure by showing within-US geographic variation.

**X-axis:** Full-sample time average of `bartik_delta` per state (1992Q3–2021Q4). Non-residualized by design — structural cross-state heterogeneity in exit propensity, not causal identification. Confirmed: residualized delta is not appropriate here (geographic variation is the signal).

**Outcome measures** (LF-weighted average across 3 recessions per state):
- **CUG** — cumulative unemployment gap: Σ max(u_{t} − pre_u, 0) over H=20 quarters from recession start (pp-qtrs)
- **CVS** — cumulative vacancy shortfall: Σ max(pre_v − v_{t}, 0) over H=20 quarters from recession start (pp-qtrs)
- **u_rec_end / v_rec_end** — quarters from NBER recession end to 50% recovery of trough gap (appendix)
- **u_rec_trough / v_rec_trough** — quarters from series trough to 50% recovery (appendix)

**Figure split:**
- **Primary** (`state_scatter_primary.png`): 1×2 — CUG (left) and CVS (right). These are the headline cross-state results.
- **Appendix** (`state_scatter_appendix.png`): 2×2 — recovery speed under both anchors (rec-end and trough). Included for completeness but not leading results.
- Outputs also: `state_scatter.csv`, `state_scatter_latex.tex` (both LaTeX snippets).

**Key weighted correlations (r_w, LF weights):**

| Outcome | r_w | Notes |
|---------|-----|-------|
| CUG | +0.647 | Headline result — correct sign, strong |
| CVS | +0.282 | Correct sign, moderate |
| u_rec_end | −0.267 | Wrong sign — see below |
| v_rec_end | +0.075 | Near zero |
| u_rec_trough | −0.439 | Still wrong sign |
| v_rec_trough | +0.059 | Near zero |

**Why recovery speed has wrong/null sign — diagnosis:**
Recovery speed is confounded by state economic structure. High-delta states are predominantly service/Sun Belt economies (CA, FL, NV, HI, TX, NJ, NY) that experienced the 2001 recession mildly and recover quickly on average. Low-delta states are Manufacturing/Rust Belt (IN, OH, KY, ME, AR, IA) that recover slowly regardless of recession. The cross-state variation in exit rates is driven by industry composition, not treatment intensity — switching anchors (trough vs. rec-end) does not fix this. CV of avg_delta_rate = 4.8%; unweighted r for u_rec_end = +0.011, confirming the negative weighted result is driven by large states, not a general pattern.

**Resolution:** CUG is the right metric because it integrates depth × duration, combining both channels through which exit shocks generate persistent damage. Recovery speed conflates structural recovery capacity with the exit shock effect. Paper leads with CUG/CVS as the cross-state motivating evidence; recovery speed panels appear only in the appendix with a brief caveat about industry composition.

**Paper narrative for fig:state_scatter:** States with structurally higher establishment exit rates accumulate significantly larger total unemployment gaps during recessions (r_w = 0.65), consistent with the model's prediction that firm destruction shocks generate persistent rather than transitory Beveridge curve dislocations. The vacancy shortfall shows the same directional pattern (r_w = 0.28) but more weakly, consistent with a partial reposting channel.

---

### 15.3 Sensitivity Check: X-Axis Definition (part10b_state_scatter_sensitivity.py)

**Purpose:** Test whether the cross-state recovery speed correlations depend on how the X-axis establishment exit rate is defined — full-sample average vs. recession-quarter average (13 quarters spanning the three NBER episodes).

**Results:**

| Outcome | Full-sample X r_w | Recession-quarter X r_w | Full-sample slope | Rec-quarter slope |
|---------|------------------|------------------------|-------------------|-------------------|
| U recovery | −0.267 | −0.284 | −1323 | −1080 |
| V recovery | +0.075 | +0.081 | +133 | +118 |

**Key facts:** Cross-state correlation between the two X measures = 0.995; mean ratio = 1.183 (recession quarters have ~18% higher exit rates uniformly — pure level rescaling, no reranking). Weighted correlations barely move and slopes rescale approximately by 1/1.18 as expected. The null/wrong-sign result for recovery speed is equally present under both X definitions, confirming the issue is economic structure, not measurement choice.

---

### 15.4 Key Files Added This Session

| File | Function |
|------|----------|
| `part9_recession_scatter.py` | Cross-recession scatter: cum. exit rate vs. outcomes (3 recessions) |
| `part10_state_scatter.py` | Cross-state scatter: avg. exit rate vs. CUG/CVS (primary) + recovery speed (appendix) |
| `part10b_state_scatter_sensitivity.py` | Sensitivity: full-sample vs. recession-quarter X for recovery speed panels |

| Output | Description |
|--------|-------------|
| `data/results/recession_scatter_primary.png` | Cross-recession primary figure (rec-end recovery + CUG/CVS) |
| `data/results/recession_scatter_appendix.png` | Cross-recession appendix figure (depth + trough recovery) |
| `data/results/recession_scatter_latex.tex` | LaTeX snippets for both recession scatter figures |
| `data/results/state_scatter_primary.png` | Cross-state primary figure (CUG + CVS) — fig:state_scatter |
| `data/results/state_scatter_appendix.png` | Cross-state appendix figure (recovery speed, 2 anchors) |
| `data/results/state_scatter.csv` | State-level outcomes: avg_delta_rate, CUG, CVS, recovery speeds |
| `data/results/state_scatter_latex.tex` | LaTeX snippets for both state scatter figures |
| `data/results/state_scatter_sensitivity.png` | Sensitivity figure: full-sample vs. recession-quarter X |
| `data/results/state_scatter_sensitivity_latex.tex` | LaTeX snippet for sensitivity figure |

### 15.5 Pending Work (Updated)

Priority order:
1. **Wild cluster bootstrap** (part5_wcrb.py) — n=50 clusters at lower bound; asymptotic SEs may be too narrow at long horizons
2. **GFC diagnostic for LD unemployment** (part7d_ld_gfc.py) — LD unemployment rises monotonically through h=20; apply GFC_{t+h} outcome-quarter control
3. **build_report_html.py** — run locally to generate HTML report (requires pandoc)
4. **SLOOS interaction** (Priority 3) — `part7b_sloos.py`: C&I net tightening as alternative to NFCI
5. **Recession-severity placebo** (Priority 4a) — extend `part5_lp.py` with Δu_nat interaction
6. **Pre-GFC sample restriction** (Priority 4c) — add `max_qt` option to `run_lp()`

---

## 16. New Work — May 5, 2026 (Session 2)

### 16.1 JF Table 2 Replication and Extension (part11_jf_table2.py)

**Purpose:** Replicate and extend Table 2 of Jaimovich & Floetotto (2008, JME) — "Firm entry, labor market dynamics, and the business cycle" — using the project's BED data. Provides empirical support for calibrating the relative importance of establishment entry and exit in aggregate gross job flows.

**Script:** `part11_jf_table2.py` in `Data/Bartek analysis/`

**Data source:** `data/cache/bed_gross_flows_correct.parquet`
- Column naming caveat (critical): `'expanding'` column = elem0001 = **total** gross gains (not just expanding establishments); `'openings'` column = elem0002 = gains at **expanding (continuing)** establishments only
- Derived series: G_O (births) = `expanding − openings`; L_C (closings/deaths) = `closings` (elem0006)
- Total private = sum of 12 Bartik supersectors; industry code `000000` for BED aggregate

**Industries:** 12 Bartik supersectors + TOTAL (aggregate). Ordered: TOTAL, 10 (Mining), 20 (Construction), 30 (Manufacturing), 41 (Wholesale), 42 (Retail), 43 (T&U), 50 (Information), 55 (Finance/RE), 60 (PBS), 65 (Edu/Health), 70 (Arts/Food), 80 (Other services).

**Sample periods:**
1. `jf_orig`: 1992Q3–2006Q3 (replicates JF original)
2. `ext_2019`: 1992Q3–2019Q4 (COVID cutoff, primary extension)
3. `ext_2024`: 1992Q3–2024Q2 (full BED coverage)

**Table columns:**

| Col | Description | Method |
|-----|-------------|--------|
| C1 | Entry share of gross gains: Σ(G_O)/Σ(G) | Raw level means |
| C2 | Exit share of gross losses: Σ(L_C)/Σ(L) | Raw level means |
| C3 | Relative std dev: sd(HP log G_O) / sd(HP log G) | HP filter λ=1600; values >1 = openings amplify cycle |
| C4 | Relative std dev: sd(HP log L_C) / sd(HP log L) | HP filter λ=1600; values >1 = closings amplify cycle |
| C3H | Hamilton relative std dev: sd(Ham log G_O) / sd(Ham log G) | Hamilton (2018) filter; robustness for C3 |
| C4H | Hamilton relative std dev: sd(Ham log L_C) / sd(Ham log L) | Hamilton (2018) filter; robustness for C4 |
| C5 | u-projected fitted var ratio (HP): Var(fitted HP G_O) / Var(fitted HP G) | HP-filter log flows and log u; project on [1, cu_t, cu_{t-1}]; ratio of fitted variances |
| C6 | u-projected fitted var ratio (HP): Var(fitted HP L_C) / Var(fitted HP L) | Same for losses |
| C5H | Same as C5 via Hamilton filter | Hamilton filter throughout |
| C6H | Same as C6 via Hamilton filter | Hamilton filter throughout |

**Key methodological notes:**
- All four cyclical columns (3-6) use HP-detrended series (λ=1600), following JF's "detrended series" language. Hamilton filter backup table produced separately.
- Cols 3-4: relative standard deviation (not R²) — standard RBC-style moment. HP applied to **raw levels** (not logs). Log-transforming before HP inflates the ratio mechanically (G_O ~ 20% of G, so log G_O has higher HP-cycle amplitude than log G even when level cycles are proportional). Raw-level HP gives C3-HP = 0.369, C4-HP = 0.336 for total private, matching JF's ~0.33.
- Cols 5-6: unemployment used as projector instead of JF's GDP; both flows and unemployment HP-filtered before projection.

**Key findings (TOTAL row, all samples):**

| Sample | C1 | C2 | C3-HP | C4-HP | C5-HP | C6-HP |
|--------|----|----|-------|-------|-------|-------|
| jf_orig (1992Q3–2006Q3) | 0.213 | 0.206 | 0.369 | 0.336 | 0.007 | 0.022 |
| ext_2019 (preferred, –2019Q4) | 0.201 | 0.195 | 0.282 | 0.224 | 0.009 | 0.013 |
| ext_2021 (incl. COVID, –2021Q4) | 0.200 | 0.194 | 0.207 | 0.133 | 0.016 | 0.015 |

- **C3-HP/C4-HP match JF (~0.33)** in jf_orig sample — validates methodology
- **Secular decline in C1/C2**: entry/exit shares fall ~1 pp post-GFC (business dynamism decline)
- **C3/C4 fall in extended samples**: GFC intensive margin (mass layoffs at continuing firms) dominates total flow variance, making the extensive margin look relatively less volatile — compositional, not a calibration problem
- **COVID distortion (ext_2021)**: the 2020Q1–Q2 spike in gross flows warps the HP trend across the whole sample; ext_2019 is the preferred extended sample for structural comparisons
- **C5-HP/C6-HP very small (~0.01) throughout**: entry/exit margin accounts for little of the unemployment-driven component of total flows — intensive margin dominates. Supports calibrating δ/τ ≈ 10% as a minority channel.
- **BED cache ends 2021Q4**: ext_2021 is the true BED endpoint; no data beyond that in cache

**Output files:** `jf_table2_{jf_orig,ext_2019,ext_2021}_latex.tex` (HP main); `_ham_latex.tex` (Hamilton backup)

**Output files in `data/results/`:**
- `jf_table2_results.parquet` — full panel (all samples × industries × statistics)
- `jf_table2_extension.txt` — human-readable ASCII table with documentation header
- `jf_table2_jf_orig_latex.tex` — LaTeX tabular for 1992Q3–2006Q3 sample
- `jf_table2_ext_2019_latex.tex` — LaTeX tabular for 1992Q3–2019Q4 sample
- `jf_table2_ext_2024_latex.tex` — LaTeX tabular for 1992Q3–2024Q2 sample

**LaTeX console print:** Script prints all three LaTeX tables to stdout with `%===` separators for direct copy-paste into .tex document.

**Technical issues resolved this session:**
- `bed_gross_flows_correct.parquet` column names are inverted vs. economic meaning (see data source caveat above) — correct mapping confirmed by summing goods-producing supersectors against BED aggregate
- UnicodeEncodeError in Spyder (cp1252 default): fixed with `encoding="utf-8"` on all `open()` calls + `.encode("ascii", errors="replace")` guard in `write_text_table`
- scipy not installed in sandbox: `pip install scipy --break-system-packages`

### 16.2 δ/τ Calibration (Completed)

- **δ_bar = 0.940%/qtr**: BDS employment-weighted exit rate (annual, ÷4), common sample 1978–2019
- **τ_bar = 9.34%/qtr**: Shimer (2012) total separation rate, common sample
- **δ/τ = 10.07%**: product-line destruction accounts for ~10% of all separations at steady state
- These are the model calibration targets; the BED-based JF Table 2 provides supporting micro-evidence on the quantitative role of entry/exit in aggregate flows

---

## 17. Revisions and Findings — May 8, 2026 (Session 3)

### 17.1 part11_jf_table2.py — Major Methodology Revision

**Column structure redesigned** (breaking change vs. prior session):

| Col | Old spec | New spec |
|-----|----------|----------|
| 3 | HP rel. std dev | Hamilton (h=8, p=4) rel. std dev — MAIN |
| 4 | HP rel. std dev | Hamilton (h=8, p=4) rel. std dev — MAIN |
| 5 | R² (filtered log u on filtered log G_O) | HP (λ=1600) rel. std dev — JF comparison |
| 6 | R² (filtered log u on filtered log L_C) | HP (λ=1600) rel. std dev — JF comparison |

**Why the R² spec was dropped (cols 5-6):** Unemployment u_t is highly autocorrelated, so any bivariate regression of filtered(log u) on filtered(log G_O) picks up the predetermined lagged component of u mechanically regardless of the actual entry/exit signal. The R² values were economically uninterpretable and sample-unstable (C5 flipped from 0.017 to 0.166 across jf_orig vs ext_2019). The parallel structure (Hamilton cols 3-4, HP cols 5-6, same statistic both) is cleaner and directly answers the filter-sensitivity question.

**Why HP rather than Hamilton for the JF comparison (cols 5-6):** HP in the jf_orig sample gives C5=0.362, C6=0.333, closely matching JF's 0.33/0.34. Hamilton gives slightly higher values (0.391/0.335) because it is less aggressive at trend removal over short windows. HP remains the appropriate filter for direct JF replication; Hamilton is preferred for extended samples containing GFC/Covid.

**Quast-Wolters modification evaluated and rejected:** QW averages Hamilton residuals over h=4..12 to reduce horizon sensitivity. Computed and compared — values differ by ≤0.015 from Hamilton h=8 across all samples. QW does not reduce the secular decline in C3/C4 across extended samples because the Covid/GFC spikes inflate sd(L) in the denominator regardless of h. Not worth adding as a separate column.

**Table format change:** Removed `tablenotes` environment entirely. Notes moved into a compact `\multicolumn{7}{p{...}}` row between the last data row and `\bottomrule`, inside the tabular. Caption is now factual description only (no findings text). Three separate `.tex` files output per sample; no HP-backup appendix file.

**Output files (current, data/results/):**
- `jf_table2_jf_orig_latex.tex` — 1992Q3–2006Q3 (JF replication)
- `jf_table2_ext_2019_latex.tex` — 1992Q3–2019Q4 (preferred extended)
- `jf_table2_ext_2024_latex.tex` — 1992Q3–2024Q4 (BED end; truncates at 2021Q4 until cache refreshed locally)
- `jf_table2_results.parquet` — columns: sample, pip, c1, c2, c3 (Ham), c4 (Ham), c5 (HP), c6 (HP), n
- `jf_table2_extension.txt` — human-readable table

**BED column naming clarification (confirmed):** `df["G_O"] = df["expanding"] - df["openings"]` is correct. The legacy cache column `"openings"` contains BED elem0002 = gains at *expanding (continuing)* establishments, not births. `"expanding"` contains elem0001 = total gross gains. The subtraction yields gains at *opening (birth)* establishments. Confirmed via `frac_open` column in cache (~0.19–0.22, consistent with C1 ≈ 20%).

### 17.2 Key Numbers for part11 (ext_2019 preferred sample)

**Total private:**
- C1 = 0.200, C2 = 0.195 (entry/exit mean shares)
- C3 = 0.269, C4 = 0.239 (Hamilton rel. std dev)
- C5 = 0.275, C6 = 0.223 (HP rel. std dev)

**Industry range (ext_2019):**
- Manufacturing: C1=12.5%, C2=13.9% (lowest entry/exit shares)
- Financial Activities: C1=23.5%, C2=24.7%; Leisure & Hospitality: C1=24.5%

**jf_orig total private (replication check):**
- C5(HP)=0.362, C6(HP)=0.333 vs. JF benchmark 0.33/0.34 — gap is ~3 pp from supersector aggregation

**Secular decline diagnosis (confirmed):** The drop in C3/C4/C5/C6 from jf_orig to ext_2019/ext_2024 is entirely driven by the denominator: sd(L) jumps from 254 (jf_orig) → 322 (ext_2019) → 1208 (ext_2024), while sd(LC) barely moves (85 → 72 → 159). The Covid spike creates max|cL| = 12,168 vs. 738 in the original sample — a 16× outlier. Hamilton is more robust than HP to this but neither filter recovers the original ratios. This is genuine economics (intensive margin dominates GFC/Covid variance), not a methodological artifact.

### 17.3 Empirical Motivation — Four Facts Structure (Paper Draft)

The empirical motivation section is organized around four interconnected facts, presented in plain "First / Second / Third / Fourth" structure (not bolded headers — bolded headers read like slides, not a journal paper):

**Fact 1 — Job flows track establishment entry/exit closely (BDS):**
- Establishment exit rate avg 10.4%/yr, firm exit 8.9%/yr, job destruction 14.1%/yr
- Establishment–firm exit correlation = 0.943

**Fact 2 — Entry/exit account for ~20% of gross flows and are cyclically relevant (BED/JF Table 2):**
- Mean shares: 20.0% of gains, 19.5% of losses (ext_2019)
- Range: 12–14% (Manufacturing) to 23–25% (Financial Activities, Leisure & Hospitality)
- Cyclical volatility ratios: 0.27/0.24 (Hamilton), 0.28/0.22 (HP) — below JF due to GFC/Covid intensive-margin dominance, but economically significant
- Key calibration contrast: monthly separation rate 3.1% >> monthly establishment exit rate 0.31% (= δ̄ 0.940%/qtr ÷ 3). Gap motivates the model: not all separations destroy a business opportunity; not all firm exits are preceded by worker separations.

**CAUTION on monthly exit rate figure in draft:** Text originally said 0.651%/month — does not match any series. Correct value is δ̄ = 0.940%/qtr = 0.313%/month (from BDS). Text should use 0.313%/month or equivalently 0.940%/qtr.

**Fact 3 — Cumulative exit predicts recession severity across recessions (recession scatter):**
- 2001: cum_delta=5.6%, CUG=30.0 pp-qtrs, CVS=25.7 pp-qtrs, u_peak=5.9%
- 2008-09: cum_delta=7.8%, CUG=76.6 pp-qtrs, CVS=18.8 pp-qtrs, u_peak=9.9%
- 2020: cum_delta=2.9%, CUG=26.6 pp-qtrs, CVS=1.7 pp-qtrs, u_peak=13.0%
- Note: 2001 CVS (25.7) > 2008-09 CVS (18.8) despite lower exit rate — prolonged dot-com vacancy freeze, not just exit
- 2020 CUG (26.6) is substantial despite tiny CVS — demand recovered fast but unemployment remained elevated briefly

**Fact 4 — High-exit states accumulate larger labor market gaps (state scatter):**
- Weighted correlations: r_w(CUG) = 0.647, r_w(CVS) = 0.282
- Recovery speed has wrong/null sign due to industry-composition confound (high-exit states are service/Sun Belt, fast structural recovery); CUG is the right metric

### 17.4 Pending Work (Updated Priority Order)

1. **Wild cluster bootstrap** (`part5_wcrb.py`) — n=50 clusters at lower bound; asymptotic SEs may be too narrow at long horizons ❌
2. **GFC diagnostic for LD unemployment** (`part7d_ld_gfc.py`) — LD unemployment rises monotonically through h=20; apply GFC_{t+h} outcome-quarter control ❌
3. **SLOOS interaction** (`part7b_sloos.py`) — C&I net tightening as alternative to NFCI; tests credit supply vs. recession severity ❌
4. **Pre-GFC sample restriction** — add `max_qt="2007Q4"` option to `run_lp()` in `part5_lp.py` ❌
5. **Refresh BED cache to 2024Q4** — run `refresh_bed_cache.py` locally (blocked in sandbox by BLS API rate limit); ext_2024 currently truncates at 2021Q4 ❌
6. **build_report_html.py** — run locally to generate HTML report (requires pandoc) ❌