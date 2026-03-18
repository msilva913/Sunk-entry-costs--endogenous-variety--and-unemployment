# Empirical LP Exercise — Project Brief

**Date:** March 18, 2026
**Author:** Mario (compiled from collaborative sessions)
**Purpose:** Handoff document for a task-oriented agent to execute remaining empirical work on the local projections exercise for the firm entry/exit DSGE paper.

---

## 1. Project Context

This is the empirical component of a DSGE paper extending the Gabrovski-Silva framework (itself building on Bilbiie, Ghironi, and Melitz 2012) with aggregate shocks to three variables: δ (firm/product line destruction), s (match separation), and z (technology). The model adds full business formation and exit to a search-and-matching labor market with finitely elastic vacancy creation (Coles and Kelishomi 2018).

**The model's central prediction:** δ shocks destroy both product lines and vacancies simultaneously (vacancies are tied to product lines), generating persistent Beveridge curve dynamics. s shocks destroy only matches — surviving firms repost vacancies at zero additional sunk cost, so the vacancy stock partially recovers. These two shocks therefore have qualitatively different transmission through the vacancy law of motion.

**The LP exercise's role:** Provide external empirical discipline on the δ vs. s transmission asymmetry before structural estimation. The LP uses Bartik (shift-share) instruments to identify shock-specific impulse responses of state-level unemployment to δ and s shocks.

---

## 2. Existing Pipeline — What Has Been Built and Run

All scripts live in a single directory. The pipeline runs sequentially:

| Script | Function | Key Output |
|--------|----------|------------|
| `construct_delta_instrument.py` | Library: QCEW shares, BED shock rates, Bartik aggregation for δ | (imported by runners) |
| `construct_s_instrument.py` | Library: JOLTS shock rates, Bartik aggregation for s | (imported by runners) |
| `part1_shares.py` | Compute base-year (2006) state×supersector employment shares ω_{s,j} from QCEW | `shares_base2006.parquet` |
| `part2_shock_rates.py` | Fetch BED establishment closing rates → LOO national industry shock rates g^δ_{-s,j,t} | `shock_rates_1992Q3_2024Q2.parquet` |
| `part2_shock_rates_s.py` | Fetch JOLTS total separation rates → LOO national industry shock rates g^s_{-s,j,t} | `s_shock_rates_2001Q1_2023Q1.parquet` |
| `part2b_shock_comovement.py` | Regress g^δ and g^s on Δlog(productivity) → residuals ν^δ, ν^s; compute cross-series correlations | `nu_delta.parquet`, `nu_s.parquet` |
| `part3_instrument.py` | Bartik aggregation for δ: B^δ_{s,t} = Σ_j ω_{s,j} · g^δ_{-s,j,t} | `delta_instrument_base2006.csv` |
| `part3_instrument_s.py` | Bartik aggregation for s | `s_instrument_base2006.csv` |
| `part3_resid_instruments.py` | Bartik aggregation using residualized ν^δ, ν^s from part2b | `delta_instrument_resid_base2006.csv`, `s_instrument_resid_base2006.csv` |
| `part4_outcomes.py` | Fetch state-level LAUS unemployment rates and labor force (1990M1–present, 50 states) | `laus_quarterly.parquet` |
| `part5_lp.py` | Panel local projections: baseline, post-2001 δ, NFCI interaction, residualized IRFs | CSVs + plots in `data/results/` |
| `part6_persistence.py` | AR(1) estimation of shock persistence: ρ_δ ≈ 0.617, ρ_s ≈ 0.751 | persistence results |
| `part7_nfci.py` | Fetch Chicago Fed NFCI (weekly → quarterly, 1990Q1+), save headline + risk subindex | `nfci_quarterly.parquet` |

**Data sources used:** BLS QCEW (annual bulk CSVs), BLS BED (establishment closings), BLS JOLTS (separations), BLS LAUS (state unemployment), Chicago Fed NFCI.

**Base year for shares:** 2006.
**12 BLS supersectors** mapped to pipeline industry codes (from QCEW agglvl=54).
**Instrument samples:** δ starts 1992Q3 (BED availability); s starts 2001Q1 (JOLTS availability). Both end ~2023Q1.
**Outcome window:** Capped at 2019Q4 to exclude COVID contamination.
**LP horizons:** h = 0, 1, 2, ..., 16 quarters.

---

## 3. LP Specification (Current Baseline)

For shock k ∈ {δ, s}, horizon h:

```
y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h · B^(k)_{s,t} + γ₁ · u_{s,t-1} + γ₂ · log(LF_{s,t-1}) + ε_{s,t,h}
```

- y_{s,t} = state unemployment rate (%)
- α_s = state fixed effects
- α_t = time fixed effects (absorb national cycle; identification from cross-state dispersion)
- B^(k)_{s,t} = Bartik instrument with LOO national shock rates and predetermined base-year shares
- u_{s,t-1} = lagged unemployment rate (predetermined control)
- log(LF_{s,t-1}) = lagged log labor force (scale control)
- Standard errors clustered at state level (n=50)

**NFCI interaction variant:**
```
... + β_h · B_{s,t} + δ_h · B_{s,t} × NFCI_risk_dm_t + ...
```
NFCI risk subindex demeaned within estimation sample so β_h = IRF at average financial conditions.

**IRF standardization:** All coefficients multiplied by the cross-sectional SD of the instrument, so raw and residualized IRFs are on the same scale (response to 1-SD instrument shock).

---

## 4. Key Empirical Findings So Far

### 4.1 δ shock (baseline)
- Significant unemployment response that plateaus around h=5–9, broadly consistent with model predictions and ρ_δ = 0.617 half-life.
- Post-2001 subsample produces similar IRF shape (sample comparability confirmed).

### 4.2 s shock (baseline)
- **Anomaly:** IRF rises monotonically through h=16 even though the shock itself dissipates by h=8–10 (ρ_s = 0.751 → shock essentially zero by h=10). This is mathematically inconsistent with linear transmission.

### 4.3 NFCI interaction results
- **δ shock:** β_h (at average financial conditions) is insignificant at all horizons beyond h=0. δ_h (interaction) is highly significant (p<0.001) throughout. **Interpretation:** The entire δ unemployment effect operates through tight financial conditions. At average conditions, δ shocks have no detectable unemployment effect. This is a substantive finding — but it is a puzzle relative to the model, which has no financial friction.
- **s shock:** β_h grows steadily and becomes significant from h=7 onward, even at average financial conditions. δ_h is also significant but smaller. The monotonically rising s IRF **survives** the NFCI interaction — financial amplification does not explain the anomaly.

### 4.4 Residualized instruments (after purging productivity comovement)
- After residualizing g^δ and g^s on Δlog(productivity), the two IRFs are on similar orders of magnitude. The productivity component was inflating raw δ estimates.
- Cross-instrument correlation (raw ≈ 0.54) is reduced by residualization but remains meaningfully positive — δ and s instruments still share substantial variation.

### 4.5 Central tension
The model predicts s shocks should partially self-correct via reposting while δ generates persistent Beveridge curve dynamics. Empirically, δ effects concentrate almost entirely in tight financial conditions (not in the model), while s generates a monotonically rising IRF even at average financial conditions (inconsistent with shock persistence). This tension needs to be addressed directly in the paper.

---

## 5. Remaining Tasks — Prioritized

### PRIORITY 1: Add Vacancies as LP Outcome Variable

**Why:** This is the most direct test of the model's core asymmetry. The model predicts:
- δ shock → vacancy stock falls (product line destroyed, vacancy destroyed)
- s shock → vacancy stock rises or is flat (surviving firms repost)

Without this result, the paper cannot claim empirical support for the mechanism in the vacancy law of motion.

**Implementation:**
- Use state-level JOLTS job openings data (available from Dec 2000, all states from 2010; earlier state estimates may exist but check coverage)
- Add to part4_outcomes.py or a new part4b script
- Run existing LP specification with Δ(vacancy rate) or Δlog(vacancies) as outcome variable y
- Produce IRF plots for both δ and s shocks with vacancy outcome alongside unemployment outcome
- The s sample (post-2001) aligns with JOLTS availability

**Expected diagnostic value:** If vacancies fall after both shocks, the reposting channel is empirically weak. If vacancies diverge (fall after δ, flat/rise after s), the model's mechanism is validated.

### PRIORITY 2: JOLTS Decomposition — Layoffs+Discharges vs. Quits

**Why:** The current s Bartik instrument conflates three types of separations with fundamentally different economic content:
1. **Quits** — procyclical, no model counterpart (no on-the-job search), unclear reposting incentive. Including them contaminates the instrument with procyclical variation.
2. **Temporary layoffs** — recall arrangements bypass the reposting channel entirely. Worst-case contamination for testing the model.
3. **Permanent discharges at continuing establishments** — the cleanest analog to the model's s_t. Involuntary, match truly dissolved, firm faces genuine reposting decision.

This is a first-order identification concern, not a minor robustness issue.

**Implementation:**
- Modify `construct_s_instrument.py` / create variant to pull JOLTS **layoffs-and-discharges** and **quits** separately by industry (series available from 2001)
- Construct separate Bartik instruments: B^{LD}_{s,t} and B^{quits}_{s,t}
- Run three LP specifications:
  1. Baseline s (total separations, for comparability)
  2. Layoffs-and-discharges only (primary — closest to model's s_t)
  3. Quits only (placebo — model predicts different IRF profile)
- Compute cross-instrument correlations: corr(B^δ, B^{LD}), corr(B^δ, B^{quits}), corr(B^{LD}, B^{quits})
- Expected: B^{LD} and B^δ moderately correlated (both countercyclical); B^{quits} and B^δ low or negative correlation (opposite cyclicality). This would solve the instrument collinearity problem.
- Also run productivity residualization on each component separately: expect γ_{quits} > 0, γ_{LD} < 0 (opposite productivity loadings)

### PRIORITY 3: SLOOS C&I Index as Alternative Interaction Variable

**Why:** The NFCI risk subindex is a composite that may proxy recession severity rather than credit supply. The Federal Reserve's Senior Loan Officer Opinion Survey (SLOOS) C&I net tightening index is a direct measure of bank credit supply to businesses. It provides a more interpretable robustness check on whether the NFCI interaction is genuinely capturing credit supply constraints.

**Implementation:**
- Fetch SLOOS C&I net tightening index from FRED (quarterly, available from ~1990)
- Create part7b_sloos.py or extend part7 to save alongside NFCI
- Re-run interaction LP replacing NFCI_risk with SLOOS C&I
- Compare δ_h significance across the two measures
- SLOOS varies only over time (like NFCI), so it does not solve the industry-level variation problem, but provides cleaner structural interpretation

### PRIORITY 4: Additional Robustness Checks

**4a. Recession-severity placebo for interaction term:**
Run three parallel interaction versions:
1. B_{s,t} × NFCI_t (current)
2. B_{s,t} × Δu^{nat}_t (national unemployment change — pure recession severity)
3. Both simultaneously

If δ_h from version 1 is absorbed by version 2, NFCI is proxying recession depth, not financial amplification. If δ_h survives version 3, there is independent financial content.

**4b. Wild cluster bootstrap p-values:**
Current state-clustered SEs are at the lower bound of reliability with n=50 clusters. Supplement with wild cluster bootstrap.

**4c. Pre-GFC sample restriction for s:**
Restrict s sample to pre-GFC quarters only to test whether the monotonically rising s IRF is driven entirely by GFC observations.

---

## 6. Items Explicitly Deferred

- **Industry-level financial conditions residualization:** Theoretically appealing but empirically infeasible. Compustat coverage is adequate for Manufacturing, Information, and Financial Activities, but poor to unusable for Construction, Leisure & Hospitality, Education & Health Services, and Other Services — precisely the sectors most financially constrained and most relevant to GFC results. Flag in paper as ideal robustness check that data limitations prevent implementing cleanly.
- **Financial frictions in the model:** The model itself has no financial friction. Financial conditions enter only through the empirical LP interaction design. Model extension is a separate project.
- **BDS annual exit rates as alternative δ instrument:** BDS classifies exit only if zero employment for the entire following year (permanent by construction), which filters out temporary COVID closures. But annual frequency limits within-recession dynamics. Useful as robustness but not primary specification.

---

## 7. File Locations and Data

**Working directory:** All scripts in one directory.
**Data subdirectories:**
- `data/instruments/` — shares, shock rates, instrument CSVs (raw and residualized)
- `data/outcomes/` — LAUS parquet
- `data/nfci/` — NFCI quarterly parquet
- `data/results/` — LP IRF CSVs and plots

**Model draft:** `Draft_21Feb2026.tex`
**Steady state code:** `steady_state_refactored.jl`

**Key parameters and conventions:**
- Base year: 2006
- Monthly model frequency (Shimer calibration); instruments at quarterly frequency; conversion to monthly in calibration step
- Bartik instruments rescaled by ×1000 to get correct quarterly rate units
- 50 states (DC excluded from LAUS in current pipeline)
- State FIPS as zero-padded 2-digit strings throughout
- NFCI risk subindex demeaned within estimation sample before constructing interaction

---

## 8. Key Principles for the Agent

1. **The δ vs. s asymmetry is the model's central prediction.** Every empirical choice should be evaluated against whether it helps or hinders testing this asymmetry through the vacancy law of motion.
2. **The vacancy LP outcome is the single most important missing test.** Without showing vacancies respond differently to δ vs. s, the paper's core mechanism is untested.
3. **The s instrument composition problem is first-order.** The JOLTS decomposition into layoffs+discharges vs. quits is the most important robustness check for identification, not just for completeness.
4. **The δ financial state-dependence finding and the s monotonic rise are puzzles, not confirmations.** They require direct engagement in the paper text.
5. **No financial frictions in the model.** Do not introduce them in the empirical specification beyond the interaction term.
6. **Pipeline principle: each file has a single responsibility.** Data fetching, instrument construction, and regression analysis are separated. New data sources get their own part files.
7. **Always use LOO (leave-one-out) construction for national shock rates** to avoid mechanical correlation between instrument and local conditions.
8. **Outcome window capped at 2019Q4** to exclude COVID contamination in LP outcomes (instruments can extend further for construction purposes).
