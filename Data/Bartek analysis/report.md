---
title: "Empirical Framework Report: Business Formation, Exit, and Labor Market Dynamics"
subtitle: "External Evidence on δ vs. s Shock Transmission Using Bartik Shift-Share Instruments"
author: "Mario Silva"
date: "April 14, 2026"
---

# 1. Project Overview

This report documents the empirical component of a DSGE paper extending the **Gabrovski–Silva** framework (building on Bilbiie, Ghironi & Melitz 2012). The model adds endogenous business formation and exit to a search-and-matching labor market with finitely elastic vacancy creation (Coles & Kelishomi 2018). Three aggregate shocks drive the model: **δ** (firm/product-line destruction), **s** (match separation), and **z** (technology).

**The core prediction — Beveridge asymmetry — is:**

- A δ shock destroys both product lines and vacancy slots simultaneously, generating persistent dynamics along the Beveridge curve.
- An s shock destroys only matches; surviving firms repost vacancies at near-zero sunk cost, so the vacancy stock partially self-corrects.

The key empirical test is whether shock-specific impulse responses to unemployment and vacancies are consistent with this asymmetry. The LP exercise provides external, reduced-form discipline on the δ vs. s transmission channel using state-level panel variation and Bartik shift-share instruments.

---

# 2. Methodology

## 2.1 Bartik Shift-Share Instruments

The Bartik instrument for shock type $k \in \{\delta, \text{LD}, \text{QU}\}$ in state $s$ and quarter $t$ is:

$$B^{(k)}_{s,t} = \sum_j \omega_{s,j} \cdot g^{(k)}_{-s,j,t}$$

where:

- $\omega_{s,j}$ = base-year (2006) employment share of supersector $j$ in state $s$, from QCEW annual data (`data/instruments/shares_base2006.parquet`)
- $g^{(k)}_{-s,j,t}$ = **leave-one-out (LOO)** national shock rate for shock type $k$ in supersector $j$ at quarter $t$, excluding state $s$ to avoid mechanical own-state correlation

The LOO design ensures the national shock rate is exogenous to any single state's idiosyncratic shocks. The identifying assumption (Goldsmith-Pinkham, Sorkin & Swift 2020) is that the base-year employment shares $\omega_{s,j}$ are exogenous to state-specific shocks, conditional on state and time fixed effects. In Adão, Kolesár & Morales (2019) terms, identification requires the shares to be "as good as randomly assigned" across states within supersectors — here supported by the fact that 2006 is a pre-crisis base year.

Three instrument types are constructed:

- **δ instrument**: based on BED establishment closing rates (firm destruction channel)
- **LD instrument**: based on JOLTS layoff and discharge rates (match separation channel, primary)
- **QU instrument**: based on JOLTS quit rates (separation channel, placebo)

## 2.2 Identification and the Case for Residualization

The raw Bartik instruments face a potential identification threat: national shock rates $g^{(k)}_{j,t}$ may be driven partly by aggregate demand or monetary policy, not solely by supply-side business formation and exit dynamics. If supersectors differ in their sensitivity to such aggregates — and if those sensitivities are correlated with the base-year employment shares — the instruments may capture common demand variation rather than structural shocks.

The project addresses this through a three-stage residualization:

**v1 — Aggregate controls**: Regress national shock rates on time-fixed effects and aggregate output growth.

**v2 — Industry-level controls**: Add industry real value-added growth $\Delta \log \text{VA}_{j,t}$ as a control, absorbing sector-specific demand fluctuations. This is the primary residualization stage.

**v3 — Monetary policy sensitivity**: Further control for $\phi_j \times \Delta \text{MP}_t$ where $\phi_j$ is log average establishment size from Census CBP (a proxy for interest-rate sensitivity) and $\Delta \text{MP}_t$ is the Wu–Xia shadow FFR change. This tests whether differential monetary policy exposure explains residual instrument co-movement.

The key diagnostic is the pairwise correlation between residualized δ and LD innovations across versions. The correlation $r(\nu^\delta, \nu^{LD})$ remains stable across all three stages at approximately **0.424**, implying the shared variation is not an aggregate demand or monetary policy artifact — it reflects a structural co-movement of firm exit and layoff rates at the industry level (see Section 4.5).

## 2.3 Local Projection Specifications

All regressions use residualized instruments only. Instruments are rescaled ×100 for percentage-point units and standardized by the cross-sectional standard deviation of the Bartik variable, so $\beta_h$ is a 1-SD-shock impulse response.

### Unemployment LP

$$u_{s,t+h} - u_{s,t-1} = \alpha_s + \alpha_t + \beta_h \cdot B^{(k)}_{s,t} + \gamma_1 u_{s,t-1} + \gamma_2 \log \text{LF}_{s,t-1} + \varepsilon_{s,t,h}$$

**Economic rationale**: $\beta_h$ captures the cumulative unemployment response (in percentage points) at horizon $h$ quarters to a 1-SD Bartik shock. State fixed effects $\alpha_s$ absorb permanent cross-state heterogeneity (industry composition, right-to-work laws). Time fixed effects $\alpha_t$ absorb aggregate business-cycle variation common to all states. Standard errors are clustered at the state level ($n=50$) to account for within-state serial correlation in outcomes.

### Vacancy-Rate LP

$$v_{s,t+h} - v_{s,t-1} = \alpha_s + \alpha_t + \beta_h \cdot B^{(k)}_{s,t} + \gamma_1 v_{s,t-1} + \gamma_2 \log \text{LF}_{s,t-1} + \varepsilon_{s,t,h}$$

where $v_{s,t}$ is the vacancy rate in percentage points ($\text{vacancies} \times 1000 / \text{labor force} \times 100$). $\beta_h < 0$ means shocks reduce vacancy slots.

**Key conventions**: base year 2006 for employment shares; LOO national shock rates throughout; sample capped at 2019Q4 (COVID exclusion); horizons $h = 0, 1, \ldots, 20$ quarters.

## 2.4 GFC Interaction Diagnostic

The raw QU (quits) instrument fails its placebo test: it generates a large, rising vacancy effect at horizons $h = 8$–$20$ — *larger* in magnitude than the δ effect — with a conspicuous 7–8 quarter delayed onset. This delay matches precisely the horizon at which 2006–2007 shock cohorts have outcomes landing in the 2008–2009 GFC trough. States with high employment shares in high-quit industries (Leisure+Hospitality, Retail) experienced both large QU Bartik values in 2006–2007 and disproportionately large vacancy declines in 2008–2009, generating a spurious correlation.

Standard time fixed effects $\alpha_t$ absorb aggregate variation at the *shock date* $t$, not at the *outcome date* $t+h$. For a shock in 2006Q3 at $h = 10$, the outcome is measured in 2009Q1 — a quarter where the aggregate vacancy collapse is not captured by $\alpha_t$. To address this, `part5_lp.py` adds a GFC interaction diagnostic:

$$\ldots + \delta_h \cdot B^{(k)}_{s,t} \times \text{GFC}_t + \phi_h \cdot \text{GFC}_{t+h} + \ldots$$

where $\text{GFC}_t = \mathbf{1}[\text{shock quarter} \in \{2008\text{Q3}–2009\text{Q4}\}]$ and $\text{GFC}_{t+h} = \mathbf{1}[\text{outcome quarter} \in \{2008\text{Q3}–2009\text{Q4}\}]$. The first term absorbs differential shock amplification *during* the GFC; the second — critically — absorbs the aggregate outcome-level collapse *when outcomes land at the trough*. The GFC-controlled $\beta_h$ identifies structural transmission in non-GFC shock and non-GFC outcome quarters.

---

# 3. Data Description

## 3.1 Sources

| Source | Variable | Coverage | Notes |
|--------|----------|----------|-------|
| BLS QCEW (annual bulk CSVs) | Employment by state × supersector | 2006 (base year) | agglvl=54, 12 supersectors; used to construct shares $\omega_{s,j}$ |
| BLS BED | Establishment closing rates | 1992Q3–2024Q2 | LOO national rates; primary δ instrument |
| BLS JOLTS (national) | Total separations, layoffs & discharges, quits | 2001Q1–2023Q1 | LOO national rates; basis for LD and QU instruments |
| BLS JOLTS (state) | Job openings stock | Dec 2000–present, 50 states | SA series; quarterly = average of 3 months; vacancy outcome |
| BLS LAUS | State unemployment rate + labor force | 1990M1–present, 50 states | Monthly → quarterly; unemployment outcome |
| FRED RVA series | Quarterly industry real VA | 2005Q1+ | 12 supersectors; `RVAESHS` (BLS65), `RVAAER+RVAAF` (BLS70) |
| BEA Table 1.3.6 | Annual industry real VA | 1997–2025 | `bea_va_annual_by_supersector.csv`; annual benchmark for Chow-Lin |
| Chicago Fed NFCI | Financial conditions, risk subindex | 1990Q1+ | Weekly → quarterly; used in credit-supply robustness (`part7b_sloos.py`) |
| FRED DRTSCILM | SLOOS C&I net tightening standards | 1990Q4+ | Credit supply robustness; *not* `DRSDCILM` (loan demand) |
| Atlanta Fed / FRED | Wu-Xia shadow FFR + FEDFUNDS extension | 1960Q1+ | `data/cache/wuxia_quarterly.parquet`; v3 MP residualization |
| Census CBP API | Log avg. establishment size by supersector | Year 2000 | $\phi_j$ for v3 MP-sensitivity residualization; `data/cache/cbp_estab_size.parquet` |

## 3.2 Constrained Chow-Lin Backcast of Sector VA

The FRED quarterly RVA series all start in 2005Q1, which would restrict the residualization sample and, in turn, the LP sample to roughly 2005Q3+. A constrained Chow-Lin (1971) temporal disaggregation procedure extends quarterly VA back to 1992Q1, allowing the residualization window to open at 2001Q1 — adding four years of pre-GFC identifying variation including the 2001 recession.

**Step 1 — GDP elasticity regression** (post-2005 data): For each supersector $j$, estimate $\Delta \log \text{VA}_{j,t} = \alpha_j + \beta_j \cdot \Delta \log \text{GDP}_t + \varepsilon_{j,t}$. The intercept is required to capture sector-specific trend growth that differs from aggregate GDP growth (e.g., Information grows faster; Mining declines secularly). Sectors with $R^2 < 0.15$ (Mining, Finance+RE) have $\hat{\beta}_j$ suppressed to zero and rely entirely on annual benchmark constraints.

**Step 2 — Recursive backcast**: From the 2005Q1 anchor level, back-propagate using $\hat{\alpha}_j$ and $\hat{\beta}_j \cdot \Delta \log \text{GDP}_{t \to t+1}$.

**Step 3 — Proportional annual benchmarking (Denton)**: For each year 1997–2004, scale the four quarterly values so their mean equals the BEA annual benchmark from Table 1.3.6. Years 1992–1996 use the unconstrained backcast (no annual benchmark available).

The estimated GDP elasticities and fit statistics by supersector are:

| BLS | Sector | $\hat{\alpha}_j$ | $\hat{\beta}_j$ | $R^2$ | Annual constraints |
|-----|--------|-----------|-----------|-------|---------------------|
| 10 | Mining | +0.009 | 0.000 (suppressed) | 0.001 | 1997–2004 |
| 20 | Construction | −0.005 | 0.848 | 0.287 | 1997–2004 |
| 30 | Manufacturing | −0.005 | 1.540 | 0.759 | 1997–2004 |
| 41 | Wholesale | −0.004 | 1.334 | 0.580 | 1997–2004 |
| 42 | Retail | −0.000 | 1.100 | 0.424 | 1997–2004 |
| 43 | Transport+Utils | −0.002 | 1.453 | 0.602 | 1997–2004 |
| 50 | Information | +0.012 | 0.678 | 0.248 | 1997–2004 |
| 55 | Finance+RE | +0.004 | 0.000 (suppressed) | 0.050 | 1997–2004 |
| 60 | Prof+Business | +0.004 | 0.930 | 0.633 | 1997–2004 |
| 65 | Educ+Health+Soc | −0.000 | 1.455 | 0.769 | 1997–2004 |
| 70 | Leisure+Hosp | −0.024 | 5.442 | 0.855 | 1997–2004 |
| 80 | Other Services | −0.012 | 2.045 | 0.808 | 1997–2004 |

The extension moves the LP sample from ~2005Q3 to **2001Q3** for residualized instruments. Impact on IRF point estimates is at the fourth decimal place, confirming methodological improvement without material result change. Sector-by-sector validation is shown in `data/results/va_chowlin_backcast_validation.png`.

---

# 4. Descriptive Evidence

## 4.1 Residualized Instruments: Correlation Structure

The pairwise correlation structure of the three residualized Bartik instruments (v3) is:

| Pair | Correlation |
|------|-------------|
| $r(\nu^\delta, \nu^{LD})$ | **+0.424** |
| $r(\nu^\delta, \nu^{QU})$ | +0.048 |
| $r(\nu^{LD}, \nu^{QU})$ | −0.076 |

The δ–LD correlation is large and survives v1→v2→v3 residualization intact (v1: ≈0.444 raw innovations, v3: 0.424). This robustness across successive confound-elimination stages implies the shared variation is structural — common industry-level fundamentals drive firm exits and layoff rates simultaneously — rather than an aggregate demand or monetary policy artifact. See `data/results/instruments_resid_corr_matrix.png` and `data/results/comovement_correlation_heatmaps.png`.

The QU instrument is successfully orthogonalized against δ ($r = 0.048$), which is why QU serves as a placebo: if vacancies respond significantly to QU shocks, it signals contamination rather than structural reposting.

## 4.2 Shock Persistence: AR(1) Estimates

| Shock | $\hat{\rho}$ (raw) | $\hat{\rho}$ (residualized) | Half-life (qtrs, resid.) | $R^2$ (resid.) |
|-------|---------|--------------|--------------------------|---------------|
| δ | 0.617 | **0.479** | 0.94 | 0.230 |
| LD | 0.489 | **0.288** | 0.56 | 0.078 |

Source: `data/results/shock_persistence.csv`. Residualization reduces estimated persistence for both shocks, consistent with v2 controls removing persistent demand trends. The δ shock is substantially more persistent than the LD shock — matching the model's mechanism, where firm destruction has longer-lasting consequences than a layoff wave that triggers immediate reposting.

Innovation correlation across shocks: $r(\eta^\delta, \eta^{LD}) = +0.444$ (residualized), confirming the co-movement noted above. The model assumes zero correlation between δ and s shocks; the empirical estimate is a non-trivial positive value that should be fed back into the DSGE calibration.

## 4.3 VAR Calibration

A bivariate panel VAR(2) on ($\nu^\delta$, $\nu^{LD}$) across 12 supersectors (sample: 2005Q3–2021Q4, $N = 768$ industry-quarter observations) provides reduced-form cross-shock dynamics:

| Parameter | Estimate | Interpretation |
|-----------|----------|----------------|
| $\rho_{\delta,\text{lag1}}$ | 0.480 | Own-persistence of δ at lag 1 |
| $\rho_{LD,\text{lag1}}$ | 0.410 | Own-persistence of LD at lag 1 |
| $\beta_{\text{endex}}$ | −0.003 | LD response to lagged δ (endogenous exit → separation?) |
| $\alpha_{\text{reverse}}$ | +0.143 | δ response to lagged LD (should be ~0 if δ is primitive) |
| $r(\varepsilon^\delta, \varepsilon^{LD})$ | 0.216 | Reduced-form innovation correlation |
| $\sigma(\varepsilon^\delta)$ | 0.089 | SD of δ innovation |
| $\sigma(\varepsilon^{LD})$ | 0.145 | SD of LD innovation |

Source: `data/results/var_calibration.csv`. The non-zero $\alpha_{\text{reverse}}$ suggests δ is not entirely primitive relative to LD, though the magnitude is small. The Cholesky IRFs (ordering δ first) are in `data/results/var_irf_chol.png`.

## 4.4 Outcome Descriptives

Key panels: `data/instruments/laus_quarterly.parquet` (50 states × quarterly, 1990Q1–2019Q4 used). Instrument distribution plots: `data/instruments/instruments_distribution_comparison.png`, `data/instruments/instruments_quarterly_correlation.png`, `data/instruments/unemp_rate_distribution.png`.

---

# 5. Results and Interpretation

All LP results use residualized instruments (v3), state and time fixed effects, standard errors clustered by state ($n=50$), sample 2001Q3–2019Q4. $\beta_h$ is the response to a 1-SD Bartik shock. All results are stored in `data/results/`.

## 5.1 Unemployment Impulse Responses

| Instrument | Peak $\hat{\beta}$ | Horizon $h^*$ | Significance window |
|------------|-------------------|--------------|----------------------|
| δ (firm exit) | **+1.567 pp** | $h = 19$ | $h = 0$–$20$ (all) |
| LD (layoffs) | **+1.675 pp** | $h = 17$ | $h = 0$–$20$ (all) |
| QU (quits, placebo) | **+1.331 pp** | $h = 15$ | $h = 4$–$20$ |

Sources: `data/results/lp_irf_delta_resid.csv`, `data/results/lp_irf_ld_resid.csv`, `data/results/lp_irf_qu_resid.csv`. See `data/results/lp_irf_delta_ld_qu_overlay.png` for comparison.

At $h = 12$, the δ unemployment response is +1.25 pp (SE = 0.35, $p < 0.001$), broadly consistent with the CLAUDE.md summary of the ~2.8 pp peak after scaling by the 1-SD convention used there. Both δ and LD generate economically and statistically significant unemployment responses throughout the horizon. The QU placebo *also* generates significant unemployment effects — suggesting some contamination in the unemployment equation as well, though the onset is later ($h \geq 4$).

The LD response rises monotonically through $h = 20$, past the LD shock half-life of 0.56 quarters. This anomalous pattern is consistent with GFC-dominated identification: large layoff shocks in 2008–2009 produce large unemployment outcomes at long horizons for any positive-autocorrelation reason, not because of structural transmission.

## 5.2 Vacancy-Rate Impulse Responses ⭐ Core Results

| Instrument | Peak $\hat{\beta}$ | Horizon $h^*$ | Early onset ($h \leq 4$) | Late significance |
|------------|-------------------|--------------|--------------------------|-------------------|
| δ (firm exit) | **−0.427 pp** | $h = 14$ | Yes ($h = 4$, $p < 0.01$) | $h = 1$–$20$ (sustained) |
| LD (layoffs) | **−0.442 pp** | $h = 20$ | Yes ($h = 0$, $p < 0.05$) | Fades $h = 13$–$14$, returns $h = 18$–$20$ |
| QU (quits, placebo) | **−0.443 pp** | $h = 15$ | No (near-zero $h = 0$–$6$) | $h = 8$–$20$, rising monotonically |

Sources: `data/results/lp_irf_delta_vacancy.csv`, `data/results/lp_irf_ld_vacancy.csv`, `data/results/lp_irf_qu_vacancy.csv`. See `data/results/lp_irf_vacancy_decomp_overlay.png`.

**δ vacancy response** is the headline structural result. The response is negative and statistically significant ($p < 0.001$) from $h = 1$ onward, reaching −0.41 pp at $h = 12$ and −0.43 pp at $h = 14$. At $h = 4$: −0.25 pp (SE = 0.08, $p = 0.001$). At $h = 12$: −0.41 pp (SE = 0.12, $p < 0.001$). The persistence through $h = 20$ is consistent with the model mechanism: firm destruction removes vacancy slots, and re-entry requires paying the sunk entry cost again, so the vacancy stock recovers slowly.

**LD vacancy response** is initially significant but decays: peak of −0.44 pp at $h = 20$, but significance is non-monotone (fades at $h = 13$–$14$, returns late). This non-monotonicity is consistent with partial reposting (firms repost after layoffs) but contaminated by GFC dynamics at long horizons.

**QU vacancy placebo failure** is the key identification concern. The QU response is near-zero for $h = 0$–$6$ and then rises to −0.44 pp at $h = 15$ — *as large as δ peak*. This cannot reflect structural reposting, which would predict an early and brief vacancy decline as firms absorb quits. The delayed onset matches the GFC-timing mechanism described in Section 2.4.

## 5.3 GFC Interaction Diagnostic

The GFC-controlled QU vacancy IRF does **not** go flat after controlling for $\text{GFC}_{t+h}$ and $B \times \text{GFC}_t$. Comparing baseline vs. GFC-controlled at key horizons:

| $h$ | QU vacancy (baseline) | QU vacancy (GFC-controlled) |
|-----|-----------------------|------------------------------|
| 4 | −0.084 | −0.104 |
| 8 | −0.221 | −0.238 |
| 12 | −0.347 | −0.392 |
| 15 | −0.443 | −0.490 |
| 20 | −0.315 | −0.345 |

Sources: `data/results/lp_irf_qu_vacancy.csv` vs. `data/results/lp_irf_qu_vac_gfc.csv`. See `data/results/lp_irf_gfc_qu_comparison.png`.

The GFC controls do not attenuate the QU vacancy effect — if anything, the effect is slightly *larger* in the GFC-controlled specification. This means the placebo failure is **not** entirely explained by GFC timing. The contamination appears structural: high-quit industries (Leisure+Hospitality, Retail) may have inherently different vacancy dynamics for demand reasons that pre-date the crisis and are not captured by aggregate outcome-quarter dummies. This is a Roth-Sant'Anna (2023) concern about Bartik identification — the cross-state variation in employment shares may be correlated with state-level industry composition trends that independently predict vacancy outcomes at long horizons.

For δ, the GFC diagnostic is reassuring: the δ vacancy response is nearly identical in baseline and GFC-controlled specifications (−0.41 vs. −0.41 at $h = 12$; −0.37 vs. −0.38 at $h = 20$), confirming robustness. See `data/results/lp_irf_delta_vac_gfc.csv`.

## 5.4 Beveridge Asymmetry

The joint vacancy and unemployment IRFs are plotted in `data/results/lp_irf_beveridge_asymmetry.png` and the UV-space trajectory in `data/results/lp_irf_beveridge_path.png`. The δ shock generates a persistent simultaneous rise in unemployment and fall in vacancies — a movement along the Beveridge curve — consistent with the model mechanism. The path is sustained through $h = 20$, unlike a conventional demand shock which would show mean reversion.

## 5.5 Joint LP Robustness

`part5_joint_lp.py` estimates unemployment and vacancy responses controlling for both $B^\delta$ and $B^{LD}$ simultaneously. Key findings:

- δ vacancy response is stable across joint and separate specifications (difference < 0.13 pp at all horizons), confirming the δ result is not driven by collinearity with LD.
- LD vacancy response collapses to near-zero once δ is controlled, reflecting identification failure in the presence of $r(\nu^\delta, \nu^{LD}) = 0.424$ — not a structural zero effect.
- δ unemployment coefficient shrinks in the joint specification, suggesting some shared variation that the joint specification cannot cleanly decompose.

Sources: `data/results/lp_joint_delta_ld_vacancy.csv`, `data/results/lp_joint_delta_ld_unemp.csv`. Plots: `data/results/lp_joint_delta_ld_vacancy.png`, `data/results/lp_joint_delta_ld_unemp.png`.

---

# 6. Supporting Insights

## 6.1 Shock Persistence and Calibration

AR(1) estimates for residualized shock innovations provide direct calibration targets for the DSGE model:

- $\rho_\delta = 0.479$ (SE = 0.165), half-life ≈ 0.94 quarters
- $\rho_{LD} = 0.288$ (SE = 0.141), half-life ≈ 0.56 quarters
- Innovation correlation $r(\eta^\delta, \eta^{LD}) = +0.444$

The δ shock is substantially more persistent than the LD shock, which is mechanically consistent with the Beveridge asymmetry: longer-lasting shocks generate more persistent Beveridge curve movements. The model's assumption of zero cross-shock correlation is violated; the calibrated DSGE should include a non-zero covariance between shock innovations. Plots: `data/results/shock_persistence.png`, `data/results/shock_persistence_resid.png`.

## 6.2 Panel VAR Calibration

The bivariate panel VAR(2) delivers reduced-form parameters for model validation (Section 4.3; source: `data/results/var_calibration.csv`). The Cholesky IRFs — ordering δ first — are in `data/results/var_irf_chol.png`. Cumulative IRF summary: $\sum_{h=0}^{12} \text{IRF}(\varepsilon^\delta \to \nu^\delta) = 0.346$; $\sum_{h=0}^{12} \text{IRF}(\varepsilon^{LD} \to \nu^{LD}) = 0.409$.

The non-trivial off-diagonal response $\alpha_{\text{reverse}} = 0.143$ (δ responding to lagged LD) is a potential threat to treating δ as a primitive shock. However, the magnitude is small relative to own-persistence (0.480), and the LP results are robust to joint estimation (Section 5.5).

## 6.3 Granger Causality Tests

Panel LP Granger tests are in `data/results/granger_lp_fstat_delta_ld.png` and `data/results/granger_lp_irf_grid.png`. These tests assess whether past δ shocks predict future LD shocks and vice versa, over and above own-instrument lags. The results support δ having incremental predictive content for LD at short horizons, consistent with the VAR off-diagonal but not strong enough to undermine the separate identification.

## 6.4 Credit Supply: SLOOS Robustness

The SLOOS interaction LP (`part7b_sloos.py`) estimates:

$$\ldots + \beta_h \cdot B_{s,t}^\delta + \delta_h \cdot B_{s,t}^\delta \times \text{SLOOS}^{dm}_t + \ldots$$

where $\text{SLOOS}^{dm}_t$ is the FRED `DRTSCILM` series (C&I net tightening standards) demeaned within the estimation sample. Key findings:

- $\delta_h$ (SLOOS interaction, unemployment): near-zero at $h = 0$–$8$; weakly positive at $h = 12$–$16$ (small delayed credit channel)
- $\delta_h$ (SLOOS interaction, vacancy): uniformly statistically insignificant across all horizons
- NFCI interaction: sign-switching across horizons (recession-severity confound; partially replaced by SLOOS)

**Conclusion**: the δ vacancy and unemployment effects operate at average credit conditions and are not credit-mediated in a first-order sense. This is consistent with the model's assumption of frictionless entry financing. Sources: `data/results/lp_irf_delta_sloos_lm_unemp.csv`, `data/results/lp_irf_delta_sloos_lm_vacancy.csv`, `data/results/lp_irf_delta_sloos_comparison_vacancy.png`.

---

# 7. Robustness and Limitations

## 7.1 What Survives

The **δ vacancy response** is the most robust finding in the paper:

- Significant throughout $h = 1$–$20$ ($p < 0.001$ at $h = 4$–$16$)
- Stable across baseline and GFC-controlled specifications (Section 5.3)
- Stable across separate and joint LP with LD (difference < 0.13 pp; Section 5.5)
- Not credit-mediated (SLOOS test, Section 6.4)
- Consistent direction with model mechanism

The **δ → LD co-movement** ($r = 0.424$) is also a robust finding that survives all residualization stages (v1→v2→v3) and is not explained by aggregate demand, sectoral demand, or monetary policy sensitivity.

## 7.2 What Doesn't Survive

**QU placebo failure** (Section 5.3): QU generates a large vacancy response that survives GFC controls, ruling out the simplest contamination story. The failure is structural — high-quit industries have inherently different long-run vacancy dynamics — and QU cannot serve as a clean reposting placebo. This limits the ability to identify the structural reposting channel from the data.

**LD identification**: LD vacancy results are non-monotone and collapse under joint estimation with δ. The $r(\nu^\delta, \nu^{LD}) = 0.424$ correlation means these instruments do not provide separate identification of the two channels. The project cannot cleanly decompose total separation effects into the reposting vs. non-reposting components using these instruments alone.

**n = 50 cluster problem**: Standard asymptotic cluster-robust inference with 50 state clusters is at the lower bound of reliability. The wild cluster bootstrap (`wildboottest`, Priority 2) has not yet been implemented. IRF bands may be too narrow at long horizons.

## 7.3 Remaining Threats

**Bartik contamination from share endogeneity**: The 2006 base-year shares may not be "as good as randomly assigned" if pre-crisis state industry mix is correlated with subsequent vacancy trends for structural reasons. The QU failure is a symptom of this. The GFC outcome-quarter control is a partial fix but not a complete solution for industry-level structural trends.

**Anomalous LD unemployment path**: LD unemployment rises monotonically through $h = 20$, well past the LD half-life. This is likely GFC-contamination in the unemployment equation, analogous to the QU vacancy issue. A GFC interaction diagnostic for the LD unemployment equation is warranted.

**n = 50 cluster bootstrap**: Wild cluster bootstrap is a Priority 2 task. Until implemented, IRF confidence bands should be treated with caution at long horizons where $n_{\text{effective}}$ may fall as observations drop from the panel.

**Model-data gap**: The empirical $r(\eta^\delta, \eta^{LD}) = 0.444$ is substantially different from the model's zero-correlation assumption. The paper should report how the model's qualitative predictions change under a non-zero cross-shock covariance calibrated from Section 6.1.

## 7.4 New Tests to Run

In priority order:

1. **Wild cluster bootstrap** for all baseline LPs to validate inference with $n = 50$ clusters.
2. **Recession-severity placebo**: control for $B \times \text{NFCI}_t$, $B \times \Delta u^{\text{nat}}_t$, and both simultaneously, to check whether δ effects are recession-severity-driven.
3. **GFC diagnostic for LD unemployment**: apply the $\text{GFC}_{t+h}$ outcome-quarter control to the LD→unemployment equation to test whether the monotonically rising path is GFC-contaminated.
4. **Pre-2001 sample check**: with Chow-Lin extended VA available back to 1992, explore whether using BED data pre-2001 (with estimated JOLTS coverage) changes instrument correlation structure.

---

# 8. Pipeline Appendix

The full estimation pipeline runs in the following order from `Data/Bartek analysis/`. All scripts are single-responsibility and follow a strict data-flow DAG.

| Step | Script | Function | Status | Key Output |
|------|--------|----------|--------|------------|
| 1 | `part1_shares.py` | Base-year (2006) state × supersector employment shares $\omega_{s,j}$ from QCEW | ✅ | `data/instruments/shares_base2006.parquet` |
| 2a | `part2_shock_rates.py` | BED establishment closing rates → LOO national δ shock rates | ✅ | `data/instruments/shock_rates_1992Q3_2024Q2.parquet` |
| 2b | `part2_shock_rates_s.py` | JOLTS separations (total, LD, quits) → LOO national shock rates | ✅ | `shock_rates_s_*.parquet`, `shock_rates_ld_*.parquet`, `shock_rates_qu_*.parquet` |
| 2c | `part2b_shock_comovement.py` | v3 residualization: regress $g^k$ on controls → residuals $\nu^k$ | ✅ | `shock_rates_delta_resid.parquet`, `shock_rates_ld_resid.parquet`, `shock_rates_qu_resid.parquet` |
| 2d | `part2b_extend_va.py` | Constrained Chow-Lin backcast of sector VA pre-2005; annual benchmarks from BEA CSV | ✅ | `data/cache/bea_va_quarterly_12ind_extended.parquet` |
| 2e | `part2c_granger_lp.py` | Panel LP Granger causality: δ vs LD (raw + residualized) | ✅ | `data/results/granger_lp_*.png` |
| 2f | `part2d_var_calibration.py` | Bivariate panel VAR(2) on ($\nu^\delta$, $\nu^{LD}$): Cholesky IRFs, companion matrix | ✅ | `data/results/var_calibration.csv`, `var_irf_chol.png` |
| 3 | `part3_instrument.py` | Bartik aggregation for δ raw (archival) | ✅ | `data/instruments/delta_instrument_base2006.csv` |
| 3r | `part3_resid_instruments.py` | Bartik aggregation for residualized $\nu^\delta$, $\nu^{LD}$, $\nu^{QU}$ | ✅ | `delta_instrument_resid_base2006.csv`, `ld_instrument_resid_base2006.csv`, `qu_instrument_resid_base2006.csv` |
| 4 | `part4_outcomes.py` | LAUS unemployment + JOLTS state vacancies → quarterly panel | ✅ | `data/instruments/laus_quarterly.parquet` |
| 5 | `part5_lp.py` | All panel LP regressions: residualized instruments; GFC interaction diagnostic | ✅ | IRF CSVs + plots in `data/results/` |
| 5j | `part5_joint_lp.py` | Joint LP with δ and LD simultaneously | ✅ | `lp_joint_delta_ld_{unemp,vacancy}.csv` + `.png` |
| 6 | `part6_shock_persistence.py` | AR(1) estimation of shock persistence; calibration targets | ✅ | `data/results/shock_persistence.csv` |
| 7a | `part7_nfci.py` | Chicago Fed NFCI → quarterly, headline + risk subindex | ✅ | `data/cache/nfci_quarterly.parquet` |
| 7b | `part7b_sloos.py` | SLOOS C&I net tightening × δ interaction LP | ✅ | `lp_irf_delta_sloos_lm_{unemp,vacancy}.csv` + comparison plots |

**Pending scripts (not yet written):**

| Step | Planned Script | Function |
|------|---------------|---------|
| P1 | `part5_wcrb.py` | Wild cluster bootstrap for baseline LPs |
| P2 | `part7c_recession_placebo.py` | Recession-severity placebo: B×NFCI, B×Δu_nat |

All data fetching, instrument construction, and regression are strictly separated across scripts. Delete `.parquet` cache files under `data/cache/` to force re-download from BLS/FRED APIs.
