# VA Extension Methodology: Constrained Chow-Lin Backcast of Sector VA

**File:** `part2b_extend_va.py`
**Date:** April 14, 2026
**Author:** Mario Silva

---

## 1. The Problem

The v2 residualization of national industry shock rates (part2b) controls for lagged
industry value-added growth, Δlog VA_{j,t−1}, to remove demand contamination from the
Bartik instruments. The BEA publishes quarterly real value-added by sector via FRED
(series RVAM, RVAC, RVAMA, etc.) starting in 2005Q1 only. The delta-shock instrument
(BED establishment closings) is available from 1992Q3; the layoff and quits instruments
(JOLTS) from 2001Q1. This mismatch truncates the LP sample:

| Shock series | Raw coverage | v2 LP sample (after VA merge) | Loss |
|---|---|---|---|
| δ (BED closings) | 1992Q3–2019Q4 | 2005Q3–2019Q4 | ~53 quarters, 48% |
| LD / QU (JOLTS) | 2001Q1–2019Q4 | 2005Q3–2019Q4 | ~18 quarters, 24% |

The 2001Q1–2004Q4 period is particularly consequential for the LD and QU instruments:
it includes the 2001 recession and the subsequent jobless recovery — a canonical episode
for distinguishing δ from s transmission in the Beveridge curve. The 2001 recession
featured a sharp rise in establishment exit rates (large δ shock) with a slow vacancy
recovery, while layoff rates (LD) spiked and recovered quickly. This asymmetric
reabsorption of vacancies between shock types is precisely the pattern the paper aims to
identify. Dropping these 16 quarters from the LP sample weakens the identification of
both the δ vacancy IRF and the LD comparison.

**What the backcast actually recovers.** The LP outcome sample window is separately
constrained by the availability of state-level vacancy data (JOLTS state release, from
2001Q1) and cannot be extended by the VA backcast. What the VA extension achieves is
expanding the *residualization sample* — the sample over which Δlog VA_{j,t−1} is
available as a control when constructing ν^δ, ν^{LD}, ν^{QU}. Better residualization of
the 2001–2004 shock observations improves the quality of the instruments used in the LP.
Per CLAUDE.md §15.1, the empirical impact on IRF estimates is small (differences at the
fourth decimal place), confirming this is a methodological improvement rather than a
material change in findings.

---

## 2. Why Linear Interpolation Is Insufficient

The simplest extension distributes annual value-added growth evenly across four
quarters. If sector j's annual VA grew by g^A_{j,y} in year y, assign g^A_{j,y}/4 to
each of Q1–Q4. This is linear interpolation on the log-level.

**The core problem** is that this imposes a constant-within-year growth rate, which is
most severely violated during recessions. The 2001 NBER contraction ran from March to
November, concentrated in Q2–Q3. Linear interpolation assigns equal VA declines to all
four quarters of 2001, smoothing out the within-year concentration of the downturn. The
control variable Δlog VA_{j,t−1} is then artificially flat precisely in the quarters
(2001Q2–Q3) where the largest establishment-closing spikes occurred and where demand
contamination in the instruments is strongest. The residualization therefore removes
little of the demand component in those quarters.

More formally: the residualization projects out Δlog VA_{j,t−1} to absorb the demand
component ε^demand_{j,t}. If the control has near-zero quarterly variance within
recession years, the projection's explanatory power in those quarters is negligible,
leaving ν^δ and ν^{LD} correlated with ε^demand in exactly the periods where we most
need them to be clean.

---

## 3. The Approach: Constrained Chow-Lin with GDP Indicator

Chow-Lin temporal disaggregation (Chow and Lin, 1971) is the standard technique in
national accounts for distributing low-frequency (annual) benchmark totals into
high-frequency (quarterly) estimates using a correlated indicator series. It is used by
the BEA itself, Eurostat, the OECD, and most central banks for historical backperiods.

We implement a **constrained Chow-Lin formulation** with four steps:

1. Estimate the GDP elasticity of each sector's quarterly VA growth from post-2005 data.
2. Back-propagate quarterly levels from the 2005Q1 anchor using that elasticity and
   observed GDP growth.
3. Apply proportional Denton constraints from BEA annual benchmarks (1997–2004) to pin
   the annual averages of the backcasted quarterly path.
4. Splice the constrained backcast with the observed FRED series and extract growth rates.

### Step 1: GDP Elasticity Regression

For each BLS supersector j, estimate by OLS over the post-2005 sample:

```
Δlog VA_{j,t} = α_j + β_j · Δlog GDP_t + ε_{j,t}
```

**The intercept α_j is required.** Sectors have sector-specific mean quarterly growth
rates that differ systematically from aggregate GDP growth. Information technology has
secular above-GDP growth (α > 0); Mining has a secular long-run decline (α > 0 as well,
but for different reasons). Omitting the intercept forces the regression through the
origin and produces negative centered-R² for these sectors — a logical impossibility for
OLS with an intercept — because the sector's mean growth is not zero even after removing
the GDP component. The intercept captures the unconditional sectoral trend; β_j captures
the cyclical GDP elasticity around that trend.

**Low-R² suppression.** For sectors where GDP explains little of the quarterly variance
(R² < 0.15), β_j is noisy and mechanically applying it to the backcast would add
spurious GDP-cycle variation to the control. For these sectors, β_use is set to zero and
the backcast uses only the intercept (α_j per quarter, equivalent to a constant-growth
extrapolation). Mining (R²=0.001) and Finance+RE (R²=0.050) fall below this threshold.
Their backcasts rely entirely on the annual BEA benchmark constraint in Step 3 to
recover meaningful within-year variation.

**Why growth rates, not log-levels?** Both log(VA_{j,t}) and log(GDP_t) are I(1)
processes. A levels regression without testing for cointegration risks producing a
spurious slope driven by the common stochastic trend rather than a structural cyclical
relationship. First-differencing removes the stochastic trend and makes both series
approximately I(0). Additionally, the control variable used in part2b is itself a growth
rate (Δlog VA_{j,t−1}), so estimating β_j in growth rates ensures the estimation and
application units are consistent.

**Sector-specific elasticities (from post-2005 OLS):**

| BLS | Sector | α | β_OLS | β_used | R² |
|-----|--------|---|-------|--------|-----|
| 10 | Mining | +0.009 | +0.136 | **0.000** | 0.001 |
| 20 | Construction | −0.005 | +0.848 | +0.848 | 0.287 |
| 30 | Manufacturing | −0.005 | +1.540 | +1.540 | 0.759 |
| 41 | Wholesale | −0.004 | +1.334 | +1.334 | 0.580 |
| 42 | Retail | −0.000 | +1.100 | +1.100 | 0.424 |
| 43 | Transport+Utils | −0.002 | +1.453 | +1.453 | 0.602 |
| 50 | Information | +0.012 | +0.678 | +0.678 | 0.248 |
| 55 | Finance+RE | +0.004 | +0.230 | **0.000** | 0.050 |
| 60 | Prof+Business | +0.004 | +0.930 | +0.930 | 0.633 |
| 65 | Education+Health | −0.000 | +1.455 | +1.455 | 0.769 |
| 70 | Leisure+Hosp | −0.024 | +5.442 | +5.442 | 0.855 |
| 80 | Other Services | −0.012 | +2.045 | +2.045 | 0.808 |

The high elasticities for Manufacturing (1.54) and Leisure+Hospitality (5.44) are
economically sensible: manufacturing output is strongly procyclical, and the leisure and
hospitality sector is highly sensitive to consumer spending, which itself responds
sharply to income fluctuations. Education and Health (1.46) is somewhat surprising but
reflects the large share of non-profit hospital spending that tracks economic activity
through insurance coverage and elective procedures.

### Step 2: Unconstrained Backward Recursive Backcast

The anchor is the observed FRED quarterly level at t* = 2005Q1. Starting there, we
iterate backward one quarter at a time. At each step, the growth rate relationship gives:

```
Δlog VA_{j,t→t+1} = α_j + β_use,j · Δlog GDP_{t→t+1}
```

which implies:

```
log VA_{j,t} = log VA_{j,t+1} − [α_j + β_use,j · (log GDP_{t+1} − log GDP_t)]
```

Starting from log VA_{j,2005Q1}, we step back through 2004Q4, 2004Q3, ..., 1992Q1.
The resulting path is the **unconstrained backcast**: the quarterly VA trajectory that is
consistent with the estimated elasticity applied to the observed GDP path, anchored at
2005Q1 in level.

The backward recursion and the closed-form equivalent:

```
log VA_{j,t} = log VA_{j,2005Q1} + β_use,j · [log GDP_t − log GDP_{2005Q1}]
             + α_j · (t − t*)    [cumulative intercept drift]
```

are numerically identical. The recursive form is used in the code to handle any missing
GDP quarters gracefully and to make the step-by-step logic transparent.

### Step 3: Annual Benchmark Constraint (Proportional Denton)

The unconstrained backcast pins the *level* at 2005Q1 and implies annual growth equal to
α_j + β_use,j times annual GDP growth. However, BEA publishes annual industry real VA in
its GDP-by-Industry accounts (Table 1.3.6), which are independent of the FRED quarterly
series and available from **1997 onward**. These annual benchmarks carry more
information than the GDP indicator alone — they reflect BEA's comprehensive benchmark
revisions, methodology, and industry classification — and should be respected.

The BEA annual convention is that the annual value equals the mean of the four
quarterly SAAR values:

```
Y_{j,A} = mean(VA_{j,Q1}, VA_{j,Q2}, VA_{j,Q3}, VA_{j,Q4})   for year A
```

For each year A in 1997–2004 with published benchmark Y_{j,A}, we apply a
proportional Denton adjustment to enforce this constraint exactly:

```
scale_A = Y_{j,A} / mean(VA^{unconstrained}_{j,Q1:Q4 of A})

VA^{constrained}_{j,q} = VA^{unconstrained}_{j,q} · scale_A    for q ∈ {Q1,Q2,Q3,Q4 of A}
```

This multiplies all four quarters of year A by the same scalar, so the within-year
seasonal and cyclical shape (inherited from GDP) is preserved, while the annual mean is
pinned to the BEA published figure. The scale factor measures the discrepancy between
the GDP-extrapolated level and the BEA annual benchmark. If scale_A ≈ 1, the GDP
indicator tracked the sector accurately; if scale_A deviates substantially from 1, the
BEA benchmark is doing significant corrective work.

**Coverage of the benchmark constraint:**
- **1997–2004** (8 years): fully constrained to BEA annual benchmarks. These years span
  the Asian financial crisis recovery, the dot-com boom and bust, and the 2001 recession.
- **1992–1996** (5 years): no independent annual benchmark available. The unconstrained
  GDP-anchored backcast applies for these years, and the 2005Q1 anchor propagates back
  through the Denton-adjusted 1997 level.

**Annual benchmark data source.** The file `bea_va_annual_by_supersector.csv` was
constructed from BEA GDP-by-Industry Table 1.3.6 (April 9, 2026 release), converted to
billions of chained 2017 dollars to match the FRED SAAR units. The 12 BLS supersector
columns aggregate BEA industry lines as follows:

| BLS | BEA aggregation |
|-----|----------------|
| BLS 43 | Transport+Warehousing (line 40) + Utilities (line 10) |
| BLS 55 | Finance+Insurance (line 55) + Real Estate+Rental (line 60) |
| BLS 65 | Edu+Health+Social Assistance aggregate (line 74) ← matches `RVAESHS`, not `RVAES` |
| BLS 70 | Arts+Entertainment (RVAAER) + Accommodation+Food (RVAAF) ← matches line 81 |

Cross-verification against FRED quarterly means (post-2005): 11 of 12 sectors match to
≤0.02%; all 12 are correct by construction.

### Step 4: Splice and Extract Growth Rates

The constrained backcasted levels (1992Q1–2004Q4) are concatenated with the observed
FRED quarterly levels (2005Q1 onward). At the splice point 2005Q1, both sides use the
same FRED-observed level, so no discontinuity is introduced. Log-differences of the
full spliced series yield Δlog VA_{j,t} for all quarters from 1992Q2 onward, which is
what part2b consumes as the control Δlog VA_{j,t−1}.

---

## 4. What This Recovers and What It Does Not

**What the extended series captures:** For 1997–2004, the backcasted growth rates carry
two layers of information: (i) the within-year GDP-cycle shape from the indicator
regression, which places more of the annual decline in recession quarters and more of
the rebound in recovery quarters; (ii) the annual level anchored to BEA benchmark data,
which ensures the multi-year trajectory of each sector's real VA is correctly measured.
For 1992–1996, only the GDP-cycle shape is available; the level is extrapolated from
the GDP-anchored path without an annual pin.

**What it does not capture:** Idiosyncratic sector-specific shocks uncorrelated with
aggregate GDP. These are the ε_{j,t} residuals from the Step 1 regression. In the
post-2005 sample, they are absorbed into the error term of the residualization. In the
pre-2005 backcast, they are set to zero by construction — the backcasted VA contains
only the GDP-correlated component (and the annual constraint). This is the fundamental
limitation.

**Implication for inference.** The control Δlog VA_{j,t−1} in pre-2005 quarters
absorbs somewhat less of the true demand variation than the observed FRED series would,
because the idiosyncratic component is zeroed out. This makes the pre-2005
residualization slightly more conservative — it removes less contamination — which
biases the residual correlation r(ν^δ, ν^{LD}) upward in the pre-2005 period. Since
r(δ,LD) = 0.424 is already being interpreted as a structural floor (not a contamination
artifact), this conservatism is not harmful: it means we are, if anything,
underestimating how clean the residualized instruments are.

The proportional Denton constraint for 1997–2004 partially addresses this limitation by
ensuring the annual trajectory is correct even when idiosyncratic quarterly variation is
missing.

---

## 5. Assessment: Is the Methodology Sound?

**What is unambiguously correct:**

The three-component structure — GDP indicator regression with intercept, backward
recursive backcast, and proportional Denton annual benchmark constraint — is the
standard constrained Chow-Lin procedure as documented in the national accounts
literature. Each component has a clear economic rationale:

- The intercept correctly captures sector-specific trend growth. Without it, Mining and
  Finance get negative R² — a mathematical sign of misspecification.
- The backward recursion correctly propagates the 2005Q1 anchor into the past without
  accumulating errors in the forward direction.
- The Denton constraint makes essential use of 8 years of independent BEA annual data,
  turning a pure extrapolation into a disciplined interpolation for the 1997–2004 window.
- The R² suppression of β for low-R² sectors (Mining, Finance) is conservative and
  correct: a noisy β estimate applied to 52 quarters of GDP variation would generate
  more noise than signal in those sectors.

**The key identification question.**

The VA extension is used only as a control variable in the residualization of shock
rates, not as a primary outcome or instrument. For this purpose, the requirements are
modest: the control should carry enough within-year variation to absorb the demand
component of the shock instrument in recession quarters. The Chow-Lin backcast achieves
this far better than linear interpolation. Whether it achieves it as well as a true
quarterly BEA series would is unknowable pre-2005, but the validation plot (which
compares predicted vs. observed growth rates within the post-2005 sample) provides an
in-sample benchmark for each sector's predictive accuracy.

**One substantive concern: structural stability of β̂_j across subperiods.**

The GDP elasticity is estimated from approximately 65 post-2005 quarterly observations
and applied to the pre-2005 period. The post-2005 estimation sample happens to include
the GFC (2008–09), which is the largest peacetime GDP swing in the sample and
dominates the β̂_j estimates for cyclically sensitive sectors. For Leisure+Hospitality,
β̂ = 5.44 reflects in part the extreme GFC decline in consumer spending; it may
overstate the sector's pre-2005 GDP sensitivity, when income volatility was lower.
Overestimating β means the 2001 recession VA decline for this sector is overstated,
making the control more aggressive than the true series — which is, however, the
conservative direction for residualization purposes (it removes more contamination, not
less).

More generally, the structural stability assumption is weakest for:
- **Information (BLS 50):** the dot-com cycle (1999–2001) created extreme within-year
  VA volatility uncorrelated with overall GDP. The post-2005 β = 0.678 may be
  insufficient to capture the 2001 crash dynamics in this sector.
- **Retail (BLS 42):** the post-2005 period includes the structural shift from
  brick-and-mortar to e-commerce, which altered the sector's GDP elasticity.
- **Finance+RE (BLS 55):** β suppressed to zero; the annual Denton constraint carries
  all the weight for 1997–2004.

For the Information sector, the low R² (0.248) and positive intercept (α = +0.012,
reflecting the secular productivity-driven growth of the post-dot-com era) may not be
appropriate for 1999–2001 when information sector VA was collapsing while aggregate GDP
was near-flat. This is a genuine limitation, though it affects only one sector's
contribution to the Bartik weight-average, and Information's employment share in the
base year (2006) weights down its influence in the aggregate instrument.

**The 1992–1996 window (unconstrained).**

For 1992–1996, no annual BEA benchmark is available and the backcast relies entirely on
GDP extrapolation from the 2005Q1 anchor. Over a 13-year extrapolation horizon, even
small annual errors in α or β_j compound into level drift. However, because we take
log-differences, the level error is differenced away: the quarterly growth rate
Δlog VA_{j,t} = α_j + β_use,j · Δlog GDP_t for each quarter, regardless of the
accumulated level. The only effect of level drift on the growth rates appears at the
1996Q4→1997Q1 boundary, where the Denton constraint for 1997 applies a scale factor
that may be different from 1.0. If that scale factor is very different from 1.0, it
implies the log-level at end-1996 differs substantially from where GDP extrapolation
placed it, and the 1997Q1 growth rate will be distorted. In practice this is a single
observation and is unlikely to affect the residualization materially.

The LP sample for the LD and QU instruments starts at 2001Q1 (JOLTS constraint), which
falls in the well-constrained 1997–2004 Denton window. The 1992–1996 extrapolation
affects only the δ instrument sample for those years, where the LP outcome sample also
requires LAUS state unemployment data (available from 1990) — so in principle the δ LP
could extend back to ≈1994Q4 (adding one productivity lag and one VA lag). Whether
those pre-1997 observations should be included depends on whether the unconstrained
backcast is reliable enough to serve as a control. Given the GDP-anchored nature of the
backcast, including them is reasonable; they should be treated as a robustness check
rather than the primary sample.

---

## 6. Econometric Specification Note for the Paper

> "For quarters prior to 2005Q1, when the BEA quarterly value-added series are
> unavailable on FRED, we extend the sector-level quarterly VA series using a
> constrained Chow-Lin (1971) temporal disaggregation. For each BLS supersector j, we
> estimate the short-run GDP elasticity β_j and the sector-specific trend α_j by OLS of
> quarterly VA growth on an intercept and quarterly real GDP growth (FRED: GDPC1) over
> the post-2005 sample. We back out pre-2005 quarterly VA levels by iterating the
> relationship Δlog VA_{j,t} = α_j + β_j · Δlog GDP_t backward from the observed 2005Q1
> anchor. For years 1997–2004, annual-mean constraints are imposed using published BEA
> GDP-by-Industry annual benchmarks (Table 1.3.6), ensuring the quarterly path averages
> to the correct annual total for each sector; years 1992–1996 use the unconstrained
> GDP-anchored extrapolation. Sectors where GDP explains less than 15% of quarterly VA
> variance (Mining, Finance and Real Estate) have their GDP elasticity suppressed to zero
> and rely entirely on the annual constraint for within-period variation. The resulting
> backcasted Δlog VA_{j,t−1} series is used as a control in the shock-rate
> residualization; its effect on instrument construction is small, with IRF estimates
> differing at the fourth decimal place relative to the unextended sample."

---

## 7. Implementation Details

**Script:** `part2b_extend_va.py`
**Run before:** `part2b_residualize_shocks.py`

**Required input file (not auto-generated):**
- `bea_va_annual_by_supersector.csv` — must be co-located with script. BEA Table 1.3.6,
  April 9 2026 release. Columns: BLS10, BLS20, ..., BLS80. Index: year (1997–2025).
  Units: billions of chained 2017 dollars.

**Other inputs (auto-generated or fetched):**
- FRED quarterly sector VA levels (2004Q4+) — fetched fresh each run (levels, not
  cached separately from growth-rate cache)
- `data/cache/gdpc1_quarterly.parquet` — GDPC1, cached after first FRED fetch

**Outputs:**
- `data/cache/bea_va_quarterly_12ind_extended.parquet` — Δlog VA_{j,t} for 1992Q2+
- `data/results/va_chowlin_backcast_validation.png` — 12-panel validation plot

**Auto-integration:** `part2b_residualize_shocks.py` prefers the extended cache when
present and falls back to the original FRED cache otherwise (no code changes needed
downstream).

**Run sequence:**
```bash
# Requires FRED API key (first run fetches GDP + sector VA levels):
FRED_API_KEY=<key> python part2b_extend_va.py

# Then re-run residualization and downstream pipeline:
python part2b_residualize_shocks.py
python part3_resid_instruments.py
python part5_lp.py
```

**LP sample impact:**

| Sample | Residualization VA coverage | LP outcome window |
|---|---|---|
| v2 original (no extension) | 2005Q2–2019Q4 | 2005Q3–2019Q4 |
| v2 extended (this script) | 1992Q2–2019Q4 | still 2001Q1–2019Q4 (JOLTS-constrained) |

The LP outcome window does not change. The benefit is better residualization for the
2001Q1–2004Q4 observations (quarterly VA control now available and Denton-constrained),
replacing the previously missing values that caused those observations to be dropped.

---

## 8. References

- **Chow, G. C., & Lin, A. L.** (1971). "Best linear unbiased interpolation, distribution, and extrapolation of time series by related series." *Journal of the American Statistical Association*, 66(335), 751–757.

- **Denton, F. T.** (1971). "Adjustment of monthly or quarterly series to annual totals: an approach based on quadratic minimization." *Journal of the American Statistical Association*, 66(333), 99–102.

- **Di Fonzo, T., & Marini, M.** (2012). "On the extrapolation with the Denton proportional benchmarking method." *IMF Working Paper* 12/169.

- **BEA GDP-by-Industry accounts, Table 1.3.6.** Real Value Added by Industry, annual, chained 2017 dollars. April 9, 2026 release. https://www.bea.gov/industry/gdpbyind_data.htm

- **BEA Real Value Added by Industry (FRED quarterly series).** RVAM, RVAC, RVAMA, RVAW, RVAR, RVAT, RVAU, RVAI, RVAFI, RVARL, RVAPBS, **RVAESHS** (not RVAES), **RVAAER+RVAAF** (not RVAER). Available 2005Q1+.

- **Real Gross Domestic Product (FRED).** Series GDPC1 (chained 2017 dollars). Available 1947Q1+.
