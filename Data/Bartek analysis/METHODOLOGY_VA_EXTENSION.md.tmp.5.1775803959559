# VA Extension Methodology: Growth-Rate Chow-Lin Temporal Disaggregation

**File:** `part2b_extend_va.py`
**Date:** April 10, 2026
**Author:** Mario Silva

---

## 1. The Problem

The v2 residualization of national industry shock rates (part2b) controls for lagged
industry value-added growth, Δlog VA_{j,t−1}, to remove demand contamination from the
Bartik instruments. The BEA publishes quarterly real value-added by sector (FRED series
RVAM, RVAC, RVAMA, etc.) starting in 2005Q1. The delta-shock instrument (BED
establishment closings) is available from 1992Q3. This mismatch truncates the LP sample:

| Shock series | Raw coverage | v2 LP sample (after VA merge) | Loss |
|---|---|---|---|
| δ (BED closings) | 1992Q3–2019Q4 | 2005Q3–2019Q4 | ~53 quarters, 48% |
| LD (JOLTS layoffs) | 2001Q1–2019Q4 | 2005Q3–2019Q4 | ~18 quarters, 24% |

The 1992Q3–2004Q4 period spans the full 1990s expansion, the 2001 recession, and the
jobless recovery. These are informative cycles for identifying Beveridge curve dynamics
— the 2001 recession in particular featured a sharp rise in establishment exit rates
with a slow vacancy recovery, a canonical episode for distinguishing δ from s
transmission. Discarding this period is empirically costly.

The goal is to construct plausible quarterly Δlog VA_{j,t} series for each of the 12
BLS supersectors back to 1992Q1 using available annual BEA data and a high-frequency
aggregate indicator.

---

## 2. Why Linear Interpolation Is Insufficient

The simplest extension distributes annual value-added growth evenly across four
quarters: if sector j's annual VA grew by g^A_j in year y, assign quarterly growth
g^A_j / 4 to each of Q1–Q4. This is exact linear interpolation on the log-level.

**The problem** is that this assumes demand conditions are constant within each calendar
year. That assumption is most severely violated precisely when the VA control matters
most — during recessions. The 2001 NBER contraction ran from March to November,
concentrated in Q2–Q3. Linear interpolation assigns equal VA declines to all four
quarters of 2001, smoothing out the within-year concentration. As a result, the
interpolated Δlog VA_{j,t−1} fails to capture the demand shock timing that generated
the large establishment-closing spikes in 2001Q2–Q3, and the residualization does not
remove the demand component in those quarters. The control variable is artificially flat
precisely where it should be large.

More formally, let ε^demand_{j,t} denote the industry-specific demand shock that is
correlated with both g^δ_{j,t} (closings) and g^LD_{j,t} (layoffs). The
residualization projects out Δlog VA_{j,t−1} to absorb ε^demand. If Δlog VA^interp
has near-zero quarterly variance within recession years, the projection removes little
of the demand variation, and the residualized shocks ν^δ and ν^LD remain correlated
with ε^demand in exactly those quarters.

---

## 3. The Approach: Growth-Rate Chow-Lin

Chow-Lin temporal disaggregation (Chow and Lin, 1971) is the standard technique for
distributing low-frequency (annual) benchmark totals into high-frequency (quarterly)
estimates using a correlated indicator series. It is the method used by the BEA itself,
Eurostat, the OECD, and most central banks when official quarterly data are unavailable
for a historical backperiod.

We implement a **growth-rate formulation** with aggregate real GDP (FRED: GDPC1) as the
high-frequency indicator. GDP is available quarterly from 1947Q1 and is highly
correlated with sector VA growth in the post-2005 sample where both are observed.

The full procedure has four steps.

### Step 1: Estimate the GDP Elasticity of Sector VA Growth

For each BLS supersector j, run OLS over the post-2005 sample where both quarterly
sector VA (FRED) and quarterly real GDP (GDPC1) are observed:

```
Δlog VA_{j,t} = β_j × Δlog GDP_t + ε_{j,t}
```

**No constant.** GDP growth rates and VA growth rates are both approximately mean-zero
after abstracting from a low-frequency trend; including a constant would absorb
differences in average sectoral growth rates that are irrelevant for the within-year
timing we are trying to capture. In part2b, the residualization includes industry fixed
effects (within-industry demeaning) that absorb the sectoral trend, so the mean of
Δlog VA_{j,t} within the estimation sample is already close to zero.

The OLS estimator is:

```
β̂_j = Cov(Δlog GDP_t, Δlog VA_{j,t}) / Var(Δlog GDP_t)
```

This identifies the **short-run GDP elasticity** of sector j's value added: the
percentage increase in sector VA associated with a one-percent increase in aggregate
GDP within a quarter. Economically:

- **Manufacturing, Construction:** β̂ > 1. Output is highly cyclical (durable goods
  demand is volatile; construction responds to interest rates and credit conditions).
  Both sectors can rapidly adjust production and hours.
- **Wholesale and Retail trade:** β̂ ≈ 0.8–1.2. Moderately cyclical; inventory
  dynamics create some amplification.
- **Professional services, Finance:** β̂ ≈ 0.5–1.0. Less cyclical; long-term
  contracts, sticky revenues.
- **Education, Health, Utilities, Other services:** β̂ < 0.5. Largely acyclical;
  demand is driven by demographics and regulation rather than the business cycle.

The R² from each sector's regression is reported and serves as a validation check: high
R² confirms that GDP is a strong predictor of sector VA timing and the backcast is
well-grounded; low R² signals sector-specific idiosyncrasy and warrants caution.

**Why growth rates rather than log-levels?**

The log-levels of both sector VA and GDP are integrated I(1) processes. Regressing
log(VA_{j,t}) on log(GDP_t) in levels without testing for cointegration risks spurious
correlation — the long-run trends in both series mechanically produce a high R² that
does not reflect a structural relationship at business-cycle frequencies. First-
differencing removes the stochastic trend and makes both series approximately I(0),
rendering OLS consistent and asymptotically normal under standard conditions.

More practically: the control variable used in part2b is Δlog VA_{j,t−1}, a growth
rate, not a level. Estimating the relationship in the same units (growth rates) as the
application avoids any inconsistency between estimation and use.

### Step 2: Backward Recursive Backcast

The key insight is that we know the anchor level — the observed quarterly VA at 2005Q1
— and want to infer earlier levels by running the growth-rate relationship in reverse.
Define the anchor as t* = 2005Q1. We iterate **backward** through time, one quarter at
a time.

At each step, suppose we know log VA_{j,t+1} (the level one quarter ahead). The
growth-rate relationship implies:

```
Δlog VA_{j,t→t+1} = β̂_j × Δlog GDP_{t→t+1}
```

which means:

```
log VA_{j,t} = log VA_{j,t+1} − β̂_j × [log GDP_{t+1} − log GDP_t]
```

Starting from log VA_{j,t*} = log VA_{j,2005Q1} (observed), we step backward through
2004Q4, 2004Q3, ..., 1992Q1, computing the implied log-level at each quarter from the
log-level one period ahead and the GDP growth in that quarter. This is equivalent to
the forward recursion:

```
log VA_{j,t} = log VA_{j,t*} − β̂_j × [log GDP_{t*} − log GDP_t]
             = log VA_{j,t*} + β̂_j × [log GDP_t − log GDP_{t*}]
```

which expresses the backcasted level as the anchor level plus β̂_j times the cumulative
GDP deviation from the anchor. The backward iteration and this closed form are
numerically identical; the backward recursion is used in the code because it naturally
handles missing GDP quarters and is less susceptible to accumulation of floating-point
errors over long horizons.

**Annual consistency.** The sum of within-year quarterly growth rates equals:

```
Σ_{q in year y} Δlog VA_{j,q} = β̂_j × Σ_{q in year y} Δlog GDP_q
                               = β̂_j × [log GDP_{y,Q4} − log GDP_{y,Q1}]
```

which is exactly β̂_j times the within-year GDP growth. This is a natural consistency
property: the backcast implies that sector j's annual VA growth tracks GDP growth
scaled by the sector-specific elasticity. No additional annual benchmark correction
(such as Denton-Cholette proportional adjustment) is needed, because we do not have
independent pre-2005 annual benchmarks for the FRED RVA series — they simply do not
exist before 2005Q1.

### Step 3: Splice and Extract Growth Rates

The backcasted log-level series (1992Q1–2004Q4) is concatenated with the observed FRED
quarterly levels (2005Q1 onward). The splice is continuous by construction: the growth
rate from 2004Q4 to 2005Q1 is β̂_j × Δlog GDP_{2004Q4→2005Q1}, which is the
prediction from the estimated relationship. No level jump is introduced.

Log-differences of the full spliced series give Δlog VA_{j,t} for all quarters from
1992Q2 onward. These are saved to the extended cache and consumed by part2b with a
one-quarter lag.

---

## 4. What This Recovers and What It Does Not

**What the extended series captures:** the GDP-correlated component of sector VA
dynamics, i.e., the part of within-year variation that moves with the aggregate
business cycle. This is precisely the channel we want to control for, because the
leading hypothesis for the residual correlation r(δ, LD) ≈ 0.42 is that both
instruments respond to industry-specific demand conditions. GDP is the best available
single-frequency indicator of those aggregate demand conditions in the pre-2005 period.

**What it does not capture:** idiosyncratic sector-specific shocks that move VA without
moving aggregate GDP. These are the ε_{j,t} residuals from Step 1. In the post-2005
sample, these residuals are absorbed into the error term of the residualization
equation in part2b. In the pre-2005 backcast, by construction, we set ε_{j,t} = 0 —
the backcasted VA is entirely driven by the GDP factor. This is the key limitation of
the approach.

**Implication for inference.** The backcasted Δlog VA_{j,t−1} in pre-2005 quarters
will be a smoother control than the observed FRED series, because the idiosyncratic
sector variation is zeroed out. This means the residualization will remove slightly
less demand contamination in the pre-2005 period than in the post-2005 period. However:

1. The GDP-correlated component is the dominant source of demand contamination in the
   instrument — it is the aggregate business cycle that creates the strongest
   industry-common demand shocks. Sector-specific idiosyncrasy is less likely to
   generate spurious cross-sector comovement in the Bartik aggregation.

2. The attenuation is systematic, not random, and does not bias the LP coefficients.
   It slightly understates the residualization's effectiveness in the pre-2005 sample,
   making our results conservative.

3. The alternative — dropping pre-2005 observations entirely — is far more costly.
   Losing 53 quarters removes the entire 2001 recession from the delta-shock LP
   sample, which is a major source of identification for the δ transmission channel.

---

## 5. Limitations and Caveats

**Structural stability of β̂_j.** The GDP elasticity is estimated from the post-2005
sample (approximately 2005Q1–2021Q4, ~65 quarters). We then apply it to the pre-2005
period. This assumes β̂_j is stable across the two subperiods. This is a strong
assumption for sectors that underwent significant structural change: retail (Amazon
disruption), information (dot-com cycle), finance (deregulation pre-2000, GFC
post-2007). Researchers should inspect the validation plot for sectors where β̂_j
estimated post-2005 seems implausible as a pre-2005 elasticity.

**No pre-2005 annual benchmark.** The Chow-Lin standard procedure verifies the
backcast against annual total benchmarks. We cannot do this for most sectors because
the FRED RVA annual data also begins in 2005. As a partial substitute, BEA GDP-by-
Industry annual data (available from 1998 via BEA website; not accessible via FRED
API) could be used to pin annual totals for 1998–2004. This is a potential extension
if the backcast is used for a central empirical result rather than a control variable.

**GDP as sole indicator.** The standard Chow-Lin approach allows multiple indicator
series. For manufacturing, industrial production (FRED: INDPRO) is a better quarterly
indicator than aggregate GDP. For construction, housing starts (FRED: HOUST) would add
information beyond GDP. The single-indicator specification is a practical simplification
that could be refined sector by sector if the R² from Step 1 is low.

**Splice artifact.** The 2005Q2 growth rate (first post-splice observation) uses
log VA observed at 2005Q1 and 2005Q2, both from FRED. This introduces no artifact.
The 2004Q4→2005Q1 growth rate uses backcasted 2004Q4 and observed 2005Q1, and equals
β̂_j × Δlog GDP_{2004Q4→2005Q1} by construction. If the true GDP elasticity was
unusually low or high in early 2005 (transition from early-expansion to mid-expansion),
this splice quarter may not perfectly match the FRED series. The validation plot flags
this by overlaying the 2005Q1 observations.

---

## 6. Econometric Specification Note for the Paper

For the paper text, the appropriate description of the VA control in the pre-2005
period is:

> "For quarters prior to 2005Q1, when the BEA quarterly value-added series (available
> via FRED) are not available, we extend the sector-level quarterly VA series using a
> growth-rate Chow-Lin temporal disaggregation with aggregate real GDP (GDPC1) as the
> indicator. For each BLS supersector j, we estimate the short-run GDP elasticity β̂_j
> from OLS of quarterly VA growth on quarterly GDP growth over the post-2005 sample.
> We then back out quarterly VA levels before 2005 by iterating the relationship
> Δlog VA_{j,t} = β̂_j × Δlog GDP_t backward from the 2005Q1 anchor. The procedure
> implies that annual sector VA growth tracks β̂_j times annual GDP growth, providing
> automatic annual consistency without requiring annual benchmark data. The backcasted
> VA control absorbs the GDP-correlated component of industry demand conditions in the
> pre-2005 period; idiosyncratic sector variation (which requires the quarterly FRED
> series to observe) is omitted, making the pre-2005 residualization slightly more
> conservative than the post-2005 residualization."

---

## 7. Implementation Details

**Script:** `part2b_extend_va.py`
**Run before:** `part2b_residualize_shocks.py`

**Inputs:**
- `data/cache/bea_va_quarterly_12ind.parquet` — FRED RVA levels, 2004Q4+ (re-fetched for levels)
- `data/cache/gdpc1_quarterly.parquet` — quarterly real GDP (GDPC1), fetched once and cached

**Outputs:**
- `data/cache/bea_va_quarterly_12ind_extended.parquet` — Δlog VA_{j,t} for 1992Q2+ (backcasted 1992Q2–2004Q4, observed 2005Q2+)
- `data/results/va_chowlin_backcast_validation.png` — 12-panel validation plot

**Auto-integration:** `part2b_residualize_shocks.py` automatically uses the extended cache when present, without any further modifications to part3 or part5.

**Run sequence:**
```bash
# First time (requires FRED API key):
FRED_API_KEY=<key> python part2b_extend_va.py

# Subsequent runs use cached GDP and extended VA; then re-residualize:
python part2b_residualize_shocks.py
python part3_resid_instruments.py
python part5_lp.py
```

**Expected LP sample extension (delta instrument):**

| Sample | Quarters | Includes |
|---|---|---|
| v2 original | 2005Q3–2019Q4 | GFC, post-GFC recovery |
| v2 extended | 1994Q4–2019Q4 | 1990s expansion, 2001 recession, jobless recovery, GFC, post-GFC |

(1994Q4 rather than 1992Q3 because productivity lag reduces the start by one quarter and the VA lag by one more.)

---

## 8. References

- **Chow, G. C., & Lin, A. L.** (1971). "Best linear unbiased interpolation, distribution, and extrapolation of time series by related series." *Journal of the American Statistical Association*, 66(335), 751–757.

- **Denton, F. T.** (1971). "Adjustment of monthly or quarterly series to annual totals: an approach based on quadratic minimization." *Journal of the American Statistical Association*, 66(333), 99–102.

- **Di Fonzo, T., & Marini, M.** (2012). "On the extrapolation with the Denton proportional benchmarking method." *IMF Working Paper* 12/169.

- **BEA Real Value Added by Industry (FRED).** Series RVAM, RVAC, RVAMA, RVAW, RVAR, RVAT, RVAU, RVAI, RVAFI, RVARL, RVAPBS, RVAES, RVAHC, RVAER, RVAAF, RVAOSEG. Available from 2005Q1.

- **Real Gross Domestic Product (FRED).** Series GDPC1 (chained 2017 dollars). Available from 1947Q1.
