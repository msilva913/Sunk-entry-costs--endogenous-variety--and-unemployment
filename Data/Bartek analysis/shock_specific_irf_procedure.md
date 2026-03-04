# Shock-Specific Impulse Responses: Implementation Procedure

## Overview

The goal is to estimate state-level impulse responses to **product destruction shocks** (δ) and **match separation shocks** (s) separately, and show that the vacancy–unemployment comovement differs across shock types — directly validating the model's central mechanism before structural estimation.

**Core predictions to test:**

| Outcome | After δ shock | After s shock |
|---|---|---|
| Unemployment | Persistent elevation | Transient spike |
| Vacancies | Persistent collapse | Rapid recovery or overshoot |
| Establishment exit | Elevated for years | Transient spike, then normalization |

---

## Step 1: Construct State-Level Shock Instruments

### Shared framework

Both instruments follow the Bartik (shift-share) design:

$$B^k_{s,t} = \sum_j \omega_{s,j,t_0} \cdot g^k_{-s,j,t}$$

where $\omega_{s,j,t_0}$ are pre-recession employment shares (state $s$, supersector $j$, base year $t_0$) and $g^k_{-s,j,t}$ is a leave-one-out national shock rate for shock type $k$. The ten supersectors are the standard BLS private nonfarm supersectors (Mining, Construction, Manufacturing, Trade/transport/utilities, Information, Financial activities, Professional and business services, Education and health, Leisure and hospitality, Other services).

**Identification**
Bartik instrument exploits interaction between predetermined local industry composition and national industry-level shocks to generate plausibly exogenous variation in local labor market conditions.

The identifying assumption following Goldsmith-Pinkham, Sorkin and Swift (2020) is that **shares are exogenous** conditional on controls. The pre-determined industry composition $\omega_{s,j,t_0}$ should be uncorrelated with error term in the outcome equation, so that states did not sort into high-closing or high-layoff indusries in anticipation of future labor market outcomes.

### Employment shares $\omega_{s,j,t_0}$

**Source:** Quarterly Census of Employment and Wages (QCEW), annual average employment by state and supersector.

**Base year:** 2006 (pre-Great Recession). A single base year is used for both instruments to ensure they operate on the same weighting scheme and to hold the composition of local economies fixed.

**Construction:** For each state $s$ and supersector $j$, the share is:

$$\omega_{s,j,t_0} = \frac{E_{s,j,t_0}}{E_{s,t_0}}$$

where $E_{s,t_0}$ is total private nonfarm employment in state $s$ in the base year. Shares sum to one within each state across the ten supersectors.






**Retrieval:** QCEW data are fetched via the BLS public API v1 (`api.bls.gov/publicAPI/v1/timeseries/data/`). Series IDs follow the format `ENU{state_fips_5d}{seasonal}{supersector_code}` — for example, `ENU010000105101` for Alabama Mining, seasonally adjusted. The full set requires 500 series (50 states × 10 supersectors) plus 10 national series. Because the API v1 unregistered limit is 10 years per request and 25 series per request, calls are batched into groups of 25 series and three year-chunks (1990–2000, 2001–2010, 2011–present). Results are cached as parquet after the first pull, which is an efficient open-source columnar storage file format.

**LOO denominator:** For each (state, supersector) cell, the leave-one-out national base employment is computed from the shares directly:

$$E^{nat}_{-s,j,t_0} = E^{nat}_{j,t_0} - E_{s,j,t_0}$$

This quantity is stored in the shares parquet and reused by both instruments. It is the only object that varies across states when computing the shock rate $g^k_{-s,j,t}$.

---

### δ instrument — Business Employment Dynamics (BED)

$$B^\delta_{s,t} = \sum_j \omega_{s,j,t_0} \cdot g^\delta_{-s,j,t}, \qquad g^\delta_{-s,j,t} = \frac{\text{closings}^{nat}_{j,t}}{E^{nat}_{-s,j,t_0}}$$

**Interpretation:** Captures the establishment-exit component of separations — the $\delta$ channel in the model where both jobs and vacancies are simultaneously destroyed. BED gross job losses from closing establishments are the natural empirical counterpart to $\delta_t$ because the closure of an establishment destroys all associated vacancies along with all filled positions.

**Numerator — data source and retrieval:**

BED gross job losses from *closing* establishments, by supersector, national level, quarterly. Fetched via the BLS public API v1. Series IDs follow the format `BDS0000000000{bed_code}110006LQ5` — for example, `BDS0000000000100010110006LQ5` for Mining. The BED industry code system differs from QCEW: Trade/transport/utilities has no single aggregate in BED and must be constructed by summing Wholesale trade, Retail trade, and Transportation/warehousing/utilities sub-components. The remaining nine supersectors map directly.

Only **national** series are fetched (area code `000000000000`). State-level BED data broken down by supersector are not published via the BLS public API. Ten national series require one batch call with three year-chunks (1992–2001, 2002–2011, 2012–present), totalling three API calls. Results are cached as parquet.

**The BED series yields quarterly closing levels in units of jobs lost** (thousands). This is divided by $E^{nat}_{-s,j,t_0}$ (thousands of workers, from QCEW base year) to produce a dimensionless quarterly rate.

**Numerator variation across states:** The closing level $\text{closings}^{nat}_{j,t}$ is a single national number for each (supersector, quarter). It does **not** vary across states. The only state-specific quantity entering the shock rate is the LOO denominator $E^{nat}_{-s,j,t_0}$, which differs across states by the amount of that state's own base employment.

**LOO correction — denominator only, applied universally:**

Because state-level BED supersector closing series are unavailable, a full numerator correction (subtracting state closings from the national total) cannot be implemented. The denominator-only LOO correction is applied universally:

$$g^\delta_{-s,j,t} = \frac{\text{closings}^{nat}_{j,t}}{E^{nat}_{j,t_0} - E_{s,j,t_0}}$$

This is the standard approach in the Bartik literature (Goldsmith-Pinkham et al. 2020, Autor et al. 2013). The denominator correction does the essential econometric work: it prevents a state's own base employment from entering the denominator and creating a mechanical correlation between the instrument and local outcomes. The numerator correction, which would additionally subtract state closings from the national total, is a second-order refinement that matters most for states that constitute a large fraction of national employment in a given supersector (e.g. Wyoming in Mining). Omitting it is the standard practice when state-level data are unavailable.

**Sample coverage:** 1992Q3–2024Q2. BED was introduced in 1992Q3 and is available for the full post-introduction period.

---

### s instrument — JOLTS Layoffs and Discharges

$$B^s_{s,t} = \sum_j \omega_{s,j,t_0} \cdot g^s_{-s,j,t}, \qquad g^s_{-s,j,t} = \frac{\text{layoffs}^{nat}_{j,t}}{E^{nat}_{-s,j,t_0}}$$

**Interpretation:** Captures within-firm separations at surviving employers — the $s$ channel in the model where separated workers' vacancies are immediately reposted by the firm. JOLTS layoffs and discharges exclude job losses from establishment closings by construction (BLS reports them as a distinct category), making them a clean proxy for the match separation component of $\tau$.

**Numerator — data source and retrieval:**

JOLTS layoffs and discharges (seasonally adjusted, level in thousands), by supersector, national level. Fetched via the BLS public API v1. Series IDs follow the 21-character JOLTS format `JTS{industry(6)}{state(2)}{area(5)}{sizeclass(2)}{dataelement(2)}{ratelevel(1)}` — for example, `JTS1100990000000LDL` for national Mining layoffs. The parameters are: seasonal = `S` (SA), state = `00` (national), area = `00000`, sizeclass = `00` (all sizes), dataelement = `LD` (layoffs and discharges), ratelevel = `L` (level). All ten supersectors have direct JOLTS codes; Trade/transport/utilities is published as a single aggregate (code `400000`), requiring no sub-component summing.

Only **national** series are fetched. JOLTS does publish some geographic breakdowns, but only at the four-region level (Northeast, Midwest, South, West) and only for total private — not by supersector. State-level JOLTS supersector data do not exist in the public API. Ten national series require one batch call with three year-chunks (2000–2008, 2009–2017, 2018–present), totalling three API calls. Results are cached as parquet.

**Monthly to quarterly aggregation:** JOLTS publishes monthly levels (thousands of layoffs occurring during the month). The three calendar months within each quarter are summed to produce a quarterly flow level. Incomplete quarters — those for which fewer than three months were returned by the API — are dropped. The resulting quarterly levels are in the same units as the QCEW annual base employment denominator (both in thousands of workers).

**Numerator variation across states:** Exactly as for the δ instrument: `layoffs_nat_{j,t}` is a single national number per (supersector, quarter) and does **not** vary across states. Cross-state variation in $g^s_{-s,j,t}$ arises exclusively from the LOO denominator.

**LOO correction — denominator only, applied universally:**

State-level JOLTS layoff series by supersector are not available via the public API. The denominator-only correction is therefore applied universally, for the same reasons as the δ instrument:

$$g^s_{-s,j,t} = \frac{\text{layoffs}^{nat}_{j,t}}{E^{nat}_{j,t_0} - E_{s,j,t_0}}$$

**Sample coverage:** 2001Q1–2024Q2. JOLTS was launched in December 2000; January–March 2001 constitutes the first complete quarter for which all three monthly observations are available. The joint sample for regressions using both instruments therefore begins in 2001Q1.

---

### Summary of what varies by state vs. nationally

| Object | Varies by state? | Source |
|---|---|---|
| Employment share $\omega_{s,j,t_0}$ | Yes | QCEW, base year 2006 |
| LOO denominator $E^{nat}_{-s,j,t_0}$ | Yes (state subtracted) | QCEW, base year 2006 |
| BED closing level $\text{closings}^{nat}_{j,t}$ | **No** — one number per (supersector, quarter) | BED via BLS API v1 |
| JOLTS layoff level $\text{layoffs}^{nat}_{j,t}$ | **No** — one number per (supersector, quarter) | JOLTS via BLS API v1 |
| Shock rate $g^k_{-s,j,t}$ | Yes, through denominator only | Computed |
| Instrument $B^k_{s,t}$ | Yes | Weighted sum |

---

## Step 2: Construct Outcome Variables

Two outcomes suffice for the core mechanism test:

1. **Unemployment rate:** BLS LAUS, monthly seasonally adjusted for all 50 states, averaged to quarterly frequency.

2. **Vacancies:** State-level job postings. Options in order of preference:
   - Indeed/Burning Glass job postings at state level (~2007 onward via NBER)
   - Conference Board Help Wanted Online index by state (2005 onward)
   - Conference Board print Help Wanted Index by metro area for pre-2005 (noisier)

Employment is **not needed** as a separate outcome. The unemployment rate captures matching market dynamics, and employment conflates these with hours, composition effects, and participation margin changes that are outside the model. The Beveridge curve is defined in $(u, v)$ space.

---

## Step 3: Local Projections

**Frequency:** Quarterly throughout. BED is only available quarterly. JOLTS monthly data are summed to quarterly. Outcome variables are averaged to quarterly. There is no gain from running monthly local projections since the binding data constraint (BED) is quarterly, and doing so would require mixing observation frequencies.

For each outcome $y \in \{\text{unemployment rate, vacancies}\}$, horizon $h$, and shock type $k \in \{\delta, s\}$, estimate:

$$\Delta y_{s,h}^{(k)} = \alpha_h^{(k)} + \beta_h^{(k)} \cdot B_s^{(k)} + \gamma_h^{(k)\prime} X_s + \varepsilon_{s,h}^{(k)}$$

where:
- $\Delta y_{s,h}$ is the change in outcome in state $s$ at horizon $h$ quarters relative to the recession peak
- $B_s^{(k)}$ is the Bartik instrument for shock type $k$
- $X_s$ is a vector of pre-recession state controls (log population, pre-recession employment growth, census region indicators)
- Standard errors are HC1 robust; supplement with permutation p-values from 1,000 random reassignments of the exposure measure

**Horizon grid:** $h \in \{1, 2, 4, 8, 12, 16, 20\}$ quarters.

**Sign convention:**
- δ instrument: exposed states should show persistently higher unemployment and persistently lower vacancies
- s instrument: exposed states should show a transient unemployment spike and rapid vacancy recovery

---

## Step 4: The Beveridge Curve Signature

Plot $\hat\beta_h^{(\delta)}$ for unemployment and $\hat\beta_h^{(\delta)}$ for vacancies on the same time axis, and repeat for $\hat\beta_h^{(s)}$. The model predicts:

- **δ shock:** the two IRFs move in *opposite* directions throughout — unemployment rises, vacancies fall, consistent with the Beveridge curve
- **s shock:** the two IRFs briefly move *together* at short horizons (surviving firms repost vacancies, pushing vacancies up alongside unemployment), then both normalize

Sharpen with the bivariate regression at each horizon:

$$\Delta v_{s,h} = \alpha_h + \gamma_h \cdot \Delta u_{s,h} + \delta_h \cdot B_s^{(k)} + \varepsilon_{s,h}$$

Test whether $\hat\gamma_h < 0$ (Beveridge-consistent negative comovement) for the δ instrument, and $\hat\gamma_h \geq 0$ at short horizons for the s instrument.

---

## Step 5: Use IRFs to Discipline Structural Estimation

### Option A — Additional SMM moments

Simulate the model under a parameter draw, impose a δ shock and an s shock separately, and compute model-implied LP coefficients by running the same regression on simulated data. Augment the moment vector with:

$$m^{LP} = \left\{\hat\beta_h^{(\delta)},\ \hat\beta_h^{(s)}\right\}_{h \in \{2,4,8,12,20\},\ y \in \{u, v\}}$$

### Option B — Indirect inference

Define the binding function $b(\theta) = \mathbb{E}_\theta[\hat\beta_h^{(\delta)}, \hat\beta_h^{(s)}]$ and estimate $\theta$ by minimizing the distance between empirical LP coefficients and the model-implied binding function.

### Option C — Just-identified moments for shock parameters

Use empirical half-lives from the LP as direct targets:
- Half-life of $\hat\beta_h^{(\delta)}$ for unemployment identifies $\rho_\delta$
- Half-life of $\hat\beta_h^{(s)}$ for vacancies identifies $\rho_s$
- Ratio of peak $\hat\beta^{(\delta)}$ to peak $\hat\beta^{(s)}$ for vacancies identifies $\bar\delta / \bar s$ conditional on structural parameters

---

## Data Sources and Pipeline Files

| Object | Source | API | Frequency | Coverage | Pipeline file |
|---|---|---|---|---|---|
| Employment shares $\omega_{s,j,t_0}$ | QCEW | BLS API v1 | Annual (base year only) | 1990–present, all states | `part1_shares.py` |
| LOO denominator $E^{nat}_{-s,j,t_0}$ | QCEW | (derived from shares) | Base year only | — | `part1_shares.py` |
| δ shock numerator | BED closings, national × supersector | BLS API v1 | Quarterly | 1992Q3–present | `part2_shock_rates.py` |
| s shock numerator | JOLTS layoffs, national × supersector | BLS API v1 | Monthly → quarterly | 2001Q1–present | `part2_shock_rates_s.py` |
| δ instrument $B^\delta_{s,t}$ | Derived | — | Quarterly | 1992Q3–2024Q2 | `part3_instrument.py` |
| s instrument $B^s_{s,t}$ | Derived | — | Quarterly | 2001Q1–2024Q2 | `part3_instrument_s.py` |
| Unemployment rate | BLS LAUS | BLS API | Monthly → quarterly | 1976–present, all states | TBD |
| Vacancies | Indeed/Burning Glass | — | Monthly → quarterly | ~2007–present | TBD |

**Joint sample for regressions using both instruments:** 2001Q1 onward (bound by JOLTS availability).

---

## Implementation Notes

**API rate limits:** The BLS public API v1 (no registration key) imposes a limit of 25 series per request and 10 years per request. All pipeline scripts split requests into batches of 25 series and three year-chunks of roughly 9 years each. A 1-second sleep is inserted between calls. All raw API results are cached as parquet files in `data/cache/` after the first pull; subsequent runs read from cache and skip all network calls.

**Stale cache warning:** If the API returns fewer quarters than expected, delete the relevant parquet cache file and rerun. Silent truncation occurs when year-chunks exceed the 10-year unregistered limit — the API returns data without an error message.

**Trade/transport/utilities in BED:** This supersector has no single national aggregate in BED. It is constructed by summing gross job losses from closing establishments across Wholesale trade (`20`), Retail trade (`44-45`), and Transportation, warehousing, and utilities (`48-49`, `22`) sub-components. All other nine supersectors map directly to a single BED series.

**Trade/transport/utilities in JOLTS:** Unlike BED, JOLTS publishes Trade/transport/utilities as a direct aggregate (industry code `400000`). No sub-component summing is required.

**Seasonality:** All BED and JOLTS series used are seasonally adjusted at source. No additional seasonal adjustment is applied in the pipeline.

**Dominant-state concern:** The denominator-only LOO means the numerator is not purged of state $s$'s own closings or layoffs. For states that constitute a large fraction of national employment in a given supersector (e.g. Wyoming in Mining, ~7% of national Mining employment), a small mechanical correlation between the instrument and the local economy remains. This is negligible for most state-supersector cells and is standard practice in the literature when state-level national data sources are unavailable.
