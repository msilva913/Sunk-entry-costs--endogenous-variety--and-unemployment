# Handoff: Vacancy LP Outcome — part4 fetch + part5 LP

## Context

This is a continuation of a Bartik instrument pipeline for a DSGE paper
on firm entry/exit and unemployment (Gabrovski-Silva extension with BGM
business formation, three aggregate shocks: δ, s, z).

The pipeline currently estimates panel local projections (LPs) with
**state unemployment rate** as the outcome variable. The task here is
to add **state-level job vacancies** as a second outcome variable.

The model's core asymmetry — δ destroys vacancies, s triggers reposting —
is most directly tested by the vacancy IRF, not the unemployment IRF. This
is the primary motivation.

---

## Working directory

```
Data/Bartek analysis/
```

All scripts run from this directory. Key subdirectories:
```
data/cache/          parquet cache files
data/instruments/    Bartik instrument CSVs and parquets
data/results/        LP output CSVs and plots
```

---

## Data source: JOLTS state-level job openings

**Source:** BLS JOLTS state estimates program  
**Coverage:** All 50 states + DC, total nonfarm, December 2000 – present  
**Nature:** Model-imputed estimates (composite synthetic model using QCEW-LDB
and regional JOLTS data) — not direct survey measurements at state level.
Flag this in any paper discussion.  
**Frequency:** Monthly → sum/average to quarterly  
**Seasonal adjustment:** Use seasonally adjusted (SA) series  

### ⚠️ CRITICAL: Verify series ID format before implementing

The exact BLS API series ID format for state-level JOLTS job openings
**must be verified** before coding the fetch. The national total nonfarm
job openings SA series is `JTSJOL`. State-level series use a different
schema that may involve state FIPS codes.

**Verification step:** Query the BLS LABSTAT flat files or use the BLS
series search at https://data.bls.gov/cgi-bin/surveymost?jt to confirm
the correct format. The state JOLTS data page is at:
https://www.bls.gov/jlt/jlt_statedata.htm

A plausible format based on the JOLTS national series structure is:
`JTS{ST}{000000}{00}{00000}{00}{JO}{L}` where ST is the 2-digit state
FIPS code — but **do not assume this without verification**. Wrong series
IDs return empty results silently.

The existing `fetch_jolts_separations_national()` function in
`construct_s_instrument.py` provides a working template for BLS API
calls including rate limiting, chunking, and caching. Use it as the
model for the vacancy fetch.

---

## Task 1: Add vacancy fetch to part4_outcomes.py

### What part4 currently does

`part4_outcomes.py` fetches LAUS (Local Area Unemployment Statistics)
state-level monthly unemployment rates and labor force, aggregates to
quarterly, and saves:
```
data/instruments/laus_quarterly.parquet
```
Columns: `state_fips`, `state`, `quarter_label`, `unemp_rate`,
`labor_force`

### What to add

Add a new function `fetch_jolts_vacancies_state()` that:

1. Constructs BLS API series IDs for state-level total nonfarm job
   openings (SA, level, thousands) for all 50 states + DC.
2. Fetches via BLS public API v1 (same infrastructure as existing
   JOLTS fetch — see `construct_s_instrument.py`).
3. Aggregates monthly levels to quarterly by **summing** the three
   months within each calendar quarter (consistent with how JOLTS
   separations are aggregated in the existing pipeline).
4. Saves cache:
   ```
   data/cache/jolts_vacancies_state_quarterly.parquet
   ```
   Columns: `state_fips`, `quarter_label`, `vacancies`
   (level in thousands of jobs)

5. Merges with the existing LAUS parquet to produce an augmented
   outcomes file:
   ```
   data/instruments/laus_quarterly.parquet
   ```
   Now includes: `state_fips`, `state`, `quarter_label`, `unemp_rate`,
   `labor_force`, `vacancies` (NaN where not available)

### Sample coverage target
- 2001Q1 – 2023Q1 matching the s instrument sample
- 50 states (DC optional — exclude if it causes issues)

---

## Task 2: Add vacancy LP outcome to part5_lp.py

### Current LP specification

```
y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h B_{s,t}^(k)
                         + γ₁ u_{s,t-1} + γ₂ log(LF_{s,t-1}) + ε
```

where y_{s,t} is the state unemployment rate (%).

### Vacancy LP specification

The vacancy outcome uses **log vacancies** rather than levels because
vacancies have a strong upward trend and vary substantially in scale
across states.

```
log(v_{s,t+h}) - log(v_{s,t-1}) = α_s + α_t + β_h B_{s,t}^(k)
                                   + γ₁ log(v_{s,t-1}) + γ₂ log(LF_{s,t-1})
                                   + ε
```

Key differences from the unemployment LP:
- Outcome: log change in vacancies from t-1 to t+h
- Control: lagged log vacancies replace lagged unemployment rate
- Interpretation: β_h = log-point change in vacancy stock per 1-SD shock

### Implementation in part5_lp.py

The cleanest implementation adds an `outcome` parameter to `run_lp()`:

```python
def run_lp(base_panel, shock_col, label,
           outcome="unemp",      # "unemp" or "vacancy"
           include_nfci=False):
```

When `outcome="vacancy"`:
- `dep_var = log_v_future - log_v_lag1`  (instead of y_future - unemp_lag1)
- Core regressors swap `unemp_lag1` for `log_v_lag1`
- `log(LF)` control retained

The `build_panel()` function needs to:
1. Accept vacancies from the augmented LAUS parquet
2. Compute `log_v = log(vacancies)` and `log_v_lag1` (lagged one quarter,
   same gap-checking logic as `unemp_lag1`)
3. Compute `log_v_future` for each horizon h (same forward-lookup logic
   as `y_future`)

### New LP calls to add

After the existing unemployment LP runs, add:

```python
# Vacancy outcome LPs (if vacancies available)
if "vacancies" in delta_panel.columns:
    irf_delta_vac      = run_lp(delta_post2001, "bartik_delta",
                                 "δ shock — vacancies", outcome="vacancy")
    irf_s_vac          = run_lp(s_panel,        "bartik_s",
                                 "s shock — vacancies",  outcome="vacancy")
    irf_ld_vac         = run_lp(ld_resid_panel, "bartik_ld",
                                 "LD shock — vacancies", outcome="vacancy")

    irf_delta_vac.to_csv(RESULTS_DIR / "lp_irf_delta_vacancy.csv",  ...)
    irf_s_vac.to_csv(    RESULTS_DIR / "lp_irf_s_vacancy.csv",      ...)
    irf_ld_vac.to_csv(   RESULTS_DIR / "lp_irf_ld_vacancy.csv",     ...)

    # Side-by-side plot: δ vacancy vs s/LD vacancy IRFs
    # This is the core asymmetry test: δ should show vacancy decline,
    # LD should show vacancy increase (reposting) or flat
    plot_irf(
        {"δ shock": irf_delta_vac, "LD shock": irf_ld_vac},
        out_path=RESULTS_DIR / "lp_irf_vacancy_delta_ld.png",
    )
```

### Why LD not TS for the vacancy test

The vacancy test is specifically about the reposting channel. LD
(layoffs and discharges at continuing establishments) is the
cleanest analog to the model's s_t because surviving firms can
repost. TS (total separations) conflates quits, which have no
clean reposting prediction. Use LD resid as the primary s-type
instrument for the vacancy LP.

---

## Model prediction being tested

From equation (lom_v) in the model:

```
v_t = (1 - δ_{e,t}) [(1 - q(θ_{t-1})) v_{t-1} + s_{t-1}(1 - u_{t-1})] + e_t
```

- **δ shock**: δ_{e,t} rises → both surviving unmatched vacancies
  (term i) AND reposting by surviving firms (term ii) collapse
  simultaneously → vacancies FALL, unemployment RISES → Beveridge
  curve dynamics confirmed

- **LD shock**: s_{t-1} rises but δ_{e,t} unchanged → term ii EXPANDS
  (surviving firms repost) → vacancies FLAT or RISE, unemployment
  also RISES → positive UV co-movement

If IRFs match these predictions, this is the single strongest
empirical confirmation of the model's core mechanism.

---

## Existing relevant files

| File | Role |
|------|------|
| `part4_outcomes.py` | Fetches LAUS; ADD vacancy fetch here |
| `part5_lp.py` | LP engine; ADD vacancy outcome option here |
| `construct_s_instrument.py` | BLS API fetch template (see `fetch_jolts_separations_national`) |
| `data/instruments/laus_quarterly.parquet` | Current outcomes file |
| `data/instruments/delta_instrument_resid_base2006.csv` | δ residualized instrument |
| `data/instruments/ld_instrument_resid_base2006.csv` | LD residualized instrument |
| `data/instruments/s_instrument_base2006.csv` | TS raw instrument |

---

## Run order after implementation

```
python part4_outcomes.py          # adds vacancies to laus_quarterly.parquet
python part5_lp.py                # runs all LPs including vacancy outcome
```

No other pipeline files need modification.

---

## Notes on unit-SD standardization

`run_lp()` already standardizes β_h by the instrument SD before
returning results. This applies automatically to the vacancy LP
as well — no separate handling needed.

## Notes on COVID cap

`MAX_OUTCOME_QUARTER = "2019Q4"` is set globally. This applies to
vacancy outcomes automatically through the existing `future_ql`
filtering logic — no change needed.

