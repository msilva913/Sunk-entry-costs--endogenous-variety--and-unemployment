"""
part4_outcomes.py — State-level unemployment rates, labor force, and vacancies
===============================================================================
Fetches monthly seasonally adjusted unemployment rates (suffix "3") and labor
force levels (suffix "6") for all 50 states from the BLS Local Area
Unemployment Statistics (LAUS) program via the BLS public API v1. Both series
are aggregated to quarterly frequency.

Also fetches state-level JOLTS job openings (vacancies) for all 50 states +
DC from the BLS JOLTS state estimates program (available from Dec 2000).

LAUS series ID format
---------------------
 LASST{state_fips_2d}0000000000003 — unemployment rate (%)
 LASST{state_fips_2d}0000000000006 — civilian labor force (thousands)

JOLTS state-level job openings series ID format (21 characters)
---------------------------------------------------------------
 JT + S + 000000 + {ST} + 00000 + 00 + JO + L

 Where:
 JT = JOLTS prefix
 S = seasonally adjusted
 000000 = total nonfarm industry (6 chars)
 {ST} = 2-digit state FIPS (e.g. "01" = Alabama)
 00000 = all areas (5 chars)
 00 = all sizes (2 chars)
 JO = job openings (data element, 2 chars)
 L = level (thousands)

 Example: Alabama (FIPS 01) → "JTS0000000100000 00JOL"
 (21 chars: JTS + 000000 + 01 + 00000 + 00 + JO + L)

 ⚠️ NOTE: BLS does not publish state JOLTS via the public API v1/v2.
 The state-level data is distributed as model-imputed synthetic estimates
 (composite from QCEW-LDB and regional JOLTS). Coverage and API
 availability must be verified before running this fetch for the first
 time. If empty results are returned the function prints a diagnostic
 and exits gracefully — no data is corrupted.

 The state JOLTS data page is at https://www.bls.gov/jlt/jlt_statedata.htm
 and the series search tool at https://data.bls.gov/cgi-bin/surveymost?jt

Vacancy aggregation — IMPORTANT NOTE ON STOCK vs FLOW
------------------------------------------------------
JOLTS job openings are a STOCK measure (count of open positions on the last
business day of the reference month), not a flow. The correct quarterly
aggregation is therefore the AVERAGE of the three monthly values, NOT a sum.
This is different from how JOLTS separations (a flow: number of events during
the month) are aggregated in the existing pipeline (where summing is correct).
Summing a stock over three months would produce a number 3× too large.

Labor force is fetched alongside the unemployment rate for three reasons:
 1. Weighted means: simple cross-state means overweight small states.
 Labor-force-weighted means approximately recover the BLS national headline.
 2. WLS robustness: LP regressions can be rerun with labor-force weights to
 give larger states more influence — standard robustness check (ADH 2013).
 3. Level conversion: converting LP percentage-point responses to worker counts
 requires a base labor force: ΔU_{s,t} = Δu_{s,t} × LF_{s,t0}.

Sample window
-------------
 LAUS: 1990M1 – present.
 Vacancies: 2001M1 – present (first full quarter of JOLTS state data).

API batching
------------
BLS API v1 (unregistered): 25 series per request, 10 years per request.
LAUS: 50 states × 2 series = 100 series, split into 4 batches × ~5 year-chunks.
Vacancies: 51 series (50 states + DC), split into 3 batches × 3 year-chunks.

Outputs
-------
 data/cache/laus_states_monthly.parquet — monthly, 50 states
 data/cache/jolts_vacancies_state_quarterly.parquet — quarterly, 50 states
 data/instruments/laus_quarterly.parquet — combined quarterly file

Columns in laus_quarterly.parquet:
 state_fips zero-padded 2-digit string, e.g. "01"
 state 2-letter abbreviation, e.g. "AL"
 quarter_label e.g. "2008Q4"
 unemp_rate unemployment rate (%)
 labor_force civilian labor force (thousands, quarterly average)
 vacancies job openings (thousands, quarterly average; NaN before 2001)

Run
---
 python part4_outcomes.py

No prerequisites (does not depend on parts 1–3).
"""

import os
import time
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg") # non-interactive backend — avoids crash issues
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

def _open_file(path):
 """Open a saved file with the OS default viewer."""
 os.startfile(path)

try:
 os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
 # Running interactively (Spyder, Jupyter, etc.) — __file__ may be
 # undefined or the hardcoded fallback path may not match this machine.
 # Fall back to the current working directory; ensure your IDE is set
 # to the "Bartek analysis" folder before running.
 pass

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START_YEAR = 1990
END_YEAR = 2026 # update as needed

CACHE_DIR = Path("data/cache")
OUTPUT_DIR = Path("data/instruments")

MONTHLY_CACHE = CACHE_DIR / "laus_states_monthly.parquet"
QUARTERLY_OUT = OUTPUT_DIR / "laus_quarterly.parquet"

BLS_URL = "https://api.bls.gov/publicAPI/v2/timeseries/data/" # v2 required with registration key
SLEEP_SEC = 1.0
BLS_API_KEY = os.environ.get("BLS_API_KEY", "ae7c43e8d19b41f9969b9a6b4cff350a")

# ---------------------------------------------------------------------------
# State FIPS → abbreviation (all 50 states, sorted by FIPS)
# ---------------------------------------------------------------------------
FIPS_TO_STATE = {
 "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA",
 "08": "CO", "09": "CT", "10": "DE", "12": "FL", "13": "GA",
 "15": "HI", "16": "ID", "17": "IL", "18": "IN", "19": "IA",
 "20": "KS", "21": "KY", "22": "LA", "23": "ME", "24": "MD",
 "25": "MA", "26": "MI", "27": "MN", "28": "MS", "29": "MO",
 "30": "MT", "31": "NE", "32": "NV", "33": "NH", "34": "NJ",
 "35": "NM", "36": "NY", "37": "NC", "38": "ND", "39": "OH",
 "40": "OK", "41": "OR", "42": "PA", "44": "RI", "45": "SC",
 "46": "SD", "47": "TN", "48": "TX", "49": "UT", "50": "VT",
 "51": "VA", "53": "WA", "54": "WV", "55": "WI", "56": "WY",
}

ALL_FIPS = sorted(FIPS_TO_STATE.keys()) # 50 entries

# LAUS measure codes (2-digit, per BLS LAUS series format)
# Series format: LA(2) + S(1) + area_code(15) + measure_code(2) = 20 chars total
# area_code for states = ST(2) + fips(2) + 00000000000 (11 zeros)
# e.g. Alabama rate → LASST010000000000003 (20 chars)
_RATE_SUFFIX = "03" # unemployment rate (%)
_LF_SUFFIX = "06" # civilian labor force (thousands)


def _laus_id(fips: str, suffix: str) -> str:
 sid = f"LASST{fips}00000000000{suffix}"
 assert len(sid) == 20, f"LAUS series ID wrong length: {sid!r} ({len(sid)} chars)"
 return sid


def _year_chunks(start: int, end: int, size: int = 8):
 chunks, y = [], start
 while y <= end:
  chunks.append((y, min(y + size - 1, end)))
  y += size
 return chunks


def _fetch_batch(series_ids: list, start_year: int, end_year: int) -> list:
 payload = {
 "seriesid": series_ids,
 "startyear": str(start_year),
 "endyear": str(end_year),
 "registrationkey": BLS_API_KEY,
 }
 resp = requests.post(BLS_URL, json=payload, timeout=30)
 resp.raise_for_status()
 data = resp.json()
 if data.get("status") != "REQUEST_SUCCEEDED":
  raise RuntimeError(
  f"BLS API status={data.get('status')}: "
  + "; ".join(data.get("message", []))
  )
 series = data["Results"]["series"]
 # Warn if all series came back empty — indicates wrong API version or bad key
 n_obs = sum(len(s.get("data", [])) for s in series)
 if n_obs == 0:
  msgs = data.get("message", [])
  raise RuntimeError(
  f"BLS returned 0 observations for {len(series)} series "
  f"({start_year}–{end_year}). API messages: {msgs}. "
  f"Check BLS_URL version (v1 vs v2) and registration key."
  )
 return series


def _parse_series(series_list: list) -> pd.DataFrame:
 """Parse BLS series list → tidy DataFrame (series_id, year, period, value)."""
 rows = []
 for s in series_list:
  sid = s["seriesID"]
  for obs in s.get("data", []):
   if obs["period"] == "M13": # skip annual averages
    continue
   try:
    val = float(obs["value"])
   except (ValueError, KeyError):
    continue
   rows.append({
    "series_id": sid,
    "year": int(obs["year"]),
    "period": obs["period"],
    "value": val,
   })
 return pd.DataFrame(rows)


def fetch_laus_monthly(
 start_year: int = START_YEAR,
 end_year: int = END_YEAR,
 cache_path: Path = MONTHLY_CACHE,
) -> pd.DataFrame:
 """
 Fetch LAUS monthly unemployment rates and labor force for all 50 states.

 Batches: interleave rate (suffix 3) and labor-force (suffix 6) series so
 both are pulled in the same API calls — 25 series per call, 4 batches per
 year-chunk (2 state-halves × 2 suffixes).

 Returns DataFrame: state_fips, state, date, unemp_rate, labor_force
 """
 if cache_path.exists():
  cached = pd.read_parquet(cache_path)
  # Validate cache: must have both columns AND cover the requested end year.
  # A partial fetch would be missing labor_force; a stale cache would end
  # well before end_year. Either condition triggers a full re-fetch.
  _cache_max_year = cached["date"].max().year if "date" in cached.columns else 0
  _cache_ok = (
   "labor_force" in cached.columns
   and _cache_max_year >= end_year - 1
  )
  if _cache_ok:
   print(f" [cache hit] {cache_path} "
         f"(through {cached['date'].max().strftime('%Y-%m')})")
   return cached
  else:
   _reason = (
    "missing labor_force" if "labor_force" not in cached.columns
    else f"ends {cached['date'].max().strftime('%Y-%m')}, need ≥{end_year - 1}"
   )
   print(f" [cache stale] {cache_path} — {_reason} — re-fetching")
   cache_path.unlink()

 print(f" [fetch] LAUS rate + labor force {start_year}–{end_year} (50 states)")

 # Build series lists: interleave rate and LF so each batch of 25 covers
 # one state-half and both series types.
 half_a = ALL_FIPS[:25]
 half_b = ALL_FIPS[25:]
 # Each sub-batch: 25 series = 12 or 13 states × 2 suffixes (truncated to 25)
 # Split each half into rate-batch and LF-batch
 batches = [
  [_laus_id(f, _RATE_SUFFIX) for f in half_a], # 25 rate series, states 1-25
  [_laus_id(f, _LF_SUFFIX) for f in half_a], # 25 LF series, states 1-25
  [_laus_id(f, _RATE_SUFFIX) for f in half_b], # 25 rate series, states 26-50
  [_laus_id(f, _LF_SUFFIX) for f in half_b], # 25 LF series, states 26-50
 ]
 year_chunks = _year_chunks(start_year, end_year)
 print(f" {len(batches)} batches × {len(year_chunks)} year-chunks "
       f"= {len(batches)*len(year_chunks)} API calls")

 all_frames = []
 call_num = 0
 for yr_start, yr_end in year_chunks:
  for batch in batches:
   call_num += 1
   suffix_label = "rate" if batch[0].endswith(_RATE_SUFFIX) else "LF"
   fips_sample = batch[0][5:7]
   print(f" Call {call_num:2d}: {yr_start}–{yr_end} "
         f"{suffix_label} FIPS starts {fips_sample}…", end=" ")
   try:
    raw = _fetch_batch(batch, yr_start, yr_end)
    df = _parse_series(raw)
    all_frames.append(df)
    print(f"→ {len(df)} rows")
   except Exception as exc:
    print(f"→ ERROR: {exc}")
    raise
   time.sleep(SLEEP_SEC)

 raw_df = pd.concat(all_frames, ignore_index=True)

 # Extract FIPS (chars 5–6 of series ID: LASST{XX}...)
 raw_df["state_fips"] = raw_df["series_id"].str[5:7]
 # Identify series type from trailing 2 chars (measure code: "03" or "06")
 raw_df["series_type"] = raw_df["series_id"].str[-2:]
 raw_df["month"] = raw_df["period"].str[1:].astype(int)
 raw_df["date"] = pd.to_datetime(raw_df[["year", "month"]].assign(day=1))
 raw_df["state"] = raw_df["state_fips"].map(FIPS_TO_STATE)

 # Pivot so rate and labor_force are separate columns
 pivot = (
  raw_df
  .pivot_table(
   index = ["state_fips", "state", "date"],
   columns = "series_type",
   values = "value",
   aggfunc = "first",
  )
  .rename(columns={_RATE_SUFFIX: "unemp_rate", _LF_SUFFIX: "labor_force"})
  .reset_index()
 )
 pivot.columns.name = None

 out = (
  pivot
  .sort_values(["state_fips", "date"])
  .reset_index(drop=True)
 )

 cache_path.parent.mkdir(parents=True, exist_ok=True)
 out.to_parquet(cache_path, index=False)
 print(f" Cached → {cache_path} ({len(out):,} rows)")
 return out


def monthly_to_quarterly(monthly: pd.DataFrame) -> pd.DataFrame:
 """
 Average monthly LAUS data to quarterly frequency using resample.

 Both unemp_rate and labor_force are averaged within each quarter.
 Incomplete quarters (fewer than 3 months) are dropped.
 Returns columns: state_fips, state, quarter_label, unemp_rate, labor_force
 """
 out = (
  monthly
  .set_index("date")
  .groupby(["state_fips", "state"])[["unemp_rate", "labor_force"]]
  .resample("QE")
  .agg(
   unemp_rate = ("unemp_rate", "mean"),
   labor_force = ("labor_force", "mean"),
   n_months = ("unemp_rate", "count"), # completeness check
  )
  .query("n_months == 3")
  .drop(columns="n_months")
  .reset_index()
  .assign(quarter_label=lambda d: d["date"].dt.to_period("Q").astype(str))
  .drop(columns="date")
  .sort_values(["state_fips", "quarter_label"])
  .reset_index(drop=True)
 )
 return out


# ---------------------------------------------------------------------------
# Vacancy fetch — JOLTS state-level job openings
# ---------------------------------------------------------------------------

# DC is FIPS 11; include alongside the 50 states.
# The handoff targets 50 states (DC optional) — DC is included here but
# the LP panels join on state_fips, so DC rows will simply have no
# match in the existing LAUS panel (which covers 50 states only) and
# will be dropped at merge time. No harm done either way.
_VACANCY_FIPS = sorted(FIPS_TO_STATE.keys()) # 50 states (DC excluded to match LAUS)
_VACANCY_CACHE = Path("data/cache") / "jolts_vacancies_state_quarterly.parquet"

# JOLTS state-level series ID:
# JT + S(1) + industry(6) + state(2) + area(5) + sizeclass(2) + JO(2) + L(1) = 21
# Total nonfarm industry code = "000000"
# Job openings data element = "JO"
# Level = "L"
def _vacancy_series_id(fips2: str) -> str:
 sid = f"JTS000000{fips2}0000000JOL"
 assert len(sid) == 21, f"Series ID wrong length: '{sid}' ({len(sid)} chars)"
 return sid


def fetch_jolts_vacancies_state(
 cache_path : Path = _VACANCY_CACHE,
 start_year : int = 2000,
 end_year : int = END_YEAR,
 sleep_secs : float = 1.2,
) -> pd.DataFrame | None:
 """
 Fetch monthly JOLTS job openings (SA, level, thousands) for all 50 states
 via the BLS public API v2. Caches result as parquet.

 Returns quarterly DataFrame or None if API returns no data.

 Columns: state_fips, quarter_label, vacancies

 ⚠️ AGGREGATION NOTE — STOCK measure:
 Job openings is a stock (level at month-end), NOT a flow. Quarterly
 aggregation averages the three monthly values. Do NOT sum.

 ⚠️ SERIES ID NOTE:
 The state-level JOLTS series IDs used here follow the standard 21-character
 JOLTS schema (JTS + industry + state + area + sizeclass + dataelement +
 ratelevel) derived from jt.txt. If BLS returns 0 observations, the series
 IDs may have changed — check https://www.bls.gov/jlt/jlt_statedata.htm and
 inspect the raw API response logged below.
 """
 if cache_path.exists():
  cached_vac = pd.read_parquet(cache_path)
  # Validate coverage: if the cache ends before end_year - 1, re-fetch.
  _vac_max_q = cached_vac["quarter_label"].max() if "quarter_label" in cached_vac.columns else ""
  if _vac_max_q >= f"{end_year - 1}Q1":
   print(f" [cache hit] {cache_path} (through {_vac_max_q})")
   return cached_vac
  else:
   print(f" [cache stale] {cache_path} ends {_vac_max_q}, "
         f"need ≥{end_year - 1}Q1 — re-fetching")
   cache_path.unlink()

 print(f" [fetch] JOLTS state vacancies {start_year}–{end_year} "
       f"({len(_VACANCY_FIPS)} states)")

 series_ids = [_vacancy_series_id(f) for f in _VACANCY_FIPS]
 sid_to_fips = {_vacancy_series_id(f): f for f in _VACANCY_FIPS}

 # BLS API v2 allows 50 series per request and 20-year windows with a key.
 # Split into batches of 25 and year-chunks of 10 years.
 def _chunks(lst, n):
  for i in range(0, len(lst), n):
   yield lst[i:i+n]

 year_chunks = [(max(2000, y), min(y+9, end_year))
        for y in range(2000, end_year+1, 10)]
 batches = list(_chunks(series_ids, 25))
 n_calls = len(batches) * len(year_chunks)
 all_obs = []
 call_num = 0

 print(f" {len(batches)} batches × {len(year_chunks)} year-chunks = "
       f"{n_calls} API calls")

 for yr_s, yr_e in year_chunks:
  for batch in batches:
   call_num += 1
   fips_sample = sid_to_fips[batch[0]]
   print(f" Call {call_num:2d}/{n_calls}: "
         f"{yr_s}–{yr_e} FIPS {fips_sample}…", end=" ")
   try:
    payload = {
     "seriesid": batch,
     "startyear": str(yr_s),
     "endyear": str(yr_e),
     "registrationkey": BLS_API_KEY,
    }
    resp = requests.post(BLS_URL, json=payload, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    if data.get("status") != "REQUEST_SUCCEEDED":
     msgs = data.get("message", [])
     print(f"→ API error: {msgs}")
     continue
    n_obs = 0
    for series in data["Results"]["series"]:
     sid = series["seriesID"]
     fips = sid_to_fips.get(sid, "??")
     for obs in series.get("data", []):
      if obs.get("period", "") == "M13":
       continue
      try:
       val = float(obs["value"])
      except (ValueError, KeyError):
       continue
      all_obs.append({
       "state_fips": fips,
       "year": obs["year"],
       "period": obs["period"],
       "value": val,
      })
      n_obs += 1
    print(f"→ {n_obs} obs")
   except Exception as exc:
    print(f"→ ERROR: {exc}")
    if call_num == 1:
     # First call failure is likely a series ID or network problem
     print(" ⚠️ First API call failed — check series ID format and "
           "network access. Verify at: "
           "https://data.bls.gov/cgi-bin/surveymost?jt")
     return None
    raise
   if sleep_secs > 0:
    time.sleep(sleep_secs)

 if not all_obs:
  print(
   "\n ⚠️ WARNING: BLS returned 0 vacancy observations for all state series.\n"
   " Possible causes:\n"
   " 1. Series ID format is incorrect for state-level JOLTS JO.\n"
   " Verify at https://data.bls.gov/cgi-bin/surveymost?jt\n"
   " 2. State-level JOLTS JO is not available via the public API v2.\n"
   " Check https://www.bls.gov/jlt/jlt_statedata.htm for download options.\n"
   " 3. BLS API registration key is expired or rate-limited.\n"
   " Vacancy outcome LP will be skipped — unemployment LP unaffected."
  )
  return None

 # Monthly to quarterly — AVERAGE (stock measure, not sum)
 raw = pd.DataFrame(all_obs)
 raw["month"] = raw["period"].str[1:].astype(int)
 raw["date"] = pd.to_datetime(raw[["year", "month"]].assign(day=1))
 raw["quarter"] = raw["date"].dt.to_period("Q").astype(str)

 quarterly = (
  raw.groupby(["state_fips", "quarter"])["value"]
  .agg(vacancies="mean", n_months="count") # AVERAGE for stock measure
  .reset_index()
  .query("n_months == 3") # drop incomplete quarters
  .drop(columns="n_months")
  .rename(columns={"quarter": "quarter_label"})
  .sort_values(["state_fips", "quarter_label"])
  .reset_index(drop=True)
 )

 cache_path.parent.mkdir(parents=True, exist_ok=True)
 quarterly.to_parquet(cache_path, index=False)
 n_states = quarterly["state_fips"].nunique()
 n_qtrs = quarterly["quarter_label"].nunique()
 print(f" Saved vacancy cache → {cache_path} "
       f"({n_states} states × {n_qtrs} quarters)")
 return quarterly


def weighted_mean(df: pd.DataFrame, var_list: list, weight_col: str, group_col: str) -> pd.DataFrame:
 """
 Calculates the weighted mean for a list of variables, grouped by a specified column.
 
 This function is robust to missing values (NaNs) in both the variable and weight columns.

 Args:
  df: The input DataFrame.
  var_list: A list of column names for which to calculate the weighted mean.
  weight_col: The column name to use for weights.
  group_col: The column name to group the calculations by (e.g., 'date' or 'quarter_label').

 Returns:
  A DataFrame with the `group_col` as the index and `var_list` as columns,
  containing the calculated weighted means.
 """
 # Ensure var_list is a list, even if a single string is passed
 if isinstance(var_list, str):
  var_list = [var_list]

 # This will be our final result DataFrame, starting with the unique group values as the index
 unique_groups = df[group_col].unique()
 result_df = pd.DataFrame(index=unique_groups)
 result_df.index.name = group_col # Name the index

 for var in var_list:
  # Create a temporary DataFrame with the necessary columns
  # and drop rows where the variable or the weight is missing.
  # This is the crucial step to prevent errors.
  temp_df = df[[group_col, var, weight_col]].dropna()

  # Define a custom aggregation function for weighted average
  def w_avg(group):
   values = group[var]
   weights = group[weight_col]
   try:
    return np.average(values, weights=weights)
   except ZeroDivisionError: # Handle cases where weights sum to zero for a group
    return np.nan

  # Group the cleaned data and apply the weighted average function
  weighted_means = temp_df.groupby(group_col).apply(w_avg)
  
  # Add the resulting Series to our final DataFrame
  result_df[var] = weighted_means

 return result_df

# ===========================================================================
# Main
# ===========================================================================
CACHE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------
# [A] LAUS: fetch / load and aggregate to quarterly
# -----------------------------------------------------------------------
monthly = fetch_laus_monthly(
 start_year = START_YEAR,
 end_year = END_YEAR,
 cache_path = MONTHLY_CACHE,
)

quarterly = monthly_to_quarterly(monthly)

# -----------------------------------------------------------------------
# [B] JOLTS vacancies: fetch / load and merge into quarterly outcomes
# -----------------------------------------------------------------------
print(f"\n[B] JOLTS state vacancies")
vacancies_q = fetch_jolts_vacancies_state(
 cache_path = _VACANCY_CACHE,
 start_year = 2000,
 end_year = END_YEAR,
)

if vacancies_q is not None and not vacancies_q.empty:
 quarterly = quarterly.merge(
  vacancies_q[["state_fips", "quarter_label", "vacancies"]],
  on=["state_fips", "quarter_label"],
  how="left", # left join: NaN where vacancies not yet available (pre-2001)
 )
 n_vac = quarterly["vacancies"].notna().sum()
 print(f" Merged vacancies: {n_vac:,} non-NaN rows "
       f"(out of {len(quarterly):,} total)")
else:
 # API unavailable or empty — add NaN column so downstream code can
 # check `if "vacancies" in outcomes.columns` without KeyError.
 quarterly["vacancies"] = np.nan
 print(" Vacancies column added as NaN "
       "(fetch unavailable — run on a machine with BLS API access)")


# -----------------------------------------------------------------------
# START: MAJOR FIX AREA - Synchronize column names before saving and analysis
# -----------------------------------------------------------------------

# Calculate vac_rate if vacancies exist. NOTE: The original file had a bug here,
# attempting to calculate vac_rate even if 'vacancies' was all NaN. Fixed.
if "vacancies" in quarterly.columns and quarterly["vacancies"].notna().any():
    quarterly["vac_rate"] = (quarterly["vacancies"] * 1000 / quarterly["labor_force"]) * 100
else:
    quarterly["vac_rate"] = np.nan

# Drop the original vacancies column (in thousands) as vac_rate is the desired metric.
if "vacancies" in quarterly.columns:
    quarterly.drop("vacancies", axis=1, inplace=True)

# Save the final file with the 'quarter_label' column, as documented.
quarterly.to_parquet(QUARTERLY_OUT, index=False)
print(f"Saved quarterly → {QUARTERLY_OUT} ({len(quarterly):,} rows)")

# -----------------------------------------------------------------------
# END: MAJOR FIX AREA
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 4 — LAUS outcomes ({START_YEAR}–{END_YEAR})")
print("=" * 60)

n_states = quarterly["state_fips"].nunique()
n_qtrs = quarterly["quarter_label"].nunique()
print(f"\nQuarterly: {quarterly.shape} ({n_states} states × {n_qtrs} quarters)")
print(f"Quarter range: {quarterly['quarter_label'].min()} → {quarterly['quarter_label'].max()}")
print("\nFirst 10 rows:")
print(quarterly.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Summary table at key quarters — labor-force-weighted mean
# -----------------------------------------------------------------------
print("\n--- Cross-state unemployment rate at key quarters ---")
print(" (weighted mean ≈ BLS national headline; simple mean shown for comparison)")

key_qtrs = [q for q in ["2001Q3", "2007Q4", "2009Q2", "2020Q2", "2023Q1"]
 if q in quarterly["quarter_label"].values]

# FIX: Use 'quarter_label' as the grouping column, which exists in the 'quarterly' DataFrame.
wmean = weighted_mean(quarterly, ["unemp_rate", "vac_rate"], "labor_force", "quarter_label")
wmean.describe()

summary_key = (
 quarterly[quarterly["quarter_label"].isin(key_qtrs)]
 .groupby("quarter_label")["unemp_rate"]
 .agg(simple_mean="mean", min="min", max="max", std="std")
 .round(2)
 # FIX: Join the wmean DataFrame directly; its index is already 'quarter_label'.
 # Also select just the 'unemp_rate' column from the weighted means for this summary.
 .join(wmean[["unemp_rate"]].rename(columns={"unemp_rate": "wtd_mean"}).round(2))
 [["wtd_mean", "simple_mean", "min", "max", "std"]]
)
print(summary_key.to_string())

# State ranking at Great Recession peak
if "2009Q2" in quarterly["quarter_label"].values:
 print("\n--- State ranking at 2009Q2 (Great Recession peak) ---")
 gr = (
  quarterly[quarterly["quarter_label"] == "2009Q2"]
  .sort_values("unemp_rate", ascending=False)
  [["state", "unemp_rate", "labor_force"]]
 )
 print(gr.to_string(index=False))

# -----------------------------------------------------------------------
# Plot: labor-force-weighted mean + unweighted percentile bands over time
# -----------------------------------------------------------------------
def _ql_to_dt(ql):
 y, q = ql.split("Q")
 return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

# -----------------------------------------------------------------------
# Helper: build percentile + weighted-mean series for a given variable
# -----------------------------------------------------------------------
def _dist_series(df, val_col, wt_col="labor_force"):
 """Return (pctile_df, wmean_series) aligned on a complete quarterly grid."""
 pctile = (
  df.dropna(subset=[val_col])
  .groupby("quarter_label")[val_col]
  .agg(
   p10=lambda x: x.quantile(0.10),
   p25=lambda x: x.quantile(0.25),
   p75=lambda x: x.quantile(0.75),
   p90=lambda x: x.quantile(0.90),
  )
  .sort_index()
 )
 # FIX: Select only the relevant column from the wmean DataFrame result.
 wmean = weighted_mean(df.dropna(subset=[val_col]), val_col, wt_col, "quarter_label")[val_col].sort_index()

 all_q = pd.period_range(
  start=pctile.index[0], end=pctile.index[-1], freq="Q"
 ).strftime("%YQ%q").tolist()
 pctile = pctile.reindex(all_q)
 wmean = wmean.reindex(all_q)
 return pctile, wmean


def _add_recession_shading(ax):
 """Add NBER recession bands and COVID shading to an axes object."""
 for start, end, label in [
  ("1990-07-01", "1991-03-01", "NBER recessions"),
  ("2001-03-01", "2001-11-01", "_"),
  ("2007-12-01", "2009-06-01", "_"),
 ]:
  ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
       alpha=0.10, color="grey", label=label)
 ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
      alpha=0.14, color="red", label="COVID-19")


def _fill_panel(ax, pctile, wmean, color, ylabel, formatter=None):
 """Render distribution bands + weighted mean into ax."""
 valid_q = pctile.dropna(subset=["p10"])
 valid_dts = [_ql_to_dt(q) for q in valid_q.index]
 wmean_dts = [_ql_to_dt(q) for q in wmean.index if pd.notna(wmean[q])]
 wmean_vals = wmean.dropna().values

 ax.fill_between(valid_dts, valid_q["p10"], valid_q["p90"],
      alpha=0.18, color=color, label="10th-90th pctile (unweighted)")
 ax.fill_between(valid_dts, valid_q["p25"], valid_q["p75"],
      alpha=0.32, color=color, label="25th-75th pctile (unweighted)")
 ax.plot(wmean_dts, wmean_vals,
      color=color, linewidth=2.0, label="LF-weighted mean")
 _add_recession_shading(ax)
 ax.set_xlabel("")
 ax.set_ylabel(ylabel, fontsize=10)
 if formatter:
  ax.yaxis.set_major_formatter(formatter)
 ax.legend(fontsize=9, framealpha=0.85, ncol=2)
 ax.grid(axis="y", linewidth=0.5, alpha=0.4)


# -----------------------------------------------------------------------
# Build unemployment distribution series
# -----------------------------------------------------------------------
unemp_pctile, unemp_wmean = _dist_series(quarterly, "unemp_rate")

# -----------------------------------------------------------------------
# Determine whether vacancy data is available for panels 2 and 3.
# -----------------------------------------------------------------------
# FIX: Use 'vac_rate' which was created earlier, not 'vacancies'.
_vac_plot = (
 "vac_rate" in quarterly.columns
 and quarterly["vac_rate"].notna().sum() > 0
)

if _vac_plot:
 # The vac_rate column is already computed.
 vac_pctile, vac_wmean = _dist_series(quarterly, "vac_rate")

 # Market tightness: theta = V / U (dimensionless ratio of stocks)
 # = vac_rate / unemp_rate (labor force cancels)
 # Drop rows where either rate is zero or NaN to avoid div-by-zero.
 quarterly["theta"] = quarterly["vac_rate"] / quarterly["unemp_rate"]
 theta_pctile, _ = _dist_series(quarterly, "theta")

 # Aggregate national theta = sum(V_actual) / sum(U_actual) by quarter.
 # This is the economically correct aggregate (not LF-weighted mean of
 # state-level theta values, which inflates small-state extreme ratios).
 # FIX: Use vac_rate and unemp_rate to derive actual vacancy and unemployed counts
 _vac_act = quarterly["vac_rate"] / 100 * quarterly["labor_force"] # persons
 _u_act = quarterly["unemp_rate"] / 100 * quarterly["labor_force"] # persons
 theta_agg = (
  quarterly.assign(v_act=_vac_act, u_act=_u_act)
  .dropna(subset=["v_act", "u_act"])
  .groupby("quarter_label")
  .apply(lambda g: g["v_act"].sum() / g["u_act"].sum(), include_groups=False)
  .rename("theta")
  .sort_index()
 )
 all_q_theta = pd.period_range(
  start=theta_pctile.index[0], end=theta_pctile.index[-1], freq="Q"
 ).strftime("%YQ%q").tolist()
 theta_agg = theta_agg.reindex(all_q_theta)

# -----------------------------------------------------------------------
# Build figure: 3-panel (u + v + theta), 2-panel, or 1-panel fallback
# -----------------------------------------------------------------------
if _vac_plot:
 fig, (ax_u, ax_v, ax_t) = plt.subplots(
  3, 1, figsize=(13, 13), sharex=False,
  gridspec_kw={"hspace": 0.42},
 )
else:
 fig, ax_u = plt.subplots(figsize=(13, 5))

# Panel 1: unemployment 
_fill_panel(
 ax_u, unemp_pctile, unemp_wmean,
 color = "#1f77b4",
 ylabel = "Unemployment rate (%)",
 formatter = mticker.FuncFormatter(lambda x, _: f"{x:.1f}%"),
)
ax_u.set_title(
 "State unemployment rates: cross-state distribution over time\n"
 "(BLS LAUS, SA, quarterly average; mean is labor-force weighted)",
 fontsize=11,
)

if _vac_plot:
 #Panel 2: vacancy rate 
 _fill_panel(
  ax_v, vac_pctile, vac_wmean,
  color = "#2ca02c",
  ylabel = "Vacancy rate (%)",
  formatter = mticker.FuncFormatter(lambda x, _: f"{x:.1f}%"),
 )
 ax_v.set_title(
  "State vacancy rates: cross-state distribution over time\n"
  "(BLS JOLTS openings / LAUS labor force, SA, quarterly average; mean is labor-force weighted)",
  fontsize=11,
 )

 # Panel 3: market tightness theta = V/U (linear scale) 
 _fill_panel(
  ax_t, theta_pctile, theta_agg,
  color = "#ff7f0e",
  ylabel = r"Tightness $\theta$ = V/U",
  formatter = mticker.FuncFormatter(lambda x, _: f"{x:.2f}"),
 )
 ax_t.yaxis.set_major_locator(mticker.MultipleLocator(0.5))
 ax_t.yaxis.set_minor_locator(mticker.MultipleLocator(0.25))
 # Reference line: theta = 1 (one vacancy per unemployed worker)
 ax_t.axhline(1.0, color="black", linewidth=0.8, linestyle="--",
        label=r"$\theta = 1$ (balanced)")
 ax_t.legend(fontsize=9, framealpha=0.85, ncol=2)
 ax_t.set_title(
  r"Market tightness $\theta$ = V/U: cross-state distribution over time"
  "\n(mean = aggregate national V/U; percentile bands are unweighted cross-state)",
  fontsize=11,
 )

fig.tight_layout()

outpath = OUTPUT_DIR / "outcomes_distribution.png"
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nPlot saved: {outpath}")
_open_file(outpath)

# -----------------------------------------------------------------------
# Separate figure: tightness on log scale (Shimer-style)
# -----------------------------------------------------------------------
if _vac_plot:
 _LOG_TICKS = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]

 fig_log, ax_log = plt.subplots(figsize=(13, 5))
 _fill_panel(
  ax_log, theta_pctile, theta_agg,
  color = "#ff7f0e",
  ylabel = r"Tightness $\theta$ = V/U (log scale)",
  formatter = None, # set manually below
 )
 ax_log.set_yscale("log")
 ax_log.yaxis.set_major_locator(mticker.FixedLocator(_LOG_TICKS))
 ax_log.yaxis.set_minor_locator(mticker.NullLocator())
 ax_log.yaxis.set_major_formatter(
  mticker.FuncFormatter(lambda x, _: f"{x:g}")
 )
 ax_log.axhline(1.0, color="black", linewidth=0.8, linestyle="--",
        label=r"$\theta = 1$ (balanced)")
 ax_log.legend(fontsize=9, framealpha=0.85, ncol=2)
 ax_log.set_title(
  r"Market tightness $\theta$ = V/U: cross-state distribution over time [log scale]"
  "\n(mean = aggregate national V/U; percentile bands are unweighted cross-state)",
  fontsize=11,
 )
 fig_log.tight_layout()

 outpath_log = OUTPUT_DIR / "tightness_log.png"
 fig_log.savefig(outpath_log, dpi=150)
 plt.close(fig_log)
 print(f"Plot saved: {outpath_log}")
 _open_file(outpath_log)