"""
part4b_vacancies.py — State-level JOLTS job openings (vacancies) from BLS JOLTS
================================================================================
Fetches state-level JOLTS job openings (JO, level in thousands) for all 50
states from the BLS JOLTS program via the BLS public API v2.  Series are
aggregated from monthly to quarterly frequency and saved to parquet files.

IMPORTANT DATA AVAILABILITY NOTE
--------------------------------
JOLTS publishes state-level data only for total nonfarm at the state level
(not broken down by industry).  State-level JOLTS series by supersector are
not available via the public API.  This script fetches total nonfarm job
openings (JO) at the state level.

Coverage and history:
    - JOLTS program launched: December 2000
    - State-level data coverage: Varies by state
    - Most states: 2001Q1 onward
    - A few states: Later (see diagnostics after fetch)
    - National series: 2000M12 onward (2001Q1 first complete quarter)

JOLTS series ID format  (21 chars)
---------------------------------
    JT + seasonal(1) + industry(6) + state(2) + area(5) + sizeclass(2)
       + dataelement(2) + ratelevel(1)

For state-level total nonfarm job openings (level, thousands):
    seasonal     = S   (seasonally adjusted)
    industry     = 000000  (total nonfarm)
    state        = SS   (state FIPS 2-digit, e.g. "01" for Alabama)
    area         = 00000  (state total, not metro)
    sizeclass    = 00   (all sizes)
    dataelement  = JO   (job openings)
    ratelevel    = L    (level, thousands of jobs)

Example series IDs:
    JTS00000001000000JOL  — Alabama total nonfarm job openings (level)
    JTS00000006000000JOL  — California total nonfarm job openings (level)

Monthly to quarterly aggregation
---------------------------------
JOLTS publishes monthly levels (thousands of job openings at mid-month).
We average the three months in each calendar quarter (not sum, since these
are stock measures from a specific date each month, not flows).
Incomplete quarters (fewer than 3 months returned) are dropped.

Outputs
-------
    data/cache/jolts_jo_states_monthly.parquet    — monthly, 50 states
    data/instruments/jolts_jo_quarterly.parquet   — quarterly averages

Columns in both parquets:
    state_fips      zero-padded 2-digit string, e.g. "01"
    state           2-letter abbreviation, e.g. "AL"
    date / quarter_label
    job_openings    job openings (level, thousands)

API batching
------------
BLS API v2 (registered): 50 series per request, 20 years per request.
We fetch 50 state total nonfarm JO series in batches of 50, with year-chunks
of 8–10 years.  Year-chunks: 2000–2008 / 2009–2017 / 2018–2030 = 3 chunks ×
1 batch = 3 API calls total.  A 1-second sleep is inserted between calls.
Results cached as parquet.

Run
---
    python part4b_vacancies.py

No prerequisites (does not depend on parts 1–3).
"""

import os
import time
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — avoids crash issues
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

def _open_file(path):
    """Open a saved file with the OS default viewer."""
    os.startfile(path)

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents"
        / "GitHub"
        / "Sunk-entry-costs--endogenous-variety--and-unemployment"
        / "Data"
        / "Bartek analysis"
    )

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START_YEAR = 2000
END_YEAR   = 2024   # update as needed

CACHE_DIR  = Path("data/cache")
OUTPUT_DIR = Path("data/instruments")

MONTHLY_CACHE = CACHE_DIR  / "jolts_jo_states_monthly.parquet"
QUARTERLY_OUT = OUTPUT_DIR / "jolts_jo_quarterly.parquet"

BLS_URL     = "https://api.bls.gov/publicAPI/v2/timeseries/data/"  # v2 required with registration key
SLEEP_SEC   = 1.0
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

ALL_FIPS = sorted(FIPS_TO_STATE.keys())   # 50 entries


def _jolts_jo_id(state_fips: str) -> str:
    """
    Construct JOLTS state-level total nonfarm job openings series ID.

    JTS + 000000 (total nonfarm industry) + SS (state FIPS 2-digit)
       + 00000 (state total area) + 00 (all sizes) + JO (job openings)
       + L (level, thousands)

    Parameters
    ----------
    state_fips : str
        2-digit zero-padded state FIPS code (e.g. "01" for Alabama)

    Returns
    -------
    str
        21-character JOLTS series ID, e.g. "JTS00000001000000JOL"
    """
    sid = f"JTS000000{state_fips}000000JOL"
    assert len(sid) == 21, f"Series ID length error: '{sid}' is {len(sid)} chars"
    return sid


def _year_chunks(start: int, end: int, size: int = 8):
    chunks, y = [], start
    while y <= end:
        chunks.append((y, min(y + size - 1, end)))
        y += size
    return chunks


def _fetch_batch(series_ids: list, start_year: int, end_year: int) -> list:
    payload = {
        "seriesid":        series_ids,
        "startyear":       str(start_year),
        "endyear":         str(end_year),
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
            if obs["period"] == "M13":   # skip annual averages
                continue
            try:
                val = float(obs["value"])
            except (ValueError, KeyError):
                continue
            rows.append({
                "series_id": sid,
                "year":   int(obs["year"]),
                "period": obs["period"],
                "value":  val,
            })
    return pd.DataFrame(rows)


def fetch_jolts_jo_monthly(
    start_year: int = START_YEAR,
    end_year:   int = END_YEAR,
    cache_path: Path = MONTHLY_CACHE,
) -> pd.DataFrame:
    """
    Fetch JOLTS state-level job openings (total nonfarm) for all 50 states.

    JOLTS publishes state-level data only for total nonfarm (not by industry).
    For each state, fetches monthly seasonally-adjusted job openings levels
    (thousands) from the most recent mid-month snapshot.

    Returns DataFrame: state_fips, state, date, job_openings
    """
    if cache_path.exists():
        cached = pd.read_parquet(cache_path)
        # Validate cache has job_openings column
        if "job_openings" in cached.columns:
            print(f"  [cache hit] {cache_path}")
            return cached
        else:
            print(f"  [cache invalid] {cache_path} is missing job_openings — re-fetching")
            cache_path.unlink()

    print(f"  [fetch] JOLTS state-level job openings  {start_year}–{end_year}  (50 states)")

    # Build series list: all 50 states × 1 series = 50 series
    # Fits in one batch (v2 API: 50 series per request)
    series_ids = [_jolts_jo_id(fips) for fips in ALL_FIPS]

    year_chunks = _year_chunks(start_year, end_year)
    print(f"  {len(year_chunks)} year-chunks × 1 batch = {len(year_chunks)} API calls")

    all_frames = []
    call_num   = 0
    for yr_start, yr_end in year_chunks:
        call_num += 1
        print(f"    Call {call_num:2d}: {yr_start}–{yr_end}  "
              f"50 states…", end=" ")
        try:
            raw = _fetch_batch(series_ids, yr_start, yr_end)
            df  = _parse_series(raw)
            all_frames.append(df)
            print(f"→ {len(df)} rows")
        except Exception as exc:
            print(f"→ ERROR: {exc}")
            raise
        time.sleep(SLEEP_SEC)

    raw_df = pd.concat(all_frames, ignore_index=True)

    # Extract FIPS (chars 7–8 of series ID: JTS000000{XX}...)
    raw_df["state_fips"] = raw_df["series_id"].str[7:9]
    raw_df["month"] = raw_df["period"].str[1:].astype(int)
    raw_df["date"]  = pd.to_datetime(raw_df[["year", "month"]].assign(day=1))
    raw_df["state"] = raw_df["state_fips"].map(FIPS_TO_STATE)

    # Rename value column
    raw_df["job_openings"] = raw_df["value"]

    out = (
        raw_df[["state_fips", "state", "date", "job_openings"]]
        .sort_values(["state_fips", "date"])
        .reset_index(drop=True)
    )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(cache_path, index=False)
    print(f"  Cached → {cache_path}  ({len(out):,} rows)")
    return out


def monthly_to_quarterly(monthly: pd.DataFrame) -> pd.DataFrame:
    """
    Average monthly JOLTS job openings to quarterly frequency using resample.

    Job openings are stock measures (snapshot at mid-month), so we average
    rather than sum the three months in each quarter.  Incomplete quarters
    (fewer than 3 months) are dropped.

    Returns columns: state_fips, state, quarter_label, job_openings
    """
    out = (
        monthly
        .set_index("date")
        .groupby(["state_fips", "state"])[["job_openings"]]
        .resample("QE")
        .agg(
            job_openings = ("job_openings",  "mean"),
            n_months     = ("job_openings",  "count"),   # completeness check
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


def weighted_mean(df: pd.DataFrame, val_col: str) -> pd.Series:
    """
    Simple unweighted mean of val_col by quarter_label.

    For now, we compute unweighted means (each state = one observation).
    Future: could weight by labor force for comparison with unemployment.
    """
    return (
        df.groupby("quarter_label")[val_col]
        .mean()
        .rename(val_col)
    )


# ===========================================================================
# Main
# ===========================================================================
CACHE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------
# Fetch / load
# -----------------------------------------------------------------------
monthly = fetch_jolts_jo_monthly(
    start_year = START_YEAR,
    end_year   = END_YEAR,
    cache_path = MONTHLY_CACHE,
)

# -----------------------------------------------------------------------
# Quarterly aggregation
# -----------------------------------------------------------------------
quarterly = monthly_to_quarterly(monthly)
quarterly.to_parquet(QUARTERLY_OUT, index=False)
print(f"Saved quarterly → {QUARTERLY_OUT}  ({len(quarterly):,} rows)")

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 4b — JOLTS job openings  ({START_YEAR}–{END_YEAR})")
print("=" * 60)

n_states = quarterly["state_fips"].nunique()
n_qtrs   = quarterly["quarter_label"].nunique()
print(f"\nQuarterly: {quarterly.shape}  ({n_states} states × {n_qtrs} quarters)")
print(f"Quarter range: {quarterly['quarter_label'].min()} → {quarterly['quarter_label'].max()}")
print("\nFirst 10 rows:")
print(quarterly.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Data coverage by state
# -----------------------------------------------------------------------
print("\n--- Data coverage by state ---")
state_coverage = (
    quarterly
    .groupby("state")
    .agg(
        n_qtrs=("quarter_label", "count"),
        q_start=("quarter_label", "min"),
        q_end=("quarter_label", "max"),
    )
    .sort_values("q_start")
    .reset_index()
)
print(state_coverage.to_string(index=False))

# -----------------------------------------------------------------------
# Summary statistics at key quarters
# -----------------------------------------------------------------------
print("\n--- Cross-state job openings at key quarters ---")

key_qtrs = [q for q in ["2001Q1", "2007Q4", "2009Q2", "2019Q4", "2023Q1"]
            if q in quarterly["quarter_label"].values]

summary_key = (
    quarterly[quarterly["quarter_label"].isin(key_qtrs)]
    .groupby("quarter_label")["job_openings"]
    .agg(mean="mean", min="min", max="max", std="std")
    .round(2)
)
print(summary_key.to_string())

# -----------------------------------------------------------------------
# Plot: cross-state mean and percentile bands over time
# -----------------------------------------------------------------------
def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

# Unweighted percentile bands (cross-state distribution)
pctile = (
    quarterly
    .groupby("quarter_label")["job_openings"]
    .agg(
        p10 = lambda x: x.quantile(0.10),
        p25 = lambda x: x.quantile(0.25),
        p75 = lambda x: x.quantile(0.75),
        p90 = lambda x: x.quantile(0.90),
    )
    .sort_index()
)

# Simple cross-state mean
mean_q = weighted_mean(quarterly, "job_openings").sort_index()

# Align on complete quarterly grid
all_q   = pd.period_range(
    start=pctile.index[0], end=pctile.index[-1], freq="Q"
).strftime("%YQ%q").tolist()
pctile  = pctile.reindex(all_q)
mean_q  = mean_q.reindex(all_q)

valid_q   = pctile.dropna(subset=["p10"])
valid_dts = [_ql_to_dt(q) for q in valid_q.index]
mean_dts  = [_ql_to_dt(q) for q in mean_q.index if pd.notna(mean_q[q])]
mean_vals = mean_q.dropna().values

fig, ax = plt.subplots(figsize=(13, 5))

# Unweighted percentile bands
ax.fill_between(valid_dts, valid_q["p10"], valid_q["p90"],
                alpha=0.18, color="#2ca02c", label="10–90th pctile (unweighted)")
ax.fill_between(valid_dts, valid_q["p25"], valid_q["p75"],
                alpha=0.32, color="#2ca02c", label="25–75th pctile (unweighted)")

# Simple cross-state mean
ax.plot(mean_dts, mean_vals,
        color="#2ca02c", linewidth=2.0, label="Cross-state mean")

# NBER recession shading
for start, end, label in [
    ("1990-07-01", "1991-03-01", "NBER recessions"),
    ("2001-03-01", "2001-11-01", "_"),
    ("2007-12-01", "2009-06-01", "_"),
]:
    ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
               alpha=0.10, color="grey", label=label)
ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
           alpha=0.14, color="red", label="COVID-19")

ax.set_title(
    "State job openings: cross-state distribution over time\n"
    "(BLS JOLTS, SA, quarterly average, thousands)",
    fontsize=12,
)
ax.set_xlabel("")
ax.set_ylabel("Job openings (thousands)", fontsize=10)
ax.legend(fontsize=9, framealpha=0.85, ncol=2)
ax.grid(axis="y", linewidth=0.5, alpha=0.4)
fig.tight_layout()

outpath = OUTPUT_DIR / "jo_distribution.png"
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nPlot saved: {outpath}")
_open_file(outpath)
