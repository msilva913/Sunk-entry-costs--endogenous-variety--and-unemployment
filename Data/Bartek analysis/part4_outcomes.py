"""
part4_outcomes.py — State-level unemployment rates and labor force from BLS LAUS
=================================================================================
Fetches monthly seasonally adjusted unemployment rates (suffix "3") and labor
force levels (suffix "6") for all 50 states from the BLS Local Area
Unemployment Statistics (LAUS) program via the BLS public API v1.  Both series
are aggregated to quarterly frequency and saved to a single parquet.

LAUS series ID format
---------------------
    LASST{state_fips_2d}0000000000003   — unemployment rate (%)
    LASST{state_fips_2d}0000000000006   — civilian labor force (thousands)

Labor force is fetched alongside the unemployment rate for three reasons:
  1. Weighted means: simple cross-state means overweight small states.
     Labor-force-weighted means approximately recover the BLS national headline.
  2. WLS robustness: LP regressions can be rerun with labor-force weights to
     give larger states more influence — standard robustness check (ADH 2013).
  3. Level conversion: converting LP percentage-point responses to worker counts
     requires a base labor force: ΔU_{s,t} = Δu_{s,t} × LF_{s,t0}.

Sample window
-------------
    1990M1 – present.  Provides a buffer before the 2001Q1 joint instrument
    sample start and covers all post-cold-war cycles (1990–91, 2001, GR, COVID).

API batching
------------
BLS API v1 (unregistered): 25 series per request, 10 years per request.
We fetch rate and labor-force series in the same batch by interleaving both
suffixes — 50 states × 2 series = 100 series total, split into 4 batches of
25 per year-chunk.  Year-chunks of 8 years → ~5 chunks × 4 batches = 20 calls.
A 1-second sleep is inserted between calls.  Results cached as parquet.

Outputs
-------
    data/cache/laus_states_monthly.parquet    — monthly, 50 states
    data/instruments/laus_quarterly.parquet   — quarterly averages

Columns in both parquets:
    state_fips      zero-padded 2-digit string, e.g. "01"
    state           2-letter abbreviation, e.g. "AL"
    date / quarter_label
    unemp_rate      unemployment rate (%)
    labor_force     civilian labor force (thousands)

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
START_YEAR = 1990
END_YEAR   = 2024   # update as needed

CACHE_DIR  = Path("data/cache")
OUTPUT_DIR = Path("data/instruments")

MONTHLY_CACHE = CACHE_DIR  / "laus_states_monthly.parquet"
QUARTERLY_OUT = OUTPUT_DIR / "laus_quarterly.parquet"

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

# LAUS suffix codes
_RATE_SUFFIX = "3"   # unemployment rate (%)
_LF_SUFFIX   = "6"   # civilian labor force (thousands)


def _laus_id(fips: str, suffix: str) -> str:
    return f"LASST{fips}0000000000{suffix}"


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


def fetch_laus_monthly(
    start_year: int = START_YEAR,
    end_year:   int = END_YEAR,
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
        # Validate cache has both columns — a partial fetch (e.g. after hitting
        # the API request limit mid-run) would have cached only unemp_rate.
        # Reject and re-fetch if labor_force is missing.
        if "labor_force" in cached.columns:
            print(f"  [cache hit] {cache_path}")
            return cached
        else:
            print(f"  [cache invalid] {cache_path} is missing labor_force — re-fetching")
            cache_path.unlink()

    print(f"  [fetch] LAUS rate + labor force  {start_year}–{end_year}  (50 states)")

    # Build series lists: interleave rate and LF so each batch of 25 covers
    # one state-half and both series types.
    half_a = ALL_FIPS[:25]
    half_b = ALL_FIPS[25:]
    # Each sub-batch: 25 series = 12 or 13 states × 2 suffixes (truncated to 25)
    # Split each half into rate-batch and LF-batch
    batches = [
        [_laus_id(f, _RATE_SUFFIX) for f in half_a],   # 25 rate series, states 1-25
        [_laus_id(f, _LF_SUFFIX)   for f in half_a],   # 25 LF series, states 1-25
        [_laus_id(f, _RATE_SUFFIX) for f in half_b],   # 25 rate series, states 26-50
        [_laus_id(f, _LF_SUFFIX)   for f in half_b],   # 25 LF series, states 26-50
    ]
    year_chunks = _year_chunks(start_year, end_year)
    print(f"  {len(batches)} batches × {len(year_chunks)} year-chunks "
          f"= {len(batches)*len(year_chunks)} API calls")

    all_frames = []
    call_num   = 0
    for yr_start, yr_end in year_chunks:
        for batch in batches:
            call_num += 1
            suffix_label = "rate" if batch[0].endswith(_RATE_SUFFIX) else "LF"
            fips_sample  = batch[0][5:7]
            print(f"    Call {call_num:2d}: {yr_start}–{yr_end}  "
                  f"{suffix_label}  FIPS starts {fips_sample}…", end=" ")
            try:
                raw = _fetch_batch(batch, yr_start, yr_end)
                df  = _parse_series(raw)
                all_frames.append(df)
                print(f"→ {len(df)} rows")
            except Exception as exc:
                print(f"→ ERROR: {exc}")
                raise
            time.sleep(SLEEP_SEC)

    raw_df = pd.concat(all_frames, ignore_index=True)

    # Extract FIPS (chars 5–6 of series ID: LASST{XX}...)
    raw_df["state_fips"] = raw_df["series_id"].str[5:7]
    # Identify series type from trailing digit
    raw_df["series_type"] = raw_df["series_id"].str[-1]
    raw_df["month"] = raw_df["period"].str[1:].astype(int)
    raw_df["date"]  = pd.to_datetime(raw_df[["year", "month"]].assign(day=1))
    raw_df["state"] = raw_df["state_fips"].map(FIPS_TO_STATE)

    # Pivot so rate and labor_force are separate columns
    pivot = (
        raw_df
        .pivot_table(
            index   = ["state_fips", "state", "date"],
            columns = "series_type",
            values  = "value",
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
    print(f"  Cached → {cache_path}  ({len(out):,} rows)")
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
            unemp_rate  = ("unemp_rate",  "mean"),
            labor_force = ("labor_force", "mean"),
            n_months    = ("unemp_rate",  "count"),   # completeness check
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


def weighted_mean(df: pd.DataFrame, val_col: str, wt_col: str) -> pd.Series:
    """
    Labor-force-weighted mean of val_col by quarter_label.

    Used in place of simple means throughout: simple means overweight
    small states (Wyoming = California).  Weighted means approximately
    recover the BLS published national headline rate.
    """
    return (
        df.groupby("quarter_label")
        .apply(lambda g: np.average(g[val_col], weights=g[wt_col]))
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
monthly = fetch_laus_monthly(
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
print(f"Part 4 — LAUS outcomes  ({START_YEAR}–{END_YEAR})")
print("=" * 60)

n_states = quarterly["state_fips"].nunique()
n_qtrs   = quarterly["quarter_label"].nunique()
print(f"\nQuarterly: {quarterly.shape}  ({n_states} states × {n_qtrs} quarters)")
print(f"Quarter range: {quarterly['quarter_label'].min()} → {quarterly['quarter_label'].max()}")
print("\nFirst 10 rows:")
print(quarterly.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Summary table at key quarters — labor-force-weighted mean
# Simple cross-state mean is also shown for reference.
# -----------------------------------------------------------------------
print("\n--- Cross-state unemployment rate at key quarters ---")
print("    (weighted mean ≈ BLS national headline; simple mean shown for comparison)")

key_qtrs = [q for q in ["2001Q3", "2007Q4", "2009Q2", "2020Q2", "2023Q1"]
            if q in quarterly["quarter_label"].values]

wmean = weighted_mean(quarterly, "unemp_rate", "labor_force")

summary_key = (
    quarterly[quarterly["quarter_label"].isin(key_qtrs)]
    .groupby("quarter_label")["unemp_rate"]
    .agg(simple_mean="mean", min="min", max="max", std="std")
    .round(2)
    .join(wmean.rename("wtd_mean").round(2))
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
#
# Mean line: labor-force-weighted — recovers BLS national headline.
# Percentile bands: unweighted — describe the cross-state distribution
#   as seen by the LP regression (each state = one observation).
# -----------------------------------------------------------------------
def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

# Unweighted percentile bands (cross-state distribution for LP context)
pctile = (
    quarterly
    .groupby("quarter_label")["unemp_rate"]
    .agg(
        p10 = lambda x: x.quantile(0.10),
        p25 = lambda x: x.quantile(0.25),
        p75 = lambda x: x.quantile(0.75),
        p90 = lambda x: x.quantile(0.90),
    )
    .sort_index()
)

# Labor-force-weighted mean (national aggregate for the mean line)
wmean_q = weighted_mean(quarterly, "unemp_rate", "labor_force").sort_index()

# Align on complete quarterly grid
all_q   = pd.period_range(
    start=pctile.index[0], end=pctile.index[-1], freq="Q"
).strftime("%YQ%q").tolist()
pctile  = pctile.reindex(all_q)
wmean_q = wmean_q.reindex(all_q)

valid_q   = pctile.dropna(subset=["p10"])
valid_dts = [_ql_to_dt(q) for q in valid_q.index]
wmean_dts = [_ql_to_dt(q) for q in wmean_q.index if pd.notna(wmean_q[q])]
wmean_vals = wmean_q.dropna().values

fig, ax = plt.subplots(figsize=(13, 5))

# Unweighted percentile bands — reflect LP cross-section
ax.fill_between(valid_dts, valid_q["p10"], valid_q["p90"],
                alpha=0.18, color="#1f77b4", label="10–90th pctile (unweighted)")
ax.fill_between(valid_dts, valid_q["p25"], valid_q["p75"],
                alpha=0.32, color="#1f77b4", label="25–75th pctile (unweighted)")

# Labor-force-weighted mean — recovers BLS national headline
ax.plot(wmean_dts, wmean_vals,
        color="#1f77b4", linewidth=2.0, label="LF-weighted mean (≈ BLS national)")

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
    "State unemployment rates: cross-state distribution over time\n"
    "(BLS LAUS, SA, quarterly average; mean is labor-force weighted)",
    fontsize=12,
)
ax.set_xlabel("")
ax.set_ylabel("Unemployment rate (%)", fontsize=10)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.1f}%"))
ax.legend(fontsize=9, framealpha=0.85, ncol=2)
ax.grid(axis="y", linewidth=0.5, alpha=0.4)
fig.tight_layout()

outpath = OUTPUT_DIR / "unemp_rate_distribution.png"
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nPlot saved: {outpath}")
_open_file(outpath)