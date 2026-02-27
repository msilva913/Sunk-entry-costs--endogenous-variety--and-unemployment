"""
part2_shock_rates_s.py — National LOO layoff rates  (g^s_{-s,j,t})
====================================================================
Loads the employment shares saved by part1_shares.py, downloads the
JOLTS Layoffs & Discharges flat file (national + state), and computes
the leave-one-out shock rate for every (state, supersector, quarter)
cell:

    g^s_{-s,j,t} = layoffs^nat_{-s,j,t} / E^nat_{-s,j,t0}

LOO correction applied as follows (option a):
    Full LOO    : where state JOLTS is published (numerator and denominator)
    Denom-only  : where state JOLTS is suppressed (small-state / small-industry)

The loo_type column records which correction was applied for each cell.

Output
------
    data/instruments/shock_rates_s_{START}_{END}.parquet

Columns in output:
    state_fips     2-digit state FIPS
    industry_code  2-digit supersector code (pipeline convention)
    quarter_label  e.g. "2008Q4"
    g_s_loo        LOO layoff rate (layoffs / base employment)
    loo_type       "full" or "denom_only"

Run
---
    python part2_shock_rates_s.py

Prerequisite
------------
    python part1_shares.py   (must have run first)

Sample window
-------------
    2003Q1 – 2024Q2  (first quarter with state JOLTS supersector coverage)
"""

import sys
import pandas as pd
from construct_s_instrument import (
    BASE_YEAR,
    START_QUARTER,
    END_QUARTER,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_S_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_national_shock_rates_s,
    fetch_jolts_layoffs,
)

# -----------------------------------------------------------------------
# Load Part 1 output  (shares are shared with the delta instrument)
# -----------------------------------------------------------------------
if not SHARES_PATH.exists():
    sys.exit(
        f"ERROR: {SHARES_PATH} not found.\n"
        f"       Run part1_shares.py first."
    )
shares = pd.read_parquet(SHARES_PATH)
print(f"Loaded shares from {SHARES_PATH}  ({len(shares):,} rows)")

# -----------------------------------------------------------------------
# Compute shock rates
# -----------------------------------------------------------------------
shock_rates = build_national_shock_rates_s(
    shares        = shares,
    start_quarter = START_QUARTER,
    end_quarter   = END_QUARTER,
    cache_dir     = DEFAULT_CACHE_DIR,
    save_output   = True,
    output_dir    = DEFAULT_OUTPUT_DIR,
)

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 2s — LOO s-shock rates  ({START_QUARTER} – {END_QUARTER})")
print("=" * 60)

n_states = shock_rates["state_fips"].nunique()
n_ss     = shock_rates["industry_code"].nunique()
n_qtrs   = shock_rates["quarter_label"].nunique()
n_valid  = shock_rates["g_s_loo"].notna().sum()
print(f"\nShape: {shock_rates.shape}  "
      f"({n_states} states × {n_ss} supersectors × {n_qtrs} quarters)")
print(f"Valid cells: {n_valid:,} / {len(shock_rates):,}")

print("\nFirst 10 rows:")
print(shock_rates.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# LOO type breakdown by supersector
# -----------------------------------------------------------------------
print("\n--- LOO coverage by supersector ---")
print("    'full' = both numerator and denominator corrected;")
print("    'denom_only' = state layoffs suppressed, denominator-only LOO\n")
loo_counts = (
    shock_rates
    .groupby(["industry_code", "loo_type"])
    .size()
    .unstack("loo_type", fill_value=0)
    .rename(index=INDUSTRY_LABELS)
)
# Add pct_full column
for col in ["full", "denom_only"]:
    if col not in loo_counts.columns:
        loo_counts[col] = 0
loo_counts["total"]    = loo_counts["full"] + loo_counts["denom_only"]
loo_counts["pct_full"] = (
    loo_counts["full"] / loo_counts["total"] * 100
).round(1)
print(loo_counts.to_string())

# -----------------------------------------------------------------------
# Mean LOO rate by supersector
# -----------------------------------------------------------------------
print("\n--- Mean LOO layoff rate by supersector (×1 000) ---")
mean_rate = (
    shock_rates.groupby("industry_code")["g_s_loo"]
    .mean()
    .rename(index=INDUSTRY_LABELS)
    .sort_values(ascending=False)
    .mul(1000)
    .round(3)
)
for name, rate in mean_rate.items():
    bar = "█" * int(rate * 6)
    print(f"  {name:<36}  {rate:6.3f}  {bar}")

# -----------------------------------------------------------------------
# Time series: cross-state mean rate by supersector (sample period)
# -----------------------------------------------------------------------
print("\n--- Time series: cross-state mean rate by supersector (×1 000) ---")
print("    COVID surge and Great Recession should be visible.\n")
pivot = (
    shock_rates
    .groupby(["quarter_label", "industry_code"])["g_s_loo"]
    .mean()
    .unstack("industry_code")
    .rename(columns=INDUSTRY_LABELS)
    .mul(1000)
    .round(3)
)
pivot.columns = [c[:12] for c in pivot.columns]
with pd.option_context("display.max_columns", 12, "display.width", 130,
                       "display.float_format", "{:6.3f}".format):
    print(pivot.loc["2007Q1":"2011Q4"].to_string())

# -----------------------------------------------------------------------
# Peak layoff quarter by supersector
# -----------------------------------------------------------------------
print("\n--- Peak layoff quarter by supersector ---")
peak_by_ss = (
    shock_rates
    .groupby(["industry_code", "quarter_label"])["g_s_loo"]
    .mean()
    .reset_index()
    .loc[lambda d: d.groupby("industry_code")["g_s_loo"]
         .transform("max") == d["g_s_loo"]]
    .assign(
        supersector = lambda d: d["industry_code"].map(INDUSTRY_LABELS),
        g_x1000     = lambda d: (d["g_s_loo"] * 1000).round(3),
    )
    [["supersector", "quarter_label", "g_x1000"]]
    .sort_values("g_x1000", ascending=False)
)
print(peak_by_ss.to_string(index=False))

# -----------------------------------------------------------------------
# Sanity check: national JOLTS levels (should be 1 000–6 000 thousands)
# -----------------------------------------------------------------------
print("\n--- Raw national JOLTS layoff levels by supersector (sanity check) ---")
print("    All levels in thousands of workers per quarter.")
print("    Plausible ranges: ~500k (Mining) to ~5 000k (Trade/transport)\n")
jolts_nat, _ = fetch_jolts_layoffs(
    cache_dir     = DEFAULT_CACHE_DIR,
    start_quarter = START_QUARTER,
    end_quarter   = END_QUARTER,
)
nat_summary = (
    jolts_nat
    .groupby("industry_code")["layoffs_nat"]
    .agg(["mean", "min", "max", "count"])
    .rename(index=INDUSTRY_LABELS)
    .rename(columns={"mean": "mean_k", "min": "min_k",
                     "max": "max_k", "count": "n_quarters"})
    .sort_values("mean_k", ascending=False)
    .round(0)
)
print(nat_summary.to_string())
