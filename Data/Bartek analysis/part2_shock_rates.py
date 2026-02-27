"""
part2_shock_rates.py — National LOO closing rates  (g^delta_{-s,j,t})
=======================================================================
Loads the employment shares saved by part1_shares.py, fetches BED
national closing employment from the BLS public API, and computes the
leave-one-out shock rate for every (state, supersector, quarter) cell:

    g^delta_{-s,j,t} = closings^nat_{j,t} / E^nat_{-s,j,t0}

The numerator is the same for all states; the denominator is state-
specific because E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}.

Output
------
    data/instruments/shock_rates_{START}_{END}.parquet

Columns in output:
    state_fips     2-digit state FIPS
    industry_code  2-digit supersector code
    quarter_label  e.g. "2008Q4"
    g_delta_loo    LOO closing rate (employment lost / base employment)

Run
---
    python part2_shock_rates.py

Prerequisite
------------
    python part1_shares.py   (must have run first)
"""

import sys
import pandas as pd
import os
os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment\Data\Bartek analysis")  # Set working directory to script's location

from construct_delta_instrument import (
    BASE_YEAR,
    START_QUARTER,
    END_QUARTER,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_national_shock_rates,
)

# -----------------------------------------------------------------------
# Load Part 1 output
# -----------------------------------------------------------------------
if not SHARES_PATH.exists():
    sys.exit(
        f"ERROR: {SHARES_PATH} not found.\n"
        f"       Run part1_shares.py first."
    )
shares = pd.read_parquet(SHARES_PATH)
print(f"Loaded shares from {SHARES_PATH}  ({len(shares):,} rows)")

# -----------------------------------------------------------------------
# Compute
# -----------------------------------------------------------------------
shock_rates = build_national_shock_rates(
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
print(f"Part 2 — LOO shock rates  ({START_QUARTER} – {END_QUARTER})")
print("=" * 60)

n_states = shock_rates["state_fips"].nunique()
n_ss     = shock_rates["industry_code"].nunique()
n_qtrs   = shock_rates["quarter_label"].nunique()
n_valid  = shock_rates["g_delta_loo"].notna().sum()
print(f"\nShape: {shock_rates.shape}  "
      f"({n_states} states × {n_ss} supersectors × {n_qtrs} quarters)")
print(f"Valid cells: {n_valid:,} / {len(shock_rates):,}")

print("\nFirst 10 rows:")
print(shock_rates.head(10).to_string(index=False))

print("\n--- Mean LOO closing rate by supersector (×1 000) ---")
mean_rate = (
    shock_rates.groupby("industry_code")["g_delta_loo"]
    .mean()
    .rename(index=INDUSTRY_LABELS)
    .sort_values(ascending=False)
    .mul(1000)
    .round(3)
)
for name, rate in mean_rate.items():
    bar = "█" * int(rate * 40)
    print(f"  {name:<36}  {rate:6.3f}  {bar}")

print("\n--- Time series: cross-state mean rate by supersector (×1 000) ---")
print("    Recession-era rows highlighted by elevated values.\n")
pivot = (
    shock_rates
    .groupby(["quarter_label", "industry_code"])["g_delta_loo"]
    .mean()
    .unstack("industry_code")
    .rename(columns=INDUSTRY_LABELS)
    .mul(1000)
    .round(3)
)
# Print abbreviated supersector names as column headers
pivot.columns = [c[:12] for c in pivot.columns]
with pd.option_context("display.max_columns", 12, "display.width", 130,
                       "display.float_format", "{:6.3f}".format):
    print(pivot.loc["2006Q1":"2011Q4"].to_string())

print("\n--- Peak closing quarter by supersector ---")
peak_by_ss = (
    shock_rates
    .groupby(["industry_code", "quarter_label"])["g_delta_loo"]
    .mean()
    .reset_index()
    .loc[lambda d: d.groupby("industry_code")["g_delta_loo"]
         .transform("max") == d["g_delta_loo"]]
    .assign(supersector = lambda d: d["industry_code"].map(INDUSTRY_LABELS),
            g_x1000     = lambda d: (d["g_delta_loo"] * 1000).round(3))
    [["supersector", "quarter_label", "g_x1000"]]
    .sort_values("g_x1000", ascending=False)
)
print(peak_by_ss.to_string(index=False))
