"""
part3_instrument_s.py — Bartik s-shock instrument  (B^s_{s,t})
===============================================================
Loads the outputs of part1_shares.py and part2_shock_rates_s.py,
aggregates them into the state-quarter Bartik s-instrument, and saves
the final CSV:

    B^s_{s,t} = sum_j  omega_{s,j,t0} * g^s_{-s,j,t}

No BLS data is fetched here — pure arithmetic on the two parquet files
produced by the earlier parts.

Output
------
    data/instruments/s_instrument_base{BASE_YEAR}.csv

Columns in output:
    state           2-letter state abbreviation
    state_fips      2-digit state FIPS
    quarter_label   e.g. "2008Q4"
    bartik_s        the instrument value B^s_{s,t}
    weight_sum      sum of share weights used (should be ≈ 1.0)
    n_supersectors  number of supersectors contributing to this cell
    n_full_loo      number of supersectors where full LOO was applied

Run
---
    python part3_instrument_s.py

Prerequisites
-------------
    python part1_shares.py          (must have run first)
    python part2_shock_rates_s.py   (must have run second)

Joint use with delta instrument
--------------------------------
The s and delta instruments share the same employment weights
(shares_base{YEAR}.parquet).  To use both in a joint IV regression:

    import pandas as pd
    delta = pd.read_csv("data/instruments/delta_instrument_base2006.csv")
    s     = pd.read_csv("data/instruments/s_instrument_base2006.csv")
    both  = delta.merge(s, on=["state", "state_fips", "quarter_label"])
    # both now has bartik_delta and bartik_s aligned on the same
    # (state, quarter) grid, ready for 2SLS or GMM estimation.
    # Note: delta runs from 1992Q3; s runs from 2003Q1.
    # Joint analysis should restrict to 2003Q1 onward.
"""

import sys
import pandas as pd
from construct_s_instrument import (
    BASE_YEAR,
    START_QUARTER,
    END_QUARTER,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_S_PATH,
    S_INSTRUMENT_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_bartik_instrument_s,
)

# -----------------------------------------------------------------------
# Load Parts 1 and 2 outputs
# -----------------------------------------------------------------------
for path in [SHARES_PATH, SHOCK_RATES_S_PATH]:
    if not path.exists():
        if "shares" in path.name:
            script = "part1_shares.py"
        else:
            script = "part2_shock_rates_s.py"
        sys.exit(f"ERROR: {path} not found.\n       Run {script} first.")

shares        = pd.read_parquet(SHARES_PATH)
shock_rates_s = pd.read_parquet(SHOCK_RATES_S_PATH)
print(f"Loaded shares:        {SHARES_PATH}  ({len(shares):,} rows)")
print(f"Loaded shock rates s: {SHOCK_RATES_S_PATH}  ({len(shock_rates_s):,} rows)")

# -----------------------------------------------------------------------
# Compute
# -----------------------------------------------------------------------
print("\n[Bartik] aggregating B^s_{s,t} ...")
instrument = build_bartik_instrument_s(shares, shock_rates_s)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
instrument.to_csv(S_INSTRUMENT_PATH, index=False)
print(f"Saved: {S_INSTRUMENT_PATH}  ({len(instrument):,} rows)")

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 3s — Bartik s-instrument  (base year {BASE_YEAR})")
print("=" * 60)

n_states = instrument["state_fips"].nunique()
n_qtrs   = instrument["quarter_label"].nunique()
print(f"\nShape: {instrument.shape}  ({n_states} states × {n_qtrs} quarters)")

print("\nFirst 10 rows:")
print(instrument.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Weight coverage check
# -----------------------------------------------------------------------
low_weight = instrument[instrument["weight_sum"] < 0.80]
if low_weight.empty:
    print("\nWeight coverage: all cells ≥ 80%  ✓")
else:
    print(f"\nWARNING: {len(low_weight)} cells have weight coverage < 80%:")
    print(low_weight.to_string(index=False))

# -----------------------------------------------------------------------
# LOO quality: distribution of n_full_loo across state-quarters
# -----------------------------------------------------------------------
print("\n--- Full-LOO supersector count per state-quarter ---")
print("    Max possible = 10 (all supersectors have state data)")
loo_dist = instrument["n_full_loo"].value_counts().sort_index()
print(loo_dist.to_string())
pct_all_full = (instrument["n_full_loo"] == instrument["n_supersectors"]).mean()
print(f"Cells with all supersectors on full LOO: {pct_all_full:.1%}")

# -----------------------------------------------------------------------
# Cross-state distribution over time
# -----------------------------------------------------------------------
print("\n--- Cross-state distribution by quarter (×1 000) ---")
summary = (
    instrument
    .groupby("quarter_label")["bartik_s"]
    .agg(["mean", "std", "min", "max"])
    .mul(1000)
    .round(4)
)
with pd.option_context("display.max_rows", 200):
    print(summary.to_string())

# -----------------------------------------------------------------------
# Peak dispersion quarter and state ranking
# -----------------------------------------------------------------------
std_series   = summary["std"]
peak_quarter = std_series.idxmax()
print(f"\nPeak cross-state dispersion: {peak_quarter}  "
      f"(std = {std_series[peak_quarter]:.4f} ×10⁻³)")
print("(Great Recession or COVID should produce the largest cross-state spread.)")

print(f"\n--- State ranking at {peak_quarter} ---")
peak_df = (
    instrument
    .query("quarter_label == @peak_quarter")
    .assign(bartik_x1000 = lambda d: d["bartik_s"].mul(1000).round(4))
    .sort_values("bartik_s", ascending=False)
    [["state", "bartik_x1000", "weight_sum", "n_supersectors", "n_full_loo"]]
)
print(peak_df.to_string(index=False))

# -----------------------------------------------------------------------
# Comparison with delta instrument (if available)
# -----------------------------------------------------------------------
delta_path = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
if delta_path.exists():
    print("\n--- Correlation with delta instrument (2003Q1 onward) ---")
    delta = pd.read_csv(delta_path)
    both = (
        instrument
        [["state_fips", "quarter_label", "bartik_s"]]
        .merge(
            delta[["state_fips", "quarter_label", "bartik_delta"]],
            on=["state_fips", "quarter_label"],
        )
        .query("quarter_label >= '2003Q1'")
    )
    corr = both[["bartik_s", "bartik_delta"]].corr().iloc[0, 1]
    print(f"  Pearson r(bartik_s, bartik_delta) = {corr:.4f}")
    print(f"  N = {len(both):,} state-quarter observations")
    print()
    if abs(corr) > 0.90:
        print("  NOTE: High correlation — instruments may be nearly collinear.")
        print("  Consider checking whether they carry independent identifying")
        print("  variation before using both in a joint IV specification.")
    else:
        print("  Correlation < 0.90: instruments appear to carry independent")
        print("  variation.  Suitable for joint IV estimation.")
else:
    print(
        f"\n(Delta instrument not found at {delta_path} — "
        f"skipping correlation check.)"
    )
