"""
part3_instrument.py — Bartik delta instrument  (B^delta_{s,t})
===============================================================
Loads the outputs of parts 1 and 2, aggregates them into the
state-quarter Bartik instrument, and saves the final CSV:

    B^delta_{s,t} = sum_j  omega_{s,j,t0} * g^delta_{-s,j,t}

No BLS data is fetched in this step — it is pure arithmetic on the
two parquet files produced by the earlier parts.

Output
------
    data/instruments/delta_instrument_base{BASE_YEAR}.csv

Columns in output:
    state           2-letter state abbreviation
    state_fips      2-digit state FIPS
    quarter_label   e.g. "2008Q4"
    bartik_delta    the instrument value B^delta_{s,t}
    weight_sum      sum of share weights used (should be ≈ 1.0)
    n_supersectors  number of supersectors contributing to this cell

Run
---
    python part3_instrument.py

Prerequisites
-------------
    python part1_shares.py      (must have run first)
    python part2_shock_rates.py (must have run second)
"""

import sys
import pandas as pd
# ---------------------------------------------------------------------------
# Working directory: set to the folder containing this script so that
# relative paths (data/cache/, data/instruments/) resolve correctly
# regardless of where Python is launched from.
#
# Path(__file__) is used when the script is run directly (e.g. python
# part3_instrument.py or F5 in VS Code with "Run Python File").
# The fallback handles interactive/REPL execution (e.g. VS Code's
# "Run Selection" or Jupyter-style terminals) where __file__ is undefined.
# ---------------------------------------------------------------------------
import os
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment\Data\Bartek analysis")

from construct_delta_instrument import (
    BASE_YEAR,
    START_QUARTER,
    END_QUARTER,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_PATH,
    INSTRUMENT_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_bartik_instrument,
)

# -----------------------------------------------------------------------
# Load Parts 1 and 2 outputs
# -----------------------------------------------------------------------
for path in [SHARES_PATH, SHOCK_RATES_PATH]:
    if not path.exists():
        script = "part1_shares.py" if "shares" in path.name else "part2_shock_rates.py"
        sys.exit(f"ERROR: {path} not found.\n       Run {script} first.")

shares      = pd.read_parquet(SHARES_PATH)
shock_rates = pd.read_parquet(SHOCK_RATES_PATH)
print(f"Loaded shares:      {SHARES_PATH}  ({len(shares):,} rows)")
print(f"Loaded shock rates: {SHOCK_RATES_PATH}  ({len(shock_rates):,} rows)")

# -----------------------------------------------------------------------
# Compute
# -----------------------------------------------------------------------
print("\n[Bartik] aggregating B^delta_{s,t} ...")
instrument = build_bartik_instrument(shares, shock_rates)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
instrument.to_csv(INSTRUMENT_PATH, index=False)
print(f"Saved: {INSTRUMENT_PATH}  ({len(instrument):,} rows)")

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 3 — Bartik delta instrument  (base year {BASE_YEAR})")
print("=" * 60)

n_states = instrument["state_fips"].nunique()
n_qtrs   = instrument["quarter_label"].nunique()
print(f"\nShape: {instrument.shape}  ({n_states} states × {n_qtrs} quarters)")

print("\nFirst 10 rows:")
print(instrument.head(10).to_string(index=False))

# Weight coverage check
low_weight = instrument[instrument["weight_sum"] < 0.80]
if low_weight.empty:
    print("\nWeight coverage: all cells ≥ 80%  ✓")
else:
    print(f"\nWARNING: {len(low_weight)} cells have weight coverage < 80%:")
    print(low_weight.to_string(index=False))

# Cross-state distribution over time
print("\n--- Cross-state distribution by quarter (×1 000) ---")
summary = (
    instrument
    .groupby("quarter_label")["bartik_delta"]
    .agg(["mean", "std", "min", "max"])
    .mul(1000)
    .round(4)
)
with pd.option_context("display.max_rows", 200):
    print(summary.to_string())

# Identify peak dispersion
std_series   = summary["std"]
peak_quarter = std_series.idxmax()
print(f"\nPeak cross-state dispersion: {peak_quarter}  "
      f"(std = {std_series[peak_quarter]:.4f} ×10⁻³)")
print("(Great Recession instrument should peak around 2008Q4–2009Q2.)")

# State ranking at peak
print(f"\n--- State ranking at {peak_quarter} ---")
peak_df = (
    instrument
    .query("quarter_label == @peak_quarter")
    .assign(bartik_x1000 = lambda d: d["bartik_delta"].mul(1000).round(4))
    .sort_values("bartik_delta", ascending=False)
    [["state", "bartik_x1000", "weight_sum", "n_supersectors"]]
)
print(peak_df.to_string(index=False))

# -----------------------------------------------------------------------
# Plot: cross-state distribution of B^delta_{s,t} over time
# Mean + 10-90 percentile band + 25-75 IQR band
# -----------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import numpy as np

    def _ql_to_dt(ql):
        y, q = ql.split("Q")
        return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

    # Build percentile panel: one row per quarter
    dist = (
        instrument
        .groupby("quarter_label")["bartik_delta"]
        .agg(
            mean  = "mean",
            p10   = lambda x: np.percentile(x, 10),
            p25   = lambda x: np.percentile(x, 25),
            p75   = lambda x: np.percentile(x, 75),
            p90   = lambda x: np.percentile(x, 90),
        )
        .mul(1000)
        .sort_index()
    )
    # Reindex to complete quarterly grid so gaps render as breaks
    all_q = pd.period_range(
        start=dist.index[0], end=dist.index[-1], freq="Q"
    ).strftime("%YQ%q").tolist()
    dist = dist.reindex(all_q)
    dates = [_ql_to_dt(q) for q in dist.index if pd.notna(dist.loc[q, "mean"])]
    dist_valid = dist.dropna(subset=["mean"])

    # Identify the two most extreme states at the peak dispersion quarter
    peak_q = summary["std"].idxmax()
    peak_vals = (
        instrument
        .query("quarter_label == @peak_q")
        .assign(b1000=lambda d: d["bartik_delta"] * 1000)
        .sort_values("b1000")
    )
    state_lo = peak_vals.iloc[0][["state", "b1000"]]
    state_hi = peak_vals.iloc[-1][["state", "b1000"]]

    # State time series for the two extreme states
    def _state_series(abbrev):
        s = (
            instrument[instrument["state"] == abbrev]
            .set_index("quarter_label")["bartik_delta"]
            .mul(1000)
            .reindex(all_q)
        )
        dates_s = [_ql_to_dt(q) for q in s.index]
        return dates_s, s.values

    fig, ax = plt.subplots(figsize=(13, 5))

    # Shaded bands
    ax.fill_between(
        [_ql_to_dt(q) for q in dist_valid.index],
        dist_valid["p10"], dist_valid["p90"],
        alpha=0.18, color="#1f77b4", label="10–90th pctile",
    )
    ax.fill_between(
        [_ql_to_dt(q) for q in dist_valid.index],
        dist_valid["p25"], dist_valid["p75"],
        alpha=0.30, color="#1f77b4", label="25–75th pctile (IQR)",
    )

    # Mean
    ax.plot(
        [_ql_to_dt(q) for q in dist_valid.index],
        dist_valid["mean"],
        color="#1f77b4", linewidth=2.0, label="Cross-state mean",
    )

    # Extreme state reference lines
    ds, vs = _state_series(state_lo["state"])
    ax.plot(ds, vs, color="firebrick", linewidth=0.9, linestyle=":",
            label=f"{state_lo['state']} (lowest at {peak_q})")
    ds, vs = _state_series(state_hi["state"])
    ax.plot(ds, vs, color="darkgreen", linewidth=0.9, linestyle=":",
            label=f"{state_hi['state']} (highest at {peak_q})")

    # Recession shading
    ax.axvspan(pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
               alpha=0.08, color="grey")
    ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
               alpha=0.10, color="red")

    ax.set_title(
        r"Bartik $\delta$-instrument $B^\delta_{s,t}$: cross-state distribution over time"
        "\n(quarterly, ×1 000 base workers; base year 2006)",
        fontsize=11,
    )
    ax.set_ylabel(r"$B^\delta_{s,t}$ (×1 000)", fontsize=10)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x:.3f}")
    )
    ax.legend(fontsize=8, framealpha=0.85, ncol=2)
    ax.grid(axis="y", linewidth=0.5, alpha=0.4)
    fig.tight_layout()

    outpath = DEFAULT_OUTPUT_DIR / "instrument_delta_distribution.png"
    fig.savefig(outpath, dpi=150)
    plt.show()
    plt.close(fig)
    print(f"\nPlot saved: {outpath}")

except ImportError:
    print("\n(matplotlib not available -- skipping plot)")