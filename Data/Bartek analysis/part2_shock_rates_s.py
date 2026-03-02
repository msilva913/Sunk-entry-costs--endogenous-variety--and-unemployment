"""
part2_shock_rates_s.py -- National LOO layoff rates  (g^s_{-s,j,t})
=====================================================================
Loads the employment shares saved by part1_shares.py, fetches national
JOLTS Layoffs & Discharges from the BLS public API, and computes the
leave-one-out shock rate for every (state, supersector, quarter) cell:

    g^s_{-s,j,t} = layoffs^nat_{j,t} / E^nat_{-s,j,t0}

Denominator LOO applied universally (JOLTS does not publish state-level
layoffs at the supersector level via the public API).

Output
------
    data/instruments/shock_rates_s_{START}_{END}.parquet

Columns in output:
    state_fips     2-digit state FIPS
    industry_code  2-digit supersector code (pipeline convention)
    quarter_label  e.g. "2008Q4"
    g_s_loo        LOO layoff rate (layoffs / base employment)

Run
---
    python part2_shock_rates_s.py

Prerequisite
------------
    python part1_shares.py   (must have run first)

Sample window
-------------
    2001Q1 - 2024Q2  (first complete quarter of national JOLTS coverage)
"""

import sys
import pandas as pd
# ---------------------------------------------------------------------------
# Working directory: set to the folder containing this script so that
# relative paths (data/cache/, data/instruments/) resolve correctly
# regardless of where Python is launched from.
#
# Path(__file__) is used when the script is run directly (e.g. python
# part2_shock_rates.py or F5 in VS Code with "Run Python File").
# The fallback handles interactive/REPL execution (e.g. VS Code's
# "Run Selection" or Jupyter-style terminals) where __file__ is undefined.
# ---------------------------------------------------------------------------
import os
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment\Data\Bartek analysis")

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
    fetch_jolts_layoffs_national,
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
print(f"Part 2s -- LOO s-shock rates  ({START_QUARTER} - {END_QUARTER})")
print("=" * 60)

n_states = shock_rates["state_fips"].nunique()
n_ss     = shock_rates["industry_code"].nunique()
n_qtrs   = shock_rates["quarter_label"].nunique()
n_valid  = shock_rates["g_s_loo"].notna().sum()
print(f"\nShape: {shock_rates.shape}  "
      f"({n_states} states x {n_ss} supersectors x {n_qtrs} quarters)")
print(f"Valid cells: {n_valid:,} / {len(shock_rates):,}")

print("\nFirst 10 rows:")
print(shock_rates.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Mean LOO rate by supersector
# -----------------------------------------------------------------------
print("\n--- Mean LOO layoff rate by supersector (x1 000) ---")
mean_rate = (
    shock_rates.groupby("industry_code")["g_s_loo"]
    .mean()
    .rename(index=INDUSTRY_LABELS)
    .sort_values(ascending=False)
    .mul(1000)
    .round(3)
)
for name, rate in mean_rate.items():
    bar = chr(9608) * int(rate * 6)
    print(f"  {name:<36}  {rate:6.3f}  {bar}")

# -----------------------------------------------------------------------
# Time series: national layoff rate by supersector
# -----------------------------------------------------------------------
print("\n--- Time series: national layoff rate by supersector (x1 000) ---")
print("    Grey = Great Recession, Red = COVID.\n")
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
# Sanity check: raw national JOLTS levels
# -----------------------------------------------------------------------
print("\n--- Raw national JOLTS layoff levels by supersector (sanity check) ---")
print("    All levels in thousands of workers per quarter.")
print("    Plausible ranges: ~150k (Mining) to ~5 000k (Trade/transport)\n")
jolts_nat = fetch_jolts_layoffs_national(
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

# -----------------------------------------------------------------------
# Plot: LOO layoff rates over time, all supersectors on one axis
# -----------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import pandas as pd

    # Build national layoff rate: layoffs_nat / emp_nat (no LOO).
    # The numerator is identical across states; showing the national rate
    # directly is the correct representation of the shock series.
    _nat = (
        jolts_nat                        # fetched in the sanity-check block above
        .merge(
            shares[["industry_code", "emp_nat_ind"]]
            .drop_duplicates("industry_code"),
            on="industry_code",
        )
        .assign(g_nat=lambda d: d["layoffs_nat"] / d["emp_nat_ind"] * 1000)
    )
    _pivot = (
        _nat
        .pivot(index="quarter_label", columns="industry_code", values="g_nat")
        .rename(columns=INDUSTRY_LABELS)
        .sort_index()
    )
    # Reindex to complete quarterly grid to expose any gaps as breaks
    all_quarters = pd.period_range(
        start=_pivot.index[0], end=_pivot.index[-1], freq="Q"
    ).strftime("%YQ%q").tolist()
    plot_data = _pivot.reindex(all_quarters)

    def _ql_to_dt(ql):
        y, q = ql.split("Q")
        return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

    plot_data.index = [_ql_to_dt(q) for q in plot_data.index]

    colors = [
        "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]
    mean_rank = plot_data.mean().sort_values(ascending=False)
    top5      = set(mean_rank.index[:5])

    fig, ax = plt.subplots(figsize=(13, 6))

    for i, col in enumerate(mean_rank.index):
        ls = "-" if col in top5 else "--"
        lw = 1.6 if col in top5 else 1.2
        ax.plot(
            plot_data.index,
            plot_data[col],
            label     = col,
            color     = colors[i % len(colors)],
            linestyle = ls,
            linewidth = lw,
        )

    # Recession shading: GR and COVID
    ax.axvspan(
        pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
        alpha=0.10, color="grey", label="_GR"
    )
    ax.axvspan(
        pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
        alpha=0.10, color="red", label="_COVID"
    )

    ax.set_title(
        "National layoff & discharge rates by supersector\n"
        "(quarterly, per 1 000 base workers)",
        fontsize=12,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Layoff rate (x1 000)", fontsize=10)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
    ax.legend(
        loc            = "upper left",
        fontsize       = 8,
        framealpha     = 0.85,
        ncol           = 2,
        title          = "Supersector",
        title_fontsize = 8,
    )
    ax.grid(axis="y", linewidth=0.5, alpha=0.4)
    fig.tight_layout()

    outpath = DEFAULT_OUTPUT_DIR / "shock_rates_s_by_supersector.png"
    DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=150)
    plt.show()
    plt.close(fig)
    print(f"\nPlot saved: {outpath}")

except ImportError:
    print("\n(matplotlib not available -- skipping plot)")