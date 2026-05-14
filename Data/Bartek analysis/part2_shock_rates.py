"""
part2_shock_rates.py — BED Deaths LOO shock rates  (g^delta_{-s,j,t})
======================================================================
Loads the employment shares saved by part1_shares.py, fetches BED Deaths
(dataclass=08: absent 4+ consecutive quarters), and computes the
leave-one-out shock rate for every (state, supersector, quarter) cell:

    g^delta_{-s,j,t} = deaths_{j,t} / E^nat_{-s,j,t0}

BED Deaths excludes temporary shutdowns by construction (an establishment
that reopens within 1–3 quarters never accumulates four consecutive
absences), so no BDS permanence calibration is required.

The numerator is the same for all states; the denominator is state-
specific because E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}.

Output
------
    data/instruments/shock_rates_{START}_{END}.parquet

Columns in shock_rates parquet:
    state_fips     2-digit state FIPS
    industry_code  2-digit supersector code
    quarter_label  e.g. "2008Q4"
    g_delta_loo    LOO death rate (employment at dying estabs / base emp.)

Run
---
    python part2_shock_rates.py

Prerequisite
------------
    python part1_shares.py   (must have run first)
"""

import sys
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — avoids all crash issues
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

def _open_file(path):
    """Open a saved plot with the OS default viewer."""
    os.startfile(path)

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(Path.home() / "Documents" / "GitHub" / "Sunk_entry_costs_endogenous_variety_unemployment" / "Data" / "Bartek analysis")

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
    build_employment_shares,
    fetch_bed_deaths_national,
    build_loo_shock_rates,
    build_bartik_instrument,
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
# Fetch BED Deaths and compute LOO shock rates
# -----------------------------------------------------------------------
print(f"\n[BED Deaths] fetching national deaths from BLS ...")
bed_deaths = fetch_bed_deaths_national(
    cache_dir     = DEFAULT_CACHE_DIR,
    start_quarter = START_QUARTER,
    end_quarter   = END_QUARTER,
)
print(
    f"  {bed_deaths['quarter_label'].nunique()} quarters  |  "
    f"{bed_deaths['industry_code'].nunique()} supersectors"
)

# Rename deaths_nat → closings_nat so build_loo_shock_rates is agnostic
bed_window = bed_deaths.rename(columns={"deaths_nat": "closings_nat"}).copy()

shock_rates = build_loo_shock_rates(shares, bed_window)

# Save
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
shock_rates.to_parquet(SHOCK_RATES_PATH, index=False)
print(f"\nSaved: {SHOCK_RATES_PATH}  ({len(shock_rates):,} rows)")

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 2 — LOO death rates  ({START_QUARTER} – {END_QUARTER})")
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

print("\n--- Mean LOO death rate by supersector (×1 000) ---")
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

print("\n--- Time series: national death rate by supersector (×1 000) ---")
print("    Grey = Great Recession, Red = COVID.\n")
pivot = (
    shock_rates
    .groupby(["quarter_label", "industry_code"])["g_delta_loo"]
    .mean()
    .unstack("industry_code")
    .rename(columns=INDUSTRY_LABELS)
    .mul(1000)
    .round(3)
)
pivot.columns = [c[:12] for c in pivot.columns]
with pd.option_context("display.max_columns", 12, "display.width", 130,
                       "display.float_format", "{:6.3f}".format):
    print(pivot.loc["2006Q1":"2011Q4"].to_string())

print("\n--- Peak death quarter by supersector ---")
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

# -----------------------------------------------------------------------
# Plot: BED death rates over time, all supersectors
# -----------------------------------------------------------------------
_nat = (
    bed_window
    .merge(
        shares[["industry_code", "emp_nat_ind"]]
        .drop_duplicates("industry_code"),
        on="industry_code",
    )
    .assign(g_nat=lambda d: d["closings_nat"] / d["emp_nat_ind"] * 1000)
)
_pivot = (
    _nat
    .pivot(index="quarter_label", columns="industry_code", values="g_nat")
    .rename(columns=INDUSTRY_LABELS)
    .sort_index()
)
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

ax.axvspan(pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
           alpha=0.10, color="grey", label="_GR")
ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
           alpha=0.10, color="red", label="_COVID")
ax.set_title(
    "National establishment death rates by supersector\n"
    "(BED Deaths, quarterly, per 1 000 base workers)",
    fontsize=12,
)
ax.set_xlabel("")
ax.set_ylabel("Death rate (×1 000)", fontsize=10)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
ax.legend(loc="upper left", fontsize=8, framealpha=0.85, ncol=2,
          title="Supersector", title_fontsize=8)
ax.grid(axis="y", linewidth=0.5, alpha=0.4)
fig.tight_layout()

outpath = DEFAULT_OUTPUT_DIR / "shock_rates_delta_by_supersector.png"
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nPlot saved: {outpath}")
