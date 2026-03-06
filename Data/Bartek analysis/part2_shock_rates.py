"""
part2_shock_rates.py — Permanence-adjusted national LOO closing rates  (g^delta_{-s,j,t})
==========================================================================================
Loads the employment shares saved by part1_shares.py, fetches BED national
closing employment, applies the BDS/BED permanence calibration, and computes
the leave-one-out shock rate for every (state, supersector, quarter) cell:

    g^delta_{-s,j,t} = closings^perm_{j,t} / E^nat_{-s,j,t0}

where closings^perm_{j,t} = π_{j,y(t)} × BED closings_{j,t} and π_{j,y} is
the fraction of BED closings in supersector j × BDS year y that proved
permanent (not temporary shutdowns).

The numerator is the same for all states; the denominator is state-
specific because E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}.

Output
------
    data/instruments/shock_rates_{START}_{END}.parquet
    data/instruments/permanence_ratios_by_supersector.parquet  (diagnostic)

Columns in shock_rates parquet:
    state_fips     2-digit state FIPS
    industry_code  2-digit supersector code
    quarter_label  e.g. "2008Q4"
    g_delta_loo    LOO adjusted-closing rate (perm. employment lost / base emp.)

Run
---
    python part2_shock_rates.py

Prerequisite
------------
    python part1_shares.py   (must have run first)

Note on Census API key
----------------------
The BDS fetch in this step calls the Census Bureau API.  Without a key,
requests are limited to ~500/day per IP.  Export your key as an environment
variable before running:

    export CENSUS_API_KEY=your_key_here    (Linux/macOS)
    set CENSUS_API_KEY=your_key_here       (Windows)

Or pass it directly: build_national_shock_rates(..., census_key="your_key")
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


from construct_delta_instrument import (
    BASE_YEAR,
    START_QUARTER,
    END_QUARTER,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_PATH,
    PERM_RATIOS_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_national_shock_rates,
    fetch_bed_closings_national,
    fetch_bds_exits_national,
    compute_permanence_ratios,
    apply_permanence_adjustment,
)

# Optional Census API key — set CENSUS_API_KEY env var or edit here
import os
CENSUS_KEY = os.environ.get("CENSUS_API_KEY", "")

# Path to locally downloaded BDS national-by-sector CSV (preferred).
# Download once from:
#   https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html
#   → National → Sector table  (e.g. bds2023_sec_nat.csv)
# Set to None to fall back to the Census API instead.
BDS_FILE = Path("data/raw/bds2023_sec_nat.csv")   # edit filename as needed

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
    census_key    = CENSUS_KEY,
    bds_file      = BDS_FILE,
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

print("\n--- Time series: national closing rate by supersector (×1 000) ---")
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

print("\n--- Raw BED closings by supersector (sanity check) ---")
print("    Closings should be in the hundreds of thousands per quarter.")
print("    Near-zero values indicate a broken or missing BED series.\n")
bed_nat = fetch_bed_closings_national(
    cache_dir     = DEFAULT_CACHE_DIR,
    start_quarter = START_QUARTER,
    end_quarter   = END_QUARTER,
)
bed_summary = (
    bed_nat
    .groupby("industry_code")["closings_nat"]
    .agg(["mean", "min", "max", "count"])
    .rename(index=INDUSTRY_LABELS)
    .rename(columns={"mean": "mean_jobs", "min": "min_jobs",
                     "max": "max_jobs", "count": "n_quarters"})
    .sort_values("mean_jobs", ascending=False)
    .round(0)
)
print(bed_summary.to_string())
print("\n--- Trade/transport raw values (first 20 quarters) ---")
tt = bed_nat[bed_nat["industry_code"] == "40"].sort_values("quarter_label")
print(tt.head(20).to_string(index=False))

# -----------------------------------------------------------------------
# Permanence ratio diagnostic
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print("Permanence ratios π_{j,y}  (BDS exits / BED closings sum)")
print("=" * 60)

if PERM_RATIOS_PATH.exists():
    perm = pd.read_parquet(PERM_RATIOS_PATH)

    # Cross-supersector mean π by year
    print("\n--- Mean π by year (all supersectors) ---")
    pi_by_year = (
        perm.groupby("bds_year")["pi"]
        .mean()
        .round(3)
    )
    with pd.option_context("display.max_rows", 50):
        print(pi_by_year.to_string())

    # π by supersector for key recession years
    key_years = [y for y in [2009, 2010, 2020, 2021, 2022] if y in perm["bds_year"].values]
    if key_years:
        print(f"\n--- π by supersector for key years {key_years} ---")
        pivot_pi = (
            perm[perm["bds_year"].isin(key_years)]
            .pivot(index="industry_code", columns="bds_year", values="pi")
            .rename(index=INDUSTRY_LABELS)
            .round(3)
        )
        print(pivot_pi.to_string())

    # Flag unusually low π (potential temporary-closure episodes)
    low_pi = perm[perm["pi"] < 0.5].sort_values("pi")
    if not low_pi.empty:
        print(f"\n--- Cells with π < 0.50 (high temporary-closure contamination) ---")
        low_pi_display = low_pi.assign(
            supersector=lambda d: d["industry_code"].map(INDUSTRY_LABELS)
        )[["supersector", "bds_year", "pi", "bds_exits", "bed_closings_sum"]].round(3)
        print(low_pi_display.to_string(index=False))
    else:
        print("\nNo cells with π < 0.50.")
else:
    print(f"  (Permanence ratios not found at {PERM_RATIOS_PATH} — run build step first)")

# -----------------------------------------------------------------------
# Plot 1: Permanence ratios π_{j,y} by supersector over time
# -----------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    if PERM_RATIOS_PATH.exists():
        perm_plot = pd.read_parquet(PERM_RATIOS_PATH)
        perm_pivot = (
            perm_plot
            .pivot(index="bds_year", columns="industry_code", values="pi")
            .rename(columns=INDUSTRY_LABELS)
            .sort_index()
        )
        colors_pi = [
            "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]
        fig_pi, ax_pi = plt.subplots(figsize=(13, 5))
        for i, col in enumerate(perm_pivot.columns):
            ax_pi.plot(
                perm_pivot.index, perm_pivot[col],
                label=col, color=colors_pi[i % len(colors_pi)],
                linewidth=1.4, marker="o", markersize=3,
            )
        ax_pi.axhline(1.0, color="black", linewidth=0.8, linestyle=":", alpha=0.5)
        ax_pi.axvspan(2008, 2010.5, alpha=0.10, color="grey", label="_GR")
        ax_pi.axvspan(2020, 2022.5, alpha=0.10, color="red",  label="_COVID")
        ax_pi.set_title(
            r"Permanence ratio $\pi_{j,y}$ = BDS exits / BED closings sum"
            "\nby supersector and BDS year  (π ≈ 1 → all closings permanent; "
            "π ↓ in 2021 → many COVID closings were temporary)",
            fontsize=10,
        )
        ax_pi.set_xlabel("BDS year")
        ax_pi.set_ylabel(r"$\pi_{j,y}$", fontsize=10)
        ax_pi.set_ylim(0, 1.15)
        ax_pi.legend(fontsize=7, framealpha=0.85, ncol=2,
                     title="Supersector", title_fontsize=7)
        ax_pi.grid(axis="y", linewidth=0.5, alpha=0.4)
        fig_pi.tight_layout()
        outpath_pi = DEFAULT_OUTPUT_DIR / "permanence_ratios_by_supersector.png"
        DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        fig_pi.savefig(outpath_pi, dpi=150)
        plt.close(fig_pi)
        print(f"\nPlot saved: {outpath_pi}")

except ImportError:
    print("\n(matplotlib not available -- skipping permanence ratio plot)")

# -----------------------------------------------------------------------
# Plot 2: LOO adjusted-closing rates over time, all supersectors
# -----------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
   

# Build national closing rate: closings_nat / emp_nat (no LOO).
# The LOO perturbation in the denominator is negligible for plotting
# purposes and the numerator is identical across states -- so the
# correct thing to show is the underlying national shock series.
_nat = (
    bed_nat                          # fetched in the sanity-check block above
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
# Reindex to a complete quarterly grid so matplotlib shows gaps as breaks
# rather than drawing misleading straight lines across missing quarters.
all_quarters = pd.period_range(
    start=_pivot.index[0], end=_pivot.index[-1], freq="Q"
).strftime("%YQ%q").tolist()
plot_data = _pivot.reindex(all_quarters)

# Parse quarter_label to datetime for a clean x-axis
def _ql_to_dt(ql):
    import pandas as pd
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

plot_data.index = [_ql_to_dt(q) for q in plot_data.index]

# Color palette: 10 visually distinct colors
colors = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
# Line styles: solid for top-5 by mean, dashed for bottom-5
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
    "National establishment closing rates by supersector\n"
    "(quarterly, per 1 000 base workers)",
    fontsize=12,
)
ax.set_xlabel("")
ax.set_ylabel("Closing rate (×1 000)", fontsize=10)
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
plt.show()
outpath = DEFAULT_OUTPUT_DIR / "shock_rates_delta_by_supersector.png"
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nPlot saved: {outpath}")

# Permanence ratios by industry 
# Load 
pi = pd.read_parquet("data/instruments/permanence_ratios_by_supersector.parquet")
print(pi.pivot(index="bds_year", columns="industry_code", values="pi").round(3).to_string())

pivot = pi.pivot(index="bds_year", columns="industry_code", values="pi")
pivot.columns = [INDUSTRY_LABELS[c] for c in pivot.columns]
pivot.plot(figsize=(14, 5), title="Permanence ratios π_{j,y} by industry")
plt.axhline(1.0, color="black", linewidth=0.8, linestyle="--")
plt.ylabel("π  (fraction of BED closings that are permanent exits)")
plt.tight_layout()
plt.show()