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
import os
import numpy as np
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
# Rescale to quarterly rates (pure dimensionless fractions)
# BED closings and JOLTS separations are in thousands of workers;
# QCEW emp_nat_ind is in raw worker counts → multiply by 1000 to correct.
# Resulting units: quarterly rate as fraction of base-year employment.
# -----------------------------------------------------------------------
instrument["bartik_delta"] = instrument["bartik_delta"] * 1000

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
print(f"Part 3 Bartik delta instrument  (base year {BASE_YEAR})")
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
print("\n--- Cross-state distribution by quarter (quarterly rate) ---")
summary = (
    instrument
    .groupby("quarter_label")["bartik_delta"]
    .agg(["mean", "std", "min", "max"])
    .round(6)
)
with pd.option_context("display.max_rows", 200):
    print(summary.to_string())

# Identify peak dispersion
std_series   = summary["std"]
peak_quarter = std_series.idxmax()
print(f"\nPeak cross-state dispersion: {peak_quarter}  "
      f"(std = {std_series[peak_quarter]:.6f})")
print("(Great Recession instrument should peak around 2008Q4–2009Q2.)")

# State ranking at peak
print(f"\n--- State ranking at {peak_quarter} ---")
peak_df = (
    instrument
    .query("quarter_label == @peak_quarter")
    .assign(bartik_delta = lambda d: d["bartik_delta"].round(6))
    .sort_values("bartik_delta", ascending=False)
    [["state", "bartik_delta", "weight_sum", "n_supersectors"]]
)
print(peak_df.to_string(index=False))

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Plot: two-panel figure
#   Left  — residualized instrument in original units (pp), with
#            10-90th and 25-75th cross-state percentile bands
#   Right — z-score overlay: residualized (blue) vs raw (red) cross-
#            state means, to show what residualization changes
# -----------------------------------------------------------------------
def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

# Load residualized instrument (full available range: 1997Q1 onward via Chow-Lin)
resid_path = DEFAULT_OUTPUT_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
if not resid_path.exists():
    print(f"WARNING: {resid_path} not found; run part2b_residualize_shocks.py first.")
    resid = None
else:
    resid = pd.read_csv(resid_path)
    # No sample restriction — show full range from Chow-Lin extension (1997Q1)
    resid = resid[resid["quarter_label"] >= "1997Q1"].copy()

# Both panels: 1997Q1-2019Q4 (full range where both raw and resid available;
# stop at 2019Q4 to exclude COVID from the comparison window)
raw_lp = instrument[
    (instrument["quarter_label"] >= "1997Q1") &
    (instrument["quarter_label"] <= "2019Q4")
].copy()
resid_lp = resid[
    (resid["quarter_label"] >= "1997Q1") &
    (resid["quarter_label"] <= "2019Q4")
].copy() if resid is not None else None

def _build_dist(df, col="bartik_delta"):
    d = (
        df.groupby("quarter_label")[col]
        .agg(
            mean = "mean",
            p10  = lambda x: np.percentile(x, 10),
            p25  = lambda x: np.percentile(x, 25),
            p75  = lambda x: np.percentile(x, 75),
            p90  = lambda x: np.percentile(x, 90),
        )
        .sort_index()
    )
    all_q = pd.period_range(
        start=d.index[0], end=d.index[-1], freq="Q"
    ).strftime("%YQ%q").tolist()
    return d.reindex(all_q)

def _zscored_mean(dist_df):
    m = dist_df["mean"].dropna()
    return (m - m.mean()) / m.std()

def _shade(ax, lp_sample=False):
    NBER = [("2007-12-01","2009-06-01"), ("2001-01-01","2001-12-01"), ("2020-01-01","2020-07-01")]
    for r0, r1 in NBER:
        ax.axvspan(pd.Timestamp(r0), pd.Timestamp(r1),
                   color="grey", alpha=0.12, zorder=0)
    if lp_sample:
        # Shade the pre-LP region (1997Q1-2000Q4) to distinguish from main sample
        ax.axvspan(pd.Timestamp("1997-01-01"), pd.Timestamp("2001-01-01"),
                   color="#f5a623", alpha=0.08, zorder=0)
        ax.axvline(pd.Timestamp("2001-01-01"), color="#f5a623",
                   lw=1.0, ls="--", alpha=0.6)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.subplots_adjust(wspace=0.28)

# ── Left panel: residualized instrument, original units ──────────────
ax1 = axes[0]
if resid is not None:
    dist_r = _build_dist(resid)
    dv     = dist_r.dropna(subset=["mean"])
    xs     = [_ql_to_dt(q) for q in dv.index]
    _shade(ax1, lp_sample=True)
    ax1.fill_between(xs, dv["p10"], dv["p90"],
                     alpha=0.18, color="#1f77b4", label="10–90th pctile")
    ax1.fill_between(xs, dv["p25"], dv["p75"],
                     alpha=0.32, color="#1f77b4", label="25–75th pctile (IQR)")
    ax1.plot(xs, dv["mean"], color="#1f77b4", lw=2.0, label="Cross-state mean")
    ax1.axhline(0, color="black", lw=0.7, ls="--", alpha=0.5)
    sigma = resid["bartik_delta"].std()
    ax1.set_ylabel(r"$\tilde{B}^\delta_{s,t}$ (pp, residualized)", fontsize=10)
    ax1.set_title(
        r"Residualized $\delta$ instrument" + "\n" +
        r"($\hat{\sigma}=" + f"{sigma:.1f}" + r"$ pp; base year 2006)",
        fontsize=10,
    )
    ax1.legend(fontsize=8, framealpha=0.85)
    ax1.grid(axis="y", lw=0.4, alpha=0.4)
    ax1.tick_params(labelsize=9)

else:
    ax1.text(0.5, 0.5, "Residualized instrument not found",
             ha="center", va="center", transform=ax1.transAxes)

# Right panel: z-score overlay raw vs residualized 
ax2 = axes[1]
dist_raw = _build_dist(raw_lp)
dv_raw   = dist_raw.dropna(subset=["mean"])
_shade(ax2, lp_sample=True)

z_raw  = _zscored_mean(dist_raw)
xs_raw = [_ql_to_dt(q) for q in z_raw.index]
ax2.plot(xs_raw, z_raw.values, color="#d62728", lw=1.8,
         label="Raw Bartik (z-score)")

if resid_lp is not None:
    dist_r_lp = _build_dist(resid_lp)
    z_res  = _zscored_mean(dist_r_lp)
    xs_res = [_ql_to_dt(q) for q in z_res.index]
    ax2.plot(xs_res, z_res.values, color="#1f77b4", lw=1.8, ls="--",
             label="Residualized (z-score)")
    common = z_raw.index.intersection(z_res.index)
    corr   = np.corrcoef(z_raw.loc[common].values, z_res.loc[common].values)[0, 1]
    ax2.annotate(f"$r = {corr:.3f}$",
                 xy=(0.05, 0.93), xycoords="axes fraction",
                 fontsize=9, color="black",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

ax2.axhline(0, color="black", lw=0.7, ls="-", alpha=0.3)
ax2.set_ylabel("Standardized units (z-score)", fontsize=10)
ax2.set_title(
    "Raw vs residualized instrument\n(z-scored means, 1997Q1--2019Q4)",
    fontsize=10,
)
ax2.legend(fontsize=8, framealpha=0.85)
ax2.grid(axis="y", lw=0.4, alpha=0.4)
ax2.tick_params(labelsize=9)

fig.suptitle(
    r"Bartik $\delta$-instrument: cross-state distribution, 1997Q1–2019Q4  (LP sample shaded: 2001Q1–2019Q4)",
    fontsize=11, y=1.01,
)
fig.tight_layout()

outpath = DEFAULT_OUTPUT_DIR / "instrument_delta_distribution.png"
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(outpath, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nPlot saved: {outpath}")