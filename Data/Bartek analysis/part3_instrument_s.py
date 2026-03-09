"""
part3_instrument_s.py -- Bartik s-shock instrument  (B^s_{s,t})
================================================================
Loads the outputs of part1_shares.py and part2_shock_rates_s.py,
aggregates them into the state-quarter Bartik s-instrument, and saves
the final CSV:

    B^s_{s,t} = sum_j  omega_{s,j,t0} * g^s_{-s,j,t}

No BLS data is fetched here -- pure arithmetic on the two parquet files
produced by the earlier parts.

Output
------
    data/instruments/s_instrument_base{BASE_YEAR}.csv

Columns in output:
    state           2-letter state abbreviation
    state_fips      2-digit state FIPS
    quarter_label   e.g. "2008Q4"
    bartik_s        the instrument value B^s_{s,t}
    weight_sum      sum of share weights used (should be approx 1.0)
    n_supersectors  number of supersectors contributing to this cell

Run
---
    python part3_instrument_s.py

Prerequisites
-------------
    python part1_shares.py          (must have run first)
    python part2_shock_rates_s.py   (must have run second)

Joint use with delta instrument
---------------------------------
The s and delta instruments share the same employment weights.
To use both in a joint IV regression:

    import pandas as pd
    delta = pd.read_csv("data/instruments/delta_instrument_base2006.csv")
    s     = pd.read_csv("data/instruments/s_instrument_base2006.csv")
    both  = delta.merge(s, on=["state", "state_fips", "quarter_label"])
    # both has bartik_delta and bartik_s aligned on the same grid.
    # delta runs from 1992Q3; s runs from 2001Q1.
    # Joint analysis should restrict to 2001Q1 onward.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

def _open_file(path):
    os.startfile(path)

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(Path.home() / "Documents" / "GitHub" / "Sunk_entry_costs_endogenous_variety_unemployment" / "Data" / "Bartek analysis")

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
        script = "part1_shares.py" if "shares" in path.name else "part2_shock_rates_s.py"
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
print(f"Part 3s -- Bartik s-instrument  (base year {BASE_YEAR})")
print("=" * 60)

n_states = instrument["state_fips"].nunique()
n_qtrs   = instrument["quarter_label"].nunique()
print(f"\nShape: {instrument.shape}  ({n_states} states x {n_qtrs} quarters)")

print("\nFirst 10 rows:")
print(instrument.head(10).to_string(index=False))

# -----------------------------------------------------------------------
# Weight coverage check
# -----------------------------------------------------------------------
low_weight = instrument[instrument["weight_sum"] < 0.80]
if low_weight.empty:
    print("\nWeight coverage: all cells >= 80%  OK")
else:
    print(f"\nWARNING: {len(low_weight)} cells have weight coverage < 80%:")
    print(low_weight.to_string(index=False))

# -----------------------------------------------------------------------
# Cross-state distribution over time
# -----------------------------------------------------------------------
print("\n--- Cross-state distribution by quarter (x1 000) ---")
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
      f"(std = {std_series[peak_quarter]:.4f} x10^-3)")
print("(Expect 2009 or 2020 given Great Recession and COVID separation spikes."
      "  Note: with total separations (TS) instead of layoffs (LD), the"
      "  COVID spike may be attenuated since quits collapsed in 2020Q2.)")

print(f"\n--- State ranking at {peak_quarter} ---")
peak_df = (
    instrument
    .query("quarter_label == @peak_quarter")
    .assign(bartik_x1000 = lambda d: d["bartik_s"].mul(1000).round(4))
    .sort_values("bartik_s", ascending=False)
    [["state", "bartik_x1000", "weight_sum", "n_supersectors"]]
)
print(peak_df.to_string(index=False))

# -----------------------------------------------------------------------
# Correlation with delta instrument (if available)
# -----------------------------------------------------------------------
delta_path = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
if delta_path.exists():
    print("\n--- Correlation with delta instrument (2001Q1 onward) ---")
    delta = pd.read_csv(delta_path)
    both = (
        instrument[["state_fips", "quarter_label", "bartik_s"]]
        .merge(
            delta[["state_fips", "quarter_label", "bartik_delta"]],
            on=["state_fips", "quarter_label"],
        )
        .query("quarter_label >= '2001Q1'")
    )
    corr = both[["bartik_s", "bartik_delta"]].corr().iloc[0, 1]
    print(f"  Pearson r(bartik_s, bartik_delta) = {corr:.4f}")
    print(f"  N = {len(both):,} state-quarter observations")
    if abs(corr) > 0.90:
        print("\n  NOTE: High correlation -- instruments may be nearly collinear.")
        print("  Verify they carry independent variation before joint IV use.")
    else:
        print("\n  Correlation < 0.90: instruments appear to carry independent")
        print("  variation.  Suitable for joint IV estimation.")
else:
    print(
        f"\n(Delta instrument not found at {delta_path} -- "
        f"skipping correlation check.)"
    )

# -----------------------------------------------------------------------
# Plot: two-panel comparison of B^delta and B^s distributions over time
# -----------------------------------------------------------------------
delta_path = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
if not delta_path.exists():
    print(f"\n(Two-panel plot skipped: {delta_path} not found)")
else:
    delta = pd.read_csv(delta_path)

    def _ql_to_dt(ql):
        y, q = ql.split("Q")
        return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

    def _build_dist(df, col):
        raw = (
            df.groupby("quarter_label")[col]
            .agg(
                mean = "mean",
                p10  = lambda x: np.percentile(x, 10),
                p25  = lambda x: np.percentile(x, 25),
                p75  = lambda x: np.percentile(x, 75),
                p90  = lambda x: np.percentile(x, 90),
            )
            .mul(1000)
            .sort_index()
        )
        all_q = pd.period_range(
            start=raw.index[0], end=raw.index[-1], freq="Q"
        ).strftime("%YQ%q").tolist()
        return raw.reindex(all_q)

    def _extreme_states(df, col):
        peak_q = df.groupby("quarter_label")[col].std().idxmax()
        vals = (
            df.query("quarter_label == @peak_q")
            .assign(b=lambda d: d[col] * 1000)
            .sort_values("b")
        )
        return vals.iloc[0]["state"], vals.iloc[-1]["state"], peak_q

    def _state_ts(df, col, abbrev, all_q):
        s = (
            df[df["state"] == abbrev]
            .set_index("quarter_label")[col]
            .mul(1000)
            .reindex(all_q)
        )
        return [_ql_to_dt(q) for q in s.index], s.values

    dist_d = _build_dist(delta,      "bartik_delta")
    dist_s = _build_dist(instrument, "bartik_s")

    lo_d, hi_d, peak_d = _extreme_states(delta,      "bartik_delta")
    lo_s, hi_s, peak_s = _extreme_states(instrument, "bartik_s")

    all_q_d = dist_d.index.tolist()
    all_q_s = dist_s.index.tolist()

    COLOR_D = "#1f77b4"
    COLOR_S = "#d62728"

    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=False)

    for ax, dist, color, col, df, lo, hi, peak_q, title, ylabel in [
        (
            axes[0], dist_d, COLOR_D, "bartik_delta", delta,
            lo_d, hi_d, peak_d,
            r"$\delta$-instrument  $B^\delta_{s,t}$"
            "\n(BED establishment closings)",
            r"$B^\delta_{s,t}$ (×1 000)",
        ),
        (
            axes[1], dist_s, COLOR_S, "bartik_s", instrument,
            lo_s, hi_s, peak_s,
            r"$s$-instrument  $B^s_{s,t}$"
            "\n(JOLTS total separations)",
            r"$B^s_{s,t}$ (×1 000)",
        ),
    ]:
        valid = dist.dropna(subset=["mean"])
        dts   = [_ql_to_dt(q) for q in valid.index]
        all_q = dist.index.tolist()

        ax.fill_between(dts, valid["p10"], valid["p90"],
                        alpha=0.18, color=color, label="10–90th pctile")
        ax.fill_between(dts, valid["p25"], valid["p75"],
                        alpha=0.32, color=color, label="25–75th pctile (IQR)")
        ax.plot(dts, valid["mean"], color=color, linewidth=2.0,
                label="Cross-state mean")

        ds, vs = _state_ts(df, col, lo, all_q)
        ax.plot(ds, vs, color="firebrick", linewidth=0.85, linestyle=":",
                label=f"{lo} (lowest at {peak_q})")
        ds, vs = _state_ts(df, col, hi, all_q)
        ax.plot(ds, vs, color="darkgreen", linewidth=0.85, linestyle=":",
                label=f"{hi} (highest at {peak_q})")

        ax.axvspan(pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
                   alpha=0.08, color="grey")
        ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
                   alpha=0.10, color="red")

        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.3f}")
        )
        ax.legend(fontsize=8, framealpha=0.85, ncol=1)
        ax.grid(axis="y", linewidth=0.5, alpha=0.4)

    fig.suptitle(
        "Bartik instrument cross-state distributions over time\n"
        "(quarterly, ×1 000 base workers; employment shares fixed at 2006)",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()

    outpath = DEFAULT_OUTPUT_DIR / "instruments_distribution_comparison.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    _open_file(outpath)
    print(f"\nTwo-panel plot saved: {outpath}")