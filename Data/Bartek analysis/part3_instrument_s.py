"""
part3_instrument_s.py -- Bartik s-shock instrument  (B^s_{s,t})
================================================================
Loads part1_shares.py and part2_shock_rates_s.py outputs, aggregates
into the state-quarter Bartik s-instrument, and saves:

    data/instruments/s_instrument_base{BASE_YEAR}.csv

Columns: state, state_fips, quarter_label, bartik_s, weight_sum,
         n_supersectors

If delta_instrument_base{BASE_YEAR}.csv is also present:
  - Reports pooled Pearson r(B^s, B^delta) over 2001Q1 onward.
  - Plots quarterly cross-sectional correlation r_t(B^s, B^delta):
    for each quarter, the Pearson r across states between the two
    instruments. This is the key diagnostic for whether the two
    instruments carry genuinely distinct cross-state variation, and
    whether that distinction holds during recession quarters.
  - Plots cross-state distributions of both instruments over time.

Run
---
    python part3_instrument_s.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates_s.py
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

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_s_instrument import (
    BASE_YEAR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    SHOCK_RATES_S_PATH,
    S_INSTRUMENT_PATH,
    build_bartik_instrument_s,
)

# ── helpers ────────────────────────────────────────────────────────────────

def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

def _open_file(path):
    try:
        os.startfile(path)          # Windows only; silently skipped elsewhere
    except AttributeError:
        pass

# ── load inputs ────────────────────────────────────────────────────────────

for path in [SHARES_PATH, SHOCK_RATES_S_PATH]:
    if not path.exists():
        script = "part1_shares.py" if "shares" in path.name else "part2_shock_rates_s.py"
        sys.exit(f"ERROR: {path} not found.  Run {script} first.")

shares        = pd.read_parquet(SHARES_PATH)
shock_rates_s = pd.read_parquet(SHOCK_RATES_S_PATH)
print(f"Loaded shares:        {SHARES_PATH}  ({len(shares):,} rows)")
print(f"Loaded shock rates s: {SHOCK_RATES_S_PATH}  ({len(shock_rates_s):,} rows)")

# ── build instrument ───────────────────────────────────────────────────────

print("\n[Bartik] aggregating B^s_{s,t} ...")
instrument = build_bartik_instrument_s(shares, shock_rates_s)

# JOLTS separations are in thousands; QCEW employment counts are in raw units.
instrument["bartik_s"] = instrument["bartik_s"] * 1000

DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
instrument.to_csv(S_INSTRUMENT_PATH, index=False)
print(f"Saved: {S_INSTRUMENT_PATH}  ({len(instrument):,} rows)")

# ── summary ────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print(f"Part 3s -- Bartik s-instrument  (base year {BASE_YEAR})")
print("=" * 60)

n_states = instrument["state_fips"].nunique()
n_qtrs   = instrument["quarter_label"].nunique()
print(f"\nShape: {instrument.shape}  ({n_states} states x {n_qtrs} quarters)")
print("\nFirst 10 rows:")
print(instrument.head(10).to_string(index=False))

# weight coverage
low_weight = instrument[instrument["weight_sum"] < 0.80]
if low_weight.empty:
    print("\nWeight coverage: all cells >= 80%  OK")
else:
    print(f"\nWARNING: {len(low_weight)} cells have weight coverage < 80%:")
    print(low_weight.to_string(index=False))

# cross-state distribution -- annual summary (avoids printing 90+ quarter rows)
print("\n--- Cross-state distribution (annual averages) ---")
summary = (
    instrument
    .assign(year=instrument["quarter_label"].str[:4])
    .groupby("year")["bartik_s"]
    .agg(mean="mean", std="std",
         p10=lambda x: np.percentile(x, 10),
         p90=lambda x: np.percentile(x, 90))
    .round(6)
)
print(summary.to_string())

# peak dispersion quarter
std_q  = instrument.groupby("quarter_label")["bartik_s"].std()
peak_q = std_q.idxmax()
print(f"\nPeak cross-state dispersion: {peak_q}  (std = {std_q[peak_q]:.6f})")
peak_df = (
    instrument.query("quarter_label == @peak_q")
    .assign(bartik_s=lambda d: d["bartik_s"].round(6))
    .sort_values("bartik_s", ascending=False)
    [["state", "bartik_s", "weight_sum", "n_supersectors"]]
)
print(peak_df.to_string(index=False))

# ── delta instrument -- correlation diagnostics and plots ──────────────────

delta_path = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"

if not delta_path.exists():
    print(f"\n(Delta instrument not found at {delta_path} -- skipping correlation checks.)")
else:
    delta = pd.read_csv(delta_path)

    # Coerce state_fips to int in both frames to avoid object/int64 merge error.
    instrument["state_fips"] = instrument["state_fips"].astype(int)
    delta["state_fips"]      = delta["state_fips"].astype(int)

    both = (
        instrument[["state_fips", "state", "quarter_label", "bartik_s"]]
        .merge(delta[["state_fips", "quarter_label", "bartik_delta"]],
               on=["state_fips", "quarter_label"])
        .query("quarter_label >= '2001Q1'")
        .sort_values(["quarter_label", "state_fips"])
    )

    # ── pooled correlation ─────────────────────────────────────────────────
    r_pooled = both[["bartik_s", "bartik_delta"]].corr().iloc[0, 1]
    print(f"\n--- Correlation with delta instrument (2001Q1 onward) ---")
    print(f"  Pooled Pearson r(B^s, B^delta) = {r_pooled:.4f}  "
          f"(N = {len(both):,} state-quarters)")
    if abs(r_pooled) > 0.80:
        print("  WARNING: high pooled correlation -- instruments may be nearly "
              "collinear in joint IV.")
    else:
        print("  Pooled correlation < 0.80: instruments carry independent "
              "variation.  Suitable for joint IV.")

    # ── quarterly cross-sectional correlation ──────────────────────────────
    # For each quarter t, compute Pearson r across the ~51 states between
    # B^s_{s,t} and B^delta_{s,t}.  Adjacent quarters share no observations.
    # Key diagnostic: if r_t spikes toward 1 during GFC quarters, the two
    # instruments collapse to the same variation precisely when financial
    # amplification results suggest the decomposition matters most.
    r_qt = (
        both
        .groupby("quarter_label")
        .apply(lambda g: g["bartik_s"].corr(g["bartik_delta"]))
        .rename("r_cross")
        .reset_index()
        .sort_values("quarter_label")
    )
    r_qt["date"] = r_qt["quarter_label"].map(_ql_to_dt)

    print(f"\n--- Quarterly cross-sectional r_t(B^s, B^delta) ---")
    print(f"  Mean = {r_qt['r_cross'].mean():.4f}  "
          f"Min = {r_qt['r_cross'].min():.4f}  "
          f"Max = {r_qt['r_cross'].max():.4f}")
    high_r = r_qt[r_qt["r_cross"].abs() > 0.70]
    if not high_r.empty:
        print(f"  Quarters with |r| > 0.70 ({len(high_r)}):")
        print(high_r[["quarter_label", "r_cross"]].to_string(index=False))
    else:
        print("  No quarters with |r| > 0.70.")

    # ── plot 1: quarterly cross-sectional correlation time series ──────────
    fig, ax = plt.subplots(figsize=(11, 4))

    ax.axhline(0,    color="black", linewidth=0.8)
    ax.axhline(0.50, color="grey",  linewidth=0.6, linestyle="--", alpha=0.6,
               label="r = 0.50")
    ax.axhline(0.70, color="grey",  linewidth=0.6, linestyle=":",  alpha=0.6,
               label="r = 0.70 (concern threshold)")

    ax.plot(r_qt["date"], r_qt["r_cross"],
            color="#1f77b4", linewidth=1.8, label=r"$r_t(B^s,\,B^\delta)$")

    # NBER recession shading
    for start, end in [("2001-03-01", "2001-11-01"),
                        ("2007-12-01", "2009-06-01"),
                        ("2020-01-01", "2020-07-01")]:
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
                   alpha=0.10, color="grey",
                   label="NBER recessions" if start == "2001-03-01" else "_nolegend_")

    ax.set_title(
        r"Quarterly cross-sectional correlation: $r_t(B^s,\,B^\delta)$"
        "\n"
        r"Each point = Pearson $r$ across ~51 states in quarter $t$  "
        r"(no rolling window)",
        fontsize=11,
    )
    ax.set_ylabel(r"$r_t(B^s,\,B^\delta)$", fontsize=10)
    ax.set_ylim(-0.35, 1.05)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.legend(fontsize=9, framealpha=0.85)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)
    fig.tight_layout()

    corr_plot = DEFAULT_OUTPUT_DIR / "instruments_quarterly_correlation.png"
    fig.savefig(corr_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    _open_file(corr_plot)
    print(f"\nQuarterly correlation plot saved: {corr_plot}")

    # ── plot 2: cross-state distributions of both instruments ──────────────

    def _build_dist(df, col):
        raw = (
            df.groupby("quarter_label")[col]
            .agg(mean="mean",
                 p10=lambda x: np.percentile(x, 10),
                 p25=lambda x: np.percentile(x, 25),
                 p75=lambda x: np.percentile(x, 75),
                 p90=lambda x: np.percentile(x, 90))
            .sort_index()
        )
        all_q = pd.period_range(raw.index[0], raw.index[-1],
                                freq="Q").strftime("%YQ%q").tolist()
        return raw.reindex(all_q)

    def _extreme_states(df, col):
        peak = df.groupby("quarter_label")[col].std().idxmax()
        vals = df.query("quarter_label == @peak").sort_values(col)
        return vals.iloc[0]["state"], vals.iloc[-1]["state"], peak

    def _state_ts(df, col, abbrev, all_q):
        s = df[df["state"] == abbrev].set_index("quarter_label")[col].reindex(all_q)
        return [_ql_to_dt(q) for q in s.index], s.values

    dist_d = _build_dist(delta,      "bartik_delta")
    dist_s = _build_dist(instrument, "bartik_s")
    lo_d, hi_d, peak_d = _extreme_states(delta,      "bartik_delta")
    lo_s, hi_s, peak_s = _extreme_states(instrument, "bartik_s")

    COLOR_D, COLOR_S = "#1f77b4", "#d62728"
    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=False)

    panels = [
        (axes[0], dist_d, COLOR_D, "bartik_delta", delta,
         lo_d, hi_d, peak_d,
         r"$\delta$-instrument  $B^\delta_{s,t}$" "\n(BED establishment closings)",
         r"$B^\delta_{s,t}$ (quarterly rate)"),
        (axes[1], dist_s, COLOR_S, "bartik_s", instrument,
         lo_s, hi_s, peak_s,
         r"$s$-instrument  $B^s_{s,t}$" "\n(JOLTS total separations)",
         r"$B^s_{s,t}$ (quarterly rate)"),
    ]

    for ax, dist, color, col, df, lo, hi, pq, title, ylabel in panels:
        valid = dist.dropna(subset=["mean"])
        dts   = [_ql_to_dt(q) for q in valid.index]
        all_q = dist.index.tolist()

        ax.fill_between(dts, valid["p10"], valid["p90"],
                        alpha=0.18, color=color, label="10–90th pctile")
        ax.fill_between(dts, valid["p25"], valid["p75"],
                        alpha=0.32, color=color, label="IQR")
        ax.plot(dts, valid["mean"], color=color, linewidth=2.0,
                label="Cross-state mean")

        xs, ys = _state_ts(df, col, lo, all_q)
        ax.plot(xs, ys, color="firebrick", linewidth=0.85, linestyle=":",
                label=f"{lo} (low at {pq})")
        xs, ys = _state_ts(df, col, hi, all_q)
        ax.plot(xs, ys, color="darkgreen", linewidth=0.85, linestyle=":",
                label=f"{hi} (high at {pq})")

        ax.axvspan(pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
                   alpha=0.08, color="grey")
        ax.axvspan(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-07-01"),
                   alpha=0.10, color="red")

        ax.set_title(title, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
        ax.legend(fontsize=8, framealpha=0.85)
        ax.grid(axis="y", linewidth=0.5, alpha=0.4)

    fig.suptitle(
        "Bartik instrument cross-state distributions over time\n"
        "(quarterly rates; employment shares fixed at 2006)",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()

    dist_plot = DEFAULT_OUTPUT_DIR / "instruments_distribution_comparison.png"
    fig.savefig(dist_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    _open_file(dist_plot)
    print(f"Distribution comparison plot saved: {dist_plot}")