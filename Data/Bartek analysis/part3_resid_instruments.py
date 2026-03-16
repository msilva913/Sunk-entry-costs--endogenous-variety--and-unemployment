"""
part3_resid_instruments.py — Bartik instruments from productivity-residualized shocks
======================================================================================
Constructs productivity-residualized versions of both Bartik instruments:

    B̃^δ_{s,t} = Σ_j  ω_{s,j,t0}  ν^δ_{j,t}
    B̃^s_{s,t} = Σ_j  ω_{s,j,t0}  ν^s_{j,t}

where ν^δ and ν^s are the residuals from part2b:

    log g^k_{j,t} = α_j + γ_k Δlog p_t + ν^k_{j,t}

The employment shares ω_{s,j,t0} are identical to those used in part3 and
part3_instrument_s.  The residualized instruments are compared with the raw
instruments to assess how much of the cross-state variation is productivity-
driven.

Outputs
-------
    data/instruments/delta_instrument_resid_base{BASE_YEAR}.csv
    data/instruments/s_instrument_resid_base{BASE_YEAR}.csv
    data/results/instruments_resid_vs_raw_scatter.png
    data/results/instruments_resid_quarterly_correlation.png

Schema of output CSVs (identical to raw instrument CSVs):
    state, state_fips, quarter_label, bartik_delta (or bartik_s),
    weight_sum, n_supersectors

Run
---
    python part3_resid_instruments.py

Prerequisites
-------------
    python part1_shares.py
    python part2b_shock_comovement.py
"""

import os
import sys
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

from construct_delta_instrument import (
    BASE_YEAR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    FIPS2D_TO_STATE,
)

def _open_file(path):
    try:
        os.startfile(path)
    except AttributeError:
        pass

def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

# ── paths ──────────────────────────────────────────────────────────────────

RESID_D_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_S_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_s_resid.parquet"
OUT_D_PATH    = DEFAULT_OUTPUT_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
OUT_S_PATH    = DEFAULT_OUTPUT_DIR / f"s_instrument_resid_base{BASE_YEAR}.csv"
RAW_D_PATH    = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
RAW_S_PATH    = DEFAULT_OUTPUT_DIR / f"s_instrument_base{BASE_YEAR}.csv"
RESULTS_DIR   = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── check inputs ───────────────────────────────────────────────────────────

print("=" * 65)
print("Part 3 (resid) — Residualized Bartik instruments")
print("=" * 65)

missing = [p for p in [SHARES_PATH, RESID_D_PATH, RESID_S_PATH]
           if not p.exists()]
if missing:
    for p in missing:
        print(f"  ERROR: {p} not found.")
    sys.exit("Run part1_shares.py and part2b_shock_comovement.py first.")

# ── load inputs ────────────────────────────────────────────────────────────

shares  = pd.read_parquet(SHARES_PATH)
resid_d = pd.read_parquet(RESID_D_PATH)   # columns: industry_code, quarter_label, nu_delta
resid_s = pd.read_parquet(RESID_S_PATH)   # columns: industry_code, quarter_label, nu_s

# Coerce industry_code to same type in all three DataFrames.
for df in [shares, resid_d, resid_s]:
    df["industry_code"] = df["industry_code"].astype(str)

print(f"\nShares:      {len(shares):,} rows  "
      f"({shares['state_fips'].nunique()} states, "
      f"{shares['industry_code'].nunique()} industries)")
print(f"ν^δ series:  {len(resid_d):,} rows  "
      f"({resid_d['quarter_label'].nunique()} quarters, "
      f"{resid_d['industry_code'].nunique()} industries)")
print(f"ν^s series:  {len(resid_s):,} rows  "
      f"({resid_s['quarter_label'].nunique()} quarters, "
      f"{resid_s['industry_code'].nunique()} industries)")

# ── Bartik aggregation ─────────────────────────────────────────────────────
# B̃^k_{s,t} = Σ_j  ω_{s,j}  ν^k_{j,t}
# Merge shares (state × industry) with national residuals (industry × quarter),
# then sum within (state, quarter).

def _build_resid_bartik(shares, resid, shock_col, bartik_col):
    """
    Aggregate residualized national industry shocks to state-quarter Bartik.

    Parameters
    ----------
    shares     : DataFrame with columns state_fips, industry_code, share
    resid      : DataFrame with columns industry_code, quarter_label, shock_col
    shock_col  : name of the residual column in resid (e.g. "nu_delta")
    bartik_col : name to give the resulting instrument column

    Returns
    -------
    DataFrame with columns: state, state_fips, quarter_label,
                            bartik_col, weight_sum, n_supersectors
    """
    merged = (
        resid.dropna(subset=[shock_col])
        .merge(shares[["state_fips", "industry_code", "share"]],
               on="industry_code", how="inner")
    )

    def _agg(grp):
        return pd.Series({
            bartik_col      : (grp["share"] * grp[shock_col]).sum(),
            "weight_sum"    : grp["share"].sum(),
            "n_supersectors": len(grp),
        })

    instrument = (
        merged
        .groupby(["state_fips", "quarter_label"], group_keys=False)
        .apply(_agg)
        .reset_index()
    )

    low = instrument["weight_sum"] < 0.80
    if low.any():
        print(f"  WARNING: {low.sum()} cells have weight coverage < 80%")

    instrument["state_fips"] = instrument["state_fips"].astype(int)
    instrument["state"]      = instrument["state_fips"].map(FIPS2D_TO_STATE)

    return instrument[["state", "state_fips", "quarter_label",
                        bartik_col, "weight_sum", "n_supersectors"]]

print("\n[Bartik] aggregating B̃^δ_{s,t} ...")
instr_d = _build_resid_bartik(shares, resid_d, "nu_delta", "bartik_delta")

print("[Bartik] aggregating B̃^s_{s,t} ...")
instr_s = _build_resid_bartik(shares, resid_s, "nu_s", "bartik_s")

# ── save ───────────────────────────────────────────────────────────────────

instr_d.to_csv(OUT_D_PATH, index=False)
instr_s.to_csv(OUT_S_PATH, index=False)
print(f"\nSaved: {OUT_D_PATH}  ({len(instr_d):,} rows)")
print(f"Saved: {OUT_S_PATH}  ({len(instr_s):,} rows)")

# ── summary ────────────────────────────────────────────────────────────────

for name, df, col in [("δ (resid)", instr_d, "bartik_delta"),
                       ("s (resid)", instr_s, "bartik_s")]:
    print(f"\n{name}:  mean={df[col].mean():.6f}  "
          f"std={df[col].std():.6f}  "
          f"min={df[col].min():.6f}  max={df[col].max():.6f}")
    low = df[df["weight_sum"] < 0.80]
    print(f"  Weight coverage < 80%: {len(low)} cells"
          if not low.empty else "  Weight coverage: all cells ≥ 80%  OK")

# ── cross-state correlation check (residualized vs. each other) ───────────

both = (
    instr_d[["state_fips", "quarter_label", "bartik_delta"]]
    .merge(instr_s[["state_fips", "quarter_label", "bartik_s"]],
           on=["state_fips", "quarter_label"])
    .query("quarter_label >= '2001Q1'")
    .sort_values(["quarter_label", "state_fips"])
)

r_pooled = both["bartik_delta"].corr(both["bartik_s"])
print(f"\n--- Residualized instruments: pooled cross-state correlation ---")
print(f"  r(B̃^δ, B̃^s) = {r_pooled:.4f}  (N = {len(both):,} state-quarters)")

r_qt = (
    both
    .groupby("quarter_label")
    .apply(lambda g: g["bartik_delta"].corr(g["bartik_s"]))
    .rename("r_cross")
    .reset_index()
    .sort_values("quarter_label")
)
r_qt["date"] = r_qt["quarter_label"].map(_ql_to_dt)
print(f"  Quarterly cross-sectional r: "
      f"mean={r_qt['r_cross'].mean():.4f}  "
      f"min={r_qt['r_cross'].min():.4f}  "
      f"max={r_qt['r_cross'].max():.4f}")

# Compare with raw instruments if available.
if RAW_D_PATH.exists() and RAW_S_PATH.exists():
    raw_d = pd.read_csv(RAW_D_PATH)
    raw_s = pd.read_csv(RAW_S_PATH)
    raw_d["state_fips"] = raw_d["state_fips"].astype(int)
    raw_s["state_fips"] = raw_s["state_fips"].astype(int)
    both_raw = (
        raw_d[["state_fips", "quarter_label", "bartik_delta"]]
        .merge(raw_s[["state_fips", "quarter_label", "bartik_s"]],
               on=["state_fips", "quarter_label"])
        .query("quarter_label >= '2001Q1'")
    )
    r_raw = both_raw["bartik_delta"].corr(both_raw["bartik_s"])
    print(f"\n  Raw instruments pooled r:         {r_raw:.4f}")
    print(f"  Residualized instruments pooled r: {r_pooled:.4f}")
    print(f"  Reduction: {r_raw - r_pooled:+.4f}")

# ── plot 1: quarterly cross-sectional correlation ─────────────────────────

fig, ax = plt.subplots(figsize=(11, 4))
ax.axhline(0,    color="black", linewidth=0.8)
ax.axhline(0.50, color="grey",  linewidth=0.6, linestyle="--", alpha=0.6,
           label="r = 0.50")
ax.axhline(0.70, color="grey",  linewidth=0.6, linestyle=":",  alpha=0.6,
           label="r = 0.70")
ax.plot(r_qt["date"], r_qt["r_cross"],
        color="#2ca02c", linewidth=1.8,
        label=r"$r_t(\tilde{B}^\delta,\,\tilde{B}^s)$  residualized")

# Overlay raw quarterly correlation if available.
if RAW_D_PATH.exists() and RAW_S_PATH.exists():
    r_qt_raw = (
        both_raw
        .groupby("quarter_label")
        .apply(lambda g: g["bartik_delta"].corr(g["bartik_s"]))
        .rename("r_cross").reset_index().sort_values("quarter_label")
    )
    r_qt_raw["date"] = r_qt_raw["quarter_label"].map(_ql_to_dt)
    ax.plot(r_qt_raw["date"], r_qt_raw["r_cross"],
            color="#1f77b4", linewidth=1.8, linestyle="--", alpha=0.7,
            label=r"$r_t(B^\delta,\,B^s)$  raw")

for start, end in [("2001-03-01", "2001-11-01"),
                    ("2007-12-01", "2009-06-01"),
                    ("2020-01-01", "2020-07-01")]:
    ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
               alpha=0.10, color="grey",
               label="NBER recessions" if start == "2001-03-01" else "_nolegend_")

ax.set_title(
    r"Quarterly cross-sectional correlation: raw vs. residualized instruments"
    "\n"
    r"Each point = Pearson $r$ across ~51 states in quarter $t$",
    fontsize=11,
)
ax.set_ylabel("Cross-sectional r", fontsize=10)
ax.set_ylim(-0.35, 1.05)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.2f}"))
ax.legend(fontsize=9, framealpha=0.85)
ax.grid(axis="y", linewidth=0.4, alpha=0.4)
fig.tight_layout()
p = RESULTS_DIR / "instruments_resid_quarterly_correlation.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot saved: {p}")

# ── plot 2: raw vs. residualized scatter (δ and s side by side) ───────────

if RAW_D_PATH.exists() and RAW_S_PATH.exists():
    comp_d = (
        raw_d[["state_fips", "quarter_label", "bartik_delta"]]
        .rename(columns={"bartik_delta": "raw"})
        .merge(instr_d[["state_fips", "quarter_label", "bartik_delta"]]
               .rename(columns={"bartik_delta": "resid"}),
               on=["state_fips", "quarter_label"])
    )
    comp_s = (
        raw_s[["state_fips", "quarter_label", "bartik_s"]]
        .rename(columns={"bartik_s": "raw"})
        .merge(instr_s[["state_fips", "quarter_label", "bartik_s"]]
               .rename(columns={"bartik_s": "resid"}),
               on=["state_fips", "quarter_label"])
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, comp, color, title in [
        (axes[0], comp_d, "#1f77b4", r"$\delta$ instrument"),
        (axes[1], comp_s, "#d62728", r"$s$ instrument"),
    ]:
        r = comp["raw"].corr(comp["resid"])
        ax.scatter(comp["raw"], comp["resid"], s=4, alpha=0.3, color=color)
        ax.set_xlabel("Raw Bartik", fontsize=10)
        ax.set_ylabel("Residualized Bartik", fontsize=10)
        ax.set_title(f"{title}  (r = {r:.3f})", fontsize=10)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.grid(linewidth=0.4, alpha=0.4)

    fig.suptitle("Raw vs. productivity-residualized Bartik instruments\n"
                 "(each point = one state-quarter)", fontsize=11)
    fig.tight_layout()
    p = RESULTS_DIR / "instruments_resid_vs_raw_scatter.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
    print(f"Plot saved: {p}")

print("\nDone.")
