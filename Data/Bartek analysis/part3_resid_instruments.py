"""
part3_resid_instruments.py — Bartik instruments from residualized shocks
=========================================================================
Constructs residualized Bartik instruments for all four shock series
produced by part2b_shock_comovement.py:

    B̃^δ_{s,t}  = Σ_j  ω_{s,j}  ν^δ_{j,t}    (product-line destruction)
    B̃^TS_{s,t} = Σ_j  ω_{s,j}  ν^TS_{j,t}   (total separations)
    B̃^LD_{s,t} = Σ_j  ω_{s,j}  ν^LD_{j,t}   (layoffs and discharges)
    B̃^QU_{s,t} = Σ_j  ω_{s,j}  ν^QU_{j,t}   (quits)

where ν^k_{j,t} are residuals from shock-specific regressions (see part2b):
    log g^δ_{j,t}  = α_j + γ_δ  Δlog p_t                    + ν^δ_{j,t}
    log g^TS_{j,t} = α_j + γ_TS Δlog p_t                    + ν^TS_{j,t}
    log g^LD_{j,t} = α_j + γ_LD Δlog p_t + λ_LD log θ_t^nat + ν^LD_{j,t}
    log g^QU_{j,t} = α_j + γ_QU Δlog p_t + λ_QU log θ_t^nat + ν^QU_{j,t}

The employment shares ω_{s,j} are identical to those used in part3 and
part3_instrument_s (base-year 2006 QCEW shares).

The key diagnostic is the cross-state instrument correlation matrix after
residualization, especially r(B̃^δ, B̃^LD) and r(B̃^δ, B̃^QU), which
should be substantially lower than the raw r(B^δ, B^TS) = 0.705 if the
decomposition restores empirical distinctness of the two shock types.

Outputs
-------
    data/instruments/delta_instrument_resid_base{BASE_YEAR}.csv
    data/instruments/s_instrument_resid_base{BASE_YEAR}.csv
    data/instruments/ld_instrument_resid_base{BASE_YEAR}.csv
    data/instruments/qu_instrument_resid_base{BASE_YEAR}.csv
    data/results/instruments_resid_quarterly_correlation.png
    data/results/instruments_resid_vs_raw_scatter.png
    data/results/instruments_resid_corr_matrix.png

Run
---
    python part3_resid_instruments.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates_s.py   (produces LD and QU parquets)
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

REC_SPANS = [("2001-03-01", "2001-11-01"),
             ("2007-12-01", "2009-06-01"),
             ("2020-01-01", "2020-07-01")]

# ── paths ──────────────────────────────────────────────────────────────────

RESID_D_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_S_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_s_resid.parquet"
RESID_LD_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_ld_resid.parquet"
RESID_QU_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_qu_resid.parquet"

OUT_D_PATH    = DEFAULT_OUTPUT_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
OUT_S_PATH    = DEFAULT_OUTPUT_DIR / f"s_instrument_resid_base{BASE_YEAR}.csv"
OUT_LD_PATH   = DEFAULT_OUTPUT_DIR / f"ld_instrument_resid_base{BASE_YEAR}.csv"
OUT_QU_PATH   = DEFAULT_OUTPUT_DIR / f"qu_instrument_resid_base{BASE_YEAR}.csv"

RAW_D_PATH    = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
RAW_S_PATH    = DEFAULT_OUTPUT_DIR / f"s_instrument_base{BASE_YEAR}.csv"

RESULTS_DIR   = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── check inputs ───────────────────────────────────────────────────────────

print("=" * 65)
print("Part 3 (resid) — Residualized Bartik instruments (δ, TS, LD, QU)")
print("=" * 65)

required = [SHARES_PATH, RESID_D_PATH, RESID_S_PATH, RESID_LD_PATH, RESID_QU_PATH]
missing  = [p for p in required if not p.exists()]
if missing:
    msg = "Missing required files:\n" + "\n".join(f"  {p}" for p in missing)
    msg += ("\nRun part2_shock_rates_s.py then part2b_shock_comovement.py "
            "to generate them.")
    raise FileNotFoundError(msg)

# ── load inputs ────────────────────────────────────────────────────────────

shares   = pd.read_parquet(SHARES_PATH)
resid_d  = pd.read_parquet(RESID_D_PATH)    # nu_delta
resid_s  = pd.read_parquet(RESID_S_PATH)    # nu_s
resid_ld = pd.read_parquet(RESID_LD_PATH)   # nu_ld
resid_qu = pd.read_parquet(RESID_QU_PATH)   # nu_qu

for df in [shares, resid_d, resid_s, resid_ld, resid_qu]:
    df["industry_code"] = df["industry_code"].astype(str)

print(f"\nShares:       {len(shares):,} rows  "
      f"({shares['state_fips'].nunique()} states, "
      f"{shares['industry_code'].nunique()} industries)")
for lbl, df, col in [("ν^δ",  resid_d,  "nu_delta"),
                      ("ν^TS", resid_s,  "nu_s"),
                      ("ν^LD", resid_ld, "nu_ld"),
                      ("ν^QU", resid_qu, "nu_qu")]:
    print(f"  {lbl:<6} {len(df):,} rows  "
          f"({df['quarter_label'].nunique()} quarters, "
          f"{df['industry_code'].nunique()} industries)")

# ── Bartik aggregation ─────────────────────────────────────────────────────
# B̃^k_{s,t} = Σ_j  ω_{s,j}  ν^k_{j,t}

def _build_resid_bartik(shares, resid, shock_col, bartik_col):
    """
    Aggregate residualized national industry shocks to state-quarter Bartik.

    Parameters
    ----------
    shares     : state_fips, industry_code, share
    resid      : industry_code, quarter_label, shock_col
    shock_col  : residual column name (e.g. "nu_delta")
    bartik_col : output instrument column name

    Returns
    -------
    DataFrame: state, state_fips, quarter_label, bartik_col,
               weight_sum, n_supersectors
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

    instr = (
        merged
        .groupby(["state_fips", "quarter_label"], group_keys=False)
        .apply(_agg)
        .reset_index()
    )

    low = instr["weight_sum"] < 0.80
    if low.any():
        print(f"    WARNING: {low.sum()} cells with weight coverage < 80%")

    instr["state_fips"] = instr["state_fips"].astype(int)
    instr["state"]      = instr["state_fips"].map(FIPS2D_TO_STATE)

    return instr[["state", "state_fips", "quarter_label",
                  bartik_col, "weight_sum", "n_supersectors"]]

print()
instruments = {}
for lbl, resid, shock_col, bartik_col, out_path in [
    ("δ  (resid)", resid_d,  "nu_delta", "bartik_delta", OUT_D_PATH),
    ("TS (resid)", resid_s,  "nu_s",     "bartik_s",     OUT_S_PATH),
    ("LD (resid)", resid_ld, "nu_ld",    "bartik_ld",    OUT_LD_PATH),
    ("QU (resid)", resid_qu, "nu_qu",    "bartik_qu",    OUT_QU_PATH),
]:
    print(f"[Bartik] aggregating B̃^{lbl} ...")
    instr = _build_resid_bartik(shares, resid, shock_col, bartik_col)
    instr.to_csv(out_path, index=False)
    instruments[lbl] = (instr, bartik_col)
    v = instr[bartik_col]
    print(f"  Saved: {out_path}  ({len(instr):,} rows)  "
          f"mean={v.mean():.6f}  std={v.std():.6f}")

instr_d  = instruments["δ  (resid)"][0]
instr_s  = instruments["TS (resid)"][0]
instr_ld = instruments["LD (resid)"][0]
instr_qu = instruments["QU (resid)"][0]

# ── cross-state correlation matrix ────────────────────────────────────────
# Merge all four residualized instruments on (state_fips, quarter_label).
# Key diagnostic: r(B̃^δ, B̃^LD) and r(B̃^δ, B̃^QU) vs. baseline r(B^δ, B^TS).

print("\n--- Residualized instrument correlation matrix (pooled, 2001Q1+) ---")

all_instr = (
    instr_d[["state_fips", "quarter_label", "bartik_delta"]]
    .merge(instr_s[["state_fips",  "quarter_label", "bartik_s"]],
           on=["state_fips", "quarter_label"])
    .merge(instr_ld[["state_fips", "quarter_label", "bartik_ld"]],
           on=["state_fips", "quarter_label"])
    .merge(instr_qu[["state_fips", "quarter_label", "bartik_qu"]],
           on=["state_fips", "quarter_label"])
    .query("quarter_label >= '2001Q1'")
)

cols  = ["bartik_delta", "bartik_s", "bartik_ld", "bartik_qu"]
labels = ["δ", "TS", "LD", "QU"]
corr_resid = all_instr[cols].corr()
corr_resid.index   = labels
corr_resid.columns = labels
print(corr_resid.round(4).to_string())

# Print the key pairs explicitly.
print("\n  Key pairs:")
for (a, ca), (b, cb) in [
    (("δ", "bartik_delta"), ("TS", "bartik_s")),
    (("δ", "bartik_delta"), ("LD", "bartik_ld")),
    (("δ", "bartik_delta"), ("QU", "bartik_qu")),
    (("LD","bartik_ld"),    ("QU", "bartik_qu")),
]:
    r = all_instr[ca].corr(all_instr[cb])
    print(f"    r(B̃^{a}, B̃^{b}) = {r:+.4f}")

# Compare with raw δ vs TS if available.
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
    r_resid_dts = all_instr["bartik_delta"].corr(all_instr["bartik_s"])
    print(f"\n  Baseline raw r(B^δ, B^TS):          {r_raw:+.4f}")
    print(f"  Residualized r(B̃^δ, B̃^TS):         {r_resid_dts:+.4f}")
    print(f"  Residualized r(B̃^δ, B̃^LD):         "
          f"{all_instr['bartik_delta'].corr(all_instr['bartik_ld']):+.4f}")
    print(f"  Residualized r(B̃^δ, B̃^QU):         "
          f"{all_instr['bartik_delta'].corr(all_instr['bartik_qu']):+.4f}")

# ── quarterly cross-sectional correlation time series ─────────────────────

def _quarterly_r(df, col1, col2):
    return (
        df.groupby("quarter_label")
        .apply(lambda g: g[col1].corr(g[col2]))
        .rename("r_cross").reset_index()
        .sort_values("quarter_label")
        .assign(date=lambda d: d["quarter_label"].map(_ql_to_dt))
    )

r_dts  = _quarterly_r(all_instr, "bartik_delta", "bartik_s")
r_dld  = _quarterly_r(all_instr, "bartik_delta", "bartik_ld")
r_dqu  = _quarterly_r(all_instr, "bartik_delta", "bartik_qu")

# ── plot 1: quarterly cross-sectional correlation ─────────────────────────

fig, ax = plt.subplots(figsize=(12, 4))
ax.axhline(0,    color="black", linewidth=0.7)
ax.axhline(0.50, color="grey",  linewidth=0.5, linestyle="--", alpha=0.5,
           label="r = 0.50")
ax.axhline(0.70, color="grey",  linewidth=0.5, linestyle=":",  alpha=0.5,
           label="r = 0.70")

colors = {"δ vs TS": "#1f77b4", "δ vs LD": "#d62728", "δ vs QU": "#2ca02c"}
for r_df, lbl in [(r_dts, "δ vs TS"), (r_dld, "δ vs LD"), (r_dqu, "δ vs QU")]:
    ax.plot(r_df["date"], r_df["r_cross"],
            linewidth=1.6, label=lbl, color=colors[lbl])

# Overlay raw δ vs TS if available.
if RAW_D_PATH.exists() and RAW_S_PATH.exists() and "both_raw" in dir():
    r_raw_qt = _quarterly_r(both_raw, "bartik_delta", "bartik_s")
    ax.plot(r_raw_qt["date"], r_raw_qt["r_cross"],
            color="#1f77b4", linewidth=1.6, linestyle="--", alpha=0.5,
            label="δ vs TS (raw)")

for start, end in REC_SPANS:
    ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
               alpha=0.10, color="grey",
               label="NBER recessions" if start == "2001-03-01" else "_nolegend_")

ax.set_title(
    r"Quarterly cross-sectional correlation with $\tilde{B}^\delta$"
    "\n"
    r"Each point = Pearson $r$ across ~51 states in quarter $t$",
    fontsize=11,
)
ax.set_ylabel("Cross-sectional r", fontsize=10)
ax.set_ylim(-0.55, 1.05)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.2f}"))
ax.legend(fontsize=9, framealpha=0.85, ncol=2)
ax.grid(axis="y", linewidth=0.4, alpha=0.4)
fig.tight_layout()
p = RESULTS_DIR / "instruments_resid_quarterly_correlation.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot 1 saved: {p}")

# ── plot 2: pooled correlation heatmap ────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

def _heatmap(ax, mat, title):
    im = ax.imshow(mat.values, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    lbls = mat.columns.tolist()
    ax.set_xticks(range(len(lbls))); ax.set_xticklabels(lbls, fontsize=10)
    ax.set_yticks(range(len(lbls))); ax.set_yticklabels(lbls, fontsize=10)
    ax.set_title(title, fontsize=10)
    for i in range(len(lbls)):
        for j in range(len(lbls)):
            v = mat.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=9, color="white" if abs(v) > 0.6 else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

_heatmap(axes[0], corr_resid, "Residualized instruments\n(pooled, 2001Q1+)")

# Raw correlation matrix for comparison (δ and TS only; LD/QU not available raw).
if RAW_D_PATH.exists() and RAW_S_PATH.exists():
    corr_raw_mat = pd.DataFrame(
        [[1.0, r_raw], [r_raw, 1.0]],
        index=["δ", "TS"], columns=["δ", "TS"]
    )
    _heatmap(axes[1], corr_raw_mat,
             f"Raw instruments (δ vs TS)\nr = {r_raw:.3f}")
else:
    axes[1].axis("off")
    axes[1].text(0.5, 0.5, "Raw instruments\nnot available",
                 ha="center", va="center", fontsize=10)

fig.suptitle("Cross-state instrument correlation: raw vs. residualized",
             fontsize=11, y=1.01)
fig.tight_layout()
p = RESULTS_DIR / "instruments_resid_corr_matrix.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot 2 saved: {p}")

# ── plot 3: raw vs. residualized scatter for δ and TS ─────────────────────

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
        (axes[1], comp_s, "#d62728", r"$s$ (TS) instrument"),
    ]:
        r = comp["raw"].corr(comp["resid"])
        ax.scatter(comp["raw"], comp["resid"], s=4, alpha=0.3, color=color)
        ax.set_xlabel("Raw Bartik", fontsize=10)
        ax.set_ylabel("Residualized Bartik", fontsize=10)
        ax.set_title(f"{title}  (r = {r:.3f})", fontsize=10)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.grid(linewidth=0.4, alpha=0.4)
    fig.suptitle("Raw vs. residualized Bartik instruments\n"
                 r"($\delta$/TS: $\Delta\log p$ only — LD/QU: $\Delta\log p + \log\theta^{nat}$)"
                 "\n"
                 "(each point = one state-quarter)", fontsize=11)
    fig.tight_layout()
    p = RESULTS_DIR / "instruments_resid_vs_raw_scatter.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
    print(f"Plot 3 saved: {p}")

print("\nDone.")

