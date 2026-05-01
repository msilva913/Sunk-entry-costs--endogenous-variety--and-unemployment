"""
part3_resid_instruments.py — Bartik instruments from residualized shocks
=========================================================================
Constructs residualized Bartik instruments for all four shock series
produced by part2b_shock_comovement.py:

    B~^δ_{s,t}  = Σ_j  ω_{s,j}  nu_delta_{j,t}    (product-line destruction)
    B~^TS_{s,t} = Σ_j  ω_{s,j}  nu_TS_{j,t}   (total separations)
    B~^LD_{s,t} = Σ_j  ω_{s,j}  nu_LD_{j,t}   (layoffs and discharges)
    B~^QU_{s,t} = Σ_j  ω_{s,j}  nu_QU_{j,t}   (quits)

where ν^k_{j,t} are residuals from the v2 enriched regression in part2b:
    log g^k_{j,t} = α_j + γ_k Δlog p_t + λ_k Δlog VA_{j,t-1} + μ_k log θ_{t-1} + ν^k_{j,t}

    Δlog p_t       = aggregate productivity growth (OPHNFB, quarterly)
    Δlog VA_{j,t-1}= lagged real value-added growth in supersector j (BEA via FRED)
    log θ_{t-1}    = lagged national market tightness log(V/U)

    The VA and tightness regressors are included for all shock types.
    Both are lagged one quarter to ensure predetermination w.r.t. g^k_{j,t}.

The employment shares ω_{s,j} are identical to those used in part3 and
part3_instrument_s (base-year 2006 QCEW shares).

The key diagnostic is the cross-state instrument correlation matrix after
residualization, especially r(B~^δ, B~^LD) and r(B~^δ, B~^QU), which
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
for lbl, df, col in [("nu_delta",  resid_d,  "nu_delta"),
                      ("nu_TS", resid_s,  "nu_s"),
                      ("nu_LD", resid_ld, "nu_ld"),
                      ("nu_QU", resid_qu, "nu_qu")]:
    print(f"  {lbl:<6} {len(df):,} rows  "
          f"({df['quarter_label'].nunique()} quarters, "
          f"{df['industry_code'].nunique()} industries)")

# ── Bartik aggregation ─────────────────────────────────────────────────────
# B~^k_{s,t} = Σ_j  ω_{s,j}  ν^k_{j,t}

def _build_resid_bartik(shares, resid, shock_col, bartik_col,
                        weight_warn_threshold=0.80):
    """
    Aggregate residualized national industry shocks to state-quarter Bartik.

    Parameters
    ----------
    shares                : state_fips, industry_code, share
    resid                 : industry_code, quarter_label, shock_col
    shock_col             : residual column name (e.g. "nu_delta")
    bartik_col            : output instrument column name
    weight_warn_threshold : print warning if weight_sum < this (default 0.80)

    Returns
    -------
    DataFrame: state, state_fips, quarter_label, bartik_col,
               weight_sum, n_supersectors

    Notes
    -----
    Weight coverage can fall below 1.0 when the residual parquet is missing
    data for some industries (e.g. JOLTS suppresses quits for thin-sample
    supersectors like Mining or Information).  The diagnostic block below
    identifies which industries are missing so the caller can decide whether
    coverage is adequate for the intended instrument.
    """
    # ── coverage diagnostic ───────────────────────────────────────────────
    # Report which industries have ANY non-null residuals vs. which are
    # completely absent.  A missing industry lowers weight_sum for every
    # state that has employment in it.
    all_inds   = set(shares["industry_code"].unique())
    resid_inds = set(resid.dropna(subset=[shock_col])["industry_code"].unique())
    missing_inds = all_inds - resid_inds
    if missing_inds:
        # Compute the average employment share of missing industries
        # (approximates how much weight_sum will be depressed on average)
        avg_missing_share = (
            shares[shares["industry_code"].isin(missing_inds)]
            .groupby("state_fips")["share"].sum()
            .mean()
        )
        print(f"    Coverage: {len(resid_inds)}/{len(all_inds)} industries "
              f"have {shock_col} data.")
        print(f"    Missing industries: {sorted(missing_inds)}")
        print(f"    Average missing employment share per state: "
              f"{avg_missing_share:.3f}  "
              f"(→ expected mean weight_sum ≈ {1 - avg_missing_share:.3f})")
    else:
        print(f"    Coverage: all {len(all_inds)} industries present.")

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
        .apply(_agg, include_groups=False)
        .reset_index()
    )

    low = instr["weight_sum"] < weight_warn_threshold
    if low.any():
        print(f"    NOTE: {low.sum()} of {len(instr):,} cells have "
              f"weight_sum < {weight_warn_threshold:.0%} "
              f"(mean weight = {instr['weight_sum'].mean():.3f}, "
              f"min = {instr['weight_sum'].min():.3f}).  "
              f"Interpret instrument with caution if missing industries "
              f"are economically large.")

    instr["state_fips"] = instr["state_fips"].astype(str).str.zfill(2)
    instr["state"]      = instr["state_fips"].map(FIPS2D_TO_STATE)
    instr["state_fips"] = instr["state_fips"].astype(int)

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
    print(f"[Bartik] aggregating B~^{lbl} ...")
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
# Key diagnostic: r(B~^δ, B~^LD) and r(B~^δ, B~^QU) vs. baseline r(B^δ, B^TS).

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
    print(f"    r(B~^{a}, B~^{b}) = {r:+.4f}")

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
    print(f"  Residualized r(B~^δ, B~^TS):         {r_resid_dts:+.4f}")
    print(f"  Residualized r(B~^δ, B~^LD):         "
          f"{all_instr['bartik_delta'].corr(all_instr['bartik_ld']):+.4f}")
    print(f"  Residualized r(B~^δ, B~^QU):         "
          f"{all_instr['bartik_delta'].corr(all_instr['bartik_qu']):+.4f}")

# ── quarterly cross-sectional correlation time series ─────────────────────

def _quarterly_r(df, col1, col2):
    return (
        df.groupby("quarter_label")
        .apply(lambda g: g[col1].corr(g[col2]), include_groups=False)
        .rename("r_cross").reset_index()
        .sort_values("quarter_label")
        .assign(date=lambda d: d["quarter_label"].map(_ql_to_dt))
    )

r_dts  = _quarterly_r(all_instr, "bartik_delta", "bartik_s")
r_dld  = _quarterly_r(all_instr, "bartik_delta", "bartik_ld")
r_dqu  = _quarterly_r(all_instr, "bartik_delta", "bartik_qu")

# ── plot 1: quarterly cross-sectional correlation ─────────────────────────

fig, ax = plt.subplots(figsize=(12, 4))
_ = ax.axhline(0,    color="black", linewidth=0.7)
_ = ax.axhline(0.50, color="grey",  linewidth=0.5, linestyle="--", alpha=0.5,
           label="r = 0.50")
_ = ax.axhline(0.70, color="grey",  linewidth=0.5, linestyle=":",  alpha=0.5,
           label="r = 0.70")

colors = {"δ vs TS": "#1f77b4", "δ vs LD": "#d62728", "δ vs QU": "#2ca02c"}
for r_df, lbl in [(r_dts, "δ vs TS"), (r_dld, "δ vs LD"), (r_dqu, "δ vs QU")]:
    _ = ax.plot(r_df["date"], r_df["r_cross"],
            linewidth=1.6, label=lbl, color=colors[lbl])

# Overlay raw δ vs TS if available.
if RAW_D_PATH.exists() and RAW_S_PATH.exists() and "both_raw" in dir():
    r_raw_qt = _quarterly_r(both_raw, "bartik_delta", "bartik_s")
    _ = ax.plot(r_raw_qt["date"], r_raw_qt["r_cross"],
            color="#1f77b4", linewidth=1.6, linestyle="--", alpha=0.5,
            label="δ vs TS (raw)")

for start, end in REC_SPANS:
    _ = ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
               alpha=0.10, color="grey",
               label="NBER recessions" if start == "2001-03-01" else "_nolegend_")

_ = ax.set_title(
    r"Quarterly cross-sectional correlation with $\tilde{B}^\delta$"
    "\n"
    r"Each point = Pearson $r$ across ~51 states in quarter $t$",
    fontsize=11,
)
_ = ax.set_ylabel("Cross-sectional r", fontsize=10)
_ = ax.set_ylim(-0.55, 1.05)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.2f}"))
_ = ax.legend(fontsize=9, framealpha=0.85, ncol=2)
_ = ax.grid(axis="y", linewidth=0.4, alpha=0.4)
fig.tight_layout()
p = RESULTS_DIR / "instruments_resid_quarterly_correlation.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot 1 saved: {p}")

# ── plot 2: correlation heatmap — delta, LD, QU (the three key instruments)
# This is the core diagnostic: shows whether residualization makes delta
# empirically distinct from the two separation-type instruments.

def _heatmap(ax, mat, title):
    im = ax.imshow(mat.values, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    lbls = mat.columns.tolist()
    _ = ax.set_xticks(range(len(lbls))); ax.set_xticklabels(lbls, fontsize=11)
    _ = ax.set_yticks(range(len(lbls))); ax.set_yticklabels(lbls, fontsize=11)
    _ = ax.set_title(title, fontsize=10)
    for i in range(len(lbls)):
        for j in range(len(lbls)):
            v = mat.values[i, j]
            _ = ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                    fontsize=10, color="white" if abs(v) > 0.6 else "black")
    _ = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# 3x3 matrix: delta, LD, QU only
cols_3  = ["bartik_delta", "bartik_ld", "bartik_qu"]
lbls_3  = [r"$\delta$", "LD", "QU"]
corr_3x3 = all_instr[cols_3].corr()
corr_3x3.index   = lbls_3
corr_3x3.columns = lbls_3

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

_heatmap(axes[0], corr_3x3,
         r"Residualized: $\delta$, LD, QU" + "\n(pooled cross-state, 2001Q1+)")

# Right panel: full 4x4 residualized matrix for reference
corr_resid.index   = [r"$\delta$", "TS", "LD", "QU"]
corr_resid.columns = [r"$\delta$", "TS", "LD", "QU"]
_heatmap(axes[1], corr_resid,
         "All residualized instruments\n(pooled cross-state, 2001Q1+)")

# Annotate with raw r(delta, TS) baseline if available
if RAW_D_PATH.exists() and RAW_S_PATH.exists():
    fig.suptitle(
        f"Residualized instrument correlations\n"
        f"Raw baseline r(B^delta, B^TS) = {r_raw:.3f}  ->  "
        f"r(B~^delta, B~^LD) = {all_instr['bartik_delta'].corr(all_instr['bartik_ld']):+.3f}  "
        f"r(B~^delta, B~^QU) = {all_instr['bartik_delta'].corr(all_instr['bartik_qu']):+.3f}",
        fontsize=10, y=1.02,
    )
else:
    fig.suptitle("Residualized instrument correlations", fontsize=11, y=1.01)

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
        _ = ax.set_xlabel("Raw Bartik", fontsize=10)
        _ = ax.set_ylabel("Residualized Bartik", fontsize=10)
        _ = ax.set_title(f"{title}  (r = {r:.3f})", fontsize=10)
        _ = ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        _ = ax.grid(linewidth=0.4, alpha=0.4)
    fig.suptitle(
        "Raw vs. residualized Bartik instruments (each point = one state-quarter)\n"
        r"v2 spec: $\log g^k_{j,t} = \alpha_j + \gamma_k \Delta\log p_t"
        r"+ \lambda_k \Delta\log VA_{j,t-1} + \mu_k \log\theta_{t-1} + \nu^k_{j,t}$",
        fontsize=10,
    )
    fig.tight_layout()
    p = RESULTS_DIR / "instruments_resid_vs_raw_scatter.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
    print(f"Plot 3 saved: {p}")

# plot 4: residualized instrument cross-state distributions over time 
# Three-panel fan chart (δ, LD, QU) matching the style of the raw comparison
# in part3_instrument_s.py (percentile bands + two extreme-state traces).

def _ql_to_dt_r(q: str):
    """'2005Q1' -> pd.Timestamp at quarter start."""
    y, qq = int(q[:4]), int(q[-1])
    return pd.Timestamp(year=y, month=(qq - 1) * 3 + 1, day=1)

def _build_dist_r(df, col):
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

def _extreme_states_r(df, col):
    """Return (low_state, high_state, peak_quarter) by cross-state SD."""
    peak = df.groupby("quarter_label")[col].std().idxmax()
    vals = df.query("quarter_label == @peak").sort_values(col)
    return vals.iloc[0]["state_fips"], vals.iloc[-1]["state_fips"], peak

def _state_ts_r(df, col, fips, all_q):
    s = df[df["state_fips"] == fips].set_index("quarter_label")[col].reindex(all_q)
    return [_ql_to_dt_r(q) for q in s.index if q is not None], \
           s.dropna().values if s.notna().any() else []

COLORS_R = {"delta": "#1f77b4", "ld": "#2ca02c", "qu": "#9467bd"}

panels_r = [
    (instr_d,  "bartik_delta", COLORS_R["delta"],
     r"$\tilde{B}^\delta_{s,t}$ � firm destruction (residualized)"),
    (instr_ld, "bartik_ld",    COLORS_R["ld"],
     r"$\tilde{B}^{LD}_{s,t}$� layoffs & discharges (residualized)"),
    (instr_qu, "bartik_qu",    COLORS_R["qu"],
     r"$\tilde{B}^{QU}_{s,t}$� quits / placebo (residualized)"),
]

fig, axes = plt.subplots(1, 3, figsize=(21, 5), sharey=False)

for ax, (df, col, color, title) in zip(axes, panels_r):
    dist   = _build_dist_r(df, col)
    lo, hi, pq = _extreme_states_r(df, col)
    valid  = dist.dropna(subset=["mean"])
    dts    = [_ql_to_dt_r(q) for q in valid.index]
    all_q  = dist.index.tolist()

    ax.fill_between(dts, valid["p10"], valid["p90"],
                    alpha=0.18, color=color, label="10-90th pctile")
    ax.fill_between(dts, valid["p25"], valid["p75"],
                    alpha=0.32, color=color, label="IQR")
    ax.plot(dts, valid["mean"], color=color, linewidth=2.0,
            label="Cross-state mean")

    # Extreme-state traces
    xs_lo = [_ql_to_dt_r(q) for q in df[df["state_fips"] == lo]
             .sort_values("quarter_label")["quarter_label"].tolist()]
    ys_lo = df[df["state_fips"] == lo].sort_values("quarter_label")[col].values
    xs_hi = [_ql_to_dt_r(q) for q in df[df["state_fips"] == hi]
             .sort_values("quarter_label")["quarter_label"].tolist()]
    ys_hi = df[df["state_fips"] == hi].sort_values("quarter_label")[col].values

    if len(xs_lo) and len(ys_lo):
        ax.plot(xs_lo, ys_lo, color="firebrick", linewidth=0.85, linestyle=":",
                label=f"FIPS {lo} (low at {pq})")
    if len(xs_hi) and len(ys_hi):
        ax.plot(xs_hi, ys_hi, color="darkgreen", linewidth=0.85, linestyle=":",
                label=f"FIPS {hi} (high at {pq})")

    ax.axvspan(pd.Timestamp("2007-12-01"), pd.Timestamp("2009-06-01"),
               alpha=0.08, color="grey", label="GFC")
    ax.axhline(0, color="black", linewidth=0.6)

    ax.set_title(title, fontsize=10)
    ax.set_ylabel("Bartik instrument value", fontsize=9)
    ax.legend(fontsize=7.5, framealpha=0.85)
    ax.grid(axis="y", linewidth=0.5, alpha=0.4)

fig.suptitle(
    "Residualized Bartik instruments: cross-state distributions over time\n"
    r"(v2 spec: $\log g^k_{j,t}$ residualized on $\Delta\log p_t$, "
    r"$\Delta\log VA_{j,t-1}$, $\log\theta_{t-1}$; shares fixed at 2006)",
    fontsize=11, y=1.02,
)
fig.tight_layout()
p4 = RESULTS_DIR / "instruments_resid_distribution_comparison.png"
fig.savefig(p4, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p4)
print(f"Plot 4 saved: {p4}")

print("\nDone.")

