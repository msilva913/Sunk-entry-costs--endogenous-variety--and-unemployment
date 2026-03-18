"""
part2b_shock_comovement.py — Shock-productivity decomposition and comovement
=============================================================================
Diagnostic step between part2 (shock-rate construction) and part3 (Bartik
aggregation).

Purpose
-------
In the model, δ and s are exogenous AR(1) processes independent of technology
z.  In the data, both series contain an endogenous component driven by
aggregate productivity: firms close and workers separate partly because
productivity falls.  This file asks two questions:

  1. How much do g^δ and g^s co-move across industries and over time?
     High co-movement implies the Bartik instruments B^δ and B^s will be
     nearly collinear, undermining separate identification.

  2. How much of that co-movement is explained by aggregate productivity?
     We follow Coles-Kelishomi (2018), who find corr(p_t, δ_t) = −0.63,
     and residualize each series on productivity to isolate the idiosyncratic
     component.  The regressor is Δlog p_t (quarterly log-growth) rather
     than log-levels, which have a trend and would produce spurious
     coefficients.

Regression:
    log g^k_{j,t} = α_j + γ_k Δlog p_t + ν^k_{j,t},   k ∈ {δ, s}

The residuals ν^δ and ν^s are the closest empirical analog to the truly
exogenous shocks in the model.

Outputs
-------
    data/instruments/shock_rates_delta_resid.parquet   residualized δ series
    data/instruments/shock_rates_s_resid.parquet       residualized s series
    data/results/comovement_raw_series.png             raw series by industry
    data/results/comovement_correlation_heatmaps.png   4-panel heatmap
    data/results/comovement_scatter.png                cross-series scatter

Run
---
    python part2b_shock_comovement.py

Prerequisites
-------------
    python part2_shock_rates.py
    python part2_shock_rates_s.py
    Internet access on first run (fetches OPHNFB from FRED; cached thereafter)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from fredapi import Fred

# ── working directory ──────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import (
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHOCK_RATES_PATH,
    INDUSTRY_LABELS,
)
from construct_s_instrument import (
    SHOCK_RATES_S_PATH,
    SHOCK_RATES_LD_PATH,
    SHOCK_RATES_QU_PATH,
)

# ── paths ──────────────────────────────────────────────────────────────────
PROD_CACHE    = DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet"
RESID_D_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_S_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_s_resid.parquet"
RESID_LD_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_ld_resid.parquet"
RESID_QU_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_qu_resid.parquet"
RESULTS_DIR  = Path("data/results")

DEFAULT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────
def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

def _dt_to_ql(dt):
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"

def _open_file(path):
    try:
        os.startfile(path)
    except AttributeError:
        pass

REC_SPANS = [("2001-03-01", "2001-11-01"),
             ("2007-12-01", "2009-06-01"),
             ("2020-01-01", "2020-07-01")]

def _add_recessions(ax):
    for start, end in REC_SPANS:
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end),
                   alpha=0.10, color="grey")

# ── 1. check inputs ────────────────────────────────────────────────────────
print("=" * 65)
print("Part 2b — Shock-productivity decomposition and comovement")
print("=" * 65)

_missing = [p for p in [SHOCK_RATES_PATH, SHOCK_RATES_S_PATH,
                        SHOCK_RATES_LD_PATH, SHOCK_RATES_QU_PATH]
            if not p.exists()]
if _missing:
    msg = "Missing required shock rate files:\n"
    msg += "\n".join(f"  {p}" for p in _missing)
    msg += "\nRun part2_shock_rates.py and part2_shock_rates_s.py first."
    raise FileNotFoundError(msg)

# ── 2. national industry series ────────────────────────────────────────────
# Average across states to recover the national series.
# The LOO correction only affects the denominator and is negligible here.

raw_d  = pd.read_parquet(SHOCK_RATES_PATH)
raw_s  = pd.read_parquet(SHOCK_RATES_S_PATH)
raw_ld = pd.read_parquet(SHOCK_RATES_LD_PATH)
raw_qu = pd.read_parquet(SHOCK_RATES_QU_PATH)

def _nat_series(raw, col, new_col):
    return (raw.groupby(["industry_code", "quarter_label"])[col]
               .mean().reset_index()
               .rename(columns={col: new_col}))

nat_d  = _nat_series(raw_d,  "g_delta_loo", "g_delta")
nat_s  = _nat_series(raw_s,  "g_s_loo",     "g_s")
nat_ld = _nat_series(raw_ld, "g_ld_loo",    "g_ld")
nat_qu = _nat_series(raw_qu, "g_qu_loo",    "g_qu")

for lbl, df in [("δ",     nat_d), ("s (TS)", nat_s),
                ("s (LD)", nat_ld), ("s (QU)", nat_qu)]:
    print(f"  {lbl}: {df['quarter_label'].nunique()} quarters, "
          f"{df['industry_code'].nunique()} industries")

# Restrict to overlapping sample.
nat = nat_d.merge(nat_s,  on=["industry_code", "quarter_label"], how="inner")
nat = nat.merge(nat_ld, on=["industry_code", "quarter_label"], how="inner")
nat = nat.merge(nat_qu, on=["industry_code", "quarter_label"], how="inner")
print(f"Overlapping sample: {nat['quarter_label'].min()} – "
      f"{nat['quarter_label'].max()}  "
      f"({nat['quarter_label'].nunique()} quarters, "
      f"{nat['industry_code'].nunique()} industries, "
      f"{len(nat):,} obs)")

# Remove zeros/negatives (log undefined); attach labels.
pos = ((nat["g_delta"] > 0) & (nat["g_s"] > 0) &
       (nat["g_ld"] > 0) & (nat["g_qu"] > 0))
nat = nat[pos].copy()
nat["log_g_delta"]    = np.log(nat["g_delta"])
nat["log_g_s"]        = np.log(nat["g_s"])
nat["log_g_ld"]       = np.log(nat["g_ld"])
nat["log_g_qu"]       = np.log(nat["g_qu"])
nat["industry_label"] = nat["industry_code"].map(INDUSTRY_LABELS).fillna(
                            nat["industry_code"].astype(str))
nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)

# ── 3. aggregate productivity ──────────────────────────────────────────────
# OPHNFB: BLS output per hour, nonfarm business, quarterly SA (index 2017=100).
# We use Δlog p_t (quarterly log-growth) rather than log-levels to avoid
# regressing a stationary shock rate on a trending series.

if PROD_CACHE.exists():
    prod = pd.read_parquet(PROD_CACHE)
    print(f"\nProductivity: loaded from cache ({PROD_CACHE})")
else:
    print("\nFetching OPHNFB from FRED via fredapi ...")
    api_key = os.environ.get("FRED_API_KEY", "")
    if not api_key:
        sys.exit("ERROR: set FRED_API_KEY environment variable before running.")
    fred  = Fred(api_key=api_key)
    s     = fred.get_series("OPHNFB")          # quarterly SA index, 2017=100
    prod  = (s.rename("productivity")
              .reset_index()
              .rename(columns={"index": "date"}))
    prod["quarter_label"] = prod["date"].map(_dt_to_ql)
    prod  = prod[["quarter_label", "productivity"]].dropna().copy()
    prod.to_parquet(PROD_CACHE, index=False)
    print(f"Saved: {PROD_CACHE}")

prod = prod.sort_values("quarter_label").reset_index(drop=True)
prod["log_p"]      = np.log(prod["productivity"])
prod["dlog_p"]     = prod["log_p"].diff()    # quarterly log-growth (primary)
prod["dlog_p_yoy"] = prod["log_p"].diff(4)   # year-over-year (robustness)
prod = prod.dropna(subset=["dlog_p"])

nat = nat.merge(prod[["quarter_label", "dlog_p", "dlog_p_yoy"]],
                on="quarter_label", how="inner")
print(f"After merging productivity: {len(nat):,} industry-quarter obs")

# ── 4. residualization ─────────────────────────────────────────────────────
# Regression:  log g^k_{j,t} = α_j + γ_k Δlog p_t + ν^k_{j,t}
# Industry FEs absorbed via within-transformation (subtract industry means).
#
# Returns a DataFrame keyed by (industry_code, quarter_label) with the
# residual as a named column.  We merge back into nat rather than assigning
# into it directly, which is the most robust approach across pandas versions.

def _residualize(df, dep_col, resid_col, prod_col="dlog_p"):
    """
    OLS of dep_col on industry FEs + prod_col.
    Prints γ and within-R².

    df must be sorted by [industry_code, quarter_label] with a clean
    0-based index (call nat.sort_values(...).reset_index(drop=True) first).
    Returns a DataFrame with the same row order as df and a column named
    resid_col containing the residuals as a numpy-backed float64 series.
    """
    # Work on a fresh 0-based copy in the same row order as df.
    work = df[["industry_code", dep_col, prod_col]].copy().reset_index(drop=True)

    # Within-transform: subtract each industry's mean.
    work["y_dm"] = (work[dep_col]
                    - work.groupby("industry_code")[dep_col].transform("mean"))
    work["p_dm"] = (work[prod_col]
                    - work.groupby("industry_code")[prod_col].transform("mean"))

    gamma  = ((work["y_dm"] * work["p_dm"]).sum()
              / (work["p_dm"] ** 2).sum())
    resid  = (work["y_dm"] - gamma * work["p_dm"]).to_numpy()   # plain numpy array
    ss_tot = (work["y_dm"] ** 2).sum()
    r2     = 1 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else np.nan

    print(f"  {resid_col:<10}  γ = {gamma:+.4f}   R²(within) = {r2:.4f}")

    out = df[["industry_code", "quarter_label"]].copy().reset_index(drop=True)
    out[resid_col] = resid
    return out

# Sort nat by [industry_code, quarter_label] and reset to a clean 0-based
# index before residualizing.  _residualize works on a copy in the same
# row order, so returned arrays align positionally with nat's rows.
nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)

print(f"\n--- Productivity regression  (\u0394log p_t, quarterly growth) ---")
resid_d  = _residualize(nat, "log_g_delta", "nu_delta")
resid_s  = _residualize(nat, "log_g_s",     "nu_s")
resid_ld = _residualize(nat, "log_g_ld",    "nu_ld")
resid_qu = _residualize(nat, "log_g_qu",    "nu_qu")

nat["nu_delta"] = resid_d["nu_delta"].to_numpy()
nat["nu_s"]     = resid_s["nu_s"].to_numpy()
nat["nu_ld"]    = resid_ld["nu_ld"].to_numpy()
nat["nu_qu"]    = resid_qu["nu_qu"].to_numpy()

for col in ["nu_delta", "nu_s", "nu_ld", "nu_qu"]:
    n_null = nat[col].isna().sum()
    status = f"WARNING: {n_null} nulls" if n_null else f"OK ({len(nat):,} obs)"
    print(f"  {col}: {status}")

# ── 5. comovement diagnostics ──────────────────────────────────────────────

def _pooled_r(df, col1, col2):
    """Pearson r across all industry-quarter observations."""
    return df[col1].corr(df[col2])

def _by_industry_r(df, col1, col2):
    """Per-industry time-series Pearson r between col1 and col2."""
    return (
        df.groupby("industry_label", group_keys=False)
          .apply(lambda g: pd.Series({"r": g[col1].corr(g[col2])}))
          .reset_index()
          .sort_values("r")
    )

# ── Pooled pairwise correlations ──────────────────────────────────────────
pairs = [
    ("δ vs TS",  "log_g_delta", "log_g_s",  "nu_delta", "nu_s"),
    ("δ vs LD",  "log_g_delta", "log_g_ld", "nu_delta", "nu_ld"),
    ("δ vs QU",  "log_g_delta", "log_g_qu", "nu_delta", "nu_qu"),
    ("LD vs QU", "log_g_ld",    "log_g_qu", "nu_ld",    "nu_qu"),
    ("TS vs LD", "log_g_s",     "log_g_ld", "nu_s",     "nu_ld"),
    ("TS vs QU", "log_g_s",     "log_g_qu", "nu_s",     "nu_qu"),
]

print(f"\n--- Pooled cross-series correlations ---")
print(f"  {'Pair':<12}  {'r raw':>8}  {'r resid':>8}  {'reduction':>10}")
print(f"  {'-'*46}")
for label, rc1, rc2, rc3, rc4 in pairs:
    r_raw_p   = _pooled_r(nat, rc1, rc2)
    r_resid_p = _pooled_r(nat, rc3, rc4)
    print(f"  {label:<12}  {r_raw_p:>+8.4f}  {r_resid_p:>+8.4f}  {r_raw_p-r_resid_p:>+10.4f}")

print(f"\n--- Per-industry corr(δ vs LD) and corr(δ vs QU) ---")
by_dld_raw   = _by_industry_r(nat, "log_g_delta", "log_g_ld").rename(columns={"r": "r_dLD_raw"})
by_dld_resid = _by_industry_r(nat, "nu_delta",    "nu_ld"   ).rename(columns={"r": "r_dLD_resid"})
by_dqu_raw   = _by_industry_r(nat, "log_g_delta", "log_g_qu").rename(columns={"r": "r_dQU_raw"})
by_dqu_resid = _by_industry_r(nat, "nu_delta",    "nu_qu"   ).rename(columns={"r": "r_dQU_resid"})
by_ind = (by_dld_raw
          .merge(by_dld_resid, on="industry_label")
          .merge(by_dqu_raw,   on="industry_label")
          .merge(by_dqu_resid, on="industry_label")
          .sort_values("r_dLD_raw"))
print(by_ind.to_string(index=False))

# Keep for plot titles
r_raw   = _pooled_r(nat, "log_g_delta", "log_g_s")
r_resid = _pooled_r(nat, "nu_delta",    "nu_s")

# ── 6. save residualized series ────────────────────────────────────────────

nat[["industry_code", "quarter_label", "nu_delta"]].to_parquet(RESID_D_PATH,  index=False)
nat[["industry_code", "quarter_label", "nu_s"]    ].to_parquet(RESID_S_PATH,  index=False)
nat[["industry_code", "quarter_label", "nu_ld"]   ].to_parquet(RESID_LD_PATH, index=False)
nat[["industry_code", "quarter_label", "nu_qu"]   ].to_parquet(RESID_QU_PATH, index=False)
for p in [RESID_D_PATH, RESID_S_PATH, RESID_LD_PATH, RESID_QU_PATH]:
    print(f"Saved: {p}")

# ── 7. plots ───────────────────────────────────────────────────────────────

industries = sorted(nat["industry_label"].unique())
quarters   = sorted(nat["quarter_label"].unique())
dates      = [_ql_to_dt(q) for q in quarters]
n_ind      = len(industries)

# ── plot A: raw log series by industry, all four shock types ──────────────
# Four columns: δ, TS, LD, QU
series_cols = [
    ("log_g_delta", r"$\log g^\delta$",  "#1f77b4"),
    ("log_g_s",     r"$\log g^{TS}$",     "#d62728"),
    ("log_g_ld",    r"$\log g^{LD}$",     "#2ca02c"),
    ("log_g_qu",    r"$\log g^{QU}$",     "#ff7f0e"),
]
fig, axes = plt.subplots(n_ind, 4,
                         figsize=(20, 2.0 * n_ind),
                         sharex=True, squeeze=False)
for i, ind in enumerate(industries):
    sub = (nat[nat["industry_label"] == ind]
           .set_index("quarter_label")
           .reindex(quarters))
    for j, (col, _, color) in enumerate(series_cols):
        ax = axes[i, j]
        ax.plot(dates, sub[col].values, color=color, linewidth=1.1)
        _add_recessions(ax)
        ax.set_ylabel(ind if j == 0 else "", fontsize=7)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
for j, (_, title, _) in enumerate(series_cols):
    axes[0, j].set_title(title, fontsize=10)
fig.suptitle("National industry shock rates by supersector (raw log)",
             fontsize=11, y=1.005)
fig.tight_layout()
p = RESULTS_DIR / "comovement_raw_series.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot A saved: {p}")

# ── plot B: cross-industry correlation heatmaps, 2×4 grid ─────────────────
# Row 1: raw series; Row 2: residualized
def _corr_matrix(df, col):
    """Cross-industry Pearson r matrix (industries × industries)."""
    wide = df.pivot_table(index="quarter_label",
                          columns="industry_label",
                          values=col, aggfunc="mean")
    return wide.corr()

hmap_panels = [
    (r"$\log g^\delta$ raw",  "log_g_delta"),
    (r"$\log g^{TS}$ raw",     "log_g_s"),
    (r"$\log g^{LD}$ raw",     "log_g_ld"),
    (r"$\log g^{QU}$ raw",     "log_g_qu"),
    (r"$\nu^\delta$ resid",   "nu_delta"),
    (r"$\nu^{TS}$ resid",      "nu_s"),
    (r"$\nu^{LD}$ resid",      "nu_ld"),
    (r"$\nu^{QU}$ resid",      "nu_qu"),
]
fig, axes = plt.subplots(2, 4, figsize=(22, 11))
for idx, (title, col) in enumerate(hmap_panels):
    ax  = axes[idx // 4, idx % 4]
    mat = _corr_matrix(nat, col)
    im  = ax.imshow(mat.values, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    labels = mat.columns.tolist()
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title(title, fontsize=9)
    for ii in range(len(labels)):
        for jj in range(len(labels)):
            v = mat.values[ii, jj]
            ax.text(jj, ii, f"{v:.2f}", ha="center", va="center",
                    fontsize=5, color="white" if abs(v) > 0.7 else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
fig.suptitle(
    "Cross-industry correlation matrices: raw vs. residualized (all four shock types)",
    fontsize=11, y=1.01,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_correlation_heatmaps.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot B saved: {p}")

# ── plot C: δ vs each s-type scatter, raw and residualized ────────────────
# Three pairs: δ vs TS, δ vs LD, δ vs QU  (2 rows × 3 cols)
scatter_pairs = [
    ("δ vs TS", "log_g_delta", "log_g_s",  "nu_delta", "nu_s",
     r"$\log g^\delta$", r"$\log g^{TS}$",
     r"$\nu^\delta$",    r"$\nu^{TS}$"),
    ("δ vs LD", "log_g_delta", "log_g_ld", "nu_delta", "nu_ld",
     r"$\log g^\delta$", r"$\log g^{LD}$",
     r"$\nu^\delta$",    r"$\nu^{LD}$"),
    ("δ vs QU", "log_g_delta", "log_g_qu", "nu_delta", "nu_qu",
     r"$\log g^\delta$", r"$\log g^{QU}$",
     r"$\nu^\delta$",    r"$\nu^{QU}$"),
]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for col_idx, (lbl, rx, ry, nx, ny, rxl, ryl, nxl, nyl) in enumerate(scatter_pairs):
    r_r = nat[rx].corr(nat[ry])
    r_n = nat[nx].corr(nat[ny])
    # top row: raw
    ax = axes[0, col_idx]
    ax.scatter(nat[rx], nat[ry], s=5, alpha=0.35, color="#555")
    ax.set_xlabel(rxl, fontsize=9); ax.set_ylabel(ryl, fontsize=9)
    ax.set_title(f"{lbl} — raw  (r = {r_r:.3f})", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.grid(linewidth=0.4, alpha=0.4)
    # bottom row: residualized
    ax = axes[1, col_idx]
    ax.scatter(nat[nx], nat[ny], s=5, alpha=0.35, color="#1f77b4")
    ax.set_xlabel(nxl, fontsize=9); ax.set_ylabel(nyl, fontsize=9)
    ax.set_title(f"{lbl} — residualized  (r = {r_n:.3f})", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.grid(linewidth=0.4, alpha=0.4)
fig.suptitle(
    r"δ vs each s-type: raw (top) and residualized on $\Delta\log p_t$ (bottom)",
    fontsize=11,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_scatter.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot C saved: {p}")

print("\nDone.")
