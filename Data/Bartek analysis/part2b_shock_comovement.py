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
    g^k_{j,t}: industry j shock rate of type k ∈ {δ, s} in quarter t
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
from construct_s_instrument import SHOCK_RATES_S_PATH

# ── paths ──────────────────────────────────────────────────────────────────
PROD_CACHE   = DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet"
RESID_D_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_S_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_s_resid.parquet"
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

for p in [SHOCK_RATES_PATH, SHOCK_RATES_S_PATH]:
    if not p.exists():
        script = ("part2_shock_rates.py" if "shock_rates" in p.name
                  else "part2_shock_rates_s.py")
        sys.exit(f"ERROR: {p} not found.  Run {script} first.")

# ── 2. national industry series ────────────────────────────────────────────
# Average across states to recover the national series.
# The LOO correction only affects the denominator and is negligible here.

raw_d = pd.read_parquet(SHOCK_RATES_PATH)
raw_s = pd.read_parquet(SHOCK_RATES_S_PATH)

nat_d = (raw_d.groupby(["industry_code", "quarter_label"])["g_delta_loo"]
               .mean().reset_index()
               .rename(columns={"g_delta_loo": "g_delta"}))

nat_s = (raw_s.groupby(["industry_code", "quarter_label"])["g_s_loo"]
               .mean().reset_index()
               .rename(columns={"g_s_loo": "g_s"}))

print(f"\nδ: {nat_d['quarter_label'].nunique()} quarters, "
      f"{nat_d['industry_code'].nunique()} industries")
print(f"s: {nat_s['quarter_label'].nunique()} quarters, "
      f"{nat_s['industry_code'].nunique()} industries")

# Restrict to overlapping sample (s starts 2001Q1).
nat = nat_d.merge(nat_s, on=["industry_code", "quarter_label"], how="inner")
print(f"Overlapping sample: {nat['quarter_label'].min()} – "
      f"{nat['quarter_label'].max()}  "
      f"({nat['quarter_label'].nunique()} quarters, "
      f"{nat['industry_code'].nunique()} industries, "
      f"{len(nat):,} obs)")

# Remove zeros/negatives (log undefined); attach labels.
nat = nat[(nat["g_delta"] > 0) & (nat["g_s"] > 0)].copy()
nat["log_g_delta"]    = np.log(nat["g_delta"])
nat["log_g_s"]        = np.log(nat["g_s"])
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
    api_key = os.environ.get("FRED_API_KEY", "d35aabd7dc07cd94481af3d1e2f0ecf3")
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
resid_d = _residualize(nat, "log_g_delta", "nu_delta")
resid_s = _residualize(nat, "log_g_s",     "nu_s")

# Merge residuals back into nat on keys — no index alignment issues.
# Assign residuals by position — nat is already sorted and 0-based indexed,
# and _residualize works on an identically-sorted copy, so the numpy arrays
# are in the same row order.  No merge or index alignment needed.
nat["nu_delta"] = resid_d["nu_delta"].to_numpy()
nat["nu_s"]     = resid_s["nu_s"].to_numpy()

# Confirm columns are present and non-null.
for col in ["nu_delta", "nu_s"]:
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

r_raw   = _pooled_r(nat, "log_g_delta", "log_g_s")
r_resid = _pooled_r(nat, "nu_delta",    "nu_s")

print(f"\n--- Pooled cross-series correlation corr(g^δ, g^s) ---")
print(f"  Raw        r = {r_raw:+.4f}")
print(f"  Residualized r = {r_resid:+.4f}  "
      f"(reduction = {r_raw - r_resid:+.4f})")
if abs(r_raw - r_resid) > 0.15:
    print("  → Productivity explains most co-movement.")
else:
    print("  → Productivity explains little; deeper co-movement present.")

print(f"\n--- Per-industry corr(g^δ, g^s) ---")
by_raw   = _by_industry_r(nat, "log_g_delta", "log_g_s").rename(columns={"r": "r_raw"})
by_resid = _by_industry_r(nat, "nu_delta",    "nu_s"   ).rename(columns={"r": "r_resid"})
by_ind   = by_raw.merge(by_resid, on="industry_label").sort_values("r_raw")
print(by_ind.to_string(index=False))

# ── 6. save residualized series ────────────────────────────────────────────

nat[["industry_code", "quarter_label", "nu_delta"]].to_parquet(RESID_D_PATH, index=False)
nat[["industry_code", "quarter_label", "nu_s"]    ].to_parquet(RESID_S_PATH, index=False)
print(f"\nSaved: {RESID_D_PATH}")
print(f"Saved: {RESID_S_PATH}")

# ── 7. plots ───────────────────────────────────────────────────────────────

industries = sorted(nat["industry_label"].unique())
quarters   = sorted(nat["quarter_label"].unique())
dates      = [_ql_to_dt(q) for q in quarters]
n_ind      = len(industries)

# ── plot A: raw log series by industry ────────────────────────────────────
fig, axes = plt.subplots(n_ind, 2,
                         figsize=(14, 2.2 * n_ind),
                         sharex=True, squeeze=False)
for i, ind in enumerate(industries):
    sub = (nat[nat["industry_label"] == ind]
           .set_index("quarter_label")
           .reindex(quarters))
    for ax, col, color in [(axes[i, 0], "log_g_delta", "#1f77b4"),
                            (axes[i, 1], "log_g_s",    "#d62728")]:
        ax.plot(dates, sub[col].values, color=color, linewidth=1.2)
        _add_recessions(ax)
        ax.set_ylabel(ind, fontsize=7)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
axes[0, 0].set_title(r"$\log g^\delta_{j,t}$  (δ shock rate)",  fontsize=10)
axes[0, 1].set_title(r"$\log g^s_{j,t}$  (s shock rate)", fontsize=10)
fig.suptitle("National industry shock rates by supersector (raw)",
             fontsize=11, y=1.005)
fig.tight_layout()
p = RESULTS_DIR / "comovement_raw_series.png"
fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot A saved: {p}")

# ── plot B: 2×2 correlation heatmaps ──────────────────────────────────────
def _corr_matrix(df, col):
    """Cross-industry Pearson r matrix (industries × industries)."""
    wide = df.pivot_table(index="quarter_label",
                          columns="industry_label",
                          values=col,
                          aggfunc="mean")   # aggfunc guards against duplicates
    return wide.corr()

cm = {
    "raw_d":   _corr_matrix(nat, "log_g_delta"),
    "raw_s":   _corr_matrix(nat, "log_g_s"),
    "resid_d": _corr_matrix(nat, "nu_delta"),
    "resid_s": _corr_matrix(nat, "nu_s"),
}

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
panels = [
    (axes[0, 0], cm["raw_d"],   r"$\log g^\delta_{j,t}$ (raw)"),
    (axes[0, 1], cm["raw_s"],   r"$\log g^s_{j,t}$ (raw)"),
    (axes[1, 0], cm["resid_d"], r"$\nu^\delta_{j,t}$ (residualized)"),
    (axes[1, 1], cm["resid_s"], r"$\nu^s_{j,t}$ (residualized)"),
]
for ax, mat, title in panels:
    im = ax.imshow(mat.values, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    labels = mat.columns.tolist()
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45,
                                                           ha="right", fontsize=7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(title, fontsize=10)
    for ii in range(len(labels)):
        for jj in range(len(labels)):
            v = mat.values[ii, jj]
            ax.text(jj, ii, f"{v:.2f}", ha="center", va="center",
                    fontsize=6, color="white" if abs(v) > 0.7 else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle(
    "Cross-industry correlation matrices: raw vs. productivity-residualized\n"
    f"Pooled corr(g^δ, g^s): raw = {r_raw:.3f}  →  residualized = {r_resid:.3f}",
    fontsize=11, y=1.01,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_correlation_heatmaps.png"
fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot B saved: {p}")

# ── plot C: cross-series scatter raw vs. residualized ─────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
axes[0].scatter(nat["log_g_delta"], nat["log_g_s"], s=6, alpha=0.4, color="#555")
axes[0].set_xlabel(r"$\log g^\delta_{j,t}$")
axes[0].set_ylabel(r"$\log g^s_{j,t}$")
axes[0].set_title(f"Raw  (r = {r_raw:.3f})", fontsize=10)

axes[1].scatter(nat["nu_delta"], nat["nu_s"], s=6, alpha=0.4, color="#1f77b4")
axes[1].set_xlabel(r"$\nu^\delta_{j,t}$")
axes[1].set_ylabel(r"$\nu^s_{j,t}$")
axes[1].set_title(f"Residualized  (r = {r_resid:.3f})", fontsize=10)

for ax in axes:
    ax.axhline(0, color="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.grid(linewidth=0.4, alpha=0.4)

fig.suptitle(
    r"$g^\delta_{j,t}$ vs. $g^s_{j,t}$: raw and residualized on $\Delta\log p_t$",
    fontsize=11,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_scatter.png"
fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot C saved: {p}")

print("\nDone.")