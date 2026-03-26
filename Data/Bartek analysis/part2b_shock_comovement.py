"""
part2b_shock_comovement.py — Shock-productivity decomposition and comovement
=============================================================================
Diagnostic step between part2 (shock-rate construction) and part3 (Bartik
aggregation).

Purpose
-------
In the model, δ and s are exogenous AR(1) processes independent of technology
z.  In the data, both series contain endogenous components driven by aggregate
productivity and/or labor-market tightness.  This file asks two questions:

  1. How much do g^δ and g^s co-move across industries and over time?
     High co-movement implies the Bartik instruments B^δ and B^s will be
     nearly collinear, undermining separate identification.

  2. How much of that co-movement is explained by aggregate productivity
     and market tightness?  We residualize each series on the appropriate
     aggregate controls to isolate the idiosyncratic component.

Residualization strategy (shock-specific):
    log g^δ_{j,t}  = α_j + γ_δ  Δlog p_t                     + ν^δ_{j,t}
    log g^TS_{j,t} = α_j + γ_TS Δlog p_t                     + ν^TS_{j,t}
    log g^LD_{j,t} = α_j + γ_LD Δlog p_t + λ_LD log θ_t^nat  + ν^LD_{j,t}
    log g^QU_{j,t} = α_j + γ_QU Δlog p_t + λ_QU log θ_t^nat  + ν^QU_{j,t}

Motivation:
    δ    — endogenous product-line destruction driven by idiosyncratic
            productivity draws; Δlog p_t is the appropriate control.
    LD   — employer-initiated layoffs respond to both productivity and
            aggregate demand (tightness raises retention costs); both
            controls are included.
    QU   — quit propensity is driven by outside options summarised by
            market tightness θ = V/U; both controls are included.
    TS   — total separations = LD + QU; productivity control only
            (dominated by LD in recessions; tightness channel addressed
            by the separate LD/QU decomposition).

θ_t^nat = V_t^nat / U_t^nat: national JOLTS job openings ÷ civilian
unemployment (both seasonally adjusted, levels).  log θ_t is stationary
so it enters in log-levels; Δlog p_t is used because log p has a trend.

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
PROD_CACHE      = DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet"
TIGHTNESS_CACHE = DEFAULT_CACHE_DIR / "national_tightness_quarterly.parquet"
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

# ── 3b. national market tightness ──────────────────────────────────────────
# θ_t^nat = V_t^nat / U_t^nat.
# Primary source: FRED (JTSJOL ÷ UNEMPLOY, both SA thousands).
# Fallback: sum JOLTS state vacancies (already cached) over sum LAUS state
#   unemployment — avoids requiring FRED_API_KEY after the first run.
# log θ is stationary so we use log-levels as the regressor.

def _build_tightness_from_local() -> pd.DataFrame:
    """Compute national log(V/U) from already-cached state-level data."""
    # National vacancies: sum state JOLTS vacancies (units: thousands)
    _vac_path = DEFAULT_CACHE_DIR / "jolts_vacancies_state_quarterly.parquet"
    _laus_path = DEFAULT_CACHE_DIR / "laus_states_monthly.parquet"
    if not _vac_path.exists() or not _laus_path.exists():
        raise FileNotFoundError(
            "Cannot build tightness from local data: missing "
            f"{_vac_path} or {_laus_path}"
        )
    vac_df = pd.read_parquet(_vac_path)
    nat_vac = (vac_df.groupby("quarter_label")["vacancies"]
                     .sum().reset_index()
                     .rename(columns={"vacancies": "vac_nat"}))
    # National unemployment: mean-within-quarter of sum-state LAUS unemployment
    laus_df = pd.read_parquet(_laus_path)
    laus_df["date"] = pd.to_datetime(laus_df["date"])
    laus_df["quarter_label"] = laus_df["date"].map(_dt_to_ql)
    laus_df["unemp"] = laus_df["labor_force"] * laus_df["unemp_rate"] / 100
    nat_unemp = (laus_df.groupby("quarter_label")["unemp"]
                        .mean().reset_index()
                        .rename(columns={"unemp": "unemp_nat"}))
    t = nat_vac.merge(nat_unemp, on="quarter_label", how="inner")
    # vac_nat in thousands, unemp_nat in persons → convert vacancies to persons
    t["theta"]     = (t["vac_nat"] * 1000) / t["unemp_nat"]
    t["log_theta"] = np.log(t["theta"])
    return t[["quarter_label", "log_theta"]].copy()

if TIGHTNESS_CACHE.exists():
    tight = pd.read_parquet(TIGHTNESS_CACHE)
    print(f"\nTightness: loaded from cache ({TIGHTNESS_CACHE})")
else:
    api_key = os.environ.get("FRED_API_KEY", "")
    if api_key:
        print("\nFetching national tightness (JTSJOL, UNEMPLOY) from FRED ...")
        fred    = Fred(api_key=api_key)
        vac_m   = fred.get_series("JTSJOL")     # job openings, thousands, SA
        unemp_m = fred.get_series("UNEMPLOY")   # unemployment, thousands, SA
        tight_m = pd.DataFrame({"vac": vac_m, "unemp": unemp_m}).dropna()
        tight_m.index = pd.to_datetime(tight_m.index)
        tight_m["quarter_label"] = tight_m.index.map(_dt_to_ql)
        tight = (tight_m.groupby("quarter_label")[["vac", "unemp"]]
                        .mean().reset_index())
        tight["theta"]     = tight["vac"] / tight["unemp"]
        tight["log_theta"] = np.log(tight["theta"])
        tight = tight[["quarter_label", "log_theta"]].copy()
    else:
        print("\nFRED_API_KEY not set — building tightness from cached "
              "state vacancy + LAUS data ...")
        tight = _build_tightness_from_local()
    tight.to_parquet(TIGHTNESS_CACHE, index=False)
    print(f"  Saved: {TIGHTNESS_CACHE}")

tight = tight.sort_values("quarter_label").reset_index(drop=True)
print(f"  Tightness: {tight['quarter_label'].min()} – "
      f"{tight['quarter_label'].max()}  ({len(tight)} quarters)")

nat = nat.merge(tight, on="quarter_label", how="inner")
print(f"After merging tightness: {len(nat):,} industry-quarter obs")

# ── 4. residualization ─────────────────────────────────────────────────────
# Regression:  log g^k_{j,t} = α_j + γ_k Δlog p_t + ν^k_{j,t}
# Industry FEs absorbed via within-transformation (subtract industry means).
#
# Returns a DataFrame keyed by (industry_code, quarter_label) with the
# residual as a named column.  We merge back into nat rather than assigning
# into it directly, which is the most robust approach across pandas versions.

def _residualize(df, dep_col, resid_col, reg_cols=("dlog_p",)):
    """
    Within-industry OLS of dep_col on industry FEs + reg_cols.

    reg_cols : tuple of column names to use as regressors alongside
               industry fixed effects.  Homogeneous coefficients across
               industries (single γ / λ per regressor).

    df must be sorted by [industry_code, quarter_label] with a clean
    0-based index (call nat.sort_values(...).reset_index(drop=True) first).
    Returns a DataFrame with the same row order as df and a column named
    resid_col containing the residuals as a numpy-backed float64 series.
    """
    work = df[["industry_code", dep_col, *reg_cols]].copy().reset_index(drop=True)

    # Within-transform: subtract each industry's time-mean.
    work["y_dm"] = (work[dep_col]
                    - work.groupby("industry_code")[dep_col].transform("mean"))
    dm_cols = []
    for rc in reg_cols:
        dc = f"_dm_{rc}"
        work[dc] = (work[rc]
                    - work.groupby("industry_code")[rc].transform("mean"))
        dm_cols.append(dc)

    X_dm   = work[dm_cols].to_numpy()          # (n, K)
    y_dm   = work["y_dm"].to_numpy()           # (n,)
    gammas = np.linalg.lstsq(X_dm, y_dm, rcond=None)[0]   # (K,)
    resid  = y_dm - X_dm @ gammas

    ss_tot = (y_dm ** 2).sum()
    r2     = 1 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else np.nan

    coef_str = "  ".join(
        f"γ({rc}) = {g:+.4f}" for rc, g in zip(reg_cols, gammas)
    )
    print(f"  {resid_col:<12}  {coef_str}   R²(within) = {r2:.4f}")

    out = df[["industry_code", "quarter_label"]].copy().reset_index(drop=True)
    out[resid_col] = resid
    return out

# Sort nat by [industry_code, quarter_label] and reset to a clean 0-based
# index before residualizing.  _residualize works on a copy in the same
# row order, so returned arrays align positionally with nat's rows.
nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)

print(f"\n--- Residualization (within-industry OLS) ---")
print(f"  δ,  TS : regressors = [Δlog p_t]")
print(f"  LD, QU : regressors = [Δlog p_t,  log θ_t^nat]")
resid_d  = _residualize(nat, "log_g_delta", "nu_delta",
                        reg_cols=("dlog_p",))
resid_s  = _residualize(nat, "log_g_s",     "nu_s",
                        reg_cols=("dlog_p",))
resid_ld = _residualize(nat, "log_g_ld",    "nu_ld",
                        reg_cols=("dlog_p", "log_theta"))
resid_qu = _residualize(nat, "log_g_qu",    "nu_qu",
                        reg_cols=("dlog_p", "log_theta"))

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
    """Cross-industry Pearson r matrix (industries × industries).

    Builds the quarter × industry matrix via direct numpy integer indexing,
    completely bypassing pandas groupby / pivot_table / pivot — all of which
    have version-specific bugs in pandas 2.x that silently drop the aggregated
    column from the result (causing KeyError on the subsequent values lookup).

    nat is already aggregated at (industry_code, quarter_label) level, so each
    (quarter_label, industry_label) cell holds at most one value; if a
    duplicate somehow exists the last row wins (harmless).
    """
    sub = (df[["quarter_label", "industry_label", col]]
           .dropna(subset=[col])
           .reset_index(drop=True))
    quarters   = sorted(sub["quarter_label"].unique())
    industries = sorted(sub["industry_label"].unique())
    if not quarters or not industries:
        return pd.DataFrame()
    q_map = {q: i for i, q in enumerate(quarters)}
    i_map = {ind: i for i, ind in enumerate(industries)}
    mat = np.full((len(quarters), len(industries)), np.nan)
    qi = sub["quarter_label"].map(q_map).to_numpy(dtype=int)
    ii = sub["industry_label"].map(i_map).to_numpy(dtype=int)
    mat[qi, ii] = sub[col].to_numpy()
    wide = pd.DataFrame(mat, index=quarters, columns=industries)
    return wide.dropna(how="all", axis=1).corr()

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
    r"δ vs each s-type: raw (top) vs residualized (bottom)"
    "\n"
    r"δ/TS: $\Delta\log p_t$ only — LD/QU: $\Delta\log p_t + \log\theta_t^{nat}$",
    fontsize=11,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_scatter.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot C saved: {p}")

print("\nDone.")
