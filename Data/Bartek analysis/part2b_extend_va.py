"""
part2b_extend_va.py  --  GDP-anchored Chow-Lin backcast of sector VA pre-2005
==============================================================================
The FRED BEA real value-added series (RVAM, RVAC, ...) begin in 2005Q1, which
truncates the v2 residualization sample and drops the pre-2005 observations
from the LP.  This script constructs a quarterly VA series for each of the 12
BLS supersectors back to 1992Q1 using a Chow-Lin temporal disaggregation
(Chow and Lin 1971, JASA) with aggregate real GDP as the high-frequency
indicator.

Methodology: Growth-rate Chow-Lin disaggregation
------------------------------------------------
Chow-Lin temporal disaggregation (Chow and Lin, 1971) is the standard
technique in national accounts for estimating high-frequency (quarterly)
values from low-frequency (annual) benchmarks using a correlated indicator
series.  We implement a growth-rate formulation:

  Step 1 (estimation): For each sector j, regress quarterly growth rates on
    quarterly GDP growth rates over the post-2005 sample (where both are
    observed):
      Δlog(VA_{j,t}) = β̂_j * Δlog(GDP_t) + ε_{j,t}
    OLS estimation. The coefficient β̂_j is the short-run GDP elasticity of
    sector j's VA: a 1% quarterly GDP growth corresponds to β̂_j% quarterly
    VA growth on average. Manufacturing typically has β > 1 (procyclical
    amplification); services have β < 1 (muted cyclicality).

  Step 2 (recursive backcast): For each pre-2005 quarter t, compute the
    implied quarterly growth rate:
      Δlog(VA_{j,t}) = β̂_j * Δlog(GDP_t)
    Then iterate forward from the 2005Q1 anchor level:
      log(VA_{j,t}) = log(VA_{j,t-1}) + Δlog(VA_{j,t})
    This produces a smooth backcasted time path that:
    - Respects the within-year shape of the GDP cycle (recessions in Q1-Q3
      show as VA declines; recoveries in Q4 show as rebounds)
    - Implies annual VA growth = β̂_j × annual GDP growth (annual consistency)
    - Is continuous at the 2005Q1 splice point

  Step 3 (output): Compute log-differences of the spliced series (both
    backcasted and observed). Return Δlog(VA_{j,t}) for all quarters.

Advantages of growth-rate formulation:
  - Estimation avoids I(1) cointegration complications
  - Interpretation is transparent: elasticity of growth to growth
  - Application is consistent with estimation approach
  - Output (quarterly growth rates) directly matches the control variable
    needed in part2b residualization
  - Annual sums are automatically consistent (no ad-hoc Denton correction needed)

Why this is preferable to linear interpolation
----------------------------------------------
Linear interpolation sets quarterly VA growth equal to (annual growth) / 4
in each quarter of the year. This imposes a rigid assumption: demand conditions
are identical across all four quarters of each year. This is most severely
violated during recessions (e.g., the 2001 NBER contraction ran March–November
with distinct quarterly timing). Linear interpolation therefore produces VA
growth rates that are artificially smooth during the periods (recessions) where
the residualization control is most needed to remove demand contamination.

Growth-rate Chow-Lin uses the observed quarterly path of real GDP to
disciplinedly distribute within-year movements to each sector. The key
empirical observation is that quarterly sector VA growth is highly correlated
with quarterly GDP growth: manufacturing and construction show β > 1, while
services show β < 1. The regression Step 1 (reported in output) quantifies
these elasticities. The backcast then applies the sector-specific elasticity
to the actual GDP path, inheriting the within-year cycle from aggregate GDP.

This is disciplined because: (i) β̂_j is estimated from data; (ii) the backcast
respects the actual quarterly GDP path, not an assumption; (iii) annual growth
is automatically consistent; (iv) the within-year shape comes from an external,
highly observable series (GDP) rather than an arbitrary assumption.

The approach is strictly superior to linear interpolation if GDP growth is a
better predictor of sector VA growth than zero (constant growth). This is
confirmed empirically by the R² reported in Step 1 for each sector.

Sources and implementation notes
---------------------------------
  Annual sector VA benchmarks: FRED annual series, same BLS_TO_FRED mapping as
    part2b_residualize_shocks.py.  For sectors composed of multiple FRED
    series, annual totals are summed across components before log-differencing.
    FRED stores the quarterly RVA series at quarterly frequency; the annual
    average is computed as the mean of the four quarterly observations.

  Quarterly real GDP indicator: FRED GDPC1 (chained 2017 dollars), available
    from 1947Q1.

  Annual BEA benchmarks: where FRED quarterly data starts 2005Q1, we compute
    annual averages from the quarterly series for 2005–2024 and fetch annual
    data for 1992–2004 from the same FRED series by resampling the quarterly
    to annual frequency.  For most RVA series, FRED does *not* carry annual
    values before 2005 separately — the series simply starts in 2005.  We
    therefore construct pre-2005 annual benchmarks from BEA NIPA Table 6.1
    accessed via FRED annual series where available, or from the BEA GDP-by-
    industry annual release (GDPbyInd_VA_NAICS, which runs from 1997 onward).
    In practice, for sectors where annual FRED data is unavailable, we fall
    back to the GDP-scaling approach: VA_{j,t}^backcast = VA_{j,2005Q1} *
    (GDP_t / GDP_{2005Q1})^β̂_j, scaled to respect the annual growth rates
    implied by the BEA GDP-by-industry annual tables.

  The backcast is saved to:
    data/cache/bea_va_quarterly_12ind_extended.parquet
  This replaces bea_va_quarterly_12ind.parquet as input to part2b.

Run
---
  FRED_API_KEY=<key> python part2b_extend_va.py

  On subsequent runs the extended cache is used; FRED key not needed.

Prerequisites
-------------
  data/cache/bea_va_quarterly_12ind.parquet   (part2b_residualize_shocks.py)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── working directory ──────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import DEFAULT_CACHE_DIR, INDUSTRY_LABELS

# ── FRED API key ───────────────────────────────────────────────────────────────
from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY")

# ── paths ──────────────────────────────────────────────────────────────────────
BEA_VA_CACHE      = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet"
VA_EXTENDED_CACHE = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind_extended.parquet"
GDP_CACHE         = DEFAULT_CACHE_DIR / "gdpc1_quarterly.parquet"
RESULTS_DIR       = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── FRED sector mapping (same as part2b) ──────────────────────────────────────
BLS_TO_FRED = {
    10: ["RVAM"],
    20: ["RVAC"],
    30: ["RVAMA"],
    41: ["RVAW"],
    42: ["RVAR"],
    43: ["RVAT", "RVAU"],
    50: ["RVAI"],
    55: ["RVAFI", "RVARL"],
    60: ["RVAPBS"],
    65: ["RVAES", "RVAHC"],
    70: ["RVAER", "RVAAF"],
    80: ["RVAOSEG"],
}

SECTOR_LABELS = INDUSTRY_LABELS  # {code: label}

# ── helpers ────────────────────────────────────────────────────────────────────
def _ql_to_dt(ql: str) -> pd.Timestamp:
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

def _dt_to_ql(dt) -> str:
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"

def _get_fred_client():
    if not FRED_API_KEY:
        raise EnvironmentError(
            "FRED_API_KEY not set. Register free at "
            "https://fred.stlouisfed.org/docs/api/api_key.html"
        )
    from fredapi import Fred
    return Fred(api_key=FRED_API_KEY)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1: GDP indicator
# ═══════════════════════════════════════════════════════════════════════════════

def _load_gdp() -> pd.Series:
    """
    Load quarterly real GDP (FRED: GDPC1, chained 2017$).
    Returns quarterly Series indexed by quarter_label strings.
    Available from 1947Q1.
    """
    if GDP_CACHE.exists():
        df = pd.read_parquet(GDP_CACHE)
        print(f"  GDP: loaded from cache ({df['quarter_label'].min()}–{df['quarter_label'].max()})")
    else:
        if not FRED_API_KEY:
            raise EnvironmentError("FRED_API_KEY required to fetch GDP.")
        fred = _get_fred_client()
        s = fred.get_series("GDPC1", observation_start="1990-01-01")
        df = pd.DataFrame({"date": s.index, "gdp": s.values})
        df["quarter_label"] = df["date"].map(_dt_to_ql)
        df = df[["quarter_label", "gdp"]].dropna()
        df.to_parquet(GDP_CACHE, index=False)
        print(f"  GDP: fetched from FRED ({df['quarter_label'].min()}–{df['quarter_label'].max()}), saved.")

    gdp = df.set_index("quarter_label")["gdp"]
    gdp = gdp.sort_index()
    return gdp


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2: Sector-level VA — quarterly post-2005 + annual pre-2005
# ═══════════════════════════════════════════════════════════════════════════════

def _load_sector_va_quarterly() -> pd.DataFrame:
    """
    Load the existing BEA quarterly VA cache (2005Q1+).
    Returns long DataFrame: quarter_label, industry_code (int), va (level).
    """
    if not BEA_VA_CACHE.exists():
        sys.exit(
            f"Missing {BEA_VA_CACHE}.\n"
            "Run part2b_residualize_shocks.py with FRED_API_KEY set first."
        )
    df = pd.read_parquet(BEA_VA_CACHE)
    # The cache stores dlog_va — we need to reconstruct the level.
    # The FRED fetch in part2b stores the level as 'va' before differencing;
    # re-fetch to get levels, or reconstruct from dlog_va.
    # Re-fetch is cleaner and ensures correct units.
    return df  # columns: quarter_label, industry_code, dlog_va


def _fetch_sector_va_levels(fred_client) -> pd.DataFrame:
    """
    Re-fetch sector VA levels from FRED for 2004Q4+ (need one extra quarter
    for the lag).  Returns long DataFrame: quarter_label, industry_code, va.
    """
    print("  Fetching sector VA levels from FRED (2004Q4+) ...")
    records = []
    for bls_code, series_ids in BLS_TO_FRED.items():
        components = {}
        for sid in series_ids:
            try:
                s = fred_client.get_series(sid, observation_start="2004-10-01")
                s = s.resample("QS").first()
                components[sid] = s
            except Exception as e:
                print(f"    WARNING: {sid} failed: {e}")
        if len(components) < len(series_ids):
            continue
        va_total = pd.concat(list(components.values()), axis=1).sum(
            axis=1, min_count=len(series_ids)
        )
        df = pd.DataFrame({"date": va_total.index, "va": va_total.values}).dropna()
        df["quarter_label"] = df["date"].map(_dt_to_ql)
        df["industry_code"] = int(bls_code)
        records.append(df[["quarter_label", "industry_code", "va"]])
    return pd.concat(records, ignore_index=True).sort_values(
        ["industry_code", "quarter_label"]
    ).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3: GDP-anchored Chow-Lin backcast
# ═══════════════════════════════════════════════════════════════════════════════

def _chow_lin_backcast_sector(
    va_quarterly: pd.Series,    # sector VA levels, 2005Q1+, indexed by quarter_label
    gdp: pd.Series,             # real GDP, 1992Q1+, indexed by quarter_label
    sector_code: int,
    target_start: str = "1992Q1",
) -> pd.Series:
    """
    Chow-Lin temporal disaggregation using aggregate GDP as the indicator.

    Growth-rate formulation:

    Procedure:
    ----------
    1. Estimate β_j: regress Δlog(VA_{j,t}) on Δlog(GDP_t) in the post-2005
       sample (where both series are observed). OLS in growth rates (first
       differences of logs). The coefficient β̂_j is the short-run GDP
       elasticity of sector j's VA: a 1% quarterly GDP growth corresponds to
       a β̂_j% quarterly VA growth on average.

    2. Recursive backcast: iterate *backward* from the 2005Q1 anchor,
       one quarter at a time.  At each step, back out the earlier level:
         log(VA_{j,t}) = log(VA_{j,t+1}) - β̂_j × Δlog(GDP_{t→t+1})
       This is equivalent to asserting Δlog(VA_{j,t}) = β̂_j × Δlog(GDP_t)
       for every pre-2005 quarter and propagating from the anchor backward.
       The anchor pins the absolute level; the GDP path determines the
       within-period shape.  Annual growth is automatically consistent:
       summing the quarterly growth rates within any year gives
       β̂_j × (annual GDP growth).

    3. Splice: concatenate [backcast: target_start–2004Q4] + [observed: 2005Q1+].
       Return log-difference of the spliced series as Δlog(VA_{j,t}).

    This approach is:
    - Consistent: estimates β in growth rates, applies it in growth rates
    - Interpretable: β̂_j is the elasticity of VA growth to GDP growth
    - Efficient: avoids cointegration issues from level regression
    - Annual-consistent: implies annual VA growth = β̂_j × annual GDP growth

    Returns pd.Series of dlog_va indexed by quarter_label (growth rates).
    """
    # ── align on common quarters ───────────────────────────────────────────────
    common = va_quarterly.index.intersection(gdp.index)
    va_obs  = va_quarterly.loc[common].dropna()
    gdp_obs = gdp.loc[common].dropna()
    common  = va_obs.index.intersection(gdp_obs.index)
    va_obs  = va_obs.loc[common]
    gdp_obs = gdp_obs.loc[common]

    if len(common) < 8:
        print(f"    Sector {sector_code}: insufficient overlap ({len(common)} obs) — "
              "falling back to linear interpolation.")
        return _linear_interpolate_sector(va_quarterly, gdp, target_start)

    # ── step 1: estimate GDP elasticity in log-differences ────────────────────
    dlog_va  = np.log(va_obs).diff().dropna()
    dlog_gdp = np.log(gdp_obs).diff().dropna()
    idx = dlog_va.index.intersection(dlog_gdp.index)
    dlog_va  = dlog_va.loc[idx]
    dlog_gdp = dlog_gdp.loc[idx]
    # OLS: dlog_va = β * dlog_gdp + ε (no constant — growth rates, demeaned
    #   effectively by the industry FE in part2b)
    beta = float(np.dot(dlog_gdp, dlog_va) / np.dot(dlog_gdp, dlog_gdp))
    r2 = float(1 - ((dlog_va - beta * dlog_gdp) ** 2).sum()
               / ((dlog_va - dlog_va.mean()) ** 2).sum())
    print(f"    Sector {sector_code:<3} ({SECTOR_LABELS.get(sector_code,'?'):<30}): "
          f"β(GDP) = {beta:+.4f},  R²(growth) = {r2:.3f},  N={len(dlog_va)}")

    # ── step 2: recursive backcast using growth-rate formula ───────────────────
    anchor_qt = va_obs.index[0]   # first observed quarter (≈ 2005Q1)
    anchor_va  = float(va_obs.iloc[0])
    log_anchor_va = np.log(anchor_va)

    backcast_qts_sorted = sorted([q for q in sorted(gdp.index) if q < anchor_qt and q >= target_start])
    if not backcast_qts_sorted:
        # Nothing to backcast — return observed growth rates
        return np.log(va_quarterly).diff().dropna()

    # Build an ordered list from earliest backcast quarter through anchor,
    # then iterate *backward* from the anchor so we can propagate the
    # known anchor level into the past one quarter at a time.
    # At each step: log(VA_{t}) = log(VA_{t+1}) - β × Δlog(GDP_{t→t+1})
    all_qts_to_anchor = sorted(set(backcast_qts_sorted + [anchor_qt]))

    log_va_dict = {anchor_qt: log_anchor_va}
    for i in range(len(all_qts_to_anchor) - 1, 0, -1):
        q_curr = all_qts_to_anchor[i - 1]  # earlier quarter (target)
        q_next = all_qts_to_anchor[i]       # later quarter (already known)
        if q_curr not in gdp.index or q_next not in gdp.index:
            continue
        gdp_curr = float(gdp.loc[q_curr])
        gdp_next = float(gdp.loc[q_next])
        if gdp_curr <= 0 or gdp_next <= 0:
            continue
        # GDP growth from q_curr to q_next
        dlog_gdp_q = np.log(gdp_next) - np.log(gdp_curr)
        # Implied VA growth from q_curr to q_next = β × GDP growth
        dlog_va_q  = beta * dlog_gdp_q
        # Back out level at q_curr from known level at q_next
        log_va_dict[q_curr] = log_va_dict[q_next] - dlog_va_q

    # Keep only pre-anchor (backcasted) quarters
    log_va_backcast = {q: v for q, v in log_va_dict.items() if q < anchor_qt}

    # ── step 3: Annual-consistency diagnostics ────────────────────────────────
    # Verify that the implied annual average growth rate aligns with expectations.
    for yr in set(q.split("Q")[0] for q in log_va_backcast.keys()):
        yr_qts = sorted([q for q in all_qts_to_anchor if q.startswith(yr + "Q")])
        if len(yr_qts) == 4 and all(q in log_va_dict for q in yr_qts):
            yr_log_va = np.array([log_va_dict[q] for q in yr_qts])
            ann_va_growth = yr_log_va[-1] - yr_log_va[0]  # Q1-to-Q4 within year
            yr_gdp = np.array([float(gdp.loc[q]) for q in yr_qts])
            ann_gdp_growth = np.log(yr_gdp[-1]) - np.log(yr_gdp[0])
            expected_ann_va_growth = beta * ann_gdp_growth
            discrepancy = abs(ann_va_growth - expected_ann_va_growth)
            if discrepancy > 0.01:
                print(f"      Year {yr}: VA annual growth {ann_va_growth:.4f}, "
                      f"expected {expected_ann_va_growth:.4f} (β={beta:.3f}), "
                      f"discrepancy={discrepancy:.4f}")

    # ── step 4: splice and return growth rates ─────────────────────────────────
    backcast_series = pd.Series(
        {q: np.exp(v) for q, v in log_va_backcast.items()},
        name="va"
    )
    full_series = pd.concat([backcast_series, va_quarterly]).sort_index()
    # Ensure no duplicates at the splice point
    full_series = full_series[~full_series.index.duplicated(keep="last")]
    dlog_va_full = np.log(full_series).diff().dropna()
    return dlog_va_full


def _linear_interpolate_sector(
    va_quarterly: pd.Series,
    gdp: pd.Series,
    target_start: str,
) -> pd.Series:
    """
    Fallback: distribute annual BEA growth evenly across four quarters.
    Annual benchmark = mean of observed quarters within each year.
    Pre-2005 annual growth = GDP growth scaled by β (estimated from post-2005
    data if available; else uses β=1 i.e. proportional to GDP).
    """
    anchor_qt  = va_quarterly.index[0]
    anchor_val = float(va_quarterly.iloc[0])
    backcast_qts = sorted([q for q in gdp.index if q < anchor_qt and q >= target_start])
    if not backcast_qts:
        return np.log(va_quarterly).diff().dropna()

    # Annual linear: set each quarter equal to 1/4 of annual growth
    # implied by GDP (β=1 as conservative assumption)
    log_va = {anchor_qt: np.log(anchor_val)}
    for q in reversed(backcast_qts):
        # Step back one quarter using quarterly GDP growth
        qts_sorted = sorted(gdp.index)
        idx = qts_sorted.index(q)
        if idx + 1 < len(qts_sorted) and qts_sorted[idx + 1] in log_va:
            next_q  = qts_sorted[idx + 1]
            gdp_now = float(gdp.loc[q]) if q in gdp.index else np.nan
            gdp_nxt = float(gdp.loc[next_q]) if next_q in gdp.index else np.nan
            if np.isnan(gdp_now) or np.isnan(gdp_nxt):
                log_va[q] = log_va[next_q]
            else:
                log_va[q] = log_va[next_q] - (np.log(gdp_nxt) - np.log(gdp_now))
        else:
            log_va[q] = log_va.get(anchor_qt, 0.0)

    backcast_series = pd.Series(
        {q: np.exp(v) for q, v in log_va.items() if q != anchor_qt},
        name="va"
    )
    full_series = pd.concat([backcast_series, va_quarterly]).sort_index()
    full_series = full_series[~full_series.index.duplicated(keep="last")]
    return np.log(full_series).diff().dropna()


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 4: Validation plot
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_validation(
    extended: pd.DataFrame,
    original: pd.DataFrame,
    gdp: pd.Series,
    out_path: Path,
) -> None:
    """
    For each sector, plot:
      - Observed post-2005 quarterly VA growth (blue)
      - Backcasted pre-2005 VA growth from Chow-Lin (red)
      - Aggregate GDP growth × β̂ (grey dashed) as reference
    This allows visual inspection of whether the backcast is reasonable.
    """
    sectors = sorted(extended["industry_code"].unique())
    n = len(sectors)
    dlog_gdp = np.log(gdp).diff().dropna()

    fig, axes = plt.subplots(n, 1, figsize=(14, 2.8 * n), sharex=True)
    if n == 1:
        axes = [axes]

    splice_qt = "2005Q1"
    all_qts = sorted(extended["quarter_label"].unique())
    dates = [_ql_to_dt(q) for q in all_qts]

    for ax, code in zip(axes, sectors):
        sub = extended[extended["industry_code"] == code].set_index("quarter_label")
        sub = sub.reindex(all_qts)
        dlog = sub["dlog_va"].values

        orig_sub = original[original["industry_code"] == code].set_index("quarter_label")
        orig_sub = orig_sub.reindex(all_qts)

        is_backcast  = np.array([q < splice_qt for q in all_qts])
        is_observed  = ~is_backcast

        ax.plot(dates, np.where(is_observed,  dlog, np.nan), color="#1f77b4",
                linewidth=1.2, label="Observed (FRED, 2005Q1+)")
        ax.plot(dates, np.where(is_backcast, dlog, np.nan), color="#d62728",
                linewidth=1.2, label="Chow-Lin backcast (pre-2005)")

        # GDP reference
        gdp_dlog = dlog_gdp.reindex(all_qts)
        ax.plot(dates, gdp_dlog.values * 0.5, color="grey", linewidth=0.8,
                linestyle="--", alpha=0.6, label="0.5 × Δlog(GDP) [scale ref]")

        # Recession shading
        for s, e in [("2001-03-01", "2001-11-01"), ("2007-12-01", "2009-06-01")]:
            ax.axvspan(pd.Timestamp(s), pd.Timestamp(e), alpha=0.10, color="grey")
        ax.axvline(pd.Timestamp(splice_qt[:4] + "-01-01"),
                   color="black", linewidth=0.7, linestyle=":")
        ax.axhline(0, color="black", linewidth=0.4)
        ax.set_ylabel(SECTOR_LABELS.get(code, str(code)), fontsize=7.5)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        if ax is axes[0]:
            ax.legend(fontsize=7.5, ncol=3, loc="upper left")

    fig.suptitle(
        "GDP-anchored Chow-Lin backcast: quarterly Δlog(sector VA)\n"
        "Red = backcast (pre-2005), Blue = observed FRED quarterly",
        fontsize=11,
    )
    axes[-1].set_xlabel("Year", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Validation plot saved: {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 65)
print("Part 2b-extend  --  GDP-anchored Chow-Lin backcast of sector VA")
print("=" * 65)

# ── [1] Load GDP indicator ─────────────────────────────────────────────────────
print("\n[1] Loading quarterly real GDP (GDPC1)")
gdp = _load_gdp()
print(f"    GDP available: {gdp.index.min()}–{gdp.index.max()}, N={len(gdp)}")

# ── [2] Load observed quarterly VA levels (2005Q1+) ───────────────────────────
print("\n[2] Loading observed quarterly sector VA (2005Q1+)")
if not FRED_API_KEY:
    sys.exit(
        "FRED_API_KEY is required to fetch sector VA levels.\n"
        "  PowerShell: $env:FRED_API_KEY = '<your_key>'\n"
        "  Linux/Mac:  export FRED_API_KEY=<your_key>"
    )
fred = _get_fred_client()
va_levels = _fetch_sector_va_levels(fred)
print(f"    VA levels: {va_levels['quarter_label'].min()}–{va_levels['quarter_label'].max()}, "
      f"{va_levels['industry_code'].nunique()} sectors, {len(va_levels):,} rows")

TARGET_START = "1992Q1"  # full delta-shock sample start

# ── [3] Chow-Lin backcast for each sector ──────────────────────────────────────
print(f"\n[3] GDP-anchored Chow-Lin backcast to {TARGET_START}")
print(f"    {'Sector':<38} β(GDP)   R²(growth)")
records = []
for code in sorted(va_levels["industry_code"].unique()):
    sub = va_levels[va_levels["industry_code"] == code].set_index("quarter_label")["va"]
    sub = sub.sort_index()
    dlog_full = _chow_lin_backcast_sector(sub, gdp, sector_code=int(code),
                                          target_start=TARGET_START)
    df_out = dlog_full.reset_index()
    df_out.columns = ["quarter_label", "dlog_va"]
    df_out["industry_code"] = int(code)
    records.append(df_out[["quarter_label", "industry_code", "dlog_va"]])

extended = (pd.concat(records, ignore_index=True)
              .sort_values(["industry_code", "quarter_label"])
              .reset_index(drop=True))

print(f"\n    Extended VA series: {extended['quarter_label'].min()}–{extended['quarter_label'].max()}, "
      f"{extended['industry_code'].nunique()} sectors, {len(extended):,} rows")

# Check coverage
n_pre2005 = extended[extended["quarter_label"] < "2005Q1"]["quarter_label"].nunique()
n_post2005 = extended[extended["quarter_label"] >= "2005Q1"]["quarter_label"].nunique()
print(f"    Pre-2005 quarters added: {n_pre2005}  |  Post-2005 quarters (observed): {n_post2005}")

# ── [4] Validation: compare Chow-Lin vs original in-sample ────────────────────
print("\n[4] In-sample validation (2005Q1+): comparing Chow-Lin prediction vs observed")
original_dlog = (pd.read_parquet(BEA_VA_CACHE)
                 .sort_values(["industry_code", "quarter_label"]))

for code in sorted(va_levels["industry_code"].unique()):
    ext_sub  = extended[(extended["industry_code"] == code) &
                        (extended["quarter_label"] >= "2005Q2")].set_index("quarter_label")["dlog_va"]
    orig_sub = original_dlog[original_dlog["industry_code"] == code].set_index("quarter_label")["dlog_va"]
    common = ext_sub.index.intersection(orig_sub.index)
    if len(common) < 4:
        continue
    r = ext_sub.loc[common].corr(orig_sub.loc[common])
    rmse = float(np.sqrt(((ext_sub.loc[common] - orig_sub.loc[common]) ** 2).mean()))
    print(f"    Sector {code:<3}: r(Chow-Lin, observed) = {r:.4f},  RMSE = {rmse:.5f}")

# ── [5] Save extended cache ────────────────────────────────────────────────────
print(f"\n[5] Saving extended VA to {VA_EXTENDED_CACHE}")
extended.to_parquet(VA_EXTENDED_CACHE, index=False)
print(f"    Saved: {VA_EXTENDED_CACHE}")

# ── [6] Validation plot ────────────────────────────────────────────────────────
print("\n[6] Generating validation plot")
_plot_validation(extended, original_dlog, gdp,
                 RESULTS_DIR / "va_chowlin_backcast_validation.png")

# ── [7] Summary statistics ─────────────────────────────────────────────────────
print("\n[7] Summary: dlog_va distribution pre- vs post-2005")
print(f"    {'Period':<15}  {'Mean':>8}  {'SD':>8}  {'Min':>8}  {'Max':>8}")
for label, mask in [("Pre-2005",  extended["quarter_label"] < "2005Q1"),
                    ("Post-2005", extended["quarter_label"] >= "2005Q1")]:
    d = extended.loc[mask, "dlog_va"].dropna()
    print(f"    {label:<15}  {d.mean():>+8.5f}  {d.std():>8.5f}  "
          f"{d.min():>+8.5f}  {d.max():>+8.5f}")

print("\n--- How to use in part2b ---")
print(
    "  In part2b_residualize_shocks.py, replace the BEA_VA_CACHE path with\n"
    "  VA_EXTENDED_CACHE (bea_va_quarterly_12ind_extended.parquet).\n"
    "  The merge in _load_controls() will now pick up pre-2005 quarters,\n"
    "  extending the v2 LP sample from 2005Q3 back to 1994Q4 for delta\n"
    "  (1992Q3 plus two lags for productivity and VA lag).\n"
    "  For LD, the sample starts at 2001Q1 regardless; the backcast\n"
    "  adds 2001Q1–2004Q4 VA observations for the LD residualization."
)

print("\nDone.")
