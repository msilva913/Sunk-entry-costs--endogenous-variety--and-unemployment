"""
part2b_residualize_shocks.py  --  Residualize national industry shock rates
==========================================================================
Constructs the residualized shock series ν^δ, ν^{TS}, ν^{LD}, ν^{QU} that
serve as inputs to the Bartik instrument aggregation in part3.

Economic motivation
-------------------
The raw shock rates g^k_{j,t} (establishment closings, layoffs, quits) are
correlated across shock types because all three series respond to common
aggregate and industry-specific demand conditions.  Simply differencing δ
from s at the Bartik stage does not remove this comovement.  We purge it
by projecting each log shock rate on predetermined controls and retaining
the residual, which is the closest empirical analog to the exogenous shock
draw in the model.

Residualization specification (v2)
------------------------------------
For δ and TS:
    log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t-1} + ν^k_{j,t}

For LD and QU (additionally controls for lagged tightness, since quits and
layoffs respond to labor market slack as well as demand conditions):
    log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t-1}
                        + μ log θ_{t-1} + ν^k_{j,t}

Regressors
----------
  Δlog p_t          aggregate nonfarm productivity growth  (FRED: OPHNFB)
  Δlog VA_{j,t-1}   real value-added growth in supersector j, lagged 1Q
                    (BEA via FRED; series listed in BLS_TO_FRED below)
  log θ_{t-1}       national log market tightness V/U, lagged 1Q
                    (FRED: JTSJOL / UNEMPLOY)

All regressors are either aggregate or lagged one quarter to ensure they
are predetermined with respect to the current-period shock draw.

The estimation uses a within-industry OLS (industry fixed effects absorbed
by demeaning) with homogeneous slope coefficients across industries.

Industry VA: FRED series IDs
------------------------------
  BLS 10 -> RVAM               Mining
  BLS 20 -> RVAC               Construction
  BLS 30 -> RVAMA              Manufacturing
  BLS 41 -> RVAW               Wholesale Trade
  BLS 42 -> RVAR               Retail Trade
  BLS 43 -> RVAT + RVAU        Transport+Warehousing + Utilities (summed)
  BLS 50 -> RVAI               Information
  BLS 55 -> RVAFI + RVARL      Finance+Insurance + Real Estate (summed)
  BLS 60 -> RVAPBS             Professional & Business Services
  BLS 65 -> RVAES      Education + Health & Social Assistance (summed)
  BLS 70 -> RVAAF      Arts/Entertainment + Accommodation/Food (summed)
  BLS 80 -> RVAOSEG            Other Services

Coverage: all FRED RVA series start 2005Q1; Chow-Lin extended cache
covers 1992Q1+. Sample for δ/TS residualization: 1997Q1+ (matching
BEA annual benchmark coverage; see part2b_extend_va.py). Sample for
LD/QU residualization: 2001Q1+ (JOLTS vacancies and tightness start
2001Q1). The two panels are built and residualized independently to
avoid JOLTS coverage silently clipping the δ sample.

Outputs
-------
  data/instruments/shock_rates_delta_resid.parquet   -- ν^δ by industry-quarter
  data/instruments/shock_rates_s_resid.parquet       -- ν^{TS}
  data/instruments/shock_rates_ld_resid.parquet      -- ν^{LD}
  data/instruments/shock_rates_qu_resid.parquet      -- ν^{QU}
  data/results/comovement_raw_series.png             -- raw log shock rates by industry
  data/results/comovement_scatter.png                -- raw vs. residualized scatter

Run
---
  FRED_API_KEY=<key> python part2b_shock_comovement.py

  On subsequent runs the caches are used and FRED_API_KEY is not needed.
  Free API key: https://fred.stlouisfed.org/docs/api/api_key.html

Prerequisites
-------------
  python part2_shock_rates.py
  python part2_shock_rates_s.py
"""

import os
import sys
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 8)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── working directory ─────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import (
    DEFAULT_CACHE_DIR, DEFAULT_OUTPUT_DIR, SHOCK_RATES_PATH, INDUSTRY_LABELS,
)
from construct_s_instrument import (
    SHOCK_RATES_S_PATH, SHOCK_RATES_LD_PATH, SHOCK_RATES_QU_PATH,
)

# ── paths ─────────────────────────────────────────────────────────────────────
PROD_CACHE      = DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet"
TIGHTNESS_CACHE = DEFAULT_CACHE_DIR / "national_tightness_quarterly.parquet"
_BEA_VA_EXTENDED = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind_extended.parquet"
_BEA_VA_ORIGINAL = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet"
# Prefer the Chow-Lin extended cache (covers 1992Q1+) when available;
# fall back to original FRED cache (2005Q1+) otherwise.
BEA_VA_CACHE    = _BEA_VA_EXTENDED if _BEA_VA_EXTENDED.exists() else _BEA_VA_ORIGINAL
RESID_D_PATH    = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_S_PATH    = DEFAULT_OUTPUT_DIR / "shock_rates_s_resid.parquet"
RESID_LD_PATH   = DEFAULT_OUTPUT_DIR / "shock_rates_ld_resid.parquet"
RESID_QU_PATH   = DEFAULT_OUTPUT_DIR / "shock_rates_qu_resid.parquet"
RESULTS_DIR     = Path("data/results")

for d in [DEFAULT_CACHE_DIR, DEFAULT_OUTPUT_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── FRED API key ──────────────────────────────────────────────────────────────
from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY")

# ── FRED series -> BLS supersector mapping ────────────────────────────────────
# Multi-component sectors are summed at the level before log-differencing.
# Corrections vs. prior version:
#   BLS 65: RVAESHS replaces RVAES (education-only subsector, ~16% of total).
#   BLS 70: RVAAER + RVAAF replaces RVAER (non-existent) + RVAAF.
BLS_TO_FRED = {
    10: ["RVAM"],            # Mining
    20: ["RVAC"],            # Construction
    30: ["RVAMA"],           # Manufacturing
    41: ["RVAW"],            # Wholesale Trade
    42: ["RVAR"],            # Retail Trade
    43: ["RVAT", "RVAU"],    # Transport+Warehousing + Utilities
    50: ["RVAI"],            # Information
    55: ["RVAFI", "RVARL"],  # Finance+Insurance + Real Estate+Rental
    60: ["RVAPBS"],          # Professional+Business Services
    65: ["RVAESHS"],         # Education+Health+Social Assistance (full aggregate)
    70: ["RVAAER", "RVAAF"], # Arts+Entertainment + Accommodation+Food
    80: ["RVAOSEG"],         # Other Services excl Govt
}

# ── small helpers ─────────────────────────────────────────────────────────────
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
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), alpha=0.10, color="grey")


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1: FRED FETCH HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _get_fred_client():
    """Return a Fred client, or raise clearly if the API key is missing."""
    if not FRED_API_KEY:
        raise EnvironmentError(
            "FRED_API_KEY environment variable is not set.\n"
            "  Register free at https://fred.stlouisfed.org/docs/api/api_key.html\n"
            "  Then: export FRED_API_KEY=<your_key>"
        )
    from fredapi import Fred
    return Fred(api_key=FRED_API_KEY)


def _fetch_fred_series(fred_client, sid: str) -> pd.Series:
    """Fetch one FRED series, resampled to quarter-start frequency."""
    try:
        s = fred_client.get_series(sid, observation_start="2004-10-01")
        s = s.resample("QS").first()
        if s.empty:
            raise ValueError(f"Series {sid} returned no data after resampling.")
        return s
    except Exception as e:
        raise ValueError(
            f"FRED series '{sid}' could not be fetched: {e}\n"
            f"  Verify at https://fred.stlouisfed.org/series/{sid}"
        ) from e


def _fetch_fred_va() -> pd.DataFrame:
    """
    Fetch quarterly real value-added by BLS supersector from FRED.
    Returns long DataFrame: quarter_label, industry_code (int), dlog_va.
    """
    fred_client = _get_fred_client()
    records, failed = [], []

    for bls_code, series_ids in BLS_TO_FRED.items():
        components = {}
        for sid in series_ids:
            try:
                components[sid] = _fetch_fred_series(fred_client, sid)
                print(f"    {sid}: {len(components[sid])} obs  OK")
            except ValueError as e:
                print(f"    {sid}: FAILED — {e}")
                failed.append((bls_code, sid, str(e)))

        if len(components) < len(series_ids):
            missing = [s for s in series_ids if s not in components]
            print(f"  BLS {bls_code}: SKIPPED (missing: {missing})")
            continue

        va_total = pd.concat(list(components.values()), axis=1).sum(
            axis=1, min_count=len(series_ids)
        )
        va_total.index = va_total.index.to_period("Q").strftime("%YQ%q")
        va_total.name  = "va"

        df = (va_total.dropna()
                      .reset_index()
                      .rename(columns={"index": "quarter_label"})
                      .sort_values("quarter_label")
                      .reset_index(drop=True))
        df["log_va"]        = np.log(df["va"].replace(0, np.nan))
        df["dlog_va"]       = df["log_va"].diff()
        df["industry_code"] = int(bls_code)
        records.append(df[["quarter_label", "industry_code", "dlog_va"]].dropna())

    if not records:
        raise ValueError(
            "FRED VA fetch returned no data for any sector.\n"
            "Check all series IDs in BLS_TO_FRED are valid FRED identifiers."
        )

    result  = pd.concat(records, ignore_index=True)
    n_sec   = result["industry_code"].nunique()
    print(f"\n  FRED VA fetch complete: {n_sec}/12 sectors, "
          f"{result['quarter_label'].min()}–{result['quarter_label'].max()}, "
          f"{len(result):,} rows")
    if n_sec < 12:
        missing_codes = set(BLS_TO_FRED) - set(result["industry_code"].unique())
        print(f"  WARNING: missing BLS codes {sorted(missing_codes)}. "
              "Enriched residualization will exclude these sectors.")
    return result


def _fetch_fred_tightness() -> pd.DataFrame:
    """
    Construct national log market tightness log(V/U) from FRED.
    V = JTSJOL (job openings, thousands); U = UNEMPLOY (persons, thousands).
    Returns DataFrame: quarter_label, log_theta.
    """
    fred_client = _get_fred_client()
    print("  Fetching JTSJOL ...")
    vac_m   = _fetch_fred_series(fred_client, "JTSJOL")
    print("  Fetching UNEMPLOY ...")
    unemp_m = _fetch_fred_series(fred_client, "UNEMPLOY")

    tight_m = pd.DataFrame({"vac": vac_m, "unemp": unemp_m}).dropna()
    tight_m.index = pd.to_datetime(tight_m.index)
    tight_m["quarter_label"] = tight_m.index.map(_dt_to_ql)
    tight = (tight_m.groupby("quarter_label")[["vac", "unemp"]]
                    .mean().reset_index())
    tight["theta"]     = tight["vac"] / tight["unemp"]
    tight["log_theta"] = np.log(tight["theta"])
    return tight[["quarter_label", "log_theta"]].copy()


def _build_tightness_from_local() -> pd.DataFrame:
    """
    Fallback: construct log(V/U) from already-cached state-level data.
    Used when FRED_API_KEY is not set but part4 caches exist.
    """
    vac_path  = DEFAULT_CACHE_DIR / "jolts_vacancies_state_quarterly.parquet"
    laus_path = DEFAULT_CACHE_DIR / "laus_states_monthly.parquet"
    if not vac_path.exists() or not laus_path.exists():
        raise FileNotFoundError(
            f"Cannot build tightness from local data: missing {vac_path} "
            f"or {laus_path}. Either set FRED_API_KEY or run part4_outcomes.py first."
        )
    vac_df  = pd.read_parquet(vac_path)
    nat_vac = (vac_df.groupby("quarter_label")["vacancies"]
                     .sum().reset_index()
                     .rename(columns={"vacancies": "vac_nat"}))
    laus_df = pd.read_parquet(laus_path)
    laus_df["date"]          = pd.to_datetime(laus_df["date"])
    laus_df["quarter_label"] = laus_df["date"].map(_dt_to_ql)
    laus_df["unemp"]         = laus_df["labor_force"] * laus_df["unemp_rate"] / 100
    nat_unemp = (laus_df.groupby("quarter_label")["unemp"]
                        .mean().reset_index()
                        .rename(columns={"unemp": "unemp_nat"}))
    t = nat_vac.merge(nat_unemp, on="quarter_label", how="inner")
    t["theta"]     = (t["vac_nat"] * 1000) / t["unemp_nat"]
    t["log_theta"] = np.log(t["theta"])
    return t[["quarter_label", "log_theta"]].copy()


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2: RESIDUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def _residualize(df, dep_col, resid_col, reg_cols):
    """
    Within-industry OLS of dep_col on industry FEs + reg_cols.
    Slopes are homogeneous across industries (pooled within estimator).
    Returns DataFrame with [industry_code, quarter_label, resid_col].
    """
    work = df[["industry_code", dep_col, *reg_cols]].copy().reset_index(drop=True)

    # Demean within industry
    work["y_dm"] = work[dep_col] - work.groupby("industry_code")[dep_col].transform("mean")
    dm_cols = []
    for rc in reg_cols:
        dc = f"_dm_{rc}"
        work[dc] = work[rc] - work.groupby("industry_code")[rc].transform("mean")
        dm_cols.append(dc)

    X      = work[dm_cols].to_numpy()
    y      = work["y_dm"].to_numpy()
    # OLS coefficients and residuals
    gammas = np.linalg.lstsq(X, y, rcond=None)[0]
    resid  = y - X @ gammas
    
    #R²
    r2 = 1 - (resid ** 2).sum() / (y ** 2).sum()
    coef_str = "  ".join(f"{rc}={g:+.4f}" for rc, g in zip(reg_cols, gammas))
    print(f"  {resid_col:<12}  {coef_str}   R²={r2:.4f}")

    out = df[["industry_code", "quarter_label"]].copy().reset_index(drop=True)
    out[resid_col] = resid
    return out


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3: DATA ASSEMBLY
#
#  Two separate panel builders with different sample ranges:
#
#    _load_panel_delta()  — δ and TS shocks
#        Sample: 1997Q1+ (BEA annual benchmarks available from 1997; BED
#        starts 1992Q3 but Chow-Lin VA extension is most credible from 1997).
#        Controls: Δlog p_t, Δlog VA_{j,t-1}.
#        NO tightness: JOLTS vacancies start 2001Q1; tightness is not in
#        the δ/TS residualization spec anyway.
#
#    _load_panel_ldqu()   — LD and QU shocks
#        Sample: 2001Q1+ (JOLTS separations and tightness both start here).
#        Controls: Δlog p_t, Δlog VA_{j,t-1}, log θ_{t-1}.
#
#  The key root cause of the prior clipping: all four shock series were
#  inner-joined in a single panel, silently limiting δ to the JOLTS start
#  date (2001Q1) even though BED data extends to 1992Q3.  Separating the
#  panels allows δ residualization to use its full pre-2001 sample.
# ═══════════════════════════════════════════════════════════════════════════════

# Shared helpers called by both panel builders
def _load_productivity() -> pd.DataFrame:
    """Aggregate nonfarm productivity growth (OPHNFB). Cached."""
    print("\n--- Aggregate productivity (OPHNFB) ---")
    if PROD_CACHE.exists():
        prod = pd.read_parquet(PROD_CACHE)
        print("  Loaded from cache")
    else:
        print("  Fetching from FRED ...")
        fred_client = _get_fred_client()
        s = fred_client.get_series("OPHNFB")
        prod = (s.rename("productivity").reset_index()
                  .rename(columns={"index": "date"}))
        prod["quarter_label"] = prod["date"].map(_dt_to_ql)
        prod = prod[["quarter_label", "productivity"]].dropna()
        prod.to_parquet(PROD_CACHE, index=False)
        print(f"  Saved: {PROD_CACHE}")
    prod = prod.sort_values("quarter_label").reset_index(drop=True)
    prod["dlog_p"] = np.log(prod["productivity"]).diff()
    return prod[["quarter_label", "dlog_p"]].dropna()


def _load_va_lag() -> pd.DataFrame:
    """Lagged industry VA growth (Δlog VA_{j,t-1}). Cached."""
    extended_label = (
        " [Chow-Lin extended, 1992Q1+]"
        if BEA_VA_CACHE == _BEA_VA_EXTENDED else " [FRED 2005Q1+ only]"
    )
    print(f"\n--- Industry real value-added{extended_label} ---")
    if BEA_VA_CACHE.exists():
        try:
            bea_va = pd.read_parquet(BEA_VA_CACHE)
            print(f"  Loaded: {bea_va['industry_code'].nunique()} sectors, "
                  f"{bea_va['quarter_label'].min()}–{bea_va['quarter_label'].max()}")
        except Exception as e:
            print(f"  Cache unreadable ({e}) — re-fetching from FRED.")
            BEA_VA_CACHE.unlink(missing_ok=True)
            bea_va = None
    else:
        bea_va = None

    if bea_va is None:
        print("  Fetching from FRED ...")
        try:
            bea_va = _fetch_fred_va()
            bea_va.to_parquet(BEA_VA_CACHE, index=False)
            print(f"  Saved: {BEA_VA_CACHE}")
        except Exception as e:
            sys.exit(f"ERROR: FRED industry VA fetch failed: {e}\n"
                     "Check your FRED_API_KEY.")

    bea_va = bea_va.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)
    bea_va["dlog_va_lag"] = bea_va.groupby("industry_code")["dlog_va"].shift(1)
    out = bea_va[["industry_code", "quarter_label", "dlog_va_lag"]].dropna().copy()
    out["industry_code"] = out["industry_code"].astype(int)
    return out


def _load_tightness_lag() -> pd.DataFrame:
    """Lagged log market tightness log(θ_{t-1}). Starts 2001Q2. Cached."""
    print("\n--- Market tightness log(V/U) ---")
    if TIGHTNESS_CACHE.exists():
        tight = pd.read_parquet(TIGHTNESS_CACHE)
        print("  Loaded from cache")
    else:
        print("  Fetching from FRED ...")
        try:
            tight = _fetch_fred_tightness()
        except Exception as e:
            sys.exit(f"ERROR: tightness fetch failed: {e}\n"
                     "Check your FRED_API_KEY.")
        tight.to_parquet(TIGHTNESS_CACHE, index=False)
        print(f"  Saved: {TIGHTNESS_CACHE}")
    tight = tight.sort_values("quarter_label").reset_index(drop=True)
    tight["log_theta_lag"] = tight["log_theta"].shift(1)
    return tight[["quarter_label", "log_theta_lag"]].dropna()


def _load_panel_delta(prod: pd.DataFrame, bea_lag: pd.DataFrame,
                      min_quarter: str = "1997Q1") -> pd.DataFrame:
    """
    Build the δ (and TS) residualization panel.

    Sample lower bound: min_quarter (default 1997Q1), matching the start of
    the BEA annual benchmark series used in the Chow-Lin backcast.  The BED
    data extends to 1992Q3 and the Chow-Lin VA is available from 1992Q2, but
    the annual-constrained VA growth is most reliable from 1997 onward where
    the annual benchmarks pin the level.  Using 1997Q1 is conservative and
    aligns the residualization sample with the annual constraint window.

    Controls: Δlog p_t, Δlog VA_{j,t-1}.
    Tightness is excluded — not in the δ/TS spec, and JOLTS vacancies
    (required to construct θ) start 2001Q1 anyway.
    """
    print("\n--- δ shock panel ---")
    if not SHOCK_RATES_PATH.exists():
        sys.exit(f"Missing: {SHOCK_RATES_PATH}\nRun part2_shock_rates.py first.")

    def _nat(raw, col, new_col):
        return (raw.groupby(["industry_code", "quarter_label"])[col]
                   .mean().reset_index()
                   .rename(columns={col: new_col}))

    nat = _nat(pd.read_parquet(SHOCK_RATES_PATH), "g_delta_loo", "g_delta")

    # TS (total separations) is JOLTS-based and starts 2001Q1 — it is NOT
    # joined here to avoid clipping the δ sample. TS is residualized in
    # _load_panel_ldqu() alongside LD and QU where it belongs.

    nat = nat[nat["g_delta"] > 0].copy()
    nat["log_g_delta"] = np.log(nat["g_delta"])
    nat["industry_label"] = (nat["industry_code"].map(INDUSTRY_LABELS)
                                                   .fillna(nat["industry_code"].astype(str)))
    nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)
    print(f"  Raw shock panel: {nat['quarter_label'].min()}–{nat['quarter_label'].max()}, "
          f"{nat['industry_code'].nunique()} industries, {len(nat):,} obs")

    # Attach controls (inner merges are safe: both prod and bea_lag cover 1992+)
    nat["industry_code"] = nat["industry_code"].astype(int)
    nat = nat.merge(prod,    on="quarter_label", how="inner")
    nat = nat.merge(bea_lag, on=["industry_code", "quarter_label"], how="inner")
    nat = nat.dropna(subset=["dlog_p", "dlog_va_lag"])

    # Apply sample lower bound
    nat = nat[nat["quarter_label"] >= min_quarter].copy()
    nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)
    print(f"  After controls + {min_quarter}+ filter: {len(nat):,} obs "
          f"({nat['quarter_label'].min()}–{nat['quarter_label'].max()})")
    return nat


def _load_panel_ldqu(prod: pd.DataFrame, bea_lag: pd.DataFrame,
                     tight_lag: pd.DataFrame) -> pd.DataFrame:
    """
    Build the LD, QU, and TS residualization panel.

    Sample: 2001Q1+ (JOLTS separations and tightness both start here).
    Controls: Δlog p_t, Δlog VA_{j,t-1}, log θ_{t-1}.

    TS (total separations) is included here — not in the δ panel — because
    it is also JOLTS-based and starts 2001Q1. Joining TS with δ would clip
    the δ sample to 2001Q1, which is the bug this split was introduced to fix.
    """
    print("\n--- TS / LD / QU shock panel ---")
    missing = [p for p in [SHOCK_RATES_S_PATH, SHOCK_RATES_LD_PATH,
                            SHOCK_RATES_QU_PATH] if not p.exists()]
    if missing:
        sys.exit("Missing: " + ", ".join(str(p) for p in missing) +
                 "\nRun part2_shock_rates_s.py first.")

    def _nat(raw, col, new_col):
        return (raw.groupby(["industry_code", "quarter_label"])[col]
                   .mean().reset_index()
                   .rename(columns={col: new_col}))

    nat_s  = _nat(pd.read_parquet(SHOCK_RATES_S_PATH),  "g_s_loo",  "g_s")
    nat_ld = _nat(pd.read_parquet(SHOCK_RATES_LD_PATH), "g_ld_loo", "g_ld")
    nat_qu = _nat(pd.read_parquet(SHOCK_RATES_QU_PATH), "g_qu_loo", "g_qu")

    nat = (nat_s
           .merge(nat_ld, on=["industry_code", "quarter_label"], how="inner")
           .merge(nat_qu, on=["industry_code", "quarter_label"], how="inner"))
    nat = nat[(nat["g_s"] > 0) & (nat["g_ld"] > 0) & (nat["g_qu"] > 0)].copy()
    nat["log_g_s"]  = np.log(nat["g_s"])
    nat["log_g_ld"] = np.log(nat["g_ld"])
    nat["log_g_qu"] = np.log(nat["g_qu"])
    nat["industry_label"] = (nat["industry_code"].map(INDUSTRY_LABELS)
                                                   .fillna(nat["industry_code"].astype(str)))
    nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)
    print(f"  Raw shock panel: {nat['quarter_label'].min()}–{nat['quarter_label'].max()}, "
          f"{nat['industry_code'].nunique()} industries, {len(nat):,} obs")

    nat["industry_code"] = nat["industry_code"].astype(int)
    nat = nat.merge(prod,      on="quarter_label", how="inner")
    nat = nat.merge(bea_lag,   on=["industry_code", "quarter_label"], how="inner")
    nat = nat.merge(tight_lag, on="quarter_label", how="inner")
    nat = nat.dropna(subset=["dlog_p", "dlog_va_lag", "log_theta_lag"])
    nat = nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)
    print(f"  After controls merge: {len(nat):,} obs "
          f"({nat['quarter_label'].min()}–{nat['quarter_label'].max()})")
    return nat


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 65)
print("Part 2b  --  Residualize national industry shock rates")
print("=" * 65)
if not FRED_API_KEY:
    sys.exit(
        "ERROR: FRED_API_KEY is not set.\n"
        "  Register free at https://fred.stlouisfed.org/docs/api/api_key.html\n"
        "  Windows PowerShell:  $env:FRED_API_KEY = '<your_key>'\n"
        "  Linux/macOS:         export FRED_API_KEY=<your_key>"
    )
print(f"FRED_API_KEY: set ({FRED_API_KEY[:8]}...)")

# ── 1. Load shared controls (each fetched/cached once) ────────────────────────
prod     = _load_productivity()
bea_lag  = _load_va_lag()
tight_lag = _load_tightness_lag()

# ── 2. Build separate panels for δ/TS vs LD/QU ───────────────────────────────
#
# δ panel: BED data from 1997Q1+ (matching annual benchmark coverage).
#   No tightness control — JOLTS vacancies start 2001Q1 and tightness is not
#   in the δ spec. This is the key fix: previously, joining δ with LD/QU
#   clipped the entire panel to 2001Q1 even though BED covers 1992Q3+.
#
# LD/QU panel: JOLTS data from 2001Q1+. Includes tightness control.
nat_d  = _load_panel_delta(prod, bea_lag, min_quarter="1997Q1")
nat_sq = _load_panel_ldqu(prod, bea_lag, tight_lag)

# ── 3. Residualize ────────────────────────────────────────────────────────────
SPEC      = "v2: Δlog(p) + Δlog(VA_lag) + log(θ_lag) [LD/QU only]"
reg_base  = ("dlog_p", "dlog_va_lag")
reg_tight = ("dlog_p", "dlog_va_lag", "log_theta_lag")

print(f"\n--- Residualization [{SPEC}] ---")
print(f"  δ    : controls = {reg_base}  sample = {nat_d['quarter_label'].min()}–{nat_d['quarter_label'].max()}")
print(f"  TS/LD/QU: controls = {reg_tight}  sample = {nat_sq['quarter_label'].min()}–{nat_sq['quarter_label'].max()}")

resid_d  = _residualize(nat_d,  "log_g_delta", "nu_delta", reg_base)
resid_s  = _residualize(nat_sq, "log_g_s",     "nu_s",     reg_base)   # TS: no tightness
resid_ld = _residualize(nat_sq, "log_g_ld",    "nu_ld",    reg_tight)
resid_qu = _residualize(nat_sq, "log_g_qu",    "nu_qu",    reg_tight)

nat_d["nu_delta"]  = resid_d["nu_delta"].to_numpy()
nat_sq["nu_s"]     = resid_s["nu_s"].to_numpy()
nat_sq["nu_ld"]    = resid_ld["nu_ld"].to_numpy()
nat_sq["nu_qu"]    = resid_qu["nu_qu"].to_numpy()

for panel, col in [(nat_d,"nu_delta"),(nat_sq,"nu_s"),(nat_sq,"nu_ld"),(nat_sq,"nu_qu")]:
    n_null = panel[col].isna().sum()
    print(f"  {col}: {'WARNING: '+str(n_null)+' nulls' if n_null else 'OK'}")

# ── 4. Save residuals ─────────────────────────────────────────────────────────
nat_d[["industry_code","quarter_label","nu_delta"]].to_parquet(RESID_D_PATH,  index=False)
nat_sq[["industry_code","quarter_label","nu_s"]   ].to_parquet(RESID_S_PATH,  index=False)
nat_sq[["industry_code","quarter_label","nu_ld"]  ].to_parquet(RESID_LD_PATH, index=False)
nat_sq[["industry_code","quarter_label","nu_qu"]  ].to_parquet(RESID_QU_PATH, index=False)
for p in [RESID_D_PATH, RESID_S_PATH, RESID_LD_PATH, RESID_QU_PATH]:
    print(f"Saved: {p}")
print(f"\nResiduals saved [{SPEC}].")

# ── 5. Cross-series correlations ─────────────────────────────────────────────
# Note: δ and LD/QU are on different sample frames. For cross-series
# comparisons we merge on the overlapping 2001Q1+ window.
nat_overlap = (
    nat_d[["industry_code","quarter_label","log_g_delta","nu_delta"]]
    .merge(
        nat_sq[["industry_code","quarter_label",
                "log_g_s","log_g_ld","log_g_qu","nu_s","nu_ld","nu_qu"]],
        on=["industry_code","quarter_label"], how="inner"
    )
)
print(f"\n  Cross-series overlap sample: "
      f"{nat_overlap['quarter_label'].min()}–{nat_overlap['quarter_label'].max()}, "
      f"{len(nat_overlap):,} obs")

pairs = [
    ("delta vs TS", "log_g_delta","log_g_s",  "nu_delta","nu_s"),
    ("delta vs LD", "log_g_delta","log_g_ld", "nu_delta","nu_ld"),
    ("delta vs QU", "log_g_delta","log_g_qu", "nu_delta","nu_qu"),
    ("LD vs QU",    "log_g_ld",   "log_g_qu", "nu_ld",   "nu_qu"),
]
print(f"\n--- Pooled cross-series correlations [{SPEC}] ---")
print(f"  {'Pair':<14}  {'r raw':>8}  {'r resid':>8}  {'reduction':>10}")
print(f"  {'-'*48}")
for lbl, rc1, rc2, rc3, rc4 in pairs:
    r_raw   = nat_overlap[rc1].corr(nat_overlap[rc2])
    r_resid = nat_overlap[rc3].corr(nat_overlap[rc4])
    print(f"  {lbl:<14}  {r_raw:>+8.4f}  {r_resid:>+8.4f}  {r_raw-r_resid:>+10.4f}")

def _by_industry_r(df, c1, c2):
    return (df.groupby("industry_label", group_keys=False)
              .apply(lambda g: pd.Series({"r": g[c1].corr(g[c2])}),
                     include_groups=False)
              .reset_index().sort_values("r"))

# Per-industry correlations computed on the 2001Q1+ overlap window
nat_overlap["industry_label"] = (
    nat_overlap["industry_code"].map(INDUSTRY_LABELS)
    .fillna(nat_overlap["industry_code"].astype(str))
)
print(f"\n--- Per-industry corr(delta vs LD) and corr(delta vs QU) ---")
by_ind = (_by_industry_r(nat_overlap,"log_g_delta","log_g_ld").rename(columns={"r":"r_dLD_raw"})
           .merge(_by_industry_r(nat_overlap,"nu_delta","nu_ld").rename(columns={"r":"r_dLD_resid"}),
                  on="industry_label")
           .merge(_by_industry_r(nat_overlap,"log_g_delta","log_g_qu").rename(columns={"r":"r_dQU_raw"}),
                  on="industry_label")
           .merge(_by_industry_r(nat_overlap,"nu_delta","nu_qu").rename(columns={"r":"r_dQU_resid"}),
                  on="industry_label")
           .sort_values("r_dLD_raw"))
print(by_ind.to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════════
#  DIAGNOSTICS  (plots — can be skipped without affecting outputs)
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_diagnostics():
    plt.ioff()   # suppress interactive display; avoids REPL repr spam

    # Plot A: raw log shock rates by industry
    # δ/TS use nat_d (1997Q1+); LD/QU use nat_sq (2001Q1+)
    industries = sorted(nat_d["industry_label"].unique())
    n_ind      = len(industries)
    series_defs = [
        (nat_d,  "log_g_delta", r"$\log g^\delta$", "#1f77b4"),
        (nat_sq, "log_g_s",     r"$\log g^{TS}$",   "#d62728"),
        (nat_sq, "log_g_ld",    r"$\log g^{LD}$",   "#2ca02c"),
        (nat_sq, "log_g_qu",    r"$\log g^{QU}$",   "#ff7f0e"),
    ]
    fig, axes = plt.subplots(n_ind, 4, figsize=(20, 2.0 * n_ind),
                             sharex=False, squeeze=False)
    for i, ind in enumerate(industries):
        for j, (src, col, _ttl, color) in enumerate(series_defs):
            ax = axes[i, j]
            sub = src[src["industry_label"] == ind].sort_values("quarter_label")
            dates_j = [_ql_to_dt(q) for q in sub["quarter_label"]]
            _ = ax.plot(dates_j, sub[col].values, color=color, linewidth=1.1)
            _add_recessions(ax)
            _ = ax.set_ylabel(ind if j == 0 else "", fontsize=7)
            _ = ax.grid(axis="y", linewidth=0.4, alpha=0.4)
    for j, (_, _, title, _c) in enumerate(series_defs):
        _ = axes[0, j].set_title(title, fontsize=10)
    _ = fig.suptitle(
        "National industry shock rates by supersector (raw log)\n"
        "δ/TS: 1997Q1+  |  LD/QU: 2001Q1+",
        fontsize=11, y=1.005
    )
    fig.tight_layout()
    p = RESULTS_DIR / "comovement_raw_series.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    _open_file(p)
    print(f"\nPlot A saved: {p}")

    # Plot B: δ vs LD and δ vs QU — raw vs residualized (2001Q1+ overlap)
    scatter_pairs = [
        ("delta vs TS", "log_g_delta","log_g_s",  "nu_delta","nu_s",
         r"$\log g^\delta$",r"$\log g^{TS}$",r"$\nu^\delta$",r"$\nu^{TS}$"),
        ("delta vs LD", "log_g_delta","log_g_ld", "nu_delta","nu_ld",
         r"$\log g^\delta$",r"$\log g^{LD}$",r"$\nu^\delta$",r"$\nu^{LD}$"),
        ("delta vs QU", "log_g_delta","log_g_qu", "nu_delta","nu_qu",
         r"$\log g^\delta$",r"$\log g^{QU}$",r"$\nu^\delta$",r"$\nu^{QU}$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ci, (lbl, rx, ry, nx, ny, rxl, ryl, nxl, nyl) in enumerate(scatter_pairs):
        r_r = nat_overlap[rx].corr(nat_overlap[ry])
        r_n = nat_overlap[nx].corr(nat_overlap[ny])
        for row, (cx, cy, xl, yl, color, r_val) in enumerate([
            (rx, ry, rxl, ryl, "#555",    r_r),
            (nx, ny, nxl, nyl, "#1f77b4", r_n),
        ]):
            ax = axes[row, ci]
            _ = ax.scatter(nat_overlap[cx], nat_overlap[cy],
                           s=5, alpha=0.35, color=color)
            _ = ax.set_xlabel(xl, fontsize=9)
            _ = ax.set_ylabel(yl, fontsize=9)
            _ = ax.set_title(
                f"{lbl}  {'raw' if row==0 else 'residualized'}  (r={r_val:.3f})",
                fontsize=9
            )
            _ = ax.axhline(0, color="black", linewidth=0.5)
            _ = ax.axvline(0, color="black", linewidth=0.5)
            _ = ax.grid(linewidth=0.4, alpha=0.4)
    _ = fig.suptitle(
        rf"$\delta$ vs s-type: raw (top) vs residualized (bottom) [{SPEC}]"
        "\n(overlap sample: 2001Q1+)",
        fontsize=11,
    )
    fig.tight_layout()
    p = RESULTS_DIR / "comovement_scatter.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    _open_file(p)
    print(f"Plot B saved: {p}")


_plot_diagnostics()
print("\nDone.")