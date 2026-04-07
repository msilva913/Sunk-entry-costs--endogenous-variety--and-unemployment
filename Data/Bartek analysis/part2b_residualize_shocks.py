"""
part2b_shock_comovement.py  --  Residualize national industry shock rates
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
  BLS 65 -> RVAES + RVAHC      Education + Health & Social Assistance (summed)
  BLS 70 -> RVAER + RVAAF      Arts/Entertainment + Accommodation/Food (summed)
  BLS 80 -> RVAOSEG            Other Services

Coverage: all FRED RVA series start 2005Q1.  Sample for this script is
2005Q2 onward (one lag needed).  LP sample in part5 will be 2005Q2–2019Q4.

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
BEA_VA_CACHE    = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet"
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

def _residualize(df, dep_col, resid_col, reg_cols, heterogeneous_slopes=True):
    """
    Residualize dep_col on industry FEs + reg_cols.

    If heterogeneous_slopes=True (default): runs separate OLS within each
    industry, allowing every slope coefficient to differ across industries.
    This absorbs industry-specific sensitivity to aggregate and demand
    controls, leaving residuals that reflect structural shock variation.

    If heterogeneous_slopes=False: pools all industries after within-demeaning
    (original behavior — one set of slope coefficients for all industries).

    Returns DataFrame with [industry_code, quarter_label, resid_col].
    """
    work = df[["industry_code", "quarter_label", dep_col, *list(reg_cols)]].copy().reset_index(drop=True)
    industries = sorted(work["industry_code"].unique())

    if heterogeneous_slopes:
        # ── separate OLS per industry ──────────────────────────────────────
        records = []
        coef_rows = []

        for ind in industries:
            sub = work[work["industry_code"] == ind].copy()

            # Drop rows with any missing values in this industry
            sub = sub.dropna(subset=[dep_col] + list(reg_cols))
            if len(sub) <= len(reg_cols) + 1:
                print(f"  WARNING: {ind} has only {len(sub)} obs — skipping")
                continue

            y = sub[dep_col].to_numpy()
            X = np.column_stack([np.ones(len(sub)), sub[list(reg_cols)].to_numpy()])

            # OLS: [intercept, γ_j, λ_j, μ_j, ...]
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            resid = y - X @ coeffs

            # R² within this industry
            r2 = 1 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum()

            coef_str = "  ".join(
                f"{rc}={c:+.4f}" for rc, c in zip(reg_cols, coeffs[1:])
            )
            print(f"  {resid_col:<12}  industry={str(ind):<6}  "
                  f"{coef_str}  R²={r2:.4f}")

            coef_rows.append({
                "resid_col":    resid_col,
                "industry_code": ind,
                **{rc: c for rc, c in zip(reg_cols, coeffs[1:])},
                "R2": r2,
            })

            sub = sub.copy()
            sub[resid_col] = resid
            records.append(sub[["industry_code", "quarter_label", resid_col]])

        out = pd.concat(records, ignore_index=True)

        # Print cross-industry coefficient summary
        coef_df = pd.DataFrame(coef_rows)
        print(f"\n  --- Cross-industry coefficient ranges [{resid_col}] ---")
        for rc in reg_cols:
            print(f"    {rc:<20}  "
                  f"mean={coef_df[rc].mean():+.4f}  "
                  f"min={coef_df[rc].min():+.4f}  "
                  f"max={coef_df[rc].max():+.4f}  "
                  f"std={coef_df[rc].std():.4f}")

    else:
        # ── original pooled within estimator ──────────────────────────────
        work["y_dm"] = (work[dep_col]
                        - work.groupby("industry_code")[dep_col].transform("mean"))
        dm_cols = []
        for rc in reg_cols:
            dc = f"_dm_{rc}"
            work[dc] = (work[rc]
                        - work.groupby("industry_code")[rc].transform("mean"))
            dm_cols.append(dc)

        X      = work[dm_cols].to_numpy()
        y      = work["y_dm"].to_numpy()
        gammas = np.linalg.lstsq(X, y, rcond=None)[0]
        resid  = y - X @ gammas

        r2 = 1 - (resid ** 2).sum() / (y ** 2).sum()
        coef_str = "  ".join(
            f"{rc}={g:+.4f}" for rc, g in zip(reg_cols, gammas)
        )
        print(f"  {resid_col:<12}  {coef_str}   R²={r2:.4f}")

        out = df[["industry_code", "quarter_label"]].copy().reset_index(drop=True)
        out[resid_col] = resid

    return out


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3: DATA ASSEMBLY
#  Loads all shock series and attaches predetermined controls.
#  Returns a single merged panel ready for residualization.
# ═══════════════════════════════════════════════════════════════════════════════

def _load_controls() -> pd.DataFrame:
    """
    Build the national industry-quarter panel with shock series and controls.
    Steps:
      (a) Load raw shock rates from part2 / part2_s parquets.
      (b) Attach aggregate productivity growth (Δlog p).
      (c) Attach lagged log market tightness (log θ_{t-1}).
      (d) Attach lagged industry VA growth (Δlog VA_{j,t-1}).
    Returns merged DataFrame ready to pass to _residualize().
    """
    # ── (a) shock rates ───────────────────────────────────────────────────────
    print("--- Loading shock rates ---")
    missing = [p for p in [SHOCK_RATES_PATH, SHOCK_RATES_S_PATH,
                            SHOCK_RATES_LD_PATH, SHOCK_RATES_QU_PATH]
               if not p.exists()]
    if missing:
        sys.exit("Missing shock rate files:\n" +
                 "\n".join(f"  {p}" for p in missing) +
                 "\nRun part2_shock_rates.py and part2_shock_rates_s.py first.")

    def _nat(raw, col, new_col):
        return (raw.groupby(["industry_code", "quarter_label"])[col]
                   .mean().reset_index()
                   .rename(columns={col: new_col}))

    nat_d  = _nat(pd.read_parquet(SHOCK_RATES_PATH),    "g_delta_loo", "g_delta")
    nat_s  = _nat(pd.read_parquet(SHOCK_RATES_S_PATH),  "g_s_loo",     "g_s")
    nat_ld = _nat(pd.read_parquet(SHOCK_RATES_LD_PATH), "g_ld_loo",    "g_ld")
    nat_qu = _nat(pd.read_parquet(SHOCK_RATES_QU_PATH), "g_qu_loo",    "g_qu")

    nat = (nat_d
           .merge(nat_s,  on=["industry_code", "quarter_label"], how="inner")
           .merge(nat_ld, on=["industry_code", "quarter_label"], how="inner")
           .merge(nat_qu, on=["industry_code", "quarter_label"], how="inner"))
    nat = nat[(nat["g_delta"] > 0) & (nat["g_s"] > 0) &
              (nat["g_ld"]   > 0) & (nat["g_qu"] > 0)].copy()
    for col, label in [("g_delta","log_g_delta"),("g_s","log_g_s"),
                       ("g_ld","log_g_ld"),("g_qu","log_g_qu")]:
        nat[label] = np.log(nat[col])
    nat["industry_label"] = (nat["industry_code"].map(INDUSTRY_LABELS)
                                                  .fillna(nat["industry_code"].astype(str)))
    nat = nat.sort_values(["industry_code","quarter_label"]).reset_index(drop=True)
    print(f"  Shock panel: {nat['quarter_label'].min()}–{nat['quarter_label'].max()}, "
          f"{nat['industry_code'].nunique()} industries, {len(nat):,} obs")

    # ── (b) aggregate productivity ────────────────────────────────────────────
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
    prod["log_p"]  = np.log(prod["productivity"])
    prod["dlog_p"] = prod["log_p"].diff()
    nat = nat.merge(prod[["quarter_label","dlog_p"]].dropna(),
                    on="quarter_label", how="inner")
    print(f"  After merge: {len(nat):,} obs")

    # ── (c) lagged market tightness ───────────────────────────────────────────
    print("\n--- Market tightness log(V/U) ---")
    if TIGHTNESS_CACHE.exists():
        tight = pd.read_parquet(TIGHTNESS_CACHE)
        print("  Loaded from cache")
    else:
        print("  Fetching from FRED ...")
        try:
            tight = _fetch_fred_tightness()
        except Exception as e:
            sys.exit(f"ERROR: FRED tightness fetch failed: {e}\n"
                     "Check your FRED_API_KEY and network connection.")
        tight.to_parquet(TIGHTNESS_CACHE, index=False)
        print(f"  Saved: {TIGHTNESS_CACHE}")

    tight = tight.sort_values("quarter_label").reset_index(drop=True)
    tight["log_theta_lag"] = tight["log_theta"].shift(1)
    nat = nat.merge(tight[["quarter_label","log_theta_lag"]].dropna(),
                    on="quarter_label", how="inner")
    print(f"  After merge: {len(nat):,} obs")

    # ── (d) lagged industry value-added growth ────────────────────────────────
    print("\n--- Industry real value-added (BEA via FRED) ---")
    if BEA_VA_CACHE.exists():
        try:
            bea_va = pd.read_parquet(BEA_VA_CACHE)
            print(f"  Loaded from cache: {bea_va['industry_code'].nunique()} sectors, "
                  f"{bea_va['quarter_label'].min()}–{bea_va['quarter_label'].max()}")
        except Exception as e:
            print(f"  Cache unreadable ({e}) — will re-fetch.")
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
                     "Check your FRED_API_KEY and network connection.")

    bea_va = bea_va.sort_values(["industry_code","quarter_label"]).reset_index(drop=True)
    bea_va["dlog_va_lag"] = bea_va.groupby("industry_code")["dlog_va"].shift(1)
    bea_lag = bea_va[["industry_code","quarter_label","dlog_va_lag"]].dropna().copy()
    bea_lag["industry_code"] = bea_lag["industry_code"].astype(int)
    nat["industry_code"]     = nat["industry_code"].astype(int)
    nat = nat.merge(bea_lag, on=["industry_code","quarter_label"], how="inner")
    nat = nat.dropna(subset=["dlog_va_lag"])
    print(f"  After merge: {len(nat):,} obs "
          f"({nat['quarter_label'].min()}–{nat['quarter_label'].max()})")

    return nat.sort_values(["industry_code","quarter_label"]).reset_index(drop=True)


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

# ── 1. assemble controls panel ────────────────────────────────────────────────
nat = _load_controls()
#Inspect
nat.groupby(['industry_code']).head(5)
# ── 2. residualize ────────────────────────────────────────────────────────────
SPEC      = "v2: Δlog(p) + Δlog(VA_lag) + log(θ_lag) [LD/QU only]"
reg_base  = ("dlog_p", "dlog_va_lag")
reg_tight = ("dlog_p", "dlog_va_lag", "log_theta_lag")

print(f"\n--- Residualization [{SPEC}] ---")
print(f"  δ, TS : controls = {reg_base}")
print(f"  LD, QU: controls = {reg_tight}")

resid_d  = _residualize(nat, "log_g_delta", "nu_delta", reg_base)
resid_s  = _residualize(nat, "log_g_s",     "nu_s",     reg_base)
resid_ld = _residualize(nat, "log_g_ld",    "nu_ld",    reg_tight)
resid_qu = _residualize(nat, "log_g_qu",    "nu_qu",    reg_tight)

nat["nu_delta"] = resid_d["nu_delta"].to_numpy()
nat["nu_s"]     = resid_s["nu_s"].to_numpy()
nat["nu_ld"]    = resid_ld["nu_ld"].to_numpy()
nat["nu_qu"]    = resid_qu["nu_qu"].to_numpy()

for col in ["nu_delta","nu_s","nu_ld","nu_qu"]:
    n_null = nat[col].isna().sum()
    print(f"  {col}: {'WARNING: '+str(n_null)+' nulls' if n_null else 'OK'}")

# ── 3. save residuals ─────────────────────────────────────────────────────────
nat[["industry_code","quarter_label","nu_delta"]].to_parquet(RESID_D_PATH,  index=False)
nat[["industry_code","quarter_label","nu_s"]    ].to_parquet(RESID_S_PATH,  index=False)
nat[["industry_code","quarter_label","nu_ld"]   ].to_parquet(RESID_LD_PATH, index=False)
nat[["industry_code","quarter_label","nu_qu"]   ].to_parquet(RESID_QU_PATH, index=False)
for p in [RESID_D_PATH, RESID_S_PATH, RESID_LD_PATH, RESID_QU_PATH]:
    print(f"Saved: {p}")
print(f"\nResiduals saved [{SPEC}].")

# ── 4. cross-series correlation summary ──────────────────────────────────────
pairs = [
    ("delta vs TS", "log_g_delta","log_g_s",  "nu_delta","nu_s"),
    ("delta vs LD", "log_g_delta","log_g_ld", "nu_delta","nu_ld"),
    ("delta vs QU", "log_g_delta","log_g_qu", "nu_delta","nu_qu"),
    ("LD vs QU",    "log_g_ld",  "log_g_qu", "nu_ld",   "nu_qu"),
    ("TS vs LD",    "log_g_s",   "log_g_ld", "nu_s",    "nu_ld"),
    ("TS vs QU",    "log_g_s",   "log_g_qu", "nu_s",    "nu_qu"),
]
print(f"\n--- Pooled cross-series correlations [{SPEC}] ---")
print(f"  {'Pair':<14}  {'r raw':>8}  {'r resid':>8}  {'reduction':>10}")
print(f"  {'-'*48}")
for lbl, rc1, rc2, rc3, rc4 in pairs:
    r_raw   = nat[rc1].corr(nat[rc2])
    r_resid = nat[rc3].corr(nat[rc4])
    print(f"  {lbl:<14}  {r_raw:>+8.4f}  {r_resid:>+8.4f}  {r_raw-r_resid:>+10.4f}")

print(f"\n--- Per-industry corr(delta vs LD) and corr(delta vs QU) ---")

def _by_industry_r(df, c1, c2):
    return (df.groupby("industry_label", group_keys=False)
              .apply(lambda g: pd.Series({"r": g[c1].corr(g[c2])}),
                     include_groups=False)
              .reset_index().sort_values("r"))

by_ind = (_by_industry_r(nat,"log_g_delta","log_g_ld").rename(columns={"r":"r_dLD_raw"})
           .merge(_by_industry_r(nat,"nu_delta","nu_ld").rename(columns={"r":"r_dLD_resid"}),
                  on="industry_label")
           .merge(_by_industry_r(nat,"log_g_delta","log_g_qu").rename(columns={"r":"r_dQU_raw"}),
                  on="industry_label")
           .merge(_by_industry_r(nat,"nu_delta","nu_qu").rename(columns={"r":"r_dQU_resid"}),
                  on="industry_label")
           .sort_values("r_dLD_raw"))
print(by_ind.to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════════
#  DIAGNOSTICS  (plots — can be skipped without affecting outputs)
# ═══════════════════════════════════════════════════════════════════════════════

industries = sorted(nat["industry_label"].unique())
quarters   = sorted(nat["quarter_label"].unique())
dates      = [_ql_to_dt(q) for q in quarters]
n_ind      = len(industries)

# Plot A: raw log shock rates by industry
series_cols = [
    ("log_g_delta", r"$\log g^\delta$", "#1f77b4"),
    ("log_g_s",     r"$\log g^{TS}$",   "#d62728"),
    ("log_g_ld",    r"$\log g^{LD}$",   "#2ca02c"),
    ("log_g_qu",    r"$\log g^{QU}$",   "#ff7f0e"),
]
fig, axes = plt.subplots(n_ind, 4, figsize=(20, 2.0 * n_ind),
                         sharex=True, squeeze=False)
for i, ind in enumerate(industries):
    sub = nat[nat["industry_label"] == ind].set_index("quarter_label").reindex(quarters)
    for j, (col, _, color) in enumerate(series_cols):
        ax = axes[i, j]
        ax.plot(dates, sub[col].values, color=color, linewidth=1.1)
        _add_recessions(ax)
        ax.set_ylabel(ind if j == 0 else "", fontsize=7)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
for j, (_, title, _) in enumerate(series_cols):
    axes[0, j].set_title(title, fontsize=10)
fig.suptitle("National industry shock rates by supersector (raw log)", fontsize=11, y=1.005)
fig.tight_layout()
p = RESULTS_DIR / "comovement_raw_series.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot A saved: {p}")

# Plot B: delta vs each s-type — raw (top row) vs residualized (bottom row)
scatter_pairs = [
    ("delta vs TS", "log_g_delta","log_g_s",  "nu_delta","nu_s",
     r"$\log g^\delta$",r"$\log g^{TS}$",r"$\nu^\delta$",r"$\nu^{TS}$"),
    ("delta vs LD", "log_g_delta","log_g_ld", "nu_delta","nu_ld",
     r"$\log g^\delta$",r"$\log g^{LD}$",r"$\nu^\delta$",r"$\nu^{LD}$"),
    ("delta vs QU", "log_g_delta","log_g_qu", "nu_delta","nu_qu",
     r"$\log g^\delta$",r"$\log g^{QU}$",r"$\nu^\delta$",r"$\nu^{QU}$"),
]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ci, (lbl, rx, ry, nx, ny, rxl, ryl, nxl, nyl) in enumerate(scatter_pairs):
    r_r = nat[rx].corr(nat[ry])
    r_n = nat[nx].corr(nat[ny])
    for row, (cx, cy, xl, yl, color, r_val) in enumerate([
        (rx, ry, rxl, ryl, "#555",    r_r),
        (nx, ny, nxl, nyl, "#1f77b4", r_n),
    ]):
        ax = axes[row, ci]
        ax.scatter(nat[cx], nat[cy], s=5, alpha=0.35, color=color)
        ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
        ax.set_title(f"{lbl}  {'raw' if row==0 else 'residualized'}  (r={r_val:.3f})",
                     fontsize=9)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.grid(linewidth=0.4, alpha=0.4)
fig.suptitle(
    rf"$\delta$ vs each s-type: raw (top) vs residualized (bottom) [{SPEC}]",
    fontsize=11,
)
fig.tight_layout()
p = RESULTS_DIR / "comovement_scatter.png"
fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot B saved: {p}")

print("\nDone.")
