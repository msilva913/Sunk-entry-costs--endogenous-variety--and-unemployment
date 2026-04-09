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

Residualization specification (v3)
------------------------------------
For δ and TS:
    log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t-1}
                        + ψ (φ_j × ΔMP_t) + ν^k_{j,t}

For LD and QU (additionally controls for lagged tightness):
    log g^k_{j,t} = α_j + γ Δlog p_t + λ Δlog VA_{j,t-1}
                        + μ log θ_{t-1} + ψ (φ_j × ΔMP_t) + ν^k_{j,t}

Regressors
----------
  Δlog p_t          aggregate nonfarm productivity growth  (FRED: OPHNFB)
  Δlog VA_{j,t-1}   real value-added growth in supersector j, lagged 1Q
                    (BEA via FRED; series listed in BLS_TO_FRED below)
  log θ_{t-1}       national log market tightness V/U, lagged 1Q
                    (FRED: JTSJOL / UNEMPLOY)
  φ_j × ΔMP_t       monetary-policy × establishment-size interaction  [v3]
                    φ_j = log avg establishment size from Census CBP 2000
                    ΔMP_t = quarterly Δ Wu-Xia shadow FFR (Atlanta Fed),
                    extended with FEDFUNDS (FRED) post-2022Q1 — theoretically
                    exact since shadow rate = FFR when FFR > 0.25%.

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
  Wu-Xia and CBP caches require no API keys.

Prerequisites
-------------
  python part2_shock_rates.py
  python part2_shock_rates_s.py
"""

import io
import os
import sys
import numpy as np
import pandas as pd
import requests
pd.set_option('display.max_columns', 8)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.ioff()   # suppress interactive display; prevents REPL from printing repr of artists
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
WUXIA_CACHE     = DEFAULT_CACHE_DIR / "wuxia_quarterly.parquet"        # v3
ESTAB_SZ_CACHE  = DEFAULT_CACHE_DIR / "cbp_estab_size.parquet"         # v3
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

# ── v3: Census CBP NAICS → BLS supersector mapping (year 2000, NAICS1997) ────
# Multi-NAICS supersectors are queried separately and summed (emp + estabs).
BLS_TO_CBP_NAICS = {
    10: ["21"],
    20: ["23"],
    30: ["31-33"],
    41: ["42"],
    42: ["44-45"],
    43: ["48-49", "22"],
    50: ["51"],
    55: ["52", "53"],
    60: ["54", "55", "56"],
    65: ["61", "62"],
    70: ["71", "72"],
    80: ["81"],
}
CBP_YEAR      = 2000          # pre-sample reference year (predetermined φ_j)
CBP_NAICS_VAR = "NAICS1997"   # Census CBP variable name for year 2000

# ── v3: Wu-Xia shadow FFR source ──────────────────────────────────────────────
WUXIA_URL = (
    "https://www.atlantafed.org/-/media/Project/Atlanta/FRBA/Documents/"
    "datafiles/cqer/research/wu-xia-shadow-federal-funds-rate/WuXiaShadowRate.xlsx"
)

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


# ── v3: Wu-Xia shadow federal funds rate ──────────────────────────────────────

def _fetch_wuxia() -> pd.DataFrame:
    """
    Download Wu-Xia (2016) shadow FFR from Atlanta Fed and extend with
    FEDFUNDS from FRED for the post-suspension period.

    Excel layout (confirmed):
      Col 0: date (monthly, last business day)
      Col 1: effective federal funds rate (% p.a.)
      Col 2: Wu-Xia shadow FFR — non-null from 1990M1; NaN before 1990
      Row 0:  column headers; data from row 1 onward.
      Last observation: 2022M2 (updates suspended when FOMC exited ZLB).

    Strategy:
      - Within Excel coverage: use Wu-Xia where non-null (1990Q1–2022Q1);
        fall back to FEDFUNDS column for pre-1990 quarters.
      - Post-2022M2: extend with FEDFUNDS from FRED. This is theoretically
        exact: Wu-Xia model states shadow rate = FFR whenever FFR > 0.25%,
        which holds throughout 2022Q2 onward.
      - Quarterly series = average of monthly values within each quarter.
      - Output first-differenced (ΔMP_t): stationary, interpretable as
        marginal change in policy stance, consistent with other flow controls.

    Returns DataFrame: [quarter_label, wuxia_level, delta_wuxia]
    """
    if WUXIA_CACHE.exists():
        df = pd.read_parquet(WUXIA_CACHE)
        print(f"  Wu-Xia: loaded from cache "
              f"({df['quarter_label'].min()}–{df['quarter_label'].max()})")
        return df

    # ── Step 1: download Atlanta Fed Excel ───────────────────────────────────
    print("  Downloading Wu-Xia shadow FFR from Atlanta Fed ...")
    try:
        resp = requests.get(WUXIA_URL, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        raise RuntimeError(
            f"Wu-Xia download failed: {e}\n"
            f"  URL: {WUXIA_URL}\n"
            f"  Download the Excel manually and save as:\n"
            f"  {WUXIA_CACHE}"
        ) from e

    raw = pd.read_excel(io.BytesIO(resp.content), sheet_name="Data", header=None)
    raw = raw.iloc[1:].reset_index(drop=True)   # skip header row
    raw["date"]     = pd.to_datetime(raw[0])
    raw["fedfunds"] = pd.to_numeric(raw[1], errors="coerce")
    raw["wuxia"]    = pd.to_numeric(raw[2], errors="coerce")

    # Prefer shadow rate; fall back to FEDFUNDS for pre-1990
    raw["rate"] = raw["wuxia"].where(raw["wuxia"].notna(), raw["fedfunds"])
    excel_end   = raw["date"].max()
    raw = raw[["date", "rate"]].dropna().copy()
    print(f"  Excel covers: {raw['date'].min().date()} – {excel_end.date()}")

    # ── Step 2: extend with FEDFUNDS from FRED past Excel end date ───────────
    # Shadow rate = FFR by construction when FFR > 0.25% (post-ZLB period).
    extension_start = excel_end + pd.DateOffset(months=1)
    print(f"  Fetching FEDFUNDS from FRED from {extension_start.date()} onward ...")
    try:
        fred_client = _get_fred_client()
        ff = fred_client.get_series(
            "FEDFUNDS",
            observation_start=extension_start.strftime("%Y-%m-%d"),
        )
        ff = ff.dropna()
        if not ff.empty:
            ext = pd.DataFrame({
                "date": ff.index.to_pydatetime(),
                "rate": ff.values,
            })
            ext["date"] = pd.to_datetime(ext["date"])
            raw = pd.concat([raw, ext], ignore_index=True)
            print(f"  FEDFUNDS extension: {ext['date'].min().date()} – "
                  f"{ext['date'].max().date()}  ({len(ext)} months appended)")
        else:
            print("  [warn] FEDFUNDS extension returned no data — "
                  "series will end at Excel coverage.")
    except Exception as e:
        print(f"  [warn] FEDFUNDS extension failed ({e}). "
              "Series will end at Excel coverage. "
              "Set FRED_API_KEY to enable extension.")

    # ── Step 3: monthly → quarterly average → first difference ───────────────
    raw = raw.sort_values("date").reset_index(drop=True)
    raw["quarter_label"] = raw["date"].apply(_dt_to_ql)
    qtr = (raw.groupby("quarter_label")["rate"]
              .mean()
              .reset_index()
              .sort_values("quarter_label")
              .reset_index(drop=True))
    qtr.rename(columns={"rate": "wuxia_level"}, inplace=True)
    qtr["delta_wuxia"] = qtr["wuxia_level"].diff()

    print(f"  Combined series: {qtr['quarter_label'].min()}–"
          f"{qtr['quarter_label'].max()}, "
          f"{len(qtr)} quarters, "
          f"{qtr['delta_wuxia'].notna().sum()} with ΔMP_t non-null")
    qtr.to_parquet(WUXIA_CACHE, index=False)
    print(f"  Saved: {WUXIA_CACHE}")
    return qtr


# ── v3: Census CBP establishment size ─────────────────────────────────────────

def _fetch_cbp_size() -> pd.DataFrame:
    """
    Compute log average establishment size φ_j by BLS supersector from Census
    County Business Patterns (CBP), year 2000 (pre-sample, predetermined).

    Identification (Gertler-Gilchrist 1994):
      φ_j = log(avg employment per establishment) in supersector j.
      Small establishments are bank-dependent and most sensitive to monetary
      tightening. Large establishments have capital-market access and are less
      affected by C&I lending standards.

      High φ_j → less bank-dependent → weaker response to ΔMP_t.
      Low  φ_j → more bank-dependent → stronger response to ΔMP_t.

    The interaction φ_j × ΔMP_t in the residualization removes the component
    of industry shock variation driven by differential credit supply tightening
    across monetary policy cycles — the leading candidate for why r(δ,LD) = 0.423
    survives v2 residualization.

    Source: US Census Bureau CBP public API (api.census.gov). No key required.

    Returns DataFrame: [industry_code (int), avg_size, log_avg_size]
    """
    if ESTAB_SZ_CACHE.exists():
        df = pd.read_parquet(ESTAB_SZ_CACHE)
        print(f"  CBP size: loaded from cache ({len(df)} supersectors)")
        return df

    print(f"  Fetching Census CBP {CBP_YEAR} establishment size ...")
    CBP_URL = f"https://api.census.gov/data/{CBP_YEAR}/cbp"
    records = []

    for bls_code, naics_list in BLS_TO_CBP_NAICS.items():
        total_emp, total_est = 0, 0
        ok = True
        for naics in naics_list:
            params = {
                "get": f"{CBP_NAICS_VAR},EMP,ESTAB",
                "for": "us:*",
                CBP_NAICS_VAR: naics,
            }
            try:
                r = requests.get(CBP_URL, params=params, timeout=20)
                r.raise_for_status()
                data = r.json()
                if len(data) < 2:
                    raise ValueError(f"Empty response for NAICS {naics}")
                row = data[1]
                total_emp += int(row[1])
                total_est += int(row[2])
            except Exception as e:
                print(f"  [warn] BLS {bls_code}, NAICS {naics}: {e}")
                ok = False
                break

        if not ok or total_est == 0:
            sys.exit(
                f"[ERROR] CBP fetch failed for BLS supersector {bls_code}.\n"
                "Check Census API availability and NAICS mapping."
            )

        avg_size     = total_emp / total_est
        log_avg_size = np.log(avg_size)
        records.append({
            "industry_code":  int(bls_code),
            "avg_size":       avg_size,
            "log_avg_size":   log_avg_size,
        })
        print(f"  BLS {bls_code:>2d}: emp={total_emp:>11,}  "
              f"estabs={total_est:>9,}  "
              f"avg_size={avg_size:>6.1f}  log_size={log_avg_size:.3f}")

    if len(records) < 12:
        sys.exit(
            f"[ERROR] CBP size fetch incomplete: {len(records)}/12 supersectors.\n"
            "Verify Census API and BLS_TO_CBP_NAICS mapping."
        )

    df = pd.DataFrame(records)[["industry_code", "avg_size", "log_avg_size"]]
    df.to_parquet(ESTAB_SZ_CACHE, index=False)
    print(f"  Saved: {ESTAB_SZ_CACHE}")
    return df


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

    # ── (e) monetary policy × establishment size interaction (v3) ─────────────
    print("\n--- Monetary policy × establishment size (v3) ---")

    # Wu-Xia shadow rate: quarterly first difference ΔMP_t
    if WUXIA_CACHE.exists():
        wuxia = pd.read_parquet(WUXIA_CACHE)
        print("  Wu-Xia: loaded from cache")
    else:
        wuxia = _fetch_wuxia()

    wuxia_delta = wuxia[["quarter_label", "delta_wuxia"]].dropna().copy()
    print(f"  ΔMP_t: {wuxia_delta['quarter_label'].min()}–"
          f"{wuxia_delta['quarter_label'].max()}  ({len(wuxia_delta)} qtrs)")

    # CBP establishment size: predetermined φ_j (time-invariant per industry)
    if ESTAB_SZ_CACHE.exists():
        size = pd.read_parquet(ESTAB_SZ_CACHE)
        print("  CBP size: loaded from cache")
    else:
        size = _fetch_cbp_size()

    size["industry_code"] = size["industry_code"].astype(int)
    nat["industry_code"]  = nat["industry_code"].astype(int)

    # Print φ_j values for transparency
    print("  φ_j (log avg establishment size), low→high credit sensitivity:")
    for _, row in size.sort_values("log_avg_size").iterrows():
        lbl = INDUSTRY_LABELS.get(int(row["industry_code"]),
                                   str(int(row["industry_code"])))
        print(f"    BLS {int(row['industry_code']):>2d}  {lbl:<35}  "
              f"φ_j={row['log_avg_size']:.3f}  "
              f"(avg {row['avg_size']:.1f} workers/estab)")

    # Merge φ_j (time-invariant) and ΔMP_t (industry-invariant) into panel
    nat = nat.merge(size[["industry_code", "log_avg_size"]],
                    on="industry_code", how="left")
    nat = nat.merge(wuxia_delta, on="quarter_label", how="left")

    # Construct interaction: φ_j × ΔMP_t (varies by both industry AND quarter)
    nat["phi_x_dmp"] = nat["log_avg_size"] * nat["delta_wuxia"]

    n_miss = nat["phi_x_dmp"].isna().sum()
    if n_miss:
        print(f"  [info] {n_miss} rows missing phi_x_dmp "
              "(first quarter of Wu-Xia series, expected) — dropped.")
    nat = nat.dropna(subset=["phi_x_dmp"])
    print(f"  After v3 merge: {len(nat):,} obs "
          f"({nat['quarter_label'].min()}–{nat['quarter_label'].max()})")

    return nat.sort_values(["industry_code", "quarter_label"]).reset_index(drop=True)


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
nat.head(5)

# ── 2. residualize ────────────────────────────────────────────────────────────
SPEC      = "v3: Δlog(p) + Δlog(VA_lag) + φ_j×ΔMP_t [+ log(θ_lag) for LD/QU]"
reg_base  = ("dlog_p", "dlog_va_lag", "phi_x_dmp")
reg_tight = ("dlog_p", "dlog_va_lag", "log_theta_lag", "phi_x_dmp")

# Full regression equation estimated by _residualize() via within-industry OLS:
#
#   For δ and total separations (TS):
#     log g^k_{j,t} = α_j
#                   + γ   · Δlog p_t           [aggregate productivity growth]
#                   + λ   · Δlog VA_{j,t-1}    [lagged industry VA growth]
#                   + ψ   · φ_j × ΔMP_t        [MP × log avg establishment size]
#                   + ν^k_{j,t}                [residual → new instrument]
#
#   For LD and QU (additionally absorbs aggregate labour-market conditions):
#     log g^k_{j,t} = α_j
#                   + γ   · Δlog p_t
#                   + λ   · Δlog VA_{j,t-1}
#                   + μ   · log θ_{t-1}        [lagged log market tightness V/U]
#                   + ψ   · φ_j × ΔMP_t
#                   + ν^k_{j,t}
#
#   α_j absorbed by within-industry demeaning (pooled OLS on demeaned data).
#   All regressors predetermined: aggregate values or lagged one quarter.
#   φ_j fixed at pre-sample CBP 2000 industry characteristics (time-invariant).
#   ΔMP_t = quarterly first difference of Wu-Xia shadow FFR (stationary).
#   After within-industry demeaning the MP interaction becomes:
#     φ_j × (ΔMP_t − mean_t[ΔMP_t])
#   identifying differential sensitivity to the monetary policy cycle by
#   industry bank-dependence, orthogonal to the industry FE.

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
#  DIAGNOSTICS  (plots — wrapped in function to suppress REPL repr output)
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_diagnostics(nat: pd.DataFrame, spec: str) -> None:
    """
    Produce two diagnostic plots and save to data/results/.
    Wrapping in a function prevents matplotlib artist repr from printing
    in interactive / REPL sessions.

    Plot A: raw log shock rates by industry × quarter (time series grid).
    Plot B: scatter of δ vs each s-type, raw (top) vs residualized (bottom).
    """
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
    fig.suptitle("National industry shock rates by supersector (raw log)",
                 fontsize=11, y=1.005)
    fig.tight_layout()
    out_a = RESULTS_DIR / "comovement_raw_series.png"
    fig.savefig(out_a, dpi=120, bbox_inches="tight")
    plt.close(fig)
    _open_file(out_a)
    print(f"\nPlot A saved: {out_a}")

    # Plot B: δ vs each s-type — raw (top row) vs residualized (bottom row)
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
            ax.set_xlabel(xl, fontsize=9)
            ax.set_ylabel(yl, fontsize=9)
            ax.set_title(
                f"{lbl}  {'raw' if row==0 else 'residualized'}  (r={r_val:.3f})",
                fontsize=9,
            )
            ax.axhline(0, color="black", linewidth=0.5)
            ax.axvline(0, color="black", linewidth=0.5)
            ax.grid(linewidth=0.4, alpha=0.4)
    fig.suptitle(
        rf"$\delta$ vs each s-type: raw (top) vs residualized (bottom) [{spec}]",
        fontsize=11,
    )
    fig.tight_layout()
    out_b = RESULTS_DIR / "comovement_scatter.png"
    fig.savefig(out_b, dpi=120, bbox_inches="tight")
    plt.close(fig)
    _open_file(out_b)
    print(f"Plot B saved: {out_b}")


# ── 5. plots ──────────────────────────────────────────────────────────────────
_plot_diagnostics(nat, SPEC)

print("\nDone.")