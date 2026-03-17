"""
construct_s_instrument.py
=========================
Constructs the match-separation (s) Bartik instrument for the shock-specific
local projection exercise.

The instrument for state s at quarter t is:

    B^s_{s,t} = sum_j  omega_{s,j,t0} * g^s_{-s,j,t}

where:
    omega_{s,j,t0}  = share of supersector j in state s total private
                      nonfarm employment in base year t0  (QCEW, identical
                      to the delta instrument — reuses shares_base{YEAR}.parquet)

    g^s_{-s,j,t}   = leave-one-out quarterly total separation rate at
                      approximately surviving employers in supersector j at
                      quarter t  (JOLTS total separations):

                          separations^nat_{j,t}
                          ---------------------
                           E^nat_{-s,j,t0}

Numerator: JOLTS total separations
-----------------------------------
The JOLTS 'total separations' series (data element TS) captures all
employer-employee separations regardless of initiation: layoffs, discharges,
voluntary quits, and retirements.  This is the appropriate empirical
counterpart to the standard DMP match dissolution rate s, which is agnostic
about initiation — a quit triggers the same vacancy repost as a layoff.

By contrast, JOLTS layoffs and discharges (LD) would exclude quits, which
constitute a large and economically relevant share of total separations and
belong in the model's s shock: when a worker quits at a surviving firm, the
firm immediately reposts the vacancy, generating exactly the same response
as an employer-initiated layoff.

Davis-Faberman-Haltiwanger (DFH) note
---------------------------------------
JOLTS significantly undercounts separations at closing establishments
because those establishments exit the survey sample before responding.
This sampling limitation is advantageous for our purpose: the published
JOLTS total separations series already approximates separations at
*continuing* establishments only — the closing-establishment component
is largely absent.  No correction for closing-establishment contamination
is therefore needed.  The primary residual DFH concern (undercounting at
rapidly contracting continuing establishments) affects levels but not the
cross-industry variation that identifies the Bartik instrument.

LOO correction
--------------
Denominator-only LOO is applied universally:

    E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}   (from QCEW base year)

JOLTS does not publish state-level total separations at the supersector
level via the public API.

Data source
-----------
    National JOLTS total separations (TS): available from 2000M12 at
    supersector level.  First complete quarter: 2001Q1.

Sample start: 2001Q1

JOLTS series ID format  (21 chars, from jt.txt)
-----------------------------------------------
    JT + seasonal(1) + industry(6) + state(2) + area(5) + sizeclass(2)
       + dataelement(2) + ratelevel(1)

Target parameters:
    seasonal     = S   (seasonally adjusted)
    state        = 00  (national)
    area         = 00000
    sizeclass    = 00  (all sizes)
    dataelement  = TS  (total separations)  ← changed from LD (layoffs)
    ratelevel    = L   (level, thousands of jobs)

Example: JT+S+110099+00+00000+00+TS+L  =  "JTS1100990000000TSL"  (21 chars)

JOLTS industry codes (confirmed from jt.industry):
    110099 = Mining and logging
    230000 = Construction
    300000 = Manufacturing
    400000 = Trade, transportation, and utilities  <- direct aggregate
    510000 = Information
    510099 = Financial activities
    540099 = Professional and business services
    600000 = Private education and health services
    700000 = Leisure and hospitality
    810000 = Other services

Monthly to quarterly aggregation
---------------------------------
JOLTS publishes monthly levels (thousands of separations during the month).
We sum the three months in each calendar quarter to obtain a quarterly flow
comparable in units to the annual QCEW base employment denominator.
Incomplete quarters (fewer than 3 months returned) are dropped.

API call count (first run only; results cached as parquet thereafter)
----------------------------------------------------------------------
    10 series, 1 batch, 3 year-chunks = 3 API calls total
    Chunks: 2000-2008 / 2009-2017 / 2018-2030

Requirements
------------
    pip install pandas numpy requests

Typical usage
-------------
    from construct_s_instrument import build_s_instrument

    df = build_s_instrument(base_year=2006)
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests


# ---------------------------------------------------------------------------
# 1.  CONSTANTS
# ---------------------------------------------------------------------------

STATE_FIPS = {
    "AL": "01000", "AK": "02000", "AZ": "04000", "AR": "05000", "CA": "06000",
    "CO": "08000", "CT": "09000", "DE": "10000", "FL": "12000", "GA": "13000",
    "HI": "15000", "ID": "16000", "IL": "17000", "IN": "18000", "IA": "19000",
    "KS": "20000", "KY": "21000", "LA": "22000", "ME": "23000", "MD": "24000",
    "MA": "25000", "MI": "26000", "MN": "27000", "MS": "28000", "MO": "29000",
    "MT": "30000", "NE": "31000", "NV": "32000", "NH": "33000", "NJ": "34000",
    "NM": "35000", "NY": "36000", "NC": "37000", "ND": "38000", "OH": "39000",
    "OK": "40000", "OR": "41000", "PA": "42000", "RI": "44000", "SC": "45000",
    "SD": "46000", "TN": "47000", "TX": "48000", "UT": "49000", "VT": "50000",
    "VA": "51000", "WA": "53000", "WV": "54000", "WI": "55000", "WY": "56000",
}
FIPS_TO_STATE   = {v: k for k, v in STATE_FIPS.items()}
STATE_FIPS_2D   = {abbrev: code[:2] for abbrev, code in STATE_FIPS.items()}
FIPS2D_TO_STATE = {v: k for k, v in STATE_FIPS_2D.items()}

# ---------------------------------------------------------------------------
# Supersector definitions
# ---------------------------------------------------------------------------
INDUSTRY_CODES = ["10", "20", "30", "41", "42", "43", "50", "55", "60", "65", "70", "80"]

INDUSTRY_LABELS = {
    "10": "Mining",
    "20": "Construction",
    "30": "Manufacturing",
    "41": "Wholesale trade",
    "42": "Retail trade",
    "43": "Transport, warehousing & utilities",
    "50": "Information",
    "55": "Financial activities",
    "60": "Professional & business services",
    "65": "Education & health",
    "70": "Leisure & hospitality",
    "80": "Other services",
}

# JOLTS industry codes from jt.industry.
# Unlike BED, JOLTS publishes the three TTU sub-industries as separate series:
#   "41" Wholesale Trade                  → JOLTS code 420000  (NAICS 42)
#   "42" Retail Trade                     → JOLTS code 440000  (NAICS 44-45)
#   "43" Transport, Warehousing & Util    → JOLTS code 480099  (NAICS 48-49+22)
#        (BLS combines 48-49 and 22 into a single published TWU series)
# The aggregate TTU series (400000) is no longer used.
JOLTS_INDUSTRY_MAP = {
    "10": "110099",   # Mining and logging
    "20": "230000",   # Construction
    "30": "300000",   # Manufacturing
    "41": "420000",   # Wholesale trade
    "42": "440000",   # Retail trade
    "43": "480099",   # Transportation, warehousing, and utilities (combined)
    "50": "510000",   # Information
    "55": "510099",   # Financial activities
    "60": "540099",   # Professional and business services
    "65": "600000",   # Private education and health services
    "70": "700000",   # Leisure and hospitality
    "80": "810000",   # Other services
}
JOLTS_TO_INDUSTRY_CODE = {v: k for k, v in JOLTS_INDUSTRY_MAP.items()}

# Thresholds (matching delta instrument)
MIN_SHARE_THRESHOLD     = 0.001
DOMINANT_CELL_THRESHOLD = 0.20

# BLS API
BLS_API_URL = "https://api.bls.gov/publicAPI/v1/timeseries/data/"

# Local cache / output
DEFAULT_CACHE_DIR  = Path("data/cache")
DEFAULT_OUTPUT_DIR = Path("data/instruments")

# Pipeline configuration
# START_QUARTER = 2001Q1: first complete quarter of national JOLTS data.
# (JOLTS launched December 2000; January-March 2001 = first full quarter.)
# END_QUARTER = 2023Q1: aligned with the delta instrument's BDS endpoint.
# Joint estimation sample therefore runs 2001Q1-2023Q1 for both instruments.
BASE_YEAR      = 2006
START_QUARTER  = "2001Q1"
END_QUARTER    = "2023Q1"

SHARES_PATH        = DEFAULT_OUTPUT_DIR / f"shares_base{BASE_YEAR}.parquet"
SHOCK_RATES_S_PATH  = DEFAULT_OUTPUT_DIR / f"shock_rates_s_{START_QUARTER}_{END_QUARTER}.parquet"
SHOCK_RATES_LD_PATH = DEFAULT_OUTPUT_DIR / f"shock_rates_ld_{START_QUARTER}_{END_QUARTER}.parquet"
SHOCK_RATES_QU_PATH = DEFAULT_OUTPUT_DIR / f"shock_rates_qu_{START_QUARTER}_{END_QUARTER}.parquet"
S_INSTRUMENT_PATH  = DEFAULT_OUTPUT_DIR / f"s_instrument_base{BASE_YEAR}.csv"


# ---------------------------------------------------------------------------
# 2.  JOLTS API FETCHING
# ---------------------------------------------------------------------------

# Valid JOLTS dataelement codes for separations decomposition:
#   TS = total separations  (layoffs + discharges + quits + other)
#   LD = layoffs and discharges
#   QU = quits
JOLTS_DATAELEMENTS = {"TS", "LD", "QU"}

# Cache file names keyed by dataelement
JOLTS_CACHE_FILES = {
    "TS": "jolts_separations_national_12ind.parquet",
    "LD": "jolts_ld_national_12ind.parquet",
    "QU": "jolts_quits_national_12ind.parquet",
}

# Shock rate parquet names keyed by dataelement
def _shock_rates_path(dataelement: str, start_q: str, end_q: str,
                      output_dir: Path) -> Path:
    tag = {"TS": "s", "LD": "ld", "QU": "qu"}[dataelement]
    return output_dir / f"shock_rates_{tag}_{start_q}_{end_q}.parquet"


def _build_series_id(industry: str, dataelement: str = "TS") -> str:
    """
    Construct a 21-character national JOLTS series ID.

    JT + S + industry(6) + 00(state) + 00000(area) + 00(size) + DE + L
    Example: "JTS1100990000000TSL"  (total separations)
             "JTS1100990000000LDL"  (layoffs and discharges)
             "JTS1100990000000QUL"  (quits)

    Parameters
    ----------
    industry    : 6-char JOLTS industry code, e.g. "110099"
    dataelement : "TS" (default), "LD", or "QU"
    """
    if dataelement not in JOLTS_DATAELEMENTS:
        raise ValueError(f"dataelement must be one of {JOLTS_DATAELEMENTS}")
    sid = f"JTS{industry}000000000{dataelement}L"
    assert len(sid) == 21, f"Series ID length error: '{sid}' is {len(sid)} chars"
    return sid


def _month_to_quarter(period: str) -> str | None:
    """
    Convert JOLTS period code (e.g. "M03") to quarter string ("1").
    Returns None for annual averages (M13).
    """
    if not period.startswith("M") or period == "M13":
        return None
    return str((int(period[1:]) - 1) // 3 + 1)


def _chunked(lst: list, n: int):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def fetch_jolts_separations_national(
    cache_dir     : Path  = DEFAULT_CACHE_DIR,
    start_quarter : str   = START_QUARTER,
    end_quarter   : str   = END_QUARTER,
    sleep_secs    : float = 1.0,
    dataelement   : str   = "TS",
) -> pd.DataFrame:
    """
    Fetch national JOLTS separation series (SA, level) via the BLS public
    API v1.  Results are cached as parquet after the first call.

    Parameters
    ----------
    cache_dir     : Directory for caching the parquet file.
    start_quarter : First quarter to return, e.g. "2001Q1".
    end_quarter   : Last quarter to return,  e.g. "2024Q2".
    sleep_secs    : Pause between API calls (rate-limit courtesy).
    dataelement   : "TS" total separations (default), "LD" layoffs and
                    discharges, or "QU" quits.

    Returns
    -------
    pd.DataFrame
        Columns: industry_code, quarter_label, separations_nat
        One row per (supersector, quarter).
    """
    if dataelement not in JOLTS_DATAELEMENTS:
        raise ValueError(f"dataelement must be one of {JOLTS_DATAELEMENTS}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / JOLTS_CACHE_FILES[dataelement]

    if cache_file.exists():
        print(f"  [JOLTS {dataelement}] loading from cache")
        df = pd.read_parquet(cache_file)
    else:
        series_ids = [_build_series_id(jt, dataelement)
                      for jt in JOLTS_INDUSTRY_MAP.values()]
        sid_to_ind = {
            _build_series_id(jt, dataelement): pip
            for pip, jt in JOLTS_INDUSTRY_MAP.items()
        }

        print(
            f"  [JOLTS {dataelement}] fetching {len(series_ids)} industry "
            f"separation series from BLS API v1 ..."
        )
        for sid in series_ids:
            print(f"    {sid}  [{INDUSTRY_LABELS[sid_to_ind[sid]]}]")

        # BLS API v1 limit: 10 years per request for unregistered users
        # (20 years with a free registered key).  Use three chunks to stay
        # safely within the limit across the full 2000-present JOLTS history.
        year_chunks = [("2000", "2008"), ("2009", "2017"), ("2018", "2030")]
        batches     = list(_chunked(series_ids, 25))
        n_batches   = len(batches)
        all_obs     = []

        for b_idx, batch in enumerate(batches):
            for c_idx, (start_yr, end_yr) in enumerate(year_chunks):
                call_num    = b_idx * len(year_chunks) + c_idx + 1
                total_calls = n_batches * len(year_chunks)
                print(
                    f"    call {call_num}/{total_calls}  "
                    f"({len(batch)} series, {start_yr}-{end_yr}) ...",
                    end=" ", flush=True,
                )
                resp = requests.post(
                    BLS_API_URL,
                    json    = {"seriesid"  : batch,
                               "startyear" : start_yr,
                               "endyear"   : end_yr},
                    headers = {"Content-Type": "application/json"},
                    timeout = 60,
                )
                resp.raise_for_status()
                result = resp.json()

                if result.get("status") != "REQUEST_SUCCEEDED":
                    raise ValueError(
                        f"BLS API error (call {call_num}): "
                        f"{result.get('message', result)}"
                    )

                n_obs = 0
                for series in result["Results"]["series"]:
                    sid = series["seriesID"]
                    for obs in series["data"]:
                        if obs.get("value", "-") == "-":
                            continue
                        all_obs.append({
                            "industry_code": sid_to_ind[sid],
                            "year"         : obs["year"],
                            "period"       : obs["period"],
                            "value"        : obs["value"],
                        })
                        n_obs += 1
                print(f"{n_obs} obs")

                is_last = (
                    b_idx == n_batches - 1 and
                    c_idx == len(year_chunks) - 1
                )
                if sleep_secs > 0 and not is_last:
                    time.sleep(sleep_secs)

        if not all_obs:
            raise ValueError(
                "BLS API returned no observations for any JOLTS national "
                "layoff series.  Check series IDs and API availability.\n"
                f"  Series tried: {series_ids}"
            )

        # Monthly -> quarterly: sum, keep only complete quarters (3 months)
        raw = pd.DataFrame(all_obs)
        raw["quarter"] = raw["period"].apply(_month_to_quarter)
        raw = raw[raw["quarter"].notna()].copy()
        raw["quarter_label"] = raw["year"] + "Q" + raw["quarter"]
        raw["layoffs"] = pd.to_numeric(
            raw["value"].str.replace(",", ""), errors="coerce"
        )

        quarterly = (
            raw
            .groupby(["industry_code", "quarter_label"])["layoffs"]
            .agg(count="count", total="sum")
            .reset_index()
        )
        quarterly = quarterly[quarterly["count"] == 3].copy()
        df = (
            quarterly
            .rename(columns={"total": "separations_nat"})
            .drop(columns="count")
            .sort_values(["industry_code", "quarter_label"])
            .reset_index(drop=True)
        )

        df.to_parquet(cache_file, index=False)
        print(
            f"  [JOLTS {dataelement}] cached {len(df):,} rows | "
            f"{df['industry_code'].nunique()} supersectors | "
            f"{df['quarter_label'].nunique()} quarters"
        )

    return df[
        (df["quarter_label"] >= start_quarter) &
        (df["quarter_label"] <= end_quarter)
    ].copy()


# ---------------------------------------------------------------------------
# 3.  LOO SHOCK RATES
# ---------------------------------------------------------------------------

def build_loo_shock_rates_s(
    shares      : pd.DataFrame,
    jolts_nat   : pd.DataFrame,
    out_col     : str = "g_s_loo",
) -> pd.DataFrame:
    """
    Compute LOO shock rates for every (state, supersector, quarter) cell.

        g^k_{-s,j,t} = separations_nat_{j,t} / E^nat_{-s,j,t0}

    Works for any separation type (TS, LD, QU) — the caller controls
    which JOLTS series is passed in `jolts_nat` and what the output
    column is named via `out_col`.

    Parameters
    ----------
    shares    : Must contain [state_fips, industry_code, emp_nat_loo].
    jolts_nat : Must contain [industry_code, quarter_label, separations_nat].
    out_col   : Name for the output rate column.  Default "g_s_loo".

    Returns
    -------
    pd.DataFrame
        Columns: state_fips, industry_code, quarter_label, {out_col}
    """
    print(f"\n[Shock rates] LOO applied — output column: {out_col}")

    state_ind = shares[["state_fips", "industry_code"]].drop_duplicates()
    panel = state_ind.merge(jolts_nat, on="industry_code", how="inner")

    panel = panel.merge(
        shares[["state_fips", "industry_code", "emp_nat_loo"]],
        on=["state_fips", "industry_code"],
        how="left",
    )

    panel[out_col] = np.where(
        panel["emp_nat_loo"] > 0,
        panel["separations_nat"] / panel["emp_nat_loo"],
        np.nan,
    )

    n_nan = panel[out_col].isna().sum()
    if n_nan > 0:
        print(f"  WARNING: {n_nan} NaN cells (emp_nat_loo <= 0) -- excluded")

    n_valid = panel[out_col].notna().sum()
    print(
        f"  {n_valid:,} valid cells  |  "
        f"{panel['state_fips'].nunique()} states x "
        f"{panel['industry_code'].nunique()} supersectors x "
        f"{panel['quarter_label'].nunique()} quarters"
    )

    return panel[["state_fips", "industry_code", "quarter_label", out_col]]


# ---------------------------------------------------------------------------
# 4.  BARTIK AGGREGATION
# ---------------------------------------------------------------------------

def build_bartik_instrument_s(
    shares        : pd.DataFrame,
    shock_rates_s : pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate LOO s-shock rates to the state-quarter Bartik instrument:

        B^s_{s,t} = sum_j  omega_{s,j,t0} * g^s_{-s,j,t}

    Parameters
    ----------
    shares        : Employment shares — must contain
                    [state_fips, industry_code, share].
    shock_rates_s : Output of build_loo_shock_rates_s().

    Returns
    -------
    pd.DataFrame
        Columns: state, state_fips, quarter_label,
                 bartik_s, weight_sum, n_supersectors
    """
    merged = (
        shock_rates_s.dropna(subset=["g_s_loo"])
        .merge(
            shares[["state_fips", "industry_code", "share"]],
            on=["state_fips", "industry_code"],
            how="inner",
        )
    )

    def _agg(grp):
        return pd.Series({
            "bartik_s"      : (grp["share"] * grp["g_s_loo"]).sum(),
            "weight_sum"    : grp["share"].sum(),
            "n_supersectors": len(grp),
        })

    instrument = (
        merged
        .groupby(["state_fips", "quarter_label"])
        .apply(_agg)
        .reset_index()
    )

    low = instrument["weight_sum"] < 0.80
    if low.any():
        print(
            f"\n  WARNING: {low.sum()} state-quarter cells have weight "
            f"coverage < 80%."
        )

    instrument["state"] = instrument["state_fips"].map(FIPS2D_TO_STATE)

    return instrument[[
        "state", "state_fips", "quarter_label",
        "bartik_s", "weight_sum", "n_supersectors",
    ]]


# ---------------------------------------------------------------------------
# 5.  TOP-LEVEL FUNCTION
# ---------------------------------------------------------------------------

def build_s_instrument(
    base_year     : int  = BASE_YEAR,
    start_quarter : str  = START_QUARTER,
    end_quarter   : str  = END_QUARTER,
    shares        : pd.DataFrame | None = None,
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    save_output   : bool = True,
    output_dir    : Path = DEFAULT_OUTPUT_DIR,
) -> pd.DataFrame:
    """
    End-to-end construction of the s-shock Bartik instrument.

    Parameters
    ----------
    base_year     : Base year for employment shares.  Default 2006.
    start_quarter : First quarter.  Default "2001Q1".
    end_quarter   : Last quarter.   Default "2024Q2".
    shares        : Pre-loaded shares DataFrame.  If None, loaded from disk.
    cache_dir     : Directory for caching BLS downloads.
    save_output   : Write CSV.  Default True.
    output_dir    : Output directory.

    Returns
    -------
    pd.DataFrame
        Columns: state, state_fips, quarter_label,
                 bartik_s, weight_sum, n_supersectors
    """
    print("=" * 60)
    print(
        f"s instrument  |  base year = {base_year}  |  "
        f"{start_quarter} - {end_quarter}"
    )
    print("=" * 60)

    if shares is None:
        shares_path = output_dir / f"shares_base{base_year}.parquet"
        if not shares_path.exists():
            raise FileNotFoundError(
                f"{shares_path} not found.  Run part1_shares.py first."
            )
        shares = pd.read_parquet(shares_path)
        print(f"\n[Shares] loaded from {shares_path}  ({len(shares):,} rows)")
    else:
        print(f"\n[Shares] {len(shares):,} rows passed directly")

    print(f"\n[JOLTS] fetching national total separations ...")
    jolts_nat = fetch_jolts_separations_national(
        cache_dir, start_quarter, end_quarter
    )
    print(
        f"  {jolts_nat['quarter_label'].nunique()} quarters x "
        f"{jolts_nat['industry_code'].nunique()} supersectors"
    )

    shock_rates = build_loo_shock_rates_s(shares, jolts_nat)

    print("\n[Bartik] aggregating B^s_{s,t} ...")
    instrument = build_bartik_instrument_s(shares, shock_rates)

    print("\n" + "=" * 60)
    print("Instrument summary (cross-state, last 8 quarters):")
    print("=" * 60)
    summary = (
        instrument
        .groupby("quarter_label")["bartik_s"]
        .agg(["mean", "std", "min", "max"])
        .round(6)
    )
    print(summary.tail(8).to_string())

    if save_output:
        output_dir.mkdir(parents=True, exist_ok=True)
        fpath = output_dir / f"s_instrument_base{base_year}.csv"
        instrument.to_csv(fpath, index=False)
        print(f"\nSaved: {fpath}  ({len(instrument):,} rows)")

    return instrument


# ---------------------------------------------------------------------------
# 6.  SHOCK RATES WRAPPER  (called by part2_shock_rates_s.py)
# ---------------------------------------------------------------------------

def build_national_shock_rates_s(
    shares        : pd.DataFrame,
    start_quarter : str  = START_QUARTER,
    end_quarter   : str  = END_QUARTER,
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    save_output   : bool = True,
    output_dir    : Path = DEFAULT_OUTPUT_DIR,
    dataelement   : str  = "TS",
) -> pd.DataFrame:
    """
    Fetch JOLTS national separation series and compute LOO shock rates
    for every (state, supersector, quarter) cell.  Part 2 of the pipeline.

    Parameters
    ----------
    dataelement : "TS" total separations (default), "LD" layoffs and
                  discharges, "QU" quits.

    Returns
    -------
    pd.DataFrame
        Columns: state_fips, industry_code, quarter_label, g_{tag}_loo
        where tag is "s" for TS, "ld" for LD, "qu" for QU.
    """
    tag     = {"TS": "s", "LD": "ld", "QU": "qu"}[dataelement]
    out_col = f"g_{tag}_loo"
    label   = {"TS": "total separations", "LD": "layoffs & discharges",
               "QU": "quits"}[dataelement]

    print(f"\n[JOLTS] fetching national {label} ...")
    jolts_nat = fetch_jolts_separations_national(
        cache_dir, start_quarter, end_quarter, dataelement=dataelement
    )
    print(
        f"  {jolts_nat['quarter_label'].nunique()} quarters x "
        f"{jolts_nat['industry_code'].nunique()} supersectors"
    )

    shock_rates = build_loo_shock_rates_s(shares, jolts_nat, out_col=out_col)

    print(f"\n  Mean LOO {label} rate by supersector (×1 000):")
    mean_by_ss = (
        shock_rates
        .groupby("industry_code")[out_col]
        .mean()
        .mul(1000)
        .round(3)
        .rename(index=INDUSTRY_LABELS)
    )
    print(mean_by_ss.to_string())

    if save_output:
        output_dir.mkdir(parents=True, exist_ok=True)
        fpath = _shock_rates_path(dataelement, start_quarter, end_quarter,
                                  output_dir)
        shock_rates.to_parquet(fpath, index=False)
        print(f"\n  Saved: {fpath}")

    return shock_rates


# Convenience wrappers for LD and QU (used by part2_shock_rates_s.py
# and part2b_shock_comovement.py).
def build_national_shock_rates_ld(shares, start_quarter=START_QUARTER,
                                   end_quarter=END_QUARTER,
                                   cache_dir=DEFAULT_CACHE_DIR,
                                   save_output=True,
                                   output_dir=DEFAULT_OUTPUT_DIR):
    """Layoffs-and-discharges LOO shock rates.  See build_national_shock_rates_s."""
    return build_national_shock_rates_s(
        shares, start_quarter, end_quarter, cache_dir,
        save_output, output_dir, dataelement="LD"
    )


def build_national_shock_rates_qu(shares, start_quarter=START_QUARTER,
                                   end_quarter=END_QUARTER,
                                   cache_dir=DEFAULT_CACHE_DIR,
                                   save_output=True,
                                   output_dir=DEFAULT_OUTPUT_DIR):
    """Quits LOO shock rates.  See build_national_shock_rates_s."""
    return build_national_shock_rates_s(
        shares, start_quarter, end_quarter, cache_dir,
        save_output, output_dir, dataelement="QU"
    )


# ---------------------------------------------------------------------------
# Convenience aliases for each separation component.
# ---------------------------------------------------------------------------
def fetch_jolts_ld_national(**kwargs) -> pd.DataFrame:
    """Fetch JOLTS Layoffs and Discharges (LD). See
    fetch_jolts_separations_national for full parameter docs."""
    return fetch_jolts_separations_national(dataelement="LD", **kwargs)

def fetch_jolts_quits_national(**kwargs) -> pd.DataFrame:
    """Fetch JOLTS Quits (QU). See
    fetch_jolts_separations_national for full parameter docs."""
    return fetch_jolts_separations_national(dataelement="QU", **kwargs)

# Backward-compat: old name preserved so existing imports do not break.
fetch_jolts_layoffs_national = fetch_jolts_separations_national