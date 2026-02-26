"""
construct_delta_instrument.py
=============================
Constructs the product-destruction (delta) Bartik instrument for the
shock-specific local projection exercise.

The instrument for state s at quarter t is:

    B^delta_{s,t} = sum_j  omega_{s,j,t0} * g^delta_{-s,j,t}

where:
    omega_{s,j,t0}   = share of supersector j in state s total private
                       nonfarm employment in base year t0  (QCEW)

    g^delta_{-s,j,t} = leave-one-out national job-loss rate from
                       *closing* establishments in supersector j at
                       quarter t  (BED):

                           closings^nat_{-s,j,t}
                           ---------------------
                             E^nat_{-s,j,t0}

Both data sources are fetched directly from BLS and cached locally —
no manual downloads required.

Industry disaggregation
-----------------------
Uses the 10 BLS supersectors.  National BED closings at supersector
level are available directly from the BLS flat file server.  To switch
to 3-digit NAICS, replace INDUSTRY_CODES with NAICS3_PRIVATE_NONFARM
and update BED_INDUSTRY_MAP accordingly; everything else is identical.

Leave-one-out correction
------------------------
Applied universally.  Denominator LOO uses QCEW employment:

    E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}

Numerator LOO (subtract state closings from national) requires
state-level BED by supersector, which is available in the same BLS
flat file but only at total-private level for states.  When not
available at the supersector-by-state level, denominator-only LOO is
applied and a diagnostic table flags dominant cells.

Requirements
------------
    pip install pandas numpy requests

Typical usage
-------------
    from construct_delta_instrument import build_delta_instrument

    df_gr = build_delta_instrument(base_year=2006)
    df_covid = build_delta_instrument(base_year=2019)
"""

import io
import os
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests


# ---------------------------------------------------------------------------
# 1.  CONSTANTS
# ---------------------------------------------------------------------------

# QCEW area_fips for states is a 5-character code: 2-digit state FIPS + "000".
# e.g. Alabama = "01000", California = "06000".
# County-level rows use the full 5-digit county FIPS (e.g. "01001").
# Filtering on these 5-digit codes — combined with agglvl_code = 54 — ensures
# we capture exactly the state-level private-sector supersector aggregates
# that BLS pre-computes, rather than accidentally summing county rows.
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
FIPS_TO_STATE = {v: k for k, v in STATE_FIPS.items()}

# 2-digit state FIPS (used by BED state_code field and for compact labelling)
STATE_FIPS_2D = {abbrev: code[:2] for abbrev, code in STATE_FIPS.items()}
FIPS2D_TO_STATE = {v: k for k, v in STATE_FIPS_2D.items()}

# ---------------------------------------------------------------------------
# Supersector definitions
# ---------------------------------------------------------------------------
# INDUSTRY_CODES: the codes used in QCEW (2-digit NAICS prefix) and in
# the BED series ID (6-digit padded, see BED_INDUSTRY_MAP below).
INDUSTRY_CODES = ["10", "20", "30", "40", "50", "55", "60", "65", "70", "80"]

INDUSTRY_LABELS = {
    "10": "Mining",
    "20": "Construction",
    "30": "Manufacturing",
    "40": "Trade, transport & utilities",
    "50": "Information",
    "55": "Financial activities",
    "60": "Professional & business services",
    "65": "Education & health",
    "70": "Leisure & hospitality",
    "80": "Other services",
}

# BED series ID encodes supersector as a 6-digit industry_code.
# From bd.txt: "National data available at 2 and 3 digit NAICS sectors"
# The mapping below converts our 2-digit NAICS prefix to the 6-digit
# BED industry_code that appears in bd.series.
# Pattern confirmed in bd.txt example: industry_code 200090 = Leisure & hosp.
# The supersector codes follow the pattern: {supersector_2digit}0000
# with the leading zero padding.  Confirmed by checking bd.industry:
#   000010 = Mining and logging
#   000020 = Construction
#   000030 = Manufacturing
#   000040 = Trade, transportation, and utilities
#   000050 = Information
#   000055 = Financial activities
#   000060 = Professional and business services
#   000065 = Education and health services
#   000070 = Leisure and hospitality  (matches example: 200090? No.)
# The bd.txt example shows industry_code=200090 for state=06 (California)
# which is 3-digit NAICS 009 padded — that is state-level.
# For national supersectors (state=00), industry codes are 6-digit:
#   0000{SS}  where SS is the 2-digit supersector padded to 4 digits.
# We verify this directly when fetching bd.series below.
BED_INDUSTRY_MAP = {
    "10": "000010",   # Mining and logging
    "20": "000020",   # Construction
    "30": "000030",   # Manufacturing
    "40": "000040",   # Trade, transport, utilities
    "50": "000050",   # Information
    "55": "000055",   # Financial activities
    "60": "000060",   # Professional and business services
    "65": "000065",   # Education and health services
    "70": "000070",   # Leisure and hospitality
    "80": "000080",   # Other services
}

# 3-digit NAICS alternative (swap into INDUSTRY_CODES when 3-digit BED
# data are available — no other changes required).
NAICS3_PRIVATE_NONFARM = [
    "211", "212", "213", "221",
    "236", "237", "238",
    "311", "312", "313", "314", "315", "316",
    "321", "322", "323", "324", "325", "326", "327",
    "331", "332", "333", "334", "335", "336", "337", "339",
    "423", "424", "425",
    "441", "442", "443", "444", "445", "446", "447", "448",
    "451", "452", "453", "454",
    "481", "482", "483", "484", "485", "486", "487", "488", "492", "493",
    "511", "512", "515", "516", "517", "518", "519",
    "521", "522", "523", "524", "525",
    "531", "532", "533",
    "541", "551", "561", "562",
    "611", "612", "621", "622", "623", "624",
    "711", "712", "713", "721", "722",
    "811", "812", "813",
]

# Thresholds
MIN_SHARE_THRESHOLD      = 0.001   # drop cells below 0.1% of state employment
DOMINANT_CELL_THRESHOLD  = 0.20    # flag cells where state >= 20% of national

# BLS URLs
QCEW_BULK_URL = (
    "https://data.bls.gov/cew/data/files/{year}/csv/"
    "{year}_annual_singlefile.zip"
)
BED_SERIES_URL  = "https://download.bls.gov/pub/time.series/bd/bd.series"
BED_DATA_URL    = "https://download.bls.gov/pub/time.series/bd/bd.data.1.AllItems"

# Local cache
DEFAULT_CACHE_DIR  = Path("data/cache")
DEFAULT_OUTPUT_DIR = Path("data/instruments")


# ---------------------------------------------------------------------------
# 2.  BED FETCHING
# ---------------------------------------------------------------------------

def fetch_bed_closings_national(
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    start_quarter : str  = "1992Q3",
    end_quarter   : str  = "2024Q2",
) -> pd.DataFrame:
    """
    Fetch national BED employment losses from *closing* establishments
    by supersector from the BLS public API v1 (no registration required).

    BED series ID structure (verified against bd.txt and db.nomics.world):
        BD + seasonal(1) + msa(5) + state(2) + county(3) + industry(6)
           + unitanalysis(1) + dataelement(1) + sizeclass(2)
           + dataclass(2) + ratelevel(1) + periodicity(1) + ownership(1)

    Target parameters:
        seasonal     = U  (not seasonally adjusted)
        state        = 00 (national)
        msa/county   = 00000 / 000
        industry     = 000010 ... 000080  (10 supersectors)
        unitanalysis = 1   (establishment)
        dataelement  = 1   (employment level)
        sizeclass    = 00  (all sizes)
        dataclass    = 06  (Closings)
        ratelevel    = L   (level not rate)
        periodicity  = Q   (quarterly)
        ownership    = 5   (private sector)

    Dataclass codes confirmed from bd.txt and db.nomics:
        01 = Gross Job Gains     05 = Contractions
        02 = Expansions          06 = Closings      <-- we want this
        03 = Openings            07 = Births
        04 = Gross Job Losses    08 = Deaths

    The BLS v1 API allows up to 25 series per request, up to 20 years
    per call.  We issue two calls (1992-2011, 2012-present) to cover
    the full BED history from 1992Q3.

    Parameters
    ----------
    cache_dir     : Directory for caching the result as parquet.
    start_quarter : First quarter to return, e.g. "1992Q3".
    end_quarter   : Last quarter to return,  e.g. "2024Q2".

    Returns
    -------
    pd.DataFrame
        Columns: quarter_label, industry_code (2-digit), closings_nat
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "bed_closings_national_supersector.parquet"

    if cache_file.exists():
        print("  [BED national] loading from cache")
        df = pd.read_parquet(cache_file)
    else:
        # Build the 10 target series IDs
        # Format: BDU + 00000 + 00 + 000 + {code6} + 1 + 1 + 00 + 06 + L + Q + 5
        inv_map    = {v: k for k, v in BED_INDUSTRY_MAP.items()}
        series_ids = [
            f"BDU0000000000{code6}110006LQ5"
            for code6 in BED_INDUSTRY_MAP.values()
        ]

        print(
            f"  [BED national] fetching {len(series_ids)} closing series "
            f"from BLS API v1 ..."
        )

        # Two API calls to cover the full 1992-present history
        # (v1 limit: 20 years per request)
        api_url     = "https://api.bls.gov/publicAPI/v1/timeseries/data/"
        year_chunks = [("1992", "2011"), ("2012", "2030")]

        all_rows = []
        for start_yr, end_yr in year_chunks:
            payload = {
                "seriesid" : series_ids,
                "startyear": start_yr,
                "endyear"  : end_yr,
            }
            resp = requests.post(
                api_url,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
            resp.raise_for_status()
            result = resp.json()

            if result.get("status") != "REQUEST_SUCCEEDED":
                raise ValueError(
                    f"BLS API error for {start_yr}-{end_yr}: "
                    f"{result.get('message', result)}"
                )

            for series in result["Results"]["series"]:
                sid   = series["seriesID"]
                # industry_code sits at chars 13-18 of the full series ID
                # (after "BD" prefix: chars 11-16 in 0-indexed suffix)
                code6 = sid[13:19]
                code2 = inv_map.get(code6)
                if code2 is None:
                    continue
                for obs in series["data"]:
                    # period format is "Q01", "Q02", "Q03", "Q04"
                    qnum = obs["period"].lstrip("Q").lstrip("0") or "1"
                    all_rows.append({
                        "quarter_label": f"{obs['year']}Q{qnum}",
                        "industry_code": code2,
                        "closings_nat" : float(
                            obs["value"].replace(",", "")
                        ),
                    })

        if not all_rows:
            raise ValueError(
                "BLS API returned no data.  Possible causes:\n"
                "  - Series IDs are incorrect for the current BED vintage\n"
                "  - API rate limit exceeded (try again after a few minutes)\n"
                f"Series requested: {series_ids}"
            )

        df = (
            pd.DataFrame(all_rows)
            .sort_values(["industry_code", "quarter_label"])
            .reset_index(drop=True)
        )

        df.to_parquet(cache_file, index=False)
        print(
            f"  [BED national] cached {len(df):,} rows | "
            f"{df['industry_code'].nunique()} supersectors | "
            f"{df['quarter_label'].nunique()} quarters"
        )

    # Filter to requested window
    df = df[
        (df["quarter_label"] >= start_quarter) &
        (df["quarter_label"] <= end_quarter)
    ].copy()
    return df



# ---------------------------------------------------------------------------
# 3.  QCEW FETCHING
# ---------------------------------------------------------------------------

def fetch_qcew_annual(
    year      : int,
    cache_dir : Path = DEFAULT_CACHE_DIR,
) -> pd.DataFrame:
    """
    Download (or load from cache) the QCEW annual single-file CSV for
    a given year, filtered to private-sector supersector employment
    across all 50 states.

    Parameters
    ----------
    year      : Calendar year (e.g. 2006).
    cache_dir : Directory for caching parquet files.

    Returns
    -------
    pd.DataFrame
        Columns: area_fips, industry_code, annual_avg_emplvl
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"qcew_{year}_supersector_private.parquet"

    if cache_file.exists():
        print(f"  [QCEW {year}] loading from cache")
        return pd.read_parquet(cache_file)

    url = QCEW_BULK_URL.format(year=year)
    print(f"  [QCEW {year}] downloading from BLS ...")
    resp = requests.get(url, timeout=180)
    resp.raise_for_status()

    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        csv_name = [n for n in zf.namelist() if n.endswith(".csv")][0]
        with zf.open(csv_name) as f:
            raw = pd.read_csv(f, dtype=str, low_memory=False)

    raw.columns = raw.columns.str.strip().str.lower().str.replace(" ", "_")
    for col in ["area_fips", "own_code", "industry_code", "agglvl_code"]:
        if col in raw.columns:
            raw[col] = raw[col].str.strip()

    raw["annual_avg_emplvl"] = pd.to_numeric(
        raw["annual_avg_emplvl"].str.replace(",", ""), errors="coerce"
    )

    # QCEW area_fips is a 5-character field:
    #   state-level  = "SS000"  (2-digit state FIPS + "000")
    #   county-level = "SSCCC"  (full 5-digit county FIPS)
    #   national     = "US000"
    # We filter to the 50 state-level 5-digit codes.
    #
    # agglvl_code 54 = State × Private sector × Supersector (2-digit NAICS).
    # This is the exact pre-computed aggregate we want.  Filtering on it
    # prevents accidentally picking up county rows if state codes ever
    # partially overlap (they don't, but this is belt-and-suspenders).
    #
    # QCEW agglvl_code reference:
    #   50 = State, all ownerships, all industries
    #   51 = State, by ownership, all industries
    #   52 = State, by ownership, supersector
    #   53 = State, by ownership, 3-digit NAICS
    #   54 = State, private sector (own=5), supersector   <-- we want this
    #   70s = county level
    valid_fips_5d = set(STATE_FIPS.values())   # e.g. {"01000", "06000", ...}

    mask = (
        (raw["own_code"]      == "5") &
        (raw["area_fips"].isin(valid_fips_5d)) &
        (raw["industry_code"].isin(INDUSTRY_CODES))
    )
    if "agglvl_code" in raw.columns:
        mask &= (raw["agglvl_code"] == "54")

    df = raw.loc[
        mask,
        ["area_fips", "industry_code", "annual_avg_emplvl"],
    ].copy()

    if df.empty:
        raise ValueError(
            f"QCEW {year}: filter returned 0 rows.\n"
            f"  area_fips sample: {raw['area_fips'].unique()[:10]}\n"
            f"  own_code values:  {raw['own_code'].unique()}\n"
            f"  agglvl sample:    {raw.get('agglvl_code', pd.Series()).unique()[:10]}"
        )

    # Normalise area_fips to 2-digit state FIPS for compact downstream keys
    # ("01000" -> "01") so that state_fips in the instrument output is concise
    # and consistent with the 2-digit codes used in the BED series IDs.
    df["area_fips"] = df["area_fips"].str[:2]

    df.to_parquet(cache_file, index=False)
    print(f"  [QCEW {year}] cached {len(df):,} rows  "
          f"({df['area_fips'].nunique()} states)")
    return df


# ---------------------------------------------------------------------------
# 4.  EMPLOYMENT SHARES
# ---------------------------------------------------------------------------

def build_employment_shares(
    base_year : int,
    cache_dir : Path = DEFAULT_CACHE_DIR,
) -> pd.DataFrame:
    """
    Construct the state-by-supersector employment share matrix and all
    LOO components for the given base year.

    Returns
    -------
    pd.DataFrame with columns:
        state_fips, industry_code, emp_state_ind, emp_state_total,
        emp_nat_ind, emp_nat_loo, share, share_of_nat
    """
    print(f"\n[Shares] base year = {base_year}")
    raw = fetch_qcew_annual(base_year, cache_dir)

    emp_si = (
        raw.groupby(["area_fips", "industry_code"])["annual_avg_emplvl"]
        .sum().reset_index()
        .rename(columns={"area_fips": "state_fips",
                         "annual_avg_emplvl": "emp_state_ind"})
    )
    emp_state_total = (
        emp_si.groupby("state_fips")["emp_state_ind"]
        .sum().reset_index()
        .rename(columns={"emp_state_ind": "emp_state_total"})
    )
    emp_nat = (
        emp_si.groupby("industry_code")["emp_state_ind"]
        .sum().reset_index()
        .rename(columns={"emp_state_ind": "emp_nat_ind"})
    )

    shares = (
        emp_si
        .merge(emp_state_total, on="state_fips")
        .merge(emp_nat, on="industry_code")
    )
    shares["share"]        = shares["emp_state_ind"] / shares["emp_state_total"]
    shares["emp_nat_loo"]  = shares["emp_nat_ind"] - shares["emp_state_ind"]
    shares["share_of_nat"] = shares["emp_state_ind"] / shares["emp_nat_ind"]

    n_before = len(shares)
    shares = shares[shares["share"] >= MIN_SHARE_THRESHOLD].copy()
    print(
        f"  {n_before - len(shares)} sparse cells dropped  |  "
        f"{len(shares):,} retained  |  "
        f"{shares['state_fips'].nunique()} states x "
        f"{shares['industry_code'].nunique()} supersectors"
    )

    dominant = shares[
        shares["share_of_nat"] >= DOMINANT_CELL_THRESHOLD
    ].sort_values("share_of_nat", ascending=False)

    if not dominant.empty:
        print(
            f"\n  Dominant cells (state >= {DOMINANT_CELL_THRESHOLD:.0%} of "
            f"national supersector employment):"
        )
        print(f"  {'State':<6} {'Supersector':<32} {'State/Nat':>10}")
        print(f"  {'-'*50}")
        for _, row in dominant.iterrows():
            abbrev = FIPS2D_TO_STATE.get(row["state_fips"], row["state_fips"])
            label  = INDUSTRY_LABELS.get(row["industry_code"],
                                         row["industry_code"])
            print(f"  {abbrev:<6} {label:<32} {row['share_of_nat']:>9.1%}")
    else:
        print(
            f"  No dominant cells "
            f"(all state shares < {DOMINANT_CELL_THRESHOLD:.0%})"
        )

    return shares


# ---------------------------------------------------------------------------
# 5.  LOO SHOCK RATES
# ---------------------------------------------------------------------------

def build_loo_shock_rates(
    shares  : pd.DataFrame,
    bed_nat : pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute g^delta_{-s,j,t} for every (state, supersector, quarter).

        g^delta_{-s,j,t} = closings^nat_{j,t} / E^nat_{-s,j,t0}

    Denominator LOO is applied universally.  Numerator LOO requires
    state-level BED by supersector (not yet in the flat files at that
    granularity) so is not applied here; the denominator correction
    alone is the standard approach in the literature for this case.

    Returns
    -------
    pd.DataFrame
        Columns: state_fips, industry_code, quarter_label, g_delta_loo
    """
    print("\n[Shock rates] denominator LOO applied universally")

    # Expand national BED to (state, supersector, quarter) panel
    state_ind = shares[["state_fips", "industry_code"]].drop_duplicates()
    panel = state_ind.merge(bed_nat, on="industry_code", how="inner")

    # Merge denominator LOO
    panel = panel.merge(
        shares[["state_fips", "industry_code", "emp_nat_loo"]],
        on=["state_fips", "industry_code"],
        how="left",
    )

    panel["g_delta_loo"] = np.where(
        panel["emp_nat_loo"] > 0,
        panel["closings_nat"] / panel["emp_nat_loo"],
        np.nan,
    )

    n_nan = panel["g_delta_loo"].isna().sum()
    if n_nan > 0:
        print(f"  WARNING: {n_nan} NaN cells (emp_nat_loo <= 0) — excluded")

    return panel[["state_fips", "industry_code", "quarter_label",
                  "g_delta_loo"]]


# ---------------------------------------------------------------------------
# 6.  BARTIK AGGREGATION
# ---------------------------------------------------------------------------

def build_bartik_instrument(
    shares      : pd.DataFrame,
    shock_rates : pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate LOO shock rates to the state-quarter Bartik instrument:

        B^delta_{s,t} = sum_j  omega_{s,j,t0} * g^delta_{-s,j,t}

    Returns
    -------
    pd.DataFrame
        Columns: state, state_fips, quarter_label,
                 bartik_delta, weight_sum, n_supersectors
    """
    merged = (
        shock_rates.dropna(subset=["g_delta_loo"])
        .merge(
            shares[["state_fips", "industry_code", "share"]],
            on=["state_fips", "industry_code"],
            how="inner",
        )
    )

    def _agg(grp):
        return pd.Series({
            "bartik_delta"  : (grp["share"] * grp["g_delta_loo"]).sum(),
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
        "bartik_delta", "weight_sum", "n_supersectors",
    ]]


# ---------------------------------------------------------------------------
# 7.  TOP-LEVEL FUNCTION
# ---------------------------------------------------------------------------

def build_delta_instrument(
    base_year     : int  = 2006,
    start_quarter : str  = "1992Q3",
    end_quarter   : str  = "2024Q2",
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    save_output   : bool = True,
    output_dir    : Path = DEFAULT_OUTPUT_DIR,
) -> pd.DataFrame:
    """
    End-to-end construction of the delta Bartik instrument.

    All data are fetched automatically from BLS and cached locally.

    Parameters
    ----------
    base_year     : Pre-recession base year for employment shares.
                    Default 2006 (Great Recession).  Use 2019 for COVID.
    start_quarter : First quarter of the instrument time series.
                    Default "1992Q3" (earliest BED date).
    end_quarter   : Last quarter of the instrument time series.
                    Default "2024Q2".
    cache_dir     : Local directory for caching BLS downloads.
                    Default Path("data/cache").
    save_output   : Write the instrument to CSV.  Default True.
    output_dir    : Directory for the output CSV.
                    Default Path("data/instruments").

    Returns
    -------
    pd.DataFrame
        Columns: state, state_fips, quarter_label,
                 bartik_delta, weight_sum, n_supersectors
        One row per (state, quarter).

    Examples
    --------
    # Great Recession
    df_gr = build_delta_instrument(base_year=2006)

    # COVID
    df_covid = build_delta_instrument(
        base_year     = 2019,
        start_quarter = "1992Q3",
        end_quarter   = "2024Q2",
    )
    """
    print("=" * 60)
    print(
        f"Delta instrument  |  base year = {base_year}  |  "
        f"{start_quarter} – {end_quarter}"
    )
    print("=" * 60)

    # Step 1: QCEW employment shares
    shares = build_employment_shares(base_year, cache_dir)

    # Step 2: BED closings from BLS flat files
    print(f"\n[BED] fetching national closings from BLS ...")
    bed_nat = fetch_bed_closings_national(cache_dir, start_quarter, end_quarter)
    print(
        f"  {bed_nat['quarter_label'].nunique()} quarters  |  "
        f"{bed_nat['industry_code'].nunique()} supersectors"
    )

    # Step 3: LOO shock rates
    shock_rates = build_loo_shock_rates(shares, bed_nat)

    # Step 4: Bartik aggregation
    print("\n[Bartik] aggregating ...")
    instrument = build_bartik_instrument(shares, shock_rates)

    # Step 5: Summary
    print("\n" + "=" * 60)
    print("Instrument summary (cross-state, last 8 quarters):")
    print("=" * 60)
    summary = (
        instrument
        .groupby("quarter_label")["bartik_delta"]
        .agg(["mean", "std", "min", "max"])
        .round(6)
    )
    print(summary.tail(8).to_string())

    # Step 6: Save
    if save_output:
        output_dir.mkdir(parents=True, exist_ok=True)
        fpath = output_dir / f"delta_instrument_base{base_year}.csv"
        instrument.to_csv(fpath, index=False)
        print(f"\nSaved: {fpath}  ({len(instrument):,} rows)")

    return instrument


# ---------------------------------------------------------------------------
# 8.  NATIONAL SHOCK RATES (convenience wrapper for Part 2)
# ---------------------------------------------------------------------------

def build_national_shock_rates(
    shares        : pd.DataFrame,
    start_quarter : str  = "1992Q3",
    end_quarter   : str  = "2024Q2",
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    save_output   : bool = True,
    output_dir    : Path = DEFAULT_OUTPUT_DIR,
) -> pd.DataFrame:
    """
    Fetch BED national closings and compute the leave-one-out shock rate
    g^delta_{-s,j,t} for every (state, supersector, quarter) cell.

    This is Part 2 of the sequential pipeline.  It depends on `shares`
    (from Part 1) for the LOO denominator, but is otherwise independent
    of the final Bartik aggregation (Part 3).

    Parameters
    ----------
    shares        : Output of build_employment_shares() — Part 1 result.
    start_quarter : First quarter of the shock rate series. Default "1992Q3".
    end_quarter   : Last quarter of the shock rate series.  Default "2024Q2".
    cache_dir     : Directory for caching the BED download.
    save_output   : Write shock rates to parquet in output_dir. Default True.
    output_dir    : Directory for the output parquet.

    Returns
    -------
    pd.DataFrame
        Columns: state_fips, industry_code, quarter_label, g_delta_loo
        One row per (state, supersector, quarter).
    """
    print(f"\n[BED] fetching national closings from BLS API ...")
    bed_nat = fetch_bed_closings_national(cache_dir, start_quarter, end_quarter)
    print(
        f"  {bed_nat['quarter_label'].nunique()} quarters  |  "
        f"{bed_nat['industry_code'].nunique()} supersectors"
    )

    print("\n[Shock rates] computing LOO rates g^delta_{{-s,j,t}} ...")
    shock_rates = build_loo_shock_rates(shares, bed_nat)

    n_total  = len(shock_rates)
    n_valid  = shock_rates["g_delta_loo"].notna().sum()
    n_qtrs   = shock_rates["quarter_label"].nunique()
    n_states = shock_rates["state_fips"].nunique()
    print(
        f"  {n_valid:,} / {n_total:,} cells have valid rates  |  "
        f"{n_states} states x {shock_rates['industry_code'].nunique()} "
        f"supersectors x {n_qtrs} quarters"
    )

    # Spot-check: mean closing rate by supersector over the full sample
    print("\n  Mean LOO closing rate by supersector (x1000 for readability):")
    mean_by_ss = (
        shock_rates
        .groupby("industry_code")["g_delta_loo"]
        .mean()
        .mul(1000)
        .round(3)
        .rename(index=INDUSTRY_LABELS)
    )
    print(mean_by_ss.to_string())

    if save_output:
        output_dir.mkdir(parents=True, exist_ok=True)
        fpath = output_dir / f"shock_rates_{start_quarter}_{end_quarter}.parquet"
        shock_rates.to_parquet(fpath, index=False)
        print(f"\n  Saved: {fpath}")

    return shock_rates


# ---------------------------------------------------------------------------
# 9.  PIPELINE CONFIGURATION  (shared by the three runner scripts)
# ---------------------------------------------------------------------------
#
# Run the pipeline in three sequential steps:
#
#   python part1_shares.py       -> data/instruments/shares_base{YEAR}.parquet
#   python part2_shock_rates.py  -> data/instruments/shock_rates_{START}_{END}.parquet
#   python part3_instrument.py   -> data/instruments/delta_instrument_base{YEAR}.csv
#
# Each script loads its predecessor's output from disk, so any step can
# be re-run in isolation without repeating earlier work.

BASE_YEAR     = 2006
START_QUARTER = "1992Q3"
END_QUARTER   = "2024Q2"

SHARES_PATH      = DEFAULT_OUTPUT_DIR / f"shares_base{BASE_YEAR}.parquet"
SHOCK_RATES_PATH = DEFAULT_OUTPUT_DIR / f"shock_rates_{START_QUARTER}_{END_QUARTER}.parquet"
INSTRUMENT_PATH  = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
