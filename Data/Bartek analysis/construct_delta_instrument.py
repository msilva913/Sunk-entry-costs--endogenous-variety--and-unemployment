"""
construct_delta_instrument.py
=============================
Constructs the permanence-adjusted product-destruction (delta) Bartik
instrument for the shock-specific local projection exercise.

The instrument for state s at quarter t is:

    B^delta_{s,t} = sum_j  omega_{s,j,t0} * g^delta_{-s,j,t}

where:
    omega_{s,j,t0}   = share of supersector j in state s total private
                       nonfarm employment in base year t0  (QCEW)

    g^delta_{-s,j,t} = leave-one-out national permanent-exit rate from
                       closing establishments in supersector j at
                       quarter t  (BED × BDS permanence calibration):

                           closings^perm_{j,t}
                           -------------------
                             E^nat_{-s,j,t0}

                       where closings^perm_{j,t} = π_{j,y(t)} × BED closings_{j,t}

Permanence calibration
----------------------
Raw BED closings count both permanent exits and temporary shutdowns
(establishments that reopen within 1–3 quarters).  In 2020Q2 this
contamination is severe — approximately 400k of the recorded closings
were lockdown-induced temporary shutdowns that later reopened.

The calibration uses Census BDS annual exit counts by NAICS sector
to compute the supersector-level permanence ratio:

    π_{j,y} = BDS exits_{j,y} / Σ_{q: y(q)=y} BED closings_{j,q}

BDS year y covers March y-1 to March y (BED quarters Q2(y-1)–Q1(y)).
BDS exits are defined identically to BED deaths: establishments absent
from March y payrolls that were present in March y-1.  π_{j,y} is
therefore the fraction of BED closings in supersector j × BDS year y
that proved permanent.

Sample endpoint
---------------
BDS data are available through 2023 (as of early 2026).  BDS year 2023
maps to BED quarters 2022Q2–2023Q1.  The permanence-adjusted delta
instrument runs from 1992Q3 to 2023Q1 by default.

Both data sources are fetched directly from BLS/Census and cached locally —
no manual downloads required.

Industry disaggregation
-----------------------
Uses the 10 BLS supersectors.  BED closings by supersector are fetched
from the BLS API; BDS exits by NAICS sector are fetched from the Census
BDS API and aggregated to supersectors via BDS_SECTOR_TO_INDUSTRY.

Leave-one-out correction
------------------------
Denominator-only LOO is applied universally:

    E^nat_{-s,j,t0} = E^nat_{j,t0} - E_{s,j,t0}

State-level BED closings by supersector are unavailable via the public
API so numerator LOO cannot be implemented.

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
# Supersector / industry definitions
# ---------------------------------------------------------------------------
# INDUSTRY_CODES: the 12 pipeline codes used throughout.  The former
# Trade/transport/utilities aggregate ("40") is split into three:
#
#   "41" = Wholesale trade                      (NAICS 42)
#   "42" = Retail trade                         (NAICS 44-45)
#   "43" = Transportation, warehousing          (NAICS 48-49 + 22)
#           & utilities
#
# Note: pipeline codes "41"/"42"/"43" are internal identifiers only.
# They do NOT correspond to NAICS sector numbers (NAICS 41 does not exist;
# NAICS 42 is Wholesale which maps to pipeline "41", etc.).  The
# INDUSTRY_LABELS dict makes every code human-readable throughout.
#
# QCEW agglvl codes:
#   agglvl=53 — State, by BLS Supersector, by ownership
#               4-digit codes 1011-1027 (BLS-internal, NOT NAICS)
#   agglvl=54 — State, by NAICS Sector, by ownership
#               2-digit NAICS codes as strings: "22","42","44","45","48","49"
#               (combined NAICS sectors like 44-45 and 48-49 appear as two
#               separate rows in the bulk CSV)
#
# For the 9 non-TTU supersectors fetch_qcew_annual uses agglvl=53 (unchanged).
# For the three TTU sub-industries it uses agglvl=54 with NAICS sector codes.
INDUSTRY_CODES = ["10", "20", "30", "41", "42", "43", "50", "55", "60", "65", "70", "80"]

# agglvl=53: BLS supersector codes → pipeline codes (TTU removed)
QCEW_TO_INDUSTRY_CODE = {
    "1011": "10",   # Mining and logging
    "1012": "20",   # Construction
    "1013": "30",   # Manufacturing
    # "1021" Trade/transport/utilities is handled at agglvl=54 below
    "1022": "50",   # Information
    "1023": "55",   # Financial activities
    "1024": "60",   # Professional and business services
    "1025": "65",   # Education and health services
    "1026": "70",   # Leisure and hospitality
    "1027": "80",   # Other services
}
QCEW_SUPERSECTOR_CODES = list(QCEW_TO_INDUSTRY_CODE.keys())  # 9 non-TTU codes

# agglvl=54: NAICS sector industry_code strings → pipeline codes.
# In the QCEW annual bulk CSV at agglvl=54 each 2-digit NAICS sector
# appears as its own row.  Retail (44, 45) and Transport/WH (48, 49)
# are NOT pre-combined: they appear as four separate rows.
# The groupby in fetch_qcew_annual sums 44+45 → "42" and 48+49+22 → "43".
# (The QCEW industry reference catalog shows "44-45" and "48-49" as
# combined labels, but those are display names — not the actual CSV values.)
QCEW_TTU_TO_INDUSTRY_CODE = {
    "42":   "41",   # Wholesale trade
    "44":   "42",   # Retail trade pt 1 (NAICS 44)
    "45":   "42",   # Retail trade pt 2 (NAICS 45)
    "48":   "43",   # Transportation (NAICS 48)
    "49":   "43",   # Warehousing (NAICS 49)
    "22":   "43",   # Utilities (NAICS 22)
}
QCEW_TTU_NAICS_CODES = list(QCEW_TTU_TO_INDUSTRY_CODE.keys())  # 6 codes

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

# BED series ID encodes supersector as a 6-digit industry_code.
# BED separates Goods-producing (1xxxxx) from Service-providing (2xxxxx).
# The 9 non-TTU supersectors each have a direct 1-to-1 BED series.
# The three TTU sub-industries are built by assigning the four existing
# BED sub-series (200010-200040) directly to pipeline codes (see BED_TTU_MAP).
#
# Direct 1-to-1 series (pipeline code → single BED industry code):
BED_INDUSTRY_MAP = {
    "10": "100010",   # Natural resources and mining  (Goods-producing)
    "20": "100020",   # Construction                  (Goods-producing)
    "30": "100030",   # Manufacturing                 (Goods-producing)
    "50": "200050",   # Information
    "55": "200060",   # Financial activities
    "60": "200070",   # Professional and business services
    "65": "200080",   # Education and health services
    "70": "200090",   # Leisure and hospitality
    "80": "200100",   # Other services
}

# TTU sub-component BED series → pipeline codes.
# 200030 (Transportation & WH) and 200040 (Utilities) both map to "43" and
# are summed by the groupby in fetch_bed_closings_national.
BED_TTU_MAP = {
    "200010": "41",   # Wholesale trade
    "200020": "42",   # Retail trade
    "200030": "43",   # Transportation and warehousing → Transport/WH/Util
    "200040": "43",   # Utilities                      → Transport/WH/Util
}

# ---------------------------------------------------------------------------
# BDS sector-to-supersector mapping
# ---------------------------------------------------------------------------
# The Census Bureau's Business Dynamics Statistics (BDS) publishes annual
# establishment exits by 2-digit NAICS sector.  These map to our ten BLS
# supersectors as follows.  Several supersectors require aggregating multiple
# NAICS sectors (Trade/transport/utilities, Financial activities, etc.).
#
# Note on Mining: BLS "Mining and logging" includes NAICS 21 (Mining) plus
# NAICS 113 (Logging, a sub-sector of NAICS 11).  BDS publishes NAICS 11
# (Agriculture, forestry, fishing) as a whole — we use only NAICS 21 as the
# BDS proxy for Mining, since logging is a negligible share and cannot be
# cleanly extracted from NAICS 11.
BDS_SECTOR_TO_INDUSTRY = {
    "11":    "10",   # Agriculture, forestry, fishing & hunting → Natural resources & mining
    "21":    "10",   # Mining, quarrying, oil & gas extraction  → Natural resources & mining
    # BED pipeline "10" = "Natural resources and mining" supersector,
    # which covers both NAICS 11 and NAICS 21.  Both must be included here
    # so that the BDS numerator and BED denominator are measured over the
    # same industry scope.  Omitting NAICS 11 causes π ≈ 0.15 for Mining
    # (BED denominator ~3-4× too large relative to BDS numerator).
    "23":    "20",   # Construction
    "31-33": "30",   # Manufacturing
    "42":    "41",   # Wholesale trade
    "44-45": "42",   # Retail trade
    "48-49": "43",   # Transportation & warehousing → Transport/WH/Util
    "22":    "43",   # Utilities                    → Transport/WH/Util
    "51":    "50",   # Information
    "52":    "55",   # Finance & insurance           → Financial activities
    "53":    "55",   # Real estate & rental          → Financial activities
    "54":    "60",   # Professional & technical svcs → Prof & business svcs
    "55":    "60",   # Management of companies       → Prof & business svcs
    "56":    "60",   # Admin & support               → Prof & business svcs
    "61":    "65",   # Educational services          → Education & health
    "62":    "65",   # Health care & social assistance → Education & health
    "71":    "70",   # Arts, entertainment & recreation → Leisure & hospitality
    "72":    "70",   # Accommodation & food services  → Leisure & hospitality
    "81":    "80",   # Other services
}
BDS_SECTORS   = list(BDS_SECTOR_TO_INDUSTRY.keys())  # 18 NAICS sectors
BDS_API_URL   = "https://api.census.gov/data/timeseries/bds"
# BDS data availability: 1978 through the most recently published year.
# As of early 2026, 2023 is the last complete year.
BDS_END_YEAR  = 2023
BDS_START_YEAR = 1992   # first year we need (aligns with BED start 1992Q3)

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

    BED series ID structure (28 chars, confirmed from bd.txt and bd.industry):
        BD + seasonal(1) + msa(5) + state(2) + county(3) + industry(6)
           + unitanalysis(1) + dataelement(1) + sizeclass(2)
           + dataclass(2) + ratelevel(1) + periodicity(1) + ownership(1)

    Target parameters:
        seasonal     = S  (seasonally adjusted)
        state        = 00 (national)
        msa/county   = 00000 / 000
        industry     = BED industry code from bd.industry (NOT NAICS codes)
                       Goods-producing: 100010, 100020, 100030
                       Service-providing: 200010–200100
        unitanalysis = 1   (establishment)
        dataelement  = 1   (employment level)
        sizeclass    = 00  (all sizes, industry breakdown)
        dataclass    = 06  (Closings)
        ratelevel    = L   (level, not rate)
        periodicity  = Q   (quarterly)
        ownership    = 5   (private sector)

    Trade/transport/utilities (pipeline "40") has no single BED aggregate.
    It is built by summing BED_TRADE_COMPONENTS:
        200010 = Wholesale trade
        200020 = Retail trade
        200030 = Transportation and warehousing
        200040 = Utilities

    Dataclass codes (from bd.txt):
        01 = Gross Job Gains     05 = Contractions
        02 = Expansions          06 = Closings      <-- we want this
        03 = Openings            07 = Births
        04 = Gross Job Losses    08 = Deaths

    Parameters
    ----------
    cache_dir     : Directory for caching the result as parquet.
    start_quarter : First quarter to return, e.g. "1992Q3".
    end_quarter   : Last quarter to return,  e.g. "2024Q2".

    Returns
    -------
    pd.DataFrame
        Columns: quarter_label, industry_code (pipeline supersector code), closings_nat
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    # Cache name reflects the 12-industry disaggregation.
    # Delete bed_closings_national_supersector.parquet (old 10-industry cache)
    # if present — it will be ignored since the filename changed.
    cache_file = cache_dir / "bed_closings_national_12ind.parquet"

    if cache_file.exists():
        print("  [BED national] loading from cache")
        df = pd.read_parquet(cache_file)
    else:
        # All BED codes to request: 9 direct + 4 TTU sub-series = 13 series
        # BED_TTU_MAP maps 200010/200020/200030/200040 directly to pipeline
        # codes "41"/"42"/"43" (200030 and 200040 both → "43", summed below).
        # Format: BDS + 00000(msa) + 00(state) + 000(county)
        #             + {bed_code}(6) + 1(unit) + 1(element)
        #             + 00(sizeclass) + 06(closings) + L(level) + Q + 5(private)
        direct_inv  = {v: k for k, v in BED_INDUSTRY_MAP.items()}
        ttu_codes   = set(BED_TTU_MAP.keys())
        all_bed_codes = list(BED_INDUSTRY_MAP.values()) + list(ttu_codes)
        series_ids = [f"BDS0000000000{c}110006LQ5" for c in all_bed_codes]

        print(f"  [BED national] fetching {len(series_ids)} series "
              f"(9 direct + 4 TTU sub-components) from BLS API v1 ...")
        for c in all_bed_codes:
            pip = direct_inv.get(c) or BED_TTU_MAP.get(c)
            label = INDUSTRY_LABELS.get(pip, c)
            print(f"    BDS0000000000{c}110006LQ5  [{label}]")

        # BLS API v1 limit: 10 years per request for unregistered users.
        # Three chunks cover the full BED history (1992Q3 to present).
        api_url     = "https://api.bls.gov/publicAPI/v1/timeseries/data/"
        year_chunks = [("1992", "2001"), ("2002", "2011"), ("2012", "2030")]

        import time
        raw_rows = []   # one row per (bed_code, quarter)
        n_chunks  = len(year_chunks)
        for c_idx, (start_yr, end_yr) in enumerate(year_chunks):
            print(
                f"    chunk {c_idx+1}/{n_chunks}  "
                f"({len(series_ids)} series, {start_yr}-{end_yr}) ...",
                end=" ", flush=True,
            )
            payload = {
                "seriesid" : series_ids,
                "startyear": start_yr,
                "endyear"  : end_yr,
            }
            resp = requests.post(
                api_url,
                json    = payload,
                headers = {"Content-Type": "application/json"},
                timeout = 60,
            )
            resp.raise_for_status()
            result = resp.json()

            if result.get("status") != "REQUEST_SUCCEEDED":
                raise ValueError(
                    f"BLS API error for {start_yr}-{end_yr}: "
                    f"{result.get('message', result)}"
                )

            n_obs = 0
            for series in result["Results"]["series"]:
                sid      = series["seriesID"]
                bed_code = sid[13:19]
                if bed_code not in direct_inv and bed_code not in ttu_codes:
                    continue
                for obs in series["data"]:
                    if obs.get("value", "-") == "-":
                        continue
                    qnum = obs["period"].lstrip("Q").lstrip("0") or "1"
                    raw_rows.append({
                        "quarter_label": f"{obs['year']}Q{qnum}",
                        "bed_code"     : bed_code,
                        "closings_nat" : float(obs["value"].replace(",", "")),
                    })
                    n_obs += 1
            print(f"{n_obs} obs")
            if n_obs == 0:
                raise ValueError(
                    f"BLS API returned 0 observations for chunk {start_yr}-{end_yr}. "
                    f"This usually means rate limiting. Wait 60 seconds and retry, "
                    f"or register a free BLS API key at https://data.bls.gov/registrationEngine/"
                )

            # Sleep between chunks to avoid silent rate-limit truncation.
            # The BLS API v1 returns REQUEST_SUCCEEDED even when throttled,
            # yielding 0 observations with no error message.
            if c_idx < n_chunks - 1:
                time.sleep(2.0)

        if not raw_rows:
            _probe = requests.post(
                api_url,
                json    = {"seriesid": series_ids[:3],
                           "startyear": "2008", "endyear": "2010"},
                headers = {"Content-Type": "application/json"},
                timeout = 60,
            ).json()
            _sers = "\n".join(
                f"  {s['seriesID']}  →  {len(s.get('data',[]))} obs  "
                f"{s.get('message','')}"
                for s in _probe.get("Results", {}).get("series", [])
            )
            raise ValueError(
                "BLS API returned no data.\n"
                f"  status : {_probe.get('status')}\n"
                f"  message: {_probe.get('message','(none)')}\n"
                f"  probe (2008-2010, first 3 series):\n{_sers}\n"
                f"  series tried: {series_ids}"
            )

        raw = pd.DataFrame(raw_rows)

        # Assign pipeline code: direct series use BED_INDUSTRY_MAP inverse;
        # TTU sub-series use BED_TTU_MAP (200030 and 200040 both → "43").
        raw["industry_code"] = raw["bed_code"].map(
            lambda c: direct_inv.get(c) or BED_TTU_MAP.get(c)
        )
        raw = raw[raw["industry_code"].notna()]

        # Groupby sums the two TWU sub-series (200030+200040) into "43"
        # and keeps all other series as single-row groups unchanged.
        df = (
            raw
            .groupby(["quarter_label", "industry_code"])["closings_nat"]
            .sum()
            .reset_index()
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
# 2b.  BDS FETCHING  (annual exits by NAICS sector, national)
# ---------------------------------------------------------------------------

def _quarter_to_bds_year(quarter_label: str) -> int:
    """
    Map a BED quarter label to the BDS fiscal year that contains it.

    BDS defines an exit as zero employment in March of year y with positive
    employment in March of year y-1.  The four BED quarters that fall within
    BDS year y are:  Q2(y-1), Q3(y-1), Q4(y-1), Q1(y).

    Mapping rule
    ------------
        Q1 of calendar year Y  →  BDS year Y
        Q2, Q3, Q4 of year Y   →  BDS year Y + 1

    Examples
    --------
        "2020Q1"  →  2020    (covered by BDS year 2020: Apr 2019 – Mar 2020)
        "2020Q2"  →  2021    (covered by BDS year 2021: Apr 2020 – Mar 2021)
        "2020Q3"  →  2021
        "2020Q4"  →  2021
    """
    cal_year = int(quarter_label[:4])
    q        = int(quarter_label[5])
    return cal_year if q == 1 else cal_year + 1


def fetch_bds_exits_national(
    cache_dir   : Path  = DEFAULT_CACHE_DIR,
    census_key  : str   = "",
    bds_file    : Path  = None,
) -> pd.DataFrame:
    """
    Load annual establishment exit counts by NAICS sector (national) from
    the Census Bureau Business Dynamics Statistics (BDS).

    Two source modes
    ----------------
    1. Local CSV file (preferred, no API key required):
       Pass bds_file = Path("path/to/bds2023_sec_nat.csv") or equivalent.
       Download once from:
         https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html
       Choose National → Sector table.  Expected columns (tab-separated):
         year, sector, estabs_exit  (plus other columns which are ignored).

    2. Census BDS API (fallback, requires internet):
       Called automatically if bds_file is None or not found.
       Optional census_key raises the rate limit from 500 to unlimited calls/day.
       Register free at https://api.census.gov/data/key_signup.html.

    BDS exit definition: establishment with positive March employment in year
    y-1 and zero employment in March of year y, with no reopening for ≥4
    consecutive quarters — identical non-reopening criterion to BED 'deaths'.

    The BDS column used is `job_destruction_deaths` — employment at permanently
    exiting establishments (raw worker count, converted to thousands to match
    BED API units).  BED `closings_nat` from the BLS API v1 is in thousands of
    workers; BDS reports raw workers.  The /1000 conversion is applied when
    reading the CSV so the π ratio is dimensionally consistent.

    Parameters
    ----------
    cache_dir  : Local directory for the cached parquet file.
    census_key : Optional Census API key (only used in API mode).
    bds_file   : Path to a locally downloaded BDS national-by-sector CSV.
                 If supplied and the file exists, the API is never called.

    Returns
    -------
    pd.DataFrame
        Columns: bds_year (int), industry_code (pipeline supersector key),
                 bds_exits (int)
        One row per (supersector, year).
    """
    import time as _time

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "bds_job_destruction_deaths_supersector.parquet"

    if cache_file.exists():
        print("  [BDS national] loading from cache")
        return pd.read_parquet(cache_file)

    # ------------------------------------------------------------------ #
    # MODE 1: local CSV                                                    #
    # ------------------------------------------------------------------ #
    if bds_file is not None and Path(bds_file).exists():
        print(f"  [BDS national] reading from local file: {bds_file}")
        raw_csv = pd.read_csv(
            bds_file,
            sep    = None,       # auto-detect tab vs comma
            engine = "python",
            dtype  = str,
        )
        # Normalise column names (lowercase, strip whitespace)
        raw_csv.columns = raw_csv.columns.str.strip().str.lower()

        required = {"year", "sector", "job_destruction_deaths"}
        missing_cols = required - set(raw_csv.columns)
        if missing_cols:
            raise ValueError(
                f"BDS CSV is missing expected columns: {missing_cols}\n"
                f"  Found columns: {list(raw_csv.columns)}"
            )

        # job_destruction_deaths = employment at establishments that permanently
        # exited (zero employment for 4+ consecutive quarters after March y).
        # BDS reports this in RAW worker counts.
        # BED closings_nat (from BLS API v1) is in THOUSANDS of workers.
        # Divide BDS by 1000 to align units before computing the π ratio.
        raw_csv = raw_csv[["year", "sector", "job_destruction_deaths"]].copy()
        raw_csv["year"]                  = pd.to_numeric(raw_csv["year"],                  errors="coerce")
        raw_csv["job_destruction_deaths"] = pd.to_numeric(raw_csv["job_destruction_deaths"], errors="coerce")
        raw_csv = raw_csv.dropna(subset=["year", "job_destruction_deaths"])
        raw_csv["year"]                  = raw_csv["year"].astype(int)
        # Convert raw workers → thousands to match BED API units
        raw_csv["job_destruction_deaths"] = raw_csv["job_destruction_deaths"] / 1000.0
        raw_csv["sector"]                 = raw_csv["sector"].str.strip()

        raw_csv = raw_csv[raw_csv["sector"].isin(BDS_SECTORS)].copy()
        if raw_csv.empty:
            raise ValueError(
                f"No matching rows after filtering BDS CSV to pipeline sectors.\n"
                f"  Expected sectors: {BDS_SECTORS}\n"
                f"  Sectors found in file: {raw_csv['sector'].unique().tolist()}"
            )

        raw = raw_csv.rename(columns={
            "year"                  : "bds_year",
            "sector"                : "naics_sector",
            "job_destruction_deaths": "bds_exits_raw",
        })

    # ------------------------------------------------------------------ #
    # MODE 2: Census BDS API                                               #
    # ------------------------------------------------------------------ #
    else:
        if bds_file is not None:
            print(f"  [BDS national] WARNING: bds_file not found at {bds_file} — falling back to API")
        print(
            f"  [BDS national] fetching {len(BDS_SECTORS)} NAICS sectors "
            f"from Census BDS API ..."
        )
        if not census_key:
            print(
                "    NOTE: no Census API key supplied.  Calls are rate-limited "
                "to 500/day.  Register free at https://api.census.gov/data/key_signup.html"
            )

        raw_rows = []
        for i, sector in enumerate(BDS_SECTORS):
            params = {
                "get"   : "YEAR,ESTABS_EXIT",
                "for"   : "us:*",
                "SECTOR": sector,
            }
            if census_key:
                params["key"] = census_key

            try:
                resp = requests.get(BDS_API_URL, params=params, timeout=30)
                resp.raise_for_status()
                data = resp.json()
            except Exception as exc:
                raise RuntimeError(
                    f"BDS API call failed for SECTOR={sector}: {exc}\n"
                    f"  URL tried: {resp.url if 'resp' in dir() else 'N/A'}"
                ) from exc

            if not data or len(data) < 2:
                raise ValueError(
                    f"BDS API returned empty data for SECTOR={sector}.\n"
                    f"  Response: {data}"
                )
            headers = data[0]
            try:
                year_col  = headers.index("YEAR")
                exits_col = headers.index("ESTABS_EXIT")
            except ValueError as exc:
                raise ValueError(
                    f"Unexpected BDS column names for SECTOR={sector}: {headers}"
                ) from exc

            for row in data[1:]:
                try:
                    yr = int(row[year_col])
                    ex = int(row[exits_col]) if row[exits_col] not in (None, "N", "") else None
                except (ValueError, TypeError):
                    continue
                if ex is None:
                    continue
                raw_rows.append({
                    "bds_year"      : yr,
                    "naics_sector"  : sector,
                    "bds_exits_raw" : ex,
                })

            if i < len(BDS_SECTORS) - 1:
                _time.sleep(0.3)

        if not raw_rows:
            raise ValueError(
                "BDS API returned no exit counts across all sectors.\n"
                "  Check Census API availability and sector codes."
            )

        raw = pd.DataFrame(raw_rows)

    # Map NAICS sectors to pipeline supersector codes and aggregate
    raw["industry_code"] = raw["naics_sector"].map(BDS_SECTOR_TO_INDUSTRY)
    raw = raw[raw["industry_code"].notna()].copy()

    df = (
        raw
        .groupby(["bds_year", "industry_code"])["bds_exits_raw"]
        .sum()
        .reset_index()
        .rename(columns={"bds_exits_raw": "bds_exits"})
        .sort_values(["industry_code", "bds_year"])
        .reset_index(drop=True)
    )

    n_years = df["bds_year"].nunique()
    n_ss    = df["industry_code"].nunique()
    print(
        f"  [BDS national] cached {len(df):,} rows | "
        f"{n_ss} supersectors | "
        f"{df['bds_year'].min()}–{df['bds_year'].max()}  ({n_years} years)"
    )
    df.to_parquet(cache_file, index=False)
    return df


def compute_permanence_ratios(
    bed_closings : pd.DataFrame,
    bds_exits    : pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute the quarterly permanence ratio π_{j,y} for each supersector j
    and BDS fiscal year y.

    Definition
    ----------
        π_{j,y} = BDS exits_{j,y}  /  Σ_{q: y(q)=y} BED closings_{j,q}

    where y(q) maps BED quarter q to its BDS year via _quarter_to_bds_year().
    The ratio lies in (0, 1) and represents the fraction of within-year BED
    closings that turned out to be permanent exits rather than temporary
    shutdowns.

    In normal years π ≈ 0.85–0.92 for most supersectors.
    In 2021 (covering BED 2020Q2–2021Q1), π is substantially below 1.0 for
    high-contact service industries (Leisure, Other services) because many
    2020Q2 COVID lockdown closures reopened before March 2021.

    Parameters
    ----------
    bed_closings : Output of fetch_bed_closings_national().
                   Columns: quarter_label, industry_code, closings_nat
    bds_exits    : Output of fetch_bds_exits_national().
                   Columns: bds_year, industry_code, bds_exits

    Returns
    -------
    pd.DataFrame
        Columns: industry_code, bds_year (int), pi (float in (0,1]),
                 bds_exits, bed_closings_sum, n_quarters
        One row per (supersector, BDS year).

    Notes
    -----
    BDS years without a corresponding BDS exit count (either before BDS
    coverage or after its endpoint) receive π = NaN and are excluded from
    the adjustment.  The first BDS year (1993) has only 3 of 4 BED quarters
    available (BED starts 1992Q3, missing 1992Q2) — the ratio is computed on
    the available quarters and is slightly upward-biased for that year only.
    """
    # Map each BED quarter to its BDS year
    bed = bed_closings.copy()
    bed["bds_year"] = bed["quarter_label"].apply(_quarter_to_bds_year)

    # Sum BED closings within each (supersector, BDS year)
    bed_sum = (
        bed.groupby(["industry_code", "bds_year"])
        .agg(bed_closings_sum=("closings_nat", "sum"),
             n_quarters=("closings_nat", "count"))
        .reset_index()
    )

    # Merge with BDS exits
    merged = bed_sum.merge(bds_exits, on=["industry_code", "bds_year"], how="left")

    # Compute ratio; cap at 1.0 to handle rounding artifacts
    merged["pi"] = np.where(
        (merged["bed_closings_sum"] > 0) & merged["bds_exits"].notna(),
        (merged["bds_exits"] / merged["bed_closings_sum"]).clip(upper=1.0),
        np.nan,
    )

    n_nan = merged["pi"].isna().sum()
    if n_nan > 0:
        print(
            f"  [π] {n_nan} (supersector, year) cells with NaN π — "
            f"outside BDS coverage or zero BED closings"
        )

    print(
        f"  [π] computed for {merged['pi'].notna().sum()} "
        f"(supersector, BDS year) cells  |  "
        f"BDS years {int(merged.loc[merged['pi'].notna(), 'bds_year'].min())}–"
        f"{int(merged.loc[merged['pi'].notna(), 'bds_year'].max())}"
    )
    print(
        f"  [π] cross-supersector mean: {merged['pi'].mean():.3f}  "
        f"min: {merged['pi'].min():.3f}  max: {merged['pi'].max():.3f}"
    )

    return merged[[
        "industry_code", "bds_year",
        "pi", "bds_exits", "bed_closings_sum", "n_quarters",
    ]].sort_values(["industry_code", "bds_year"]).reset_index(drop=True)


def apply_permanence_adjustment(
    bed_closings     : pd.DataFrame,
    permanence_ratios: pd.DataFrame,
) -> pd.DataFrame:
    """
    Multiply quarterly BED closings by the permanence ratio π_{j,y(q)} to
    produce a permanent-exit-only quarterly closing series.

    The operation preserves within-year quarterly variation (the relative
    magnitudes of closings across the four quarters of a given BDS year are
    unchanged) while scaling the annual level to confirmed permanent exits.

    Parameters
    ----------
    bed_closings      : Columns: quarter_label, industry_code, closings_nat
    permanence_ratios : Output of compute_permanence_ratios().
                        Columns: industry_code, bds_year, pi

    Returns
    -------
    pd.DataFrame
        Same columns as bed_closings plus:
            bds_year      : BDS fiscal year assigned to each quarter
            pi            : permanence ratio used
            closings_raw  : original unadjusted BED closings (for diagnostics)
        The 'closings_nat' column now contains the permanence-adjusted values.

    Notes
    -----
    Quarters for which π is NaN (no BDS coverage) retain their raw BED
    closing values and a warning is printed.  For the primary estimation
    sample (2001Q1–2023Q1) all quarters have BDS coverage.
    """
    adj = bed_closings.copy()
    adj["bds_year"]     = adj["quarter_label"].apply(_quarter_to_bds_year)
    adj["closings_raw"] = adj["closings_nat"].copy()

    adj = adj.merge(
        permanence_ratios[["industry_code", "bds_year", "pi"]],
        on=["industry_code", "bds_year"],
        how="left",
    )

    n_no_pi = adj["pi"].isna().sum()
    if n_no_pi > 0:
        missing_qtrs = (
            adj.loc[adj["pi"].isna(), "quarter_label"].unique().tolist()[:6]
        )
        print(
            f"  [adj] WARNING: {n_no_pi} quarter-supersector cells have no "
            f"π (no BDS coverage) — raw closings used.  "
            f"Sample quarters: {missing_qtrs}"
        )
        adj["pi"] = adj["pi"].fillna(1.0)   # raw = adjusted when no pi

    adj["closings_nat"] = adj["closings_raw"] * adj["pi"]

    mean_raw = adj["closings_raw"].mean()
    mean_adj = adj["closings_nat"].mean()
    print(
        f"  [adj] mean quarterly closings: raw={mean_raw:,.0f}k  "
        f"adjusted={mean_adj:,.0f}k  "
        f"(ratio={mean_adj/mean_raw:.3f})"
    )

    return adj[[
        "quarter_label", "industry_code",
        "closings_nat", "closings_raw", "bds_year", "pi",
    ]]


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
    # Cache name reflects the 12-industry disaggregation.
    # Delete qcew_{year}_supersector_private.parquet (old 10-industry cache)
    # if present — it will be ignored since the filename changed.
    cache_file = cache_dir / f"qcew_{year}_12ind_private.parquet"

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
    # Two agglvl passes:
    #
    #   agglvl=53 — State, by BLS Supersector, by ownership.
    #               9 non-TTU supersectors (BLS codes 1011-1027 excl. 1021).
    #
    #   agglvl=54 — State, by NAICS Sector, by ownership.
    #               QCEW pre-combines NAICS 44+45 → "44-45" and 48+49 → "48-49"
    #               (same convention as BDS).  Utilities is code "22".
    #               Four industry_code strings: "42", "44-45", "48-49", "22".
    #               "44-45" and "48-49" each map 1-to-1 to pipeline codes,
    #               so no within-pipeline summing is required for Retail or
    #               Transport/Warehousing.  Only "22" and "48-49" share a
    #               target ("43") and are summed by the groupby below.
    valid_fips_5d = set(STATE_FIPS.values())

    # Pass 1: 9 non-TTU supersectors at agglvl=53
    mask53 = (
        (raw["own_code"]      == "5") &
        (raw["area_fips"].isin(valid_fips_5d)) &
        (raw["industry_code"].isin(QCEW_SUPERSECTOR_CODES))
    )
    if "agglvl_code" in raw.columns:
        mask53 &= (raw["agglvl_code"] == "53")

    df53 = raw.loc[mask53, ["area_fips", "industry_code", "annual_avg_emplvl"]].copy()
    df53["industry_code"] = df53["industry_code"].map(QCEW_TO_INDUSTRY_CODE)

    # Pass 2: TTU NAICS sectors at agglvl=54
    mask54 = (
        (raw["own_code"]      == "5") &
        (raw["area_fips"].isin(valid_fips_5d)) &
        (raw["industry_code"].isin(QCEW_TTU_NAICS_CODES))
    )
    if "agglvl_code" in raw.columns:
        mask54 &= (raw["agglvl_code"] == "54")

    df54 = raw.loc[mask54, ["area_fips", "industry_code", "annual_avg_emplvl"]].copy()
    df54["industry_code"] = df54["industry_code"].map(QCEW_TTU_TO_INDUSTRY_CODE)

    df = pd.concat([df53, df54], ignore_index=True)

    if df.empty:
        raise ValueError(
            f"QCEW {year}: both agglvl=53 and agglvl=54 filters returned 0 rows.\n"
            f"  agglvl=53 codes sought: {QCEW_SUPERSECTOR_CODES}\n"
            f"  agglvl=54 codes sought: {QCEW_TTU_NAICS_CODES}\n"
            f"  area_fips sample: {raw['area_fips'].unique()[:8]}\n"
            f"  own_code values:  {raw['own_code'].unique()}\n"
            f"  agglvl sample:    {raw.get('agglvl_code', pd.Series()).unique()[:10]}"
        )

    # Groupby sums: NAICS 44+45 → pipeline "42" (Retail),
    # NAICS 48+49+22 → pipeline "43" (Transport/WH/Util).
    # Wholesale ("41") passes through as a single row.
    df = (
        df.groupby(["area_fips", "industry_code"])["annual_avg_emplvl"]
        .sum().reset_index()
    )

    # Normalise area_fips to 2-digit state FIPS
    df["area_fips"] = df["area_fips"].str[:2]

    df.to_parquet(cache_file, index=False)
    print(f"  [QCEW {year}] cached {len(df):,} rows  "
          f"({df['area_fips'].nunique()} states x "
          f"{df['industry_code'].nunique()} industries)")
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


def shares_to_grid(shares: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot the shares DataFrame into a (states x supersectors) grid.

    Parameters
    ----------
    shares : output of build_employment_shares()

    Returns
    -------
    pd.DataFrame
        Rows   : state abbreviations (AL, AK, ...), sorted alphabetically
        Columns: supersector names (Mining, Construction, ...)
        Values : share  (omega_{s,j} = emp_state_ind / emp_state_total)
    """
    grid = (
        shares
        .pivot(index="state_fips", columns="industry_code", values="share")
        .rename(index=FIPS2D_TO_STATE, columns=INDUSTRY_LABELS)
        .sort_index()
    )
    grid.index.name   = None
    grid.columns.name = None
    return grid


# ---------------------------------------------------------------------------
# 5.  LOO SHOCK RATES
# ---------------------------------------------------------------------------

def build_loo_shock_rates(
    shares  : pd.DataFrame,
    bed_nat : pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute g^delta_{-s,j,t} for every (state, supersector, quarter).

        g^delta_{-s,j,t} = closings^perm_{j,t} / E^nat_{-s,j,t0}

    The numerator is taken from the 'closings_nat' column of bed_nat.
    When bed_nat is the output of apply_permanence_adjustment(), this column
    contains permanence-adjusted values (BED closings × π_{j,y}).  When
    bed_nat is the raw output of fetch_bed_closings_national(), it contains
    unadjusted values — the function is agnostic to which is supplied.

    Denominator LOO is applied universally.  Numerator LOO requires
    state-level BED by supersector (not published via the public API) so is
    not applied here; the denominator correction alone is the standard
    approach in the Bartik literature for this case.

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
    census_key    : str  = "",
    bds_file      : Path = None,
) -> pd.DataFrame:
    """
    End-to-end construction of the delta Bartik instrument using BED Deaths.

    The numerator uses BED Deaths (dataclass=08): establishments absent from
    the UI payroll base for four or more consecutive quarters.  By construction
    Deaths ⊆ Closings, so temporary shutdowns (establishments that reopen
    within 1–3 quarters) are excluded without requiring a BDS permanence ratio
    calibration.  This is the conceptually clean measure of permanent
    product-line destruction — the structural δ shock in the model.

    All data are fetched automatically from BLS and cached locally.

    Parameters
    ----------
    base_year     : Pre-recession base year for employment shares.
                    Default 2006 (Great Recession).  Use 2019 for COVID.
    start_quarter : First quarter of the instrument time series.
                    Default "1992Q3" (earliest BED Deaths date).
    end_quarter   : Last quarter.  Default "2024Q2".
    cache_dir     : Local directory for caching BLS downloads.
    save_output   : Write the instrument to CSV.  Default True.
    output_dir    : Directory for the output CSV.
    census_key    : Unused; retained for call-site compatibility.
    bds_file      : Unused; retained for call-site compatibility.

    Returns
    -------
    pd.DataFrame
        Columns: state, state_fips, quarter_label,
                 bartik_delta, weight_sum, n_supersectors
        One row per (state, quarter).

    Examples
    --------
    # Great Recession identification episode
    df_gr = build_delta_instrument(base_year=2006)

    # COVID episode
    df_covid = build_delta_instrument(base_year=2019, end_quarter="2024Q2")
    """
    print("=" * 60)
    print(
        f"Delta instrument (BED Deaths)  |  "
        f"base year = {base_year}  |  {start_quarter} – {end_quarter}"
    )
    print("=" * 60)

    # Step 1: QCEW employment shares
    shares = build_employment_shares(base_year, cache_dir)

    # Step 2: BED Deaths (dataclass=08) — no BDS permanence calibration needed
    print(f"\n[BED Deaths] fetching national deaths from BLS ...")
    bed_deaths = fetch_bed_deaths_national(
        cache_dir, start_quarter=start_quarter, end_quarter=end_quarter
    )
    print(
        f"  {bed_deaths['quarter_label'].nunique()} quarters  |  "
        f"{bed_deaths['industry_code'].nunique()} supersectors"
    )

    # Rename deaths_nat → closings_nat so build_loo_shock_rates is agnostic
    bed_window = bed_deaths.rename(columns={"deaths_nat": "closings_nat"}).copy()

    # Step 3: LOO shock rates
    shock_rates = build_loo_shock_rates(shares, bed_window)

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

    # Step 8: Save
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
    census_key    : str  = "",
    bds_file      : Path = None,
) -> pd.DataFrame:
    """
    Fetch BED national closings, apply the BDS/BED permanence calibration,
    and compute the leave-one-out shock rate g^delta_{-s,j,t} for every
    (state, supersector, quarter) cell.

    This is Part 2 of the sequential delta-instrument pipeline.

    Permanence calibration
    ----------------------
    Raw BED closings include both permanent exits and temporary shutdowns
    (establishments that reopen within one to three quarters).  The
    calibration multiplies quarterly BED closings by the permanence ratio
    π_{j,y}, computed as:

        π_{j,y} = BDS exits_{j,y} / Σ_{q: y(q)=y} BED closings_{j,q}

    where BDS exits are the Census Bureau's annual establishment exit counts
    by NAICS sector — establishments absent from March y payrolls that were
    present in March y-1.  π_{j,y} is the fraction of BED closings in
    supersector j during BDS fiscal year y that proved to be permanent.

    The adjustment is computed using the full BED history (1992Q3 onward) as
    the denominator for each BDS year, then applied to the requested window.

    Sample endpoint
    ---------------
    BDS data are published through {BDS_END_YEAR} (as of early 2026).
    BDS year {BDS_END_YEAR} maps to BED quarters
    {BDS_END_YEAR-1}Q2 – {BDS_END_YEAR}Q1.  The permanence-adjusted δ
    instrument therefore has a practical endpoint of {BDS_END_YEAR}Q1.
    Quarters beyond that lack a BDS permanence ratio; raw BED closings
    are used as a fallback with a warning.

    Parameters
    ----------
    shares        : Output of build_employment_shares() — Part 1 result.
    start_quarter : First quarter of the shock rate series. Default "1992Q3".
    end_quarter   : Last quarter.  Default "2023Q1" (last BDS-covered quarter).
    cache_dir     : Directory for caching BED and BDS downloads.
    save_output   : Write shock rates to parquet in output_dir.
    output_dir    : Directory for the output parquet.
    census_key    : Optional Census API key for BDS fetch (API mode only).
    bds_file      : Path to locally downloaded BDS national-by-sector CSV.
                    If supplied and the file exists, the Census API is not
                    called.  Download from:
                    https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html

    Returns
    -------
    pd.DataFrame
        Columns: state_fips, industry_code, quarter_label, g_delta_loo
        One row per (state, supersector, quarter).
    """
    # ------------------------------------------------------------------
    # Step 1: Fetch raw BED closings (full history for denominator use)
    # ------------------------------------------------------------------
    print(f"\n[BED] fetching national closings (raw) from BLS ...")
    bed_raw = fetch_bed_closings_national(
        cache_dir, start_quarter="1992Q3", end_quarter="2024Q2"
    )
    print(
        f"  {bed_raw['quarter_label'].nunique()} quarters  |  "
        f"{bed_raw['industry_code'].nunique()} supersectors"
    )

    # ------------------------------------------------------------------
    # Step 2: Fetch BDS annual exits by NAICS sector
    # ------------------------------------------------------------------
    print(f"\n[BDS] fetching annual exits by NAICS sector ...")
    bds_exits = fetch_bds_exits_national(
        cache_dir, census_key=census_key, bds_file=bds_file
    )

    # ------------------------------------------------------------------
    # Step 3: Compute permanence ratios π_{j,y}
    # ------------------------------------------------------------------
    print(f"\n[π] computing permanence ratios ...")
    perm_ratios = compute_permanence_ratios(bed_raw, bds_exits)

    # Cache permanence ratios as a standalone diagnostic file
    perm_path = output_dir / "permanence_ratios_by_supersector.parquet"
    output_dir.mkdir(parents=True, exist_ok=True)
    perm_ratios.to_parquet(perm_path, index=False)
    print(f"  [π] saved: {perm_path}")

    # ------------------------------------------------------------------
    # Step 4: Apply permanence adjustment
    # ------------------------------------------------------------------
    print(f"\n[adj] applying permanence adjustment to BED closings ...")
    bed_adj = apply_permanence_adjustment(bed_raw, perm_ratios)

    # Filter to requested window
    bed_adj_window = bed_adj[
        (bed_adj["quarter_label"] >= start_quarter) &
        (bed_adj["quarter_label"] <= end_quarter)
    ].copy()

    print(
        f"  {bed_adj_window['quarter_label'].nunique()} quarters in window  |  "
        f"{bed_adj_window['industry_code'].nunique()} supersectors"
    )

    # ------------------------------------------------------------------
    # Step 5: Build LOO shock rates on adjusted closings
    # ------------------------------------------------------------------
    print("\n[Shock rates] computing LOO rates g^delta_{{-s,j,t}} "
          "(permanence-adjusted numerator) ...")
    shock_rates = build_loo_shock_rates(shares, bed_adj_window)

    n_total  = len(shock_rates)
    n_valid  = shock_rates["g_delta_loo"].notna().sum()
    n_qtrs   = shock_rates["quarter_label"].nunique()
    print(
        f"  {n_valid:,} / {n_total:,} cells valid  |  "
        f"{shock_rates['state_fips'].nunique()} states × "
        f"{shock_rates['industry_code'].nunique()} supersectors × "
        f"{n_qtrs} quarters"
    )

    # Spot-check: mean adjusted closing rate by supersector
    print("\n  Mean LOO adjusted closing rate by supersector (×1 000):")
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
# END_QUARTER: BED Deaths series available through at least 2024Q2 via BLS API.
# No BDS coverage limit applies — Deaths requires no permanence calibration.
END_QUARTER   = "2024Q2"

SHARES_PATH      = DEFAULT_OUTPUT_DIR / f"shares_base{BASE_YEAR}.parquet"
SHOCK_RATES_PATH = DEFAULT_OUTPUT_DIR / f"shock_rates_{START_QUARTER}_{END_QUARTER}.parquet"
INSTRUMENT_PATH  = DEFAULT_OUTPUT_DIR / f"delta_instrument_base{BASE_YEAR}.csv"


# ---------------------------------------------------------------------------
# 10.  BED DEATHS FETCH AND CLOSINGS-VS-DEATHS DIAGNOSTIC
# ---------------------------------------------------------------------------

def fetch_bed_deaths_national(
    cache_dir     : Path = DEFAULT_CACHE_DIR,
    start_quarter : str  = "1992Q3",
    end_quarter   : str  = "2024Q2",
) -> pd.DataFrame:
    """
    Fetch national BED employment losses from *dying* establishments
    (dataclass=08: absent for 4+ consecutive quarters) by supersector.

    Deaths/Closings <= 1 by construction in every industry x quarter because
    Deaths is a strict subset of Closings.  This series is the conceptually
    clean measure of permanent establishment exits and eliminates the need
    for the BDS/BED permanence ratio calibration.

    Series ID structure: identical to closings but dataclass=08 instead of 06.
        BDS0000000000{bed_code}110008LQ5  (28 chars)
    bed_code at positions [13:19] -- same as fetch_bed_closings_national().
    """
    import time

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "bed_deaths_national_12ind.parquet"

    if cache_file.exists():
        print("  [BED Deaths] loading from cache")
        df = pd.read_parquet(cache_file)
    else:
        direct_inv    = {v: k for k, v in BED_INDUSTRY_MAP.items()}
        ttu_codes     = set(BED_TTU_MAP.keys())
        all_bed_codes = list(BED_INDUSTRY_MAP.values()) + list(ttu_codes)
        # dataclass=08 (Deaths) instead of 06 (Closings)
        series_ids = [f"BDS0000000000{c}110008LQ5" for c in all_bed_codes]

        print(f"  [BED Deaths] fetching {len(series_ids)} series from BLS API v1 ...")

        api_url     = "https://api.bls.gov/publicAPI/v1/timeseries/data/"
        year_chunks = [("1992", "2001"), ("2002", "2011"), ("2012", "2030")]
        raw_rows    = []

        for c_idx, (start_yr, end_yr) in enumerate(year_chunks):
            print(
                f"    chunk {c_idx+1}/{len(year_chunks)}  "
                f"({start_yr}-{end_yr}) ...",
                end=" ", flush=True,
            )
            payload = {
                "seriesid" : series_ids,
                "startyear": start_yr,
                "endyear"  : end_yr,
            }
            resp = requests.post(
                api_url,
                json    = payload,
                headers = {"Content-Type": "application/json"},
                timeout = 60,
            )
            resp.raise_for_status()
            result = resp.json()

            if result.get("status") != "REQUEST_SUCCEEDED":
                raise ValueError(
                    f"BLS API error for {start_yr}-{end_yr}: "
                    f"{result.get('message', result)}"
                )

            n_obs = 0
            for series in result["Results"]["series"]:
                sid      = series["seriesID"]
                bed_code = sid[13:19]    # positions 13-18 inclusive; same as closings
                if bed_code not in direct_inv and bed_code not in ttu_codes:
                    continue
                for obs in series["data"]:
                    if obs.get("value", "-") == "-":
                        continue
                    # Strip leading zero: "Q04" -> quarter "4"
                    qnum = str(int(obs["period"].lstrip("Q")))
                    raw_rows.append({
                        "quarter_label": f"{obs['year']}Q{qnum}",
                        "bed_code"     : bed_code,
                        "deaths_nat"   : float(obs["value"].replace(",", "")),
                    })
                    n_obs += 1
            print(f"{n_obs} obs")
            if c_idx < len(year_chunks) - 1:
                time.sleep(2.0)

        if not raw_rows:
            raise ValueError(
                "BLS API returned no Deaths data.  Check daily quota or register "
                "a free API key at https://data.bls.gov/registrationEngine/"
            )

        raw = pd.DataFrame(raw_rows)
        raw["industry_code"] = raw["bed_code"].map(
            lambda c: direct_inv.get(c) or BED_TTU_MAP.get(c)
        )
        raw = raw[raw["industry_code"].notna()]

        df = (
            raw
            .groupby(["quarter_label", "industry_code"])["deaths_nat"]
            .sum()
            .reset_index()
            .sort_values(["industry_code", "quarter_label"])
            .reset_index(drop=True)
        )
        df.to_parquet(cache_file, index=False)
        print(
            f"  [BED Deaths] cached {len(df):,} rows | "
            f"{df['industry_code'].nunique()} supersectors | "
            f"{df['quarter_label'].nunique()} quarters"
        )

    df = df[
        (df["quarter_label"] >= start_quarter) &
        (df["quarter_label"] <= end_quarter)
    ].copy()
    return df


def plot_closings_vs_deaths(
    output_dir    : Path  = Path("data/results"),
    cache_dir     : Path  = DEFAULT_CACHE_DIR,
    start_quarter : str   = "1992Q3",
    end_quarter   : str   = "2019Q4",
    percentiles   : tuple = (25, 75),
) -> None:
    """
    Two-panel time-series figure comparing BED Closings vs Deaths rates.

    Panel 1 (top): Employment-weighted mean closing rate (blue) and death rate
      (red) across 12 supersectors, with shaded bands for the lower-upper
      industry percentile range.  Units: percent of employment per quarter.

    Panel 2 (bottom): Deaths/Closings ratio -- employment-weighted mean
      (green) with percentile band and a horizontal reference line at 1.0.

    Employment weights: base-year (2006) national industry employment from
      data/instruments/shares_base2006.parquet.

    NBER recession shading: 2001Q1-2001Q4, 2007Q4-2009Q2, 2020Q1-2020Q2.

    Saves to: output_dir / "bed_closings_vs_deaths.png"
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("[plot] loading closings and deaths ...")
    closings = fetch_bed_closings_national(
        cache_dir, start_quarter=start_quarter, end_quarter=end_quarter
    )
    deaths = fetch_bed_deaths_national(
        cache_dir, start_quarter=start_quarter, end_quarter=end_quarter
    )

    # Employment weights: base-year national employment per industry (persons)
    shares_path = DEFAULT_OUTPUT_DIR / "shares_base2006.parquet"
    if shares_path.exists():
        shares = pd.read_parquet(shares_path)
        emp_weights = shares.groupby("industry_code")["emp_nat_ind"].first().to_dict()
    else:
        print("  [plot] shares_base2006.parquet not found; using equal weights")
        emp_weights = {c: 1.0 for c in closings["industry_code"].unique()}

    merged = pd.merge(
        closings, deaths,
        on=["quarter_label", "industry_code"], how="inner"
    )
    # BED levels in thousands of workers; emp_weights from QCEW are in persons.
    # Divide by 1000 so rates are in percent-of-employment units.
    merged["emp_w"]        = merged["industry_code"].map(emp_weights).fillna(1.0) / 1000.0
    merged["closing_rate"] = merged["closings_nat"] / merged["emp_w"] * 100
    merged["death_rate"]   = merged["deaths_nat"]   / merged["emp_w"] * 100
    merged["ratio"]        = merged["deaths_nat"] / merged["closings_nat"].replace(0, np.nan)

    all_quarters = sorted(merged["quarter_label"].unique())
    q_to_idx     = {q: i for i, q in enumerate(all_quarters)}
    merged["q_idx"] = merged["quarter_label"].map(q_to_idx)

    def _wpercentile(values, weights, q):
        arr = np.asarray(values, dtype=float)
        wt  = np.asarray(weights, dtype=float)
        mask = np.isfinite(arr) & np.isfinite(wt) & (wt > 0)
        arr, wt = arr[mask], wt[mask]
        if len(arr) == 0:
            return np.nan
        idx = np.argsort(arr)
        arr, wt = arr[idx], wt[idx]
        cdf = np.cumsum(wt) / wt.sum()
        return float(np.interp(q / 100.0, cdf, arr))

    plo, phi = percentiles
    records = []
    for q_label, grp in merged.groupby("quarter_label"):
        w  = grp["emp_w"].values
        wt = w.sum()
        records.append({
            "q_idx"      : q_to_idx[q_label],
            "close_mean" : (grp["closing_rate"] * w).sum() / wt,
            "death_mean" : (grp["death_rate"]   * w).sum() / wt,
            "ratio_mean" : (grp["ratio"].fillna(0) * w).sum() / wt,
            "close_plo"  : _wpercentile(grp["closing_rate"].values, w, plo),
            "close_phi"  : _wpercentile(grp["closing_rate"].values, w, phi),
            "death_plo"  : _wpercentile(grp["death_rate"].values,   w, plo),
            "death_phi"  : _wpercentile(grp["death_rate"].values,   w, phi),
            "ratio_plo"  : _wpercentile(grp["ratio"].values, w, plo),
            "ratio_phi"  : _wpercentile(grp["ratio"].values, w, phi),
        })

    agg = pd.DataFrame(records).sort_values("q_idx").reset_index(drop=True)
    xs  = agg["q_idx"].values
    xtick_idx    = [i for i, q in enumerate(all_quarters) if q.endswith("Q1")]
    xtick_labels = [all_quarters[i][:4] for i in xtick_idx]

    NBER_RECESSIONS = [
        ("2001Q1", "2001Q4"),
        ("2007Q4", "2009Q2"),
        ("2020Q1", "2020Q2"),
    ]

    def _rec_shade(ax, alpha=0.15):
        for r_start, r_end in NBER_RECESSIONS:
            x0 = q_to_idx.get(r_start)
            x1 = q_to_idx.get(r_end)
            if x0 is not None and x1 is not None:
                ax.axvspan(x0, x1, color="grey", alpha=alpha, zorder=0)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    fig.subplots_adjust(hspace=0.12)

    # Panel 1: Closings and Deaths rates
    ax1 = axes[0]
    _rec_shade(ax1)
    ax1.fill_between(xs, agg["close_plo"], agg["close_phi"],
                     color="#2166ac", alpha=0.15,
                     label=f"Closings {plo}-{phi}th pctile (industry)")
    ax1.fill_between(xs, agg["death_plo"], agg["death_phi"],
                     color="#d6604d", alpha=0.20,
                     label=f"Deaths {plo}-{phi}th pctile (industry)")
    ax1.plot(xs, agg["close_mean"], color="#2166ac", lw=1.8, label="Closings (wtd mean)")
    ax1.plot(xs, agg["death_mean"], color="#d6604d", lw=1.8, label="Deaths (wtd mean)")
    ax1.set_ylabel("Pct of employment per quarter", fontsize=10)
    ax1.set_title(
        "BED Closings vs. Deaths: employment-weighted rates, 12 supersectors",
        fontsize=11,
    )
    ax1.legend(fontsize=8, ncol=2, loc="upper right")
    ax1.tick_params(axis="y", labelsize=9)
    ax1.grid(axis="y", lw=0.4, alpha=0.5)

    # Panel 2: Deaths/Closings ratio
    ax2 = axes[1]
    _rec_shade(ax2)
    ax2.fill_between(xs, agg["ratio_plo"], agg["ratio_phi"],
                     color="#4dac26", alpha=0.20,
                     label=f"Deaths/Closings {plo}-{phi}th pctile")
    ax2.plot(xs, agg["ratio_mean"], color="#4dac26", lw=1.8,
             label="Deaths/Closings (wtd mean)")
    ax2.axhline(1.0, color="black", lw=0.9, ls="--", label="Ratio = 1.0 (upper bound)")
    ax2.set_ylabel("Deaths / Closings", fontsize=10)
    ax2.set_xlabel("Year", fontsize=10)
    ax2.legend(fontsize=8, loc="lower right")
    ax2.tick_params(axis="both", labelsize=9)
    ax2.grid(axis="y", lw=0.4, alpha=0.5)
    ax2.set_ylim(0, 1.15)
    ax2.set_xticks(xtick_idx)
    ax2.set_xticklabels(xtick_labels, rotation=45, ha="right")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "bed_closings_vs_deaths.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] saved: {out_path}")

    print(f"\nDeaths/Closings ratio (wtd mean):")
    print(f"  Full-sample mean : {agg['ratio_mean'].mean():.3f}")
    print(f"  Min              : {agg['ratio_mean'].min():.3f}")
    print(f"  Max              : {agg['ratio_mean'].max():.3f}")
    print(f"  Mean closing rate: {agg['close_mean'].mean():.4f}% of emp/qtr")
    print(f"  Mean deaths  rate: {agg['death_mean'].mean():.4f}% of emp/qtr")


# ---------------------------------------------------------------------------
# 11.  CLI ENTRY POINTS
# ---------------------------------------------------------------------------
# Run from Data/Bartek analysis/:
#
#   python construct_delta_instrument.py              -- plot using cached data
#   python construct_delta_instrument.py --refresh    -- delete Deaths cache and re-fetch
#
# To force a re-fetch without the flag, delete the cache file manually:
#   data/cache/bed_deaths_national_12ind.parquet

if __name__ == "__main__":
    import sys

    refresh = "--refresh" in sys.argv
    if refresh:
        deaths_cache = DEFAULT_CACHE_DIR / "bed_deaths_national_12ind.parquet"
        if deaths_cache.exists():
            deaths_cache.unlink()
            print(f"[cache] deleted {deaths_cache} -- will re-fetch from BLS API")
        else:
            print(f"[cache] {deaths_cache} not found -- will fetch fresh")

    plot_closings_vs_deaths(
        output_dir    = Path("data/results"),
        cache_dir     = DEFAULT_CACHE_DIR,
        start_quarter = START_QUARTER,
        end_quarter   = "2019Q4",
        percentiles   = (25, 75),
    )
