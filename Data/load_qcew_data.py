"""
QCEW Establishment Data Loader
================================

Downloads establishment counts for six industries from the BLS QCEW
annual bulk ZIP files.

ZIP structure (all years use this format)
-----------------------------------------
Each annual ZIP contains one CSV per industry, named:

    {year}.q1-q4.by_industry/{year}.q1-q4 {code} {title}.csv

e.g.:
    2010.q1-q4.by_industry/2010.q1-q4 1021 Trade, transportation, and utilities.csv
    2010.q1-q4.by_industry/2010.q1-q4 31 Food manufacturing.csv

Each CSV covers all four quarters of the year. Columns include
area_fips, own_code, qtr, qtrly_estabs (and others).

We find each needed file by matching the industry code at the start of
the filename (after the date prefix), read only those files, and filter
to area_fips='US000' and own_code=5 (private sector).

Nondurable Manufacturing
------------------------
No QCEW supersector exists for nondurable manufacturing. Constructed by
summing NAICS 31 + NAICS 32 (both have individual files in the ZIP).

Bulk ZIP URL:
    https://data.bls.gov/cew/data/files/{year}/csv/{year}_qtrly_by_industry.zip

Usage
-----
    from load_qcew_data import load_qcew_data
    df = load_qcew_data(start_year=2001, end_year=2023)
    df.to_csv("qcew_establishments_industry.csv")
"""

import io
import re
import zipfile
import time
import requests
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================

BULK_ZIP_URL = (
    "https://data.bls.gov/cew/data/files/{year}/csv/{year}_qtrly_by_industry.zip"
)

AREA_FIPS_NATIONAL = "US000"
OWN_CODE_PRIVATE   = "5"     # string for CSV comparison
REQUEST_DELAY      = 1.0     # seconds between annual ZIP downloads

# Industry codes to extract. Lists are summed to form the BED category.
QCEW_INDUSTRY_CODES = {
    "NDR": ["31", "32"],   # Nondurable Mfg: NAICS 31 + 32 (no supersector)
    "DUR": ["33"],
    "CON": ["1012"],       # Construction supersector (not NAICS 23 in these files)
    "TTU": ["1021"],
    "PBS": ["1024"],
    "LHS": ["1026"],
}

BED_CODE_TO_INDUSTRY = {
    "NDR": "Nondurable Goods Manufacturing",
    "DUR": "Durable Goods Manufacturing",
    "CON": "Construction",
    "TTU": "Trade, Transportation, and Utilities",
    "PBS": "Professional and Business Services",
    "LHS": "Leisure and Hospitality",
}

QTR_END_MONTH = {1: 3, 2: 6, 3: 9, 4: 12}

# =============================================================================
# HELPERS
# =============================================================================

def _quarter_end_date(year: int, qtr: int) -> pd.Timestamp:
    month = QTR_END_MONTH[qtr]
    return pd.Timestamp(year=year, month=month, day=1) + pd.offsets.MonthEnd(0)


def _build_code_to_filename(names: list[str], year: int) -> dict[str, str]:
    """
    Build a mapping {industry_code: zip_entry_name} by parsing ZIP filenames.

    Filename pattern:
        {year}.q1-q4.by_industry/{year}.q1-q4 {CODE} {Title}.csv

    The industry code is the token immediately after the date prefix,
    e.g. '1021' in '2010.q1-q4 1021 Trade, transportation...csv'.
    """
    code_map = {}
    prefix = f"{year}.q1-q4 "   # space-terminated so '31' won't match '310'

    for name in names:
        # Get just the filename part (after the last '/')
        basename = name.split("/")[-1]
        if not basename.lower().endswith(".csv"):
            continue
        if not basename.startswith(prefix):
            continue

        # Extract the code: everything between prefix and the next space
        remainder = basename[len(prefix):]          # e.g. "1021 Trade...csv"
        code = remainder.split(" ")[0].strip()      # e.g. "1021"
        code_map[code] = name

    return code_map


def _read_industry_csv(zf: zipfile.ZipFile, entry: str) -> pd.DataFrame:
    """
    Read one industry CSV from an open ZipFile and return a DataFrame
    filtered to US national private-sector rows.

    Columns returned: qtr (int), qtrly_estabs (int)
    """
    with zf.open(entry) as f:
        df = pd.read_csv(f, dtype=str)

    df.columns = df.columns.str.strip().str.lower()

    required = {"area_fips", "own_code", "qtr", "qtrly_estabs"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"Missing columns {missing} in '{entry}'.\n"
            f"Available: {list(df.columns)}"
        )

    # Filter to national private sector
    mask = (
        (df["area_fips"].str.strip() == AREA_FIPS_NATIONAL) &
        (df["own_code"].str.strip()   == OWN_CODE_PRIVATE)
    )
    df = df.loc[mask, ["qtr", "qtrly_estabs"]].copy()

    df["qtr"] = pd.to_numeric(df["qtr"], errors="coerce").astype("Int64")
    df["qtrly_estabs"] = pd.to_numeric(
        df["qtrly_estabs"].str.replace(",", "", regex=False),
        errors="coerce"
    ).astype("Int64")

    # Keep only quarterly rows (qtr 1-4); some files include annual summary (qtr=5)
    df = df.loc[df["qtr"].between(1, 4)]

    return df.reset_index(drop=True)


# =============================================================================
# MAIN LOADER
# =============================================================================

def load_qcew_data(
    start_year: int = 2001,
    end_year:   int = 2023,
    verbose:    bool = True,
) -> pd.DataFrame:
    """
    Download QCEW quarterly establishment counts for six industries.

    Downloads one ZIP per year. Each ZIP contains one CSV per industry;
    we read only the files for the industry codes we need.

    Parameters
    ----------
    start_year : int   First year (inclusive). NAICS data from 2001 onward.
    end_year   : int   Last year (inclusive).
    verbose    : bool  Print progress messages.

    Returns
    -------
    pd.DataFrame
        Long-format panel indexed by (industry, date).
        Column: establishments (private-sector quarterly reporting units).
    """
    all_records = []
    n_years = end_year - start_year + 1

    for i, year in enumerate(range(start_year, end_year + 1), 1):

        if verbose:
            print(f"[{i}/{n_years}] {year}: downloading...", end=" ", flush=True)

        url = BULK_ZIP_URL.format(year=year)
        response = requests.get(url, timeout=120)

        if response.status_code != 200:
            print(f"ERROR HTTP {response.status_code} — skipping {year}.")
            continue

        zf = zipfile.ZipFile(io.BytesIO(response.content))
        code_map = _build_code_to_filename(zf.namelist(), year)

        if verbose:
            print(f"ZIP has {len(code_map)} industry files. Extracting...", end=" ", flush=True)

        # Read each needed industry code and store by code -> {qtr: estabs}
        code_data: dict[str, dict[int, int]] = {}

        for bed_code, ind_codes in QCEW_INDUSTRY_CODES.items():
            for code in ind_codes:
                if code not in code_map:
                    print(f"\n  WARNING: code '{code}' not in {year} ZIP. "
                          f"Available codes (sample): {list(code_map.keys())[:20]}")
                    code_data[code] = {}
                    continue

                try:
                    df_ind = _read_industry_csv(zf, code_map[code])
                except RuntimeError as err:
                    print(f"\n  WARNING: {err}")
                    code_data[code] = {}
                    continue

                # Store as {qtr: estabs}
                code_data[code] = dict(
                    zip(df_ind["qtr"].astype(int),
                        df_ind["qtrly_estabs"].astype(int))
                )

        # Aggregate into records
        for qtr in range(1, 5):
            date = _quarter_end_date(year, qtr)

            for bed_code, ind_codes in QCEW_INDUSTRY_CODES.items():
                total = 0
                ok = True
                for code in ind_codes:
                    val = code_data.get(code, {}).get(qtr)
                    if val is None:
                        print(f"\n  WARNING: no Q{qtr} value for code '{code}' ({year})")
                        ok = False
                        break
                    total += val

                all_records.append({
                    "date":           date,
                    "industry":       BED_CODE_TO_INDUSTRY[bed_code],
                    "establishments": total if ok else pd.NA,
                })

        if verbose:
            print("done.")

        time.sleep(REQUEST_DELAY)

    df = (
        pd.DataFrame(all_records)
        .assign(date=lambda x: pd.to_datetime(x["date"]))
        .sort_values(["industry", "date"])
        .set_index(["industry", "date"])
    )

    if verbose:
        print(f"\nFinished. {len(df)} records total.")

    return df


# =============================================================================
# SCRIPT ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import os

    OUTPUT_DIR  = r"C:\Users\msilva913\Documents\GitHub\Sunk-entry-costs--endogenous-variety--and-unemployment\Data"
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "qcew_establishments_industry.csv")

    print("=" * 60)
    print("QCEW Establishment Data Download")
    print("=" * 60)
    
    " Test "
    df = load_qcew_data(start_year=2010, end_year=2010, verbose=True)
    
    df = load_qcew_data(start_year=2001, end_year=2023, verbose=True)

    print(f"\nShape: {df.shape}")
    print(f"\nSample:\n{df.head(12)}")
    print(f"\nMissing:\n{df.isnull().sum()}")

    df.to_csv(OUTPUT_FILE)
    print(f"\nSaved -> {OUTPUT_FILE}")
