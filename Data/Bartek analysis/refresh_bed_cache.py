# -*- coding: utf-8 -*-
"""
refresh_bed_cache.py
====================
Fetch BED gross job flow data from the BLS public API and write (or overwrite)
the cache file used by part11_jf_table2.py.

    data/cache/bed_gross_flows_correct.parquet

Run from:  Data/Bartek analysis/
Requires:  requests, pandas, pyarrow

WHY THIS SCRIPT EXISTS
----------------------
The BLS API caps unregistered users at 25 requests/day and 10 years/request.
The sandbox environment exhausts the daily quota quickly.  Run this locally
once; the parquet is cached and part11 will not re-fetch on subsequent runs.

SERIES STRUCTURE
----------------
Format: BDS0000000000{ind6}11{dc}LQ5
  ind6 : 6-digit BED industry code (see BED_INDUSTRY_MAP / BED_TTU_MAP)
  dc   : 2-digit dataclass code
    01 = Gross Job Gains       (G,  "expanding" in legacy cache)
    02 = Expansions            (G - G_O, "openings" in legacy cache -- NB misleading name)
    04 = Gross Job Losses      (L,  "contracting")
    06 = Closings / Deaths     (L_C, "closings")

BED data availability: 1992Q3 through latest quarter (2024Q4 as of May 2026).
API rate limits (v1, no key): 25 series/request, 10 years/request, 25 req/day.
Three 10-year windows cover 1992-2024: (1992,2001), (2002,2011), (2012,2024).
Batches of 25 series per call.
"""

import time
import requests
import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CACHE_PATH = Path("data/cache/bed_gross_flows_correct.parquet")
API_URL    = "https://api.bls.gov/publicAPI/v1/timeseries/data/"

BED_INDUSTRY_MAP = {
    "10": "100010",  # Natural resources & mining
    "20": "100020",  # Construction
    "30": "100030",  # Manufacturing
    "50": "200050",  # Information
    "55": "200060",  # Financial activities
    "60": "200070",  # Professional & business services
    "65": "200080",  # Education & health services
    "70": "200090",  # Leisure & hospitality
    "80": "200100",  # Other services
}
BED_TTU_MAP = {
    "200010": "41",  # Wholesale trade
    "200020": "42",  # Retail trade
    "200030": "43",  # Transportation & warehousing  \_ both sum into pip 43
    "200040": "43",  # Utilities                     /
}
ALL_BED_CODES = list(BED_INDUSTRY_MAP.values()) + list(BED_TTU_MAP.keys())
# 9 direct + 4 TTU sub-components = 13 BED codes

# Dataclass codes and their column names in the output parquet
# NOTE: legacy column names are kept for backward compatibility with
#       part11_jf_table2.py's load_bed_panel() remapping.
#   "expanding"   = elem 01 = total gross GAINS   (NOT just expanding estabs)
#   "openings"    = elem 02 = gains at EXPANDING (continuing) establishments
#   "contracting" = elem 04 = total gross LOSSES
#   "closings"    = elem 06 = losses at CLOSING establishments
DATACLS = {"01": "expanding", "02": "openings", "04": "contracting", "06": "closings"}

# Three 10-year windows to satisfy the BLS API year-range limit
WINDOWS = [("1992", "2001"), ("2002", "2011"), ("2012", "2024")]

BATCH_SIZE = 25      # series per API call (v1 unregistered limit)
SLEEP_SEC  = 0.8    # polite pause between calls

# ---------------------------------------------------------------------------
# Build series ID list
# ---------------------------------------------------------------------------
series_meta = []
for dc, col in DATACLS.items():
    for bed_code in ALL_BED_CODES:
        sid = f"BDS0000000000{bed_code}11{dc:0>2}LQ5"
        pip = ({v: k for k, v in BED_INDUSTRY_MAP.items()}.get(bed_code)
               or BED_TTU_MAP.get(bed_code))
        series_meta.append({"sid": sid, "bed_code": bed_code,
                             "pip": pip, "col": col})

sid_to_meta = {m["sid"]: m for m in series_meta}
all_sids    = [m["sid"] for m in series_meta]

# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------

def fetch_window(sids: list, start_yr: str, end_yr: str) -> list:
    """Fetch one time window for a batch of series; return raw observation rows."""
    rows = []
    for i in range(0, len(sids), BATCH_SIZE):
        batch   = sids[i : i + BATCH_SIZE]
        payload = {"seriesid": batch, "startyear": start_yr, "endyear": end_yr}
        resp    = requests.post(
            API_URL,
            json    = payload,
            headers = {"Content-Type": "application/json"},
            timeout = 60,
        )
        resp.raise_for_status()
        result = resp.json()
        if result.get("status") != "REQUEST_SUCCEEDED":
            msgs = result.get("message", [])
            print(f"  WARNING ({start_yr}-{end_yr}): {msgs}")
        for series in result.get("Results", {}).get("series", []):
            sid  = series["seriesID"]
            meta = sid_to_meta.get(sid)
            if not meta:
                continue
            for obs in series["data"]:
                if not obs["period"].startswith("Q"):
                    continue
                qnum   = obs["period"].replace("Q0", "").replace("Q", "")
                qlabel = f"{obs['year']}Q{qnum}"
                try:
                    val = float(obs["value"].replace(",", ""))
                except (ValueError, AttributeError):
                    val = np.nan
                rows.append({
                    "pip": meta["pip"], "col": meta["col"],
                    "quarter_label": qlabel, "value": val,
                })
        time.sleep(SLEEP_SEC)
    return rows


def build_panel(raw_rows: list) -> pd.DataFrame:
    """Pivot raw rows to wide format, merge T&U components, add TOTAL."""
    df = pd.DataFrame(raw_rows).drop_duplicates(["pip", "col", "quarter_label"])

    # Pivot: index=(pip, quarter_label), columns=col
    # aggfunc='sum' correctly accumulates 43a and 43b into pip '43'
    wide = (df.pivot_table(
                index   = ["pip", "quarter_label"],
                columns = "col",
                values  = "value",
                aggfunc = "sum",
            )
            .reset_index())
    wide.columns.name = None

    # Add TOTAL row: sum of 12 supersectors (pip != "00")
    total = (wide.groupby("quarter_label")
                 [["expanding", "openings", "contracting", "closings"]]
                 .sum()
                 .reset_index())
    total["pip"] = "TOTAL"
    out = pd.concat([wide, total], ignore_index=True)
    out = out.rename(columns={"pip": "pip_code"})

    # Derived convenience columns (not strictly needed by part11 but kept for
    # backward compatibility with any other scripts reading the cache)
    out["year"]         = out["quarter_label"].str[:4].astype(int)
    out["gains_total"]  = out["expanding"]
    out["losses_total"] = out["contracting"]
    out["frac_open"]    = (out["expanding"] - out["openings"]) / out["expanding"]
    out["frac_close"]   = out["closings"] / out["contracting"]

    return out.sort_values(["pip_code", "quarter_label"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)

    print(f"Fetching {len(all_sids)} BED series across {len(WINDOWS)} time windows...")
    all_rows = []
    for start_yr, end_yr in WINDOWS:
        print(f"  Window {start_yr}-{end_yr}...", end=" ", flush=True)
        rows = fetch_window(all_sids, start_yr, end_yr)
        all_rows.extend(rows)
        print(f"{len(rows)} observations")

    print(f"Total raw observations: {len(all_rows)}")
    panel = build_panel(all_rows)

    q_min = panel["quarter_label"].min()
    q_max = panel["quarter_label"].max()
    print(f"Quarter range: {q_min} - {q_max}")
    print(f"pip_codes: {sorted(panel['pip_code'].unique())}")

    # Spot-checks
    for q in ["2001Q1", "2019Q4", "2021Q4", "2024Q4"]:
        sub = panel[
            (panel["quarter_label"] == q) &
            (~panel["pip_code"].isin(["00", "TOTAL"]))
        ]
        total_exp = panel[
            (panel["quarter_label"] == q) &
            (panel["pip_code"] == "TOTAL")
        ]["expanding"]
        total_val = total_exp.values[0] if len(total_exp) else float("nan")
        print(f"  {q}: TOTAL expanding={total_val:.0f}, N supersectors={len(sub)}")

    panel.to_parquet(CACHE_PATH, index=False)
    print(f"\nSaved: {CACHE_PATH}")
    print("Re-run part11_jf_table2.py to regenerate tables with updated data.")
