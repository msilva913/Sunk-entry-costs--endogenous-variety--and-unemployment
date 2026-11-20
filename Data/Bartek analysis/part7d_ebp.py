"""
part7d_ebp.py -- fetch the Gilchrist-Zakrajsek credit spread and excess bond premium
====================================================================================
Downloads the monthly GZ credit spread and excess bond premium (EBP), updated by the
Federal Reserve Board (Favara, Gilchrist, Lewis and Zakrajsek, FEDS Notes 2016; citation
unverified), and saves quarterly means for q3_delta_finance.py.

Source
------
    https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv
File format (verified Nov 10, 2026): columns date (YYYY-MM-01), gz_spread, ebp, est_prob;
monthly from 1973-01. gz_spread and ebp are in percentage points.

Quarterly aggregation: the mean of the three months of the quarter; a quarter with fewer
than three months is dropped.

Output
------
    data/cache/ebp_quarterly.parquet   quarter_label, gz_spread, ebp (percentage points)

Run: python part7d_ebp.py (standalone data fetch, no prerequisites).
"""

import os
from pathlib import Path

import pandas as pd

try:
    os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
    pass                                   # interactive: run from "Bartek analysis"

EBP_URL  = "https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv"
REQUIRED = ["date", "gz_spread", "ebp"]
OUT      = Path("data/cache/ebp_quarterly.parquet")


def quarterly(raw):
    """Quarterly means of gz_spread and ebp; complete quarters only."""
    df = raw[REQUIRED].copy()
    df["quarter"] = pd.to_datetime(df["date"]).dt.to_period("Q")
    g = df.groupby("quarter")
    out = g[["gz_spread", "ebp"]].mean()
    out = out[g.size() == 3]
    out.insert(0, "quarter_label", out.index.astype(str))
    return out.reset_index(drop=True)


if __name__ == "__main__":
    print(f"  Fetching: {EBP_URL}")
    raw = pd.read_csv(EBP_URL, storage_options={"User-Agent": "Mozilla/5.0"})   # the default agent gets 403
    missing = [c for c in REQUIRED if c not in raw.columns]
    if missing:
        raise ValueError(f"EBP file format changed; missing columns {missing}")
    q = quarterly(raw)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    q.to_parquet(OUT, index=False)
    print(f"  {q.quarter_label.iloc[0]}-{q.quarter_label.iloc[-1]}, {len(q)} quarters")
    print(f"Saved: {OUT}")
