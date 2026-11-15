"""
part7e_tcu.py -- fetch total capacity utilization (Federal Reserve G.17) and save quarterly means
=================================================================================================
Total capacity utilization (manufacturing, mining, and electric and gas utilities), seasonally
adjusted, percent of capacity; the series FRED calls TCU. FRED refuses this container's requests,
so the source is the Board's G.17 text file of utilization rates.

Source
------
    https://www.federalreserve.gov/releases/g17/ipdisk/utl_sa.txt
File format (verified Nov 15, 2026): a quoted label line per series, then one line per year:
"<series code>" <year> <12 monthly values>, whitespace separated. Total index: code "B50001"
(label "B50001: Total index"); monthly from 1967.

Quarterly aggregation: the mean of the three months; a quarter with fewer than three months is
dropped.

Output
------
    data/cache/tcu_quarterly.parquet   quarter_label, tcu (percent of capacity)

Run: python part7e_tcu.py (standalone data fetch, no prerequisites).
"""

import io
import os
from pathlib import Path

import pandas as pd
import requests

try:
    os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
    pass                                   # interactive: run from "Bartek analysis"

URL    = "https://www.federalreserve.gov/releases/g17/ipdisk/utl_sa.txt"
SERIES = "B50001"
OUT    = Path("data/cache/tcu_quarterly.parquet")


def parse(text, code=SERIES):
    """Monthly series `code` from the G.17 text layout, as a Series indexed by month."""
    rows = []
    for line in io.StringIO(text):
        parts = line.split()
        if len(parts) >= 3 and parts[0] == f'"{code}"':
            year = int(parts[1])
            for m, val in enumerate(parts[2:14], start=1):
                rows.append((pd.Period(year=year, month=m, freq="M"), float(val)))
    if not rows:
        raise ValueError(f"series {code} not found: the G.17 file format may have changed")
    return pd.Series(dict(rows)).sort_index()


def quarterly(monthly):
    """Quarterly means; complete quarters only."""
    q = monthly.groupby(monthly.index.asfreq("Q"))
    out = q.mean()[q.size() == 3]
    return pd.DataFrame({"quarter_label": out.index.astype(str), "tcu": out.to_numpy()})


if __name__ == "__main__":
    print(f"  Fetching: {URL}")
    resp = requests.get(URL, headers={"User-Agent": "Mozilla/5.0"}, timeout=60)
    resp.raise_for_status()
    q = quarterly(parse(resp.text))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    q.to_parquet(OUT, index=False)
    print(f"  {q.quarter_label.iloc[0]}-{q.quarter_label.iloc[-1]}, {len(q)} quarters; "
          f"mean {q.tcu.mean():.1f} %")
    print(f"Saved: {OUT}")
