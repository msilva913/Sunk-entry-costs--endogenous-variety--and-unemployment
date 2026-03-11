"""
part7_nfci.py — Fetch and prepare NFCI data for LP interaction term
====================================================================
Downloads the Chicago Fed National Financial Conditions Index (NFCI)
and its risk subindex, aggregates from weekly to quarterly, and saves
a parquet for use by part5_lp.py.

Source
------
Chicago Fed direct CSV download:
    https://www.chicagofed.org/-/media/publications/nfci/nfci-data-series-csv.csv

File format (verified against Chicago Fed download 2026-03-11):
    - Column "Friday_of_Week"       : week-ending Friday date, format YYYY-MM-DD
    - Column "NFCI"                 : headline National Financial Conditions Index
    - Column "ANFCI"                : adjusted NFCI (controls for macro conditions)
    - Column "Risk"                 : risk subindex (volatility + funding risk)
    - Column "Credit"               : credit subindex
    - Column "Leverage"             : leverage subindex
    - Column "Nonfinancial_Leverage": nonfinancial leverage subindex (present, not used)

    All index values are mean-zero, SD-one over 1971–present.
    Positive = tighter-than-average financial conditions.

Why the risk subindex?
-----------------------
The headline NFCI is highly collinear with recession severity: it
tightens precisely when unemployment rises for aggregate reasons.
Including NFCI as an interaction term in the LP would therefore make
it difficult to distinguish financial amplification of the separation
shock from a generic deep-recession effect.

The risk subindex captures volatility and funding risk in the financial
sector — components that can tighten for purely financial reasons
(e.g., LTCM 1998, money market freeze 2008) independently of real
activity. Using the risk subindex as the interaction variable gives
cleaner identification of a financial amplification channel.

Quarterly aggregation
----------------------
Each quarterly value is the simple mean of the weekly values whose
week-ending Friday falls within the quarter. Weeks that straddle
quarter boundaries are assigned to the quarter containing the Friday.

ANFCI note
-----------
The adjusted NFCI (ANFCI) removes variation attributable to economic
activity and inflation before computing the index. It is saved to the
parquet for convenience — it is an alternative to the risk subindex
for isolating financial conditions from macro conditions. However, the
ANFCI adjustment is estimated with noise and subject to larger
revisions than the raw subindices. We recommend the risk subindex as
the primary interaction variable.

Outputs
-------
    data/cache/nfci_quarterly.parquet
        Columns:
            quarter_label  str    e.g. "2008Q4"
            nfci           float  headline NFCI, quarterly mean
            nfci_risk      float  risk subindex, quarterly mean
            nfci_credit    float  credit subindex, quarterly mean
            nfci_leverage  float  leverage subindex, quarterly mean
            anfci          float  adjusted NFCI, quarterly mean

Sample
------
    1990Q1 – latest available (currently 2025Q4 or later)
    START_QUARTER = "1990Q1" matches instrument sample lower bound.

Run
---
    python part7_nfci.py

No prerequisites — standalone data fetch.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents"
        / "GitHub"
        / "Sunk-entry-costs--endogenous-variety--and-unemployment"
        / "Data"
        / "Bartek analysis"
    )

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NFCI_URL = (
    "https://www.chicagofed.org"
    "/-/media/publications/nfci/nfci-data-series-csv.csv"
)

# Expected column names in the Chicago Fed CSV (as of 2026-03-11).
# If the file format changes, assertions below will catch it immediately.
COL_DATE            = "Friday_of_Week"
COL_NFCI            = "NFCI"
COL_ANFCI           = "ANFCI"
COL_RISK            = "Risk"
COL_CREDIT          = "Credit"
COL_LEVERAGE        = "Leverage"
COL_NONFIN_LEVERAGE = "Nonfinancial_Leverage"   # present but not used in LP

REQUIRED_COLS = [COL_DATE, COL_NFCI, COL_ANFCI, COL_RISK, COL_CREDIT, COL_LEVERAGE]

START_QUARTER = "1990Q1"

CACHE_DIR   = Path("data/cache")
RESULTS_DIR = Path("data/results")
CACHE_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

OUT_PARQUET = CACHE_DIR / "nfci_quarterly.parquet"
OUT_PLOT    = RESULTS_DIR / "nfci_quarterly.png"


# ---------------------------------------------------------------------------
# Step 1: Fetch
# ---------------------------------------------------------------------------
def fetch_nfci(url: str) -> pd.DataFrame:
    """
    Download the Chicago Fed NFCI CSV and return as a raw DataFrame.
    Prints the first rows and column names for verification.
    """
    print(f"  Fetching: {url}")
    try:
        raw = pd.read_csv(url)
    except Exception as e:
        print(f"\n  ERROR: Could not fetch NFCI data.\n  {e}")
        print(
            "\n  Fallback: download the CSV manually from:\n"
            f"    {url}\n"
            "  and save it to data/raw/nfci-data-series-csv.csv\n"
            "  then re-run this script."
        )
        local = Path("data/raw/nfci-data-series-csv.csv")
        if local.exists():
            print(f"  Found local file: {local}. Loading.")
            raw = pd.read_csv(local)
        else:
            sys.exit(1)

    print(f"\n  Raw shape: {raw.shape}")
    print(f"  Columns:   {raw.columns.tolist()}")
    print(f"\n  First 3 rows:")
    print(raw.head(3).to_string())
    print(f"\n  Last 3 rows:")
    print(raw.tail(3).to_string())
    return raw


# ---------------------------------------------------------------------------
# Step 2: Validate and parse
# ---------------------------------------------------------------------------
def validate_and_parse(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Assert expected columns exist, parse dates, return clean DataFrame.
    Raises AssertionError with a helpful message if format has changed.
    """
    # Check all required columns are present
    missing = [c for c in REQUIRED_COLS if c not in raw.columns]
    if missing:
        raise AssertionError(
            f"\n  FORMAT ERROR: Expected columns not found in NFCI CSV.\n"
            f"  Missing: {missing}\n"
            f"  Found:   {raw.columns.tolist()}\n"
            f"  The Chicago Fed may have changed the CSV format.\n"
            f"  Update REQUIRED_COLS and column name constants in this file."
        )

    df = raw[REQUIRED_COLS].copy()
    df = df.rename(columns={
        COL_DATE:     "date",
        COL_NFCI:     "nfci",
        COL_ANFCI:    "anfci",
        COL_RISK:     "nfci_risk",
        COL_CREDIT:   "nfci_credit",
        COL_LEVERAGE: "nfci_leverage",
    })

    # Parse dates
    df["date"] = pd.to_datetime(df["date"], format="%m/%d/%Y")

    # Numeric conversion — coerce to float, warn on NaNs
    for col in ["nfci", "anfci", "nfci_risk", "nfci_credit", "nfci_leverage"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    n_na = df[["nfci", "nfci_risk"]].isna().sum()
    if n_na.any():
        print(f"\n  WARNING: NaN values found:\n{n_na}")

    df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)

    print(f"\n  Parsed date range: {df['date'].min().date()} – "
          f"{df['date'].max().date()}")
    print(f"  Total weekly obs: {len(df)}")

    # Sanity check: NFCI should be mean ~0, SD ~1 over long sample
    print(f"\n  NFCI  — mean={df['nfci'].mean():.3f}  "
          f"SD={df['nfci'].std():.3f}  "
          f"min={df['nfci'].min():.3f}  "
          f"max={df['nfci'].max():.3f}")
    print(f"  Risk  — mean={df['nfci_risk'].mean():.3f}  "
          f"SD={df['nfci_risk'].std():.3f}  "
          f"min={df['nfci_risk'].min():.3f}  "
          f"max={df['nfci_risk'].max():.3f}")

    return df


# ---------------------------------------------------------------------------
# Step 3: Weekly → quarterly
# ---------------------------------------------------------------------------
def to_quarterly(df: pd.DataFrame, start_quarter: str) -> pd.DataFrame:
    """
    Assign each weekly observation to its quarter (by Friday date),
    compute simple means within quarter, restrict to start_quarter onward.

    Returns DataFrame with columns:
        quarter_label, nfci, anfci, nfci_risk, nfci_credit, nfci_leverage
    """
    # Assign quarter using the Friday date
    df["period"] = df["date"].dt.to_period("Q")
    df["quarter_label"] = df["period"].astype(str)

    # Count weeks per quarter for QA
    weeks_per_q = df.groupby("quarter_label")["date"].count().rename("n_weeks")

    # Quarterly means
    value_cols = ["nfci", "anfci", "nfci_risk", "nfci_credit", "nfci_leverage"]
    qdf = (
        df.groupby("quarter_label")[value_cols]
        .mean()
        .reset_index()
    )
    qdf = qdf.merge(weeks_per_q, on="quarter_label")

    # Warn on quarters with fewer than 12 weeks (incomplete)
    incomplete = qdf[qdf["n_weeks"] < 12]
    if len(incomplete):
        print(f"\n  WARNING: {len(incomplete)} quarters with < 12 weeks "
              f"(may be incomplete):")
        print(incomplete[["quarter_label", "n_weeks"]].to_string())

    # Sort and restrict to start_quarter
    qdf = qdf.sort_values("quarter_label")
    qdf = qdf[qdf["quarter_label"] >= start_quarter].reset_index(drop=True)

    print(f"\n  Quarterly obs after {start_quarter} cutoff: {len(qdf)}")
    print(f"  Quarter range: {qdf['quarter_label'].iloc[0]} – "
          f"{qdf['quarter_label'].iloc[-1]}")
    print(f"\n  Quarterly NFCI summary:")
    print(qdf[["quarter_label"] + value_cols].describe().round(3).to_string())

    return qdf[["quarter_label"] + value_cols]


# ---------------------------------------------------------------------------
# Step 4: Plot
# ---------------------------------------------------------------------------
def plot_nfci(qdf: pd.DataFrame, out_path: Path) -> None:
    """
    Two-panel figure: headline NFCI and risk subindex over time,
    quarterly frequency, with NBER recession shading.
    """
    dts = [pd.Period(q, freq="Q").to_timestamp() for q in qdf["quarter_label"]]

    recessions = [
        ("1990-07-01", "1991-03-01"),
        ("2001-03-01", "2001-11-01"),
        ("2007-12-01", "2009-06-01"),
        ("2020-02-01", "2020-04-01"),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    for ax, col, label, color in [
        (axes[0], "nfci",     "NFCI (headline)",    "#1f77b4"),
        (axes[1], "nfci_risk","NFCI risk subindex", "#d62728"),
    ]:
        ax.plot(dts, qdf[col].values, color=color, linewidth=1.5, label=label)
        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.fill_between(dts, qdf[col].values, 0,
                        where=(qdf[col].values > 0),
                        color=color, alpha=0.15, label="Tighter than avg")
        for rs, re in recessions:
            ax.axvspan(pd.Timestamp(rs), pd.Timestamp(re),
                       alpha=0.10, color="grey")
        ax.set_ylabel("Index value (SD units)", fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)

    axes[1].set_xlabel("Quarter", fontsize=9)

    caption = (
        "Notes: Quarterly averages of weekly Chicago Fed NFCI values. "
        "The headline NFCI (top) captures overall financial conditions across "
        "money markets, debt and equity markets, and the banking system. "
        "The risk subindex (bottom) captures volatility and funding risk in "
        "the financial sector — the preferred interaction variable for the LP "
        "because it isolates financial amplification from recession severity. "
        "Both indices are mean-zero, SD-one over their full sample (1971–present). "
        "Positive values indicate tighter-than-average conditions. "
        "Grey shading marks NBER recessions."
    )
    fig.text(0.5, -0.03, caption, ha="center", va="top", fontsize=8,
             wrap=True, transform=fig.transFigure,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                       edgecolor="#cccccc", linewidth=0.8),
             multialignment="left")

    fig.suptitle(
        "NFCI and Risk Subindex — Quarterly averages\n"
        f"({qdf['quarter_label'].iloc[0]} – {qdf['quarter_label'].iloc[-1]})",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 7 — NFCI Data Preparation")
print("=" * 60)

print("\n[1] Fetching NFCI data")
raw = fetch_nfci(NFCI_URL)

print("\n[2] Validating and parsing")
weekly = validate_and_parse(raw)

print("\n[3] Aggregating to quarterly")
quarterly = to_quarterly(weekly, START_QUARTER)

print("\n[4] Saving parquet")
quarterly.to_parquet(OUT_PARQUET, index=False)
print(f"  Saved: {OUT_PARQUET}  ({len(quarterly)} rows)")
print(f"  Columns: {quarterly.columns.tolist()}")

print("\n[5] Plotting")
plot_nfci(quarterly, OUT_PLOT)

print("\n" + "=" * 60)
print("Output ready for part5_lp.py")
print(f"  Load with: pd.read_parquet('{OUT_PARQUET}')")
print(f"  Key column for LP interaction: 'nfci_risk'")
print("=" * 60)