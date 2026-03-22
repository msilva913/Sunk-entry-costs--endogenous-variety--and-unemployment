"""
part7b_sloos.py — Fetch and prepare SLOOS data for LP interaction term
========================================================================
Downloads the Senior Loan Officer Opinion Survey (SLOOS) C&I net tightening
index from FRED and aggregates from quarterly survey responses to a clean
quarterly panel for use by part5d_lp.py as an alternative financial
conditions measure to NFCI.

Source
------
Federal Reserve Economic Data (FRED):
    DRTSCILM — Net Percentage of Domestic Respondents Tightening Standards
               for Commercial and Industrial Loans to Large and Middle-Market Firms
    DRTSCIS  — Net Percentage of Domestic Respondents Tightening Standards
               for Commercial and Industrial Loans to Small Firms

Both series are quarterly (survey conducted quarterly), representing the
net percentage of respondents reporting tighter lending standards.

Why SLOOS as alternative to NFCI?
-----------------------------------
While NFCI captures broad financial conditions across markets and institutions,
SLOOS isolates supply-side behavior in the core banking system: how loan officers
actually adjust lending standards in response to risk, funding availability, or
prudential constraints. SLOOS is directly measured rather than constructed from
market prices, making it less contaminated by demand-side movements.

The C&I lending standard measure (especially for large/middle-market firms)
is preferred for LP identification because:
  1. It reflects bank capital and funding constraints more directly
  2. C&I lending is procyclical — tightening during weak demand
  3. Large/middle-market firms are less affected by information asymmetry,
     isolating the pure credit supply channel

Quarterly data
---------------
SLOOS is already quarterly (no aggregation needed). The survey is conducted
once per quarter, and FRED stores the values with quarter-end or quarter-start
dates. The script parses these dates and assigns to quarter_label format
(e.g. "1990Q2").

Outputs
-------
    data/cache/sloos_quarterly.parquet
        Columns:
            quarter_label        str    e.g. "1990Q2"
            sloos_ci             float  C&I net tightening, large/middle-market (primary)
            sloos_ci_small       float  C&I net tightening, small firms
            sloos_ci_demeaned    float  (demeaned version; computed in LP, stored here for reference)

Sample
------
    1990Q2 – latest available (currently ~2025Q4)
    START_QUARTER = "1990Q1" matches instrument sample lower bound;
    first actual SLOOS obs is typically 1990Q2.

Run
---
    python part7b_sloos.py

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
# FRED series IDs
FRED_SERIES_LARGE  = "DRTSCILM"   # Large and middle-market firms
FRED_SERIES_SMALL  = "DRTSCIS"    # Small firms

# Direct CSV download URLs from FRED
FRED_BASE_URL = "https://fred.stlouisfed.org/graph/fredgraph.csv"
SLOOS_URL_LARGE = (
    f"{FRED_BASE_URL}?id={FRED_SERIES_LARGE}"
    f"&cosd=1990-01-01&coed=2025-12-31"
)
SLOOS_URL_SMALL = (
    f"{FRED_BASE_URL}?id={FRED_SERIES_SMALL}"
    f"&cosd=1990-01-01&coed=2025-12-31"
)

START_QUARTER = "1990Q1"

CACHE_DIR   = Path("data/cache")
RESULTS_DIR = Path("data/results")
CACHE_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

OUT_PARQUET = CACHE_DIR / "sloos_quarterly.parquet"
OUT_PLOT    = RESULTS_DIR / "sloos_quarterly.png"


# ---------------------------------------------------------------------------
# Step 1: Fetch
# ---------------------------------------------------------------------------
def fetch_sloos_series(url: str, series_id: str) -> pd.DataFrame:
    """
    Download a SLOOS series CSV from FRED.
    Returns DataFrame with columns: date, value.
    """
    print(f"  Fetching {series_id}: {url}")
    try:
        raw = pd.read_csv(url)
    except Exception as e:
        print(f"\n  ERROR: Could not fetch {series_id}.\n  {e}")
        print(
            f"\n  Fallback: download the CSV manually from:\n"
            f"    {url}\n"
            f"  and save it to data/raw/{series_id.lower()}.csv\n"
            "  then re-run this script."
        )
        local = Path(f"data/raw/{series_id.lower()}.csv")
        if local.exists():
            print(f"  Found local file: {local}. Loading.")
            raw = pd.read_csv(local)
        else:
            return None

    # FRED CSVs have one column 'DATE' and one column with the series ID
    if "DATE" in raw.columns and series_id in raw.columns:
        raw = raw[["DATE", series_id]].copy()
        raw.columns = ["date", "value"]
    else:
        print(f"  Unexpected FRED CSV format for {series_id}")
        print(f"  Columns found: {raw.columns.tolist()}")
        return None

    print(f"  Raw shape: {raw.shape}")
    print(f"  First 3 rows:")
    print(raw.head(3).to_string())
    print(f"  Last 3 rows:")
    print(raw.tail(3).to_string())
    return raw


# ---------------------------------------------------------------------------
# Step 2: Validate and parse
# ---------------------------------------------------------------------------
def validate_and_parse(raw_large: pd.DataFrame, raw_small: pd.DataFrame) -> pd.DataFrame:
    """
    Parse dates, convert to numeric, align two series, return clean DataFrame.
    """
    # Parse dates and convert values
    df_large = raw_large.copy()
    df_small = raw_small.copy()

    df_large["date"] = pd.to_datetime(df_large["date"])
    df_small["date"] = pd.to_datetime(df_small["date"])

    df_large["value"] = pd.to_numeric(df_large["value"], errors="coerce")
    df_small["value"] = pd.to_numeric(df_small["value"], errors="coerce")

    # Rename for clarity
    df_large = df_large.rename(columns={"value": "sloos_ci"})
    df_small = df_small.rename(columns={"value": "sloos_ci_small"})

    # Merge on date
    df = df_large[["date", "sloos_ci"]].merge(
        df_small[["date", "sloos_ci_small"]],
        on="date",
        how="outer",
    ).sort_values("date").reset_index(drop=True)

    # Check for NaN
    n_na = df[["sloos_ci", "sloos_ci_small"]].isna().sum()
    if n_na.any():
        print(f"\n  WARNING: NaN values found:\n{n_na}")

    df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)

    print(f"\n  Parsed date range: {df['date'].min().date()} – "
          f"{df['date'].max().date()}")
    print(f"  Total quarterly obs: {len(df)}")

    # Summary statistics
    print(f"\n  SLOOS C&I (large/middle):")
    print(f"    mean={df['sloos_ci'].mean():.3f}  "
          f"SD={df['sloos_ci'].std():.3f}  "
          f"min={df['sloos_ci'].min():.3f}  "
          f"max={df['sloos_ci'].max():.3f}")
    print(f"  SLOOS C&I (small):")
    print(f"    mean={df['sloos_ci_small'].mean():.3f}  "
          f"SD={df['sloos_ci_small'].std():.3f}  "
          f"min={df['sloos_ci_small'].min():.3f}  "
          f"max={df['sloos_ci_small'].max():.3f}")

    return df


# ---------------------------------------------------------------------------
# Step 3: Assign quarter labels
# ---------------------------------------------------------------------------
def assign_quarters(df: pd.DataFrame, start_quarter: str) -> pd.DataFrame:
    """
    Assign each date to its quarter, create quarter_label (e.g. "2008Q4"),
    restrict to start_quarter onward, drop duplicates (keep first if any).

    Returns DataFrame with columns:
        quarter_label, sloos_ci, sloos_ci_small
    """
    # Assign quarter from date
    df["period"] = df["date"].dt.to_period("Q")
    df["quarter_label"] = df["period"].astype(str)

    # Drop date and period columns
    df = df[["quarter_label", "sloos_ci", "sloos_ci_small"]].copy()

    # If there are multiple dates per quarter (unlikely but possible),
    # keep first occurrence
    df = df.drop_duplicates(subset=["quarter_label"], keep="first")

    # Sort and restrict to start_quarter
    df = df.sort_values("quarter_label")
    df = df[df["quarter_label"] >= start_quarter].reset_index(drop=True)

    print(f"\n  Quarterly obs after {start_quarter} cutoff: {len(df)}")
    if len(df) > 0:
        print(f"  Quarter range: {df['quarter_label'].iloc[0]} – "
              f"{df['quarter_label'].iloc[-1]}")
        print(f"\n  SLOOS quarterly summary:")
        print(df[["quarter_label", "sloos_ci", "sloos_ci_small"]].describe().round(3).to_string())

    return df


# ---------------------------------------------------------------------------
# Step 4: Plot
# ---------------------------------------------------------------------------
def plot_sloos(qdf: pd.DataFrame, out_path: Path) -> None:
    """
    Two-panel figure: SLOOS C&I net tightening for large/middle-market
    and small firms over time, with NBER recession shading.
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
        (axes[0], "sloos_ci",       "C&I Tightening — Large/Middle-Market", "#1f77b4"),
        (axes[1], "sloos_ci_small", "C&I Tightening — Small Firms",        "#d62728"),
    ]:
        ax.plot(dts, qdf[col].values, color=color, linewidth=1.5, label=label)
        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.fill_between(dts, qdf[col].values, 0,
                        where=(qdf[col].values > 0),
                        color=color, alpha=0.15, label="Tightening (>0)")
        for rs, re in recessions:
            ax.axvspan(pd.Timestamp(rs), pd.Timestamp(re),
                       alpha=0.10, color="grey")
        ax.set_ylabel("Net % tightening", fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)

    axes[1].set_xlabel("Quarter", fontsize=9)

    caption = (
        "Notes: Quarterly values from the Federal Reserve's Senior Loan Officer "
        "Opinion Survey (SLOOS) on net percentage of domestic respondents tightening "
        "credit standards for commercial and industrial (C&I) loans. Top panel: "
        "Large and middle-market firms; bottom panel: small firms. A positive value "
        "indicates that tightening respondents outnumber easing respondents. The SLOOS "
        "is conducted once per quarter (early in the quarter) and reflects bank "
        "supply-side behavior rather than market prices. Positive values are associated "
        "with constrained credit supply. Grey shading marks NBER recessions."
    )
    fig.text(0.5, -0.03, caption, ha="center", va="top", fontsize=8,
             wrap=True, transform=fig.transFigure,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                       edgecolor="#cccccc", linewidth=0.8),
             multialignment="left")

    fig.suptitle(
        "SLOOS C&I Net Tightening — Quarterly\n"
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
print("Part 7b — SLOOS Data Preparation")
print("=" * 60)

print("\n[1] Fetching SLOOS data from FRED")
raw_large = fetch_sloos_series(SLOOS_URL_LARGE, FRED_SERIES_LARGE)
raw_small = fetch_sloos_series(SLOOS_URL_SMALL, FRED_SERIES_SMALL)

if raw_large is None or raw_small is None:
    print("\n  FATAL: Could not fetch one or both SLOOS series.")
    sys.exit(1)

print("\n[2] Validating and parsing")
weekly = validate_and_parse(raw_large, raw_small)

print("\n[3] Assigning quarters and restricting to sample")
quarterly = assign_quarters(weekly, START_QUARTER)

if quarterly.empty:
    print("\n  FATAL: No data after filtering to start quarter.")
    sys.exit(1)

print("\n[4] Saving parquet")
quarterly.to_parquet(OUT_PARQUET, index=False)
print(f"  Saved: {OUT_PARQUET}  ({len(quarterly)} rows)")
print(f"  Columns: {quarterly.columns.tolist()}")

print("\n[5] Plotting")
plot_sloos(quarterly, OUT_PLOT)

print("\n" + "=" * 60)
print("Output ready for part5d_lp.py")
print(f"  Load with: pd.read_parquet('{OUT_PARQUET}')")
print(f"  Key column for LP interaction: 'sloos_ci'")
print("=" * 60)
