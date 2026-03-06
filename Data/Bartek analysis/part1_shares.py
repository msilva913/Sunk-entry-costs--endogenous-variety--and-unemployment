"""
part1_shares.py — Employment shares  (weights  omega_{s,j,t0})
===============================================================
Fetches QCEW annual data for BASE_YEAR from the BLS bulk download,
builds the state x supersector employment share matrix, and saves it
to disk.

Output
------
    data/instruments/shares_base{BASE_YEAR}.parquet

Columns in output:
    state_fips       2-digit state FIPS code (e.g. "06" = California)
    industry_code    2-digit supersector code (e.g. "30" = Manufacturing)
    emp_state_ind    state employment in the supersector
    emp_state_total  total state private employment (denominator of share)
    emp_nat_ind      national employment in the supersector
    emp_nat_loo      national employment minus this state's (LOO denominator)
    share            omega_{s,j} = emp_state_ind / emp_state_total
    share_of_nat     emp_state_ind / emp_nat_ind  (used for LOO diagnostics)

Run
---
    python part1_shares.py
"""

import pandas as pd
import os      
#os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment\Data\Bartek analysis"
#)
from pathlib import Path
try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:
    BASE_DIR = Path.cwd()   # whatever directory VS Code opened

DEFAULT_CACHE_DIR  = BASE_DIR / "data" / "cache"
DEFAULT_OUTPUT_DIR = BASE_DIR / "data" / "instruments"
BDS_FILE = BASE_DIR / "data" / "raw" / "bds2023_sec_nat.csv"

from construct_delta_instrument import (
    BASE_YEAR,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_employment_shares,
    shares_to_grid,
)

# -----------------------------------------------------------------------
# Compute
# -----------------------------------------------------------------------
shares = build_employment_shares(
    base_year = BASE_YEAR,
    cache_dir = DEFAULT_CACHE_DIR,
)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
shares.to_parquet(SHARES_PATH, index=False)
print(f"\nSaved: {SHARES_PATH}  ({len(shares):,} rows)")

# -----------------------------------------------------------------------
# Inspect
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Part 1 — Employment shares  (base year {BASE_YEAR})")
print("=" * 60)

print(f"\nShape: {shares.shape}  "
      f"({shares['state_fips'].nunique()} states × "
      f"{shares['industry_code'].nunique()} supersectors)")

# ── Grid view ────────────────────────────────────────────────────────────
print("\n--- Share grid (omega_{s,j})  [states x supersectors] ---")
grid = shares_to_grid(shares)
with pd.option_context("display.float_format", "{:.4f}".format,
                       "display.max_columns", 12,
                       "display.width", 120):
    print(grid.to_string())

# ── Summary statistics on the grid ───────────────────────────────────────
print("\n--- Column means (average share by supersector across states) ---")
col_means = grid.mean().sort_values(ascending=False)
for name, s in col_means.items():
    bar = "█" * int(s * 200)
    print(f"  {name:<36}  {s:.4f}  {bar}")

print("\n--- Row sums (should all be ≈ 1.0) ---")
row_sums = grid.sum(axis=1)
ok = (row_sums - 1.0).abs().max() < 0.001
print(f"  Max deviation from 1.0: {(row_sums - 1.0).abs().max():.6f}  {'✓' if ok else '✗ WARNING'}")

# ── LOO coverage ─────────────────────────────────────────────────────────
print("\n--- LOO coverage: states that dominate a supersector nationally ---")
dominant = shares[shares["share_of_nat"] >= 0.15].sort_values(
    "share_of_nat", ascending=False
)
if dominant.empty:
    print("  No state accounts for ≥15% of any national supersector.")
else:
    dominant_display = dominant.assign(
        state       = lambda d: d["state_fips"].map(FIPS2D_TO_STATE),
        supersector = lambda d: d["industry_code"].map(INDUSTRY_LABELS),
    )[["state", "supersector", "share_of_nat", "emp_state_ind", "emp_nat_loo"]]
    print(dominant_display.round(4).to_string(index=False))
    print("  (These cells receive the largest LOO correction.)")
    
import matplotlib.pyplot as plt
import seaborn as sns      
plt.figure(figsize=(20, 5))
g = sns.clustermap(grid.T, cmap="viridis", figsize=(22, 6),
                   linewidths=0, cbar_pos=(0.02, 0.8, 0.03, 0.15),
                   dendrogram_ratio=0.1)
g.ax_heatmap.set_xlabel("State")
g.ax_heatmap.set_ylabel("Industry")
#plt.show()
shares



