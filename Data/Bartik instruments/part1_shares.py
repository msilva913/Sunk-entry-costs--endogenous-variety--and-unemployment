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
import os 
import pandas as pd

os.chdir(r"C:\Users\msilv\Documents\GitHub\Sunk-entry-costs--endogenous-variety--and-unemployment\Data\Bartik instruments")

from construct_delta_instrument import (
    BASE_YEAR,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    SHARES_PATH,
    INDUSTRY_LABELS,
    FIPS2D_TO_STATE,
    build_employment_shares,
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

print("\nColumns:\n  " + "  ".join(shares.columns))

print("\n--- National employment by supersector ---")
nat_emp = (
    shares.groupby("industry_code")["emp_nat_ind"]
    .first()
    .rename(index=INDUSTRY_LABELS)
    .sort_values(ascending=False)
)
for name, emp in nat_emp.items():
    print(f"  {name:<36}  {emp:>12,.0f}")

print("\n--- Top 15 state-supersector weights (share of state employment) ---")
top = (
    shares
    .nlargest(15, "share")
    .assign(
        state       = lambda d: d["state_fips"].map(FIPS2D_TO_STATE),
        supersector = lambda d: d["industry_code"].map(INDUSTRY_LABELS),
    )
    [["state", "supersector", "share", "share_of_nat",
      "emp_state_ind", "emp_nat_loo"]]
    .round({"share": 4, "share_of_nat": 4})
)
print(top.to_string(index=False))

print("\n--- Share matrix: mean weight by supersector across states ---")
mean_share = (
    shares.groupby("industry_code")["share"]
    .mean()
    .rename(index=INDUSTRY_LABELS)
    .sort_values(ascending=False)
)


for name, s in mean_share.items():
    bar = "█" * int(s * 200)
    print(f"  {name:<36}  {s:.4f}  {bar}")

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

shares.to_csv("shares_base_2006.csv")

def shares_to_grid(shares: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot the shares DataFrame into a 50 x 10 grid.

    Parameters
    ----------
    shares : output of build_employment_shares()

    Returns
    -------
    pd.DataFrame
        Rows   : states, indexed by abbreviation (AL, AK, ...), sorted
        Columns: supersectors, labelled by name (Mining, Construction, ...)
        Values : share  (omega_{s,j})
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

shares_grid = shares_to_grid(shares)
print(shares_grid)
shares.to_csv("shares_base_2006.csv")
shares_grid.to_csv("shares_grid_base_2006.csv")