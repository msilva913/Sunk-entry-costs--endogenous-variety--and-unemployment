"""
run_v2_locally.py
=================
Run this script on your LOCAL machine (not in the Cowork sandbox) to:

  1. Fetch BEA quarterly industry value-added data (requires BEA_API_KEY env var)
  2. Re-run part2b with the enriched v2 residualization
  3. Re-run part3 to rebuild Bartik instruments
  4. Re-run part5 to get updated IRFs

Prerequisites:
  - pip install pandas pyarrow requests matplotlib statsmodels linearmodels
  - export BEA_API_KEY=<your key from apps.bea.gov>  (or set in Windows env vars)
  - Run from the project working directory:
      cd "Data/Bartek analysis"
      python run_v2_locally.py

This script does NOT modify any code — it calls the existing pipeline scripts
in sequence and prints a summary of cross-instrument correlations and IRF peaks.
"""

import subprocess
import sys
import os
import pandas as pd
from pathlib import Path

SCRIPTS = ["part2b_shock_comovement.py", "part3_resid_instruments.py", "part5_lp.py"]
key = "7878D119-FFC4-411E-969D-3B168F53A1E0"
#key = os.environ["BEA_API_KEY"]
def check_prereqs():
    #key = os.environ.get("BEA_API_KEY", "")
    if not key:
        print("ERROR: BEA_API_KEY environment variable not set.")
        print("  Register at https://apps.bea.gov/API/signup/ (free, instant approval)")
        print("  Then: export BEA_API_KEY=<your-key>   (or set in Windows env vars)")
        sys.exit(1)
    print(f"BEA_API_KEY: {key[:8]}... (set)")

    bea_cache = Path("data/cache/bea_va_quarterly_12ind.parquet")
    if bea_cache.exists():
        try:
            df = pd.read_parquet(bea_cache)
            coverage_start = df["quarter_label"].min()
            print(f"BEA cache exists, starts {coverage_start}.")
            if coverage_start > "2002Q1":
                print("  WARNING: BEA cache starts after 2002Q1 — this is the 2005-start limitation.")
                print("  The v2 spec will run but cover only 2005Q2 onward.")
                print("  v1 spec (full sample) will be run as primary robustness check.")
            else:
                print(f"  Coverage is sufficient (starts {coverage_start}).")
        except Exception as e:
            print(f"  BEA cache is corrupted ({e}). Will delete and re-fetch.")
            bea_cache.unlink()

def run_script(script):
    print(f"\n{'='*60}")
    print(f"Running {script} ...")
    print('='*60)
    result = subprocess.run(
        [sys.executable, script],
        capture_output=False,
        text=True
    )
    if result.returncode != 0:
        print(f"\nERROR: {script} exited with code {result.returncode}")
        print("Fix the error above before continuing.")
        sys.exit(result.returncode)
    print(f"\n{script} completed successfully.")

def print_summary():
    print("\n" + "="*60)
    print("SUMMARY: Key results from v2 run")
    print("="*60)

    # Instrument correlations from part3
    corr_plot = Path("data/results/instruments_resid_corr_matrix.png")
    if corr_plot.exists():
        print(f"\nInstrument correlation matrix plot: {corr_plot}")

    # Key IRF peaks from part5
    results_dir = Path("data/results")
    for fname, label in [
        ("lp_irf_delta_vacancy.csv", "δ → vacancies"),
        ("lp_irf_ld_vacancy.csv",    "LD → vacancies"),
        ("lp_irf_qu_vacancy.csv",    "QU → vacancies (placebo)"),
    ]:
        fpath = results_dir / fname
        if fpath.exists():
            df = pd.read_csv(fpath)
            peak_idx = df["beta"].abs().idxmax()
            peak = df.loc[peak_idx]
            sig = "***" if peak["pval"] < 0.01 else "**" if peak["pval"] < 0.05 else "*" if peak["pval"] < 0.1 else ""
            std_beta = peak["beta"] * peak["instr_sd"]
            print(f"  {label}: peak β={std_beta:+.3f} pp at h={int(peak['h'])} (p={peak['pval']:.3f}) {sig}")

    # Check if v2 or v1 was used
    bea_cache = Path("data/cache/bea_va_quarterly_12ind.parquet")
    if bea_cache.exists():
        df = pd.read_parquet(bea_cache)
        print(f"\nResidual spec: v2 (BEA industry VA, {df['quarter_label'].min()}–{df['quarter_label'].max()})")
        print("  Compare with v1 (Dlog p only) for robustness.")
        print("  Key test: did r(δ,LD) drop further from the v1 value of 0.39?")
    else:
        print("\nResidual spec: v1 fallback (BEA fetch failed — check API key/network)")

    print(f"\nAll plots saved to: data/results/")
    print("Key plot: data/results/lp_irf_vacancy_delta_ld.png")


if __name__ == "__main__":
    print("Track B: v2 BEA-enriched residualization pipeline")
    print("="*60)
    check_prereqs()

    for script in SCRIPTS:
        if not Path(script).exists():
            print(f"ERROR: {script} not found in current directory.")
            print("Make sure you are running from: Data/Bartek analysis/")
            sys.exit(1)

    for script in SCRIPTS:
        run_script(script)

    print_summary()
    print("\nDone. Check data/results/ for updated IRF plots and CSVs.")
