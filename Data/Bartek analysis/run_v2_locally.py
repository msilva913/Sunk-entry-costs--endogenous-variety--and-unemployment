"""
run_v2_locally.py
=================
Run this script on your LOCAL machine to execute the v2 enriched residualization:

  1. Fetch BEA Real VA by Industry from FRED (requires FRED_API_KEY env var)
  2. part2b: enriched residualization (dlog_va_lag + log_theta_lag)
  3. part3: rebuild Bartik instruments from v2 residuals
  4. part5: updated IRFs

FRED API keys are FREE: https://fred.stlouisfed.org/docs/api/api_key.html


After first run the industry VA is cached; FRED_API_KEY not needed again.
"""

import subprocess
import sys
import os
import pandas as pd
from pathlib import Path

SCRIPTS = ["part2b_shock_comovement.py", "part3_resid_instruments.py", "part5_lp.py"]

from dotenv import load_dotenv
load_dotenv()  # Loads variables from .env into os.environ
key = os.getenv('FRED_API_KEY')

def check_prereqs():
    if not key:
        print("ERROR: FRED_API_KEY not set.")
        sys.exit(1)
    print(f"FRED_API_KEY: {key[:8]}... (set)")

    cache = Path("data/cache/bea_va_quarterly_12ind.parquet")
    if cache.exists():
        try:
            df = pd.read_parquet(cache)
            start = df["quarter_label"].min()
            end   = df["quarter_label"].max()
            nsec  = df["industry_code"].nunique()
            print(f"Industry VA cache: {start}-{end}, {nsec} sectors.")
            print("  v2 residualization sample: 2005Q2 onward (one lag).")
        except Exception as e:
            print(f"  Cache corrupted ({e}). Will re-fetch.")
            cache.unlink()
    else:
        print("No cache found -- will fetch from FRED.")

def run_script(script):
    print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
    result = subprocess.run([sys.executable, script])
    if result.returncode != 0:
        print(f"\nERROR: {script} failed (code {result.returncode})")
        sys.exit(result.returncode)
    print(f"\n{script} OK.")

def print_summary():
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    results_dir = Path("data/results")
    for fname, label in [
        ("lp_irf_delta_vacancy.csv",  "delta -> vacancies"),
        ("lp_irf_ld_vacancy.csv",     "LD    -> vacancies"),
        ("lp_irf_qu_vacancy.csv",     "QU    -> vacancies (placebo)"),
        ("lp_irf_delta_resid.csv",    "delta -> unemployment"),
        ("lp_irf_ld_resid.csv",       "LD    -> unemployment"),
        ("lp_irf_qu_resid.csv",       "QU    -> unemployment (placebo)"),
    ]:
        fpath = results_dir / fname
        if not fpath.exists():
            continue
        df = pd.read_csv(fpath)
        pk = df.loc[df["beta"].abs().idxmax()]
        sig = "***" if pk["pval"] < 0.01 else "**" if pk["pval"] < 0.05 else "*" if pk["pval"] < 0.1 else ""
        print(f"  {label}: peak={pk['beta']:+.3f} pp at h={int(pk['h'])} (p={pk['pval']:.3f}) {sig}")

    cache = Path("data/cache/bea_va_quarterly_12ind.parquet")
    if cache.exists():
        df = pd.read_parquet(cache)
        print(f"\nSpec: v2 (FRED industry VA, {df['quarter_label'].min()}-{df['quarter_label'].max()})")
        print("Key question: did r(delta,LD) drop below 0.39?")
        print("Key question: is QU vacancy IRF now flat/insignificant?")
    else:
        print("\nSpec: v1 fallback (FRED fetch failed)")

    print("\nPlots: data/results/lp_irf_vacancy_decomp_overlay.png")

if __name__ == "__main__":
    print("v2 enriched residualization pipeline")
    print("="*60)
    check_prereqs()
    for s in SCRIPTS:
        if not Path(s).exists():
            print(f"ERROR: {s} not found. Run from Data/Bartek analysis/")
            sys.exit(1)
    for s in SCRIPTS:
        run_script(s)
    print_summary()
    print("\nDone.")
