"""
part2b_extend_va.py  --  Constrained Chow-Lin backcast of sector VA pre-2005
=============================================================================
Extends the quarterly BEA sector value-added series from 2005Q1 back to
1992Q1 using a constrained Chow-Lin (1971, JASA) temporal disaggregation.

The resulting extended cache (bea_va_quarterly_12ind_extended.parquet) replaces
bea_va_quarterly_12ind.parquet as the VA control series in part2b_residualize_
shocks.py, expanding the v2/v3 residualization sample from 2005Q2+ to 1992Q2+.

BLS-to-FRED mapping corrections vs. prior version
--------------------------------------------------
  BLS 65: RVAESHS    (Education+Health+Social Assistance, full aggregate)
          Previously used RVAES (Education only, ~16% of total) — corrected.
  BLS 70: RVAAER + RVAAF   (Arts+Entertainment + Accommodation+Food)
          Previously used RVAER (does not exist on FRED) — corrected.

Constrained Chow-Lin formulation
---------------------------------
Let VA_{j,t} be quarterly sector VA (unobserved pre-2005), and Y_{j,A} be
the BEA annual benchmark for year A (observed 1997+).

Step 1 — Indicator regression (post-2005 quarterly data):
    Δlog(VA_{j,t}) = α_j + β_j · Δlog(GDP_t) + ε_{j,t}

    OLS with intercept. α_j captures sector-specific mean growth that differs
    from GDP (secular tech growth, mining decline, etc.); β_j is the short-run
    cyclical GDP elasticity. Omitting the intercept produces negative centered-R²
    for sectors where mean growth ≠ mean GDP growth. R² is reported for each
    sector; sectors with R² < R2_THRESHOLD receive β_use = 0 (intercept-only
    backcast) to avoid adding noisy GDP variation.

Step 2 — Unconstrained recursive backcast:
    For each pre-2005 quarter, apply implied growth = α_j + β_use·Δlog(GDP_t)
    recursively backward from the 2005Q1 anchor level.

Step 3 — Annual benchmark constraint (proportional Denton adjustment):
    BEA annual VA = average of 4 quarterly SAAR values. For each year A in
    1997–2004 with annual benchmark Y_{j,A}:

        scale_A = Y_{j,A} / mean(VA_{j,q}^unconstrained, q in {Q1..Q4 of A})
        VA_{j,q}^constrained = VA_{j,q}^unconstrained · scale_A

    This guarantees mean(VA_{j,q}^constrained) = Y_{j,A} exactly, making the
    quarterly path consistent with published BEA annual data.
    Years 1992–1996 (no annual benchmark available) use the unconstrained
    backcast only.

Annual benchmark source:
    bea_va_annual_by_supersector.csv — produced from BEA GDP-by-Industry
    Table 1.3.6, April 9 2026 release. One column per BLS supersector (BLS10
    through BLS80), row index = year, units = billions of chained 2017 dollars.
    File must be co-located with this script (or edit ANN_BENCH_CSV below).

GDP indicator:
    FRED GDPC1 (Real GDP, chained 2017 dollars). Cached after first fetch.

Outputs:
    data/cache/bea_va_quarterly_12ind_extended.parquet
    data/results/va_chowlin_backcast_validation.png

Run:
    FRED_API_KEY=<key> python part2b_extend_va.py

Prerequisites:
    part2b_residualize_shocks.py  (produces bea_va_quarterly_12ind.parquet)
    bea_va_annual_by_supersector.csv  (annual benchmark CSV)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── working directory ──────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import DEFAULT_CACHE_DIR, INDUSTRY_LABELS

from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY")

# ── paths ──────────────────────────────────────────────────────────────────────
BEA_QTR_CACHE = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet"
BEA_EXT_CACHE = DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind_extended.parquet"
GDP_CACHE     = DEFAULT_CACHE_DIR / "gdpc1_quarterly.parquet"
RESULTS_DIR   = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Annual benchmark CSV: must be in same directory as this script.
ANN_BENCH_CSV = Path("bea_va_annual_by_supersector.csv")

# ── constants ──────────────────────────────────────────────────────────────────
TARGET_START  = "1992Q1"   # desired start of extended quarterly series
SPLICE_QT     = "2005Q1"   # first observed FRED quarterly value (anchor)
R2_THRESHOLD  = 0.15       # R² below this → suppress β (intercept-only backcast)

# ── FRED quarterly series for each BLS supersector ────────────────────────────
# Used to estimate α_j and β_j from post-2005 data.
# Corrections vs. prior version:
#   BLS 65: RVAESHS replaces RVAES (education-only subsector, ~16% of total).
#   BLS 70: RVAAER + RVAAF replaces RVAER (non-existent) + RVAAF.
BLS_TO_FRED = {
    10: ["RVAM"],            # Mining
    20: ["RVAC"],            # Construction
    30: ["RVAMA"],           # Manufacturing
    41: ["RVAW"],            # Wholesale Trade
    42: ["RVAR"],            # Retail Trade
    43: ["RVAT", "RVAU"],    # Transport+Warehousing + Utilities
    50: ["RVAI"],            # Information
    55: ["RVAFI", "RVARL"],  # Finance+Insurance + Real Estate+Rental
    60: ["RVAPBS"],          # Professional+Business Services
    65: ["RVAESHS"],         # Education+Health+Social Assistance (full aggregate)
    70: ["RVAAER", "RVAAF"], # Arts+Entertainment + Accommodation+Food
    80: ["RVAOSEG"],         # Other Services excl Govt
}


# ── helpers ────────────────────────────────────────────────────────────────────
def _dt_to_ql(dt) -> str:
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"

def _ql_to_dt(ql: str) -> pd.Timestamp:
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

def _get_fred_client():
    if not FRED_API_KEY:
        raise EnvironmentError(
            "FRED_API_KEY not set.\n"
            "  Register free: https://fred.stlouisfed.org/docs/api/api_key.html\n"
            "  Windows PowerShell:  $env:FRED_API_KEY = '<key>'\n"
            "  Linux/macOS:         export FRED_API_KEY=<key>"
        )
    from fredapi import Fred
    return Fred(api_key=FRED_API_KEY)

def _sector_label(bls_code: int) -> str:
    return (INDUSTRY_LABELS.get(bls_code)
            or INDUSTRY_LABELS.get(str(bls_code), "?"))


# ── Section 1: GDP indicator ───────────────────────────────────────────────────
def _load_gdp() -> pd.Series:
    """Real GDP (GDPC1) as quarterly Series indexed by quarter_label."""
    if GDP_CACHE.exists():
        df = pd.read_parquet(GDP_CACHE)
        print(f"  GDP: cache ({df['quarter_label'].min()}–{df['quarter_label'].max()})")
    else:
        print("  GDP: fetching GDPC1 from FRED ...")
        fred = _get_fred_client()
        s    = fred.get_series("GDPC1", observation_start="1990-01-01")
        df   = (pd.DataFrame({"date": s.index, "gdp": s.values})
                  .assign(quarter_label=lambda d: d["date"].map(_dt_to_ql))
                  [["quarter_label", "gdp"]].dropna())
        df.to_parquet(GDP_CACHE, index=False)
        print(f"  GDP: saved ({df['quarter_label'].min()}–{df['quarter_label'].max()})")
    return df.set_index("quarter_label")["gdp"].sort_index()


# ── Section 2: Post-2005 quarterly VA levels ──────────────────────────────────
def _fetch_quarterly_va() -> pd.DataFrame:
    """
    Fetch post-2005 quarterly VA levels from FRED (2004Q4+ to get the lag).
    Returns long DataFrame: quarter_label, industry_code (int), va.
    """
    print("  Fetching quarterly VA from FRED ...")
    fred    = _get_fred_client()
    records = []
    for bls_code, sids in BLS_TO_FRED.items():
        components = {}
        for sid in sids:
            try:
                s = fred.get_series(sid, observation_start="2004-10-01")
                components[sid] = s.resample("QS").first()
                print(f"    {sid}: {len(components[sid])} quarters")
            except Exception as e:
                sys.exit(f"[ERROR] FRED series '{sid}' failed: {e}\n"
                         f"  Check BLS_TO_FRED mapping above.")
        va_total = pd.concat(list(components.values()), axis=1).sum(
            axis=1, min_count=len(sids)
        )
        df = (pd.DataFrame({"date": va_total.index, "va": va_total.values})
                .dropna()
                .assign(quarter_label=lambda d: d["date"].map(_dt_to_ql),
                        industry_code=int(bls_code))
                [["quarter_label", "industry_code", "va"]])
        records.append(df)

    result = (pd.concat(records, ignore_index=True)
                .sort_values(["industry_code", "quarter_label"])
                .reset_index(drop=True))
    print(f"  VA levels: {result['quarter_label'].min()}–{result['quarter_label'].max()}, "
          f"{result['industry_code'].nunique()} sectors")
    return result


# ── Section 3: Annual benchmarks ──────────────────────────────────────────────
def _load_annual_benchmarks() -> pd.DataFrame:
    """
    Load BEA annual VA benchmarks from bea_va_annual_by_supersector.csv.

    CSV schema:
      - Index column: 'year' (int, 1997–2025)
      - Value columns: BLS10, BLS20, ..., BLS80
      - Units: billions of chained 2017 dollars (same as FRED SAAR series)

    BEA convention: annual value = mean of 4 quarterly SAAR values.
    Annual constraint in Step 3 therefore requires:
        mean(VA_{Q1..Q4}^constrained) = annual_benchmark

    Coverage 1997–2004 is used for constraints; 1992–1996 has no benchmark.
    """
    if not ANN_BENCH_CSV.exists():
        sys.exit(
            f"[ERROR] Annual benchmark CSV not found: {ANN_BENCH_CSV.resolve()}\n"
            f"  Produce it from BEA GDP-by-Industry Table 1.3.6 using the\n"
            f"  companion Excel workbook, then save as:\n"
            f"  {ANN_BENCH_CSV.resolve()}"
        )
    df = pd.read_csv(ANN_BENCH_CSV, index_col="year")
    df.columns = [int(c.replace("BLS", "")) for c in df.columns]
    df.index   = df.index.astype(int)
    print(f"  Annual benchmarks: {df.index.min()}–{df.index.max()}, "
          f"{len(df.columns)} sectors")
    return df   # shape (n_years, 12), index=year, columns=bls_code (int)


# ── Section 4: Constrained Chow-Lin backcast (one sector) ─────────────────────
def _chow_lin_sector(
    va_quarterly: pd.Series,  # quarterly VA levels post-2005Q1, index=quarter_label
    gdp         : pd.Series,  # real GDP quarterly, index=quarter_label
    ann_bench   : pd.Series,  # annual benchmarks, index=year (int)
    bls_code    : int,
) -> tuple:
    """
    Constrained Chow-Lin for one BLS supersector.

    Returns:
        dlog_full   : pd.Series of Δlog(VA), quarter_label index, TARGET_START onward
        meta        : dict with alpha, beta, beta_use, r2, n_constrained
    """
    # Align quarterly VA and GDP on common post-2005 quarters
    common  = va_quarterly.index.intersection(gdp.index)
    va_obs  = va_quarterly.loc[common].dropna()
    gdp_obs = gdp.loc[va_obs.index.intersection(gdp.index)]
    va_obs  = va_obs.loc[gdp_obs.index]

    if len(va_obs) < 8:
        sys.exit(f"[ERROR] BLS {bls_code}: only {len(va_obs)} common observations. "
                 "Check FRED series.")

    # ── Step 1: OLS regression Δlog(VA) = α + β·Δlog(GDP) + ε ───────────────
    dlog_va  = np.log(va_obs).diff().dropna()
    dlog_gdp = np.log(gdp_obs).diff().reindex(dlog_va.index).dropna()
    dlog_va  = dlog_va.reindex(dlog_gdp.index)

    X      = np.column_stack([np.ones(len(dlog_gdp)), dlog_gdp.values])
    y      = dlog_va.values
    coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
    alpha, beta = float(coeffs[0]), float(coeffs[1])

    resid  = y - X @ coeffs
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2     = float(1.0 - (resid**2).sum() / ss_tot) if ss_tot > 0 else 0.0

    # Suppress β for low-R² sectors to avoid adding noisy GDP variation
    beta_use = beta if r2 >= R2_THRESHOLD else 0.0
    suppressed = r2 < R2_THRESHOLD

    lbl = f"{_sector_label(bls_code):<35}"
    print(f"    BLS {bls_code:<3} ({lbl}): "
          f"α={alpha:+.4f}  β={beta:+.5f}  R²={r2:.3f}  N={len(dlog_va)}"
          + ("  [β suppressed]" if suppressed else ""))

    # ── Step 2: unconstrained recursive backcast ──────────────────────────────
    anchor_qt    = va_obs.index[0]   # ≈ 2005Q1
    backcast_qts = sorted(q for q in gdp.index if TARGET_START <= q < anchor_qt)

    if not backcast_qts:
        dlog_obs = np.log(va_quarterly).diff().dropna()
        meta = dict(alpha=alpha, beta=beta, beta_use=beta_use, r2=r2, n_constrained=0)
        return dlog_obs, meta

    all_qts = sorted(set(backcast_qts + [anchor_qt]))
    log_va  = {anchor_qt: float(np.log(va_obs.iloc[0]))}

    for i in range(len(all_qts) - 1, 0, -1):
        q_next = all_qts[i]
        q_curr = all_qts[i - 1]
        if q_curr not in gdp.index or q_next not in gdp.index:
            continue
        g_curr, g_next = float(gdp[q_curr]), float(gdp[q_next])
        if g_curr <= 0 or g_next <= 0:
            continue
        dlog_gdp_q     = np.log(g_next) - np.log(g_curr)
        dlog_va_q      = alpha + beta_use * dlog_gdp_q
        log_va[q_curr] = log_va[q_next] - dlog_va_q

    va_backcast = pd.Series(
        {q: np.exp(v) for q, v in log_va.items() if q < anchor_qt},
        name="va"
    ).sort_index()

    # ── Step 3: annual benchmark constraint (proportional Denton) ─────────────
    # For year A in backcast window with benchmark Y_{j,A}:
    #   scale_A = Y_{j,A} / mean(VA_backcast[Q1:Q4 of A])
    #   VA_backcast[Q1:Q4 of A] *= scale_A
    # Guarantees mean(VA_constrained[Q1:Q4]) = Y_{j,A} exactly.
    n_constrained = 0
    bench_years   = set(ann_bench.dropna().index) if len(ann_bench) > 0 else set()
    for yr in sorted(set(int(q[:4]) for q in va_backcast.index)):
        if yr not in bench_years:
            continue
        yr_qts = [q for q in va_backcast.index if q.startswith(f"{yr}Q")]
        if len(yr_qts) < 4:
            continue
        predicted_mean = float(va_backcast[yr_qts].mean())
        if predicted_mean <= 0:
            continue
        scale = float(ann_bench[yr]) / predicted_mean
        va_backcast[yr_qts] = va_backcast[yr_qts] * scale
        n_constrained += 1

    # ── Step 4: splice with observed quarterly series and compute Δlog ─────────
    full_va   = pd.concat([va_backcast, va_quarterly]).sort_index()
    full_va   = full_va[~full_va.index.duplicated(keep="last")]
    dlog_full = np.log(full_va).diff().dropna()

    meta = dict(alpha=alpha, beta=beta, beta_use=beta_use,
                r2=r2, n_constrained=n_constrained)
    return dlog_full, meta


# ── Section 5: Validation plot ─────────────────────────────────────────────────
def _plot_validation(extended: pd.DataFrame, gdp: pd.Series,
                     results_meta: dict) -> None:
    """
    Per-sector plot: observed (blue) vs backcasted (red) quarterly Δlog(VA),
    with GDP×β reference (grey dashed). Marks splice point and NBER recessions.
    """
    sectors  = sorted(extended["industry_code"].unique())
    all_qts  = sorted(extended["quarter_label"].unique())
    dates    = [_ql_to_dt(q) for q in all_qts]
    dlog_gdp = np.log(gdp).diff()
    rec_spans = [("2001-03-01","2001-11-01"), ("2007-12-01","2009-06-01")]

    fig, axes = plt.subplots(len(sectors), 1,
                             figsize=(14, 2.8 * len(sectors)), sharex=True)
    if len(sectors) == 1:
        axes = [axes]

    for ax, code in zip(axes, sectors):
        sub  = (extended[extended["industry_code"] == code]
                .set_index("quarter_label")["dlog_va"]
                .reindex(all_qts))
        dlog = sub.values
        is_back = np.array([q < SPLICE_QT for q in all_qts])

        ax.plot(dates, np.where(~is_back, dlog, np.nan),
                color="#1f77b4", lw=1.2, label="Observed (FRED)")
        ax.plot(dates, np.where( is_back, dlog, np.nan),
                color="#d62728", lw=1.2, label="Chow-Lin backcast")

        m = results_meta.get(code, {})
        if m.get("beta_use", 0) != 0:
            ref = (dlog_gdp.reindex(all_qts) * m["beta_use"] + m["alpha"]).values
            ax.plot(dates, ref, color="grey", lw=0.8, ls="--", alpha=0.5,
                    label=f"α+β·ΔGDP  (β={m['beta_use']:.2f})")

        for s, e in rec_spans:
            ax.axvspan(pd.Timestamp(s), pd.Timestamp(e), alpha=0.10, color="grey")
        ax.axvline(_ql_to_dt(SPLICE_QT), color="black", lw=0.7, ls=":")
        ax.axhline(0, color="black", lw=0.4)

        r2  = m.get("r2", float("nan"))
        nc  = m.get("n_constrained", 0)
        sup = "  β=0" if m.get("beta_use", 1) == 0 else ""
        ax.set_ylabel(f"BLS {code} | R²={r2:.2f}{sup} | {nc} yrs constrained",
                      fontsize=7.5)
        ax.grid(axis="y", lw=0.4, alpha=0.4)
        if ax is axes[0]:
            ax.legend(fontsize=7, ncol=3, loc="upper left")

    axes[-1].set_xlabel("Year", fontsize=9)
    fig.suptitle(
        "Constrained Chow-Lin: quarterly Δlog(sector VA)\n"
        "Red = backcast | Blue = observed FRED | "
        "Dotted = 2005Q1 splice | Annual constraints applied 1997–2004",
        fontsize=10
    )
    fig.tight_layout()
    out = RESULTS_DIR / "va_chowlin_backcast_validation.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 70)
print("part2b_extend_va.py  --  Constrained Chow-Lin sector VA backcast")
print("=" * 70)

if not FRED_API_KEY:
    sys.exit(
        "ERROR: FRED_API_KEY not set.\n"
        "  Register free at https://fred.stlouisfed.org/docs/api/api_key.html"
    )

# [1] GDP indicator
print("\n[1] Real GDP indicator")
gdp = _load_gdp()

# [2] Post-2005 quarterly VA levels (estimation sample for α, β)
print("\n[2] Post-2005 quarterly sector VA levels")
va_levels = _fetch_quarterly_va()

# [3] Annual benchmarks (1997–2025; 1997–2004 used as constraints)
print("\n[3] Annual VA benchmarks")
ann_bench_df = _load_annual_benchmarks()

# [4] Constrained Chow-Lin for each sector
print(f"\n[4] Constrained Chow-Lin backcast  (target: {TARGET_START}, "
      f"R² threshold: {R2_THRESHOLD})")
print(f"    {'Sector':<42}  α         β        R²")

records      = []
results_meta = {}

for bls_code in sorted(va_levels["industry_code"].unique()):
    code_int = int(bls_code)
    va_sub   = (va_levels[va_levels["industry_code"] == bls_code]
                .set_index("quarter_label")["va"]
                .sort_index())
    ann_sub  = (ann_bench_df[code_int]
                if code_int in ann_bench_df.columns
                else pd.Series(dtype=float))

    dlog_full, meta = _chow_lin_sector(va_sub, gdp, ann_sub, code_int)

    df_out = dlog_full.reset_index()
    df_out.columns = ["quarter_label", "dlog_va"]
    df_out["industry_code"] = code_int
    records.append(df_out[["quarter_label", "industry_code", "dlog_va"]])
    results_meta[code_int] = meta

extended = (pd.concat(records, ignore_index=True)
              .sort_values(["industry_code", "quarter_label"])
              .reset_index(drop=True))

n_pre  = extended[extended["quarter_label"] <  SPLICE_QT]["quarter_label"].nunique()
n_post = extended[extended["quarter_label"] >= SPLICE_QT]["quarter_label"].nunique()
print(f"\n  Extended: {extended['quarter_label'].min()}–{extended['quarter_label'].max()}, "
      f"{extended['industry_code'].nunique()} sectors, {len(extended):,} rows")
print(f"  Backcasted quarters: {n_pre}  |  Observed (FRED): {n_post}")

# [5] Summary table
print("\n[5] Estimation summary")
print()
HDR = (f"  {'BLS':<5}  {'Sector / description':<35}  "
       f"{'alpha':>9}  {'beta_OLS':>9}  {'beta_use':>9}  "
       f"{'R-sq':>7}  {'Ann.':>5}")
SEP = ("  " + "-"*5 + "  " + "-"*35 + "  " + "-"*9 + "  " + "-"*9 +
       "  " + "-"*9 + "  " + "-"*7 + "  " + "-"*5)
print(HDR)
print(SEP)
for code in sorted(results_meta.keys()):
    m    = results_meta[code]
    lbl  = _sector_label(code)
    flag = " *" if m["beta_use"] == 0.0 else ""
    print(f"  {code:<5}  {lbl:<35}  "
          f"{m['alpha']:>+9.4f}  {m['beta']:>+9.4f}  {m['beta_use']:>+9.4f}  "
          f"{m['r2']:>7.3f}  {m['n_constrained']:>4}{flag}")
print(SEP)
print()
print("  Column guide:")
print("    alpha     : sector-specific mean quarterly growth (intercept)")
print("    beta_OLS  : OLS estimate of GDP elasticity (post-2005 data)")
print("    beta_used : beta applied in backcast (0 if R-sq below threshold)")
print("    R-sq      : in-sample R² of Δlog(VA) ~ α + β·Δlog(GDP)")
print("    Ann.      : years with annual BEA benchmark applied as constraint")
print(f"  * beta_used suppressed to 0 (R-sq < {R2_THRESHOLD}); backcast uses α only.")

# [6] Save
print(f"\n[6] Saving {BEA_EXT_CACHE}")
extended.to_parquet(BEA_EXT_CACHE, index=False)
print(f"  Saved. part2b_residualize_shocks.py will use this automatically.")

# [7] Validation plot
print("\n[7] Validation plot")
_plot_validation(extended, gdp, results_meta)

print("\nDone.")
