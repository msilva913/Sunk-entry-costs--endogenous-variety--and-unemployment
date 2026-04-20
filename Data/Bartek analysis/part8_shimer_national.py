"""
part8_shimer_national.py — National Shimer separation rate vs. g^LD and g^δ
============================================================================

PURPOSE
-------
Step 1 of the practical check: establish whether the Shimer-derived national
separation rate s^nat_t primarily co-moves with the LD (layoffs-and-discharges)
channel or the δ (firm destruction) channel at the aggregate level.

This is a purely diagnostic time-series exercise — it does not feed into
the LP pipeline.

WHY NOT BARTIK AGGREGATION AT NATIONAL LEVEL
---------------------------------------------
The Bartik instrument's identifying power is entirely cross-sectional: it
creates variation *across states* by interacting national industry shock rates
with different state industry compositions. At the national level there is no
cross-section to exploit.

The employment-weighted average ∑_j ω̄_j g^LD_{j,t} is by construction
essentially identical to the published national JOLTS LD rate, which BLS
already computes as an employment-weighted average of industry-level rates.
The only difference — fixed 2006 vs. current-year weights — is immaterial
for this diagnostic. The LOO correction also adds nothing here.

We therefore fetch published national rates directly from FRED:
  - g^LD_t  : FRED series JTSLDL (Layoffs & Discharges Level, thousands)
               divided by civilian employment (CE16OV) to get the rate.
               Alternatively, BLS publishes the LD rate directly.
  - g^δ_t   : national BED closing rate — NOT on FRED in usable form,
               so we still use the pipeline BED parquet and take the
               employment-weighted national average. This is fine: we are
               not imposing the Bartik cross-sectional structure, just using
               the industry-level BED data to compute one national number.

UNIT COMPARABILITY
-------------------
All three series are expressed as decimal fractions per quarter:

  s^nat_t   : Shimer E→U hazard, monthly → quarterly via 1-(1-s_m)^3
              Typical: 0.05–0.07 per quarter
              Captures ALL E→U flows (layoffs + quits→U + other)

  g^LD_t    : JOLTS LD monthly rate → quarterly via 1-(1-g_m)^3
              Typical: 0.03–0.06 per quarter
              s^nat > g^LD mechanically: quits→U and other flows in s^nat
              but not in g^LD

  g^δ_t     : BED quarterly closing rate, already quarterly
              Typical: 0.005–0.015 per quarter
              Firm-level concept; much smaller than worker-level hazards

ANALYSIS STRATEGY
-----------------
Primary: OLS regression of s^nat on g^LD and g^δ (joint and separate).
         The joint regression gives partial coefficients — the cleanest
         diagnostic — plus R² and formal tests.
         Run in demeaned levels and first-differences.

Secondary: Pairwise correlations (standard robustness check).

WHY REGRESSION OVER CORRELATION ALONE
--------------------------------------
Pairwise correlations r(s^nat, g^LD) and r(s^nat, g^δ) are contaminated
by r(g^LD, g^δ) > 0 — both are procyclical, so both correlate with s^nat
partly through their shared business cycle component.

The joint regression gives partial β coefficients: how much does s^nat
respond to g^LD holding g^δ constant, and vice versa. Additional outputs:

  β^LD / β^δ : coefficient ratio — the single most informative statistic.
               β^LD >> β^δ  → s^nat primarily tracks s-channel (LD).
               β^δ >> β^LD  → s^nat primarily tracks δ-channel (firm exit).
               β^LD ≈ β^δ   → blend; state-level Shimer rates add little.

  R²          : fraction of s^nat variation jointly explained. Low R² means
                quits→U and other flows dominate — independent dynamics that
                neither LD nor δ captures.

  Regression residual: component of s^nat not explained by either channel.
                Plotting this shows what Shimer picks up beyond LD and δ.

  F-test β^LD = β^δ: formal test of differential loading.

OUTPUTS
-------
  data/results/shimer_national_series.csv
  data/results/shimer_national_comparison.png   (time series + scatters)
  data/results/shimer_national_rolling_corr.png (rolling 8-quarter r)
  Console: regression table (demeaned + first-diff) + correlation table

PREREQUISITES
-------------
  part2_shock_rates.py   (BED parquet for g^δ)
  FRED_API_KEY in .env   (for UEMPLT5, UNEMPLOY, CLF16OV, JTSLDR)

RUN
---
  python part8_shimer_national.py
"""

import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
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

from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY", "")


def _open_file(path: Path) -> None:
    """Open a file with the default OS viewer (Windows/macOS/Linux)."""
    try:
        os.startfile(path)          # Windows
    except AttributeError:
        import subprocess, sys as _sys
        opener = "open" if _sys.platform == "darwin" else "xdg-open"
        subprocess.Popen([opener, str(path)])

# ── paths ──────────────────────────────────────────────────────────────────────
INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
CACHE_DIR   = Path("data/cache")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_START = "2001Q1"
SAMPLE_END   = "2019Q4"
PREGFC_END   = "2007Q3"

SHIMER_CACHE = CACHE_DIR / "shimer_national_monthly.parquet"
LD_CACHE     = CACHE_DIR / "jolts_ld_national_monthly.parquet"

NW_LAGS = 6   # Newey-West lag truncation for quarterly series


# ── quarter helper ─────────────────────────────────────────────────────────────
def month_to_quarter(dt: pd.Timestamp) -> str:
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"


# =============================================================================
# 1. FRED fetches
# =============================================================================
def _fetch_fred(series_ids: dict) -> pd.DataFrame:
    """
    Fetch multiple FRED series. series_ids = {fred_id: col_name}.
    Returns monthly DataFrame indexed by date.
    """
    if not FRED_API_KEY:
        sys.exit("ERROR: FRED_API_KEY not set.")
    from fredapi import Fred
    fred = Fred(api_key=FRED_API_KEY)
    frames = {}
    for sid, col in series_ids.items():
        s = fred.get_series(sid)
        s.name = col
        frames[col] = s
        print(f"    {sid}: {s.index[0].date()} – {s.index[-1].date()}")
    df = pd.DataFrame(frames)
    df.index.name = "date"
    return df.reset_index()


def fetch_shimer_inputs() -> pd.DataFrame:
    """
    UEMPLT5 : unemployed < 5 weeks, SA (thousands)
    UNEMPLOY : total unemployed, SA (thousands)
    CLF16OV  : civilian labor force, SA (thousands)
    """
    if SHIMER_CACHE.exists():
        df = pd.read_parquet(SHIMER_CACHE)
        print(f"  CPS (Shimer inputs): cached ({len(df)} months)")
        return df
    print("  Fetching CPS series from FRED ...")
    df = _fetch_fred({
        "UEMPLT5": "u_short",
        "UNEMPLOY": "u_total",
        "CLF16OV":  "lf",
    })
    df["date"] = pd.to_datetime(df["date"])
    df = df.dropna()
    # CPS redesign correction (January 1994): scale UEMPLT5 by 1.16
    # per Elsby, Michaels and Solon (2009). Our sample is post-1994 so
    # this applies uniformly and does not affect correlations.
    df.loc[df["date"] >= "1994-01-01", "u_short"] *= 1.16
    df.to_parquet(SHIMER_CACHE, index=False)
    print(f"  Cached: {SHIMER_CACHE}")
    return df


def fetch_jolts_ld() -> pd.DataFrame:
    """
    Fetch JOLTS national layoffs-and-discharges rate directly from FRED.

    We use the published BLS rate series (percent of employment per month)
    rather than constructing it via Bartik aggregation, because:
      (a) The Bartik structure adds no information at the national level.
      (b) The published rate is the exact national benchmark we want to
          compare against — it is what the Bartik construction approximates.

    FRED series used:
      JTSLDL  : Layoffs and Discharges: Total Nonfarm (Level, thousands, SA)
      CE16OV  : Civilian Employment Level (thousands, SA)

    Rate = JTSLDL / CE16OV  (decimal fraction per month)

    Alternative: BLS also publishes JTSLDR (the rate directly as percent).
    We use JTSLDL / CE16OV to construct the decimal rate, which is easier
    to convert to quarterly units consistently with s^nat.
    """
    if LD_CACHE.exists():
        df = pd.read_parquet(LD_CACHE)
        print(f"  JOLTS LD national: cached ({len(df)} months)")
        return df
    print("  Fetching JOLTS LD from FRED ...")
    raw = _fetch_fred({
        "JTSLDL": "ld_level",   # thousands of workers
        "CE16OV": "emp_level",   # thousands of employed
    })
    raw["date"] = pd.to_datetime(raw["date"])
    raw = raw.dropna()
    raw["g_ld_monthly"] = raw["ld_level"] / raw["emp_level"]  # decimal/month
    raw = raw[["date", "g_ld_monthly"]]
    raw.to_parquet(LD_CACHE, index=False)
    print(f"  Cached: {LD_CACHE}")
    print(f"  g_ld_monthly: mean={raw['g_ld_monthly'].mean():.4f}  "
          f"SD={raw['g_ld_monthly'].std():.4f}  (decimal/month)")
    return raw


# =============================================================================
# 2. Shimer (2007) national separation rate
# =============================================================================
def compute_shimer_monthly(cps: pd.DataFrame) -> pd.DataFrame:
    """
    Compute monthly f_t and s_t using Shimer (2007) continuous-time formula.

    f_t = 1 - (u_{t+1} - u^s_{t+1}) / u_t

    s_t = u^s_{t+1} * f_t / (e_t * (1 - exp(-f_t)))

    where u^s = short-term unemployed (< 5 weeks), u = total unemployed,
    e = employment = LF - u.

    Both output series are decimal fractions per month.
    """
    df = cps.sort_values("date").copy()
    df["e"] = df["lf"] - df["u_total"]
    df["u_total_lead"] = df["u_total"].shift(-1)
    df["u_short_lead"] = df["u_short"].shift(-1)
    df = df.dropna().copy()

    df["f_t"] = 1.0 - (df["u_total_lead"] - df["u_short_lead"]) / df["u_total"]
    df["f_t"] = df["f_t"].clip(0.0, 1.0)

    safe_f = df["f_t"].clip(lower=1e-6)
    df["s_t"] = (df["u_short_lead"] * df["f_t"]) / \
                (df["e"] * (1.0 - np.exp(-safe_f)))
    df["s_t"] = df["s_t"].clip(0.0, 0.20)

    df["quarter_label"] = df["date"].apply(month_to_quarter)
    print(f"  Shimer monthly: n={len(df)}  "
          f"s_t mean={df['s_t'].mean():.4f}  SD={df['s_t'].std():.4f}")
    return df[["date", "quarter_label", "f_t", "s_t"]]


def monthly_to_quarterly(monthly: pd.DataFrame,
                         rate_col: str, out_col: str) -> pd.DataFrame:
    """
    Aggregate a monthly rate to quarterly using exact compounding:
        r_quarterly = 1 - (1 - r_monthly_avg)^3

    Step 1: Average the three monthly values within each quarter.
    Step 2: Apply the compounding formula.

    This converts a monthly hazard (decimal/month) to a quarterly rate
    (decimal/quarter) in a way that is consistent with the BED quarterly
    rate g^δ, which is already expressed per quarter.

    Note: for small rates, r_q ≈ 3 * r_m, but we use exact compounding
    for precision and to avoid ambiguity when rates differ across months.
    """
    q = (monthly.groupby("quarter_label")[[rate_col]]
                .mean()
                .reset_index())
    q[out_col] = 1.0 - (1.0 - q[rate_col]) ** 3
    return q[["quarter_label", out_col]]


# =============================================================================
# 3. National BED closing rate from pipeline parquet
# =============================================================================
def load_delta_national() -> pd.DataFrame:
    """
    Load BED establishment closing rates and compute the national
    employment-weighted average across industries.

    The pipeline parquet has shape (state × industry × quarter) with column
    g_delta_loo — the leave-one-out national closing rate for each
    state × industry × quarter cell. Crucially, g_delta_loo is a NATIONAL
    rate (excluding the own state), not a state-level rate. Its value is
    the same for every state within the same industry × quarter cell, up
    to the small LOO correction. We therefore:

      Step 1: Collapse across states by taking the simple mean of g_delta_loo
              within each industry × quarter. This recovers the national
              industry-level closing rate (the LOO correction averages out
              across 50 states and becomes negligible).

      Step 2: Aggregate across industries using 2006 employment shares as
              weights to get a single national quarterly closing rate.

    Output: quarterly decimal fraction per quarter.
    The parquet stores rates as per-thousand (divide by 1000 for decimal).
    """
    candidates = [
        INSTR_DIR / "shock_rates_1992Q3_2024Q2.parquet",
        INSTR_DIR / "shock_rates_delta_base2006.parquet",
    ]
    df = None
    for p in candidates:
        if p.exists():
            df = pd.read_parquet(p)
            print(f"  BED parquet: {p}  shape={df.shape}  "
                  f"cols={df.columns.tolist()}")
            break
    if df is None:
        found = list(Path("data").rglob("*shock_rates*.parquet"))
        if not found:
            print("  WARNING: No BED parquet found. Skipping g^δ.")
            return pd.DataFrame(columns=["quarter_label", "g_delta_nat"])
        df = pd.read_parquet(found[0])
        print(f"  BED parquet (fallback): {found[0]}")

    # Identify columns
    rate_col = _find_col(df, ["g_delta_loo", "g_delta", "closing_rate",
                               "rate", "delta_rate", "exit_rate"])
    ind_col  = _find_col(df, ["industry_code", "supersector_code",
                               "supersector", "naics", "sector"])
    qt_col   = _find_col(df, ["quarter_label", "quarter"])

    if rate_col is None or qt_col is None:
        print(f"  WARNING: Cannot identify required columns. "
              f"Cols: {df.columns.tolist()}")
        return pd.DataFrame(columns=["quarter_label", "g_delta_nat"])

    if qt_col != "quarter_label":
        df = df.rename(columns={qt_col: "quarter_label"})

    # Step 1: collapse across states → national industry × quarter panel
    # g_delta_loo is a national rate so this is just averaging out the
    # small LOO correction across the 50 state rows per industry × quarter
    if ind_col is not None:
        ind_q = (df.groupby(["quarter_label", ind_col])[rate_col]
                   .mean()
                   .reset_index())
    else:
        # No industry dimension — just aggregate directly
        agg = (df.groupby("quarter_label")[rate_col]
                  .mean().reset_index()
                  .rename(columns={rate_col: "g_delta_nat"}))
        agg["g_delta_nat"] *= 1000.0   # same ×1000 as part3
        print(f"  g^δ_nat (no industry dim): {len(agg)} quarters  "
              f"mean={agg['g_delta_nat'].mean():.4f}")
        return agg[["quarter_label", "g_delta_nat"]]

    # Step 2: employment-weighted average across industries
    shares_p = INSTR_DIR / "shares_base2006.parquet"
    if shares_p.exists():
        shares  = pd.read_parquet(shares_p)
        s_ind   = _find_col(shares, [ind_col, "industry_code",
                                      "supersector_code"])
        emp_col = _find_col(shares, ["emp_share", "employment",
                                      "emp", "share"])
        if s_ind and emp_col:
            if s_ind != ind_col:
                shares = shares.rename(columns={s_ind: ind_col})
            # National share = sum of state shares per industry
            nat_wt = (shares.groupby(ind_col)[emp_col]
                            .sum().reset_index()
                            .rename(columns={emp_col: "wt"}))
            nat_wt["wt"] /= nat_wt["wt"].sum()   # normalize to sum to 1
            ind_q  = ind_q.merge(nat_wt, on=ind_col, how="left")
            ind_q["wt"] = ind_q["wt"].fillna(ind_q["wt"].mean())
            ind_q["wr"] = ind_q[rate_col] * ind_q["wt"]
            agg = (ind_q.groupby("quarter_label")["wr"]
                        .sum().reset_index()
                        .rename(columns={"wr": "g_delta_nat"}))

            # Diagnose the raw magnitude before any scaling
            raw_mean = agg["g_delta_nat"].mean()
            print(f"  g^δ raw mean after weighting: {raw_mean:.6e}")

            # The parquet stores g_delta_loo as job losses at closing
            # establishments as a fraction of *total US employment* (not
            # industry employment), making values very small (~1e-5 per
            # industry per quarter). part3 multiplies by ×1000 to bring
            # instruments to interpretable units. We apply the same here
            # so that g_delta_nat is comparable to g_ld_nat (~0.01–0.05
            # per quarter as a decimal fraction of employment).
            agg["g_delta_nat"] *= 1000.0

            print(f"  g^δ_nat (emp-weighted, ×1000): {len(agg)} quarters  "
                  f"mean={agg['g_delta_nat'].mean():.4f}  "
                  f"SD={agg['g_delta_nat'].std():.4f}  "
                  f"(decimal fraction per QUARTER after ×1000 scaling)")
            return agg[["quarter_label", "g_delta_nat"]]

    # Fallback: unweighted mean across industries
    agg = (ind_q.groupby("quarter_label")[rate_col]
                .mean().reset_index()
                .rename(columns={rate_col: "g_delta_nat"}))
    agg["g_delta_nat"] *= 1000.0   # same ×1000 as part3
    print(f"  g^δ_nat (unweighted mean): {len(agg)} quarters  "
          f"mean={agg['g_delta_nat'].mean():.4f}")
    return agg[["quarter_label", "g_delta_nat"]]


def _find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


# =============================================================================
# 4. Regression analysis (primary)
# =============================================================================
def run_regression(panel: pd.DataFrame,
                   sample_start: str, sample_end: str,
                   label: str) -> None:
    """
    Regress s^nat on g^LD and g^δ, separately and jointly.
    Run under two transformations: demeaned levels and first differences.
    Use Newey-West HAC standard errors (NW_LAGS quarters).

    Reports:
      - β^LD, β^δ (separate and joint regressions)
      - R² for each specification
      - Coefficient ratio β^LD / β^δ (joint regression)
      - F-test of joint significance
      - Implied share of s^nat variation attributable to each channel

    The regression s^nat_t = α + β^LD g^LD_t + β^δ g^δ_t + ε
    gives partial coefficients: how much does s^nat respond to each
    channel holding the other constant — the cleanest diagnostic.
    """
    sub = panel[(panel["quarter_label"] >= sample_start) &
                (panel["quarter_label"] <= sample_end)].dropna(
                subset=["s_nat", "g_ld_nat", "g_delta_nat"]).copy()

    print(f"\n{'─'*65}")
    print(f"REGRESSION RESULTS — {label}  (n={len(sub)})")
    print(f"{'─'*65}")
    print(f"{'Spec':<30} {'β^LD':>8} {'(SE)':>7} {'β^δ':>8} {'(SE)':>7}"
          f"  {'R²':>5}  {'β^LD/β^δ':>9}")
    print(f"{'─'*65}")

    for transform, tlabel in [
        ("demean",  "Demeaned levels"),
        ("fdiff",   "First differences"),
    ]:
        s  = sub["s_nat"].copy()
        gl = sub["g_ld_nat"].copy()
        gd = sub["g_delta_nat"].copy()

        if transform == "demean":
            s  = s  - s.mean()
            gl = gl - gl.mean()
            gd = gd - gd.mean()
        else:
            s  = s.diff().dropna()
            gl = gl.reindex(s.index).diff().dropna()
            gd = gd.reindex(s.index).diff().dropna()
            idx = s.index.intersection(gl.index).intersection(gd.index)
            s, gl, gd = s[idx], gl[idx], gd[idx]

        n = len(s)

        def _nw_ols(y, X_df):
            X = sm.add_constant(X_df, has_constant="add")
            m = sm.OLS(y.values, X.values).fit(
                cov_type="HAC",
                cov_kwds={"maxlags": NW_LAGS, "use_correction": True}
            )
            return m

        # Separate: s on g^LD only
        m_ld = _nw_ols(s, gl.rename("g_ld"))
        b_ld_sep = m_ld.params[1];  se_ld_sep = m_ld.bse[1]
        r2_ld    = m_ld.rsquared

        # Separate: s on g^δ only
        m_d = _nw_ols(s, gd.rename("g_d"))
        b_d_sep = m_d.params[1];  se_d_sep = m_d.bse[1]
        r2_d    = m_d.rsquared

        # Joint: s on g^LD and g^δ
        m_joint = _nw_ols(s, pd.DataFrame({"g_ld": gl, "g_d": gd}))
        b_ld_j  = m_joint.params[1];  se_ld_j = m_joint.bse[1]
        b_d_j   = m_joint.params[2];  se_d_j  = m_joint.bse[2]
        r2_j    = m_joint.rsquared
        ratio   = b_ld_j / b_d_j if abs(b_d_j) > 1e-8 else float("nan")

        print(f"\n  [{tlabel}]  n={n}")
        print(f"  {'s ~ g^LD (sep.)':<28} {b_ld_sep:>8.4f} ({se_ld_sep:.4f})"
              f"           {'':>8}          R²={r2_ld:.3f}")
        print(f"  {'s ~ g^δ (sep.)':<28}           {'':>8}  "
              f"{b_d_sep:>8.4f} ({se_d_sep:.4f})  R²={r2_d:.3f}")
        print(f"  {'s ~ g^LD + g^δ (joint)':<28} {b_ld_j:>8.4f} ({se_ld_j:.4f})"
              f"  {b_d_j:>8.4f} ({se_d_j:.4f})  R²={r2_j:.3f}  "
              f"ratio={ratio:.2f}")

        # Significance stars on joint
        def star(p): return "***" if p<.01 else "**" if p<.05 else "*" if p<.10 else ""
        p_ld = m_joint.pvalues[1];  p_d = m_joint.pvalues[2]
        print(f"  {'p-values (joint)':<28} {p_ld:>8.3f}{star(p_ld):<3}       "
              f"  {p_d:>8.3f}{star(p_d):<3}")

        # Share of s^nat variation attributable to each channel (OLS decomp)
        var_s   = np.var(s.values)
        cov_ld  = b_ld_j * np.cov(gl.values, s.values)[0, 1]
        cov_d   = b_d_j  * np.cov(gd.values, s.values)[0, 1]
        sh_ld   = cov_ld / var_s if var_s > 0 else float("nan")
        sh_d    = cov_d  / var_s if var_s > 0 else float("nan")
        sh_res  = 1.0 - sh_ld - sh_d
        print(f"  Variance decomposition (joint): "
              f"g^LD={sh_ld:.2%}  g^δ={sh_d:.2%}  residual={sh_res:.2%}")

    print()


# =============================================================================
# 5. Correlation table (secondary)
# =============================================================================
def correlation_table(panel: pd.DataFrame,
                      sample_start: str, sample_end: str,
                      label: str) -> pd.DataFrame:
    """
    Pairwise correlations under demeaned levels and first differences.
    Secondary to the regression — included for conventional reporting.
    """
    sub = panel[(panel["quarter_label"] >= sample_start) &
                (panel["quarter_label"] <= sample_end)].dropna(
                subset=["s_nat", "g_ld_nat", "g_delta_nat"]).copy()
    rows = []
    for transform, tlabel in [
        ("demean",  "Demeaned levels"),
        ("fdiff",   "First differences"),
    ]:
        s  = sub["s_nat"].copy()
        gl = sub["g_ld_nat"].copy()
        gd = sub["g_delta_nat"].copy()
        if transform == "demean":
            s -= s.mean(); gl -= gl.mean(); gd -= gd.mean()
        else:
            s = s.diff().dropna()
            gl = gl.reindex(s.index).diff().dropna()
            gd = gd.reindex(s.index).diff().dropna()
            idx = s.index.intersection(gl.index).intersection(gd.index)
            s, gl, gd = s[idx], gl[idx], gd[idx]
        rows.append({
            "sample": label, "transform": tlabel,
            "r(s,g_LD)":  round(s.corr(gl), 3),
            "r(s,g_δ)":   round(s.corr(gd), 3),
            "r(g_LD,g_δ)": round(gl.corr(gd), 3),
            "n": len(s),
        })
    return pd.DataFrame(rows)


# =============================================================================
# 6. Plots
# =============================================================================
def plot_comparison(panel: pd.DataFrame, out_path: Path) -> None:
    plt.ioff()
    sub = panel[(panel["quarter_label"] >= SAMPLE_START) &
                (panel["quarter_label"] <= SAMPLE_END)].copy()

    def zs(x): return (x - x.mean()) / x.std()
    sub = sub.copy()
    sub["s_z"]  = zs(sub["s_nat"])
    sub["ld_z"] = zs(sub["g_ld_nat"])
    sub["d_z"]  = zs(sub["g_delta_nat"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Time series
    ax = axes[0]
    x  = range(len(sub))
    xt = [i for i, q in enumerate(sub["quarter_label"])
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]
    xl = [q[:4] for i, q in enumerate(sub["quarter_label"])
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]
    _ = ax.plot(x, sub["s_z"],  "#1f77b4", lw=1.8,
                label=r"$s^{nat}_t$ (Shimer)")
    _ = ax.plot(x, sub["ld_z"], "#d62728", lw=1.6, ls="--",
                label=r"$g^{LD}_t$ (JOLTS — published national rate)")
    _ = ax.plot(x, sub["d_z"],  "#2ca02c", lw=1.6, ls=":",
                label=r"$g^{\delta}_t$ (BED closing)")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8)
    ax.set_title("Z-scored quarterly series\n2001Q1–2019Q4", fontsize=10)
    ax.set_ylabel("SD from mean", fontsize=9)
    ax.legend(fontsize=7); ax.grid(axis="y", lw=0.4, alpha=0.4)

    # Scatter: s vs g^LD
    for ax, col, color, lbl in [
        (axes[1], "ld_z", "#d62728", r"$g^{LD}_t$"),
        (axes[2], "d_z",  "#2ca02c", r"$g^{\delta}_t$"),
    ]:
        r = sub["s_z"].corr(sub[col])
        _ = ax.scatter(sub[col], sub["s_z"], alpha=0.45, s=18,
                       color=color, edgecolors="none")
        xl_ = np.linspace(sub[col].min(), sub[col].max(), 50)
        m, b = np.polyfit(sub[col].dropna(), sub["s_z"].dropna(), 1)
        _ = ax.plot(xl_, m*xl_+b, color=color, lw=1.5)
        _ = ax.set_title(fr"$s^{{nat}}$ vs {lbl}" + f"\nr = {r:.3f}",
                         fontsize=10)
        _ = ax.set_xlabel(f"{lbl} (z-score)", fontsize=9)
        _ = ax.set_ylabel(r"$s^{nat}$ (z-score)", fontsize=9)
        ax.grid(lw=0.4, alpha=0.3)

    fig.suptitle(
        r"National Shimer $s^{nat}_t$ vs. JOLTS LD and BED closing rates"
        "\n(using published national JOLTS rate, not Bartik aggregation)",
        fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


def plot_rolling_corr(panel: pd.DataFrame, out_path: Path,
                      window: int = 8) -> None:
    plt.ioff()
    sub = panel[(panel["quarter_label"] >= SAMPLE_START) &
                (panel["quarter_label"] <= SAMPLE_END)].copy().reset_index(drop=True)
    r_ld = sub["s_nat"].rolling(window).corr(sub["g_ld_nat"])
    r_d  = sub["s_nat"].rolling(window).corr(sub["g_delta_nat"])

    fig, ax = plt.subplots(figsize=(11, 4))
    x = range(len(sub))
    _ = ax.plot(x, r_ld, "#d62728", lw=1.8, label=r"$r(s^{nat},\ g^{LD})$")
    _ = ax.plot(x, r_d,  "#2ca02c", lw=1.8, ls="--",
                label=r"$r(s^{nat},\ g^{\delta})$")
    ax.axhline(0, color="black", lw=0.6)
    ax.axhline(0.5, color="grey", lw=0.4, ls=":")
    gfc_idx = next((i for i, q in enumerate(sub["quarter_label"])
                    if q == "2008Q3"), None)
    if gfc_idx:
        ax.axvline(gfc_idx, color="orange", lw=1.2, ls="--",
                   alpha=0.7, label="GFC start (2008Q3)")
    xt = [i for i, q in enumerate(sub["quarter_label"])
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]
    xl = [q[:4] for q in sub["quarter_label"]
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8)
    _ = ax.set_title(f"Rolling {window}-quarter correlation with "
                     r"$s^{nat}_t$", fontsize=11)
    _ = ax.set_ylabel("Correlation coefficient", fontsize=9)
    _ = ax.legend(fontsize=9); ax.grid(lw=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print("=" * 65)
    print("Part 8 — National Shimer separation rate diagnostic")
    print("=" * 65)

    # ── 1. Shimer rate ────────────────────────────────────────────────────────
    print("\n[1] Computing national Shimer separation rate")
    cps     = fetch_shimer_inputs()
    monthly = compute_shimer_monthly(cps)
    s_q     = monthly_to_quarterly(monthly, "s_t", "s_nat")

    # ── 2. JOLTS LD national rate (published, not Bartik) ─────────────────────
    print("\n[2] Fetching published national JOLTS LD rate from FRED")
    print("    (using JTSLDL / CE16OV — not Bartik aggregation)")
    ld_m  = fetch_jolts_ld()
    ld_m["quarter_label"] = ld_m["date"].apply(month_to_quarter)
    g_ld  = monthly_to_quarterly(ld_m, "g_ld_monthly", "g_ld_nat")

    # ── 3. BED national closing rate ──────────────────────────────────────────
    print("\n[3] Loading national BED closing rate from pipeline parquet")
    print("    (employment-weighted average of industry rates — "
          "not Bartik cross-sectional structure)")
    g_delta = load_delta_national()

    # ── 4. Merge ──────────────────────────────────────────────────────────────
    print("\n[4] Merging and restricting to estimation sample")
    panel = (s_q
             .merge(g_ld,    on="quarter_label", how="outer")
             .merge(g_delta, on="quarter_label", how="outer")
             .sort_values("quarter_label").reset_index(drop=True))

    print("\n  Unit summary (all in decimal fraction per QUARTER):")
    for col, desc in [
        ("s_nat",      "Shimer E→U rate [all separations→U]"),
        ("g_ld_nat",   "JOLTS LD rate   [layoffs→U only]"),
        ("g_delta_nat","BED closing rate[firm exit, quarterly]"),
    ]:
        sub = panel[(panel["quarter_label"] >= SAMPLE_START) &
                    (panel["quarter_label"] <= SAMPLE_END)]
        s = sub[col].dropna()
        if len(s):
            print(f"    {col:<14}: mean={s.mean():.4f}  "
                  f"SD={s.std():.4f}  n={len(s)}")

    # Save
    out_csv = RESULTS_DIR / "shimer_national_series.csv"
    panel.to_csv(out_csv, index=False)
    print(f"\n  Saved: {out_csv}")

    # ── 5. Regression (primary) ───────────────────────────────────────────────
    print("\n[5] Regression analysis (primary diagnostic)")
    print("    s^nat_t = α + β^LD·g^LD_t + β^δ·g^δ_t + ε")
    print("    Newey-West HAC SEs, lag truncation =", NW_LAGS)
    run_regression(panel, SAMPLE_START, SAMPLE_END,
                   f"Full sample {SAMPLE_START}–{SAMPLE_END}")
    run_regression(panel, SAMPLE_START, PREGFC_END,
                   f"Pre-GFC {SAMPLE_START}–{PREGFC_END}")

    # ── 6. Correlation table (secondary) ─────────────────────────────────────
    print("\n[6] Pairwise correlations (secondary)")
    c1 = correlation_table(panel, SAMPLE_START, SAMPLE_END,
                           f"Full {SAMPLE_START}–{SAMPLE_END}")
    c2 = correlation_table(panel, SAMPLE_START, PREGFC_END,
                           f"Pre-GFC {SAMPLE_START}–{PREGFC_END}")
    print(pd.concat([c1, c2]).to_string(index=False))

    # ── 7. Plots ──────────────────────────────────────────────────────────────
    print("\n[7] Plotting")
    plot_comparison(panel,
                    RESULTS_DIR / "shimer_national_comparison.png")
    plot_rolling_corr(panel,
                      RESULTS_DIR / "shimer_national_rolling_corr.png")

    print("\nDone.")