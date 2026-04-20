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
SAMPLE_END   = "2019Q4"   # excludes COVID — see documentation below
PREGFC_END   = "2007Q3"

# GFC episode window only — COVID excluded (see WHY NOT COVID below)
GFC_START    = "2008Q3"
GFC_END      = "2009Q4"

# WHY THE SAMPLE ENDS AT 2019Q4 — COVID MEASUREMENT BREAKDOWN
# ─────────────────────────────────────────────────────────────
# The Shimer (2007) separation rate formula requires:
#   s_t = u^s_{t+1} * f_t / (e_t * (1 - exp(-f_t)))
# where u^s_{t+1} = CPS short-term unemployed (< 5 weeks).
#
# During 2020Q1–2021Q4 this formula breaks down due to a well-documented
# CPS classification error. BLS survey interviewers classified millions of
# temporarily laid-off workers as "employed but absent from work" rather
# than unemployed (Bick and Blandin 2020; Abraham and Kearney 2020).
# This caused:
#   (1) u_total to be severely understated — the denominator of f_t
#       collapses, causing f_t to be mechanically mis-estimated
#   (2) u^s_{t+1} to be understated — reducing the numerator of s_t
#   (3) Paradoxical result: in 2020Q2, when g_delta spiked to 2.4%
#       (vs normal 1.1%) and g_LD spiked to 9.7% (vs normal 4%),
#       s_nat *fell* from 15.3% to 9.4% — the opposite of economics.
#
# The COVID interactions in a full-sample regression absorb this
# measurement artifact rather than identifying structural economics.
# Including them adds noise to the normal-times β estimates without
# adding interpretable information. The correct treatment is sample
# exclusion, not dummy-variable correction, because the misclassification
# is correlated with both g_LD (which spiked) and g_delta (which spiked)
# in a way that the interaction terms cannot cleanly separate from
# genuine structural amplification.
#
# The BED parquet extends to 2024Q2, and JOLTS/CPS extend further, so
# the data *exists* — the exclusion is a methodological choice, not a
# data constraint.

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
    Regress s^nat on g^LD and g^delta with GFC episode dummy and interactions.

    Specification (primary — spec C)
    ---------------------------------
    s^nat_t = alpha
              + beta^LD  g^LD_t + beta^delta g^delta_t  <- normal-times
              + gamma  D^GFC_t                          <- GFC level shift
              + delta^LD  g^LD_t x D^GFC_t              <- GFC amplification
              + delta^delta  g^delta_t x D^GFC_t
              + eps_t

    All continuous variables demeaned. beta^LD and beta^delta identify
    the structural relationship in non-GFC quarters at average conditions.

    GFC window: GFC_START - GFC_END (NBER peak-trough). n=6 quarters.

    Why only GFC, not COVID
    -----------------------
    The COVID period (2020Q1-2021Q4) is excluded from the sample entirely
    because the Shimer formula breaks down due to CPS misclassification:
    millions of temporarily laid-off workers were coded as "employed but
    absent," causing s_nat to fall in 2020Q2 precisely when g_delta and
    g_LD spiked (see SAMPLE_END documentation in constants block above).
    Dummy-variable correction cannot separate this measurement artifact
    from genuine structural amplification because both are correlated with
    the same episode indicators.

    Reports three nested specifications for stability assessment:
      [A] Baseline — no episode controls
      [B] GFC level dummy only
      [C] GFC dummy + interactions (primary)

    Newey-West HAC standard errors, lag = NW_LAGS quarters.
    """
    sub = (panel[(panel["quarter_label"] >= sample_start) &
                 (panel["quarter_label"] <= sample_end)]
           .dropna(subset=["s_nat", "g_ld_nat", "g_delta_nat"])
           .copy()
           .reset_index(drop=True))

    sub["d_gfc"] = ((sub["quarter_label"] >= GFC_START) &
                    (sub["quarter_label"] <= GFC_END)).astype(float)
    n_gfc  = int(sub["d_gfc"].sum())
    n_norm = len(sub) - n_gfc

    s  = sub["s_nat"]       - sub["s_nat"].mean()
    gl = sub["g_ld_nat"]    - sub["g_ld_nat"].mean()
    gd = sub["g_delta_nat"] - sub["g_delta_nat"].mean()
    gl_gfc = gl * sub["d_gfc"]
    gd_gfc = gd * sub["d_gfc"]

    def _nw_ols(y, X_df):
        X = sm.add_constant(X_df, has_constant="add")
        return sm.OLS(y.values, X.values).fit(
            cov_type="HAC",
            cov_kwds={"maxlags": NW_LAGS, "use_correction": True}
        )

    def star(p):
        return "***" if p < .01 else "**" if p < .05 else "*" if p < .10 else ""

    print(f"\n{chr(0x2550)*68}")
    print(f"REGRESSION — {label}")
    print(f"  n={len(sub)}  normal={n_norm}  GFC={n_gfc}"
          f"  |  GFC window: {GFC_START}-{GFC_END}")
    print(f"  NW HAC SEs lag={NW_LAGS}")
    print(chr(0x2550)*68)

    # [A] Baseline
    m_a = _nw_ols(s, pd.DataFrame({"g_ld": gl, "g_d": gd}))
    b_ld_a = m_a.params[1]; b_d_a = m_a.params[2]
    ratio_a = b_ld_a / b_d_a if abs(b_d_a) > 1e-8 else float("nan")
    print(f"\n  [A] Baseline — no episode controls  R2={m_a.rsquared:.3f}")
    print(f"    {'beta^LD':<28} {b_ld_a:>9.4f}  ({m_a.bse[1]:.4f})"
          f"  p={m_a.pvalues[1]:.3f}{star(m_a.pvalues[1])}")
    print(f"    {'beta^delta':<28} {b_d_a:>9.4f}  ({m_a.bse[2]:.4f})"
          f"  p={m_a.pvalues[2]:.3f}{star(m_a.pvalues[2])}")
    print(f"    beta^LD/beta^delta = {ratio_a:.2f}")

    # [B] GFC level dummy only
    m_b = _nw_ols(s, pd.DataFrame(
        {"g_ld": gl, "g_d": gd, "d_gfc": sub["d_gfc"]}))
    b_ld_b = m_b.params[1]; b_d_b = m_b.params[2]
    ratio_b = b_ld_b / b_d_b if abs(b_d_b) > 1e-8 else float("nan")
    print(f"\n  [B] GFC level dummy only  R2={m_b.rsquared:.3f}")
    print(f"    {'beta^LD':<28} {b_ld_b:>9.4f}  ({m_b.bse[1]:.4f})"
          f"  p={m_b.pvalues[1]:.3f}{star(m_b.pvalues[1])}")
    print(f"    {'beta^delta':<28} {b_d_b:>9.4f}  ({m_b.bse[2]:.4f})"
          f"  p={m_b.pvalues[2]:.3f}{star(m_b.pvalues[2])}")
    print(f"    {'gamma(GFC level)':<28} {m_b.params[3]:>9.4f}  ({m_b.bse[3]:.4f})"
          f"  p={m_b.pvalues[3]:.3f}{star(m_b.pvalues[3])}")
    print(f"    beta^LD/beta^delta = {ratio_b:.2f}")

    # [C] GFC dummy + interactions (primary)
    m_c = _nw_ols(s, pd.DataFrame({
        "g_ld": gl, "g_d": gd, "d_gfc": sub["d_gfc"],
        "gl_x_gfc": gl_gfc, "gd_x_gfc": gd_gfc,
    }))
    b_ld_c = m_c.params[1]; b_d_c = m_c.params[2]
    ratio_c = b_ld_c / b_d_c if abs(b_d_c) > 1e-8 else float("nan")
    var_s  = np.var(s.values)
    sh_ld  = b_ld_c * np.cov(gl.values, s.values)[0, 1] / var_s
    sh_d   = b_d_c  * np.cov(gd.values, s.values)[0, 1] / var_s
    sh_res = 1.0 - sh_ld - sh_d

    print(f"\n  [C] GFC dummy + interactions (primary)  R2={m_c.rsquared:.3f}")
    print(f"    {chr(0x2500)*52}")
    print(f"    Normal-times coefficients:")
    print(f"    {'beta^LD':<28} {b_ld_c:>9.4f}  ({m_c.bse[1]:.4f})"
          f"  p={m_c.pvalues[1]:.3f}{star(m_c.pvalues[1])}")
    print(f"    {'beta^delta':<28} {b_d_c:>9.4f}  ({m_c.bse[2]:.4f})"
          f"  p={m_c.pvalues[2]:.3f}{star(m_c.pvalues[2])}")
    print(f"    beta^LD/beta^delta = {ratio_c:.2f}")
    print(f"    {chr(0x2500)*52}")
    print(f"    GFC controls:")
    print(f"    {'gamma(GFC level)':<28} {m_c.params[3]:>9.4f}  ({m_c.bse[3]:.4f})"
          f"  p={m_c.pvalues[3]:.3f}{star(m_c.pvalues[3])}")
    print(f"    {'delta^LD (g^LD x GFC)':<28} {m_c.params[4]:>9.4f}  ({m_c.bse[4]:.4f})"
          f"  p={m_c.pvalues[4]:.3f}{star(m_c.pvalues[4])}")
    print(f"    {'delta^d  (g^d x GFC)':<28} {m_c.params[5]:>9.4f}  ({m_c.bse[5]:.4f})"
          f"  p={m_c.pvalues[5]:.3f}{star(m_c.pvalues[5])}")
    print(f"    {chr(0x2500)*52}")
    print(f"    Variance decomp (approximate, normal-times cov):")
    print(f"    g^LD={sh_ld:.1%}  g^delta={sh_d:.1%}  residual={sh_res:.1%}")
    print(f"    {chr(0x2500)*52}")
    print(f"    Coefficient stability A to C:")
    print(f"    beta^LD:    {b_ld_a:.4f} -> {b_ld_b:.4f} -> {b_ld_c:.4f}"
          f"  (delta from baseline: {b_ld_c - b_ld_a:+.4f})")
    print(f"    beta^delta: {b_d_a:.4f} -> {b_d_b:.4f} -> {b_d_c:.4f}"
          f"  (delta from baseline: {b_d_c - b_d_a:+.4f})")
    stable_ld = "STABLE" if abs(b_ld_c - b_ld_a) < 0.3 else "SENSITIVE"
    stable_d  = "STABLE" if abs(b_d_c  - b_d_a)  < 0.5 else "SENSITIVE"
    print(f"    Stability: {stable_ld} for beta^LD, {stable_d} for beta^delta")
    print()


# 5. Correlation table (secondary)
# =============================================================================
def correlation_table(panel: pd.DataFrame,
                      sample_start: str, sample_end: str,
                      label: str) -> pd.DataFrame:
    """
    Pairwise correlations, demeaned levels only.
    First differences excluded — see run_regression docstring.
    """
    sub = panel[(panel["quarter_label"] >= sample_start) &
                (panel["quarter_label"] <= sample_end)].dropna(
                subset=["s_nat", "g_ld_nat", "g_delta_nat"]).copy()
    s  = sub["s_nat"]      - sub["s_nat"].mean()
    gl = sub["g_ld_nat"]   - sub["g_ld_nat"].mean()
    gd = sub["g_delta_nat"]- sub["g_delta_nat"].mean()
    return pd.DataFrame([{
        "sample":      label,
        "transform":   "Demeaned levels",
        "r(s,g_LD)":   round(s.corr(gl), 3),
        "r(s,g_δ)":    round(s.corr(gd), 3),
        "r(g_LD,g_δ)": round(gl.corr(gd), 3),
        "n":           len(s),
    }])


# =============================================================================
# Filter utilities (HP and Hamilton) + comparison regression
# =============================================================================
def apply_hp_filter(series: pd.Series, lam: float = 1600) -> pd.Series:
    """
    HP filter — extract cyclical component.

    Applied to the FULL available series (not just the estimation window)
    so that the estimation window 2001Q1–2019Q4 sits well inside the
    filter's reliable interior, minimising endpoint bias.

    Buffer assessment:
      s_nat    : 1948Q1–2026Q1 — enormous buffer on both sides
      g_ld_nat : 2000Q4–2026Q1 — 25Q post-sample buffer after 2019Q4  ✓
      g_delta  : 1992Q3–2021Q4 — 8Q post-sample buffer after 2019Q4
                 Borderline but acceptable (HP distortion zone ≈ 5–6Q)

    Data loss within estimation window: ZERO for all three series.

    λ=1600: standard choice for quarterly business cycle analysis.
    Shimer (2007) uses HP throughout for comparability.
    """
    from statsmodels.tsa.filters.hp_filter import hpfilter
    vals = series.dropna()
    cycle, _ = hpfilter(vals, lamb=lam)
    return pd.Series(cycle.values, index=vals.index, name=series.name)


def apply_hamilton_filter(series: pd.Series,
                          h: int = 4, p: int = 4) -> pd.Series:
    """
    Hamilton (2018) filter — project y_t onto y_{t-h}, …, y_{t-h-p+1}.
    Cyclical component = OLS residual.

    Why h=4 rather than the canonical h=8:
      The JOLTS series (g_ld_nat) starts 2000Q4. Hamilton h=8 requires
      t-h through t-h-p+1 = t-11, so the first filtered quarter would be
      2000Q4 + 11Q = 2003Q3, losing 2001Q1–2002Q2 (10 quarters = 13%
      of the estimation sample). Unacceptable.

      h=4 needs t-4 through t-7 (7 lags), so the first filtered quarter
      is 2000Q4 + 7Q = 2002Q3. This loses only 2001Q1–2002Q2 (6 quarters)
      from g_ld_nat. The other two series have long pre-samples and lose
      nothing within the estimation window.

      h=4 targets cycles of at least 4 quarters (1 year), which still
      removes the secular trend while retaining business cycle variation.

    Advantages over HP:
      - No endpoint bias by construction (residual at t uses only past data)
      - No spurious cycles (Hamilton 2018 critique of HP)
      - Stationary residuals by construction

    Reference: Hamilton, J.D. (2018). "Why You Should Never Use the
    Hodrick-Prescott Filter." Review of Economics and Statistics,
    100(5), 831–843.
    """
    vals = series.dropna().copy()
    df_h = pd.DataFrame({"y": vals})
    for i in range(p):
        df_h[f"y_lag{h+i}"] = df_h["y"].shift(h + i)
    df_h = df_h.dropna()
    X = sm.add_constant(df_h[[f"y_lag{h+i}" for i in range(p)]])
    resid = sm.OLS(df_h["y"], X).fit().resid
    return pd.Series(resid.values, index=df_h.index, name=series.name)


def filter_panel(panel: pd.DataFrame, method: str,
                 lam: float = 1600,
                 h: int = 4, p: int = 4) -> pd.DataFrame:
    """
    Apply trend filter to all three series using the FULL available data
    (entire panel, not just the estimation window) to maximise pre/post
    sample buffers and minimise endpoint or initialisation bias.

    Returns a copy of the panel with filtered (cyclical) series replacing
    the originals. The estimation window restriction (SAMPLE_START–
    SAMPLE_END) is applied in run_regression as usual.

    Data loss summary:
      HP     : zero data loss for all three series ✓
      Hamilton: zero loss for s_nat and g_delta; g_ld_nat loses
               2001Q1–2002Q2 (6 quarters) due to JOLTS start date.
    """
    if method not in ("hp", "hamilton"):
        raise ValueError(f"method must be 'hp' or 'hamilton', got {method!r}")

    filtered = panel.copy()
    tag = (f"HP λ={lam}" if method == "hp"
           else f"Hamilton h={h} p={p}")
    losses = {}

    for col in ["s_nat", "g_ld_nat", "g_delta_nat"]:
        full_series = (panel.set_index("quarter_label")[col]
                            .dropna()
                            .sort_index())
        if method == "hp":
            cyc = apply_hp_filter(full_series, lam=lam)
        else:
            cyc = apply_hamilton_filter(full_series, h=h, p=p)

        # Track data loss within estimation window
        est_index = [q for q in full_series.index
                     if SAMPLE_START <= q <= SAMPLE_END]
        filtered_est = [q for q in cyc.index
                        if SAMPLE_START <= q <= SAMPLE_END]
        lost = len(est_index) - len(filtered_est)
        losses[col] = lost

        cyc_df = cyc.reset_index()
        cyc_df.columns = ["quarter_label", col]
        filtered = (filtered.drop(columns=[col])
                             .merge(cyc_df, on="quarter_label", how="left"))

    print(f"\n  Filter: {tag}")
    for col, lost in losses.items():
        status = "✓" if lost == 0 else f"⚠  {lost} quarters lost"
        print(f"    {col:<16}: {status}")
    return filtered


# =============================================================================
# F2. Run filtered regressions — full comparison table
# =============================================================================
def run_filtered_comparison(panel: pd.DataFrame) -> None:
    """
    Run and print regressions under three specifications side by side:
      (1) Demeaned levels with GFC dummy + interactions  [existing]
      (2) HP filter λ=1600                               [new]
      (3) Hamilton filter h=4, p=4                       [new]

    For the filtered specs the GFC dummy/interactions are retained.
    Even in cyclically-filtered data the GFC episode is an outlier —
    the cyclical component of g_delta spiked sharply in 2008–2009 —
    and omitting the controls would push the interaction variation
    into β^δ as in spec A of the demeaned regression.

    The run_regression() function demeans continuous variables before
    estimation. For filtered series this demean is nearly a no-op
    (HP and Hamilton residuals have mean ≈ 0 by construction) but
    is harmless and ensures a consistent intercept interpretation.
    """
    print("\n" + "=" * 68)
    print("FILTERED REGRESSION COMPARISON — Full sample 2001Q1–2019Q4")
    print("=" * 68)
    print("Purpose: assess sensitivity of normal-times β^LD and β^δ to")
    print("trend-removal method.")
    print()
    print("Spec [C] (GFC dummy + interactions) reported for each filter.")
    print()

    specs = [
        ("Demeaned levels\n(GFC dummy + interactions)",
         panel, None),
        (f"HP filter λ=1600\n(GFC dummy + interactions)",
         filter_panel(panel, "hp", lam=1600), "hp"),
        (f"Hamilton h=4 p=4\n(GFC dummy + interactions)",
         filter_panel(panel, "hamilton", h=4, p=4), "hamilton"),
    ]

    # Collect spec-C results for comparison table
    rows = []
    for label, pnl, ftype in specs:
        sub = (pnl[(pnl["quarter_label"] >= SAMPLE_START) &
                   (pnl["quarter_label"] <= SAMPLE_END)]
               .dropna(subset=["s_nat", "g_ld_nat", "g_delta_nat"])
               .copy()
               .reset_index(drop=True))

        sub["d_gfc"] = ((sub["quarter_label"] >= GFC_START) &
                        (sub["quarter_label"] <= GFC_END)).astype(float)

        s  = sub["s_nat"]       - sub["s_nat"].mean()
        gl = sub["g_ld_nat"]    - sub["g_ld_nat"].mean()
        gd = sub["g_delta_nat"] - sub["g_delta_nat"].mean()

        X = sm.add_constant(pd.DataFrame({
            "g_ld":     gl,
            "g_d":      gd,
            "d_gfc":    sub["d_gfc"],
            "gl_x_gfc": gl * sub["d_gfc"],
            "gd_x_gfc": gd * sub["d_gfc"],
        }))
        m = sm.OLS(s.values, X.values).fit(
            cov_type="HAC",
            cov_kwds={"maxlags": NW_LAGS, "use_correction": True}
        )

        def star(p):
            return "***" if p<.01 else "**" if p<.05 else "*" if p<.10 else ""

        b_ld = m.params[1]; p_ld = m.pvalues[1]
        b_d  = m.params[2]; p_d  = m.pvalues[2]
        ratio = b_ld / b_d if abs(b_d) > 1e-8 else float("nan")
        rows.append({
            "label":   label.split("\n")[0],
            "n":       int(m.nobs),
            "β^LD":    b_ld,
            "p_LD":    p_ld,
            "s_LD":    star(p_ld),
            "β^δ":     b_d,
            "p_δ":     p_d,
            "s_δ":     star(p_d),
            "ratio":   ratio,
            "R²":      m.rsquared,
        })

    # Print comparison table
    print(f"  {'Specification':<22}  {'n':>4}  {'β^LD':>8}  {'p':>6}  "
          f"{'β^δ':>8}  {'p':>6}  {'ratio':>6}  {'R²':>5}")
    print(f"  {'─'*72}")
    for r in rows:
        print(f"  {r['label']:<22}  {r['n']:>4}  {r['β^LD']:>8.4f}  "
              f"{r['p_LD']:>6.3f}{r['s_LD']:<3}  {r['β^δ']:>8.4f}  "
              f"{r['p_δ']:>6.3f}{r['s_δ']:<3}  {r['ratio']:>6.2f}  "
              f"{r['R²']:>5.3f}")
    print()
    print("  Interpretation:")
    print("  Stable β^LD across rows → LD-channel loading is not a trend artifact")
    print("  Stable β^δ across rows → δ-channel loading robust to detrending")
    print("  Ratio >> 1 across all rows → s^nat primarily tracks LD in normal times")


# =============================================================================
# 6. Plots
# =============================================================================
def plot_comparison(panel: pd.DataFrame, out_path: Path) -> None:
    """
    Three-panel figure using raw quarterly rates (no z-scoring):
      Left   : time series of all three series in levels.
                s_nat and g_ld_nat share the left y-axis (both ~0–12%).
                g_delta_nat uses the right y-axis (smaller range, ~0–2%).
      Middle : scatter s_nat vs g_ld_nat (demeaned), with regression line.
      Right  : scatter s_nat vs g_delta_nat (demeaned), with regression line.

    Using raw rates rather than z-scores preserves the economically
    meaningful level information: s_nat ≈ 7%, g_ld ≈ 4%, g_delta ≈ 1%,
    and their differences reflect real structural distinctions between
    E→U hazards, layoff rates, and firm destruction rates.
    """
    plt.ioff()
    sub = panel[(panel["quarter_label"] >= SAMPLE_START) &
                (panel["quarter_label"] <= SAMPLE_END)].copy().reset_index(drop=True)

    pct = lambda x: x * 100   # decimal fraction → percent for y-axis labels

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # ── Left: time series overlay ──────────────────────────────────────────
    ax  = axes[0]
    ax2 = ax.twinx()   # right axis for g_delta (smaller scale)

    x  = range(len(sub))
    xt = [i for i, q in enumerate(sub["quarter_label"])
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]
    xl = [q[:4] for q in sub["quarter_label"]
          if q.endswith("Q1") and int(q[:4]) % 4 == 0]

    l1, = ax.plot(x, pct(sub["s_nat"]),    "#1f77b4", lw=2.0,
                  label=r"$s^{nat}_t$ (Shimer, left)")
    l2, = ax.plot(x, pct(sub["g_ld_nat"]), "#d62728", lw=1.8, ls="--",
                  label=r"$g^{LD}_t$ (JOLTS, left)")
    l3, = ax2.plot(x, pct(sub["g_delta_nat"]), "#2ca02c", lw=1.8, ls=":",
                   label=r"$g^{\delta}_t$ (BED, right)")

    # Shade GFC and COVID episode windows
    def _shade(ax, start_q, end_q, color, label):
        qs = list(sub["quarter_label"])
        i0 = next((i for i, q in enumerate(qs) if q >= start_q), None)
        i1 = next((i for i, q in enumerate(qs) if q > end_q),   len(qs)-1)
        if i0 is not None:
            ax.axvspan(i0, i1, alpha=0.12, color=color, label=label)

    _shade(ax, GFC_START, GFC_END, "orange", "GFC (2008Q3–2009Q4)")

    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8)
    ax.set_ylabel("Quarterly rate (% of employment)", fontsize=9)
    ax2.set_ylabel(r"$g^{\delta}$ quarterly rate (%)", fontsize=9,
                   color="#2ca02c")
    ax2.tick_params(axis="y", colors="#2ca02c")
    ax.set_title("Quarterly separation and closing rates\n"
                 r"(all in % of employment per quarter — shaded = episode windows)",
                 fontsize=10)
    ax.legend(handles=[l1, l2, l3] + ax.patches[:2],
              fontsize=7, loc="upper right")
    ax.grid(axis="y", lw=0.4, alpha=0.4)

    # ── Middle and right: scatter plots (demeaned) ─────────────────────────
    dm = lambda x: x - x.mean()
    s_dm  = dm(sub["s_nat"])
    ld_dm = dm(sub["g_ld_nat"])
    d_dm  = dm(sub["g_delta_nat"])

    for ax, xvals, color, xlabel in [
        (axes[1], ld_dm, "#d62728", r"$g^{LD}_t$ demeaned (% of empl.)"),
        (axes[2], d_dm,  "#2ca02c", r"$g^{\delta}_t$ demeaned (% of empl.)"),
    ]:
        r = s_dm.corr(xvals)
        xp = pct(xvals); yp = pct(s_dm)
        _ = ax.scatter(xp, yp, alpha=0.45, s=20, color=color,
                       edgecolors="none")
        xl_ = np.linspace(xp.min(), xp.max(), 50)
        m, b = np.polyfit(xp.dropna(), yp.dropna(), 1)
        _ = ax.plot(xl_, m*xl_+b, color=color, lw=1.6)
        _ = ax.set_title(fr"$s^{{nat}}$ vs {xlabel.split()[0]}"
                         f"\nr = {r:.3f}", fontsize=10)
        _ = ax.set_xlabel(xlabel, fontsize=9)
        _ = ax.set_ylabel(r"$s^{nat}$ demeaned (% of empl.)", fontsize=9)
        ax.axhline(0, color="black", lw=0.5)
        ax.axvline(0, color="black", lw=0.5)
        ax.grid(lw=0.4, alpha=0.3)

    fig.suptitle(
        r"National Shimer $s^{nat}_t$ vs. JOLTS LD and BED closing rates"
        "\n(raw quarterly rates in comparable units — no z-scoring)",
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
def main():
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
    print("\n[5] Regression analysis")
    print(f"    Sample: {SAMPLE_START}–{SAMPLE_END}")
    print(f"    GFC episode: {GFC_START}–{GFC_END}  (n=6 quarters)")
    print(f"    COVID excluded — CPS measurement breakdown (see constants block)")
    print()
    print("    Three nested specs per run:")
    print("    [A] Baseline — no episode controls")
    print("    [B] GFC level dummy only")
    print("    [C] GFC dummy + interactions (primary)")
    print()

    run_regression(panel, SAMPLE_START, SAMPLE_END,
                   f"Full sample {SAMPLE_START}–{SAMPLE_END}")
    run_regression(panel, SAMPLE_START, PREGFC_END,
                   f"Pre-GFC only {SAMPLE_START}–{PREGFC_END}")

    # ── 6. Correlation table (secondary) ─────────────────────────────────────
    print("\n[6] Pairwise correlations (secondary)")
    c1 = correlation_table(panel, SAMPLE_START, SAMPLE_END,
                           f"Full {SAMPLE_START}–{SAMPLE_END}")
    c2 = correlation_table(panel, SAMPLE_START, PREGFC_END,
                           f"Pre-GFC {SAMPLE_START}–{PREGFC_END}")
    print(pd.concat([c1, c2]).to_string(index=False))

    # ── 6b. Filtered regression comparison ──────────────────────────────────
    print("\n[6b] Filtered regression comparison (HP and Hamilton)")
    run_filtered_comparison(panel)

    # ── 7. Plots ──────────────────────────────────────────────────────────────
    print("\n[7] Plotting")
    plot_comparison(panel,
                    RESULTS_DIR / "shimer_national_comparison.png")
    plot_rolling_corr(panel,
                      RESULTS_DIR / "shimer_national_rolling_corr.png")

    print("\nDone.")


if __name__ == "__main__":
    main()