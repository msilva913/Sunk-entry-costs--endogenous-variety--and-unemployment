"""
part5_lp.py — Panel Local Projections for δ and s Bartik instruments
=====================================================================
Estimates impulse response functions (IRFs) of the state unemployment
rate to the job-destruction (δ) and job-finding (s) Bartik shocks via
Jordà (2005) local projections on the state-quarter panel.

Specification
-------------
For shock k ∈ {δ, s}, horizon h ∈ {0, ..., 16}, estimate separately:

    y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h · B_{s,t}^(k)
                             + γ₁ · u_{s,t-1} + γ₂ · log(LF_{s,t-1})
                             + ε_{s,t,h}

where:
    y_{s,t}       unemployment rate (%) in state s, quarter t
    α_s           state fixed effects — absorb permanent cross-state
                  unemployment level differences
    α_t           time fixed effects — absorb the common national cycle
                  so identification comes from cross-state dispersion
                  within each quarter, not aggregate fluctuations
    B_{s,t}^(k)  Bartik instrument: Σ_j ω_{s,j} · g_{-s,j,t}^(k)
                  (leave-one-out national shock rates, pre-determined
                  base-year industry shares as weights)
    u_{s,t-1}    lagged unemployment rate — controls for state's position
                  in Beveridge space entering the quarter
    log(LF_{s,t-1}) lagged log labor force — controls for labor market
                  scale; used instead of log population because (a) it is
                  the correct denominator for the unemployment rate, and
                  (b) it is already available from part4 LAUS fetch

Identification
--------------
The time FE removes the national aggregate shock each quarter. What
remains is cross-state dispersion: states with higher ω_{s,j} for
industries experiencing a large national shock receive a larger
B_{s,t} that quarter. β_h measures how much more the unemployment
rate rises (relative to t-1) in a high-exposure state than a
low-exposure state h quarters after the shock quarter.

Standard errors
---------------
Clustered at the state level throughout. With n=50 clusters this is
at the lower bound of reliability; supplement with wild cluster
bootstrap p-values (future work).

Sample periods
--------------
    δ instrument: 1992Q3 – 2023Q1  (~122 quarters)
    s instrument: 2001Q1 – 2023Q1  (~89 quarters)
At horizon h, the last usable t is (instrument end) – h quarters.
The δ IRF is also estimated on the post-2001Q1 subsample to ensure
comparability with the s IRF.

Outputs
-------
    data/results/lp_irf_delta.csv              — β_h, SE, CIs for δ shock (full)
    data/results/lp_irf_s.csv                 — β_h, SE, CIs for s-TS shock
    data/results/lp_irf_delta_post2001.csv    — δ on post-2001 sample
    data/results/lp_irf_delta_nfci.csv        — δ with NFCI interaction
    data/results/lp_irf_s_nfci.csv            — s with NFCI interaction
    data/results/lp_irf_delta_resid.csv       — δ residualized (post-2001)
    data/results/lp_irf_s_resid.csv           — s-TS residualized
    data/results/lp_irf_ld_resid.csv          — s-LD (layoffs & discharges) residualized
    data/results/lp_irf_qu_resid.csv          — s-QU (quits) residualized
    data/results/lp_irf_delta_resid_comparison.png   — δ raw vs. residualized
    data/results/lp_irf_s_decomp_comparison.png      — s raw vs TS/LD/QU residualized
    data/results/lp_irf_delta_ld_qu_comparison.png   — δ vs LD vs QU (core asymmetry)
    data/results/lp_irf_delta_vacancy.csv            — δ → vacancy IRF
    data/results/lp_irf_ld_vacancy.csv               — LD → vacancy IRF
    data/results/lp_irf_ts_vacancy.csv               — TS → vacancy IRF
    data/results/lp_irf_qu_vacancy.csv               — QU → vacancy IRF
    data/results/lp_irf_vacancy_delta_ld.png         — δ vs LD vacancy side-by-side
    data/results/lp_irf_vacancy_decomp_overlay.png   — δ vs LD vs QU vs TS vacancy overlay
    data/results/lp_irf_beveridge_asymmetry.png      — u vs v for δ and LD (4-panel)

Run
---
    python part5_lp.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates.py       (for δ)
    python part2_shock_rates_s.py     (for s)
    python part3_instrument.py        (produces delta_instrument_base2006.csv)
    python part3_instrument_s.py      (produces s_instrument_base2006.csv)
    python part4_outcomes.py          (produces laus_quarterly.parquet)
    python part7_nfci.py              (produces nfci_quarterly.parquet)
    python part2b_shock_comovement.py (produces residualized shock parquets)
    python part3_resid_instruments.py (produces LD/QU residualized instruments)
"""

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

def _open_file(path):
    """Open a saved file with the OS default viewer."""
    os.startfile(path)

try:
    os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
    # Running interactively (Spyder, Jupyter, etc.) — fall back to cwd.
    # Ensure your IDE working directory is set to "Bartek analysis".
    pass

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HORIZONS   = list(range(21))          # h = 0, 1, ..., 20 quarters
BASE_YEAR  = 2006                     # must match part1/part3 base year
CI_LEVEL   = 0.90                     # confidence band width for main plot
Z90        = 1.645
Z95        = 1.960

# Cap outcome quarter at 2019Q4 to exclude COVID (2020Q1–2021Q4) from the
# h-quarter-ahead outcome window. Without this, COVID enters the outcome at
# different horizons for δ (h=6+) and s (h=12+) depending on sample end
# dates, producing sharp spurious jumps in β_h at exactly those horizons.
# Set to None to use the full sample (useful as a robustness check).
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR  = Path("data/instruments")
RESULTS_DIR = Path("data/results")

DELTA_INSTR_FILE       = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
S_INSTR_FILE           = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
DELTA_RESID_INSTR_FILE = INSTR_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
S_RESID_INSTR_FILE     = INSTR_DIR / f"s_instrument_resid_base{BASE_YEAR}.csv"
LD_RESID_INSTR_FILE    = INSTR_DIR / f"ld_instrument_resid_base{BASE_YEAR}.csv"
QU_RESID_INSTR_FILE    = INSTR_DIR / f"qu_instrument_resid_base{BASE_YEAR}.csv"
LAUS_FILE        = INSTR_DIR / "laus_quarterly.parquet"
NFCI_FILE        = Path("data/cache") / "nfci_quarterly.parquet"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Quarter arithmetic
# ---------------------------------------------------------------------------
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    """
    Shift a Series of quarter_label strings (e.g. '2008Q4') by h quarters.
    h can be negative (lag) or positive (lead).
    """
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_instruments() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load δ and s Bartik instruments produced by part3 scripts."""
    delta = pd.read_csv(DELTA_INSTR_FILE, dtype={"state_fips": str})
    s     = pd.read_csv(S_INSTR_FILE,     dtype={"state_fips": str})

    # Ensure FIPS is zero-padded 2-digit string
    delta["state_fips"] = delta["state_fips"].str.zfill(2)
    s["state_fips"]     = s["state_fips"].str.zfill(2)

    print(f"  δ instrument: {delta.shape}  "
          f"({delta['quarter_label'].min()} – {delta['quarter_label'].max()})")
    print(f"  s instrument: {s.shape}  "
          f"({s['quarter_label'].min()} – {s['quarter_label'].max()})")
    return delta, s


def load_outcomes() -> pd.DataFrame:
    """
    Load quarterly LAUS unemployment rates, labor force, and vacancies from part4.

    Vacancies column is present when part4_outcomes.py has been run with BLS
    API access.  It is NaN for quarters before 2001 and absent entirely if the
    vacancy fetch failed; in the latter case downstream code falls back to
    unemployment-only LPs (guarded by `if "vacancies" in outcomes.columns`).
    """
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    print(f"  LAUS: {laus.shape}  "
          f"({laus['quarter_label'].min()} – {laus['quarter_label'].max()})")
    # Pass through vacancies if present; downstream functions check for column.
    keep_cols = ["state_fips", "state", "quarter_label", "unemp_rate", "labor_force"]
    if "vacancies" in laus.columns:
        keep_cols.append("vacancies")
        n_vac = laus["vacancies"].notna().sum()
        pct   = 100 * n_vac / len(laus)
        print(f"  Vacancies: {n_vac:,} non-NaN rows ({pct:.0f}% coverage)")
    else:
        print("  Vacancies: column not present — run part4_outcomes.py to add")
    return laus[keep_cols]


def load_nfci() -> pd.DataFrame:
    """
    Load quarterly NFCI data produced by part7_nfci.py.
    Returns DataFrame with columns: quarter_label, nfci, nfci_risk,
    nfci_credit, nfci_leverage, anfci.
    """
    nfci = pd.read_parquet(NFCI_FILE)
    print(f"  NFCI: {nfci.shape}  "
          f"({nfci['quarter_label'].min()} – {nfci['quarter_label'].max()})")
    print(f"  nfci_risk — mean={nfci['nfci_risk'].mean():.3f}  "
          f"SD={nfci['nfci_risk'].std():.3f}  "
          f"min={nfci['nfci_risk'].min():.3f}  "
          f"max={nfci['nfci_risk'].max():.3f}")
    return nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]]


# ---------------------------------------------------------------------------
# Panel construction
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame) -> pd.DataFrame:
    """
    Merge instrument with outcomes and construct pre-determined controls.

    Controls are lagged to t-1 (one quarter before the shock quarter t)
    to ensure they are predetermined with respect to the shock:
        u_{s,t-1}         — lagged unemployment rate
        log(LF_{s,t-1})   — lagged log labor force
        log(v_{s,t-1})    — lagged log vacancies (control only)
        vac_rate_{s,t-1}  — lagged vacancy rate in % (LP outcome)

    Returns a panel DataFrame with columns:
        state, instr_col, unemp_rate, unemp_lag1, lf_log_lag1
        [ + vacancies, log_v, log_v_lag1, vac_rate, vac_rate_lag1
          if outcomes contains vacancies ]
    """
    # Merge instrument and outcomes on state × quarter
    merge_cols = ["state_fips", "quarter_label", "unemp_rate", "labor_force"]
    _has_vacancies = "vacancies" in outcomes.columns
    if _has_vacancies:
        merge_cols.append("vacancies")

    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[merge_cols],
        on=["state_fips", "quarter_label"],
        how="inner",
    )

    # Rescale instrument from decimal fraction to percentage points.
    # Instruments are stored as e.g. 0.011 (= 1.1% job destruction rate).
    # After ×100, β_h is interpretable as: pp change in unemployment rate
    # per 1 pp change in the Bartik-predicted shock rate — the natural
    # structural unit for comparison with model predictions.
    panel[instr_col] = panel[instr_col] * 100

    # Sort for lag construction
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)

    # Construct lags within state — shift by 1 quarter
    # (groupby + shift respects the state boundary so lags don't bleed
    # across states, but we also verify the quarter continuity below)
    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )

    # Vacancy lags — only when vacancies are present
    # vac_rate (%) = vacancies * 1000 / labor_force * 100
    #   vacancies is in thousands; labor_force is in persons.
    #   Symmetric to unemp_rate: both are rates in percentage points.
    # log_v is kept as a lagged *control* variable (parallel to lf_log_lag1).
    if _has_vacancies:
        panel["log_v"] = np.log(panel["vacancies"].clip(lower=0.001))
        panel["log_v_lag1"] = panel.groupby("state_fips")["log_v"].shift(1)
        panel["vac_rate"] = (
            panel["vacancies"] * 1000 / panel["labor_force"] * 100
        )
        panel["vac_rate_lag1"] = (
            panel.groupby("state_fips")["vac_rate"].shift(1)
        )

    # Verify lag is truly t-1 (not a gap across non-contiguous quarters)
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    expected_prev    = quarter_shift(panel["quarter_label"], -1)
    bad_lag          = panel["prev_ql"].notna() & (panel["prev_ql"] != expected_prev)
    if bad_lag.any():
        n_bad = bad_lag.sum()
        print(f"  [warn] {n_bad} rows with non-contiguous quarter sequence — "
              f"lagged controls set to NaN for those rows")
        lag_cols = ["unemp_lag1", "lf_log_lag1"]
        if _has_vacancies:
            lag_cols += ["log_v_lag1", "vac_rate_lag1"]
        panel.loc[bad_lag, lag_cols] = np.nan

    panel = panel.drop(columns="prev_ql")

    n_states = panel["state_fips"].nunique()
    n_qtrs   = panel["quarter_label"].nunique()
    vac_note = ""
    if _has_vacancies:
        n_vac = panel["vac_rate_lag1"].notna().sum()
        vac_note = f"  vacancy lags: {n_vac:,} non-NaN"
    print(f"  Panel: {len(panel):,} rows  ({n_states} states × {n_qtrs} quarters)"
          + vac_note)
    return panel


# ---------------------------------------------------------------------------
# LP estimation — one horizon
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int,
                   shock_col: str,
                   include_nfci: bool = False,
                   outcome: str = "unemp") -> dict | None:
    """
    Estimate the LP for a single horizon h.

    Parameters
    ----------
    base_panel   : panel from build_panel() + attach_nfci_interaction()
    h            : forecast horizon (quarters)
    shock_col    : name of the instrument column (e.g. 'bartik_delta')
    include_nfci : if True, add B_{s,t} × NFCI_risk_dm interaction
    outcome      : "unemp" (default) or "vacancy"

        "unemp"   — dep_var = u_{s,t+h} - u_{s,t-1}
                    baseline control: unemp_lag1
                    interpretation: pp change in unemployment rate

        "vacancy" — dep_var = vac_rate_{s,t+h} - vac_rate_{s,t-1}
                    baseline control: vac_rate_lag1
                    interpretation: pp change in vacancy rate (symmetric to unemp)
                    Requires 'vac_rate' and 'vac_rate_lag1' columns in base_panel
                    (present when outcomes contained a 'vacancies' column).

    Returns a dict with point estimate, SEs, CIs, and diagnostics,
    or None if there are too few observations.
    """
    if outcome not in ("unemp", "vacancy"):
        raise ValueError(f"outcome must be 'unemp' or 'vacancy', got '{outcome}'")

    df = base_panel.copy()

    # --- Construct the h-quarter-ahead outcome ---
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql

    if outcome == "unemp":
        # Look up unemployment rate at t+h
        outcomes_lookup = (
            df[["state_fips", "quarter_label", "unemp_rate"]]
            .rename(columns={"quarter_label": "future_ql",
                             "unemp_rate":    "y_future"})
        )
        df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")

        if MAX_OUTCOME_QUARTER is not None:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

        # Outcome: u_{s,t+h} - u_{s,t-1}
        # Baseline t-1 is predetermined; see original note in prior version.
        df["dep_var"]    = df["y_future"] - df["unemp_lag1"]
        lag_control      = "unemp_lag1"

    else:  # outcome == "vacancy"
        # Require vac_rate columns — present only when vacancies were fetched.
        if "vac_rate" not in df.columns or "vac_rate_lag1" not in df.columns:
            print(f"    h={h:2d}: vacancy-rate columns not in panel — skipping")
            return None

        # Look up vac_rate (%) at t+h — symmetric to unemp_rate lookup above.
        vac_lookup = (
            df[["state_fips", "quarter_label", "vac_rate"]]
            .rename(columns={"quarter_label": "future_ql",
                             "vac_rate":      "vac_rate_future"})
        )
        df = df.merge(vac_lookup, on=["state_fips", "future_ql"], how="left")

        if MAX_OUTCOME_QUARTER is not None:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

        # Outcome: vac_rate_{s,t+h} - vac_rate_{s,t-1}   (pp change)
        # Exactly symmetric to the unemployment outcome.
        df["dep_var"] = df["vac_rate_future"] - df["vac_rate_lag1"]
        lag_control   = "vac_rate_lag1"

    # Drop rows with missing dep_var or any regressor
    required = ["dep_var", shock_col, lag_control, "lf_log_lag1"]
    if include_nfci:
        required.append("instr_x_nfci")
    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        print(f"    h={h:2d}: too few obs ({len(df)}) — skipping")
        return None

    # --- Build regressor matrix ---
    state_dummies = pd.get_dummies(df["state_fips"],    prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)

    core_regressors = [shock_col, lag_control, "lf_log_lag1"]
    if include_nfci:
        core_regressors.append("instr_x_nfci")

    X = pd.concat(
        [df[core_regressors], state_dummies, time_dummies],
        axis=1,
    )
    X = sm.add_constant(X, has_constant="add")
    y = df["dep_var"]

    # --- OLS with state-clustered SEs ---
    model = sm.OLS(y, X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values},
    )

    beta  = float(model.params[shock_col])
    se    = float(model.bse[shock_col])
    tstat = float(model.tvalues[shock_col])
    pval  = float(model.pvalues[shock_col])

    if include_nfci:
        delta_h    = float(model.params["instr_x_nfci"])
        delta_h_se = float(model.bse["instr_x_nfci"])
        delta_h_t  = float(model.tvalues["instr_x_nfci"])
        delta_h_p  = float(model.pvalues["instr_x_nfci"])
    else:
        delta_h = delta_h_se = delta_h_t = delta_h_p = None

    partial_f = tstat ** 2

    return {
        "h":           h,
        "beta":        beta,
        "se":          se,
        "tstat":       tstat,
        "pval":        pval,
        "partial_f":   partial_f,
        "ci90_lo":     beta - Z90 * se,
        "ci90_hi":     beta + Z90 * se,
        "ci95_lo":     beta - Z95 * se,
        "ci95_hi":     beta + Z95 * se,
        "delta_h":     delta_h,
        "delta_h_se":  delta_h_se,
        "delta_h_t":   delta_h_t,
        "delta_h_p":   delta_h_p,
        "nobs":        int(model.nobs),
        "r2":          float(model.rsquared),
        "n_clusters":  df["state_fips"].nunique(),
        "qt_range":    f"{df['quarter_label'].min()}–{df['quarter_label'].max()}",
        "outcome":     outcome,
    }


# ---------------------------------------------------------------------------
# LP estimation — all horizons
# ---------------------------------------------------------------------------
def run_lp(base_panel: pd.DataFrame, shock_col: str,
           label: str,
           include_nfci: bool = False,
           outcome: str = "unemp") -> pd.DataFrame:
    """
    Run LP for all horizons 0..20 and return results as a DataFrame.

    Args:
        base_panel:   output of build_panel() with attach_nfci_interaction()
        shock_col:    'bartik_delta', 'bartik_s', 'bartik_ld', etc.
        label:        human-readable label for printing (e.g. 'δ shock')
        include_nfci: if True, include B_{s,t} × NFCI_risk_dm interaction.
        outcome:      "unemp" (default) or "vacancy"
                      "vacancy" requires log_v/log_v_lag1 columns in panel.
    """
    nfci_tag    = " [+NFCI interaction]" if include_nfci else ""
    outcome_tag = " → log(vacancies)" if outcome == "vacancy" else " → unemp. rate"
    print(f"\n{'='*60}")
    print(f"  LP — {label}{nfci_tag}{outcome_tag}  (h = 0 … {max(HORIZONS)})")
    print(f"{'='*60}")
    if include_nfci:
        print(f"  {'h':>3}  {'β_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}  "
              f"{'δ_h':>10}  {'δ_h SE':>8}  {'δ_h p':>6}  {'N':>6}")
    else:
        print(f"  {'h':>3}  {'β_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}  "
              f"{'partial-F':>10}  {'N':>6}")
    print(f"  {'-'*60}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(base_panel, h, shock_col,
                             include_nfci=include_nfci,
                             outcome=outcome)
        if res is None:
            continue
        rows.append(res)
        sig = ("***" if res["pval"] < 0.01 else
               "**"  if res["pval"] < 0.05 else
               "*"   if res["pval"] < 0.10 else "")
        if include_nfci:
            dsig = ("***" if res["delta_h_p"] < 0.01 else
                    "**"  if res["delta_h_p"] < 0.05 else
                    "*"   if res["delta_h_p"] < 0.10 else "")
            print(f"  {h:3d}  {res['beta']:10.4f}  {res['se']:8.4f}  "
                  f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
                  f"{res['delta_h']:10.4f}  {res['delta_h_se']:8.4f}  "
                  f"{res['delta_h_p']:6.3f}{dsig}  {res['nobs']:6d}  {sig}")
        else:
            print(f"  {h:3d}  {res['beta']:10.4f}  {res['se']:8.4f}  "
                  f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
                  f"{res['partial_f']:10.2f}  {res['nobs']:6d}  {sig}")

    # ── Standardize to unit-SD scale ───────────────────────────────────────
    # Multiply β_h, SE, and CIs by the instrument's cross-sectional SD so
    # coefficients are in per-1-SD-shock units.  Works identically for both
    # unemployment and vacancy outcomes.
    instr_sd = float(base_panel[shock_col].std())
    scale_cols = ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]
    irf = pd.DataFrame(rows)
    if irf.empty:
        return irf
    irf[scale_cols] = irf[scale_cols] * instr_sd
    if include_nfci:
        irf["delta_h"]    = irf["delta_h"]    * instr_sd
        irf["delta_h_se"] = irf["delta_h_se"] * instr_sd
    irf["instr_sd"] = instr_sd
    irf["outcome"]  = outcome    # tag for plot labelling
    return irf


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_irf(results: dict[str, pd.DataFrame],
             out_path: Path,
             instr_sds: dict[str, float] | None = None) -> None:
    """
    Side-by-side IRF plots.  Each β_h is already scaled to 1-SD units by
    run_lp(), so all panels share the same y-axis interpretation:
    pp change in unemployment rate per 1-SD shock.

    results:   dict mapping label → DataFrame from run_lp()
    instr_sds: dict mapping label → instrument SD (stored in irf["instr_sd"]);
               kept as parameter for backward compatibility but no longer
               used for scaling (scaling already done in run_lp).
    """
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6), sharey=False)
    if n == 1:
        axes = [axes]

    color_cycle = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
                   "#9467bd", "#8c564b"]

    peak_info = {}

    for (ax, (label, df)), color in zip(
            zip(axes, results.items()), color_cycle):
        if df.empty:
            ax.set_visible(False)
            continue

        h    = df["h"].values
        beta = df["beta"].values
        lo90 = df["ci90_lo"].values
        hi90 = df["ci90_hi"].values
        lo95 = df["ci95_lo"].values
        hi95 = df["ci95_hi"].values

        # Retrieve per-IRF instrument SD stored by run_lp (for caption).
        isd = float(df["instr_sd"].iloc[0]) if "instr_sd" in df.columns else np.nan
        max_h = int(h.max())

        # 95% CI contour
        ax.plot(h, lo95, color=color, linewidth=0.7, linestyle="--", alpha=0.6)
        ax.plot(h, hi95, color=color, linewidth=0.7, linestyle="--", alpha=0.6)
        # 90% CI band
        ax.fill_between(h, lo90, hi90, color=color, alpha=0.20,
                        label="90% CI")
        # Point estimate
        ax.plot(h, beta, color=color, linewidth=2.0,
                marker="o", markersize=3.5, label=label)

        ax.axhline(0, color="black", linewidth=0.8)
        for hh in range(4, max_h + 1, 4):
            ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)

        # Detect outcome type from results DataFrame if tagged
        _outcome_type = (df["outcome"].iloc[0]
                         if "outcome" in df.columns else "unemp")
        _ylabel = ("pp change in vacancy rate\nper 1-SD shock"
                   if _outcome_type == "vacancy"
                   else "pp change in unemp. rate\nper 1-SD shock")

        ax.set_title(f"IRF — {label}", fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel(_ylabel, fontsize=10)
        ax.set_xticks(h[::2] if max_h > 16 else h)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: str(int(x))))
        ax.legend(fontsize=9, framealpha=0.85)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)

        idx_peak = int(np.argmax(np.abs(beta)))
        peak_info[label] = {
            "h": int(h[idx_peak]), "beta": float(beta[idx_peak]),
            "isd": isd, "max_h": max_h,
        }

    # Classify instruments: raw (rate units) vs. residualized (log-residual units).
    labels_list = list(results.keys())
    n_resid = sum("resid" in lbl.lower() for lbl in labels_list)
    n_raw   = len(labels_list) - n_resid
    if n_resid == 0:
        instr_type  = "raw Bartik"
        units_note  = ("All instruments are raw Bartik rates; "
                       "1-SD units are directly comparable across panels.")
    elif n_raw == 0:
        instr_type  = "productivity-residualized Bartik"
        units_note  = ("All instruments are residualized on $\\Delta\\log p_t$; "
                       "1-SD units are directly comparable across panels.")
    else:
        instr_type  = "Bartik (mixed raw and residualized)"
        units_note  = ("Raw instruments are in quarterly rate units; "
                       "residualized instruments are in log-residual units. "
                       "1-SD units are NOT directly comparable across raw and "
                       "residualized panels — compare shapes, not magnitudes.")

    labels_str = ", ".join(labels_list)
    max_h_all  = max(v["max_h"] for v in peak_info.values()) if peak_info else 20

    # Detect whether panels are vacancy or unemployment outcomes (or a mix).
    _outcome_types = set()
    for df in results.values():
        if not df.empty and "outcome" in df.columns:
            _outcome_types.add(df["outcome"].iloc[0])
    _has_vac   = "vacancy" in _outcome_types
    _has_unemp = "unemp" in _outcome_types
    if _has_vac and _has_unemp:
        _outcome_label = "unemployment rate & vacancy rate"
        _controls_label = "$u_{{t-1}}$ / $v_{{t-1}}$, $\\log LF_{{t-1}}$"
        _beta_interp = ("pp change in unemployment rate (or pp change in "
                        "vacancy rate) per 1-SD increase in the Bartik instrument")
        _peak_unit = ""   # mixed — defer to per-panel ylabel
    elif _has_vac:
        _outcome_label = "vacancy rate"
        _controls_label = "$v_{{t-1}}$, $\\log LF_{{t-1}}$"
        _beta_interp = ("pp change in vacancy rate per 1-SD increase "
                        "in the Bartik instrument")
        _peak_unit = " pp"
    else:
        _outcome_label = "unemployment rate"
        _controls_label = "$u_{{t-1}}$, $\\log LF_{{t-1}}$"
        _beta_interp = ("pp change in the unemployment rate per 1-SD increase "
                        "in the Bartik instrument")
        _peak_unit = " pp"

    fig.suptitle(
        f"Panel LP — {instr_type} instruments → {_outcome_label}\n"
        f"Shocks: {labels_str}  |  "
        f"State + time FEs; controls: {_controls_label}; "
        "SEs clustered by state; outcomes capped 2019Q4",
        fontsize=10, y=1.01,
    )

    # Peak-effect sentences — omit raw instrument SD (non-interpretable units).
    peak_sentences = []
    for lbl, info in peak_info.items():
        peak_sentences.append(
            f"{lbl}: peak {info['beta']:+.2f}{_peak_unit} at h={info['h']}."
        )

    caption = (
        f"Notes: Each panel plots $\\hat{{\\beta}}_h$, $h=0,\\ldots,{max_h_all}$, "
        "from separate OLS regressions of $y_{{s,t+h}} - y_{{s,t-1}}$ on the "
        f"{instr_type} instrument, state FEs, time FEs, {_controls_label}.  "
        "Coefficients are scaled to 1-SD units (multiplied by the pooled "
        f"cross-sectional SD of the instrument), so $\\hat{{\\beta}}_h$ measures the "
        f"{_beta_interp}.  "
        + units_note + "  "
        "Baseline: $t-1$ (predetermined w.r.t. shock); "
        "$\\hat{{\\beta}}_0$ measures the within-quarter impact.  "
        "Outcome quarters capped at 2019Q4 to exclude COVID distortions.  "
        + "  ".join(peak_sentences) + "  "
        "Shaded bands: 90\\% CI; dashed lines: 95\\% CI.  "
        "SEs clustered by state (50 clusters)."
    )

    fig.text(
        0.5, -0.04, caption,
        ha="center", va="top", fontsize=7.5,
        wrap=True, transform=fig.transFigure,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                  edgecolor="#cccccc", linewidth=0.8),
        multialignment="left",
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


# ---------------------------------------------------------------------------
# Overlay IRF plot — multiple series on a single panel
# ---------------------------------------------------------------------------
def overlay_plot_irf(
    results  : dict,
    out_path : Path,
    title    : str = "",
    ylabel   : str = "pp change in unemp. rate per 1-SD shock",
) -> None:
    """
    Overlay multiple IRFs on one panel for direct shape comparison.
    results: dict mapping label -> DataFrame from run_lp() (already 1-SD scaled).
    """
    PALETTE = [
        ("#1f77b4", "-"),   # δ
        ("#d62728", "-"),   # LD
        ("#2ca02c", "-"),   # QU
        ("#ff7f0e", "--"),  # TS
        ("#9467bd", ":"),   # extra
        ("#8c564b", "-."),  # extra
    ]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(0, color="black", linewidth=0.8)
    for hh in [4, 8, 12, 16, 20]:
        ax.axvline(hh, color="grey", linewidth=0.4, linestyle=":", alpha=0.5)

    for idx, (label, df) in enumerate(results.items()):
        if df is None or df.empty:
            continue
        color, ls = PALETTE[idx % len(PALETTE)]
        h    = df["h"].values
        beta = df["beta"].values
        lo90 = df["ci90_lo"].values
        hi90 = df["ci90_hi"].values
        ax.fill_between(h, lo90, hi90, color=color, alpha=0.12)
        ax.plot(h, beta, color=color, linewidth=2.0, linestyle=ls,
                marker="o", markersize=3.0, label=label)

    ax.set_xlabel("Horizon h (quarters)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title if title else "IRF comparison", fontsize=12)
    ax.legend(fontsize=10, framealpha=0.85)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")
    _open_file(out_path)


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5 — Panel Local Projections")
print("=" * 60)

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------
print("\n[1] Loading instruments, outcomes, and NFCI")
delta_instr, s_instr = load_instruments()
outcomes = load_outcomes()
nfci     = load_nfci()

# -----------------------------------------------------------------------
# Build panels
# -----------------------------------------------------------------------
print("\n[2] Building panels")
print("  δ panel:")
delta_panel = build_panel(delta_instr, "bartik_delta", outcomes)

print("  s panel:")
s_panel     = build_panel(s_instr,     "bartik_s",     outcomes)

# δ panel restricted to post-2001Q1 for comparability with s instrument
print("  δ panel (post-2001):")
delta_post2001 = delta_panel[
    delta_panel["quarter_label"] >= "2001Q1"
].copy().reset_index(drop=True)

# -----------------------------------------------------------------------
# Merge NFCI and construct demeaned interaction term
# -----------------------------------------------------------------------
print("\n[2b] Merging NFCI and constructing interaction term")

def attach_nfci_interaction(panel: pd.DataFrame,
                             nfci: pd.DataFrame,
                             instr_col: str) -> pd.DataFrame:
    """
    Merge quarterly NFCI into panel and construct the demeaned interaction
    term for the financial conditions robustness specification.

    Demeaning is done within the panel's own estimation sample — i.e. the
    mean of nfci_risk is computed over the quarters actually present in the
    panel after merging, not over the full NFCI history. This ensures that
    the interaction term equals zero at average financial conditions
    prevailing during the estimation window, making β_h interpretable as
    the IRF at typical sample-period financial conditions.

    Two interaction columns are added:
        nfci_risk_dm     : demeaned risk subindex  (nfci_risk − mean)
        instr_x_nfci     : B_{s,t} × nfci_risk_dm  (the interaction regressor)

    The interaction is NOT included in the baseline LP regressions —
    it is available in the panel for use in robustness specifications.
    """
    panel = panel.merge(
        nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]],
        on="quarter_label",
        how="left",
    )

    # Warn if any quarters in the panel have no NFCI match
    n_missing = panel["nfci_risk"].isna().sum()
    if n_missing > 0:
        missing_qtrs = (panel.loc[panel["nfci_risk"].isna(), "quarter_label"]
                        .unique().tolist())
        print(f"  [warn] {n_missing} panel rows have no NFCI match "
              f"({len(missing_qtrs)} quarters): {missing_qtrs[:5]} ...")

    # Demean within the panel's own sample
    nfci_mean = panel["nfci_risk"].mean()
    nfci_sd   = panel["nfci_risk"].std()
    panel["nfci_risk_dm"] = panel["nfci_risk"] - nfci_mean

    print(f"  {instr_col}: nfci_risk sample mean = {nfci_mean:.4f}  "
          f"SD = {nfci_sd:.4f}  "
          f"(demeaned range: [{panel['nfci_risk_dm'].min():.3f}, "
          f"{panel['nfci_risk_dm'].max():.3f}])")

    # Interaction: instrument (already in pp) × demeaned risk index
    panel["instr_x_nfci"] = panel[instr_col] * panel["nfci_risk_dm"]

    return panel

delta_panel    = attach_nfci_interaction(delta_panel,    nfci, "bartik_delta")
s_panel        = attach_nfci_interaction(s_panel,        nfci, "bartik_s")
delta_post2001 = attach_nfci_interaction(delta_post2001, nfci, "bartik_delta")
n_post = delta_post2001["quarter_label"].nunique()
print(f"    Restricted to 2001Q1+: {len(delta_post2001):,} rows "
      f"({delta_post2001['state_fips'].nunique()} states × {n_post} quarters)")

# -----------------------------------------------------------------------
# Instrument summary statistics (after FE removal — within variation)
# -----------------------------------------------------------------------
print("\n[3] Instrument variation summary (cross-sectional SD by quarter)")
for lbl, pnl, col in [("δ", delta_panel, "bartik_delta"),
                       ("s", s_panel,     "bartik_s")]:
    within_sd = (
        pnl.groupby("quarter_label")[col].std()
        .describe()[["mean", "min", "max"]]
    )
    print(f"  {lbl}: cross-state SD  mean={within_sd['mean']:.4f}  "
          f"min={within_sd['min']:.4f}  max={within_sd['max']:.4f}")

# -----------------------------------------------------------------------
# Run LPs
# -----------------------------------------------------------------------
print("\n[4] Running panel LPs")

irf_delta       = run_lp(delta_panel,    "bartik_delta", "δ shock")
irf_s           = run_lp(s_panel,        "bartik_s",     "s shock")
irf_delta_post  = run_lp(delta_post2001, "bartik_delta", "δ shock (post-2001)")

# Interaction robustness: baseline + B × NFCI_risk_dm
irf_delta_nfci  = run_lp(delta_panel,    "bartik_delta", "δ shock",
                          include_nfci=True)
irf_s_nfci      = run_lp(s_panel,        "bartik_s",     "s shock",
                          include_nfci=True)

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
print("\n[5] Saving results")

irf_delta.to_csv(RESULTS_DIR / "lp_irf_delta.csv", index=False)
irf_s.to_csv(RESULTS_DIR / "lp_irf_s.csv", index=False)
irf_delta_post.to_csv(RESULTS_DIR / "lp_irf_delta_post2001.csv", index=False)
irf_delta_nfci.to_csv(RESULTS_DIR / "lp_irf_delta_nfci.csv", index=False)
irf_s_nfci.to_csv(RESULTS_DIR / "lp_irf_s_nfci.csv", index=False)

print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_post2001.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_nfci.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_nfci.csv'}")

# -----------------------------------------------------------------------
# Plot IRFs
# -----------------------------------------------------------------------
print("\n[6] Plotting IRFs")

# Instrument SDs after ×100 rescaling — used in caption for 1-SD effect
instr_sds_main = {
    "δ shock": float(delta_panel["bartik_delta"].std()),
    "s shock": float(s_panel["bartik_s"].std()),
}

# Main comparison: δ (full sample) and s side by side
plot_irf(
    {"δ shock": irf_delta, "s shock": irf_s},
    out_path=RESULTS_DIR / "lp_irf_combined.png",
    instr_sds=instr_sds_main,
)

# Sample robustness: δ full vs δ post-2001
if not irf_delta_post.empty:
    instr_sds_sample = {
        "δ shock":           float(delta_panel["bartik_delta"].std()),
        "δ shock (post-2001)": float(delta_post2001["bartik_delta"].std()),
    }
    plot_irf(
        {"δ shock": irf_delta, "δ shock (post-2001)": irf_delta_post},
        out_path=RESULTS_DIR / "lp_irf_delta_sample_check.png",
        instr_sds=instr_sds_sample,
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_sample_check.png'}")

_open_file(RESULTS_DIR / "lp_irf_combined.png")

# -----------------------------------------------------------------------
# Residualized instruments — LP (δ, TS, LD, QU — if available)
# -----------------------------------------------------------------------
_resid_files_present = all(p.exists() for p in [
    DELTA_RESID_INSTR_FILE, S_RESID_INSTR_FILE,
    LD_RESID_INSTR_FILE, QU_RESID_INSTR_FILE,
])

if _resid_files_present:
    print("\n[7] Running LPs with productivity-residualized instruments "
          "(δ, TS, LD, QU)")

    def _load_resid(path):
        df = pd.read_csv(path, dtype={"state_fips": str})
        df["state_fips"] = df["state_fips"].str.zfill(2)
        return df

    delta_resid_raw = _load_resid(DELTA_RESID_INSTR_FILE)
    s_resid_raw     = _load_resid(S_RESID_INSTR_FILE)
    ld_resid_raw    = _load_resid(LD_RESID_INSTR_FILE)
    qu_resid_raw    = _load_resid(QU_RESID_INSTR_FILE)

    # LD and QU CSVs use bartik_ld / bartik_qu column names.
    # build_panel is generic — pass the correct instr_col for each.
    outcomes = load_outcomes()
    delta_resid_panel = build_panel(delta_resid_raw, "bartik_delta", outcomes)
    s_resid_panel     = build_panel(s_resid_raw,     "bartik_s",     outcomes)
    ld_resid_panel    = build_panel(ld_resid_raw,    "bartik_ld",    outcomes)
    qu_resid_panel    = build_panel(qu_resid_raw,    "bartik_qu",    outcomes)

    # Restrict δ to post-2001 for direct comparability with s/LD/QU.
    delta_resid_post = delta_resid_panel[
        delta_resid_panel["quarter_label"] >= "2001Q1"
    ].copy()

    # ── run LPs ──────────────────────────────────────────────────────────
    irf_delta_resid = run_lp(delta_resid_post, "bartik_delta",
                              "δ shock (resid, post-2001)")
    #irf_s_resid     = run_lp(s_resid_panel,    "bartik_s",
     #                         "s-TS shock (resid)")
    irf_ld_resid    = run_lp(ld_resid_panel,   "bartik_ld",
                              "s-LD shock (resid)")
    irf_qu_resid    = run_lp(qu_resid_panel,   "bartik_qu",
                              "s-QU shock (resid)")

    # ── save ─────────────────────────────────────────────────────────────
    for irf, fname in [
        (irf_delta_resid, "lp_irf_delta_resid.csv"),
        (irf_s_resid,     "lp_irf_s_resid.csv"),
        (irf_ld_resid,    "lp_irf_ld_resid.csv"),
        (irf_qu_resid,    "lp_irf_qu_resid.csv"),
    ]:
        irf.to_csv(RESULTS_DIR / fname, index=False)
        print(f"  Saved: {RESULTS_DIR / fname}")

    # ── plots ─────────────────────────────────────────────────────────────
    # Plot 1: δ raw (post-2001) vs δ residualized — side-by-side
    plot_irf(
        {"δ raw (post-2001)": irf_delta_post,
         "δ residualized":    irf_delta_resid},
        out_path=RESULTS_DIR / "lp_irf_delta_resid_comparison.png",
        instr_sds={
            "δ raw (post-2001)": float(delta_post2001["bartik_delta"].std()),
            "δ residualized":    float(delta_resid_post["bartik_delta"].std()),
        },
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_resid_comparison.png'}")

    # Plot 2 (CORE): δ vs LD vs QU residualized overlaid — asymmetry test
    # Model prediction: δ destroys vacancies (sharper initial response),
    # LD triggers reposting (more muted), QU is procyclical (different shape).
    overlay_plot_irf(
        {"δ residualized":  irf_delta_resid,
         "LD residualized": irf_ld_resid,
         "QU residualized": irf_qu_resid},
        out_path=RESULTS_DIR / "lp_irf_delta_ld_qu_overlay.png",
        title=(
            r"IRF: $\delta$ vs Layoffs+Discharges vs Quits"
            "\n(productivity-residualized, 1-SD scale, post-2001)"
        ),
    )

    # Plot 3: LD vs QU vs TS residualized — s decomposition
    overlay_plot_irf(
        {"LD residualized": irf_ld_resid,
         "QU residualized": irf_qu_resid},
         # "TS residualized": irf_s_resid},
        out_path=RESULTS_DIR / "lp_irf_s_decomp_overlay.png",
        title=(
            "IRF: Layoffs+Discharges vs Quits vs Total Separations\n"
            "(productivity-residualized, 1-SD scale)"
        ),
    )

else:
    print("\n[7] Residualized instruments not found — skipping resid LP.")
    print("    Run part2_shock_rates_s.py then part3_resid_instruments.py first.")

# -----------------------------------------------------------------------
# Summary table — raw and residualized peaks
# -----------------------------------------------------------------------
print("\n[8] Summary — peak unemployment response (1-SD standardized)")
print(f"  {'Series':<32}  {'h_peak':>6}  {'β_peak':>8}  "
      f"{'SE':>8}  {'p':>6}  {'partial-F':>10}")
print(f"  {'-'*74}")

raw_series = [
    ("δ raw (full)",          irf_delta),
    ("δ raw (post-2001)",     irf_delta_post),
    ("s raw (TS)",            irf_s),
]
for lbl, df in raw_series:
    if df is None or df.empty:
        continue
    pk = df.loc[df["beta"].idxmax()]
    print(f"  {lbl:<32}  {int(pk['h']):>6}  {pk['beta']:>8.4f}  "
          f"{pk['se']:>8.4f}  {pk['pval']:>6.3f}  {pk['partial_f']:>10.2f}")

try:
    print(f"  {'':32}  {'':>6}")
    resid_series = [
        ("δ residualized (post-2001)",  irf_delta_resid),
        ("LD residualized",             irf_ld_resid),
        ("QU residualized",             irf_qu_resid),
       # ("TS residualized",             irf_s_resid),
    ]
    for lbl, df in resid_series:
        if df is None or df.empty:
            continue
        pk = df.loc[df["beta"].idxmax()]
        print(f"  {lbl:<32}  {int(pk['h']):>6}  {pk['beta']:>8.4f}  "
              f"{pk['se']:>8.4f}  {pk['pval']:>6.3f}  {pk['partial_f']:>10.2f}")
except NameError:
    print("  (residualized LPs not run — run part3_resid_instruments.py first)")

# -----------------------------------------------------------------------
# [9] Vacancy outcome LPs — the core asymmetry test
#
# Model prediction (from vacancy law of motion):
#   δ shock → δ_{e,t} rises → surviving unmatched vacancies + reposting
#             BOTH collapse → vacancies FALL, unemployment RISES
#             → Beveridge curve dynamics: UV move together, not counter
#
#   LD shock → s_{t-1} rises but δ_{e,t} unchanged → surviving firms
#              repost → vacancies FLAT or RISE, unemployment RISES
#              → positive UV co-movement (off the Beveridge curve)
#
# This is the single strongest direct test of the model's mechanism.
# Vacancies require part4_outcomes.py to have been run with BLS API access.
# -----------------------------------------------------------------------
print("\n[9] Vacancy outcome LPs — core asymmetry test")

_vac_available = "vacancies" in outcomes.columns and outcomes["vacancies"].notna().any()
if not _vac_available:
    print("  Vacancies not available — run part4_outcomes.py with BLS API access")
    print("  Skipping vacancy LPs.")
else:
    n_vac_obs = outcomes["vacancies"].notna().sum()
    vac_start = outcomes.loc[outcomes["vacancies"].notna(), "quarter_label"].min()
    vac_end   = outcomes.loc[outcomes["vacancies"].notna(), "quarter_label"].max()
    print(f"  Vacancies available: {n_vac_obs:,} obs  ({vac_start} – {vac_end})")

    # --- Panels (re-merge to pick up vacancies; build_panel now includes them) ---
    print("\n  Rebuilding panels with vacancy columns:")
    print("    δ panel (post-2001, for comparability with s/LD):")
    _delta_vac_raw   = _load_resid(DELTA_RESID_INSTR_FILE) if _resid_files_present \
                       else delta_instr
    _delta_col       = "bartik_delta"
    _delta_vac_panel = build_panel(_delta_vac_raw, _delta_col, outcomes)
    _delta_vac_post  = _delta_vac_panel[
        _delta_vac_panel["quarter_label"] >= "2001Q1"
    ].copy().reset_index(drop=True)
    _delta_vac_post  = attach_nfci_interaction(_delta_vac_post, nfci, _delta_col)

    print("    LD panel (primary s-type for vacancy test — reposting channel):")
    _ld_col = "bartik_ld"
    _ld_src = ld_resid_raw if _resid_files_present else \
              pd.read_csv(INSTR_DIR / f"ld_instrument_base{BASE_YEAR}.csv",
                          dtype={"state_fips": str}).assign(
                              state_fips=lambda d: d["state_fips"].str.zfill(2))
    _ld_vac_panel = build_panel(_ld_src, _ld_col, outcomes)
    _ld_vac_panel = attach_nfci_interaction(_ld_vac_panel, nfci, _ld_col)

    # print("    TS panel (total separations — vacancy outcome):")
    # _ts_col = "bartik_s"
    # _ts_src = s_resid_raw if _resid_files_present else s_instr
    # _ts_vac_panel = build_panel(_ts_src, _ts_col, outcomes)
    # _ts_vac_panel = attach_nfci_interaction(_ts_vac_panel, nfci, _ts_col)

    print("    QU panel (quits — vacancy outcome):")
    _qu_col = "bartik_qu"
    _qu_src = qu_resid_raw if _resid_files_present else \
              pd.read_csv(INSTR_DIR / f"qu_instrument_base{BASE_YEAR}.csv",
                          dtype={"state_fips": str}).assign(
                              state_fips=lambda d: d["state_fips"].str.zfill(2))
    _qu_vac_panel = build_panel(_qu_src, _qu_col, outcomes)
    _qu_vac_panel = attach_nfci_interaction(_qu_vac_panel, nfci, _qu_col)

    # Use the residualized instruments if available (preferred), raw otherwise
    _instr_note = ("residualized" if _resid_files_present else "raw (run part3_resid_instruments.py for residualized)")
    print(f"    Using {_instr_note} instruments for vacancy LPs")

    # --- Run vacancy LPs ---
    print()
    irf_delta_vac = run_lp(_delta_vac_post, _delta_col,
                            "δ shock — vacancies", outcome="vacancy")
    irf_ld_vac    = run_lp(_ld_vac_panel,   _ld_col,
                            "LD shock — vacancies", outcome="vacancy")
    #irf_ts_vac    = run_lp(_ts_vac_panel,   _ts_col,
               #             "TS shock — vacancies", outcome="vacancy")
    irf_qu_vac    = run_lp(_qu_vac_panel,   _qu_col,
                            "QU shock — vacancies", outcome="vacancy")

    # Save
    irf_delta_vac.to_csv(RESULTS_DIR / "lp_irf_delta_vacancy.csv", index=False)
    irf_ld_vac.to_csv(   RESULTS_DIR / "lp_irf_ld_vacancy.csv",    index=False)
    #irf_ts_vac.to_csv(   RESULTS_DIR / "lp_irf_ts_vacancy.csv",    index=False)
    irf_qu_vac.to_csv(   RESULTS_DIR / "lp_irf_qu_vacancy.csv",    index=False)
    print(f"\n  Saved: {RESULTS_DIR / 'lp_irf_delta_vacancy.csv'}")
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_ld_vacancy.csv'}")
    #print(f"  Saved: {RESULTS_DIR / 'lp_irf_ts_vacancy.csv'}")
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_qu_vacancy.csv'}")

    # --- Core asymmetry plot: δ vs LD vacancy IRFs side by side ---
    # This is the key figure for the paper: if δ vacancies fall and LD vacancies
    # are flat/rise, the model's reposting mechanism is empirically validated.
    plot_irf(
        {"δ shock → vacancies": irf_delta_vac,
         "LD shock → vacancies": irf_ld_vac},
        out_path=RESULTS_DIR / "lp_irf_vacancy_delta_ld.png",
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_vacancy_delta_ld.png'}")

    # --- Vacancy IRF overlay: δ vs LD vs QU vs TS on one panel ---
    # Full separation decomposition for vacancy outcome — parallels the
    # unemployment decomposition in Plot 3 of section [7].
    # Model predictions:
    #   δ  → vacancies FALL (firm destruction collapses both unmatched stock + reposting)
    #   LD → vacancies FLAT/RISE (surviving firms repost after layoffs)
    #   QU → ambiguous (quits may not trigger reposting; procyclical)
    #   TS → blend of LD + QU
    overlay_plot_irf(
        {"δ residualized":  irf_delta_vac,
         "LD residualized": irf_ld_vac,
         "QU residualized": irf_qu_vac},
        #  "TS residualized": irf_ts_vac},
        out_path=RESULTS_DIR / "lp_irf_vacancy_decomp_overlay.png",
        title=(
            r"Vacancy IRF: $\delta$ vs LD vs QU vs TS"
            "\n(productivity-residualized, 1-SD scale, post-2001)"
        ),
        ylabel="pp change in vacancy rate per 1-SD shock",
    )

    # --- Unemployment and vacancy IRFs side by side for δ (4-panel figure) ---
    # Combines unemployment (from [4]) and vacancy IRFs to show Beveridge dynamics
    _delta_unemp_irf = irf_delta_resid if _resid_files_present else irf_delta_post
    _ld_unemp_irf    = irf_ld_resid    if _resid_files_present else irf_s
    plot_irf(
        {"δ — unemployment": _delta_unemp_irf,
         "δ — vacancies":    irf_delta_vac,
         "LD — unemployment": _ld_unemp_irf,
         "LD — vacancies":   irf_ld_vac},
        out_path=RESULTS_DIR / "lp_irf_beveridge_asymmetry.png",
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_beveridge_asymmetry.png'}")

    # --- Quick diagnostic: vacancy IRF peaks ---
    print("\n  Vacancy LP summary:")
    for lbl, vdf in [("δ shock → vacancies",  irf_delta_vac),
                     ("LD shock → vacancies", irf_ld_vac),
                   #  ("TS shock → vacancies", irf_ts_vac),
                     ("QU shock → vacancies", irf_qu_vac)]:
        if vdf.empty:
            print(f"    {lbl}: no results (insufficient vacancy obs)")
            continue
        pk  = vdf.loc[vdf["beta"].abs().idxmax()]
        sig = ("***" if pk["pval"] < 0.01 else
               "**"  if pk["pval"] < 0.05 else
               "*"   if pk["pval"] < 0.10 else "(insig)")
        print(f"    {lbl}: peak {pk['beta']:+.4f} log-pts at h={int(pk['h'])}  "
              f"SE={pk['se']:.4f}  p={pk['pval']:.3f}  {sig}")

              