"""
part5c_lp_vacancies.py — Panel Local Projections for δ and s Bartik instruments
              with job openings (vacancies) as outcome variable
===============================================================================
Estimates impulse response functions (IRFs) of state-level job openings
(vacancies) to the job-destruction (δ) and job-finding (s) Bartik shocks via
Jordà (2005) local projections on the state-quarter panel.

This script parallels part5_lp.py but uses job openings (JO, from JOLTS) as
the dependent variable instead of unemployment rates. The model predicts
opposite outcomes across shock types:

    δ shock (establishment exit):
      - Unemployment: Persistent elevation (negative vacancy response)
      - Vacancies: Persistent collapse (reposted vacancies destroyed)

    s shock (match separations):
      - Unemployment: Transient spike (surviving firms repost vacancies)
      - Vacancies: Rapid recovery or overshoot (immediate vacancy repost)

Specification
-------------
For shock k ∈ {δ, s}, horizon h ∈ {0, ..., 16}, estimate separately:

    log(JO_{s,t+h}) - log(JO_{s,t-1}) = α_s + α_t + β_h · B_{s,t}^(k)
                                         + γ₁ · log(JO_{s,t-1})
                                         + γ₂ · u_{s,t-1}
                                         + ε_{s,t,h}

where:
    JO_{s,t}      job openings (level, thousands) in state s, quarter t
    α_s           state fixed effects — absorb permanent cross-state
                  vacancy level differences
    α_t           time fixed effects — absorb the common national cycle
    B_{s,t}^(k)  Bartik instrument for shock type k (already in pp)
    log(JO_{s,t-1})  lagged log job openings — controls for vacancy stock
    u_{s,t-1}    lagged unemployment rate — controls for labor market slack

Identification
--------------
Identical to part5_lp.py: time FE removes the national aggregate shock each
quarter, leaving only cross-state dispersion in Bartik exposure.  β_h measures
how much more (or less) the log job opening changes h quarters after the shock
in a high-exposure state than a low-exposure state.

Standard errors
---------------
Clustered at the state level throughout. With n=50 clusters this is at the
lower bound of reliability; supplement with wild cluster bootstrap p-values
(future work).

Sample periods
--------------
    δ instrument: 1992Q3 – 2019Q4 (outcomes capped to exclude COVID)
    s instrument: 2001Q1 – 2019Q4 (outcomes capped to exclude COVID)
    Joint sample for comparisons: 2001Q1 – 2019Q4
    Both instruments: JOLTS job openings available 2001Q1 onward

Outcomes
-------
    data/results/lp_irf_vacancies_delta.csv       — β_h, SE, CIs for δ shock
    data/results/lp_irf_vacancies_s.csv           — β_h, SE, CIs for s shock
    data/results/lp_irf_vacancies_delta_post2001.csv  — δ on post-2001 sample
    data/results/lp_irf_vacancies_combined.png    — δ and s side by side

Run
---
    python part5c_lp_vacancies.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates.py       (for δ)
    python part2_shock_rates_s.py     (for s)
    python part3_instrument.py        (produces delta_instrument_base2006.csv)
    python part3_instrument_s.py      (produces s_instrument_base2006.csv)
    python part4b_vacancies.py        (produces jolts_jo_quarterly.parquet)
    python part7_nfci.py              (produces nfci_quarterly.parquet)
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
T          = 20
HORIZONS   = list(range(T+1))          # h = 0, 1, ..., 16 quarters
BASE_YEAR  = 2006                     # must match part1/part3 base year
CI_LEVEL   = 0.90                     # confidence band width for main plot
Z90        = 1.645
Z95        = 1.960

# Cap outcome quarter at 2019Q4 to exclude COVID (2020Q1–2021Q4) from the
# h-quarter-ahead outcome window.
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR  = Path("data/instruments")
RESULTS_DIR = Path("data/results")

DELTA_INSTR_FILE       = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
S_INSTR_FILE           = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
JO_OUTCOMES_FILE       = INSTR_DIR / "jolts_jo_quarterly.parquet"
NFCI_FILE              = Path("data/cache") / "nfci_quarterly.parquet"

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
    """Load quarterly JOLTS job openings (vacancies) from part4b."""
    jo = pd.read_parquet(JO_OUTCOMES_FILE)
    jo["state_fips"] = jo["state_fips"].str.zfill(2)
    print(f"  Job openings: {jo.shape}  "
          f"({jo['quarter_label'].min()} – {jo['quarter_label'].max()})")
    return jo[["state_fips", "state", "quarter_label", "job_openings"]]


def load_unemployment() -> pd.DataFrame:
    """Load quarterly LAUS unemployment rates from part4 (for use as control)."""
    laus = pd.read_parquet(INSTR_DIR / "laus_quarterly.parquet")
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    print(f"  Unemployment (control): {laus.shape}  "
          f"({laus['quarter_label'].min()} – {laus['quarter_label'].max()})")
    return laus[["state_fips", "quarter_label", "unemp_rate"]]


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
                outcomes_jo: pd.DataFrame,
                outcomes_u: pd.DataFrame) -> pd.DataFrame:
    """
    Merge instrument with job openings and unemployment, construct controls.

    Controls are lagged to t-1 (one quarter before the shock quarter t):
        log(JO_{s,t-1})   — lagged log job openings
        u_{s,t-1}        — lagged unemployment rate

    Returns a panel DataFrame indexed by (state_fips, quarter_label) with
    columns: state, instr_col, job_openings, unemp_lag1, jo_log_lag1.
    """
    # Merge instrument and outcomes on state × quarter
    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes_jo[["state_fips", "quarter_label", "job_openings"]],
        on=["state_fips", "quarter_label"],
        how="inner",
    )

    # Merge unemployment as control
    panel = panel.merge(
        outcomes_u[["state_fips", "quarter_label", "unemp_rate"]],
        on=["state_fips", "quarter_label"],
        how="left",
    )

    # Rescale instrument from decimal fraction to percentage points.
    panel[instr_col] = panel[instr_col] * 100

    # Sort for lag construction
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)

    # Construct lags within state — shift by 1 quarter
    panel["jo_log_lag1"] = panel.groupby("state_fips")["job_openings"].transform(
        lambda x: np.log(x.shift(1))
    )
    panel["unemp_lag1"] = panel.groupby("state_fips")["unemp_rate"].shift(1)

    # Verify lag is truly t-1 (not a gap across non-contiguous quarters)
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    expected_prev    = quarter_shift(panel["quarter_label"], -1)
    bad_lag          = panel["prev_ql"].notna() & (panel["prev_ql"] != expected_prev)
    if bad_lag.any():
        n_bad = bad_lag.sum()
        print(f"  [warn] {n_bad} rows with non-contiguous quarter sequence — "
              f"lagged controls set to NaN for those rows")
        panel.loc[bad_lag, ["jo_log_lag1", "unemp_lag1"]] = np.nan

    panel = panel.drop(columns="prev_ql")

    n_states = panel["state_fips"].nunique()
    n_qtrs   = panel["quarter_label"].nunique()
    print(f"  Panel: {len(panel):,} rows  ({n_states} states × {n_qtrs} quarters)")
    return panel


# ---------------------------------------------------------------------------
# LP estimation — one horizon
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int,
                   shock_col: str,
                   include_nfci: bool = False) -> dict | None:
    """
    Estimate the LP for a single horizon h.

    Constructs Δlog(JO)_{s,t,h} = log(JO_{s,t+h}) - log(JO_{s,t-1}), merges
    with the instrument at time t, and runs OLS with state + time FEs and
    state-clustered standard errors.

    Returns a dict with point estimate, SEs, CIs, and diagnostics,
    or None if there are too few observations.
    """
    df = base_panel.copy()

    # --- Construct the h-quarter-ahead outcome ---
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql

    outcomes_lookup = (
        df[["state_fips", "quarter_label", "job_openings"]]
        .rename(columns={"quarter_label": "future_ql",
                         "job_openings":  "jo_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")

    # Drop observations where the outcome quarter falls in or after COVID
    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

    # Outcome: log(JO_{s,t+h}) - log(JO_{s,t-1})
    # (baseline is t-1, parallel to unemployment specification)
    df["jo_future_log"] = np.log(df["jo_future"])
    df["dep_var"] = df["jo_future_log"] - df["jo_log_lag1"]

    # Drop rows with missing dep_var or any regressor
    required = ["dep_var", shock_col, "jo_log_lag1", "unemp_lag1"]
    if include_nfci:
        required.append("instr_x_nfci")
    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        print(f"    h={h:2d}: too few obs ({len(df)}) — skipping")
        return None

    # --- Build regressor matrix ---
    state_dummies = pd.get_dummies(df["state_fips"],     prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"],  prefix="qt",
                                   drop_first=True, dtype=float)

    core_regressors = [shock_col, "jo_log_lag1", "unemp_lag1"]
    if include_nfci:
        core_regressors.append("instr_x_nfci")

    X = pd.concat(
        [df[core_regressors],
         state_dummies,
         time_dummies],
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

    # Interaction coefficient δ_h (None when include_nfci=False)
    if include_nfci:
        delta_h    = float(model.params["instr_x_nfci"])
        delta_h_se = float(model.bse["instr_x_nfci"])
        delta_h_t  = float(model.tvalues["instr_x_nfci"])
        delta_h_p  = float(model.pvalues["instr_x_nfci"])
    else:
        delta_h = delta_h_se = delta_h_t = delta_h_p = None

    # Partial F
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
    }


# ---------------------------------------------------------------------------
# LP estimation — all horizons
# ---------------------------------------------------------------------------
def run_lp(base_panel: pd.DataFrame, shock_col: str,
           label: str,
           include_nfci: bool = False) -> pd.DataFrame:
    """
    Run LP for all horizons 0..16 and return results as a DataFrame.
    """
    nfci_tag = " [+NFCI interaction]" if include_nfci else ""
    print(f"\n{'='*60}")
    print(f"  LP — {label}{nfci_tag}  (h = 0 … {max(HORIZONS)})")
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
                             include_nfci=include_nfci)
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

    # --- Standardize to unit-SD scale ---
    instr_sd = float(base_panel[shock_col].std())
    scale_cols = ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]
    irf = pd.DataFrame(rows)
    irf[scale_cols] = irf[scale_cols] * instr_sd
    if include_nfci:
        irf["delta_h"]    = irf["delta_h"]    * instr_sd
        irf["delta_h_se"] = irf["delta_h_se"] * instr_sd
    irf["instr_sd"] = instr_sd
    return irf


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_irf(results: dict[str, pd.DataFrame],
             out_path: Path,
             instr_sds: dict[str, float] | None = None) -> None:
    """
    Side-by-side IRF plots for δ and s shocks with detailed figure caption.
    Outcome variable: log job openings.
    """
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6), sharey=False)
    if n == 1:
        axes = [axes]

    colors = {"δ shock": "#1f77b4", "s shock": "#d62728",
              "δ shock (post-2001)": "#aec7e8"}

    peak_info = {}

    for ax, (label, df) in zip(axes, results.items()):
        if df.empty:
            ax.set_visible(False)
            continue

        h      = df["h"].values
        beta   = df["beta"].values
        lo90   = df["ci90_lo"].values
        hi90   = df["ci90_hi"].values
        lo95   = df["ci95_lo"].values
        hi95   = df["ci95_hi"].values

        color = colors.get(label, "#2ca02c")

        # 95% CI — thin dashed contour
        ax.plot(h, lo95, color=color, linewidth=0.7,
                linestyle="--", alpha=0.6)
        ax.plot(h, hi95, color=color, linewidth=0.7,
                linestyle="--", alpha=0.6)

        # 90% CI — filled band
        ax.fill_between(h, lo90, hi90, color=color, alpha=0.20,
                        label="90% CI (clustered)")

        # Point estimate
        ax.plot(h, beta, color=color, linewidth=2.0,
                marker="o", markersize=3.5, label=label)

        # Zero line
        ax.axhline(0, color="black", linewidth=0.8, linestyle="-")

        # Vertical grid at 4-quarter multiples
        for hh in [4, 8, 12, 16]:
            ax.axvline(hh, color="grey", linewidth=0.5,
                       linestyle=":", alpha=0.6)

        ax.set_title(f"IRF — {label}\n"
                     f"(cumulative log change in job openings from shock quarter)",
                     fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel("log point change in JO per 1 pp shock", fontsize=10)
        ax.set_xticks(h)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: str(int(x))))
        ax.legend(fontsize=9, framealpha=0.85)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)

        # Store peak for caption
        idx_peak = int(np.argmax(np.abs(beta)))
        peak_info[label] = {
            "h":    int(h[idx_peak]),
            "beta": float(beta[idx_peak]),
            "sd":   instr_sds.get(label, np.nan) if instr_sds else np.nan,
        }

    fig.suptitle(
        "Panel LP — Bartik δ and s shocks → log job openings  [Test A: baseline = t−1]\n"
        "(state + time FEs; controls: log(JO_{t-1}), u_{t-1}; "
        "SEs clustered by state; outcomes capped 2019Q4)",
        fontsize=11, y=1.01,
    )

    # -------------------------------------------------------------------
    # Detailed caption
    # -------------------------------------------------------------------
    sd_sentences = []
    for label, info in peak_info.items():
        if not np.isnan(info["sd"]):
            effect_1sd = info["beta"] * info["sd"]
            sd_sentences.append(
                f"A 1-SD move in the {label} instrument ({info['sd']:.2f} pp) "
                f"implies a peak effect of {effect_1sd:.3f} log points at h={info['h']}."
            )

    sd_text = "  ".join(sd_sentences) if sd_sentences else ""

    caption = (
        "Notes: Each panel plots 17 coefficients $\\hat{{\\beta}}_h$, $h=0,\\ldots,16$, "
        "from separate OLS regressions of $\\log(\\mathrm{{JO}}_{{s,t+h}}) - "
        "\\log(\\mathrm{{JO}}_{{s,t-1}})$ (cumulative log change in job openings from "
        "one quarter before the shock) on the Bartik instrument $B_{{s,t}}^{{(k)}}$, "
        "state fixed effects $\\alpha_s$, time fixed effects $\\alpha_t$, "
        "lagged log job openings $\\log(\\mathrm{{JO}}_{{s,t-1}})$, and lagged "
        "unemployment $u_{{s,t-1}}$.  "
        "The baseline is $t-1$ rather than $t$: differencing from the shock quarter "
        "$t$ itself would confound the contemporaneous shock impact with the outcome "
        "variable. Using $t-1$ as the baseline ensures the outcome is predetermined "
        "with respect to $B_{{s,t}}$.  "
        "Outcome quarters are capped at 2019Q4 to exclude COVID (2020Q1–2021Q4), "
        "which produces spurious jumps in $\\hat{{\\beta}}_h$ at the horizons where "
        "COVID outcomes first enter the outcome window.  "
        "The instrument is rescaled to percentage points (×100), so $\\hat{{\\beta}}_h$ "
        "measures the cumulative log change in job openings per "
        "1 pp increase in the Bartik-predicted shock rate.  "
        "The $\\delta$ and $s$ instruments have different cross-sectional standard "
        "deviations (after rescaling): "
        f"SD($B^{{\\delta}}$) $\\approx$ {instr_sds.get('δ shock', float('nan')):.2f} pp, "
        f"SD($B^{{s}}$) $\\approx$ {instr_sds.get('s shock', float('nan')):.2f} pp, "
        "so coefficients are not directly comparable in magnitude across panels.  "
        + (sd_text + "  " if sd_text else "")
        + "Shaded bands: 90\\% CI; dashed lines: 95\\% CI.  "
        "Standard errors clustered by state (50 clusters).  "
        "$\\delta$ instrument sample: 2001Q1–2019Q4 (outcomes, restricted for comparability); "
        "$s$ instrument sample: 2001Q1–2019Q4 (outcomes)."
    )

    fig.text(
        0.5, -0.04, caption,
        ha="center", va="top",
        fontsize=7.5,
        wrap=True,
        transform=fig.transFigure,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                  edgecolor="#cccccc", linewidth=0.8),
        multialignment="left",
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5c — Panel Local Projections (Job Openings Outcome)")
print("=" * 60)

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------
print("\n[1] Loading instruments, outcomes, and NFCI")
delta_instr, s_instr = load_instruments()
outcomes_jo = load_outcomes()
outcomes_u  = load_unemployment()
nfci        = load_nfci()

# -----------------------------------------------------------------------
# Build panels
# -----------------------------------------------------------------------
print("\n[2] Building panels")
print("  δ panel:")
delta_panel = build_panel(delta_instr, "bartik_delta", outcomes_jo, outcomes_u)

print("  s panel:")
s_panel     = build_panel(s_instr,     "bartik_s",     outcomes_jo, outcomes_u)

# δ panel restricted to post-2001Q1 for comparability with s instrument
print("  δ panel (post-2001):")
delta_post2001 = delta_panel[
    delta_panel["quarter_label"] >= "2001Q1"
].copy().reset_index(drop=True)

# -----------------------------------------------------------------------
# Instrument summary statistics
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
print("\n[4] Running panel LPs (job openings outcome)")

irf_delta       = run_lp(delta_post2001, "bartik_delta", "δ shock (post-2001)")
irf_s           = run_lp(s_panel,        "bartik_s",     "s shock")

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
print("\n[5] Saving results")

irf_delta.to_csv(RESULTS_DIR / "lp_irf_vacancies_delta.csv", index=False)
irf_s.to_csv(RESULTS_DIR / "lp_irf_vacancies_s.csv", index=False)

print(f"  Saved: {RESULTS_DIR / 'lp_irf_vacancies_delta.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_vacancies_s.csv'}")

# -----------------------------------------------------------------------
# Plot IRFs
# -----------------------------------------------------------------------
print("\n[6] Plotting IRFs")

# Instrument SDs after ×100 rescaling — used in caption for 1-SD effect
instr_sds_main = {
    "δ shock (post-2001)": float(delta_post2001["bartik_delta"].std()),
    "s shock": float(s_panel["bartik_s"].std()),
}

# Main comparison: δ and s side by side
plot_irf(
    {"δ shock (post-2001)": irf_delta, "s shock": irf_s},
    out_path=RESULTS_DIR / "lp_irf_vacancies_combined.png",
    instr_sds=instr_sds_main,
)

_open_file(RESULTS_DIR / "lp_irf_vacancies_combined.png")

# -----------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------
print("\n[7] Summary — peak vacancy response")
for lbl, df in [("δ (post-2001)", irf_delta), ("s", irf_s)]:
    if df.empty:
        continue
    # Find peak by absolute value
    peak_idx = df["beta"].abs().idxmax()
    peak = df.loc[peak_idx]
    print(f"\n  {lbl} shock peak:")
    print(f"    h = {int(peak['h'])}  β = {peak['beta']:.4f}  "
          f"SE = {peak['se']:.4f}  p = {peak['pval']:.3f}  "
          f"N = {int(peak['nobs'])}")
    print(f"    90% CI: [{peak['ci90_lo']:.4f}, {peak['ci90_hi']:.4f}]")

print("\n" + "=" * 60)
print("Part 5c — Complete")
print("=" * 60)
