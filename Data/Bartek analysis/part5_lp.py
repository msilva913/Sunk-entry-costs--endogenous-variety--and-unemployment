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
    data/results/lp_irf_s.csv                 — β_h, SE, CIs for s shock
    data/results/lp_irf_delta_post2001.csv    — δ on post-2001 sample
    data/results/lp_irf_delta_nfci.csv        — δ with NFCI interaction
    data/results/lp_irf_s_nfci.csv            — s with NFCI interaction
    data/results/lp_irf_delta_resid.csv       — δ residualized (post-2001)
    data/results/lp_irf_s_resid.csv           — s residualized
    data/results/lp_irf_combined.png          — raw δ and s side by side
    data/results/lp_irf_delta_resid_comparison.png  — δ raw vs. residualized
    data/results/lp_irf_s_resid_comparison.png      — s raw vs. residualized

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
HORIZONS   = list(range(17))          # h = 0, 1, ..., 16 quarters
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
    """Load quarterly LAUS unemployment rates and labor force from part4."""
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    print(f"  LAUS: {laus.shape}  "
          f"({laus['quarter_label'].min()} – {laus['quarter_label'].max()})")
    return laus[["state_fips", "state", "quarter_label",
                 "unemp_rate", "labor_force"]]


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

    Returns a panel DataFrame indexed by (state_fips, quarter_label) with
    columns: state, instr_col, unemp_rate, unemp_lag1, lf_log_lag1.
    """
    # Merge instrument and outcomes on state × quarter
    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[["state_fips", "quarter_label", "unemp_rate", "labor_force"]],
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

    # Verify lag is truly t-1 (not a gap across non-contiguous quarters)
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    expected_prev    = quarter_shift(panel["quarter_label"], -1)
    bad_lag          = panel["prev_ql"].notna() & (panel["prev_ql"] != expected_prev)
    if bad_lag.any():
        n_bad = bad_lag.sum()
        print(f"  [warn] {n_bad} rows with non-contiguous quarter sequence — "
              f"lagged controls set to NaN for those rows")
        panel.loc[bad_lag, ["unemp_lag1", "lf_log_lag1"]] = np.nan

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

    Constructs Δy_{s,t,h} = y_{s,t+h} - y_{s,t-1}, merges with the
    instrument at time t, and runs OLS with state + time FEs and
    state-clustered standard errors.

    Returns a dict with point estimate, SEs, CIs, and diagnostics,
    or None if there are too few observations.
    """
    df = base_panel.copy()

    # --- Construct the h-quarter-ahead outcome ---
    # For each (state, t), look up y_{s,t+h} by shifting the quarter label
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql

    outcomes_lookup = (
        df[["state_fips", "quarter_label", "unemp_rate"]]
        .rename(columns={"quarter_label": "future_ql",
                         "unemp_rate":    "y_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")

    # Drop observations where the outcome quarter falls in or after the COVID
    # window. COVID (2020Q1–2021Q4) caused an aggregate unemployment spike
    # completely unrelated to the Bartik instrument, contaminating β_h at
    # the horizons where those quarters first enter the outcome window.
    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

    # Outcome: y_{s,t+h} - y_{s,t-1}  (Test A — baseline is t-1, not t)
    #
    # Motivation: y_{s,t} is contaminated by the shock itself — states with
    # a large B_{s,t} already have elevated unemployment within quarter t.
    # Differencing from t therefore understates the early response and
    # mechanically produces upward drift in β_h as the cumulative effect
    # accumulates above an already-elevated starting point.
    #
    # Using t-1 as baseline avoids this: unemployment one quarter before
    # the shock is predetermined with respect to B_{s,t}. The cost is
    # near-collinearity between the baseline and u_{s,t-1} (they differ
    # only by one quarter), but this affects efficiency not identification.
    # β_0 is no longer anchored at zero by construction — it now measures
    # the within-quarter impact of the shock on unemployment.
    df["dep_var"] = df["y_future"] - df["unemp_lag1"]

    # Drop rows with missing dep_var or any regressor
    required = ["dep_var", shock_col, "unemp_lag1", "lf_log_lag1"]
    if include_nfci:
        required.append("instr_x_nfci")
    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        print(f"    h={h:2d}: too few obs ({len(df)}) — skipping")
        return None

    # --- Build regressor matrix ---
    # State and time dummies (drop_first avoids perfect multicollinearity)
    state_dummies = pd.get_dummies(df["state_fips"],     prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"],  prefix="qt",
                                   drop_first=True, dtype=float)

    core_regressors = [shock_col, "unemp_lag1", "lf_log_lag1"]
    if include_nfci:
        # Note: the main effect of nfci_risk_dm is absorbed by time FEs
        # (it varies only over time, not across states within a quarter),
        # so only the interaction term enters as a cross-sectional regressor.
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
    # Clustering allows arbitrary within-state serial correlation, which
    # is expected at long horizons due to overlapping observations
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

    # Partial F: squared t-stat on the instrument (equals partial F
    # with one instrument under clustered SEs)
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
        # Interaction coefficient δ_h (None when include_nfci=False)
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

    Args:
        base_panel:   output of build_panel() with attach_nfci_interaction()
        shock_col:    'bartik_delta' or 'bartik_s'
        label:        human-readable label for printing (e.g. 'δ shock')
        include_nfci: if True, include B_{s,t} × NFCI_risk_dm interaction.
                      The printed table gains a δ_h column.
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

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_irf(results: dict[str, pd.DataFrame],
             out_path: Path,
             instr_sds: dict[str, float] | None = None) -> None:
    """
    Side-by-side IRF plots for δ and s shocks with detailed figure caption.

    Shaded bands show 90% confidence interval (clustered SEs).
    A dashed 95% CI contour is added as a thin outer line.
    Zero line and h=4/8/12/16 grid lines aid readability.

    results:    dict mapping label → DataFrame from run_lp()
    instr_sds:  dict mapping label → cross-sectional SD of instrument (in pp,
                i.e. after ×100 rescaling) — used in caption for 1-SD scaling
    """
    n = len(results)
    # Extra bottom margin for caption
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6), sharey=False)
    if n == 1:
        axes = [axes]

    colors = {"δ shock": "#1f77b4", "s shock": "#d62728",
              "δ shock (post-2001)": "#aec7e8"}

    peak_info = {}   # collect peak β per label for caption

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
                     f"(cumulative change in unemp. rate from shock quarter, pp)",
                     fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel("pp change in unemp. rate per 1 pp shock", fontsize=10)
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
        "Panel LP — Bartik δ and s shocks → unemployment rate  [Test A: baseline = t−1]\n"
        "(state + time FEs; controls: u_{t-1}, log LF_{t-1}; "
        "SEs clustered by state; outcomes capped 2019Q4)",
        fontsize=11, y=1.01,
    )

    # -------------------------------------------------------------------
    # Detailed caption
    # -------------------------------------------------------------------
    # Build 1-SD effect sentences for each shock if SD info available
    sd_sentences = []
    for label, info in peak_info.items():
        if not np.isnan(info["sd"]):
            effect_1sd = info["beta"] * info["sd"]
            sd_sentences.append(
                f"A 1-SD move in the {label} instrument ({info['sd']:.2f} pp) "
                f"implies a peak effect of {effect_1sd:.2f} pp at h={info['h']}."
            )

    sd_text = "  ".join(sd_sentences) if sd_sentences else ""

    caption = (
        "Notes: Each panel plots 17 coefficients $\\hat{{\\beta}}_h$, $h=0,\\ldots,16$, "
        "from separate OLS regressions of $y_{{s,t+h}} - y_{{s,t-1}}$ (cumulative "
        "change in the unemployment rate from one quarter before the shock) on the "
        "Bartik instrument $B_{{s,t}}^{{(k)}}$, state fixed effects $\\alpha_s$, time "
        "fixed effects $\\alpha_t$, lagged unemployment $u_{{s,t-1}}$, and lagged log "
        "labor force $\\log(\\mathrm{{LF}}_{{s,t-1}})$.  "
        "The baseline is $t-1$ rather than $t$: differencing from the shock quarter "
        "$t$ itself is problematic because $y_{{s,t}}$ is contaminated by the "
        "contemporaneous shock — states with larger $B_{{s,t}}$ already have elevated "
        "unemployment within quarter $t$, causing the difference $y_{{s,t+h}}-y_{{s,t}}$ "
        "to understate the early response and drift upward mechanically at long horizons. "
        "Using $t-1$ as the baseline avoids this: unemployment one quarter before the "
        "shock is predetermined with respect to $B_{{s,t}}$. Unlike the $y_{{s,t}}$ "
        "baseline, $\\hat{{\\beta}}_0$ is not anchored at zero — it measures the "
        "within-quarter impact of the shock.  "
        "Outcome quarters are capped at 2019Q4 to exclude COVID (2020Q1–2021Q4), "
        "which produces spurious jumps in $\\hat{{\\beta}}_h$ at the horizons where "
        "COVID outcomes first enter the outcome window.  "
        "The instrument is rescaled to percentage points (×100), so $\\hat{{\\beta}}_h$ "
        "measures the cumulative change in the unemployment rate (in pp) per "
        "1 pp increase in the Bartik-predicted shock rate.  "
        "The $\\delta$ and $s$ instruments have different cross-sectional standard "
        "deviations (after rescaling: "
        f"SD($B^{{\\delta}}$) $\\approx$ {instr_sds.get('δ shock', float('nan')):.2f} pp, "
        f"SD($B^{{s}}$) $\\approx$ {instr_sds.get('s shock', float('nan')):.2f} pp), "
        "so coefficients are not directly comparable in magnitude across panels.  "
        + (sd_text + "  " if sd_text else "")
        + "Shaded bands: 90\\% CI; dashed lines: 95\\% CI.  "
        "Standard errors clustered by state (50 clusters).  "
        "$\\delta$ instrument sample: 1992Q3–2019Q4 (outcomes); "
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
        # Use a wide text width so caption wraps naturally
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
# Residualized instruments — LP (if available)
# -----------------------------------------------------------------------
if DELTA_RESID_INSTR_FILE.exists() and S_RESID_INSTR_FILE.exists():
    print("\n[7] Running LPs with productivity-residualized instruments")

    delta_resid_raw = pd.read_csv(DELTA_RESID_INSTR_FILE,
                                   dtype={"state_fips": str})
    s_resid_raw     = pd.read_csv(S_RESID_INSTR_FILE,
                                   dtype={"state_fips": str})
    delta_resid_raw["state_fips"] = delta_resid_raw["state_fips"].str.zfill(2)
    s_resid_raw["state_fips"]     = s_resid_raw["state_fips"].str.zfill(2)

    # Residualized series are demeaned so the column name is the same;
    # build_panel expects the instrument column to be named bartik_delta/bartik_s.
    outcomes = load_outcomes()
    delta_resid_panel = build_panel(delta_resid_raw, "bartik_delta", outcomes)
    s_resid_panel     = build_panel(s_resid_raw,     "bartik_s",     outcomes)

    # Restrict δ residualized to post-2001 for direct comparability with s.
    delta_resid_post = delta_resid_panel[
        delta_resid_panel["quarter_label"] >= "2001Q1"
    ].copy()

    irf_delta_resid = run_lp(delta_resid_post, "bartik_delta",
                              "δ shock (resid, post-2001)")
    irf_s_resid     = run_lp(s_resid_panel,    "bartik_s",
                              "s shock (resid)")

    irf_delta_resid.to_csv(RESULTS_DIR / "lp_irf_delta_resid.csv", index=False)
    irf_s_resid.to_csv(RESULTS_DIR     / "lp_irf_s_resid.csv",     index=False)
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_resid.csv'}")
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_resid.csv'}")

    # Comparison plot: raw (post-2001) vs. residualized, δ and s side by side
    instr_sds_resid = {
        "δ raw (post-2001)":  float(delta_post2001["bartik_delta"].std()),
        "δ residualized":     float(delta_resid_post["bartik_delta"].std()),
        "s raw":              float(s_panel["bartik_s"].std()),
        "s residualized":     float(s_resid_panel["bartik_s"].std()),
    }

    # δ comparison
    plot_irf(
        {"δ raw (post-2001)": irf_delta_post,
         "δ residualized":    irf_delta_resid},
        out_path=RESULTS_DIR / "lp_irf_delta_resid_comparison.png",
        instr_sds=instr_sds_resid,
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_resid_comparison.png'}")

    # s comparison
    plot_irf(
        {"s raw":         irf_s,
         "s residualized": irf_s_resid},
        out_path=RESULTS_DIR / "lp_irf_s_resid_comparison.png",
        instr_sds=instr_sds_resid,
    )
    print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_resid_comparison.png'}")

    _open_file(RESULTS_DIR / "lp_irf_delta_resid_comparison.png")
    _open_file(RESULTS_DIR / "lp_irf_s_resid_comparison.png")

else:
    print("\n[7] Residualized instruments not found — skipping resid LP.")
    print("    Run part3_resid_instruments.py first.")

# -----------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------
print("\n[8] Summary — peak unemployment response")
for lbl, df in [("δ", irf_delta), ("s", irf_s)]:
    if df.empty:
        continue
    peak = df.loc[df["beta"].idxmax()]
    print(f"\n  {lbl} shock peak:")
    print(f"    h = {int(peak['h'])}  β = {peak['beta']:.4f}  "
          f"SE = {peak['se']:.4f}  p = {peak['pval']:.3f}  "
          f"N = {int(peak['nobs'])}")
    print(f"    90% CI: [{peak['ci90_lo']:.4f}, {peak['ci90_hi']:.4f}]")