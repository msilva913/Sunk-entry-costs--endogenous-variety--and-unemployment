"""
part5d_lp_sloos.py — Panel LP with SLOOS C&I interaction as alternative to NFCI
================================================================================
Extends the baseline panel LP from part5_lp.py by running the same estimation
with B_{s,t} × SLOOS_CI_dm interaction (demeaned within sample) in place of
or alongside the NFCI risk subindex interaction.

Specification
-------------
For shock k ∈ {δ, s}, horizon h ∈ {0, ..., 16}, estimate separately:

    y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h · B_{s,t}^(k)
                             + δ_h · [B_{s,t}^(k) × SLOOS_CI_dm]
                             + γ₁ · u_{s,t-1} + γ₂ · log(LF_{s,t-1})
                             + ε_{s,t,h}

Interaction variable
---------------------
SLOOS_CI_dm is the net percentage of respondents tightening C&I lending
standards (large/middle-market firms) demeaned within the estimation sample.
Demeaning ensures β_h is interpretable at average sample financial conditions,
and δ_h isolates the marginal effect of worse financial conditions.

Comparison with NFCI
--------------------
Both NFCI risk subindex and SLOOS C&I tightening measure financial conditions,
but from different angles:

  NFCI risk subindex
  - Constructed from market prices (volatility, spreads, funding costs)
  - Capturing market participants' risk assessment
  - Includes multiple financial sectors (equity, corporate debt, etc.)
  - Published weekly, aggregated to quarterly here

  SLOOS C&I net tightening
  - Direct survey of bank loan officers
  - Isolating bank credit supply behavior
  - Focused specifically on C&I lending (core to non-financial firms)
  - Published quarterly
  - More direct but potentially subject to survey measurement error

Outputs
-------
    data/results/lp_irf_delta_sloos.csv        — β_h, δ_h, CIs for δ shock
    data/results/lp_irf_s_sloos.csv            — β_h, δ_h, CIs for s shock
    data/results/lp_irf_delta_nfci_vs_sloos.png  — side-by-side δ shock comparison
    data/results/lp_irf_s_nfci_vs_sloos.png      — side-by-side s shock comparison
    data/results/lp_sloos_interaction_table.csv  — δ_h coefficients across all h

Run
---
    python part5d_lp_sloos.py

Prerequisites
-------------
    python part5_lp.py              (produces delta and s panels, NFCI IRFs)
    python part7b_sloos.py          (produces sloos_quarterly.parquet)
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
BASE_YEAR  = 2006
CI_LEVEL   = 0.90
Z90        = 1.645
Z95        = 1.960

MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR  = Path("data/instruments")
RESULTS_DIR = Path("data/results")

DELTA_INSTR_FILE       = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
S_INSTR_FILE           = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
LAUS_FILE              = INSTR_DIR / "laus_quarterly.parquet"
SLOOS_FILE             = Path("data/cache") / "sloos_quarterly.parquet"

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


def load_sloos() -> pd.DataFrame:
    """
    Load quarterly SLOOS data produced by part7b_sloos.py.
    Returns DataFrame with columns: quarter_label, sloos_ci, sloos_ci_small.
    """
    sloos = pd.read_parquet(SLOOS_FILE)
    print(f"  SLOOS: {sloos.shape}  "
          f"({sloos['quarter_label'].min()} – {sloos['quarter_label'].max()})")
    print(f"  sloos_ci (large/middle) — mean={sloos['sloos_ci'].mean():.3f}  "
          f"SD={sloos['sloos_ci'].std():.3f}  "
          f"min={sloos['sloos_ci'].min():.3f}  "
          f"max={sloos['sloos_ci'].max():.3f}")
    return sloos[["quarter_label", "sloos_ci", "sloos_ci_small"]]


# ---------------------------------------------------------------------------
# Panel construction
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame) -> pd.DataFrame:
    """
    Merge instrument with outcomes and construct pre-determined controls.
    (Same as part5_lp.py)
    """
    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[["state_fips", "quarter_label", "unemp_rate", "labor_force"]],
        on=["state_fips", "quarter_label"],
        how="inner",
    )

    panel[instr_col] = panel[instr_col] * 100

    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)

    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )

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
# Attach SLOOS interaction
# ---------------------------------------------------------------------------
def attach_sloos_interaction(panel: pd.DataFrame,
                             sloos: pd.DataFrame,
                             instr_col: str) -> pd.DataFrame:
    """
    Merge quarterly SLOOS into panel and construct the demeaned interaction
    term B_{s,t} × SLOOS_CI_dm.

    Demeaning is done within the panel's own estimation sample — i.e. the
    mean of sloos_ci is computed over the quarters actually present in the
    panel after merging. This ensures that β_h is interpretable at average
    financial conditions prevailing during the estimation window.

    Two columns are added:
        sloos_ci_dm      : demeaned SLOOS C&I (sloos_ci − mean)
        instr_x_sloos    : B_{s,t} × sloos_ci_dm  (the interaction regressor)
    """
    panel = panel.merge(
        sloos[["quarter_label", "sloos_ci", "sloos_ci_small"]],
        on="quarter_label",
        how="left",
    )

    # Warn if any quarters in the panel have no SLOOS match
    n_missing = panel["sloos_ci"].isna().sum()
    if n_missing > 0:
        missing_qtrs = (panel.loc[panel["sloos_ci"].isna(), "quarter_label"]
                        .unique().tolist())
        print(f"  [warn] {n_missing} panel rows have no SLOOS match "
              f"({len(missing_qtrs)} quarters): {missing_qtrs[:5]} ...")

    # Demean within the panel's own sample
    sloos_mean = panel["sloos_ci"].mean()
    sloos_sd   = panel["sloos_ci"].std()
    panel["sloos_ci_dm"] = panel["sloos_ci"] - sloos_mean

    print(f"  {instr_col}: sloos_ci sample mean = {sloos_mean:.4f}  "
          f"SD = {sloos_sd:.4f}  "
          f"(demeaned range: [{panel['sloos_ci_dm'].min():.3f}, "
          f"{panel['sloos_ci_dm'].max():.3f}])")

    # Interaction: instrument (already in pp) × demeaned SLOOS
    panel["instr_x_sloos"] = panel[instr_col] * panel["sloos_ci_dm"]

    return panel


# ---------------------------------------------------------------------------
# LP estimation — one horizon
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int,
                   shock_col: str) -> dict | None:
    """
    Estimate the LP for a single horizon h with SLOOS interaction.
    (Parallel to run_lp_horizon in part5_lp.py but with instr_x_sloos.)
    """
    df = base_panel.copy()

    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql

    outcomes_lookup = (
        df[["state_fips", "quarter_label", "unemp_rate"]]
        .rename(columns={"quarter_label": "future_ql",
                         "unemp_rate":    "y_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")

    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

    df["dep_var"] = df["y_future"] - df["unemp_lag1"]

    required = ["dep_var", shock_col, "unemp_lag1", "lf_log_lag1", "instr_x_sloos"]
    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        print(f"    h={h:2d}: too few obs ({len(df)}) — skipping")
        return None

    # State and time dummies
    state_dummies = pd.get_dummies(df["state_fips"],     prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"],  prefix="qt",
                                   drop_first=True, dtype=float)

    core_regressors = [shock_col, "unemp_lag1", "lf_log_lag1", "instr_x_sloos"]

    X = pd.concat(
        [df[core_regressors],
         state_dummies,
         time_dummies],
        axis=1,
    )
    X = sm.add_constant(X, has_constant="add")
    y = df["dep_var"]

    model = sm.OLS(y, X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values},
    )

    beta  = float(model.params[shock_col])
    se    = float(model.bse[shock_col])
    tstat = float(model.tvalues[shock_col])
    pval  = float(model.pvalues[shock_col])

    # SLOOS interaction coefficient
    gamma_h    = float(model.params["instr_x_sloos"])
    gamma_h_se = float(model.bse["instr_x_sloos"])
    gamma_h_t  = float(model.tvalues["instr_x_sloos"])
    gamma_h_p  = float(model.pvalues["instr_x_sloos"])

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
        # SLOOS interaction coefficient
        "gamma_h":     gamma_h,
        "gamma_h_se":  gamma_h_se,
        "gamma_h_t":   gamma_h_t,
        "gamma_h_p":   gamma_h_p,
        "nobs":        int(model.nobs),
        "r2":          float(model.rsquared),
        "n_clusters":  df["state_fips"].nunique(),
        "qt_range":    f"{df['quarter_label'].min()}–{df['quarter_label'].max()}",
    }


# ---------------------------------------------------------------------------
# LP estimation — all horizons
# ---------------------------------------------------------------------------
def run_lp(base_panel: pd.DataFrame, shock_col: str,
           label: str) -> pd.DataFrame:
    """
    Run LP for all horizons 0..16 with SLOOS interaction.
    """
    print(f"\n{'='*60}")
    print(f"  LP — {label}  (h = 0 … {max(HORIZONS)}) [with SLOOS interaction]")
    print(f"{'='*60}")
    print(f"  {'h':>3}  {'β_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}  "
          f"{'γ_h':>10}  {'γ_h SE':>8}  {'γ_h p':>6}  {'N':>6}")
    print(f"  {'-'*60}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(base_panel, h, shock_col)
        if res is None:
            continue
        rows.append(res)
        sig = ("***" if res["pval"] < 0.01 else
               "**"  if res["pval"] < 0.05 else
               "*"   if res["pval"] < 0.10 else "")
        gsig = ("***" if res["gamma_h_p"] < 0.01 else
                "**"  if res["gamma_h_p"] < 0.05 else
                "*"   if res["gamma_h_p"] < 0.10 else "")
        print(f"  {h:3d}  {res['beta']:10.4f}  {res['se']:8.4f}  "
              f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
              f"{res['gamma_h']:10.4f}  {res['gamma_h_se']:8.4f}  "
              f"{res['gamma_h_p']:6.3f}{gsig}  {res['nobs']:6d}  {sig}")

    # ── Standardize to unit-SD scale ───────────────────────────────────────
    instr_sd = float(base_panel[shock_col].std())
    scale_cols = ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]
    irf = pd.DataFrame(rows)
    irf[scale_cols] = irf[scale_cols] * instr_sd
    # Also scale the interaction coefficient
    irf["gamma_h"]    = irf["gamma_h"]    * instr_sd
    irf["gamma_h_se"] = irf["gamma_h_se"] * instr_sd
    irf["instr_sd"] = instr_sd
    return irf


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_irf_comparison(results_nfci: pd.DataFrame,
                        results_sloos: pd.DataFrame,
                        label: str,
                        out_path: Path) -> None:
    """
    Side-by-side plot: NFCI interaction (from part5_lp.py) vs SLOOS interaction.
    Shows main effect β_h and interaction effect δ_h (NFCI) vs γ_h (SLOOS).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    colors_nfci  = {"main": "#1f77b4", "int": "#aec7e8"}
    colors_sloos = {"main": "#d62728", "int": "#ff9999"}

    # Left: NFCI
    h_nfci = results_nfci["h"].values
    beta_nfci = results_nfci["beta"].values
    lo90_nfci = results_nfci["ci90_lo"].values
    hi90_nfci = results_nfci["ci90_hi"].values
    delta_h_nfci = results_nfci["delta_h"].values

    ax = axes[0]
    ax.fill_between(h_nfci, lo90_nfci, hi90_nfci, color=colors_nfci["main"],
                   alpha=0.20, label="90% CI main effect")
    ax.plot(h_nfci, beta_nfci, color=colors_nfci["main"], linewidth=2.0,
           marker="o", markersize=3.5, label="β_h (main effect)")
    ax.plot(h_nfci, delta_h_nfci, color=colors_nfci["int"], linewidth=2.0,
           marker="s", markersize=3.5, linestyle="--",
           label="δ_h (NFCI interaction)")
    ax.axhline(0, color="black", linewidth=0.8, linestyle="-")
    for hh in [4, 8, 12, 16]:
        ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.set_title(f"NFCI Risk Subindex Interaction\n({label})", fontsize=10)
    ax.set_xlabel("Horizon h (quarters)", fontsize=9)
    ax.set_ylabel("pp change in unemp. rate per 1 pp shock", fontsize=9)
    ax.set_xticks(h_nfci)
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)

    # Right: SLOOS
    h_sloos = results_sloos["h"].values
    beta_sloos = results_sloos["beta"].values
    lo90_sloos = results_sloos["ci90_lo"].values
    hi90_sloos = results_sloos["ci90_hi"].values
    gamma_h_sloos = results_sloos["gamma_h"].values

    ax = axes[1]
    ax.fill_between(h_sloos, lo90_sloos, hi90_sloos, color=colors_sloos["main"],
                   alpha=0.20, label="90% CI main effect")
    ax.plot(h_sloos, beta_sloos, color=colors_sloos["main"], linewidth=2.0,
           marker="o", markersize=3.5, label="β_h (main effect)")
    ax.plot(h_sloos, gamma_h_sloos, color=colors_sloos["int"], linewidth=2.0,
           marker="s", markersize=3.5, linestyle="--",
           label="γ_h (SLOOS interaction)")
    ax.axhline(0, color="black", linewidth=0.8, linestyle="-")
    for hh in [4, 8, 12, 16]:
        ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.set_title(f"SLOOS C&I Tightening Interaction\n({label})", fontsize=10)
    ax.set_xlabel("Horizon h (quarters)", fontsize=9)
    ax.set_xticks(h_sloos)
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)

    caption = (
        "Notes: Left panel shows results with NFCI risk subindex interaction (from part5_lp.py). "
        "Right panel shows results with SLOOS C&I net tightening interaction. Both measure "
        "financial conditions but from different angles: NFCI is constructed from market prices "
        "(volatility, spreads); SLOOS is direct survey of bank lending behavior. β_h (main effect) "
        "is the impulse response of unemployment to the Bartik instrument at typical financial "
        "conditions (demeaned). δ_h (NFCI) and γ_h (SLOOS) capture how the response changes when "
        "financial conditions are tighter (positive = amplification of unemployment response). "
        "Blue (NFCI) and red (SLOOS) lines allow visual comparison of the two interaction specifications."
    )
    fig.text(0.5, -0.12, caption, ha="center", va="top", fontsize=7.5,
             wrap=True, transform=fig.transFigure,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                       edgecolor="#cccccc", linewidth=0.8),
             multialignment="left")

    fig.suptitle(
        f"Financial Conditions Interaction: NFCI vs. SLOOS\n{label}",
        fontsize=11, y=1.00,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5d — Panel LP with SLOOS Interaction")
print("=" * 60)

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------
print("\n[1] Loading instruments, outcomes, and SLOOS")
delta_instr, s_instr = load_instruments()
outcomes = load_outcomes()
sloos    = load_sloos()

# -----------------------------------------------------------------------
# Load baseline NFCI results for comparison
# -----------------------------------------------------------------------
print("\n[1b] Loading baseline NFCI results from part5_lp.py")
nfci_delta_file = RESULTS_DIR / "lp_irf_delta_nfci.csv"
nfci_s_file     = RESULTS_DIR / "lp_irf_s_nfci.csv"

if not nfci_delta_file.exists() or not nfci_s_file.exists():
    print("  ERROR: NFCI results not found. Run part5_lp.py first.")
    print(f"    Expected: {nfci_delta_file}")
    print(f"    Expected: {nfci_s_file}")
    import sys
    sys.exit(1)

irf_delta_nfci = pd.read_csv(nfci_delta_file)
irf_s_nfci     = pd.read_csv(nfci_s_file)
print(f"  Loaded NFCI δ results: {len(irf_delta_nfci)} horizons")
print(f"  Loaded NFCI s results: {len(irf_s_nfci)} horizons")

# -----------------------------------------------------------------------
# Build panels
# -----------------------------------------------------------------------
print("\n[2] Building panels")
print("  δ panel:")
delta_panel = build_panel(delta_instr, "bartik_delta", outcomes)

print("  s panel:")
s_panel     = build_panel(s_instr,     "bartik_s",     outcomes)

# δ panel restricted to post-2001Q1 for comparability with s
print("  δ panel (post-2001):")
delta_post2001 = delta_panel[
    delta_panel["quarter_label"] >= "2001Q1"
].copy().reset_index(drop=True)

# -----------------------------------------------------------------------
# Merge SLOOS and construct demeaned interaction term
# -----------------------------------------------------------------------
print("\n[2b] Merging SLOOS and constructing interaction term")

delta_panel    = attach_sloos_interaction(delta_panel,    sloos, "bartik_delta")
s_panel        = attach_sloos_interaction(s_panel,        sloos, "bartik_s")
delta_post2001 = attach_sloos_interaction(delta_post2001, sloos, "bartik_delta")
n_post = delta_post2001["quarter_label"].nunique()
print(f"    Restricted to 2001Q1+: {len(delta_post2001):,} rows "
      f"({delta_post2001['state_fips'].nunique()} states × {n_post} quarters)")

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
# Run LPs with SLOOS interaction
# -----------------------------------------------------------------------
print("\n[4] Running panel LPs with SLOOS interaction")

irf_delta_sloos = run_lp(delta_panel,    "bartik_delta", "δ shock")
irf_s_sloos     = run_lp(s_panel,        "bartik_s",     "s shock")

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
print("\n[5] Saving results")

irf_delta_sloos.to_csv(RESULTS_DIR / "lp_irf_delta_sloos.csv", index=False)
irf_s_sloos.to_csv(RESULTS_DIR / "lp_irf_s_sloos.csv", index=False)

print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_sloos.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_sloos.csv'}")

# -----------------------------------------------------------------------
# Create comparison table: interaction coefficients (NFCI vs SLOOS)
# -----------------------------------------------------------------------
print("\n[6] Creating interaction coefficient comparison table")

comparison_delta = pd.DataFrame({
    "h":            irf_delta_sloos["h"],
    "delta_h_nfci": irf_delta_nfci["delta_h"],
    "delta_h_nfci_se": irf_delta_nfci["delta_h_se"],
    "delta_h_nfci_p":  irf_delta_nfci["delta_h_p"],
    "gamma_h_sloos": irf_delta_sloos["gamma_h"],
    "gamma_h_sloos_se": irf_delta_sloos["gamma_h_se"],
    "gamma_h_sloos_p":  irf_delta_sloos["gamma_h_p"],
})

comparison_s = pd.DataFrame({
    "h":            irf_s_sloos["h"],
    "delta_h_nfci": irf_s_nfci["delta_h"],
    "delta_h_nfci_se": irf_s_nfci["delta_h_se"],
    "delta_h_nfci_p":  irf_s_nfci["delta_h_p"],
    "gamma_h_sloos": irf_s_sloos["gamma_h"],
    "gamma_h_sloos_se": irf_s_sloos["gamma_h_se"],
    "gamma_h_sloos_p":  irf_s_sloos["gamma_h_p"],
})

comparison_delta.to_csv(RESULTS_DIR / "lp_comparison_delta_nfci_vs_sloos.csv",
                       index=False)
comparison_s.to_csv(RESULTS_DIR / "lp_comparison_s_nfci_vs_sloos.csv",
                   index=False)

print(f"  Saved: {RESULTS_DIR / 'lp_comparison_delta_nfci_vs_sloos.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_comparison_s_nfci_vs_sloos.csv'}")

# -----------------------------------------------------------------------
# Plotting — comparison
# -----------------------------------------------------------------------
print("\n[7] Plotting NFCI vs SLOOS comparisons")

plot_irf_comparison(irf_delta_nfci, irf_delta_sloos,
                   "δ shock",
                   out_path=RESULTS_DIR / "lp_irf_delta_nfci_vs_sloos.png")

plot_irf_comparison(irf_s_nfci, irf_s_sloos,
                   "s shock",
                   out_path=RESULTS_DIR / "lp_irf_s_nfci_vs_sloos.png")

_open_file(RESULTS_DIR / "lp_irf_delta_nfci_vs_sloos.png")

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
print("\n[8] Summary — peak interaction effects")
for lbl, df_nfci, df_sloos in [("δ shock", irf_delta_nfci, irf_delta_sloos),
                                 ("s shock", irf_s_nfci,     irf_s_sloos)]:
    if df_nfci.empty or df_sloos.empty:
        continue

    peak_nfci = df_nfci.loc[df_nfci["delta_h"].abs().idxmax()]
    peak_sloos = df_sloos.loc[df_sloos["gamma_h"].abs().idxmax()]

    print(f"\n  {lbl}")
    print(f"    NFCI:  h={int(peak_nfci['h'])} δ_h={peak_nfci['delta_h']:8.4f} "
          f"SE={peak_nfci['delta_h_se']:8.4f} p={peak_nfci['delta_h_p']:.3f}")
    print(f"    SLOOS: h={int(peak_sloos['h'])} γ_h={peak_sloos['gamma_h']:8.4f} "
          f"SE={peak_sloos['gamma_h_se']:8.4f} p={peak_sloos['gamma_h_p']:.3f}")

print("\n" + "=" * 60)
print("Part 5d complete — Results ready for interpretation")
print("=" * 60)
