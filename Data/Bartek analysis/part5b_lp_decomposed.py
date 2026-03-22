"""
part5b_lp_decomposed.py — JOLTS Decomposition LP: LD vs Quits vs Total Separations
====================================================================================
PRIORITY 2 from the project brief.

Constructs separate Bartik instruments from:
  1. Total separations (TS) — baseline, for comparability
  2. Layoffs and discharges (LD) — primary, closest to model's s_t
  3. Quits (QU) — placebo, model predicts different IRF profile

Then runs the LP specification from part5_lp.py with each instrument and
produces comparison IRF plots.

Rationale
---------
The baseline s Bartik instrument conflates three types of separations:
  - Quits: procyclical, no model counterpart (no on-the-job search)
  - Temporary layoffs: recall arrangements bypass reposting channel
  - Permanent discharges: cleanest analog to model's s_t

This decomposition is a first-order identification concern.  Expected:
  - B^LD and B^δ moderately correlated (both countercyclical)
  - B^QU and B^δ low/negative correlation (opposite cyclicality)
  - LD IRF should show a more moderate/declining profile than total s
  - QU IRF should differ qualitatively from LD (procyclical shock)

Outputs
-------
  data/instruments/ld_instrument_base2006.csv
  data/instruments/qu_instrument_base2006.csv
  data/results/lp_irf_ld.csv
  data/results/lp_irf_qu.csv
  data/results/lp_irf_decomposition_comparison.png
  data/results/instrument_correlations_decomposed.csv

Run
---
    python part5b_lp_decomposed.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates_s.py   (produces LD and QU shock rates)
    python part4_outcomes.py
    python part7_nfci.py
"""

import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    pass

# ---------------------------------------------------------------------------
# Configuration  (mirrors part5_lp.py)
# ---------------------------------------------------------------------------
T          = 20
HORIZONS   = list(range(T + 1))
BASE_YEAR  = 2006
CI_LEVEL   = 0.90
Z90        = 1.645
Z95        = 1.960

MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
CACHE_DIR   = Path("data/cache")

SHARES_PATH       = INSTR_DIR / f"shares_base{BASE_YEAR}.parquet"
S_INSTR_FILE      = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
DELTA_INSTR_FILE  = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
LAUS_FILE         = INSTR_DIR / "laus_quarterly.parquet"
NFCI_FILE         = CACHE_DIR / "nfci_quarterly.parquet"

SHOCK_RATES_LD    = INSTR_DIR / "shock_rates_ld_2001Q1_2023Q1.parquet"
SHOCK_RATES_QU    = INSTR_DIR / "shock_rates_qu_2001Q1_2023Q1.parquet"

LD_INSTR_FILE     = INSTR_DIR / f"ld_instrument_base{BASE_YEAR}.csv"
QU_INSTR_FILE     = INSTR_DIR / f"qu_instrument_base{BASE_YEAR}.csv"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FIPS2D_TO_STATE = {
    "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA",
    "08": "CO", "09": "CT", "10": "DE", "12": "FL", "13": "GA",
    "15": "HI", "16": "ID", "17": "IL", "18": "IN", "19": "IA",
    "20": "KS", "21": "KY", "22": "LA", "23": "ME", "24": "MD",
    "25": "MA", "26": "MI", "27": "MN", "28": "MS", "29": "MO",
    "30": "MT", "31": "NE", "32": "NV", "33": "NH", "34": "NJ",
    "35": "NM", "36": "NY", "37": "NC", "38": "ND", "39": "OH",
    "40": "OK", "41": "OR", "42": "PA", "44": "RI", "45": "SC",
    "46": "SD", "47": "TN", "48": "TX", "49": "UT", "50": "VT",
    "51": "VA", "53": "WA", "54": "WV", "55": "WI", "56": "WY",
}


# ---------------------------------------------------------------------------
# Quarter arithmetic  (identical to part5_lp.py)
# ---------------------------------------------------------------------------
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


# ---------------------------------------------------------------------------
# Bartik aggregation for LD and QU
# ---------------------------------------------------------------------------
def build_bartik(shares: pd.DataFrame,
                 shock_rates: pd.DataFrame,
                 rate_col: str,
                 out_col: str) -> pd.DataFrame:
    """
    Aggregate LOO shock rates to state-quarter Bartik instrument.

        B_{s,t} = sum_j omega_{s,j} * g_{-s,j,t}

    Parameters
    ----------
    shares      : Must contain [state_fips, industry_code, share]
    shock_rates : Must contain [state_fips, industry_code, quarter_label, rate_col]
    rate_col    : Name of the shock rate column (e.g. 'g_ld_loo')
    out_col     : Name for the output Bartik column (e.g. 'bartik_ld')
    """
    merged = (
        shock_rates.dropna(subset=[rate_col])
        .merge(
            shares[["state_fips", "industry_code", "share"]],
            on=["state_fips", "industry_code"],
            how="inner",
        )
    )

    def _agg(grp):
        return pd.Series({
            out_col         : (grp["share"] * grp[rate_col]).sum(),
            "weight_sum"    : grp["share"].sum(),
            "n_supersectors": len(grp),
        })

    instrument = (
        merged
        .groupby(["state_fips", "quarter_label"])
        .apply(_agg)
        .reset_index()
    )

    instrument["state"] = instrument["state_fips"].map(FIPS2D_TO_STATE)

    # JOLTS levels are in thousands; QCEW employment in raw units → multiply by 1000
    instrument[out_col] = instrument[out_col] * 1000

    return instrument[["state", "state_fips", "quarter_label",
                        out_col, "weight_sum", "n_supersectors"]]


# ---------------------------------------------------------------------------
# Panel construction  (mirrors part5_lp.py)
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame) -> pd.DataFrame:
    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[["state_fips", "quarter_label", "unemp_rate", "labor_force"]],
        on=["state_fips", "quarter_label"],
        how="inner",
    )
    panel[instr_col] = panel[instr_col] * 100
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)
    panel["unemp_lag1"] = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    expected_prev = quarter_shift(panel["quarter_label"], -1)
    bad_lag = panel["prev_ql"].notna() & (panel["prev_ql"] != expected_prev)
    if bad_lag.any():
        panel.loc[bad_lag, ["unemp_lag1", "lf_log_lag1"]] = np.nan
    panel = panel.drop(columns="prev_ql")
    n_states = panel["state_fips"].nunique()
    n_qtrs = panel["quarter_label"].nunique()
    print(f"  Panel: {len(panel):,} rows  ({n_states} states x {n_qtrs} quarters)")
    return panel


# ---------------------------------------------------------------------------
# LP estimation  (identical to part5_lp.py)
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int,
                   shock_col: str,
                   include_nfci: bool = False) -> dict | None:
    df = base_panel.copy()
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql
    outcomes_lookup = (
        df[["state_fips", "quarter_label", "unemp_rate"]]
        .rename(columns={"quarter_label": "future_ql",
                         "unemp_rate": "y_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")
    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
    df["dep_var"] = df["y_future"] - df["unemp_lag1"]
    required = ["dep_var", shock_col, "unemp_lag1", "lf_log_lag1"]
    if include_nfci:
        required.append("instr_x_nfci")
    df = df.dropna(subset=required).copy()
    if len(df) < 100:
        return None

    state_dummies = pd.get_dummies(df["state_fips"], prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies = pd.get_dummies(df["quarter_label"], prefix="qt",
                                  drop_first=True, dtype=float)
    core_regressors = [shock_col, "unemp_lag1", "lf_log_lag1"]
    if include_nfci:
        core_regressors.append("instr_x_nfci")
    X = pd.concat([df[core_regressors], state_dummies, time_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    y = df["dep_var"]
    model = sm.OLS(y, X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values},
    )
    beta = float(model.params[shock_col])
    se = float(model.bse[shock_col])
    tstat = float(model.tvalues[shock_col])
    pval = float(model.pvalues[shock_col])

    delta_h = delta_h_se = delta_h_t = delta_h_p = None
    if include_nfci:
        delta_h = float(model.params["instr_x_nfci"])
        delta_h_se = float(model.bse["instr_x_nfci"])
        delta_h_t = float(model.tvalues["instr_x_nfci"])
        delta_h_p = float(model.pvalues["instr_x_nfci"])

    return {
        "h": h, "beta": beta, "se": se, "tstat": tstat, "pval": pval,
        "partial_f": tstat ** 2,
        "ci90_lo": beta - Z90 * se, "ci90_hi": beta + Z90 * se,
        "ci95_lo": beta - Z95 * se, "ci95_hi": beta + Z95 * se,
        "delta_h": delta_h, "delta_h_se": delta_h_se,
        "delta_h_t": delta_h_t, "delta_h_p": delta_h_p,
        "nobs": int(model.nobs), "r2": float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range": f"{df['quarter_label'].min()}-{df['quarter_label'].max()}",
    }


def run_lp(base_panel: pd.DataFrame, shock_col: str,
           label: str, include_nfci: bool = False) -> pd.DataFrame:
    nfci_tag = " [+NFCI interaction]" if include_nfci else ""
    print(f"\n{'='*60}")
    print(f"  LP - {label}{nfci_tag}  (h = 0 ... {max(HORIZONS)})")
    print(f"{'='*60}")
    print(f"  {'h':>3}  {'beta_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}  "
          f"{'partial-F':>10}  {'N':>6}")
    print(f"  {'-'*60}")
    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(base_panel, h, shock_col, include_nfci=include_nfci)
        if res is None:
            continue
        rows.append(res)
        sig = ("***" if res["pval"] < 0.01 else
               "**" if res["pval"] < 0.05 else
               "*" if res["pval"] < 0.10 else "")
        print(f"  {h:3d}  {res['beta']:10.4f}  {res['se']:8.4f}  "
              f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
              f"{res['partial_f']:10.2f}  {res['nobs']:6d}  {sig}")

    # Standardize to unit-SD scale
    instr_sd = float(base_panel[shock_col].std())
    scale_cols = ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]
    irf = pd.DataFrame(rows)
    if irf.empty:
        return irf
    irf[scale_cols] = irf[scale_cols] * instr_sd
    if include_nfci:
        irf["delta_h"] = irf["delta_h"] * instr_sd
        irf["delta_h_se"] = irf["delta_h_se"] * instr_sd
    irf["instr_sd"] = instr_sd
    return irf


# ---------------------------------------------------------------------------
# NFCI interaction attachment  (mirrors part5_lp.py)
# ---------------------------------------------------------------------------
def attach_nfci_interaction(panel: pd.DataFrame,
                             nfci: pd.DataFrame,
                             instr_col: str) -> pd.DataFrame:
    panel = panel.merge(
        nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]],
        on="quarter_label", how="left",
    )
    nfci_mean = panel["nfci_risk"].mean()
    panel["nfci_risk_dm"] = panel["nfci_risk"] - nfci_mean
    panel["instr_x_nfci"] = panel[instr_col] * panel["nfci_risk_dm"]
    print(f"  {instr_col}: nfci_risk sample mean = {nfci_mean:.4f}")
    return panel


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_decomposition(results: dict, out_path: Path) -> None:
    """
    Three-panel comparison: Total s, LD, QU IRFs side by side.
    """
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5.5), sharey=True)
    if n == 1:
        axes = [axes]

    colors = {
        "s (total sep.)": "#d62728",
        "LD (layoffs+disch.)": "#1f77b4",
        "QU (quits)": "#2ca02c",
    }

    for ax, (label, df) in zip(axes, results.items()):
        if df.empty:
            ax.set_visible(False)
            continue
        h = df["h"].values
        beta = df["beta"].values
        lo90 = df["ci90_lo"].values
        hi90 = df["ci90_hi"].values
        lo95 = df["ci95_lo"].values
        hi95 = df["ci95_hi"].values
        color = colors.get(label, "#ff7f0e")

        ax.plot(h, lo95, color=color, linewidth=0.7, linestyle="--", alpha=0.6)
        ax.plot(h, hi95, color=color, linewidth=0.7, linestyle="--", alpha=0.6)
        ax.fill_between(h, lo90, hi90, color=color, alpha=0.20,
                        label="90% CI (clustered)")
        ax.plot(h, beta, color=color, linewidth=2.0,
                marker="o", markersize=3.5, label=label)
        ax.axhline(0, color="black", linewidth=0.8)
        for hh in [4, 8, 12, 16]:
            ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)
        ax.set_title(f"IRF - {label}", fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel("pp change in unemp. rate\n(per 1-SD instrument shock)",
                       fontsize=9)
        ax.set_xticks(h)
        ax.legend(fontsize=8, framealpha=0.85)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)

    fig.suptitle(
        "JOLTS Decomposition: Total Separations vs Layoffs+Discharges vs Quits\n"
        "(unemployment IRFs, 1-SD standardized, state+time FEs, outcomes capped 2019Q4)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


def plot_overlay(results: dict, out_path: Path) -> None:
    """
    Overlay all three IRFs on a single panel for direct comparison.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {
        "s (total sep.)": ("#d62728", "-"),
        "LD (layoffs+disch.)": ("#1f77b4", "-"),
        "QU (quits)": ("#2ca02c", "--"),
    }

    for label, df in results.items():
        if df.empty:
            continue
        h = df["h"].values
        beta = df["beta"].values
        lo90 = df["ci90_lo"].values
        hi90 = df["ci90_hi"].values
        color, ls = colors.get(label, ("#ff7f0e", "-"))

        ax.fill_between(h, lo90, hi90, color=color, alpha=0.12)
        ax.plot(h, beta, color=color, linewidth=2.0, linestyle=ls,
                marker="o", markersize=3.5, label=label)

    ax.axhline(0, color="black", linewidth=0.8)
    for hh in [4, 8, 12, 16]:
        ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.set_title(
        "JOLTS Decomposition: Unemployment IRFs Overlaid\n"
        "(1-SD standardized, 90% CI bands, state+time FEs, outcomes capped 2019Q4)",
        fontsize=11,
    )
    ax.set_xlabel("Horizon h (quarters)", fontsize=10)
    ax.set_ylabel("pp change in unemp. rate (per 1-SD instrument shock)", fontsize=10)
    ax.legend(fontsize=10, framealpha=0.85)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5b - JOLTS Decomposition LP: LD vs QU vs Total s")
print("=" * 60)

# -----------------------------------------------------------------------
# [1] Load cached data
# -----------------------------------------------------------------------
print("\n[1] Loading cached data")

for path in [SHARES_PATH, SHOCK_RATES_LD, SHOCK_RATES_QU, LAUS_FILE, NFCI_FILE]:
    if not path.exists():
        sys.exit(f"ERROR: {path} not found. Run prerequisite scripts first.")

shares = pd.read_parquet(SHARES_PATH)
shock_ld = pd.read_parquet(SHOCK_RATES_LD)
shock_qu = pd.read_parquet(SHOCK_RATES_QU)
outcomes = pd.read_parquet(LAUS_FILE)
nfci = pd.read_parquet(NFCI_FILE)
outcomes["state_fips"] = outcomes["state_fips"].str.zfill(2)

print(f"  Shares: {shares.shape}")
print(f"  LD shock rates: {shock_ld.shape}")
print(f"  QU shock rates: {shock_qu.shape}")
print(f"  LAUS: {outcomes.shape}")

# -----------------------------------------------------------------------
# [2] Build LD and QU Bartik instruments
# -----------------------------------------------------------------------
print("\n[2] Building LD and QU Bartik instruments")

ld_instr = build_bartik(shares, shock_ld, "g_ld_loo", "bartik_ld")
qu_instr = build_bartik(shares, shock_qu, "g_qu_loo", "bartik_qu")

INSTR_DIR.mkdir(parents=True, exist_ok=True)
ld_instr.to_csv(LD_INSTR_FILE, index=False)
qu_instr.to_csv(QU_INSTR_FILE, index=False)
print(f"  Saved: {LD_INSTR_FILE}  ({len(ld_instr):,} rows)")
print(f"  Saved: {QU_INSTR_FILE}  ({len(qu_instr):,} rows)")

# Also load baseline s and delta instruments for correlation analysis
s_instr = pd.read_csv(S_INSTR_FILE, dtype={"state_fips": str})
s_instr["state_fips"] = s_instr["state_fips"].str.zfill(2)
delta_instr = pd.read_csv(DELTA_INSTR_FILE, dtype={"state_fips": str})
delta_instr["state_fips"] = delta_instr["state_fips"].str.zfill(2)

# -----------------------------------------------------------------------
# [3] Cross-instrument correlations
# -----------------------------------------------------------------------
print("\n[3] Cross-instrument correlations (2001Q1+)")

# Merge all four instruments on state x quarter
all_instr = (
    s_instr[["state_fips", "quarter_label", "bartik_s"]]
    .merge(ld_instr[["state_fips", "quarter_label", "bartik_ld"]],
           on=["state_fips", "quarter_label"])
    .merge(qu_instr[["state_fips", "quarter_label", "bartik_qu"]],
           on=["state_fips", "quarter_label"])
    .merge(delta_instr[["state_fips", "quarter_label", "bartik_delta"]],
           on=["state_fips", "quarter_label"],
           how="left")
    .query("quarter_label >= '2001Q1'")
)

corr_cols = ["bartik_delta", "bartik_s", "bartik_ld", "bartik_qu"]
corr_matrix = all_instr[corr_cols].corr().round(4)
print("\n  Pooled correlation matrix:")
print(corr_matrix.to_string())

# Save correlation matrix
corr_matrix.to_csv(RESULTS_DIR / "instrument_correlations_decomposed.csv")
print(f"\n  Saved: {RESULTS_DIR / 'instrument_correlations_decomposed.csv'}")

# Key diagnostics
print(f"\n  Key correlations:")
print(f"    corr(B^delta, B^LD)    = {corr_matrix.loc['bartik_delta', 'bartik_ld']:.4f}  "
      f"(expected: moderate positive)")
print(f"    corr(B^delta, B^QU)    = {corr_matrix.loc['bartik_delta', 'bartik_qu']:.4f}  "
      f"(expected: low or negative)")
print(f"    corr(B^LD, B^QU)       = {corr_matrix.loc['bartik_ld', 'bartik_qu']:.4f}")
print(f"    corr(B^s_total, B^LD)  = {corr_matrix.loc['bartik_s', 'bartik_ld']:.4f}")
print(f"    corr(B^s_total, B^QU)  = {corr_matrix.loc['bartik_s', 'bartik_qu']:.4f}")

# -----------------------------------------------------------------------
# [4] Build LP panels
# -----------------------------------------------------------------------
print("\n[4] Building LP panels")
print("  s (total) panel:")
s_panel = build_panel(s_instr, "bartik_s", outcomes)
print("  LD panel:")
ld_panel = build_panel(ld_instr, "bartik_ld", outcomes)
print("  QU panel:")
qu_panel = build_panel(qu_instr, "bartik_qu", outcomes)

# Attach NFCI
s_panel = attach_nfci_interaction(s_panel, nfci, "bartik_s")
ld_panel = attach_nfci_interaction(ld_panel, nfci, "bartik_ld")
qu_panel = attach_nfci_interaction(qu_panel, nfci, "bartik_qu")

# -----------------------------------------------------------------------
# [5] Run LPs
# -----------------------------------------------------------------------
print("\n[5] Running LPs")

irf_s  = run_lp(s_panel,  "bartik_s",  "s (total sep.)")
irf_ld = run_lp(ld_panel, "bartik_ld", "LD (layoffs+disch.)")
irf_qu = run_lp(qu_panel, "bartik_qu", "QU (quits)")

# Save
irf_s.to_csv(RESULTS_DIR / "lp_irf_s_decomp_baseline.csv", index=False)
irf_ld.to_csv(RESULTS_DIR / "lp_irf_ld.csv", index=False)
irf_qu.to_csv(RESULTS_DIR / "lp_irf_qu.csv", index=False)
print(f"\n  Saved: lp_irf_s_decomp_baseline.csv")
print(f"  Saved: lp_irf_ld.csv")
print(f"  Saved: lp_irf_qu.csv")

# -----------------------------------------------------------------------
# [6] Plots
# -----------------------------------------------------------------------
print("\n[6] Plotting")

decomp_results = {
    "s (total sep.)": irf_s,
    "LD (layoffs+disch.)": irf_ld,
    "QU (quits)": irf_qu,
}

plot_decomposition(
    decomp_results,
    RESULTS_DIR / "lp_irf_decomposition_comparison.png",
)

plot_overlay(
    decomp_results,
    RESULTS_DIR / "lp_irf_decomposition_overlay.png",
)

# -----------------------------------------------------------------------
# [7] Summary
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print("Summary - JOLTS Decomposition LP Results")
print("=" * 60)

for lbl, df in [("s (total)", irf_s), ("LD", irf_ld), ("QU", irf_qu)]:
    if df.empty:
        print(f"\n  {lbl}: no results")
        continue
    peak = df.loc[df["beta"].abs().idxmax()]
    sig_horizons = df[df["pval"] < 0.10]["h"].tolist()
    print(f"\n  {lbl} shock:")
    print(f"    Peak |beta|: h={int(peak['h'])}  beta={peak['beta']:.4f}  "
          f"SE={peak['se']:.4f}  p={peak['pval']:.3f}")
    print(f"    Significant at 10%: horizons {sig_horizons}")
    print(f"    Instrument SD: {peak['instr_sd']:.4f}")

print(f"\n  Cross-instrument correlations:")
print(f"    corr(delta, LD) = {corr_matrix.loc['bartik_delta', 'bartik_ld']:.4f}")
print(f"    corr(delta, QU) = {corr_matrix.loc['bartik_delta', 'bartik_qu']:.4f}")
print(f"    corr(LD, QU)    = {corr_matrix.loc['bartik_ld', 'bartik_qu']:.4f}")

print("\n" + "=" * 60)
print("Done. Check data/results/ for outputs.")
print("=" * 60)
