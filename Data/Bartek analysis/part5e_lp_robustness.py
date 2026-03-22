"""
part5e_lp_robustness.py — Priority 4 Robustness Checks for LP Results
========================================================================
PRIORITY 4 from the project brief: three robustness checks for the baseline
panel LP results on δ and s shocks.

Implements:
  4a. Recession-severity placebo: B_{s,t} × NFCI_t vs B_{s,t} × Δu^nat_t
      Tests whether NFCI interaction can be explained by recession depth.

  4b. Wild cluster bootstrap p-values
      Cameron-Gelbach-Miller (2008) with Rademacher weights at state level
      to supplement lower-bound state-clustered SEs (n=50 clusters).

  4c. Pre-GFC sample restriction for s
      Restrict s instrument sample to 2001Q1-2007Q3 to test whether the
      monotonically rising s IRF is driven entirely by GFC observations.

Outputs
-------
  data/results/lp_irf_delta_nfci_vs_recsev.csv          — δ: NFCI vs recession severity
  data/results/lp_irf_s_nfci_vs_recsev.csv              — s: NFCI vs recession severity
  data/results/lp_irf_delta_wcb.csv                     — δ: wild cluster bootstrap
  data/results/lp_irf_s_wcb.csv                         — s: wild cluster bootstrap
  data/results/lp_irf_s_preGFC.csv                      — s: pre-GFC restriction
  data/results/lp_robustness_comparison.png             — comparison plot

Run
---
    python part5e_lp_robustness.py

Prerequisites
--------------
    python part5_lp.py  (for baseline results)
    python part4_outcomes.py  (produces laus_quarterly.parquet)
    python part7_nfci.py  (produces nfci_quarterly.parquet)
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
# Configuration (matches part5_lp.py)
# ---------------------------------------------------------------------------
T          = 20
HORIZONS   = list(range(T + 1))
BASE_YEAR  = 2006
CI_LEVEL   = 0.90
Z90        = 1.645
Z95        = 1.960

MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR  = Path("data/instruments")
RESULTS_DIR = Path("data/results")
CACHE_DIR   = Path("data/cache")

DELTA_INSTR_FILE = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
S_INSTR_FILE     = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
LAUS_FILE        = INSTR_DIR / "laus_quarterly.parquet"
NFCI_FILE        = CACHE_DIR / "nfci_quarterly.parquet"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Wild cluster bootstrap config
WCB_REPS = 999  # Standard for publication-quality inference
WCB_HORIZONS = [0, 1, 2, 4, 8, 12, 16]  # Key horizons only to save time


# ---------------------------------------------------------------------------
# Quarter arithmetic (identical to part5_lp.py)
# ---------------------------------------------------------------------------
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    """Shift quarter_label strings by h quarters."""
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_instruments() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load δ and s Bartik instruments."""
    delta = pd.read_csv(DELTA_INSTR_FILE, dtype={"state_fips": str})
    s     = pd.read_csv(S_INSTR_FILE,     dtype={"state_fips": str})
    delta["state_fips"] = delta["state_fips"].str.zfill(2)
    s["state_fips"]     = s["state_fips"].str.zfill(2)
    print(f"  δ instrument: {delta.shape}")
    print(f"  s instrument: {s.shape}")
    return delta, s


def load_outcomes() -> pd.DataFrame:
    """Load quarterly LAUS unemployment rates and labor force."""
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    print(f"  LAUS: {laus.shape}")
    return laus[["state_fips", "state", "quarter_label", "unemp_rate", "labor_force"]]


def load_nfci() -> pd.DataFrame:
    """Load quarterly NFCI data."""
    nfci = pd.read_parquet(NFCI_FILE)
    print(f"  NFCI: {nfci.shape}")
    return nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]]


# ---------------------------------------------------------------------------
# Recession severity: Compute Δu^nat_t (national unemployment change)
# ---------------------------------------------------------------------------
def compute_national_unemployment_change(outcomes: pd.DataFrame) -> pd.DataFrame:
    """
    Compute national unemployment rate change from t-1 to t.

    Uses labor-force-weighted mean of state unemployment rates:
        u^nat_t = Σ_s (LF_s,t / Σ_s' LF_s',t) × u_s,t

    Then computes:
        Δu^nat_t = u^nat_t - u^nat_{t-1}

    Returns DataFrame with columns: quarter_label, unemp_nat, unemp_nat_change
    """
    # Compute labor-force weighted national unemployment each quarter
    outcomes_sorted = outcomes.sort_values(["quarter_label", "state_fips"]).reset_index(drop=True)

    nat = (
        outcomes_sorted
        .groupby("quarter_label")
        .agg({
            "unemp_rate": lambda u: (u * outcomes_sorted.loc[u.index, "labor_force"]).sum() / outcomes_sorted.loc[u.index, "labor_force"].sum(),
            "labor_force": "sum"
        })
        .rename(columns={"unemp_rate": "unemp_nat", "labor_force": "total_lf"})
        .reset_index()
    )

    # Compute change as lagged difference
    nat = nat.sort_values("quarter_label").reset_index(drop=True)
    nat["unemp_nat_change"] = nat["unemp_nat"].diff()

    print(f"\n  National unemployment:")
    print(f"    Mean u^nat = {nat['unemp_nat'].mean():.3f}")
    print(f"    Mean Δu^nat = {nat['unemp_nat_change'].mean():.4f}  "
          f"(SD = {nat['unemp_nat_change'].std():.4f})")

    return nat[["quarter_label", "unemp_nat", "unemp_nat_change"]]


# ---------------------------------------------------------------------------
# Panel construction (mirrors part5_lp.py)
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame,
                sample_end: str | None = None) -> pd.DataFrame:
    """Build LP panel with pre-determined controls."""

    if sample_end is not None:
        instr = instr[instr["quarter_label"] <= sample_end].copy()

    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[["state_fips", "quarter_label", "unemp_rate", "labor_force"]],
        on=["state_fips", "quarter_label"],
        how="inner",
    )

    panel[instr_col] = panel[instr_col] * 100
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)

    # Construct lags within state
    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )

    # Verify lag continuity
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    expected_prev    = quarter_shift(panel["quarter_label"], -1)
    bad_lag          = panel["prev_ql"].notna() & (panel["prev_ql"] != expected_prev)
    if bad_lag.any():
        panel.loc[bad_lag, ["unemp_lag1", "lf_log_lag1"]] = np.nan

    panel = panel.drop(columns="prev_ql")

    n_states = panel["state_fips"].nunique()
    n_qtrs   = panel["quarter_label"].nunique()
    print(f"  Panel: {len(panel):,} rows  ({n_states} states × {n_qtrs} quarters)")

    return panel


# ---------------------------------------------------------------------------
# Attach NFCI and recession severity interactions
# ---------------------------------------------------------------------------
def attach_interactions(panel: pd.DataFrame, nfci: pd.DataFrame,
                       nat_unemp: pd.DataFrame, instr_col: str) -> pd.DataFrame:
    """
    Merge NFCI and national unemployment change.
    Construct demeaned interaction terms:
      - instr_x_nfci:    B_{s,t} × (NFCI_risk_dm)
      - instr_x_recsev:  B_{s,t} × (Δu^nat_t_dm)
    """
    # Merge NFCI
    panel = panel.merge(
        nfci[["quarter_label", "nfci_risk"]],
        on="quarter_label",
        how="left",
    )

    # Merge national unemployment change
    panel = panel.merge(
        nat_unemp[["quarter_label", "unemp_nat_change"]],
        on="quarter_label",
        how="left",
    )

    # Demean within sample
    nfci_mean = panel["nfci_risk"].mean()
    recsev_mean = panel["unemp_nat_change"].mean()

    panel["nfci_risk_dm"] = panel["nfci_risk"] - nfci_mean
    panel["recsev_dm"] = panel["unemp_nat_change"] - recsev_mean

    # Interactions (instrument already in pp)
    panel["instr_x_nfci"] = panel[instr_col] * panel["nfci_risk_dm"]
    panel["instr_x_recsev"] = panel[instr_col] * panel["recsev_dm"]

    print(f"  {instr_col}:")
    print(f"    NFCI_risk sample mean = {nfci_mean:.4f}")
    print(f"    Δu^nat sample mean = {recsev_mean:.4f}  (SD = {panel['recsev_dm'].std():.4f})")

    return panel


# ---------------------------------------------------------------------------
# LP estimation — single horizon
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int, shock_col: str,
                   interaction_cols: list[str] | None = None) -> dict | None:
    """
    Estimate LP for single horizon h with optional interactions.

    interaction_cols: list of column names to include as additional regressors
                      (e.g., ["instr_x_nfci"] or ["instr_x_nfci", "instr_x_recsev"])
    """
    df = base_panel.copy()

    # Construct h-quarter-ahead outcome
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql

    outcomes_lookup = (
        df[["state_fips", "quarter_label", "unemp_rate"]]
        .rename(columns={"quarter_label": "future_ql", "unemp_rate": "y_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")

    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()

    df["dep_var"] = df["y_future"] - df["unemp_lag1"]

    # Check required variables
    required = ["dep_var", shock_col, "unemp_lag1", "lf_log_lag1"]
    if interaction_cols:
        required.extend(interaction_cols)

    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        return None

    # Build regressor matrix
    state_dummies = pd.get_dummies(df["state_fips"], prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)

    core_regressors = [shock_col, "unemp_lag1", "lf_log_lag1"]
    if interaction_cols:
        core_regressors.extend(interaction_cols)

    X = pd.concat([df[core_regressors], state_dummies, time_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    y = df["dep_var"]

    # OLS with state-clustered SEs
    model = sm.OLS(y, X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values},
    )

    beta  = float(model.params[shock_col])
    se    = float(model.bse[shock_col])
    tstat = float(model.tvalues[shock_col])
    pval  = float(model.pvalues[shock_col])

    result = {
        "h": h,
        "beta": beta,
        "se": se,
        "tstat": tstat,
        "pval": pval,
        "partial_f": tstat ** 2,
        "ci90_lo": beta - Z90 * se,
        "ci90_hi": beta + Z90 * se,
        "ci95_lo": beta - Z95 * se,
        "ci95_hi": beta + Z95 * se,
        "nobs": int(model.nobs),
        "r2": float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range": f"{df['quarter_label'].min()}–{df['quarter_label'].max()}",
    }

    # Store interaction coefficients if present
    if interaction_cols:
        for col in interaction_cols:
            result[f"{col}_beta"] = float(model.params[col])
            result[f"{col}_se"] = float(model.bse[col])
            result[f"{col}_t"] = float(model.tvalues[col])
            result[f"{col}_p"] = float(model.pvalues[col])

    return result


# ---------------------------------------------------------------------------
# LP estimation — all horizons
# ---------------------------------------------------------------------------
def run_lp(base_panel: pd.DataFrame, shock_col: str, label: str,
           interaction_cols: list[str] | None = None) -> pd.DataFrame:
    """Run LP for all horizons with optional interactions."""

    interaction_tag = ""
    if interaction_cols:
        interaction_tag = f" [{', '.join(interaction_cols)}]"

    print(f"\n{'='*70}")
    print(f"  LP — {label}{interaction_tag}  (h = 0 … {max(HORIZONS)})")
    print(f"{'='*70}")

    if interaction_cols:
        header = f"  {'h':>3}  {'β_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}"
        for col in interaction_cols:
            header += f"  {col}β:>8"
        header += f"  {'N':>6}"
        print(header)
    else:
        print(f"  {'h':>3}  {'β_h':>10}  {'SE':>8}  {'t':>7}  {'p':>6}  {'N':>6}")

    print(f"  {'-'*70}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(base_panel, h, shock_col, interaction_cols=interaction_cols)
        if res is None:
            continue
        rows.append(res)

        sig = ("***" if res["pval"] < 0.01 else
               "**"  if res["pval"] < 0.05 else
               "*"   if res["pval"] < 0.10 else "")

        line = f"  {h:3d}  {res['beta']:10.4f}  {res['se']:8.4f}  {res['tstat']:7.3f}  {res['pval']:6.3f}"
        if interaction_cols:
            for col in interaction_cols:
                line += f"  {res[f'{col}_beta']:8.4f}"
        line += f"  {res['nobs']:6d}  {sig}"
        print(line)

    # Standardize to unit-SD scale
    instr_sd = float(base_panel[shock_col].std())
    irf = pd.DataFrame(rows)
    if irf.empty:
        return irf

    scale_cols = ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]
    irf[scale_cols] = irf[scale_cols] * instr_sd

    if interaction_cols:
        for col in interaction_cols:
            irf[f"{col}_beta"] = irf[f"{col}_beta"] * instr_sd
            irf[f"{col}_se"] = irf[f"{col}_se"] * instr_sd

    irf["instr_sd"] = instr_sd
    return irf


# ---------------------------------------------------------------------------
# Wild cluster bootstrap (Cameron-Gelbach-Miller 2008)
# ---------------------------------------------------------------------------
def wild_cluster_bootstrap(base_panel: pd.DataFrame, h: int, shock_col: str,
                          n_reps: int = 999,
                          seed: int | None = None) -> tuple[float, float]:
    """
    Wild cluster bootstrap p-value for H0: beta_{shock_col} = 0 at horizon h.

    Implements the Cameron-Gelbach-Miller (2008) WCB-RE procedure:
      1. Estimate RESTRICTED model (excluding shock_col) → restricted residuals
      2. Estimate UNRESTRICTED model → original t-stat (with state-clustered SEs)
      3. For each bootstrap rep:
         a. Rademacher weights w_g ∈ {-1, +1} at cluster (state) level
         b. y*_b = X_R * beta_R + w_g * e_R  (null-imposed DGP)
         c. Estimate unrestricted model on y*_b → t*_b
      4. p-value = fraction of |t*_b| >= |t_orig|

    Using RESTRICTED residuals is essential: unrestricted residuals bake
    the shock_col signal into fitted values, making the bootstrap t-stats
    centered near t_orig rather than near 0 — which produces meaningless
    p-values near 0.5.

    Returns (bootstrap_p, original_t_stat_clustered).
    """
    if seed is not None:
        np.random.seed(seed)

    df = base_panel.copy()

    # Construct outcome
    future_ql = quarter_shift(df["quarter_label"], h)
    df["future_ql"] = future_ql
    outcomes_lookup = (
        df[["state_fips", "quarter_label", "unemp_rate"]]
        .rename(columns={"quarter_label": "future_ql", "unemp_rate": "y_future"})
    )
    df = df.merge(outcomes_lookup, on=["state_fips", "future_ql"], how="left")
    if MAX_OUTCOME_QUARTER is not None:
        df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
    df["dep_var"] = df["y_future"] - df["unemp_lag1"]
    df = df.dropna(subset=["dep_var", shock_col, "unemp_lag1", "lf_log_lag1"]).copy()
    if len(df) < 100:
        return np.nan, np.nan

    # Build regressor matrices
    state_dummies = pd.get_dummies(df["state_fips"], prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)
    controls = df[["unemp_lag1", "lf_log_lag1"]]

    # UNRESTRICTED model (includes shock_col)
    X_u = pd.concat([df[[shock_col]], controls, state_dummies, time_dummies], axis=1)
    X_u = sm.add_constant(X_u, has_constant="add")
    y = df["dep_var"].values

    model_u = sm.OLS(y, X_u).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values},
    )
    t_orig = float(model_u.tvalues[shock_col])

    # RESTRICTED model (excludes shock_col — imposes H0: beta = 0)
    X_r = pd.concat([controls, state_dummies, time_dummies], axis=1)
    X_r = sm.add_constant(X_r, has_constant="add")
    model_r = sm.OLS(y, X_r).fit()
    fitted_r = model_r.fittedvalues.values
    resid_r  = model_r.resid.values

    # Bootstrap
    states = df["state_fips"].values
    unique_states = np.unique(states)
    n_clusters = len(unique_states)
    # Pre-compute cluster membership for speed
    cluster_idx = {st: np.where(states == st)[0] for st in unique_states}

    t_boot = np.empty(n_reps)
    t_boot[:] = np.nan
    X_u_arr = X_u.values
    shock_col_idx = list(X_u.columns).index(shock_col)

    for b in range(n_reps):
        # Rademacher weights at cluster level
        weights = np.random.choice([-1.0, 1.0], size=n_clusters)

        # Construct y*_b under the null
        resid_boot = resid_r.copy()
        for i, st in enumerate(unique_states):
            resid_boot[cluster_idx[st]] *= weights[i]
        y_boot = fitted_r + resid_boot

        # Re-estimate unrestricted model (no clustering — just need t-stat)
        try:
            model_b = sm.OLS(y_boot, X_u_arr).fit()
            t_boot[b] = model_b.tvalues[shock_col_idx]
        except Exception:
            pass

    valid = t_boot[~np.isnan(t_boot)]
    if len(valid) == 0:
        return np.nan, t_orig

    # Two-tailed p-value
    p_boot = float(np.mean(np.abs(valid) >= np.abs(t_orig)))
    return p_boot, t_orig


# ---------------------------------------------------------------------------
# ROBUSTNESS CHECK 4a: NFCI vs Recession Severity (Placebo)
# ---------------------------------------------------------------------------
def robustness_4a(delta_panel: pd.DataFrame, s_panel: pd.DataFrame,
                  delta_instr: pd.DataFrame, s_instr: pd.DataFrame,
                  nfci: pd.DataFrame, nat_unemp: pd.DataFrame,
                  outcomes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    4a. Recession severity placebo for interaction term.

    Run three parallel interaction versions:
      1. B_{s,t} × NFCI_t (baseline)
      2. B_{s,t} × Δu^nat_t (pure recession severity)
      3. Both simultaneously (horse race)
    """
    print("\n" + "="*70)
    print("ROBUSTNESS CHECK 4a: NFCI vs Recession Severity Placebo")
    print("="*70)

    # Baseline: B × NFCI
    print("\n[4a.1] δ shock: baseline B × NFCI")
    irf_delta_nfci = run_lp(delta_panel, "bartik_delta", "δ shock",
                            interaction_cols=["instr_x_nfci"])

    print("\n[4a.2] δ shock: recession severity B × Δu^nat")
    irf_delta_recsev = run_lp(delta_panel, "bartik_delta", "δ shock",
                              interaction_cols=["instr_x_recsev"])

    print("\n[4a.3] δ shock: horse race (both interactions)")
    irf_delta_horsrace = run_lp(delta_panel, "bartik_delta", "δ shock",
                                interaction_cols=["instr_x_nfci", "instr_x_recsev"])

    # s shock
    print("\n[4a.4] s shock: baseline B × NFCI")
    irf_s_nfci = run_lp(s_panel, "bartik_s", "s shock",
                        interaction_cols=["instr_x_nfci"])

    print("\n[4a.5] s shock: recession severity B × Δu^nat")
    irf_s_recsev = run_lp(s_panel, "bartik_s", "s shock",
                          interaction_cols=["instr_x_recsev"])

    print("\n[4a.6] s shock: horse race (both interactions)")
    irf_s_horsrace = run_lp(s_panel, "bartik_s", "s shock",
                            interaction_cols=["instr_x_nfci", "instr_x_recsev"])

    # Combine for output: NFCI vs recession severity comparison
    irf_delta_comp = irf_delta_nfci.copy()
    irf_delta_comp["recsev_beta"] = irf_delta_recsev["instr_x_recsev_beta"].values
    irf_delta_comp["recsev_se"] = irf_delta_recsev["instr_x_recsev_se"].values
    irf_delta_comp["recsev_p"] = irf_delta_recsev["instr_x_recsev_p"].values
    irf_delta_comp["horsrace_nfci_beta"] = irf_delta_horsrace["instr_x_nfci_beta"].values
    irf_delta_comp["horsrace_nfci_se"] = irf_delta_horsrace["instr_x_nfci_se"].values
    irf_delta_comp["horsrace_recsev_beta"] = irf_delta_horsrace["instr_x_recsev_beta"].values
    irf_delta_comp["horsrace_recsev_se"] = irf_delta_horsrace["instr_x_recsev_se"].values

    irf_s_comp = irf_s_nfci.copy()
    irf_s_comp["recsev_beta"] = irf_s_recsev["instr_x_recsev_beta"].values
    irf_s_comp["recsev_se"] = irf_s_recsev["instr_x_recsev_se"].values
    irf_s_comp["recsev_p"] = irf_s_recsev["instr_x_recsev_p"].values
    irf_s_comp["horsrace_nfci_beta"] = irf_s_horsrace["instr_x_nfci_beta"].values
    irf_s_comp["horsrace_nfci_se"] = irf_s_horsrace["instr_x_nfci_se"].values
    irf_s_comp["horsrace_recsev_beta"] = irf_s_horsrace["instr_x_recsev_beta"].values
    irf_s_comp["horsrace_recsev_se"] = irf_s_horsrace["instr_x_recsev_se"].values

    return irf_delta_comp, irf_s_comp


# ---------------------------------------------------------------------------
# ROBUSTNESS CHECK 4b: Wild Cluster Bootstrap
# ---------------------------------------------------------------------------
def robustness_4b(delta_panel: pd.DataFrame, s_panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    4b. Wild cluster bootstrap p-values (Cameron-Gelbach-Miller 2008).

    Complement state-clustered SEs with Rademacher weights at cluster level.
    B = 49 bootstrap replications at key horizons only.
    """
    print("\n" + "="*70)
    print("ROBUSTNESS CHECK 4b: Wild Cluster Bootstrap P-values")
    print("="*70)
    print(f"\n  Running {WCB_REPS} bootstrap replications at horizons {WCB_HORIZONS}...")

    wcb_delta = []
    wcb_s = []

    for h in WCB_HORIZONS:
        # δ shock
        p_delta, t_delta = wild_cluster_bootstrap(
            delta_panel, h, "bartik_delta", n_reps=WCB_REPS, seed=42+h
        )
        wcb_delta.append({
            "h": h,
            "wcb_pval": p_delta,
            "orig_tstat": t_delta,
        })

        # s shock
        p_s, t_s = wild_cluster_bootstrap(
            s_panel, h, "bartik_s", n_reps=WCB_REPS, seed=100+h
        )
        wcb_s.append({
            "h": h,
            "wcb_pval": p_s,
            "orig_tstat": t_s,
        })

        print(f"    Completed h={h}")

    irf_delta_wcb = pd.DataFrame(wcb_delta)
    irf_s_wcb = pd.DataFrame(wcb_s)

    return irf_delta_wcb, irf_s_wcb


# ---------------------------------------------------------------------------
# ROBUSTNESS CHECK 4c: Pre-GFC Sample for s
# ---------------------------------------------------------------------------
def robustness_4c(s_instr: pd.DataFrame, outcomes: pd.DataFrame,
                  nfci: pd.DataFrame) -> pd.DataFrame:
    """
    4c. Pre-GFC sample restriction for s.

    Restrict the s instrument sample to 2001Q1-2007Q3 (pre-GFC) to test
    whether the monotonically rising s IRF is driven entirely by GFC.
    """
    print("\n" + "="*70)
    print("ROBUSTNESS CHECK 4c: Pre-GFC Sample Restriction for s")
    print("="*70)

    print("\n  Building s panel with pre-GFC restriction (2001Q1-2007Q3)")
    s_pregfc = build_panel(s_instr, "bartik_s", outcomes, sample_end="2007Q3")

    # Attach NFCI for consistency
    s_pregfc = s_pregfc.merge(
        nfci[["quarter_label", "nfci_risk"]],
        on="quarter_label",
        how="left",
    )

    nfci_mean = s_pregfc["nfci_risk"].mean()
    s_pregfc["nfci_risk_dm"] = s_pregfc["nfci_risk"] - nfci_mean
    s_pregfc["instr_x_nfci"] = s_pregfc["bartik_s"] * s_pregfc["nfci_risk_dm"]

    print(f"  NFCI_risk sample mean (pre-GFC) = {nfci_mean:.4f}")

    irf_s_pregfc = run_lp(s_pregfc, "bartik_s", "s shock (pre-GFC restriction)")

    return irf_s_pregfc


# ---------------------------------------------------------------------------
# Plotting: comparison of robustness checks
# ---------------------------------------------------------------------------
def plot_robustness_comparison(irf_delta: dict[str, pd.DataFrame],
                              irf_s: dict[str, pd.DataFrame],
                              out_path: Path) -> None:
    """
    Create comparison plots for robustness checks.

    irf_delta/irf_s: dict mapping label → DataFrame
    """
    print(f"\nPlotting robustness comparison...")

    n_rows = 2
    n_cols = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 10))
    axes = axes.flatten()

    # Color scheme
    colors = {
        "δ baseline": "#1f77b4",
        "δ NFCI": "#1f77b4",
        "δ Recession severity": "#ff7f0e",
        "δ Horse race (NFCI)": "#1f77b4",
        "δ Horse race (RecSev)": "#ff7f0e",
        "s baseline": "#d62728",
        "s NFCI": "#d62728",
        "s Recession severity": "#ff7f0e",
        "s Horse race (NFCI)": "#d62728",
        "s Horse race (RecSev)": "#ff7f0e",
        "s pre-GFC": "#2ca02c",
    }

    plot_idx = 0

    # [0] δ: NFCI vs RecSev interaction (4a)
    # The irf_delta dict contains a single DataFrame with all results
    if not irf_delta["nfci"].empty:
        ax = axes[plot_idx]
        df = irf_delta["nfci"]
        h = df["h"].values

        ax.plot(h, df["instr_x_nfci_beta"].values,
               color=colors["δ NFCI"], linewidth=2, marker="o", markersize=3,
               label="B × NFCI")
        ax.plot(h, df["recsev_beta"].values,
               color=colors["δ Recession severity"], linewidth=2, marker="s", markersize=3,
               label="B × ΔU^nat")

        ax.axhline(0, color="black", linewidth=0.8)
        for hh in [4, 8, 12, 16]:
            ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)

        ax.set_title("δ shock: NFCI vs Recession Severity Interaction", fontsize=10)
        ax.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax.set_ylabel("Interaction coeff.", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        plot_idx += 1

    # [1] s: NFCI vs RecSev interaction (4a)
    if not irf_s["nfci"].empty:
        ax = axes[plot_idx]
        df = irf_s["nfci"]
        h = df["h"].values

        ax.plot(h, df["instr_x_nfci_beta"].values,
               color=colors["s NFCI"], linewidth=2, marker="o", markersize=3,
               label="B × NFCI")
        ax.plot(h, df["recsev_beta"].values,
               color=colors["s Recession severity"], linewidth=2, marker="s", markersize=3,
               label="B × ΔU^nat")

        ax.axhline(0, color="black", linewidth=0.8)
        for hh in [4, 8, 12, 16]:
            ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)

        ax.set_title("s shock: NFCI vs Recession Severity Interaction", fontsize=10)
        ax.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax.set_ylabel("Interaction coeff.", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        plot_idx += 1

    # [2] s: Full sample vs Pre-GFC (4c)
    if not irf_s["full"].empty and not irf_s["pregfc"].empty:
        ax = axes[plot_idx]
        df_full = irf_s["full"]
        df_pregfc = irf_s["pregfc"]

        h_full = df_full["h"].values
        h_pregfc = df_pregfc["h"].values

        ax.plot(h_full, df_full["beta"].values,
               color=colors["s baseline"], linewidth=2, marker="o", markersize=3,
               label="Full sample (2001-2019Q4)")
        ax.plot(h_pregfc, df_pregfc["beta"].values,
               color=colors["s pre-GFC"], linewidth=2, marker="s", markersize=3,
               label="Pre-GFC (2001-2007Q3)")

        ax.axhline(0, color="black", linewidth=0.8)
        for hh in [4, 8, 12, 16]:
            ax.axvline(hh, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)

        ax.set_title("s shock: Full vs Pre-GFC Sample", fontsize=10)
        ax.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax.set_ylabel("IRF coefficient (1-SD shock)", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        plot_idx += 1

    # [3] Summary note
    if plot_idx < len(axes):
        ax = axes[plot_idx]
        ax.axis("off")

        summary_text = (
            "PRIORITY 4 ROBUSTNESS CHECKS\n\n"
            "4a. Recession severity placebo\n"
            "    Tests if NFCI interaction is absorbed by ΔU^nat\n"
            "    If δ_h survives horse race, there is independent financial content\n\n"
            "4b. Wild cluster bootstrap p-values\n"
            "    Complements state-clustered SEs (n=50, lower bound)\n"
            "    Uses Rademacher weights, 999 reps per horizon\n\n"
            "4c. Pre-GFC sample restriction for s\n"
            "    Tests if monotonic s IRF driven by GFC observations\n"
            "    Restricts to 2001Q1-2007Q3 instrumental sample"
        )

        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
               fontsize=8, verticalalignment="top", family="monospace",
               bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle(
        "Priority 4 Robustness Checks: LP Results\n"
        "(NFCI vs RecSev Placebo | WCB | Pre-GFC Restriction)",
        fontsize=12, y=0.98
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 70)
print("Part 5e — Priority 4 LP Robustness Checks")
print("=" * 70)

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------
print("\n[0] Loading data")
delta_instr, s_instr = load_instruments()
outcomes = load_outcomes()
nfci = load_nfci()

# Compute national unemployment change
print("\n[0b] Computing national unemployment change (recession severity)")
nat_unemp = compute_national_unemployment_change(outcomes)

# -----------------------------------------------------------------------
# Build baseline panels
# -----------------------------------------------------------------------
print("\n[1] Building baseline panels")
print("  δ panel:")
delta_panel = build_panel(delta_instr, "bartik_delta", outcomes)

print("  s panel:")
s_panel = build_panel(s_instr, "bartik_s", outcomes)

# Attach NFCI and recession severity interactions to both
print("\n[1b] Attaching NFCI and recession severity interactions")
delta_panel = attach_interactions(delta_panel, nfci, nat_unemp, "bartik_delta")
s_panel = attach_interactions(s_panel, nfci, nat_unemp, "bartik_s")

# -----------------------------------------------------------------------
# Robustness Check 4a: Recession Severity Placebo
# -----------------------------------------------------------------------
irf_delta_4a, irf_s_4a = robustness_4a(
    delta_panel, s_panel,
    delta_instr, s_instr,
    nfci, nat_unemp, outcomes
)

# Save 4a results
print("\n[2] Saving 4a results")
irf_delta_4a.to_csv(RESULTS_DIR / "lp_irf_delta_nfci_vs_recsev.csv", index=False)
irf_s_4a.to_csv(RESULTS_DIR / "lp_irf_s_nfci_vs_recsev.csv", index=False)
print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_nfci_vs_recsev.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_nfci_vs_recsev.csv'}")

# -----------------------------------------------------------------------
# Robustness Check 4b: Wild Cluster Bootstrap
# -----------------------------------------------------------------------
irf_delta_4b, irf_s_4b = robustness_4b(delta_panel, s_panel)

# Save 4b results
print("\n[3] Saving 4b results")
irf_delta_4b.to_csv(RESULTS_DIR / "lp_irf_delta_wcb.csv", index=False)
irf_s_4b.to_csv(RESULTS_DIR / "lp_irf_s_wcb.csv", index=False)
print(f"  Saved: {RESULTS_DIR / 'lp_irf_delta_wcb.csv'}")
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_wcb.csv'}")

# -----------------------------------------------------------------------
# Robustness Check 4c: Pre-GFC Sample for s
# -----------------------------------------------------------------------
irf_s_pregfc = robustness_4c(s_instr, outcomes, nfci)

# Save 4c results
print("\n[4] Saving 4c results")
irf_s_pregfc.to_csv(RESULTS_DIR / "lp_irf_s_preGFC.csv", index=False)
print(f"  Saved: {RESULTS_DIR / 'lp_irf_s_preGFC.csv'}")

# -----------------------------------------------------------------------
# Load baseline results for comparison
# -----------------------------------------------------------------------
print("\n[5] Loading baseline results for comparison")
try:
    irf_delta_baseline = pd.read_csv(RESULTS_DIR / "lp_irf_delta_post2001.csv")
    irf_s_baseline = pd.read_csv(RESULTS_DIR / "lp_irf_s.csv")
    print("  Baseline δ and s IRFs loaded")
except FileNotFoundError:
    print("  [warn] Baseline results not found — plotting will skip baseline comparison")
    irf_delta_baseline = None
    irf_s_baseline = None

# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------
print("\n[6] Creating robustness comparison plots")

# Prepare plot data — just the 4a results since they contain all info
irf_delta_dict = {
    "nfci": irf_delta_4a,
}

irf_s_dict = {
    "nfci": irf_s_4a,
    "full": irf_s_baseline if irf_s_baseline is not None else irf_s_4a,
    "pregfc": irf_s_pregfc,
}

plot_robustness_comparison(
    irf_delta_dict,
    irf_s_dict,
    RESULTS_DIR / "lp_robustness_comparison.png"
)

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
print("\n" + "="*70)
print("SUMMARY — Priority 4 Robustness Checks")
print("="*70)

print("\n[4a] Recession Severity Placebo:")
print("  δ shock: NFCI interaction coefficient")
if not irf_delta_4a.empty:
    peak_nfci = irf_delta_4a.loc[irf_delta_4a["instr_x_nfci_beta"].abs().idxmax()]
    print(f"    Peak at h={int(peak_nfci['h'])}: {peak_nfci['instr_x_nfci_beta']:.4f}  "
          f"(p={peak_nfci['instr_x_nfci_p']:.3f})")

print("  s shock: NFCI interaction coefficient")
if not irf_s_4a.empty:
    peak_nfci = irf_s_4a.loc[irf_s_4a["instr_x_nfci_beta"].abs().idxmax()]
    print(f"    Peak at h={int(peak_nfci['h'])}: {peak_nfci['instr_x_nfci_beta']:.4f}  "
          f"(p={peak_nfci['instr_x_nfci_p']:.3f})")

print("\n[4b] Wild Cluster Bootstrap:")
print(f"  Completed {WCB_REPS} replications per horizon")
print(f"  δ shock: {len(irf_delta_4b)} horizons with WCB p-values")
print(f"  s shock: {len(irf_s_4b)} horizons with WCB p-values")

print("\n[4c] Pre-GFC Sample Restriction for s:")
print(f"  Sample: 2001Q1 - 2007Q3 (pre-GFC quarters)")
if not irf_s_pregfc.empty:
    peak = irf_s_pregfc.loc[irf_s_pregfc["beta"].abs().idxmax()]
    print(f"  Peak at h={int(peak['h'])}: {peak['beta']:.4f}  SE={peak['se']:.4f}  "
          f"p={peak['pval']:.3f}")

print("\n" + "="*70)
print("Done. Check data/results/ for outputs.")
print("="*70)
