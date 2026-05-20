"""
part6b_shock_persistence_cyclical.py
=====================================
Estimates the stochastic process for z (labor productivity) and δ (permanent
establishment exit) at business-cycle frequencies, closely following the
approach of Coles and Kelishomi (2018, AEJ Macro), Table 2.

Overview
--------
CK calibrate a bivariate AR(1) with correlated innovations:

    ln p_t = ρ_p ln p_{t-1} + ε_t
    ln δ_t = ρ_δ ln δ_{t-1} + (1 - ρ_δ) ln δ̄ + η_t
    (ε_t, η_t) ~ N(0, Σ),   Σ not necessarily diagonal

They estimate the process from log deviations from HP trend (λ = 10^5 at
quarterly frequency) and choose *monthly* autocorrelation parameters so that
the implied process matches the first-order autocorrelation and cross-
correlation at quarterly intervals.

Our approach differs from CK in one key respect: we impose orthogonal
innovations via Cholesky decomposition with z ordered first. The economic
motivation is that the endogenous exit mechanism in our model already captures
the pathway z → x_c → δ_e through equilibrium dynamics; we do not want the
shock innovations themselves to also carry this correlation. The Cholesky
ordering z-first means that η^z is the pure productivity innovation, and η^δ
is the component of exit variation orthogonal to that contemporaneous
productivity shock — i.e., the structural exit shock the model is meant to
capture.

This is distinct from CK who allow correlated innovations (ρ_{pδ} = −0.63)
because their model does not have an endogenous exit mechanism to absorb the
z → δ pathway internally.

Procedure
---------
Step 1: Load raw_data.pkl (built by observables.py). Extract national
        aggregate lp (PRS85006163, output per hour, nonfarm business) and
        delta (BED Deaths employment-weighted exit rate).

Step 2: Apply log transformation then HP filter with λ = 1,600 to both
        series. This is the standard RBC convention and is consistent with
        the SMM moment construction in observables_moments.py (which was
        updated from λ = 100,000 to λ = 1,600 on May 19, 2026 — see
        context/data_and_files.md). The resulting cycle_{z,t} and cycle_{δ,t}
        are log-deviations from trend.

        Note: CK use λ = 10^5, which is close to the Shimer/CK convention for
        monthly data converted to quarterly. We use λ = 1,600 (standard
        quarterly) for consistency with our SMM moments. See
        app:filter_robustness in the draft for robustness across filter choices.

Step 3: Restrict to common sample. BED Deaths begin 1992Q3; lp is available
        from 1947Q1. Common sample is 1992Q3–2019Q4 (pre-pandemic; our LP
        estimation sample ends 2019Q4).

Step 4: Estimate a bivariate VAR(1) on [cycle_z, cycle_δ]:

        [cycle_z_t ]   [a_11  a_12] [cycle_z_{t-1}]   [u^z_t ]
        [cycle_δ_t ] = [a_21  a_22] [cycle_δ_{t-1}] + [u^δ_t ]

        At quarterly frequency. We use OLS equation-by-equation (equivalent
        to GLS since the regressors are identical). HC3 robust SEs reported.

Step 5: Extract the reduced-form innovation covariance matrix Σ_u, then
        Cholesky-decompose: Σ_u = P P', P lower triangular.

        Structural shocks (z ordered first):
            ε^z_t = u^z_t / P[0,0]                    (pure productivity)
            ε^δ_t = (u^δ_t - P[1,0] ε^z_t) / P[1,1]  (exit ⊥ productivity)

        By construction: Cov(ε^z, ε^δ) = 0 at quarterly frequency.

Step 6: Extract quarterly calibration targets:
        - ρ_z^Q  = a_11  (own-persistence of z at quarterly frequency)
        - ρ_δ^Q  = a_22  (own-persistence of δ at quarterly frequency)
        - β(δ←z) = a_21  (endogenous exit response to z; reported, not
                          structural — this pathway is captured by model eqns)
        - σ_z^Q  = std(ε^z_t)
        - σ_δ^Q  = std(ε^δ_t)

Step 7: Convert quarterly AR(1) parameters to monthly, following CK.
        A monthly AR(1) with parameter ρ_m implies a quarterly AR(1) with
        parameter ρ_q = ρ_m^3. Therefore:

            ρ_m = ρ_q^(1/3)

        Innovation SDs convert as:

            For an AR(1) x_t = ρ x_{t-1} + σ ε_t:
            Var(x_t) = σ^2 / (1 - ρ^2)   [unconditional]

            The unconditional variance is the same at monthly and quarterly
            frequency (it's a property of the stationary distribution).
            So σ_m = σ_q * sqrt((1 - ρ_q^2) / (1 - ρ_m^2)).

            Alternatively: the quarterly aggregate of 3 monthly shocks has
            variance 3 σ_m^2 under the approximation that shocks are
            contemporaneous — a rougher but common approach. We use the
            exact unconditional-variance-matching formula above.

Step 8: Save results to data/results/var_calibration_zd.csv and produce
        diagnostic plots.

Cross-check
-----------
Compare ρ_δ^Q obtained here against part6's AR(1) on the residualized Bartik
aggregate (ρ_δ = 0.617). If VAR-orthogonalized ρ_δ^Q is close, both are
measuring the same structural object. If the VAR estimate is lower, it
suggests part6's estimate includes the endogenous z→δ channel.

Output
------
    data/results/var_calibration_zd.csv   — main calibration table
    data/results/var_zd_cycles.png        — HP-filtered cycles time series
    data/results/var_zd_scatter.png       — scatter: z vs δ cycles
    data/results/var_zd_irfs.png          — VAR IRFs: z shock → z, δ

Run
---
    python part6b_shock_persistence_cyclical.py

Prerequisites
-------------
    python observables.py    (produces raw_data.pkl)
"""

import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import linalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── working directory ─────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── parameters ────────────────────────────────────────────────────────────────
HP_LAMBDA   = 1_600          # λ = 1,600: standard quarterly RBC convention
                             # (CK use 10^5; see docstring for discussion)
SAMPLE_END  = "2019Q4"       # Pre-pandemic; consistent with LP sample
RAW_DATA    = Path("../raw_data.pkl")   # built by observables.py


# =============================================================================
#  STEP 1 — Load data
# =============================================================================

def load_raw_data(path: Path) -> pd.DataFrame:
    """
    Load raw_data.pkl produced by observables.py.
    Returns a DataFrame with columns ['lp', 'delta', ...] at quarterly freq.
    """
    dat = pd.read_pickle(path)
    # Ensure quarterly period index
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    return dat


# =============================================================================
#  STEP 2 — Log + HP filter
# =============================================================================

def hp_cycle(series: pd.Series, lamb: int = HP_LAMBDA) -> pd.Series:
    """
    Apply log transformation then HP filter.
    Returns the cyclical component (log-deviation from HP trend).
    Drops leading/trailing NaNs introduced by log(0) or missing data.
    """
    log_s = np.log(series.dropna())
    cycle, _ = sm.tsa.filters.hpfilter(log_s, lamb=lamb)
    return cycle


# =============================================================================
#  STEP 3 — Common sample
# =============================================================================

def common_sample(z_cycle: pd.Series, d_cycle: pd.Series,
                  end: str = SAMPLE_END) -> tuple[pd.Series, pd.Series]:
    """
    Restrict both cycles to their common non-missing sample, up to end.
    BED Deaths begin 1992Q3; lp is available earlier.
    """
    df = pd.DataFrame({"z": z_cycle, "delta": d_cycle}).dropna()
    df = df[df.index <= end]
    return df["z"], df["delta"]


# =============================================================================
#  STEP 4 — Bivariate VAR(1) by OLS
# =============================================================================

def estimate_var1(z: pd.Series, d: pd.Series) -> dict:
    """
    Estimate a bivariate VAR(1) on [z_t, δ_t] by OLS (equation by equation).

    Model:
        z_t     = a_11 z_{t-1} + a_12 δ_{t-1} + u^z_t
        δ_t     = a_21 z_{t-1} + a_22 δ_{t-1} + u^δ_t

    We include a constant in each equation. HC3 robust SEs.

    Returns a dict with:
        A1         : 2×2 companion matrix (top-left block)
        Sigma_u    : 2×2 reduced-form innovation covariance
        resid      : T×2 array of reduced-form residuals [u^z, u^δ]
        params_z   : OLS result for z equation
        params_d   : OLS result for δ equation
        y          : T×2 aligned LHS matrix
        y_lag      : T×2 aligned lagged matrix
    """
    df = pd.DataFrame({"z": z, "d": d}).dropna()
    df_lag = df.shift(1)
    df_est = pd.concat([df, df_lag.add_suffix("_lag")], axis=1).dropna()

    Y   = df_est[["z", "d"]].values          # T × 2
    X   = sm.add_constant(df_est[["z_lag", "d_lag"]].values)  # T × 3

    # z equation
    res_z = sm.OLS(Y[:, 0], X).fit(cov_type="HC3")
    # δ equation
    res_d = sm.OLS(Y[:, 1], X).fit(cov_type="HC3")

    # Companion matrix A1: rows = equations, cols = [z_lag, d_lag]
    A1 = np.array([
        [res_z.params[1], res_z.params[2]],   # z equation lags
        [res_d.params[1], res_d.params[2]],   # δ equation lags
    ])

    # Reduced-form residuals and covariance
    resid = np.column_stack([np.asarray(res_z.resid), np.asarray(res_d.resid)])
    T     = resid.shape[0]
    Sigma_u = resid.T @ resid / (T - 3)   # bias-corrected (3 params per eq)

    return {
        "A1":       A1,
        "Sigma_u":  Sigma_u,
        "resid":    resid,
        "params_z": res_z,
        "params_d": res_d,
        "Y":        Y,
        "X":        X,
    }


# =============================================================================
#  STEP 5 — Cholesky identification (z ordered first)
# =============================================================================

def cholesky_identify(var_result: dict) -> dict:
    """
    Cholesky decomposition of Σ_u with z ordered first.

    Σ_u = P P',   P lower triangular.

    Structural shocks:
        ε^z_t = u^z_t / P[0,0]
        ε^δ_t = (u^δ_t - P[1,0] * ε^z_t) / P[1,1]

    By construction: Cov(ε^z, ε^δ) = 0.

    The impulse-response matrix for a one-SD structural shock is P (each
    column of P is the contemporaneous response to one unit of that structural
    shock).  The structural VAR impact matrix is B = P (lower triangular):
        B[:, 0] = response to z shock      (P[:, 0])
        B[:, 1] = response to δ shock      (P[:, 1])

    Returns structural shock series and their standard deviations.
    """
    Sigma_u = var_result["Sigma_u"]
    resid   = var_result["resid"]          # T × 2: [u^z, u^δ]

    P = linalg.cholesky(Sigma_u, lower=True)   # lower triangular

    # Structural shocks
    eps_z = resid[:, 0] / P[0, 0]
    eps_d = (resid[:, 1] - P[1, 0] * eps_z) / P[1, 1]

    # Sanity check: orthogonality
    corr_check = np.corrcoef(eps_z, eps_d)[0, 1]
    assert abs(corr_check) < 1e-10, (
        f"Structural shocks not orthogonal: corr = {corr_check:.6f}")

    # The structural shocks eps_z, eps_d are standardized (unit SD by
    # construction of the Cholesky normalization).  The innovation SDs in
    # original log-cycle units are the diagonal entries of P:
    #   σ_z^Q = P[0,0]  — SD of the z structural shock in log-cycle units
    #   σ_δ^Q = P[1,1]  — SD of the δ structural shock in log-cycle units
    # (P[1,0] captures how much of u^δ is explained by ε^z contemporaneously)
    return {
        "P":         P,          # 2×2 Cholesky factor (impact matrix)
        "eps_z":     eps_z,      # structural z shock series (standardized)
        "eps_d":     eps_d,      # structural δ shock series (standardized, ⊥ ε^z)
        "sigma_z_Q": float(P[0, 0]),   # log-cycle units
        "sigma_d_Q": float(P[1, 1]),   # log-cycle units
    }


# =============================================================================
#  STEP 6 — Quarterly calibration targets
# =============================================================================

def quarterly_targets(var_result: dict, chol_result: dict) -> dict:
    """
    Extract quarterly-frequency calibration targets from VAR and Cholesky.

    Returns a dict suitable for printing and CSV export.
    """
    A1 = var_result["A1"]
    return {
        # Own-persistence (diagonal of A1)
        "rho_z_Q":    float(A1[0, 0]),
        "rho_d_Q":    float(A1[1, 1]),
        # Cross-lag: δ response to lagged z
        # — this is the endogenous exit channel; reported but NOT fed into
        #   the model's shock process (model equations capture this pathway)
        "beta_dz_Q":  float(A1[1, 0]),
        # Cross-lag: z response to lagged δ (should be near zero)
        "alpha_zd_Q": float(A1[0, 1]),
        # Structural shock SDs at quarterly frequency
        "sigma_z_Q":  chol_result["sigma_z_Q"],
        "sigma_d_Q":  chol_result["sigma_d_Q"],
        # Reduced-form innovation correlation (before Cholesky)
        "corr_uu":    float(
            var_result["Sigma_u"][0, 1]
            / np.sqrt(var_result["Sigma_u"][0, 0] * var_result["Sigma_u"][1, 1])
        ),
    }


# =============================================================================
#  STEP 7 — Monthly parameter conversion (following CK)
# =============================================================================

def quarterly_to_monthly(rho_Q: float, sigma_Q: float) -> tuple[float, float]:
    """
    Convert quarterly AR(1) parameters to monthly, following Coles and
    Kelishomi (2018), Section III.

    A monthly AR(1):   x_t = ρ_m x_{t-1} + σ_m ε_t,   ε_t ~ iid N(0,1)

    implies at quarterly frequency (temporal aggregation of 3 periods):
        ρ_q ≈ ρ_m^3
    so:
        ρ_m = ρ_q^(1/3)

    For the innovation SD: the unconditional variance of the process is
        Var(x) = σ_m^2 / (1 - ρ_m^2)   [monthly]
               = σ_Q^2 / (1 - ρ_Q^2)   [quarterly, same stationary dist.]

    So:   σ_m = σ_Q * sqrt((1 - ρ_Q^2) / (1 - ρ_m^2))

    Note: this holds exactly for the AR(1) itself. The quarterly observable
    is an average of 3 monthly values, which has a more complex covariance
    structure — but CK adopt this simpler unconditional-variance-matching
    approach (p. 124), and we follow them for comparability.

    Parameters
    ----------
    rho_Q   : quarterly AR(1) coefficient
    sigma_Q : quarterly structural shock SD (from Cholesky)

    Returns
    -------
    rho_m, sigma_m : monthly AR(1) coefficient and innovation SD
    """
    if rho_Q <= 0:
        raise ValueError(f"rho_Q = {rho_Q:.4f} ≤ 0; monthly conversion undefined.")
    rho_m   = rho_Q ** (1.0 / 3.0)
    sigma_m = sigma_Q * np.sqrt((1.0 - rho_Q**2) / (1.0 - rho_m**2))
    return float(rho_m), float(sigma_m)


# =============================================================================
#  VAR IMPULSE RESPONSES
# =============================================================================

def var_irfs(A1: np.ndarray, P: np.ndarray, H: int = 20) -> dict:
    """
    Compute VAR IRFs to a one-SD structural shock (Cholesky identified).

    irf_z[:, h] = response of [z, δ] at horizon h to a 1-SD z shock
    irf_d[:, h] = response of [z, δ] at horizon h to a 1-SD δ shock

    Returns dict keyed by shock name, each a (H+1) × 2 array.
    """
    I = np.eye(2)
    # Structural impact matrix: B = P (lower triangular, 1-SD normalization)
    # Each column of P = contemporaneous response to that structural shock
    irfs = {}
    for j, name in enumerate(["z", "delta"]):
        irf = np.zeros((H + 1, 2))
        irf[0, :] = P[:, j]           # impact period: column j of P
        for h in range(1, H + 1):
            irf[h, :] = A1 @ irf[h - 1, :]
        irfs[name] = irf
    return irfs


# =============================================================================
#  PLOTTING
# =============================================================================

REC_SPANS = [
    ("1990-07-01", "1991-03-01"),
    ("2001-03-01", "2001-11-01"),
    ("2007-12-01", "2009-06-01"),
]


def _period_to_ts(idx) -> list:
    """Convert PeriodIndex or string period labels to Timestamps."""
    if hasattr(idx, "to_timestamp"):
        return list(idx.to_timestamp())
    return [pd.Period(q, freq="Q").to_timestamp() for q in idx]


def plot_cycles(z: pd.Series, d: pd.Series, out_path: Path) -> None:
    """Time series of HP-filtered log cycles for z and δ."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    dts = _period_to_ts(z.index)

    for ax, series, label, color in [
        (axes[0], z, "Labor productivity (z)", "#1f77b4"),
        (axes[1], d, "Exit rate (δ)",          "#d62728"),
    ]:
        ax.plot(dts, series.values * 100, color=color, lw=1.2)
        ax.axhline(0, color="black", lw=0.7)
        for rs, re in REC_SPANS:
            ax.axvspan(pd.Timestamp(rs), pd.Timestamp(re),
                       alpha=0.12, color="grey")
        ax.set_ylabel("Log-cycle (×100)", fontsize=9)
        ax.set_title(f"HP-filtered log cycle — {label}  (λ={HP_LAMBDA:,})",
                     fontsize=10)
        ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle("Business-cycle components: z and δ (common sample)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


def plot_scatter(z: pd.Series, d: pd.Series, out_path: Path) -> None:
    """Scatter of cycle_z vs cycle_δ to visualise raw correlation."""
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(z.values * 100, d.values * 100, alpha=0.4, s=14,
               color="#555555")
    # OLS line
    X   = sm.add_constant(z.values)
    fit = sm.OLS(d.values, X).fit()
    xl  = np.linspace(z.min(), z.max(), 100)
    ax.plot(xl * 100, (fit.params[0] + fit.params[1] * xl) * 100,
            color="red", lw=1.5, label=f"OLS slope={fit.params[1]:.2f}")
    corr = float(z.corr(d))
    ax.set_xlabel("cycle_z (×100)", fontsize=9)
    ax.set_ylabel("cycle_δ (×100)", fontsize=9)
    ax.set_title(f"cycle_z vs cycle_δ  (corr={corr:+.3f})", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(lw=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


def plot_irfs(irfs: dict, out_path: Path) -> None:
    """
    VAR IRFs: responses of z and δ to each structural shock.
    2 rows (z shock, δ shock) × 2 columns (response of z, response of δ).
    """
    H    = irfs["z"].shape[0] - 1
    h_ax = np.arange(H + 1)
    var_labels = ["z (labor productivity)", "δ (exit rate)"]
    shk_labels = ["z shock (ε^z)", "δ shock (ε^δ, ⊥ z)"]
    colors = ["#1f77b4", "#d62728"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for row, (shock, skl) in enumerate(zip(["z", "delta"], shk_labels)):
        irf = irfs[shock]          # (H+1) × 2
        for col, (vl, color) in enumerate(zip(var_labels, colors)):
            ax = axes[row, col]
            ax.plot(h_ax, irf[:, col] * 100, color=color, lw=2.0,
                    marker="o", ms=3)
            ax.axhline(0, color="black", lw=0.8)
            ax.set_title(f"{skl}  →  {vl}", fontsize=9)
            ax.set_xlabel("Horizon h (quarters)", fontsize=8)
            ax.set_ylabel("Response (×100)", fontsize=8)
            ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle(
        "VAR(1) IRFs — Cholesky structural shocks (z ordered first)\n"
        "ε^δ is exit shock orthogonal to contemporaneous productivity",
        fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


# =============================================================================
#  PRINTING
# =============================================================================

def print_var_results(var_result: dict) -> None:
    res_z = var_result["params_z"]
    res_d = var_result["params_d"]
    A1    = var_result["A1"]
    Sigma = var_result["Sigma_u"]

    print("\n  VAR(1) coefficient matrix A1:")
    print(f"         z_lag      δ_lag")
    print(f"  z  :  {A1[0,0]:+.4f}    {A1[0,1]:+.4f}   "
          f"  (SE: {res_z.bse[1]:.4f}, {res_z.bse[2]:.4f})")
    print(f"  δ  :  {A1[1,0]:+.4f}    {A1[1,1]:+.4f}   "
          f"  (SE: {res_d.bse[1]:.4f}, {res_d.bse[2]:.4f})")
    print(f"\n  Reduced-form innovation covariance Σ_u (×10^4):")
    print(f"    [{Sigma[0,0]*1e4:.4f}  {Sigma[0,1]*1e4:.4f}]")
    print(f"    [{Sigma[1,0]*1e4:.4f}  {Sigma[1,1]*1e4:.4f}]")
    corr_uu = Sigma[0,1] / np.sqrt(Sigma[0,0] * Sigma[1,1])
    print(f"    corr(u^z, u^δ) = {corr_uu:+.4f}")
    print(f"  N obs: {int(res_z.nobs)},  R²_z = {res_z.rsquared:.4f},  "
          f"R²_δ = {res_d.rsquared:.4f}")
    if res_d.rsquared < 0.15:
        print(f"  WARNING: R²_δ = {res_d.rsquared:.4f} is low — BED Deaths is a "
              f"noisy quarterly series (N=110). Own-persistence ρ_δ may be "
              f"imprecisely estimated; interpret with caution.")


def print_calibration_summary(qt: dict, rho_z_m: float, sigma_z_m: float,
                               rho_d_m: float, sigma_d_m: float) -> None:
    sep = "=" * 65
    print(f"\n{sep}")
    print("CALIBRATION SUMMARY")
    print(sep)
    print(f"\n  {'Parameter':<35}  {'Quarterly':>10}  {'Monthly':>10}")
    print(f"  {'-'*35}  {'-'*10}  {'-'*10}")
    print(f"  {'ρ_z  (own-persistence of z)':35}  {qt['rho_z_Q']:>10.4f}  "
          f"{rho_z_m:>10.4f}")
    print(f"  {'ρ_δ  (own-persistence of δ)':35}  {qt['rho_d_Q']:>10.4f}  "
          f"{rho_d_m:>10.4f}")
    print(f"  {'σ_z  (structural z shock SD)':35}  {qt['sigma_z_Q']:>10.4f}  "
          f"{sigma_z_m:>10.4f}")
    print(f"  {'σ_δ  (structural δ shock SD)':35}  {qt['sigma_d_Q']:>10.4f}  "
          f"{sigma_d_m:>10.4f}")
    print(f"\n  {'Reported (not in shock process)':35}")
    print(f"  {'β(δ←z)  (endogenous exit resp.)':35}  {qt['beta_dz_Q']:>10.4f}")
    print(f"  {'α(z←δ)  (near-zero expected)':35}  {qt['alpha_zd_Q']:>10.4f}")
    print(f"  {'corr(u^z, u^δ)  (pre-Cholesky)':35}  {qt['corr_uu']:>+10.4f}")
    print(f"\n  Cholesky ordering: z first → ε^δ ⊥ ε^z by construction")
    print(f"  corr(ε^z, ε^δ) = 0.0000  (imposed)")
    print(f"\n  Cross-check (part6 Bartik-residualized ρ_δ = 0.617):")
    print(f"    VAR-orthogonalized ρ_δ^Q = {qt['rho_d_Q']:.4f}")
    gap = qt['rho_d_Q'] - 0.617
    print(f"    Difference = {gap:+.4f}")
    if abs(gap) >= 0.08:
        print("    DIAGNOSTIC: Large gap likely reflects two effects:")
        print("      (a) VAR conditions on z_{t-1}, absorbing the portion of")
        print("          delta persistence driven by the z->delta_e endogenous channel.")
        print("      (b) BED Deaths is noisy at quarterly frequency -- part6's")
        print("          Bartik aggregation across industries reduces noise.")
        print("    The VAR estimate is the preferred model calibration target")
        print("    because it strips the endogenous channel. The gap documents")
        print("    how much of raw delta persistence is endogenous vs. structural.")


# =============================================================================
#  MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 65)
    print("Part 6b — Shock Persistence: Cyclical VAR (z and δ)")
    print(f"HP filter λ = {HP_LAMBDA:,}  |  Sample end: {SAMPLE_END}")
    print("=" * 65)

    # ── 1. Load ───────────────────────────────────────────────────────────────
    print("\n[1] Loading raw_data.pkl")
    dat = load_raw_data(RAW_DATA)
    print(f"  Loaded {len(dat)} obs, columns: {list(dat.columns)}")
    print(f"  lp   : {dat['lp'].first_valid_index()} – "
          f"{dat['lp'].last_valid_index()}")
    print(f"  delta: {dat['delta'].first_valid_index()} – "
          f"{dat['delta'].last_valid_index()}")

    # ── 2. Log + HP filter ────────────────────────────────────────────────────
    print(f"\n[2] Log + HP filter (λ = {HP_LAMBDA:,})")
    z_cycle = hp_cycle(dat["lp"])
    d_cycle = hp_cycle(dat["delta"])
    print(f"  cycle_z: {len(z_cycle)} obs, SD = {z_cycle.std()*100:.3f}×10⁻²")
    print(f"  cycle_δ: {len(d_cycle)} obs, SD = {d_cycle.std()*100:.3f}×10⁻²")

    # ── 3. Common sample ──────────────────────────────────────────────────────
    print(f"\n[3] Common sample (BED Deaths start ∩ {SAMPLE_END})")
    z, d = common_sample(z_cycle, d_cycle)
    print(f"  Sample: {z.index[0]} – {z.index[-1]},  N = {len(z)}")

    # ── 4. VAR(1) estimation ──────────────────────────────────────────────────
    print("\n[4] Bivariate VAR(1) estimation (OLS, HC3 SEs)")
    var_result = estimate_var1(z, d)
    print_var_results(var_result)

    # ── 5. Cholesky identification ────────────────────────────────────────────
    print("\n[5] Cholesky identification (z ordered first)")
    chol = cholesky_identify(var_result)
    print(f"  Cholesky impact matrix P:")
    P = chol["P"]
    print(f"    [{P[0,0]:.6f}   0.000000]")
    print(f"    [{P[1,0]:.6f}  {P[1,1]:.6f}]")
    print(f"  σ_z^Q = {chol['sigma_z_Q']:.6f}  (SD of structural z shocks)")
    print(f"  σ_δ^Q = {chol['sigma_d_Q']:.6f}  (SD of structural δ shocks ⊥ z)")
    print(f"  corr(ε^z, ε^δ) ≈ 0  [imposed by construction]")

    # ── 6. Quarterly targets ──────────────────────────────────────────────────
    print("\n[6] Quarterly calibration targets")
    qt = quarterly_targets(var_result, chol)

    # ── 7. Monthly conversion ─────────────────────────────────────────────────
    print("\n[7] Monthly parameter conversion (CK approach: ρ_m = ρ_q^(1/3))")
    rho_z_m, sigma_z_m = quarterly_to_monthly(qt["rho_z_Q"], qt["sigma_z_Q"])
    rho_d_m, sigma_d_m = quarterly_to_monthly(qt["rho_d_Q"], qt["sigma_d_Q"])
    print(f"  z: ρ_m = {rho_z_m:.4f},  σ_m = {sigma_z_m:.6f}")
    print(f"  δ: ρ_m = {rho_d_m:.4f},  σ_m = {sigma_d_m:.6f}")

    # ── 8. Summary ────────────────────────────────────────────────────────────
    print_calibration_summary(qt, rho_z_m, sigma_z_m, rho_d_m, sigma_d_m)

    # ── 9. Save CSV ───────────────────────────────────────────────────────────
    print("\n[8] Saving results")
    A1 = var_result["A1"]
    rows = [
        # Quarterly
        {"param": "rho_z",       "freq": "quarterly", "value": qt["rho_z_Q"],
         "note": "own-persistence of z (diagonal A1[0,0])"},
        {"param": "rho_delta",   "freq": "quarterly", "value": qt["rho_d_Q"],
         "note": "own-persistence of delta (diagonal A1[1,1])"},
        {"param": "sigma_z",     "freq": "quarterly", "value": qt["sigma_z_Q"],
         "note": "structural z shock SD in log-cycle units (P[0,0])"},
        {"param": "sigma_delta", "freq": "quarterly", "value": qt["sigma_d_Q"],
         "note": "structural delta shock SD orthogonal to z (P[1,1])"},
        {"param": "beta_dz",     "freq": "quarterly", "value": qt["beta_dz_Q"],
         "note": "endogenous exit response to z (A1[1,0]); not in shock process"},
        {"param": "alpha_zd",    "freq": "quarterly", "value": qt["alpha_zd_Q"],
         "note": "z response to delta (A1[0,1]); expected near zero"},
        {"param": "corr_uu",     "freq": "quarterly", "value": qt["corr_uu"],
         "note": "reduced-form innovation corr (pre-Cholesky)"},
        # Monthly
        {"param": "rho_z",       "freq": "monthly",   "value": rho_z_m,
         "note": "rho_q^(1/3)"},
        {"param": "rho_delta",   "freq": "monthly",   "value": rho_d_m,
         "note": "rho_q^(1/3)"},
        {"param": "sigma_z",     "freq": "monthly",   "value": sigma_z_m,
         "note": "unconditional-variance-matched from quarterly"},
        {"param": "sigma_delta", "freq": "monthly",   "value": sigma_d_m,
         "note": "unconditional-variance-matched from quarterly"},
        # Metadata
        {"param": "hp_lambda",   "freq": "n/a",       "value": HP_LAMBDA,
         "note": "HP filter smoothing parameter"},
        {"param": "sample_start","freq": "n/a",       "value": str(z.index[0]),
         "note": "common sample start"},
        {"param": "sample_end",  "freq": "n/a",       "value": str(z.index[-1]),
         "note": "common sample end"},
        {"param": "nobs",        "freq": "n/a",       "value": len(z),
         "note": "number of quarterly observations"},
    ]
    out_csv = RESULTS_DIR / "var_calibration_zd.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"  Saved: {out_csv}")

    # -- 10. Plots ----------------------------------------------------------------
    print("\n[9] Plotting")
    plot_cycles(z, d, RESULTS_DIR / "var_zd_cycles.png")
    plot_scatter(z, d, RESULTS_DIR / "var_zd_scatter.png")
    irfs = var_irfs(var_result["A1"], chol["P"], H=20)
    plot_irfs(irfs, RESULTS_DIR / "var_zd_irfs.png")

    print("\nDone.")
