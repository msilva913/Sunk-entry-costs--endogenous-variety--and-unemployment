"""
part6b_shock_persistence_cyclical.py
=====================================
Business-cycle processes for the model's forcing variables, at the monthly
frequency the model is solved at.

Produces one table covering four series:

    z      labor productivity        (PRS85006163, index, normalized to 1)
    delta  permanent exit rate       (BED Deaths, employment-weighted)
    tau    total separation rate     (Shimer 2012 continuous-time)
    s      match separation rate     (constructed, see below)

with, for each: the sample mean, the cyclical standard deviation, and an
AR(1) persistence and conditional standard deviation converted to monthly.
Only z, delta and s are model shocks; tau is reported because it is what the
data measure directly and what the moment table targets.

Why s is not tau
----------------
Shimer's series is the TOTAL separation rate: every worker who leaves
employment, whether the job vanished with the firm or the match simply broke
on a surviving product line. The model's s is the second kind only. The two
are linked by the accounting identity behind tab:calib_targets,

    tau = delta_e + (1 - delta_e) s    =>    s = (tau - delta_e)/(1 - delta_e)

Feeding raw tau to the s process would hand the s shock variation that belongs
to delta. In practice the correction turns out to be close to a pure level
shift -- delta_e is under a tenth of tau and the two cycles correlate at 0.99 --
but the identity is what the calibration uses, so this is the right input.

Identification of the z and delta processes
-------------------------------------------
For z and delta we estimate a bivariate VAR(1) on the HP cycles and orthogonalize
by Cholesky with z ordered first, rather than fitting two univariate AR(1)s.
The reason is that the model already generates the z -> chi^c -> delta_e pathway
endogenously through the job creation condition and the continuation-cost
threshold. Letting the exogenous innovations also carry that correlation would
double-count the mechanism under study. This is where we depart from Coles and
Kelishomi (2018), who allow correlated innovations (rho_{p,delta} = -0.63)
because their model has no endogenous exit margin to absorb the pathway.

The two pieces do separate jobs, and it is worth keeping them apart. The VAR
supplies persistence: rho_z and rho_delta are the diagonal of the estimated
companion matrix, formed by OLS before the innovation covariance is even
computed, so the Cholesky has nothing to do with them. The Cholesky supplies
one number, sigma_delta, by stripping the part of the delta innovation that
the contemporaneous z innovation explains. Since z is ordered first, sigma_z
is simply its reduced-form residual SD. The orthogonalization is therefore
doing calibration work here, not only identifying impulse responses: the model
assumes Cov(eps^z, eps^delta) = 0, and the reduced-form innovations correlate
at -0.29, so the raw delta residual SD would overstate the exogenous shock.

The off-diagonal VAR terms are reported but excluded from the shock processes
for the same reason: they are equilibrium outcomes, not primitive forcing.
The same argument applies to s, which is therefore a univariate AR(1) rather
than a third VAR equation.

Filtering and frequency
-----------------------
Log, then HP filter with lambda = 1,600 -- the standard quarterly convention,
and the one observables_moments.py uses, so the shock processes and the SMM
moments see the same cyclical object. CK use lambda = 10^5; see
app:filter_robustness for the comparison. Each series is filtered over its own
full sample and then restricted to the common window, so the trend is not
estimated off a truncated series.

Quarterly estimates convert to monthly following CK: rho_m = rho_Q^(1/3), and
sigma_m matches the unconditional variance, sigma_m = sigma_Q *
sqrt((1 - rho_Q^2)/(1 - rho_m^2)).

Output
------
    data/results/shock_calibration.csv          the table
    data/results/cycles_panel.png               HP cycles, all four series
    data/results/separation_decomposition.png   tau, delta, s: levels and cycles
    data/results/var_zd_scatter.png             cycle_z vs cycle_delta
    data/results/var_zd_irfs.png                VAR IRFs to structural shocks

Run
---
    python part6b_shock_persistence_cyclical.py

Prerequisites
-------------
    python observables.py    (produces raw_data.pkl)
"""

import os
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant
from statsmodels.tsa.filters.hp_filter import hpfilter
import matplotlib
from pathlib import Path


# ── working directory ─────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    # __file__ is undefined when code is exec'd without a file, e.g. running a
    # selection rather than the whole file in Spyder. Fail loudly instead of
    # chdir'ing somewhere that may not exist.
    if not Path("part6b_shock_persistence_cyclical.py").exists():
        raise RuntimeError(
            "Cannot locate the script directory: __file__ is undefined and the "
            "working directory is not 'Data/Bartek analysis'."
        )


# ── matplotlib backend ────────────────────────────────────────────────────────
def _interactive_session() -> bool:
    """True inside an IPython kernel (Spyder, Jupyter), False as a script."""
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False


# Agg writes files but has no display, so plt.show() on it is a silent no-op.
# matplotlib.use() is sticky for the life of a kernel, so a Spyder console that
# once ran an Agg-setting script stays on Agg through every later run; switching
# back here makes the program self-healing. Must precede the pyplot import.
INTERACTIVE = _interactive_session()
if INTERACTIVE:
    if matplotlib.get_backend().lower() == "agg":
        try:
            matplotlib.use("module://matplotlib_inline.backend_inline")
        except Exception as exc:                       # pragma: no cover
            print(f"[backend] could not select the inline backend: {exc}")
else:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ── parameters ────────────────────────────────────────────────────────────────
HP_LAMBDA   = 1_600
SAMPLE_END  = "2019Q4"                   # pre-pandemic; matches the LP sample
RAW_DATA    = Path("../raw_data.pkl")    # built by observables.py
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Display labels and colors, shared by the table and every plot.
SERIES = {
    "z":     ("Labor productivity (z)",        "#1f77b4"),
    "delta": ("Permanent exit rate (delta)",   "#d62728"),
    "tau":   ("Total separation rate (tau)",   "#7f7f7f"),
    "s":     ("Match separation rate (s)",     "#2ca02c"),
}

# Plausibility bands for the units guard. Both rate series must be MONTHLY
# fractions: tau sits near 0.031 and delta near 0.0027. Wide enough to admit a
# sensible revision, narrow enough to catch a percent/fraction or
# monthly/quarterly mix-up.
TAU_BOUNDS   = (0.010, 0.060)
DELTA_BOUNDS = (0.0005, 0.010)

# Recessions inside the sample. The 1990-91 recession is omitted: BED Deaths
# begin in 1992Q3, so shading it would mark an empty region.
REC_SPANS = [("2001-03-01", "2001-11-01"),
             ("2007-12-01", "2009-06-01")]


# =============================================================================
#  SERIES
# =============================================================================

def normalized_z(dat):
    """
    Labor productivity as an index with unit mean.

    Defined once and used everywhere, because z must be filtered over its whole
    1951- history rather than the table window: truncating it to 1992 moves
    rho_z and would silently change the calibration in the draft.
    """
    return dat["lp"] / dat["lp"].mean()


def build_levels(dat):
    """
    Assemble the three rate series: delta, tau and the constructed s.

    From observables.py both inputs are already monthly fractions observed
    quarterly: tau = 1 - exp(-s_h) is a monthly separation probability, and
    delta is the quarterly BED deaths employment share divided by 3. The
    identity below is a within-period accounting statement at the model's
    frequency, so it needs both on that footing -- hence the guard.

    z is not included. It lives on a different sample and carries no
    informative level; see normalized_z.
    """
    tau, delta = dat["s"], dat["delta"]     # note: observables' 's' IS tau

    df = pd.DataFrame({"tau": tau, "delta": delta}).dropna()
    for name, series, bounds in [("tau", df["tau"], TAU_BOUNDS),
                                 ("delta", df["delta"], DELTA_BOUNDS)]:
        if not bounds[0] <= series.mean() <= bounds[1]:
            raise ValueError(
                f"{name} averages {series.mean():.5f}, outside the "
                f"monthly-fraction band {bounds}. Check that observables.py "
                f"still delivers a monthly rate as a fraction."
            )

    df["s"] = (df["tau"] - df["delta"]) / (1.0 - df["delta"])
    if (df["s"] <= 0).any():
        first = df.index[df["s"] <= 0][0]
        raise ValueError(
            f"Match separation rate is non-positive, first at {first}: delta "
            f"exceeds tau there, so the two series are measured on "
            f"inconsistent bases."
        )

    return df[["delta", "tau", "s"]]


def hp_cycles(dat, levels):
    """
    Log-HP cycles, restricted to the common window.

    tau and s are both filtered over the s window rather than tau's own longer
    history. s cannot begin before BED Deaths do, and the decomposition
    compares the two cycles, which requires one trend window for both.
    """
    sources = {"z": normalized_z(dat)} | {c: levels[c] for c in levels}
    out = {}
    for name, series in sources.items():
        out[name], _ = hpfilter(np.log(series.dropna()), lamb=HP_LAMBDA)
    cyc = pd.DataFrame(out).dropna()
    return cyc[cyc.index <= SAMPLE_END]


# =============================================================================
#  PROCESSES
# =============================================================================

def fit_ar1(cycle: pd.Series):
    """
    AR(1) with a constant on a quarterly HP cycle, HC3 robust SEs.

    The constant is kept even though the cycle is mean-zero by construction,
    so the estimate is not forced through the origin on a subsample.
    """
    df = pd.DataFrame({"y": cycle, "x": cycle.shift(1)}).dropna()
    fit = OLS(df["y"], add_constant(df["x"])).fit(cov_type="HC3")
    return {"rho_Q":  float(fit.params["x"]),
            "rho_se": float(fit.bse["x"]),
            "sigma_Q": float(np.std(fit.resid, ddof=2)),
            "r2":     float(fit.rsquared),
            "nobs":   int(fit.nobs)}


def fit_var_cholesky(cyc):
    """
    Bivariate VAR(1) on [cycle_z, cycle_delta], Cholesky with z ordered first.

    Returns the companion matrix, the reduced-form innovation covariance, the
    impact matrix P, and the own-persistence and orthogonalized innovation SD
    for each of z and delta. The diagonal of P is the structural shock SD in
    log-cycle units; the off-diagonal absorbs the contemporaneous z->delta
    pathway that the model generates endogenously.
    """
    df = pd.concat([cyc[["z", "delta"]],
                    cyc[["z", "delta"]].shift(1).add_suffix("_lag")],
                   axis=1).dropna()
    X = add_constant(df[["z_lag", "delta_lag"]].values)
    fits = [OLS(df[c].values, X).fit(cov_type="HC3") for c in ("z", "delta")]

    A1 = np.array([f.params[1:] for f in fits])          # rows = equations
    resid = np.column_stack([f.resid for f in fits])
    Sigma_u = resid.T @ resid / (resid.shape[0] - 3)     # 3 params per equation
    P = np.linalg.cholesky(Sigma_u)                      # lower triangular

    # The own-lag coefficient sits at position 1 in the z equation and 2 in the
    # delta equation, since X = [const, z_lag, delta_lag].
    return {"A1": A1, "P": P,
            "corr_uu": float(Sigma_u[0, 1] /
                             np.sqrt(Sigma_u[0, 0] * Sigma_u[1, 1])),
            # Persistence is the diagonal of A1 -- an OLS object, formed before
            # Sigma_u exists. The Cholesky plays no part in it.
            "rho":    {"z": float(A1[0, 0]), "delta": float(A1[1, 1])},
            "rho_se": {"z": float(fits[0].bse[1]),
                       "delta": float(fits[1].bse[2])},
            "r2":     {"z": float(fits[0].rsquared),
                       "delta": float(fits[1].rsquared)},
            # Innovation SDs are the diagonal of P. For z, ordered first, that
            # is just the reduced-form residual SD. For delta it is the part
            # orthogonal to the contemporaneous z innovation, which is the only
            # number in the table the Cholesky actually changes.
            "sigma": {"z": float(P[0, 0]), "delta": float(P[1, 1])},
            "sigma_rf": {"z": float(np.sqrt(Sigma_u[0, 0])),
                         "delta": float(np.sqrt(Sigma_u[1, 1]))},
            "nobs": int(fits[0].nobs)}


def to_monthly(rho_Q, sigma_Q):
    """
    Quarterly AR(1) to monthly, following CK. A monthly AR(1) aggregates to
    rho_Q = rho_m^3, and the unconditional variance sigma^2/(1-rho^2) is a
    property of the stationary distribution, so it is the same at both
    frequencies. Exact for the AR(1) itself; the quarterly observable is an
    average of three monthly values, which CK also abstract from.
    """
    if rho_Q <= 0:
        raise ValueError(f"rho_Q = {rho_Q:.4f} <= 0; monthly conversion undefined.")
    rho_m = rho_Q ** (1 / 3)
    return rho_m, sigma_Q * np.sqrt((1 - rho_Q**2) / (1 - rho_m**2))


def build_table(levels, cyc, var):
    """
    One row per series: mean, cyclical SD, and the AR(1) process at both
    frequencies.

    The `source` column refers to the estimator, and the two parameters do not
    share one:

      rho    z, delta  diagonal of the VAR companion matrix (OLS)
             tau, s    univariate AR(1)
      sigma  z         reduced-form residual SD -- z is ordered first, so the
                       Cholesky leaves it untouched
             delta     residual SD orthogonalized against the contemporaneous
                       z innovation. The only Cholesky-dependent number here
             tau, s    AR(1) residual SD

    sigma_Q_rf keeps the un-orthogonalized SD alongside so the size of that
    single adjustment is visible, and rho_Q_uni does the same for the VAR's
    effect on persistence.

    rho_se and r2 always come from whichever fit supplied rho_Q, so they
    describe the estimate printed beside them rather than the cross-check.
    """
    rows = []
    for name in SERIES:
        uni = fit_ar1(cyc[name])
        if name in var["rho"]:                       # z and delta
            est = {k: var[k][name] for k in
                   ("rho", "rho_se", "r2", "sigma", "sigma_rf")}
            est["nobs"] = var["nobs"]
        else:                                        # tau and s
            est = {"rho":   uni["rho_Q"],  "rho_se":   uni["rho_se"],
                   "r2":    uni["r2"],     "sigma":    uni["sigma_Q"],
                   "nobs":  uni["nobs"],   "sigma_rf": uni["sigma_Q"]}

        rho_m, sigma_m = to_monthly(est["rho"], est["sigma"])
        rows.append({
            "series":     name,
            "mean":       1.0 if name == "z" else float(levels[name].mean()),
            "units":      "index (normalized)" if name == "z" else "fraction/month",
            "sd_cycle":   float(cyc[name].std()),
            "rho_Q":      est["rho"],
            "sigma_Q":    est["sigma"],
            "rho_m":      rho_m,
            "sigma_m":    sigma_m,
            "source":     {"z": "VAR(1)", "delta": "VAR(1), sigma perp z"}
                          .get(name, "AR(1)"),
            "rho_se":     est["rho_se"],
            "r2":         est["r2"],
            "nobs":       est["nobs"],
            "rho_Q_uni":  uni["rho_Q"],
            "sigma_Q_rf": est["sigma_rf"],
        })
    return pd.DataFrame(rows).set_index("series")


def report(table, var, cyc):
    """Print the headline table, then the estimation detail behind it."""
    fmt = {"mean": "{:.4f}".format, "sd_cycle": "{:.4f}".format,
           "rho_m": "{:.4f}".format, "sigma_m": "{:.6f}".format}
    head = table[["mean", "units", "sd_cycle", "rho_m", "sigma_m", "source"]]

    print(f"\n{'='*78}\nSHOCK CALIBRATION — monthly, HP({HP_LAMBDA:,}) log cycles, "
          f"{cyc.index[0]}–{cyc.index[-1]}\n{'='*78}")
    print(head.to_string(formatters=fmt))
    print("\n  tau is an outcome of delta and s in the model, not a shock; it is "
          "shown\n  because it is the series the moment table targets directly.")

    print(f"\n{'-'*78}\nQuarterly estimates and cross-checks\n{'-'*78}")
    detail = table[["rho_Q", "rho_se", "rho_Q_uni", "sigma_Q", "sigma_Q_rf",
                    "r2", "nobs"]]
    print(detail.to_string(float_format=lambda x: f"{x:.4f}"))
    print("\n  rho_Q_uni: univariate AR(1). For z and delta the table takes the "
          "VAR value\n  instead; the gap is the lagged z->delta channel the VAR "
          "conditions on.")
    d = table.loc["delta"]
    print(f"\n  sigma_Q_rf: reduced-form residual SD, before orthogonalization."
          f"\n  Only delta differs, because only it is ordered second: "
          f"{d['sigma_Q_rf']:.4f} -> {d['sigma_Q']:.4f},"
          f"\n  a {100*(1 - d['sigma_Q']/d['sigma_Q_rf']):.1f}% reduction. That "
          f"is the entire footprint of the Cholesky here."
          f"\n  Persistence is an OLS object and is untouched by it.")

    A1 = var["A1"]
    print("\n  VAR(1) off-diagonals (reported, not fed to the model):")
    print(f"    beta(delta<-z)     = {A1[1,0]:+.4f}  endogenous exit response to z")
    print(f"    alpha(z<-delta)    = {A1[0,1]:+.4f}  expected near zero")
    print(f"    corr(u^z, u^delta) = {var['corr_uu']:+.4f}  pre-Cholesky")

    if d["r2"] < 0.15:
        print(f"\n  WARNING: R^2 for the delta equation is {d['r2']:.4f}. BED "
              f"Deaths is noisy at quarterly\n  frequency (N={int(d['nobs'])}); "
              f"read rho_delta with caution.")


# =============================================================================
#  PLOTS
# =============================================================================

def _finish(fig, out_path):
    """Save, then display (Spyder) or release (script). Always saves, so the
    file outputs do not depend on how the program was launched."""
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Plot saved: {out_path}")
    plt.show() if INTERACTIVE else plt.close(fig)


def _ts(idx):
    """PeriodIndex or period labels to Timestamps."""
    return list(idx.to_timestamp()) if hasattr(idx, "to_timestamp") else \
        [pd.Period(q, freq="Q").to_timestamp() for q in idx]


def _decorate(ax):
    ax.axhline(0, color="black", lw=0.7)
    for lo, hi in REC_SPANS:
        ax.axvspan(pd.Timestamp(lo), pd.Timestamp(hi), alpha=0.12, color="grey")
    ax.grid(axis="y", lw=0.4, alpha=0.4)


def plot_cycle_panel(cyc, out_path):
    """One stacked panel per series, common sample and common scale."""
    fig, axes = plt.subplots(len(SERIES), 1, figsize=(11, 9), sharex=True)
    dts = _ts(cyc.index)
    for ax, (name, (label, color)) in zip(axes, SERIES.items()):
        ax.plot(dts, cyc[name].values * 100, color=color, lw=1.2)
        ax.set_ylabel("Log-cycle (x100)", fontsize=8)
        ax.set_title(label, fontsize=9)
        _decorate(ax)
    fig.suptitle(f"HP-filtered log cycles (lambda = {HP_LAMBDA:,})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


def plot_separation_decomposition(levels, cyc, out_path):
    """
    Overlay of tau, delta and s. The top panel shows the size of the split in
    levels, the bottom shows the cycles the AR(1) is fitted to.
    """
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    names = ["tau", "delta", "s"]
    styles = {"tau": "-", "delta": "-", "s": "--"}

    for ax, frame, scale, ylab, title in [
        (axes[0], levels, 100, "Percent per month", "Levels"),
        (axes[1], cyc,    100, "Log-cycle (x100)",
         f"HP-filtered log cycles (lambda = {HP_LAMBDA:,})"),
    ]:
        dts = _ts(frame.index)
        for name in names:
            ax.plot(dts, frame[name].values * scale, color=SERIES[name][1],
                    ls=styles[name], lw=1.3, label=SERIES[name][0])
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_title(title, fontsize=10)
        _decorate(ax)
        ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

    fig.suptitle("Separation rate decomposition: tau = delta + (1 - delta) s",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


def plot_scatter(cyc, out_path):
    """cycle_z against cycle_delta, with the OLS fit."""
    z, d = cyc["z"], cyc["delta"]
    fit = OLS(d.values, add_constant(z.values)).fit()
    grid = np.linspace(z.min(), z.max(), 100)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(z * 100, d * 100, alpha=0.4, s=14, color="#555555")
    ax.plot(grid * 100, (fit.params[0] + fit.params[1] * grid) * 100,
            color="red", lw=1.5, label=f"OLS slope = {fit.params[1]:.2f}")
    ax.set_xlabel("cycle_z (x100)", fontsize=9)
    ax.set_ylabel("cycle_delta (x100)", fontsize=9)
    ax.set_title(f"cycle_z vs cycle_delta (corr = {z.corr(d):+.3f})", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(lw=0.4, alpha=0.4)
    fig.tight_layout()
    _finish(fig, out_path)


def plot_irfs(var, out_path, H=20):
    """Responses of z and delta to each one-SD Cholesky structural shock."""
    A1, P = var["A1"], var["P"]
    names = ["z", "delta"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for row, shock in enumerate(names):
        irf = np.zeros((H + 1, 2))
        irf[0] = P[:, row]                       # impact = column of P
        for h in range(1, H + 1):
            irf[h] = A1 @ irf[h - 1]
        for col, resp in enumerate(names):
            ax = axes[row, col]
            ax.plot(np.arange(H + 1), irf[:, col] * 100,
                    color=SERIES[resp][1], lw=2.0, marker="o", ms=3)
            ax.axhline(0, color="black", lw=0.8)
            ax.set_title(f"{shock} shock  ->  {SERIES[resp][0]}", fontsize=9)
            ax.set_xlabel("Horizon h (quarters)", fontsize=8)
            ax.set_ylabel("Response (x100)", fontsize=8)
            ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle("VAR(1) IRFs — Cholesky structural shocks, z ordered first",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


# =============================================================================
#  MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("Part 6b — shock processes for z, delta, tau, s")
    print("=" * 78)

    dat = pd.read_pickle(RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    print(f"\n[1] raw_data.pkl: {len(dat)} obs, "
          f"lp {dat['lp'].first_valid_index()}–{dat['lp'].last_valid_index()}, "
          f"delta {dat['delta'].first_valid_index()}–"
          f"{dat['delta'].last_valid_index()}")

    levels = build_levels(dat)
    print(f"[2] Levels built, {levels.index[0]}–{levels.index[-1]}. "
          f"delta/tau = {levels['delta'].mean()/levels['tau'].mean():.4f}")

    cyc = hp_cycles(dat, levels)
    print(f"[3] Cycles: {len(cyc)} obs, {cyc.index[0]}–{cyc.index[-1]}. "
          f"cor(s, tau) = {cyc['s'].corr(cyc['tau']):.4f}")

    var = fit_var_cholesky(cyc)
    table = build_table(levels, cyc, var)
    report(table, var, cyc)

    out_csv = RESULTS_DIR / "shock_calibration.csv"
    table.assign(hp_lambda=HP_LAMBDA,
                 sample_start=str(cyc.index[0]),
                 sample_end=str(cyc.index[-1])).to_csv(out_csv)
    print(f"\n[4] Saved: {out_csv}")

    print(f"\n[5] Plotting  (interactive = {INTERACTIVE}, "
          f"backend = {matplotlib.get_backend()})")
    plot_cycle_panel(cyc, RESULTS_DIR / "cycles_panel.png")
    plot_separation_decomposition(levels[levels.index <= SAMPLE_END], cyc,
                                  RESULTS_DIR / "separation_decomposition.png")
    plot_scatter(cyc, RESULTS_DIR / "var_zd_scatter.png")
    plot_irfs(var, RESULTS_DIR / "var_zd_irfs.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
