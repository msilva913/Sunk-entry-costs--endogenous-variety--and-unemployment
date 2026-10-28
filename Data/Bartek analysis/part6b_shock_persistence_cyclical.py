"""
part6b_shock_persistence_cyclical.py
=====================================
Business-cycle processes for the model's three shocks, z, delta and s, at the
monthly frequency the model is solved at. Implements the sample windows of
decisions S13 and S14 (Oct 20, 2026) and the D14 choice of Oct 28, 2026 (s is
not orthogonalized), subject to the independence check in section 5 of the
report.

    z      labor productivity        (PRS85006163, index, normalized to 1)
    delta  permanent exit rate       (BED Deaths, employment-weighted)
    tau    total separation rate     (Shimer 2012 method, from CPS)
    s      match separation rate     (constructed, see SPLICE below)

The model takes three independent AR(1) processes with uncorrelated
innovations. "Independent" describes the model's processes, not the
estimator: delta comes out of a VAR, z and s out of univariate AR(1)s. The
text between METHOD and OUTPUT is printed at the top of the report, so the
.txt file carries its own logic.

METHOD
------
Windows (S13, S14). Every series ends at 2019Q4, cut BEFORE filtering: the HP
trend is two-sided, so filtering through 2020+ would leak the pandemic into the
2017-19 cycles. z, tau and s start in 1951Q1, the window of Shimer (2005),
Hagedorn-Manovskii (2008) and Coles-Kelishomi (2018). delta starts when BED
Deaths do, 1992Q3. A 1992-2019 window alone would sit entirely in the
post-1985 era, when labor productivity is acyclical, and make z look weak by
construction.

Step 1. z: univariate AR(1) on 1951Q1-2019Q4 cycles.
    z is ordered first wherever it enters a system, so its innovation is not
    purged of anything, and its cross-lag terms are insignificant. A
    univariate AR(1) is then the z equation of any of the VARs below, run on
    the longer window. It supplies rho_z and sigma_z.

Step 2. delta: bivariate (z, delta) VAR(1) on 1992Q3-2019Q4, Cholesky with z
    ordered first. Only rho_delta and sigma_delta are kept.
    - sigma_delta. The model generates z -> chi^c -> delta_e itself, through
      the continuation-cost threshold. If the exogenous delta innovation also
      carried its correlation with the z innovation, that channel would be
      counted twice. The Cholesky step strips the part of the delta
      innovation that the same-quarter z innovation explains. This is where
      we depart from Coles and Kelishomi (2018), who allow correlated
      innovations because their model has no endogenous exit margin.
    - rho_delta is the own-lag coefficient of the VAR, an OLS object formed
      before the Cholesky step, which therefore does not touch it.
    - The cross-lag terms are reported, not fed to the model: they are
      equilibrium outcomes, not primitive forcing.
    - The VAR's z equation only supplies the z innovation for the purge.
      rho_z and sigma_z come from Step 1, on the longer window.

Step 3. s: univariate AR(1) on the spliced 1951Q1-2019Q4 series. sigma_s is
    the AR(1) residual SD, NOT orthogonalized (D14).
    The double-counting argument of Step 2 does not carry over. The model has
    no endogenous s, so z never moves s in the model. Purging s of its
    correlation with z would remove separation variation the model can never
    regenerate. Persistence is the univariate AR(1) coefficient, not the
    diagonal of a larger VAR with the cross-lags dropped: the model's s is an
    AR(1) on its own, so its autocorrelation should match that of s in the
    data. The model also assumes the s and z innovations are uncorrelated.
    Section 5 of the report tests that, against a rule fixed in advance:
    revisit D14 if the z innovation explains at least 10 % of the s
    innovation variance and |t| >= 2.

Why s is not tau. Shimer's tau is the TOTAL separation rate: it counts workers
whose job vanished with the product line as well as broken matches on
surviving lines. The model's s is the second kind only. The accounting identity
behind tab:calib_targets links them:

    tau = delta_e + (1 - delta_e) s    =>    s = (tau - delta_e)/(1 - delta_e)

SPLICE (s before 1992Q3). BED Deaths begin in 1992Q3; tau begins in 1951Q1.

    t >= 1992Q3:  s_t = (tau_t - delta_t)/(1 - delta_t)    measured identity
    t <  1992Q3:  s_t = (tau_t - dbar)/(1 - dbar)          dbar a constant

    dbar = mean of the first SPLICE_ANCHOR_Q BED quarters (1992Q3-1994Q4).

The spliced series is cut at 2019Q4, logged and HP-filtered once, over the
whole 1951-2019 window. Log-linearizing at the means shows what each piece
does to the cycle:

    before 1992Q3:  d ln s = k_pre d ln tau,            k_pre = tau/(tau - dbar)
    from 1992Q3:    d ln s = k_post d ln tau - w d ln delta,
                    k_post = tau/(tau - delta),
                    w = delta (1 - tau)/((1 - delta)(tau - delta))

The report prints k_pre, k_post and w (section 4a). Three consequences:
    1. Before 1992Q3 the s cycle is tau's cycle times a near-constant factor.
    2. The level of dbar hardly matters for the cycle. A constant log shift
       goes into the HP trend, and so does a slow drift in the true delta
       (the secular fall in business dynamism). dbar only moves k_pre.
    3. The one thing lost is delta's own cycle before 1992, which enters s
       with weight w. No method built on tau alone can recover it, because tau
       does not contain it. Only a separate pre-1992 exit series could.

Why the early-quarter anchor, not the 1992-2019 mean. BED delta starts high
and falls. With the BED-window mean, s would jump at 1992Q3, and the HP
filter would turn that step into a spurious cycle on both sides of the
splice. Anchoring dbar at the start of BED nearly removes the step. The cost
is a slightly larger k_pre, which overstates sigma_Q a little on the BED
window, where the true s is known. Sections 4b and 4c show both anchors.

Rejected alternative: regress s on tau after 1992 and impute the fitted value.
In cycles, the fit is k * (tau cycle), the same object as a constant dbar up
to the scale k, so it leaves rho unchanged and only rescales sigma. On the BED
window it understates the innovation SD. In levels it extrapolates: tau
averages about 0.044 before 1992 and 0.024 after, so most pre-1992 values lie
outside the range the regression is fitted on, and the fit opens a larger
step at 1992Q3.

Checks reported:
    (a) k_pre, k_post and w.
    (b) The step in ln s at 1992Q3 under both anchors.
    (c) On the BED window, AR(1)s of the identity s against the constant-dbar
        versions and tau. This measures the pre-1992 approximation where the
        truth is known, under the untestable assumption that delta's cyclical
        role was similar before 1992.
    (d) AR(1)s of z and s on sub-windows of the full-sample cycles, including
        1985-2019 (S14 iv).
    (e) Independence (section 5): the correlation of the s and z innovations,
        against the D14 rule.
    (f) Ljung-Box tests on the residuals (section 6).

Other key points:
    - Pre-1992 LEVELS of s are imputed. Never use them for calibration targets
      (s-bar, delta_e/tau). The table reports the means of delta, tau and s
      over the BED window only.
    - tau fell by about half between the two eras. That is irrelevant for the
      cycles, but it matters for the tau-bar calibration target.
    - After 1992, BED measurement noise enters s with weight w. It lowers
      rho_s on that window relative to tau.
    - BED Deaths are noisy at quarterly frequency: the delta equation has a low
      R^2. Read rho_delta with caution.
    - The 1951-2019 window spans the Great Moderation. A single AR(1) per
      series averages two volatility regimes; section 4d shows the
      sub-windows.

Filtering and frequency. Log, then HP filter with lambda = 1,600, the quarterly
convention that observables_moments.py also uses, so the shock processes and
the moments see the same cyclical object. Each series is filtered over its own
window. Quarterly estimates convert to monthly following CK:
rho_m = rho_Q^(1/3), and sigma_m matches the unconditional variance,
sigma_m = sigma_Q * sqrt((1 - rho_Q^2)/(1 - rho_m^2)).

OUTPUT
------
    data/results/shock_calibration.csv          the table, one row per series
    data/results/shock_calibration.txt          the full report (N15)
    data/results/cycles_panel.png               HP cycles, all four series
    data/results/separation_decomposition.png   tau, delta, s: levels and cycles
    data/results/var_zd_scatter.png             cycle_z vs cycle_delta
    data/results/var_zd_irfs.png                VAR IRFs to structural shocks

Run
---
    python part6b_shock_persistence_cyclical.py   (environment: Data/environment.yml)

Prerequisites
-------------
    python observables.py    (produces raw_data.pkl)
"""

import os
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.diagnostic import acorr_ljungbox
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
HP_LAMBDA       = 1_600
SAMPLE_START    = "1951Q1"               # S14: z, tau, s and labor-market moments
SAMPLE_END      = "2019Q4"               # S13: pre-pandemic, cut before filtering
ROBUST_START    = "1985Q1"               # S14 (iv): robustness window
SPLICE_ANCHOR_Q = 10                     # dbar = mean of the first 10 BED quarters
RAW_DATA        = Path("../raw_data.pkl")    # built by observables.py
RESULTS_DIR     = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# D14 rule: revisit the independence of s if the z innovation explains at least
# this share of the s innovation variance AND the loading has |t| >= 2.
INDEP_R2_MAX = 0.10

# SHOCKS_BASE in run_solution_core.jl (pre-S13 part6b, Sept 21), for comparison.
SHOCKS_BASE = {"z": (0.902, 0.0092), "delta": (0.592, 0.0669), "s": (0.8741, 0.0854)}

# Display labels and colors, shared by the table and every plot.
SERIES = {
    "z":     ("Labor productivity (z)",        "#1f77b4"),
    "delta": ("Permanent exit rate (delta)",   "#d62728"),
    "tau":   ("Total separation rate (tau)",   "#7f7f7f"),
    "s":     ("Match separation rate (s)",     "#2ca02c"),
}

# Plausibility bands for the units guard. Both rate series must be MONTHLY
# fractions: tau sits near 0.03 and delta near 0.0023. Wide enough to admit a
# sensible revision, narrow enough to catch a percent/fraction or
# monthly/quarterly mix-up.
TAU_BOUNDS   = (0.010, 0.060)
DELTA_BOUNDS = (0.0005, 0.010)

# NBER recessions inside 1951-2019 (peak month, trough month).
REC_SPANS = [("1953-07-01", "1954-05-01"), ("1957-08-01", "1958-04-01"),
             ("1960-04-01", "1961-02-01"), ("1969-12-01", "1970-11-01"),
             ("1973-11-01", "1975-03-01"), ("1980-01-01", "1980-07-01"),
             ("1981-07-01", "1982-11-01"), ("1990-07-01", "1991-03-01"),
             ("2001-03-01", "2001-11-01"), ("2007-12-01", "2009-06-01")]

# Every line of the report goes through say(), which prints it and keeps it
# for shock_calibration.txt.
REPORT_LINES = []


def say(text=""):
    print(text)
    REPORT_LINES.append(str(text))


# =============================================================================
#  SERIES
# =============================================================================

def normalized_z(dat):
    """
    Labor productivity as an index with unit mean.

    Defined once and used everywhere, because z must be filtered over its whole
    1951- history: truncating it moves rho_z and would silently change the
    calibration in the draft.
    """
    return dat["lp"] / dat["lp"].mean()


def splice_point(delta):
    """
    The splice quarter (first BED observation) and the anchor dbar, the mean
    of the first SPLICE_ANCHOR_Q BED quarters. See SPLICE in the docstring.
    """
    delta = delta.dropna()
    return delta.index[0], float(delta.iloc[:SPLICE_ANCHOR_Q].mean())


def build_levels(dat):
    """
    Assemble the rate series on tau's window: tau, delta (BED window only),
    delta_used (dbar before the splice, measured delta from it) and the
    spliced s.

    From observables.py both inputs are already monthly fractions observed
    quarterly: tau = 1 - exp(-s_h) is a monthly separation probability, and
    delta is the quarterly BED deaths employment share divided by 3. The
    identity is a within-period accounting statement at the model's
    frequency, so it needs both on that footing -- hence the guard.

    z is not included. It carries no informative level; see normalized_z.
    """
    tau, delta = dat["s"].dropna(), dat["delta"].dropna()   # observables' 's' IS tau
    for name, series, bounds in [("tau", tau, TAU_BOUNDS),
                                 ("delta", delta, DELTA_BOUNDS)]:
        if not bounds[0] <= series.mean() <= bounds[1]:
            raise ValueError(
                f"{name} averages {series.mean():.5f}, outside the "
                f"monthly-fraction band {bounds}. Check that observables.py "
                f"still delivers a monthly rate as a fraction."
            )

    splice, dbar = splice_point(delta)
    df = pd.DataFrame({"tau": tau}).join(delta.rename("delta"))
    bed = (df.index >= splice) & (df.index <= SAMPLE_END)
    if df.loc[bed, "delta"].isna().any():
        raise ValueError("BED delta has gaps inside its window; the identity "
                         "would silently switch to dbar there.")

    # Before the splice dbar stands in for delta. From the splice on, the
    # measured delta is used, and it stays NaN after BED ends (past 2019Q4).
    df["delta_used"] = df["delta"].where(df.index >= splice, dbar)
    df["s"] = (df["tau"] - df["delta_used"]) / (1.0 - df["delta_used"])
    if (df["s"] <= 0).any():
        first = df.index[df["s"] <= 0][0]
        raise ValueError(
            f"Match separation rate is non-positive, first at {first}: delta "
            f"exceeds tau there, so the two series are measured on "
            f"inconsistent bases."
        )

    return df[["delta", "delta_used", "tau", "s"]]


def window(series, start=SAMPLE_START, end=SAMPLE_END):
    series = series.dropna()
    return series[(series.index >= start) & (series.index <= end)]


def hp_log_cycle(series):
    """Log-HP cycle of a level series, over the series' own (cut) window."""
    cycle, _ = hpfilter(np.log(series), lamb=HP_LAMBDA)
    return cycle


def hp_cycles(dat, levels):
    """
    Log-HP cycles, one column per series, each filtered over its own window.

    z, tau and s are filtered over SAMPLE_START-SAMPLE_END; delta over its BED
    window, so it is NaN before 1992Q3. Every series is cut at SAMPLE_END
    *before* filtering (S13).
    """
    sources = {"z": normalized_z(dat)} | {c: levels[c] for c in ("delta", "tau", "s")}
    return pd.DataFrame({name: hp_log_cycle(window(series))
                         for name, series in sources.items()})


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
            "nobs":   int(fit.nobs),
            "resid":  fit.resid,
            "start":  df.index[0], "end": df.index[-1]}


def fit_var_cholesky(cyc):
    """
    Bivariate VAR(1) on [cycle_z, cycle_delta] over the BED window, Cholesky
    with z ordered first (Step 2).

    Returns the companion matrix, the impact matrix P, and the own persistence
    and innovation SDs. The diagonal of P is the structural shock SD in
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
            "A1_t": np.array([f.tvalues[1:] for f in fits]),
            "corr_uu": float(Sigma_u[0, 1] /
                             np.sqrt(Sigma_u[0, 0] * Sigma_u[1, 1])),
            "rho":    {"z": float(A1[0, 0]), "delta": float(A1[1, 1])},
            "rho_se": {"z": float(fits[0].bse[1]),
                       "delta": float(fits[1].bse[2])},
            "r2":     {"z": float(fits[0].rsquared),
                       "delta": float(fits[1].rsquared)},
            # Innovation SDs are the diagonal of P. For delta it is the part
            # orthogonal to the contemporaneous z innovation, the only number
            # in the table the Cholesky changes.
            "sigma": {"z": float(P[0, 0]), "delta": float(P[1, 1])},
            "sigma_rf": {"z": float(np.sqrt(Sigma_u[0, 0])),
                         "delta": float(np.sqrt(Sigma_u[1, 1]))},
            "resid": pd.DataFrame(resid, index=df.index, columns=["z", "delta"]),
            "nobs": int(fits[0].nobs),
            "start": df.index[0], "end": df.index[-1]}


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


def build_table(levels, cyc, fits, var, bed_start):
    """
    One row per series: mean, cyclical SD, and the AR(1) process at both
    frequencies.

      rho    z, tau, s  univariate AR(1), 1951-2019 (Steps 1 and 3)
             delta      own lag of the (z, delta) VAR, BED window (Step 2)
      sigma  z, tau, s  AR(1) residual SD, not orthogonalized
             delta      VAR residual SD orthogonalized against the z innovation

    rho_Q_uni is the univariate AR(1) for every series, so the VAR's effect on
    rho_delta is visible; sigma_Q_rf keeps delta's SD before the Cholesky step.
    mean is over the BED window for the rate series: pre-1992 levels of s are
    imputed and must not feed a calibration target.
    """
    rows = []
    for name in SERIES:
        uni = fits[name]
        if name == "delta":
            est = {k: var[k]["delta"] for k in ("rho", "rho_se", "r2", "sigma", "sigma_rf")}
            est |= {"nobs": var["nobs"], "start": var["start"], "end": var["end"]}
            source = "VAR(1) z->delta, sigma perp z"
        else:
            est = {"rho": uni["rho_Q"], "rho_se": uni["rho_se"], "r2": uni["r2"],
                   "sigma": uni["sigma_Q"], "sigma_rf": uni["sigma_Q"],
                   "nobs": uni["nobs"], "start": uni["start"], "end": uni["end"]}
            source = "AR(1), spliced s" if name == "s" else "AR(1)"

        rho_m, sigma_m = to_monthly(est["rho"], est["sigma"])
        rows.append({
            "series":       name,
            "mean":         1.0 if name == "z" else
                            float(window(levels[name], start=bed_start).mean()),
            "units":        "index (normalized)" if name == "z" else "fraction/month",
            "sd_cycle":     float(window(cyc[name], est["start"], est["end"]).std()),
            "rho_Q":        est["rho"],
            "sigma_Q":      est["sigma"],
            "rho_m":        rho_m,
            "sigma_m":      sigma_m,
            "source":       source,
            "rho_se":       est["rho_se"],
            "r2":           est["r2"],
            "nobs":         est["nobs"],
            "rho_Q_uni":    uni["rho_Q"],
            "sigma_Q_rf":   est["sigma_rf"],
            "hp_lambda":    HP_LAMBDA,
            "sample_start": str(est["start"]),
            "sample_end":   str(est["end"]),
        })
    return pd.DataFrame(rows).set_index("series")


# =============================================================================
#  SPLICE AND INDEPENDENCE CHECKS
# =============================================================================

def loglin_weights(levels, splice, dbar):
    """k_pre, k_post and w of the SPLICE section, evaluated at window means."""
    pre = window(levels["tau"], end=str(splice - 1))
    bed = levels.loc[(levels.index >= splice) & (levels.index <= SAMPLE_END)]
    t, d = bed["tau"].mean(), bed["delta"].mean()
    return {"k_pre":  float(pre.mean() / (pre.mean() - dbar)),
            "k_post": float(t / (t - d)),
            "w":      float(d * (1 - t) / ((1 - d) * (t - d))),
            "tau_pre": float(pre.mean()), "tau_post": float(t)}


def splice_step(levels, splice, dbar_anchor):
    """
    The jump in ln s at the splice quarter if dbar were used there instead of
    the measured delta: ln s(identity) - ln s(constant dbar).
    """
    t, d = levels.loc[splice, "tau"], levels.loc[splice, "delta"]
    return float(np.log((t - d) / (1 - d)) - np.log((t - dbar_anchor) / (1 - dbar_anchor)))


def splice_test(levels, splice, dbar_early, dbar_full):
    """
    Check (c): on the BED window, where the true s is known, compare the
    identity s with what the pre-1992 construction would have produced.
    Each series is filtered over the BED window.
    """
    bed = window(levels[["tau", "delta"]], start=splice)
    candidates = {
        "identity s (truth)":            (bed["tau"] - bed["delta"]) / (1 - bed["delta"]),
        f"const dbar, first {SPLICE_ANCHOR_Q} q (used)":
                                         (bed["tau"] - dbar_early) / (1 - dbar_early),
        "const dbar, BED-window mean":   (bed["tau"] - dbar_full) / (1 - dbar_full),
        "tau itself":                    bed["tau"],
    }
    truth = hp_log_cycle(candidates["identity s (truth)"])
    rows = {}
    for label, series in candidates.items():
        cycle = hp_log_cycle(series)
        fit = fit_ar1(cycle)
        rows[label] = {"sd_cycle": cycle.std(), "rho_Q": fit["rho_Q"],
                       "sigma_Q": fit["sigma_Q"], "cor_truth": cycle.corr(truth)}
    return pd.DataFrame(rows).T


def subwindow_ar1(cyc, splice):
    """Check (d): AR(1)s on sub-windows of the full-sample cycles."""
    spans = {f"{SAMPLE_START}-{splice - 1}": (SAMPLE_START, str(splice - 1)),
             f"{splice}-{SAMPLE_END}":       (str(splice), SAMPLE_END),
             f"{ROBUST_START}-{SAMPLE_END}": (ROBUST_START, SAMPLE_END),
             f"{SAMPLE_START}-{SAMPLE_END}": (SAMPLE_START, SAMPLE_END)}
    rows = {}
    for name in ("z", "s"):
        for label, (lo, hi) in spans.items():
            fit = fit_ar1(window(cyc[name], lo, hi))
            rows[(name, label)] = {"rho_Q": fit["rho_Q"], "sigma_Q": fit["sigma_Q"],
                                   "nobs": fit["nobs"]}
    return pd.DataFrame(rows).T


def independence(fits, var):
    """
    Check (e): do the s and z innovations correlate? Regress the s AR(1)
    residual on the z AR(1) residual over 1951-2019, and on the z and delta
    VAR residuals over the BED window.
    """
    out = {}
    df = pd.concat([fits["s"]["resid"].rename("s"), fits["z"]["resid"].rename("z"),
                    fits["tau"]["resid"].rename("tau")], axis=1).dropna()

    def on_z(d):
        fit = OLS(d["s"], add_constant(d["z"])).fit(cov_type="HC3")
        return {"corr": float(d["s"].corr(d["z"])), "r2": float(fit.rsquared),
                "t": float(fit.tvalues["z"]), "b": float(fit.params["z"]),
                "corr_tau": float(d["tau"].corr(d["z"])), "nobs": int(fit.nobs),
                "start": d.index[0], "end": d.index[-1]}

    out["z"] = on_z(df)
    # Era split at ROBUST_START: the same residuals, from the full-window fits.
    out["eras"] = [on_z(df[df.index < ROBUST_START]), on_z(df[df.index >= ROBUST_START])]

    df = pd.concat([fits["s"]["resid"].rename("s"), var["resid"]], axis=1).dropna()
    fit = OLS(df["s"], add_constant(df[["z", "delta"]])).fit(cov_type="HC3")
    out["z_delta"] = {"r2": float(fit.rsquared), "t_z": float(fit.tvalues["z"]),
                      "t_delta": float(fit.tvalues["delta"]), "nobs": int(fit.nobs),
                      "start": df.index[0], "end": df.index[-1]}
    return out


def ljung_box(resid, lags=(4, 8)):
    lb = acorr_ljungbox(np.asarray(resid), lags=list(lags))
    return [float(p) for p in lb["lb_pvalue"]]


# =============================================================================
#  REPORT
# =============================================================================

def report_method():
    method = __doc__.split("METHOD\n------\n")[1].split("OUTPUT\n------")[0]
    say(f"\n{'-'*78}\n0. Method\n{'-'*78}")
    for line in method.rstrip().splitlines():
        say(line)


def report_table(table):
    fmt = {"mean": "{:.5f}".format, "sd_cycle": "{:.4f}".format,
           "rho_m": "{:.4f}".format, "sigma_m": "{:.6f}".format}
    say(f"\n{'='*78}\n2. Shock processes, monthly (HP({HP_LAMBDA:,}) log cycles)\n{'='*78}")
    say(table[["mean", "sd_cycle", "rho_m", "sigma_m", "sample_start", "sample_end",
               "source"]].to_string(formatters=fmt))
    say("\n  tau is an outcome of delta and s in the model, not a shock; it is shown"
        "\n  because it is the series the moment table targets directly."
        "\n  mean: delta, tau and s over the BED window only (pre-1992 s is imputed).")

    say(f"\n  Quarterly estimates")
    say(table[["rho_Q", "rho_se", "rho_Q_uni", "sigma_Q", "sigma_Q_rf", "r2",
               "nobs"]].to_string(float_format=lambda x: f"{x:.4f}"))

    say("\n  Against SHOCKS_BASE (run_solution_core.jl; pre-S13 part6b, 1992-2019 window)")
    say(f"  {'':8s}{'rho_m old':>11s}{'rho_m new':>11s}{'sigma_m old':>13s}{'sigma_m new':>13s}")
    for name, (rho_old, sig_old) in SHOCKS_BASE.items():
        r = table.loc[name]
        say(f"  {name:8s}{rho_old:11.4f}{r['rho_m']:11.4f}{sig_old:13.4f}{r['sigma_m']:13.4f}")


def report_var(table, var):
    d, A1, T = table.loc["delta"], var["A1"], var["A1_t"]
    say(f"\n{'-'*78}\n3. Step 2 detail: (z, delta) VAR(1), {var['start']}-{var['end']}, "
        f"T = {var['nobs']}\n{'-'*78}")
    say(f"  sigma_Q of delta before -> after the Cholesky step: "
        f"{d['sigma_Q_rf']:.4f} -> {d['sigma_Q']:.4f} "
        f"({100*(1 - d['sigma_Q']/d['sigma_Q_rf']):.1f}% smaller)."
        f"\n  That is the entire footprint of the Cholesky step; persistence is "
        f"untouched by it.")
    say(f"  rho_Q of delta: VAR own lag {d['rho_Q']:.4f}, univariate AR(1) "
        f"{d['rho_Q_uni']:.4f}.")
    say("\n  Off-diagonals (reported, not fed to the model; HC3 t in parentheses):")
    say(f"    delta <- L.z     = {A1[1,0]:+.4f} ({T[1,0]:+.2f})  endogenous exit response to z")
    say(f"    z <- L.delta     = {A1[0,1]:+.4f} ({T[0,1]:+.2f})  expected near zero")
    say(f"    corr(u^z, u^delta) = {var['corr_uu']:+.4f}  before the Cholesky step")
    say(f"  VAR z equation on this window (not used): rho_Q {var['rho']['z']:.4f}, "
        f"sigma_Q {var['sigma']['z']:.4f}.")
    if d["r2"] < 0.15:
        say(f"\n  WARNING: R^2 for the delta equation is {d['r2']:.4f}. BED Deaths is "
            f"noisy at quarterly\n  frequency (N = {int(d['nobs'])}); read rho_delta "
            f"with caution.")


def report_splice(levels, splice, dbar, dbar_full, cyc):
    say(f"\n{'-'*78}\n4. Splice checks\n{'-'*78}")
    say(f"  Splice quarter {splice}. dbar (used) = {dbar:.5f}, mean of "
        f"{splice}-{splice + SPLICE_ANCHOR_Q - 1};"
        f"\n  BED-window mean {dbar_full:.5f} (shown for comparison).")

    lw = loglin_weights(levels, splice, dbar)
    say(f"\n  (a) Log-linear weights. tau averages {lw['tau_pre']:.4f} before the splice "
        f"and {lw['tau_post']:.4f} after.")
    say(f"      k_pre  = {lw['k_pre']:.4f}   s cycle = k_pre * tau cycle before {splice}")
    say(f"      k_post = {lw['k_post']:.4f}   weight on the tau cycle from {splice}")
    say(f"      w      = {lw['w']:.4f}   weight on the delta cycle from {splice}; "
        f"the part lost before it")

    say(f"\n  (b) Step in ln s at {splice} if dbar replaced the measured delta there:")
    say(f"      dbar = first-{SPLICE_ANCHOR_Q}-quarter mean (used): "
        f"{100*splice_step(levels, splice, dbar):+.2f}%")
    say(f"      dbar = BED-window mean:              "
        f"{100*splice_step(levels, splice, dbar_full):+.2f}%")
    say(f"      For scale: SD of the s cycle = {100*cyc['s'].std():.2f}%.")

    say(f"\n  (c) BED window ({splice}-{SAMPLE_END}), where the true s is known. "
        f"Each series filtered\n      over that window; quarterly AR(1).")
    say(splice_test(levels, splice, dbar, dbar_full)
        .to_string(float_format=lambda x: f"{x:.4f}"))
    say("      The constant-dbar rows show what the pre-1992 construction does when"
        "\n      the truth is known. The gap to the identity row in rho_Q is delta's own"
        "\n      cycle (including BED noise), which no tau-based construction recovers.")

    say("\n  (d) Sub-windows of the full-sample cycles, quarterly AR(1)")
    say(subwindow_ar1(cyc, splice).to_string(
        float_format=lambda x: f"{x:.4f}" if abs(x) < 50 else f"{x:.0f}"))


def report_independence(fits, var):
    ind = independence(fits, var)
    z, zd = ind["z"], ind["z_delta"]
    say(f"\n{'-'*78}\n5. Independence of the s innovation (D14)\n{'-'*78}")
    say(f"  s AR(1) residual on z AR(1) residual, {z['start']}-{z['end']}, "
        f"T = {z['nobs']}:")
    say(f"      corr {z['corr']:+.4f}, share explained R^2 = {100*z['r2']:.1f}%, "
        f"HC3 t = {z['t']:+.2f}")
    say("  Same residuals, split by era (corr_tau: the tau residual on z, which shows"
        "\n  whether the splice drives the result; before 1992 s is k_pre * tau):")
    for e in ind["eras"]:
        say(f"      {e['start']}-{e['end']}  T = {e['nobs']:3d}  corr {e['corr']:+.3f}  "
            f"R^2 {100*e['r2']:5.1f}%  t {e['t']:+.2f}  loading {e['b']:+.2f}  "
            f"corr_tau {e['corr_tau']:+.3f}")
    say(f"  s AR(1) residual on the (z, delta) VAR residuals, {zd['start']}-{zd['end']}, "
        f"T = {zd['nobs']}:")
    say(f"      share explained R^2 = {100*zd['r2']:.1f}%, t(z) = {zd['t_z']:+.2f}, "
        f"t(delta) = {zd['t_delta']:+.2f}")
    revisit = z["r2"] >= INDEP_R2_MAX and abs(z["t"]) >= 2
    say(f"  D14 rule (revisit if R^2 >= {100*INDEP_R2_MAX:.0f}% and |t| >= 2): "
        f"{'REVISIT: z loads on s' if revisit else 'independence holds; keep sigma_s unpurged'}.")


def report_residuals(fits, var):
    say(f"\n{'-'*78}\n6. Ljung-Box p-values on residuals (lags 4, 8)\n{'-'*78}")
    for label, resid in [("z     AR(1)", fits["z"]["resid"]),
                         ("delta VAR  ", var["resid"]["delta"]),
                         ("s     AR(1)", fits["s"]["resid"])]:
        p4, p8 = ljung_box(resid)
        say(f"  {label}  {p4:.3f}  {p8:.3f}")


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


def _decorate(ax, splice=None):
    ax.axhline(0, color="black", lw=0.7)
    for lo, hi in REC_SPANS:
        ax.axvspan(pd.Timestamp(lo), pd.Timestamp(hi), alpha=0.12, color="grey")
    if splice is not None:
        ax.axvline(splice.to_timestamp(), color="black", ls=":", lw=1.0)
    ax.grid(axis="y", lw=0.4, alpha=0.4)


def plot_cycle_panel(cyc, splice, out_path):
    """One stacked panel per series, each on its own estimation window."""
    fig, axes = plt.subplots(len(SERIES), 1, figsize=(11, 9), sharex=True)
    dts = _ts(cyc.index)
    for ax, (name, (label, color)) in zip(axes, SERIES.items()):
        ax.plot(dts, cyc[name].values * 100, color=color, lw=1.0)
        ax.set_ylabel("Log-cycle (x100)", fontsize=8)
        ax.set_title(label, fontsize=9)
        _decorate(ax, splice)
    fig.suptitle(f"HP-filtered log cycles (lambda = {HP_LAMBDA:,}); "
                 f"dotted line: BED start / s splice", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


def plot_separation_decomposition(levels, cyc, splice, out_path):
    """
    Overlay of tau, delta and s. The top panel shows the split in levels, with
    the dbar stand-in before the splice; the bottom shows the cycles the AR(1)
    is fitted to.
    """
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    lv = levels[(levels.index >= SAMPLE_START) & (levels.index <= SAMPLE_END)]

    ax = axes[0]
    dts = _ts(lv.index)
    for name, ls in [("tau", "-"), ("delta", "-"), ("s", "--")]:
        ax.plot(dts, lv[name].values * 100, color=SERIES[name][1], ls=ls, lw=1.2,
                label=SERIES[name][0])
    pre = lv["delta"].isna()
    ax.plot(_ts(lv.index[pre]), lv.loc[pre, "delta_used"].values * 100,
            color=SERIES["delta"][1], ls=":", lw=1.2, label="dbar (splice stand-in)")
    ax.set_ylabel("Percent per month", fontsize=9)
    ax.set_title("Levels", fontsize=10)

    ax = axes[1]
    dts = _ts(cyc.index)
    for name, ls in [("tau", "-"), ("delta", "-"), ("s", "--")]:
        ax.plot(dts, cyc[name].values * 100, color=SERIES[name][1], ls=ls, lw=1.0,
                label=SERIES[name][0])
    ax.set_ylabel("Log-cycle (x100)", fontsize=9)
    ax.set_title(f"HP-filtered log cycles (lambda = {HP_LAMBDA:,})", fontsize=10)

    for ax in axes:
        _decorate(ax, splice)
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    fig.suptitle("Separation rate decomposition: tau = delta + (1 - delta) s; "
                 f"s spliced at {splice}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


def plot_scatter(cyc, out_path):
    """cycle_z against cycle_delta on the BED window, with the OLS fit."""
    df = cyc[["z", "delta"]].dropna()
    z, d = df["z"], df["delta"]
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

    fig.suptitle("Step 2 VAR(1) IRFs, Cholesky structural shocks, z ordered first",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _finish(fig, out_path)


# =============================================================================
#  RUN
# =============================================================================

if __name__ == "__main__":
    REPORT_LINES.clear()
    say("=" * 78)
    say("Part 6b — shock processes for z, delta and s  (part6b_shock_persistence_cyclical.py)")
    say("=" * 78)
    report_method()

    # ── 1. data ──────────────────────────────────────────────────────────────
    dat = pd.read_pickle(RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    levels = build_levels(dat)
    splice, dbar = splice_point(dat["delta"])
    dbar_full = float(window(levels["delta"], start=splice).mean())
    cyc = hp_cycles(dat, levels)

    say(f"\n{'-'*78}\n1. Data\n{'-'*78}")
    say(f"  raw_data.pkl: lp {dat['lp'].first_valid_index()}-{dat['lp'].last_valid_index()}, "
        f"tau {dat['s'].first_valid_index()}-{dat['s'].last_valid_index()}, "
        f"delta {dat['delta'].first_valid_index()}-{dat['delta'].last_valid_index()}")
    say(f"  Windows: z, tau, s {SAMPLE_START}-{SAMPLE_END}; delta {splice}-{SAMPLE_END}. "
        f"All cut at {SAMPLE_END} before filtering.")
    say(f"  s spliced at {splice} with dbar = {dbar:.5f}. On the BED window: "
        f"delta/tau = {dbar_full / window(levels['tau'], start=splice).mean():.4f}, "
        f"cor(s, tau) cycles = {window(cyc['s'], start=splice).corr(window(cyc['tau'], start=splice)):.4f}.")

    # ── Step 1: z, univariate AR(1) on 1951-2019 ────────────────────────────
    # ── Step 2: delta, (z, delta) VAR(1) on the BED window, z first ─────────
    # ── Step 3: s, univariate AR(1) on the spliced 1951-2019 series ─────────
    # delta also gets a univariate fit, only as the rho_Q_uni cross-check.
    fits = {name: fit_ar1(cyc[name]) for name in SERIES}
    var = fit_var_cholesky(cyc)
    table = build_table(levels, cyc, fits, var, splice)

    report_table(table)
    report_var(table, var)
    report_splice(levels, splice, dbar, dbar_full, cyc)
    report_independence(fits, var)
    report_residuals(fits, var)

    # ── outputs ─────────────────────────────────────────────────────────────
    out_csv = RESULTS_DIR / "shock_calibration.csv"
    out_txt = RESULTS_DIR / "shock_calibration.txt"
    table.to_csv(out_csv)
    out_txt.write_text("\n".join(REPORT_LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {out_csv}\nSaved: {out_txt}")

    print(f"\nPlotting  (interactive = {INTERACTIVE}, backend = {matplotlib.get_backend()})")
    plot_cycle_panel(cyc, splice, RESULTS_DIR / "cycles_panel.png")
    plot_separation_decomposition(levels, cyc, splice,
                                  RESULTS_DIR / "separation_decomposition.png")
    plot_scatter(cyc, RESULTS_DIR / "var_zd_scatter.png")
    plot_irfs(var, RESULTS_DIR / "var_zd_irfs.png")

    print("\nDone.")
