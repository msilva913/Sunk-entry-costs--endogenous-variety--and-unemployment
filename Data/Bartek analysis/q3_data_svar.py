"""
q3_data_svar.py
===============
Data side of Block B' (Q3, decision D13): recursive SVARs on the model's three
shocks (z, delta, s) and the labor-market outcomes (u, v). Reports how much of
the variance of u and v each shock explains, the impulse responses of u and v
with bootstrap bands, and a short list of key takeaways (section 6). The text between METHOD and OUTPUT is printed at the top of the report.

METHOD
------
Series. HP(1,600) log cycles, quarterly, every series cut at 2019Q4 BEFORE
filtering (S13) and filtered over its own window (S14):
    z      labor productivity (output per worker)      1951Q1-2019Q4
    delta  permanent exit rate, BED Deaths              1992Q3-2019Q4
    s      match separation rate, spliced before 1992Q3 1951Q1-2019Q4
    u      unemployment rate                            1951Q1-2019Q4
    v      vacancy rate (HWI spliced with JOLTS)        1951Q1-2019Q4
z, delta and s come from part6b_shock_persistence_cyclical.py, so the SVAR sees
exactly the cycles behind the model's shock processes.

Three SVARs:
    A  (z, s, delta, u, v)  1992Q3-2019Q4  delta's share; the full B' target
    B  (z, s, u, v)         1951Q1-2019Q4  z and s shares on the S14 window
    C  (z, s, u, v)         1992Q3-2019Q4  A against C measures what omitting
                                           delta does to the z and s shares
delta cannot enter B: before 1992Q3 it has no cycle (a constant dbar has an HP
cycle of exactly zero), so it would distort every coefficient on delta. In B
and C, delta's contribution is spread over the other shocks. A against C on the
same window measures that shift; carrying it back to 1951-2019 assumes delta's
role was similar before 1992, which cannot be tested. Before 1992Q3 the s cycle
is tau's cycle scaled by k_pre (part6b), so in B the s shock also carries a
small piece of delta (weight w of about 0.09).

Financial conditions (section 4b). Two more versions of A add a financial
conditions index F after z and before s and delta: (z, F, s, delta, u, v),
1992Q3-2019Q4, F fixed in second place. delta's shock is then exit net of
same-quarter and lagged financial conditions. 
Weighed average of about 100 financial indicators (risk spreads, volatilities, funding
conditions, leverage measures), with weights from dynamic factor model.
    A_F   F = NFCI, the Chicago Fed's headline index. It also tightens when the
          economy weakens, so it strips the broad financial-and-demand bundle
          from delta: what remains is a lower bound on pure destruction.
    A_Fa  F = ANFCI, the index adjusted for current macro conditions: delta net
          of financial shocks proper.
F is a mean-zero index, so it is HP-filtered in levels (not logs), cut at
2019Q4 first, like the other series. Source: data/cache/nfci_quarterly.parquet
(part7_nfci.py; quarterly averages of the weekly series, 1990Q1 on). The model
has no financial variable, so these SVARs inform what the model's delta should
be asked to match; they are not Block B' targets.

Mediation through delta (section 4c). How much of a shock's effect on u and
v runs through exit? Counterfactual impulse responses in the spirit of
Bernanke, Gertler and Watson (1997) and Sims and Zha (2006) (citations
unverified): delta's response is held at zero by adding, each quarter, the
delta shock that offsets it. delta is ordered after F and s, so the
offsetting shocks move only delta, u and v within the quarter. The indirect
effect through delta is the full response minus the counterfactual; the
direct effect (all other channels) is the counterfactual. Reported for the F
shock (A_F, A_Fa) and the z shock (A, A_F), with s as an alternative mediator,
as the peak and the cumulative (h = 0-12) response of u and v and the share
running through the mediator, with 68 % bootstrap bands. A reduced-form
counterfactual: it assumes the VAR's dynamics are unchanged when delta is
held fixed.

Filter robustness (section 7). Every SVAR is re-estimated on Hamilton
cycles (h = 8, p = 4; part6b's hamilton_cycle) instead of HP cycles. The
Hamilton filter is one-sided, so it does not blur the recursive timing
assumption, but it drops the first 11 quarters of each series (A starts in
1995Q2, T = 99). HP remains the baseline of sections 2-6.

Estimation. Reduced-form VAR(p) with a constant, equation by equation OLS on a
common sample. p is chosen by BIC over 1-4 (common sample of the VAR(4)); p = 1
and p = 2 are always reported as well.

Identification. Recursive (Cholesky), the shock block before (u, v), z always
first: technology is the most plausibly exogenous series, and the model
generates z -> delta_e within the period, so movements in delta caused by z
count toward z's share.
    Preferred (A): z, s, delta, u, v. Model-consistent: s is purely exogenous
        in the model, while measured exit delta_e responds within the period
        to aggregate conditions, so delta goes after s.
    Alternative (A): z, delta, s, u, v, for completeness. It assigns the
        same-quarter delta-s comovement to delta instead of s (measured s is
        built from delta, so BED noise is common to both). Only the split
        between delta and s changes; z, delta + s and the residual do not.
    B and C contain only z and s, so z first leaves one ordering.
The shares of the shock-block innovations do not depend on how u and v are
ordered among themselves, so the report combines them as "Residual".

Variance shares. Forecast-error variance decomposition at h = 4 and 8 quarters, and
at h = 200 as the unconditional share. The headline shows h = 4 and h = 200: with
HP cycles, h = 8 is nearly identical to h = 200, while h = 4 shows the short-run
composition (s acts fastest). h = 200 is (the model-side object:
share of the variance of HP-filtered u by shock). The residual u and v
innovations collect everything the three shocks do not.

Impulse responses. Preferred ordering, one-SD structural shocks, responses of u
and v (log points x 100) at h = 0-12 to the delta and s shocks, the Block B'
targets. Also z. Bands: residual bootstrap (recursive design, initial values
fixed, centered residuals resampled with replacement), N_BOOT draws, 68 % and
90 % percentile bands.

Y_t = XB + U
Y: T x k
Y_t: (T-p) x k
X: (T-p) x (1+kp)
B: (1+kp) x k

Residuals: U = Y_t − X B, of size (T − p) × k.
Covariance: Σ = UᵀU / (T − p − (1 + kp)), using the degrees-of-freedom correction.
Identification: the Cholesky factor of Σ, in the chosen ordering, gives the impact matrix P.
Downstream: the A_j matrices give the moving-average weights Ψh, and the impulse responses are Θh = Ψ_h P.
Variance shares: sums of squared entries of the Θ_h.


Caveats.
    - Recursive identification credits the earlier-ordered shocks with all
      same-quarter comovement. The ordering bounds show how much that matters.
    - Before 1992 the s shock carries a small part of delta (see above).
    - BED Deaths are noisy at quarterly frequency, so delta's equation fits
      poorly (part6b); its share is estimated imprecisely.
    - The model counterpart is the same VAR on simulated data (indirect
      inference). This program is the data side only.

OUTPUT
------
    data/results/q3_data_svar.txt        the full report (N15)
    data/results/q3_data_svar_headline.csv  headline shares of var(u), var(v)
    data/results/q3_data_svar_headline.txt  the same as plain tables, ready to paste
    data/results/q3_data_svar_fevd.csv   shares by SVAR, lag, ordering, horizon (incl. A_F, A_Fa)
    data/results/q3_data_svar_irf.csv    IRFs with bootstrap bands (filter column: hp, hamilton)
    data/results/q3_data_svar_irf.png    IRFs of u and v to the z, delta, s shocks
    data/results/q3_data_svar_levels.png the series in levels, with HP trends

Run
---
    python q3_data_svar.py      (environment: Data/environment.yml)
"""

import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    if not Path("q3_data_svar.py").exists():
        raise RuntimeError("Run from 'Data/Bartek analysis' or as a whole file.")

import part6b_shock_persistence_cyclical as p6b

# After part6b on purpose: part6b selects the backend (Agg as a script, inline in
# Spyder) before pyplot is first imported, and pyplot keeps that backend.
import matplotlib.pyplot as plt


# ── parameters ────────────────────────────────────────────────────────────────
BED_START  = "1992Q3"
HORIZONS   = (4, 8, 200)                  # 200 = unconditional
HEADLINE_H = (4, 200)                     # headline: short run and long run; 8 in the detail
H_IRF      = 12
P_MAX      = 4
N_BOOT     = 1_000
SEED       = 20261029
SHOCKS     = ("z", "delta", "s")
OUTCOMES   = ("u", "v")
SVARS = {
    "A": {"vars": ("z", "s", "delta", "u", "v"), "start": BED_START},  # preferred order
    "B": {"vars": ("z", "s", "u", "v"),          "start": p6b.SAMPLE_START},
    "C": {"vars": ("z", "s", "u", "v"),          "start": BED_START, "align": "A"},  # A's window
}
# A with a financial conditions index fixed after z (section 4b; not in the headline).
SVARS_F = {
    "A_F":  {"vars": ("z", "nfci", "s", "delta", "u", "v"),  "start": BED_START, "fixed": ("z", "nfci")},
    "A_Fa": {"vars": ("z", "anfci", "s", "delta", "u", "v"), "start": BED_START, "fixed": ("z", "anfci")},
}
NFCI_CACHE = Path("data/cache/nfci_quarterly.parquet")
OUT_TXT  = p6b.RESULTS_DIR / "q3_data_svar.txt"
OUT_FEVD = p6b.RESULTS_DIR / "q3_data_svar_fevd.csv"
OUT_HEAD = p6b.RESULTS_DIR / "q3_data_svar_headline.csv"
OUT_HEAD_TXT = p6b.RESULTS_DIR / "q3_data_svar_headline.txt"
OUT_IRF  = p6b.RESULTS_DIR / "q3_data_svar_irf.csv"
OUT_PNG  = p6b.RESULTS_DIR / "q3_data_svar_irf.png"
OUT_LEV  = p6b.RESULTS_DIR / "q3_data_svar_levels.png"

ECHO_METHOD = False     # True also prints the method text to the console (always in the .txt)

REPORT_LINES = []


def say(text="", echo=True):
    if echo:
        print(text)
    REPORT_LINES.append(str(text))


# =============================================================================
#  DATA
# =============================================================================

def build_levels_all():
    """
    The five series in levels on 1951Q1-2019Q4, before any filtering: log
    labor productivity, and delta, tau, s, u, v in percent (delta, tau and s
    monthly). s is the spliced series of part6b.
    """
    dat = pd.read_pickle(p6b.RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    lev = p6b.build_levels(dat)
    out = pd.DataFrame({"log z": np.log(p6b.normalized_z(dat)),
                        "delta": 100 * lev["delta"], "tau": 100 * lev["tau"],
                        "s": 100 * lev["s"], "u": 100 * dat["u"], "v": 100 * dat["v"]})
    return out[(out.index >= p6b.SAMPLE_START) & (out.index <= p6b.SAMPLE_END)]


def plot_levels(lev, out_path):
    """
    One panel per series in levels, with its HP(1,600) trend dashed: what the
    filter removes. Dotted line at the start of BED (1992Q3); NBER recessions
    shaded.
    """
    panels = [("log z", ["log z"], "Labor productivity (log, normalized)"),
              ("tau", ["tau", "s"], "Total separation tau and match separation s (% per month)"),
              ("delta", ["delta"], "Permanent exit delta (% per month, BED)"),
              ("u", ["u"], "Unemployment rate (%)"),
              ("v", ["v"], "Vacancy rate (%)")]
    colors = {"log z": "#1f77b4", "tau": "#7f7f7f", "s": "#2ca02c", "delta": "#d62728",
              "u": "#9467bd", "v": "#ff7f0e"}
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 12), sharex=True)
    for ax, (_, cols, title) in zip(axes, panels):
        for c in cols:
            x = lev[c].dropna()
            dts = list(x.index.to_timestamp())
            ax.plot(dts, x.values, color=colors[c], lw=1.2, label=c)
            _, trend = p6b.hpfilter(x, lamb=p6b.HP_LAMBDA)
            ax.plot(dts, np.asarray(trend), color=colors[c], lw=1.0, ls="--",
                    label=f"{c}, HP trend")
        for lo, hi in p6b.REC_SPANS:
            ax.axvspan(pd.Timestamp(lo), pd.Timestamp(hi), alpha=0.12, color="grey")
        ax.axvline(pd.Period(BED_START, freq="Q").to_timestamp(), color="black", ls=":", lw=1.0)
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", lw=0.4, alpha=0.4)
        ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
    fig.suptitle(f"Series in levels, {p6b.SAMPLE_START}-{p6b.SAMPLE_END}, with HP({p6b.HP_LAMBDA:,}) "
                 f"trends; dotted line: BED start {BED_START}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p6b._finish(fig, out_path)


def build_cycles(method="hp"):
    """z, delta, s from part6b; u and v filtered the same way ("hp" or "hamilton")."""
    dat = pd.read_pickle(p6b.RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    cyc = p6b.hp_cycles(dat, p6b.build_levels(dat), method)[["z", "delta", "s"]]
    for name in OUTCOMES:
        cyc[name] = p6b.log_cycle(p6b.window(dat[name]), method)
    fin = pd.read_parquet(NFCI_CACHE)
    fin.index = pd.PeriodIndex(fin["quarter_label"], freq="Q")
    for name in ("nfci", "anfci"):                   # mean-zero indices: filter levels
        x = p6b.window(fin[name])
        cyc[name] = (p6b.hpfilter(x, lamb=p6b.HP_LAMBDA)[0] if method == "hp"
                     else p6b.hamilton_cycle(x))
    return cyc[(cyc.index >= p6b.SAMPLE_START) & (cyc.index <= p6b.SAMPLE_END)]


def estimate_all(cyc, report=True):
    """Every SVAR on one set of cycles: fits at each lag and the FEVD rows."""
    svars, rows = {}, []
    for key, spec in SVARS.items():
        sv = prepare_svar(key, spec, cyc)
        sv["fits"] = estimate_svar(sv)
        if report:
            report_lags(sv)
        rows += [r for p, fit in sv["fits"].items() for r in fevd_records(sv, p, fit)]
        svars[key] = sv
    return svars, rows


def irf_frame(svA, filt):
    """SVAR A at its BIC lag, preferred ordering: point IRFs and bootstrap bands."""
    Y, names, p_bic = svA["Y"], svA["names"], svA["p_bic"]
    fitA = svA["fits"][p_bic]
    order = list(SVARS["A"]["vars"])
    point = irf(fitA, order, names, H_IRF)
    draws = bootstrap_irf(Y, fitA, order, names, H_IRF, N_BOOT)
    q = {lvl: np.percentile(draws, lvl, axis=0) for lvl in (5, 16, 84, 95)}
    rows = []
    for shk in SHOCKS:
        j = names.index(shk)
        for resp in list(OUTCOMES) + [shk]:          # u, v, and the shocked series itself
            i = names.index(resp)
            for h in range(H_IRF + 1):
                rows.append({"filter": filt, "svar": "A", "p": p_bic, "shock": shk, "response": resp,
                             "h": h, "irf": 100 * point[h, i, j],
                             "q05": 100 * q[5][h, i, j], "q16": 100 * q[16][h, i, j],
                             "q84": 100 * q[84][h, i, j], "q95": 100 * q[95][h, i, j]})
    return pd.DataFrame(rows)


# =============================================================================
#  VAR
# =============================================================================

def lag_matrix(Y, p):
    """Regressors [1, y_{t-1}, ..., y_{t-p}] and the matching left-hand side."""
    T = Y.shape[0]
    X = np.hstack([np.ones((T - p, 1))] + [Y[p - j:T - j] for j in range(1, p + 1)])
    return X, Y[p:]


def fit_var(Y, p):
    """OLS VAR(p) with a constant. A[j] multiplies y_{t-1-j}."""
    X, Yt = lag_matrix(Y, p)
    B = np.linalg.lstsq(X, Yt, rcond=None)[0]
    U = Yt - X @ B
    k = Y.shape[1]
    Sigma = U.T @ U / (U.shape[0] - X.shape[1])
    A = [B[1 + j * k:1 + (j + 1) * k].T for j in range(p)]
    return {"c": B[0], "A": A, "U": U, "Sigma": Sigma, "p": p, "T": U.shape[0]}


def bic_lag(Y):
    """BIC over p = 1..P_MAX on the common sample of the VAR(P_MAX)."""
    k, out = Y.shape[1], {}
    for p in range(1, P_MAX + 1):
        X, Yt = lag_matrix(Y[P_MAX - p:], p)
        U = Yt - X @ np.linalg.lstsq(X, Yt, rcond=None)[0]
        T = U.shape[0]
        out[p] = np.log(np.linalg.det(U.T @ U / T)) + (k * k * p + k) * np.log(T) / T
    return min(out, key=out.get), out


def max_root(A):
    k, p = A[0].shape[0], len(A)
    comp = np.zeros((k * p, k * p))
    comp[:k] = np.hstack(A)
    comp[k:, :-k] = np.eye(k * (p - 1))
    return float(np.max(np.abs(np.linalg.eigvals(comp))))


def ma_coefficients(A, H):
    """Psi_0..Psi_H of the reduced-form moving average."""
    k = A[0].shape[0]
    Psi = [np.eye(k)]
    for h in range(1, H + 1):
        Psi.append(sum(A[j] @ Psi[h - 1 - j] for j in range(min(h, len(A)))))
    return Psi


def structural(fit, order, names):
    """Impact matrix P (columns = shocks in `order`), mapped back to `names`."""
    idx = [names.index(v) for v in order]
    P_ord = np.linalg.cholesky(fit["Sigma"][np.ix_(idx, idx)])
    P = np.zeros_like(P_ord)
    P[np.ix_(idx, range(len(idx)))] = P_ord        # rows back in the VAR's order
    return P


def fevd(fit, order, names, horizons):
    """Share of the h-step forecast-error variance of each variable by shock."""
    Psi = ma_coefficients(fit["A"], max(horizons))
    P = structural(fit, order, names)
    Theta2 = np.array([(M @ P) ** 2 for M in Psi])            # (H+1, var, shock)
    out = {}
    for h in horizons:
        contrib = Theta2[:h].sum(axis=0)
        out[h] = contrib / contrib.sum(axis=1, keepdims=True)  # rows sum to 1
    return out


def irf(fit, order, names, H):
    P = structural(fit, order, names)
    return np.array([M @ P for M in ma_coefficients(fit["A"], H)])  # (H+1, var, shock)


def mediated(Theta, j, m):
    """
    Counterfactual responses to shock j with variable m held at zero, by
    offsetting shocks to m (Theta: horizons x variables x shocks, columns in
    the ordering; m is both a variable and a shock index). Returns
    (horizons x variables).
    """
    H = Theta.shape[0] - 1
    full = Theta[:, :, j]
    eps = np.zeros(H + 1)
    cf = full.copy()
    for h in range(H + 1):
        eps[h] = -cf[h, m] / Theta[0, m, m]
        cf[h:] += eps[h] * Theta[:H + 1 - h, :, m]
    return cf


def mediation_rows(svar, Theta, draws, names, shock, mediator):
    """Peak and cumulative u, v responses: full, counterfactual, share through mediator."""
    j, m = names.index(shock), names.index(mediator)
    out = []
    for resp in OUTCOMES:
        i = names.index(resp)
        full, cf = Theta[:, i, j], mediated(Theta, j, m)[:, i]
        share = lambda f, c: 100 * (1 - c.sum() / f.sum())
        bs = np.array([share(D[:, i, j], mediated(D, j, m)[:, i]) for D in draws])
        k = int(np.abs(full).argmax())
        out.append([svar, LABEL[shock], LABEL[mediator], resp,
                    f"{100 * full[k]:+.2f} (h {k})", f"{100 * cf[k]:+.2f}",
                    f"{100 * full.sum():+.2f}", f"{100 * cf.sum():+.2f}",
                    f"{share(full, cf):.0f} [{np.percentile(bs, 16):.0f}, {np.percentile(bs, 84):.0f}]"])
    return out


def bootstrap_irf(Y, fit, order, names, H, n_boot, seed=SEED):
    """
    Residual bootstrap, recursive design, initial values fixed. The generator is
    seeded here, not passed in, so the bands do not depend on what else ran in
    the session (a Spyder selection or a re-run would otherwise reuse an
    advanced generator).
    """
    rng = np.random.default_rng(seed)
    p, c, A = fit["p"], fit["c"], fit["A"]
    U = fit["U"] - fit["U"].mean(axis=0)
    T, k = U.shape
    draws = np.empty((n_boot, H + 1, k, k))
    for b in range(n_boot):
        e = U[rng.integers(0, T, T)]
        Ys = np.empty_like(Y)
        Ys[:p] = Y[:p]
        for t in range(p, p + T):
            Ys[t] = c + sum(A[j] @ Ys[t - 1 - j] for j in range(p)) + e[t - p]
        draws[b] = irf(fit_var(Ys, p), order, names, H)
    return draws


# =============================================================================
#  ONE SVAR: DATA, FITS, VARIANCE SHARES
# =============================================================================

def prepare_svar(key, spec, cyc):
    """The SVAR's window and data matrix, its shock block and its lags. With
    "align", the window starts where that SVAR's data start (C uses A's
    window, which the Hamilton filter moves to 1995Q2)."""
    names = list(spec["vars"])
    start = spec["start"]
    if "align" in spec:
        start = cyc[list(SVARS[spec["align"]]["vars"])].dropna().index[0]
    #Initialize according to start date
    frame = cyc.loc[cyc.index >= start, names].dropna()
    
    Y = frame.to_numpy()
    p_bic, bic = bic_lag(Y)
    return {"key": key, "names": names, "frame": frame, "Y": Y,
            "block": [v for v in names if v not in OUTCOMES],
            "fixed": list(spec.get("fixed", ("z",))),
            "p_bic": p_bic, "bic": bic, "lags": sorted({1, 2, p_bic})}


def estimate_svar(sv):
    """VAR fits at each reported lag, keyed by p."""
    return {p: fit_var(sv["Y"], p) for p in sv["lags"]}


def block_orderings(block, fixed=("z",)):
    """The fixed head (z, and F in A_F) first, the other shocks in every order."""
    fixed = list(fixed)
    assert block[:len(fixed)] == fixed, (block, fixed)
    return [fixed + list(perm) for perm in itertools.permutations(block[len(fixed):])]


def fevd_records(sv, p, fit):
    """
    One row per (ordering of the shock block, horizon, response, shock). z is
    first in every ordering; u and v always follow the block. "baseline" marks
    the preferred ordering, the SVAR's own variable order. In fevd's output,
    column j belongs to shock order[j].
    """
    rows = []
    for perm in block_orderings(sv["block"], sv["fixed"]):
        order = list(perm) + list(OUTCOMES)
        for h, M in fevd(fit, order, sv["names"], HORIZONS).items():
            rows += [{"svar": sv["key"], "p": p, "p_bic": p == sv["p_bic"],
                      "ordering": "-".join(order), "baseline": perm == sv["block"],
                      "horizon": h, "response": resp, "shock": shk, "share": float(M[i, j])}
                     for i, resp in enumerate(sv["names"]) for j, shk in enumerate(order)]
    return rows


# =============================================================================
#  REPORT
# =============================================================================

def report_method():
    method = __doc__.split("METHOD\n------\n")[1].split("OUTPUT\n------")[0]
    say(f"\n{'-'*78}\n0. Method\n{'-'*78}", echo=ECHO_METHOD)
    for line in method.rstrip().splitlines():
        say(line, echo=ECHO_METHOD)
    if not ECHO_METHOD:
        print("  (method text: see the top of data/results/q3_data_svar.txt)")


def report_lags(sv):
    frame = sv["frame"]
    say(f"  SVAR {sv['key']} {tuple(sv['names'])}  {frame.index[0]}-{frame.index[-1]}, "
        f"T = {len(frame)}; BIC picks p = {sv['p_bic']}  (BIC: "
        + ", ".join(f"{q}: {v:.3f}" for q, v in sv["bic"].items()) + ")")
    for p, fit in sv["fits"].items():
        say(f"      p = {p}: T = {fit['T']}, max |root| = {max_root(fit['A']):.3f}")


def headline_frames(rows):
    """
    Headline variance decomposition, percent: one column group per SVAR, rows
    (variable, horizon) for var(u) and var(v) at HEADLINE_H (h = 4 and the
    unconditional h = 200), each SVAR at its BIC lag.
      head  preferred ordering; "Residual" = u-own + v-own. The order of u and
            v between themselves is arbitrary, so only their sum is identified
      alt   the alternative z-first ordering, for SVARs that have one (A).
            Columns in the same order as head; "Residual" is the same, since
            the block's innovations span the same space in any order
    """
    df = pd.DataFrame(rows)
    df = df[df.response.isin(OUTCOMES) & df.p_bic & df.horizon.isin(HEADLINE_H)]
    df = df.assign(share=lambda d: 100 * d.share)
    idx = ["response", "horizon"]
    head, alt = {}, {}
    for key, spec in SVARS.items():
        d = df[df.svar == key]
        shocks = [v for v in spec["vars"] if v in SHOCKS]

        def table(sel):
            b = sel.pivot_table(index=idx, columns="shock", values="share")
            b["Residual"] = b[list(OUTCOMES)].sum(axis=1)
            return b[shocks + ["Residual"]]

        head[key] = table(d[d.baseline])
        others = sorted(set(d.ordering[~d.baseline]))
        if others:
            alt[key] = table(d[d.ordering == others[0]])
    names = {"response": "var", "horizon": "h"}
    return (pd.concat(head, axis=1).rename_axis(index=names, columns=["SVAR", "shock"]),
            pd.concat(alt, axis=1).rename_axis(index=names, columns=["SVAR", "shock"]))


LABEL   = {"z": "z", "delta": "δ", "s": "s", "u": "u", "v": "v", "nfci": "F", "anfci": "Fa"}
H_LABEL = {4: "4q", 8: "8q", 200: "LR"}


def plain_table(frame):
    """A (variable, horizon) x column frame as a narrow monospace table."""
    cols = list(frame.columns)
    rows = [[LABEL[v], H_LABEL[h]] + [f"{x:.1f}" for x in frame.loc[(v, h)]]
            for v, h in frame.index]
    return tabulate(rows, ["", "h"] + [LABEL.get(c, c) for c in cols], tablefmt="simple",
                    disable_numparse=True, colalign=("left", "left") + ("right",) * len(cols))


def headline_text(head, alt, svars):
    """
    The headline as plain monospace tables, one SVAR per block and narrow
    enough for a phone (about 40 characters). Wrapped in ``` so that WhatsApp
    and Slack show it in a fixed-width font.
    """
    blocks = ["Variance shares of u and v, % (data)"]
    for key, sv in svars.items():
        frame = sv["frame"]
        names = ", ".join(LABEL[v] for v in sv["names"])
        block = (f"SVAR {key}: {names}\n{frame.index[0]}-{frame.index[-1]}, p = {sv['p_bic']}\n"
                 + plain_table(head[key]))
        if key in alt.columns.get_level_values(0):
            shocks = [c for c in alt[key].columns if c != "Residual"]
            order = ", ".join(LABEL[v] for v in block_orderings(sv["block"])[1])
            block += f"\n\nAlternative ordering {order}\n" + plain_table(alt[key][shocks])
        blocks.append(block)
    blocks.append("4q = 4 quarters; LR = long run (h = 200).\n"
                  "z first. Preferred: z, s, δ (s is\n"
                  "exogenous in the model; exit responds\n"
                  "to the cycle). Residual = u and v\n"
                  "innovations not explained by z, δ, s;\n"
                  "the same in both orderings.")
    return "```\n" + "\n\n".join(blocks) + "\n```"


def takeaways(head, alt, irf_df, svars, rows):
    """Key takeaways, each with its numbers filled in from this run."""
    LR = HORIZONS[-1]
    A, B, C, Aa = head["A"], head["B"], head["C"], alt["A"]
    g = lambda f, var, shk, h=LR: f.loc[(var, h), shk]
    d = pd.DataFrame(rows)
    d = d[d.baseline & d.p_bic & d.response.isin(OUTCOMES)].pivot_table(
        index=["svar", "response", "shock"], columns="horizon", values="share")
    gap = 100 * (d[8] - d[LR]).abs().max()
    root = max(max_root(sv["fits"][sv["p_bic"]]["A"]) for sv in svars.values())

    def peak(shk, resp, sign=1):
        r = irf_df[(irf_df.shock == shk) & (irf_df.response == resp)]
        x = r.loc[(sign * r.irf).idxmax()]
        return int(x.h), float(x.irf)

    cmp = lambda a, b: "more than" if a > b + 0.05 else ("less than" if a < b - 0.05 else "about as much as")
    hz, hd, hs = peak("z", "u", -1)[0], peak("delta", "u")[0], peak("s", "u")[0]
    (hvd, vd), (hvs, vs) = peak("delta", "v", -1), peak("s", "v", -1)
    items = [
        f"z explains {g(A,'u','z'):.0f} % of var(u) and {g(A,'v','z'):.0f} % of var(v)\n"
        f"   (A, LR): the largest share, with z ordered first.",
        f"Preferred ordering (z, s, delta): delta explains\n"
        f"   {cmp(g(A,'u','delta'), g(A,'u','s'))} s of var(u) ({g(A,'u','delta'):.1f} vs {g(A,'u','s'):.1f} %) and\n"
        f"   {cmp(g(A,'v','delta'), g(A,'v','s'))} s of var(v) ({g(A,'v','delta'):.1f} vs {g(A,'v','s'):.1f} %).\n"
        f"   Alternative (z, delta, s): {g(Aa,'u','delta'):.1f} vs {g(Aa,'u','s'):.1f} % for u,\n"
        f"   {g(Aa,'v','delta'):.1f} vs {g(Aa,'v','s'):.1f} % for v. The ordering moves only the split:\n"
        f"   delta + s is {g(A,'u','delta') + g(A,'u','s'):.1f} % of var(u) in both.",
        f"The three shocks explain at most {100 - g(A,'u','Residual'):.0f} % of var(u) and\n"
        f"   {100 - g(A,'v','Residual'):.0f} % of var(v) (A, LR). An upper bound: they are ordered\n"
        f"   before u and v. A three-shock model should match the\n"
        f"   composition of this share, not total volatility.",
        f"s matters more in the short run: {g(A,'u','s',4):.1f} % of var(u) at 4q\n"
        f"   against {g(A,'u','s'):.1f} % in the LR (A). Shares at 8q and LR differ by at\n"
        f"   most {gap:.1f} points: HP cycles are short-lived (largest VAR\n"
        f"   root {root:.2f}), so LR is a business-cycle share.",
        f"s acts fastest: the u response peaks at h = {hs} after s,\n"
        f"   h = {hz} after z and h = {hd} after delta (A, preferred).",
        f"Leaving delta out inflates s: on the same window, s's share\n"
        f"   of var(u) is {g(A,'u','s'):.1f} % in A and {g(C,'u','s'):.1f} % in C. B's s share\n"
        f"   ({g(B,'u','s'):.1f} %) is an upper bound.",
        f"z is not weaker after 1992: its share of var(u) is\n"
        f"   {g(B,'u','z'):.1f} % in B (1951-2019) and {g(C,'u','z'):.1f} % in C (1992-2019).",
        f"v falls after delta and after s shocks: troughs {vd:+.2f} at\n"
        f"   h = {hvd} and {vs:+.2f} at h = {hvs} (log points x 100). delta destroys\n"
        f"   vacancies; v falling after s argues against costless reposting.",
    ]
    return "\n".join(f"{i}. {x}" for i, x in enumerate(items, 1))


def share_frames(rows):
    """
    Two flat tables of the shares of var(u) and var(v), in percent, indexed by
    (SVAR, lags, variable, horizon):
      base    baseline ordering; one column per shock, u and v being the
              residual innovations
      bounds  "min-max" over orderings of the shock block, one column per shock
    delta is blank for B and C, which do not contain it.
    """
    df = pd.DataFrame(rows)
    df = df[df.response.isin(OUTCOMES)].assign(share=lambda d: 100 * d.share)
    idx = ["svar", "p", "response", "horizon"]
    cols = list(SHOCKS) + list(OUTCOMES)
    base = (df[df.baseline].pivot_table(index=idx, columns="shock", values="share")
            .reindex(columns=cols))
    rng_ = (df[df.shock.isin(SHOCKS)].groupby(idx + ["shock"])["share"]
            .agg(["min", "max"]).unstack("shock"))
    bounds = pd.DataFrame({shk: [f"{lo:4.1f}-{hi:4.1f}" if pd.notna(lo) else "-"
                                 for lo, hi in zip(rng_[("min", shk)], rng_[("max", shk)])]
                           for shk in SHOCKS}, index=rng_.index)
    names = {"svar": "SVAR", "p": "lags", "response": "var", "horizon": "h"}
    return base.rename_axis(index=names, columns=None), bounds.rename_axis(index=names)


# =============================================================================
#  RUN
# =============================================================================

if __name__ == "__main__":
    REPORT_LINES.clear()
    say("=" * 78)
    say("Q3 data SVAR — Block B' data side  (q3_data_svar.py)")
    say("=" * 78)
    report_method()

    cyc = build_cycles()
    print((cyc.std()*100).map(lambda  x: f"{x:6.2f}"))
    say(f"\n{'-'*78}\n1. Data\n{'-'*78}")
    for name in cyc.columns:
        s_ = cyc[name].dropna()
        say(f"  {name:6s} {s_.index[0]}-{s_.index[-1]}  T = {len(s_):3d}  SD of cycle {100*s_.std():6.2f}")

    # ── estimation: every SVAR, BIC lag plus p = 1 and 2, every ordering ─────
    say(f"\n{'-'*78}\n2. Lag selection and stability\n{'-'*78}")
    svars, fevd_rows = estimate_all(cyc)

    # ── variance shares ──────────────────────────────────────────────────────
    say(f"\n{'-'*78}\n3. Shares of the variance of u and v, percent\n{'-'*78}")
    say(f"  Headline: h = 4 and h = {HORIZONS[-1]}, the unconditional share (the model-side object);")
    say("  the detail tables add h = 8.")
    say("  BIC picks the lag in section 2.")
    head, alt = headline_frames(fevd_rows)
    headline = headline_text(head, alt, svars)
    say("\n  Headline: preferred ordering, BIC lag; each row sums to 100")
    say(headline, echo=False)
    shares, bounds = share_frames(fevd_rows)
    say("\n  Detail 3a. Preferred ordering, p = 1 and 2, u-own and v-own separately; rows sum to 100")
    say(shares.to_string(float_format=lambda x: f"{x:5.1f}", na_rep="    -"))
    say("\n  Detail 3b. Min-max of each shock's share over the z-first orderings (A: two; B, C: one)")
    say(bounds.to_string())

    # ── A against C: what omitting delta does ────────────────────────────────
    say(f"\n{'-'*78}\n4. Omitting delta: A against C on {BED_START}-{p6b.SAMPLE_END}, "
        f"unconditional, preferred ordering\n{'-'*78}")
    df = pd.DataFrame(fevd_rows)
    for p in (1, 2):
        say(f"  p = {p}")
        for resp in OUTCOMES:
            def sh(svar, shk):
                r = df[(df.svar == svar) & (df.p == p) & (df.horizon == HORIZONS[-1]) &
                       df.baseline & (df.response == resp) & (df.shock == shk)]
                return 100 * r["share"].iloc[0]
            a_tot = sh("A", "z") + sh("A", "delta") + sh("A", "s")
            c_tot = sh("C", "z") + sh("C", "s")
            say(f"    var({resp}):  z  A {sh('A','z'):5.1f}  C {sh('C','z'):5.1f}   "
                f"s  A {sh('A','s'):5.1f}  C {sh('C','s'):5.1f}   delta A {sh('A','delta'):5.1f}   "
                f"shock block  A {a_tot:5.1f}  C {c_tot:5.1f}")
    say("  C's z and s shares minus A's show how much of delta's share the other shocks absorb"
        "\n  when delta is left out, as it must be in B.")

    # ── financial conditions: A with NFCI or ANFCI after z ──────────────────
    from tabulate import tabulate
    say(f"\n{'-'*78}\n4b. Financial conditions: A against A with F after z, BIC lag, preferred ordering\n{'-'*78}")
    fin_rows, svars_f = [], {}
    for key, spec in SVARS_F.items():
        sv = prepare_svar(key, spec, cyc)
        sv["fits"] = estimate_svar(sv)
        report_lags(sv)
        fin_rows += [r for p, fit in sv["fits"].items() for r in fevd_records(sv, p, fit)]
        svars_f[key] = sv
    allf = pd.DataFrame(fevd_rows + fin_rows)
    allf = allf[allf.p_bic & allf.response.isin(OUTCOMES)]

    def fshare(svar, resp, shk, h, base=True):
        r = allf[(allf.svar == svar) & (allf.response == resp) & (allf.shock == shk)
                 & (allf.horizon == h) & (allf.baseline == base)]
        return 100 * r["share"].iloc[0] if len(r) else np.nan

    rows = []
    for resp in OUTCOMES:
        for h in HEADLINE_H:
            for key, fcol in (("A", None), ("A_F", "nfci"), ("A_Fa", "anfci")):
                shares = [fshare(key, resp, s_, h) for s_ in ("z",)] + \
                         [fshare(key, resp, fcol, h) if fcol else np.nan] + \
                         [fshare(key, resp, s_, h) for s_ in ("s", "delta")]
                resid = 100 - np.nansum(shares)
                d_alt = fshare(key, resp, "delta", h, base=False)
                rows.append([resp, H_LABEL[h], key] + [("-" if np.isnan(x) else f"{x:.1f}") for x in shares]
                            + [f"{resid:.1f}", f"{d_alt:.1f}"])
    say("\n  Shares of var(u) and var(v), percent; F = NFCI (A_F) or ANFCI (A_Fa); last column:")
    say("  delta's share with delta ordered before s (the alternative ordering)")
    say(tabulate(rows, ["var", "h", "SVAR", "z", "F", "s", "δ", "Residual", "δ (δ before s)"],
                 tablefmt="simple", disable_numparse=True, colalign=("left",) * 3 + ("right",) * 6))
    fin_take = (f"Net of financial conditions, delta's share of var(u) (LR) is\n"
                f"   {fshare('A_F','u','delta',HORIZONS[-1]):.1f} % with NFCI and "
                f"{fshare('A_Fa','u','delta',HORIZONS[-1]):.1f} % with ANFCI, against "
                f"{fshare('A','u','delta',HORIZONS[-1]):.1f} % in A.\n"
                f"   F itself explains {fshare('A_F','u','nfci',HORIZONS[-1]):.1f} % (NFCI) and "
                f"{fshare('A_Fa','u','anfci',HORIZONS[-1]):.1f} % (ANFCI).")

    # δ → u and v, point IRFs at the BIC lag, A against A_F and A_Fa
    rows = []
    for key, sv in [("A", svars["A"])] + list(svars_f.items()):
        th = irf(sv["fits"][sv["p_bic"]], list(sv["names"]), sv["names"], H_IRF)
        j = sv["names"].index("delta")
        cells = [key]
        for resp in OUTCOMES:
            path = 100 * th[:, sv["names"].index(resp), j]
            k = int(np.abs(path).argmax())
            cells.append(f"{path[k]:+.2f} (h {k})")
        rows.append(cells)
    say("\n  delta shock, peak responses (one SD, log points x 100), preferred ordering")
    say(tabulate(rows, ["SVAR", "u", "v"], tablefmt="simple", disable_numparse=True))

    # ── mediation: how much of F's and z's effect runs through delta ───────
    say(f"\n{'-'*78}\n4c. Mediation: share of a shock's effect on u and v that runs through delta\n{'-'*78}")
    say("  delta's response held at zero by offsetting delta shocks; BIC lag, preferred ordering.")
    say("  Peak: response at the full path's peak (one SD, log points x 100); cumulative: sum over")
    say(f"  h = 0-{H_IRF}. Share = 1 - counterfactual/full (cumulative), [68 % bootstrap band].\n")
    med = []
    for key, sv, cases in (("A",    svars["A"],      (("z", "delta"), ("z", "s"))),
                           ("A_F",  svars_f["A_F"],  (("nfci", "delta"), ("nfci", "s"), ("z", "delta"))),
                           ("A_Fa", svars_f["A_Fa"], (("anfci", "delta"), ("anfci", "s")))):
        fit = sv["fits"][sv["p_bic"]]
        Theta = irf(fit, list(sv["names"]), sv["names"], H_IRF)
        draws = bootstrap_irf(sv["Y"], fit, list(sv["names"]), sv["names"], H_IRF, N_BOOT)
        for shock, mediator in cases:
            med += mediation_rows(key, Theta, draws, sv["names"], shock, mediator)
    say(tabulate(med, ["SVAR", "shock", "via", "resp", "peak full", "peak cf", "cum full", "cum cf",
                       "share via, %"], tablefmt="simple", disable_numparse=True))
    say("  'via s' rows: the same exercise with s as the mediator, for comparison.")
    med_df = pd.DataFrame(med, columns=["svar", "shock", "via", "resp", "peak_full", "peak_cf",
                                        "cum_full", "cum_cf", "share_via"])
    pick = lambda sv_, sh, via, r: med_df[(med_df.svar == sv_) & (med_df.shock == sh) &
                                          (med_df.via == via) & (med_df.resp == r)].share_via.iloc[0]
    med_take = (f"Exit as a channel: of the financial shock's cumulative effect on\n"
                f"   u, {pick('A_F', 'F', 'δ', 'u')} % runs through delta (NFCI; ANFCI {pick('A_Fa', 'Fa', 'δ', 'u')} %);\n"
                f"   for v, {pick('A_F', 'F', 'δ', 'v')} %. For z: {pick('A', 'z', 'δ', 'u')} % of its effect on u (A).\n"
                f"   [68 % bootstrap bands in section 4c.]")

    # ── impulse responses, SVAR A at the BIC lag ─────────────────────────────
    svA = svars["A"]
    p_bic = svA["p_bic"]
    irf_df = irf_frame(svA, "hp")

    say(f"\n{'-'*78}\n5. Impulse responses, SVAR A, p = {p_bic} (BIC), preferred ordering; "
        f"one-SD shocks,\n   log points x 100; [68 % band] from {N_BOOT} bootstrap draws\n{'-'*78}")
    for shk in SHOCKS:
        for resp in OUTCOMES:
            r = irf_df[(irf_df.shock == shk) & (irf_df.response == resp)]
            cells = [f"h{int(x.h)} {x.irf:+5.2f} [{x.q16:+5.2f},{x.q84:+5.2f}]"
                     for x in r.itertuples() if x.h in (0, 1, 2, 4, 6, 8, 12)]
            say(f"  {shk:5s} -> {resp}:  " + "  ".join(cells))
    say("\n  Scale: each shock's impact on its own series, and the peak u and v responses per unit")
    say("  of that impact (log points per log point; signs as estimated)")
    for shk in SHOCKS:
        own = irf_df[(irf_df.shock == shk) & (irf_df.response == shk) & (irf_df.h == 0)].irf.iloc[0]
        cells = []
        for resp in OUTCOMES:
            r = irf_df[(irf_df.shock == shk) & (irf_df.response == resp)]
            x = r.loc[r.irf.abs().idxmax()]
            cells.append(f"peak {resp} {x.irf:+5.2f} at h = {int(x.h)}, per unit {x.irf / own:+5.2f}")
        say(f"  {shk:5s} own impact {own:+5.2f}:  " + ";  ".join(cells))

    # ── filter robustness: everything again on Hamilton cycles ──────────────
    cyc_h = build_cycles("hamilton")
    svars_h, fevd_h = estimate_all(cyc_h, report=False)
    head_h, alt_h = headline_frames(fevd_h)
    irf_h = irf_frame(svars_h["A"], "hamilton")
    say(f"\n{'-'*78}\n7. Filter robustness: Hamilton (h = {p6b.HAMILTON_H}, p = {p6b.HAMILTON_P}) "
        f"against HP, BIC lag, preferred ordering\n{'-'*78}")
    for key, sv in svars_h.items():
        fr = sv["frame"]
        say(f"  SVAR {key} (Hamilton): {fr.index[0]}-{fr.index[-1]}, T = {len(fr)}, BIC p = {sv['p_bic']}, "
            f"max |root| = {max_root(sv['fits'][sv['p_bic']]['A']):.3f}")
    from tabulate import tabulate
    rows = []
    for key in SVARS:
        for v, h in head[key].index:
            for c in head[key].columns:
                rows.append([key, LABEL[v], H_LABEL[h], LABEL.get(c, c),
                             f"{head[key].loc[(v, h), c]:.1f}", f"{head_h[key].loc[(v, h), c]:.1f}"])
    say("\n  Shares of var(u) and var(v), percent: HP against Hamilton")
    say(tabulate(rows, ["SVAR", "var", "h", "shock", "HP", "Hamilton"], tablefmt="simple",
                 disable_numparse=True, colalign=("left",) * 4 + ("right",) * 2))
    say("\n  Peak responses of u and v, SVAR A, one-SD shocks, log points x 100 (h of the peak)")
    rows = []
    for shk in SHOCKS:
        for resp in OUTCOMES:
            cells = [LABEL[shk], resp]
            for d in (irf_df, irf_h):
                r = d[(d.shock == shk) & (d.response == resp)]
                x = r.loc[r.irf.abs().idxmax()]
                cells.append(f"{x.irf:+.2f} (h {int(x.h)})")
            rows.append(cells)
    say(tabulate(rows, ["shock", "resp", "HP", "Hamilton"], tablefmt="simple", disable_numparse=True))
    irf_df = pd.concat([irf_df, irf_h], ignore_index=True)

    # ── key takeaways ────────────────────────────────────────────────────────
    key_points = takeaways(head, alt, irf_df[irf_df["filter"] == "hp"], svars, fevd_rows)
    key_points += f"\n9. {fin_take}\n10. {med_take}"
    say(f"\n{'-'*78}\n6. Key takeaways (numbers from this run)\n{'-'*78}")
    say(key_points, echo=False)

    # ── outputs ──────────────────────────────────────────────────────────────
    pd.DataFrame(fevd_rows + fin_rows).to_csv(OUT_FEVD, index=False)
    head.round(2).to_csv(OUT_HEAD)          # full precision in OUT_FEVD
    OUT_HEAD_TXT.write_text(headline + "\n", encoding="utf-8")
    irf_df.to_csv(OUT_IRF, index=False)
    OUT_TXT.write_text("\n".join(REPORT_LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_HEAD}\nSaved: {OUT_HEAD_TXT}\nSaved: {OUT_FEVD}\nSaved: {OUT_IRF}")

    plot_levels(build_levels_all(), OUT_LEV)

    fig, axes = plt.subplots(len(SHOCKS), len(OUTCOMES), figsize=(10, 9), sharex=True)
    for r_, shk in enumerate(SHOCKS):
        for c_, resp in enumerate(OUTCOMES):
            ax = axes[r_, c_]
            d = irf_df[(irf_df["filter"] == "hp") & (irf_df.shock == shk) & (irf_df.response == resp)]
            dh = irf_df[(irf_df["filter"] == "hamilton") & (irf_df.shock == shk) & (irf_df.response == resp)]
            ax.plot(dh.h, dh.irf, color="#d62728", lw=1.4, ls="--", label="Hamilton")
            ax.fill_between(d.h, d.q05, d.q95, color="#1f77b4", alpha=0.15, lw=0)
            ax.fill_between(d.h, d.q16, d.q84, color="#1f77b4", alpha=0.30, lw=0)
            ax.plot(d.h, d.irf, color="#1f77b4", lw=1.8)
            ax.axhline(0, color="black", lw=0.7)
            ax.set_title(f"{shk} shock -> {resp}", fontsize=9)
            ax.grid(axis="y", lw=0.4, alpha=0.4)
    for ax in axes[-1]:
        ax.set_xlabel("Quarters", fontsize=8)
    fig.suptitle(f"SVAR A ({BED_START}-{p6b.SAMPLE_END}, p = {p_bic}, order z, s, delta): one-SD shocks, "
                 f"68 % and 90 % bootstrap bands", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p6b._finish(fig, OUT_PNG)
    # ── console summary: the headline, last thing on screen ─────────────────
    print("\nHEADLINE (copy from ``` to ```; also in " + str(OUT_HEAD_TXT) + ")\n")
    print(headline)
    print("\nKEY TAKEAWAYS\n" + key_points)
    print("\nIn the namespace: head, alt (DataFrames).")
    print("\nDone.")
