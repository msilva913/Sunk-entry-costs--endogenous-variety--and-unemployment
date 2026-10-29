"""
q3_data_svar.py
===============
Data side of Block B' (Q3, decision D13): recursive SVARs on the model's three
shocks (z, delta, s) and the labor-market outcomes (u, v). Reports how much of
the variance of u and v each shock explains, with bounds across orderings, and
the impulse responses of u and v to the delta and s shocks with bootstrap
bands. The text between METHOD and OUTPUT is printed at the top of the report.

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
    A  (z, delta, s, u, v)  1992Q3-2019Q4  delta's share; the full B' target
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

Estimation. Reduced-form VAR(p) with a constant, equation by equation OLS on a
common sample. p is chosen by BIC over 1-4 (common sample of the VAR(4)); p = 1
and p = 2 are always reported as well.

Identification. Recursive (Cholesky), the shock block (z, delta, s) ordered
before (u, v), as in the shock calibration. The shares of the shock-block
innovations do not depend on how u and v are ordered among themselves. They do
depend on the ordering inside the shock block, so the report gives the share
of each shock at its minimum and maximum over all orderings of that block
(6 for A, 2 for B and C). Baseline ordering: z, delta, s, u, v.

Variance shares. Forecast-error variance decomposition at h = 4, 8 and 20
quarters, and at h = 200 as the unconditional share (the model-side object:
share of the variance of HP-filtered u by shock). The residual u and v
innovations collect everything the three shocks do not.

Impulse responses. Baseline ordering, one-SD structural shocks, responses of u
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
    data/results/q3_data_svar_fevd.csv   shares by SVAR, lag, ordering, horizon
    data/results/q3_data_svar_irf.csv    baseline IRFs with bootstrap bands
    data/results/q3_data_svar_irf.png    IRFs of u and v to the z, delta, s shocks

Run
---
    python q3_data_svar.py      (environment: Data/environment.yml)
"""

import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd

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
HORIZONS   = (4, 8, 20, 200)              # 200 = unconditional
H_IRF      = 12
P_MAX      = 4
N_BOOT     = 1_000
SEED       = 20261029
SHOCKS     = ("z", "delta", "s")
OUTCOMES   = ("u", "v")
SVARS = {
    "A": {"vars": ("z", "delta", "s", "u", "v"), "start": BED_START},
    "B": {"vars": ("z", "s", "u", "v"),          "start": p6b.SAMPLE_START},
    "C": {"vars": ("z", "s", "u", "v"),          "start": BED_START},
}
OUT_TXT  = p6b.RESULTS_DIR / "q3_data_svar.txt"
OUT_FEVD = p6b.RESULTS_DIR / "q3_data_svar_fevd.csv"
OUT_IRF  = p6b.RESULTS_DIR / "q3_data_svar_irf.csv"
OUT_PNG  = p6b.RESULTS_DIR / "q3_data_svar_irf.png"

REPORT_LINES = []


def say(text=""):
    print(text)
    REPORT_LINES.append(str(text))


# =============================================================================
#  DATA
# =============================================================================

def build_cycles():
    """z, delta, s from part6b; u and v filtered the same way."""
    dat = pd.read_pickle(p6b.RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    cyc = p6b.hp_cycles(dat, p6b.build_levels(dat))[["z", "delta", "s"]]
    for name in OUTCOMES:
        cyc[name] = p6b.hp_log_cycle(p6b.window(dat[name]))
    return cyc[(cyc.index >= p6b.SAMPLE_START) & (cyc.index <= p6b.SAMPLE_END)]


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


def bootstrap_irf(Y, fit, order, names, H, n_boot, rng):
    """Residual bootstrap, recursive design, initial values fixed."""
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
    """The SVAR's window and data matrix, its shock block and its lags."""
    names = list(spec["vars"])
    #Initialize according to start date
    frame = cyc.loc[cyc.index >= spec["start"], names].dropna()
    
    Y = frame.to_numpy()
    p_bic, bic = bic_lag(Y)
    return {"key": key, "names": names, "frame": frame, "Y": Y,
            "block": [v for v in names if v in SHOCKS],
            "p_bic": p_bic, "bic": bic, "lags": sorted({1, 2, p_bic})}


def estimate_svar(sv):
    """VAR fits at each reported lag, keyed by p."""
    return {p: fit_var(sv["Y"], p) for p in sv["lags"]}


def fevd_records(sv, p, fit):
    """
    One row per (ordering of the shock block, horizon, response, shock). The
    shock block is permuted; u and v always follow it. In fevd's output,
    column j belongs to shock order[j].
    """
    rows = []
    for perm in itertools.permutations(sv["block"]):
        order = list(perm) + list(OUTCOMES)
        for h, M in fevd(fit, order, sv["names"], HORIZONS).items():
            rows += [{"svar": sv["key"], "p": p, "p_bic": p == sv["p_bic"],
                      "ordering": "-".join(order), "baseline": list(perm) == sv["block"],
                      "horizon": h, "response": resp, "shock": shk, "share": float(M[i, j])}
                     for i, resp in enumerate(sv["names"]) for j, shk in enumerate(order)]
    return rows


# =============================================================================
#  REPORT
# =============================================================================

def report_method():
    method = __doc__.split("METHOD\n------\n")[1].split("OUTPUT\n------")[0]
    say(f"\n{'-'*78}\n0. Method\n{'-'*78}")
    for line in method.rstrip().splitlines():
        say(line)


def report_lags(sv):
    frame = sv["frame"]
    say(f"  SVAR {sv['key']} {tuple(sv['names'])}  {frame.index[0]}-{frame.index[-1]}, "
        f"T = {len(frame)}; BIC picks p = {sv['p_bic']}  (BIC: "
        + ", ".join(f"{q}: {v:.3f}" for q, v in sv["bic"].items()) + ")")
    for p, fit in sv["fits"].items():
        say(f"      p = {p}: T = {fit['T']}, max |root| = {max_root(fit['A']):.3f}")


def share_table(rows, svar, p, horizon):
    """Baseline-ordering shares at one horizon, plus min-max over orderings."""
    df = pd.DataFrame(rows)
    df = df[(df.svar == svar) & (df.p == p) & (df.horizon == horizon)]
    base = df[df.baseline]
    lines = []
    for resp in OUTCOMES:
        cells = []
        for shk in SVARS[svar]["vars"]:
            b = base[(base.response == resp) & (base.shock == shk)]["share"].iloc[0]
            if shk in SHOCKS:
                rng_ = df[(df.response == resp) & (df.shock == shk)]["share"]
                cells.append(f"{shk} {100*b:5.1f} [{100*rng_.min():4.1f}-{100*rng_.max():4.1f}]")
            else:
                cells.append(f"{shk}-own {100*b:5.1f}")
        lines.append(f"      var({resp}):  " + "   ".join(cells))
    return lines


# =============================================================================
#  RUN
# =============================================================================

if __name__ == "__main__":
    REPORT_LINES.clear()
    rng = np.random.default_rng(SEED)
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
    svars, fevd_rows = {}, []
    for key, spec in SVARS.items():
        sv = prepare_svar(key, spec, cyc)
        sv["fits"] = estimate_svar(sv)
        report_lags(sv)
        fevd_rows += [r for p, fit in sv["fits"].items() for r in fevd_records(sv, p, fit)]
        svars[key] = sv

    # ── variance shares ──────────────────────────────────────────────────────
    say(f"\n{'-'*78}\n3. Shares of the variance of u and v, percent\n{'-'*78}")
    say("  Baseline ordering (z, delta, s, u, v); [min-max] over orderings of the shock block.")
    say("  h = 200 is the unconditional share, the model-side object. '-own' = residual innovation.")
    for key, sv in svars.items():
        for p in sv["lags"]:
            say(f"\n  SVAR {key}, p = {p}{' (BIC)' if p == sv['p_bic'] else ''}")
            for h in HORIZONS:
                say(f"    h = {h}" + (" (unconditional)" if h == HORIZONS[-1] else ""))
                for line in share_table(fevd_rows, key, p, h):
                    say(line)

    # ── A against C: what omitting delta does ────────────────────────────────
    say(f"\n{'-'*78}\n4. Omitting delta: A against C on {BED_START}-{p6b.SAMPLE_END}, "
        f"unconditional, baseline ordering\n{'-'*78}")
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

    # ── impulse responses, SVAR A at the BIC lag ─────────────────────────────
    svA = svars["A"]
    Y, names, p_bic = svA["Y"], svA["names"], svA["p_bic"]
    fitA = svA["fits"][p_bic]
    order = list(SVARS["A"]["vars"])
    point = irf(fitA, order, names, H_IRF)
    draws = bootstrap_irf(Y, fitA, order, names, H_IRF, N_BOOT, rng)
    q = {lvl: np.percentile(draws, lvl, axis=0) for lvl in (5, 16, 84, 95)}
    irf_rows = []
    for shk in SHOCKS:
        j = names.index(shk)
        for resp in OUTCOMES:
            i = names.index(resp)
            for h in range(H_IRF + 1):
                irf_rows.append({"svar": "A", "p": p_bic, "shock": shk, "response": resp,
                                 "h": h, "irf": 100 * point[h, i, j],
                                 "q05": 100 * q[5][h, i, j], "q16": 100 * q[16][h, i, j],
                                 "q84": 100 * q[84][h, i, j], "q95": 100 * q[95][h, i, j]})
    irf_df = pd.DataFrame(irf_rows)

    say(f"\n{'-'*78}\n5. Impulse responses, SVAR A, p = {p_bic} (BIC), baseline ordering; "
        f"one-SD shocks,\n   log points x 100; [68 % band] from {N_BOOT} bootstrap draws\n{'-'*78}")
    for shk in SHOCKS:
        for resp in OUTCOMES:
            r = irf_df[(irf_df.shock == shk) & (irf_df.response == resp)]
            cells = [f"h{int(x.h)} {x.irf:+5.2f} [{x.q16:+5.2f},{x.q84:+5.2f}]"
                     for x in r.itertuples() if x.h in (0, 1, 2, 4, 6, 8, 12)]
            say(f"  {shk:5s} -> {resp}:  " + "  ".join(cells))

    # ── outputs ──────────────────────────────────────────────────────────────
    pd.DataFrame(fevd_rows).to_csv(OUT_FEVD, index=False)
    irf_df.to_csv(OUT_IRF, index=False)
    OUT_TXT.write_text("\n".join(REPORT_LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_FEVD}\nSaved: {OUT_IRF}")

    fig, axes = plt.subplots(len(SHOCKS), len(OUTCOMES), figsize=(10, 9), sharex=True)
    for r_, shk in enumerate(SHOCKS):
        for c_, resp in enumerate(OUTCOMES):
            ax = axes[r_, c_]
            d = irf_df[(irf_df.shock == shk) & (irf_df.response == resp)]
            ax.fill_between(d.h, d.q05, d.q95, color="#1f77b4", alpha=0.15, lw=0)
            ax.fill_between(d.h, d.q16, d.q84, color="#1f77b4", alpha=0.30, lw=0)
            ax.plot(d.h, d.irf, color="#1f77b4", lw=1.8)
            ax.axhline(0, color="black", lw=0.7)
            ax.set_title(f"{shk} shock -> {resp}", fontsize=9)
            ax.grid(axis="y", lw=0.4, alpha=0.4)
    for ax in axes[-1]:
        ax.set_xlabel("Quarters", fontsize=8)
    fig.suptitle(f"SVAR A ({BED_START}-{p6b.SAMPLE_END}, p = {p_bic}): one-SD shocks, "
                 f"68 % and 90 % bootstrap bands", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p6b._finish(fig, OUT_PNG)
    print("\nDone.")
