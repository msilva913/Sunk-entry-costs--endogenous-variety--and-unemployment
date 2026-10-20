"""
part6c_shock_var_trivariate.py
==============================
Q2: shock processes for z, delta and s from one trivariate VAR(1).

part6b fits a bivariate (z, delta) VAR and a separate univariate AR(1) for the
constructed match separation rate s = (tau - delta)/(1 - delta). That leaves
sigma_s as a reduced-form residual SD, never orthogonalized, while the model
treats the three innovations as uncorrelated. This program puts s in the same
system and reports what that changes.

Same data as part6b by construction: the series, the HP(1,600) log cycles,
the window and the quarterly-to-monthly conversion are imported from it.

What it reports
---------------
1. VAR(1) on [z, delta, s], constant, OLS with HC3 standard errors: the full
   companion matrix with t-statistics and its eigenvalues.
2. The reduced-form residual correlation matrix.
3. Cholesky with z -> delta -> s: the share of each innovation's variance
   explained by the innovations ordered before it.
4. Candidate monthly shock processes (rho from the diagonal of A1, as in part6b):
     (a) orthogonalized    sigma = diag of the Cholesky factor
     (b) correlated        sigma = reduced-form residual SD, plus the residual
                           correlation matrix
   against the part6b values that SHOCKS_BASE uses now.
5. Robustness: ordering z -> s -> delta, lag-length criteria, residual
   autocorrelation.

Which of (a) or (b) the model uses is decision D14 (context/decisions.md).
The z -> delta argument of part6b (the model generates that pathway, so the
exogenous process must not carry it) does not transfer to s: the model has no
endogenous s, so orthogonalizing s removes variation the model cannot
reproduce. This program only measures; it decides nothing.

Output
------
    data/results/shock_var_trivariate.txt    report (quoted by §5.2, N15)
    data/results/shock_var_trivariate.csv    candidate processes, machine-readable

Run
---
    python part6c_shock_var_trivariate.py      (environment: Data/environment.yml)
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tools import add_constant

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    if not Path("part6c_shock_var_trivariate.py").exists():
        raise RuntimeError("Run from 'Data/Bartek analysis' or as a whole file.")

import part6b_shock_persistence_cyclical as p6b

VARS     = ["z", "delta", "s"]           # Cholesky order z -> delta -> s
ALT      = ["z", "s", "delta"]           # robustness ordering
OUT_TXT  = p6b.RESULTS_DIR / "shock_var_trivariate.txt"
OUT_CSV  = p6b.RESULTS_DIR / "shock_var_trivariate.csv"

# part6b values, as in SHOCKS_BASE (run_solution_core.jl), for comparison.
PART6B = {"z": (0.902, 0.0092), "delta": (0.592, 0.0669), "s": (0.8741, 0.0854)}


# =============================================================================
#  ESTIMATION
# =============================================================================

def lagged_frame(cyc, names, p):
    """Cycles with p lags each, rows with complete data only."""
    lags = [cyc[names].shift(k).add_suffix(f"_l{k}") for k in range(1, p + 1)]
    return pd.concat([cyc[names]] + lags, axis=1).dropna()


def fit_var(cyc, names, p=1, sample=None):
    """
    VAR(p) by equation-wise OLS with a constant, HC3 SEs (as part6b).
    `sample` fixes the estimation rows, so lag-length criteria compare models
    on the same observations.
    """
    df = lagged_frame(cyc, names, p)
    if sample is not None:
        df = df.loc[sample]
    regs = [c for c in df.columns if c not in names]
    X = add_constant(df[regs].values)
    fits = [OLS(df[c].values, X).fit(cov_type="HC3") for c in names]
    U = np.column_stack([f.resid for f in fits])
    T, k = U.shape[0], X.shape[1]
    Sigma = U.T @ U / (T - k)
    return {"fits": fits, "names": names, "p": p, "U": U, "T": T, "k": k,
            "Sigma": Sigma, "index": df.index,
            "A": [np.array([f.params[1 + j * len(names):1 + (j + 1) * len(names)]
                            for f in fits]) for j in range(p)],
            "tA": [np.array([f.tvalues[1 + j * len(names):1 + (j + 1) * len(names)]
                             for f in fits]) for j in range(p)]}


def info_criteria(cyc, names, pmax=4):
    """AIC/BIC for p = 1..pmax on the common sample of the VAR(pmax)."""
    common = lagged_frame(cyc, names, pmax).index
    rows = []
    for p in range(1, pmax + 1):
        v = fit_var(cyc, names, p, sample=common)
        S_ml = v["U"].T @ v["U"] / v["T"]
        ld = np.log(np.linalg.det(S_ml))
        npar = len(names) * v["k"]
        rows.append({"p": p, "AIC": ld + 2 * npar / v["T"],
                     "BIC": ld + np.log(v["T"]) * npar / v["T"]})
    return pd.DataFrame(rows).set_index("p"), len(common)


def cholesky_summary(Sigma, names):
    """Impact matrix, orthogonalized SDs, and the share of each innovation's
    variance explained by those ordered before it."""
    P = np.linalg.cholesky(Sigma)
    sd_rf = np.sqrt(np.diag(Sigma))
    sd_or = np.diag(P)
    return {"P": P, "sd_rf": dict(zip(names, sd_rf)),
            "sd_or": dict(zip(names, sd_or)),
            "share": dict(zip(names, 1 - sd_or**2 / sd_rf**2))}


def corr(Sigma):
    d = np.sqrt(np.diag(Sigma))
    return Sigma / np.outer(d, d)


# =============================================================================
#  REPORT
# =============================================================================

def main():
    dat = pd.read_pickle(p6b.RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    levels = p6b.build_levels(dat)
    cyc = p6b.hp_cycles(dat, levels)

    v = fit_var(cyc, VARS, p=1)
    A1, tA1, Sigma = v["A"][0], v["tA"][0], v["Sigma"]
    C = corr(Sigma)
    ch = cholesky_summary(Sigma, VARS)
    eig = np.abs(np.linalg.eigvals(A1))
    s_uni = p6b.fit_ar1(cyc["s"])

    # Robustness ordering: same Sigma, s before delta.
    iperm = [VARS.index(n) for n in ALT]
    ch_alt = cholesky_summary(Sigma[np.ix_(iperm, iperm)], ALT)

    ic, n_ic = info_criteria(cyc, VARS, pmax=4)
    lb = {n: acorr_ljungbox(v["U"][:, i], lags=[4, 8], return_df=True)["lb_pvalue"]
          for i, n in enumerate(VARS)}

    # Candidate monthly processes. rho is the VAR diagonal in every case, so
    # the candidates differ only in sigma (and, for (b), the correlations).
    rows = []
    for i, n in enumerate(VARS):
        rho_Q = A1[i, i]
        for spec, sQ in (("a_orthogonalized", ch["sd_or"][n]),
                         ("b_correlated", ch["sd_rf"][n]),
                         ("alt_order_z_s_delta", ch_alt["sd_or"][n])):
            rho_m, sig_m = p6b.to_monthly(rho_Q, sQ)
            rows.append({"spec": spec, "shock": n, "rho_Q": rho_Q, "sigma_Q": sQ,
                         "rho_m": rho_m, "sigma_m": sig_m})
    cand = pd.DataFrame(rows)
    cand.assign(sample_start=str(v["index"][0]), sample_end=str(v["index"][-1]),
                nobs=v["T"], hp_lambda=p6b.HP_LAMBDA).to_csv(OUT_CSV, index=False)

    L = []
    w = L.append
    w("Q2 — trivariate VAR(1) shock processes   (part6c_shock_var_trivariate.py)")
    w(f"HP({p6b.HP_LAMBDA:,}) log cycles; estimation sample {v['index'][0]}–{v['index'][-1]}, "
      f"T = {v['T']} quarters; order z -> delta -> s")
    w("Data, cycles and monthly conversion as in part6b (imported).")
    w("")
    w("1. Companion matrix A1 (rows: equation; columns: lag of z, delta, s), HC3 t-stats")
    w(f"{'':8s}" + "".join(f"{'L.'+n:>18s}" for n in VARS) + f"{'R2':>8s}")
    for i, n in enumerate(VARS):
        w(f"{n:8s}" + "".join(f"{A1[i, j]:+9.3f} ({tA1[i, j]:+5.2f})" for j in range(3))
          + f"{v['fits'][i].rsquared:8.3f}")
    w(f"   |eigenvalues| of A1: " + "  ".join(f"{e:.3f}" for e in sorted(eig, reverse=True)))
    w("")
    w("2. Reduced-form residual correlations")
    w(f"{'':8s}" + "".join(f"{n:>9s}" for n in VARS))
    for i, n in enumerate(VARS):
        w(f"{n:8s}" + "".join(f"{C[i, j]:+9.3f}" for j in range(3)))
    w("")
    w("3. Cholesky z -> delta -> s: quarterly innovation SDs")
    w(f"{'':8s}{'reduced form':>14s}{'orthogonal':>12s}{'share explained by earlier shocks':>36s}")
    for n in VARS:
        w(f"{n:8s}{ch['sd_rf'][n]:14.4f}{ch['sd_or'][n]:12.4f}{100*ch['share'][n]:35.1f}%")
    w(f"   ordering z -> s -> delta: sigma_s {ch_alt['sd_or']['s']:.4f} "
      f"({100*ch_alt['share']['s']:.1f}% explained), sigma_delta {ch_alt['sd_or']['delta']:.4f} "
      f"({100*ch_alt['share']['delta']:.1f}% explained)")
    w("")
    w("4. Candidate monthly processes (rho = VAR diagonal, cubed root; sigma matches the")
    w("   unconditional variance, as in part6b)")
    w(f"{'':8s}{'part6b (now)':>18s}{'(a) orthogonal':>20s}{'(b) correlated':>20s}{'alt order':>20s}")
    w(f"{'':8s}" + f"{'rho_m   sigma_m':>18s}" * 1 + f"{'rho_m   sigma_m':>20s}" * 3)
    for n in VARS:
        r = cand[cand.shock == n].set_index("spec")
        cell = lambda s: f"{r.loc[s, 'rho_m']:.4f}  {r.loc[s, 'sigma_m']:.4f}"
        w(f"{n:8s}{PART6B[n][0]:>10.4f}  {PART6B[n][1]:.4f}" + "".join(
            f"{cell(s):>20s}" for s in ("a_orthogonalized", "b_correlated", "alt_order_z_s_delta")))
    w(f"   (b) also needs the innovation correlations of section 2 (quarterly; carried to")
    w(f"   monthly unchanged, an approximation).")
    w(f"   Univariate AR(1) for s (part6b): rho_Q {s_uni['rho_Q']:.4f}, sigma_Q {s_uni['sigma_Q']:.4f};"
      f" VAR diagonal: rho_Q {A1[2, 2]:.4f}.")
    w("")
    w(f"5. Diagnostics")
    w(f"   Lag length (common sample of VAR(4), T = {n_ic}):")
    for p, r in ic.iterrows():
        w(f"     p = {p}:  AIC {r.AIC:9.3f}   BIC {r.BIC:9.3f}")
    w(f"     minimum: AIC p = {ic.AIC.idxmin()}, BIC p = {ic.BIC.idxmin()}")
    w("   Ljung-Box p-values on VAR(1) residuals (lags 4, 8):")
    for n in VARS:
        w(f"     {n:6s} {lb[n].iloc[0]:.3f}  {lb[n].iloc[1]:.3f}")
    txt = "\n".join(L)
    OUT_TXT.write_text(txt + "\n", encoding="utf-8")
    print(txt)
    print(f"\nSaved: {OUT_TXT}\n       {OUT_CSV}")


if __name__ == "__main__":
    main()
