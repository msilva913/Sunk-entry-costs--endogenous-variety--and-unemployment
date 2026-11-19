"""
part6e_window_processes.py -- monthly shock processes from SVAR A's own equations, by filter.

Rule (amends S14/D14, Nov 19, 2026; decisions.md D20). Each exogenous process in the model is
calibrated to the corresponding shock of SVAR A (z, s, delta, u, v; p = 1, recursive), the system
the model is compared through:
    persistence  the own-lag coefficient of the shock's equation, conditional on the lags of every
                 variable in the system;
    size         the structural shock's impact on its own series (Cholesky diagonal), i.e. the
                 innovation net of the shocks ordered before it.
The univariate AR(1) used before overstates the persistence of a series that responds to the cycle:
with gross flows, s has a univariate rho_Q of 0.73 but an own-lag of 0.15 in SVAR A, because most of
its cyclical variation follows productivity and the labor market (FEVD: own shock 30 % in the long
run). Coles and Kelishomi (2018) also use univariate persistence, with correlated innovations; here
the innovations are orthogonal, as in the SVAR.

Each quarterly target is then converted to a monthly AR(1) with part6b's simulated matching (Step 4b):
the monthly process is simulated, averaged to quarters, filtered over the series' own span as the
data are, and fitted on SVAR A's window, until the quarterly persistence and innovation size match.
    HP        window 1992Q3-2019Q4 (110 quarters); z filtered over 1951-2019, s over 1990-2019,
              delta over 1992-2019.
    Hamilton  window 1995Q2-2019Q4 (99 quarters); same spans, first 11 quarters lost.

Output: data/results/shock_calibration_window.csv, one row per series (z, delta, s): columns
<f>_rho_m, <f>_sigma_m (f = hp, ham; the processes used, read by Programs AGS via shock_processes),
the quarterly targets <f>_rho_Q, <f>_sigma_Q, <f>_attained, and the previous univariate processes
<f>_uni_rho_m, <f>_uni_sigma_m for comparison.
"""

import numpy as np
import pandas as pd

import part6b_shock_persistence_cyclical as p6b
import q3_data_svar as q3

OUT_CSV  = p6b.RESULTS_DIR / "shock_calibration_window.csv"
NAMES    = list(q3.SVARS["A"]["vars"])                    # z, s, delta, u, v
FILTERS  = {"hp": ("hp", q3.BED_START), "ham": ("hamilton", "1995Q2")}


def svar_targets(method, start):
    """Own-lag and own-impact SD of each shock in SVAR A (p = 1) under the given filter."""
    cyc = q3.build_cycles(method)
    frame = cyc.loc[cyc.index >= start, NAMES].dropna()
    fit = q3.fit_var(frame.to_numpy(), 1)
    P = q3.structural(fit, NAMES, NAMES)
    A = fit["A"][0]
    out = {}
    for x in ("z", "delta", "s"):
        i = NAMES.index(x)
        out[x] = {"rho_Q": float(A[i, i]), "sigma_Q": float(P[i, i])}
    return cyc, frame, out


def uni_col(old, fkey, what):
    """The univariate process: the uni column once this program has run, else the original column."""
    if f"{fkey}_uni_{what}" in old.columns:
        return f"{fkey}_uni_{what}"
    else:
        return f"{fkey}_{what}"


if __name__ == "__main__":
    old = pd.read_csv(OUT_CSV).set_index("series")
    rng = np.random.default_rng(p6b.MATCH_SEED)
    rows = {x: {"series": x} for x in ("z", "delta", "s")}
    for fkey, (method, start) in FILTERS.items():
        cyc, frame, tg = svar_targets(method, start)
        n_win = len(frame)
        lost = 0
        if method == "hamilton":
            lost = p6b.HAMILTON_H + p6b.HAMILTON_P - 1
        for x in ("z", "delta", "s"):
            nq = len(cyc[x].dropna()) + lost                # quarters before filtering
            rho_m, sigma_m, eps, C, attained = p6b.match_process(tg[x]["rho_Q"], tg[x]["sigma_Q"], nq, rng,
                                                                 method=method, n_keep=n_win)
            rows[x].update({f"{fkey}_rho_m": rho_m, f"{fkey}_sigma_m": sigma_m, f"{fkey}_attained": attained,
                            f"{fkey}_rho_Q": tg[x]["rho_Q"], f"{fkey}_sigma_Q": tg[x]["sigma_Q"],
                            f"{fkey}_uni_rho_m": old.loc[x, uni_col(old, fkey, "rho_m")],
                            f"{fkey}_uni_sigma_m": old.loc[x, uni_col(old, fkey, "sigma_m")]})
        print(f"{fkey}: SVAR A {frame.index[0]}-{frame.index[-1]}, T = {n_win}")
    out = pd.DataFrame(rows.values()).set_index("series")
    cols = [c for f in FILTERS for c in (f"{f}_rho_Q", f"{f}_sigma_Q", f"{f}_rho_m", f"{f}_sigma_m",
                                         f"{f}_attained", f"{f}_uni_rho_m", f"{f}_uni_sigma_m")]
    print(out[cols].round(4).to_string())
    out.to_csv(OUT_CSV)
    print(f"\nSaved: {OUT_CSV}")
