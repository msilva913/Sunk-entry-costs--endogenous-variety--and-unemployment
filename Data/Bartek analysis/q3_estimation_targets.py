"""
q3_estimation_targets.py -- data targets for the AGS estimation (decisions.md D18).

Auxiliary system: SVAR (z, s, delta, u, v, jf), recursive in that order, on HP and on Hamilton
cycles (q3_data_svar.build_cycles), from BED's start (HP: 1992Q3) or the first complete
Hamilton quarter (1995Q2), lag by BIC (q3_data_svar.bic_lag). Responses to one-SD shocks,
h = 0-12, in log points (not x 100).

Estimation moments: u, v, jf responses to the z and s shocks (78). Test moments: u, v, jf
responses to delta (39). Own responses of z, s, delta (block "own"): validation only. Data SDs from a 1,000-draw residual bootstrap (q3_data_svar.bootstrap_irf).

Outputs:
    data/results/estimation_targets.csv            filter, p, T, shock, response, h, data, sd, block
    data/results/estimation_draws_<filter>.npz     bootstrap draws of the moment vectors, for the
                                                   full covariance Omega (sandwich, D18)
"""

import numpy as np
import pandas as pd

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

NAMES   = ["z", "s", "delta", "u", "v", "jf"]
RESP    = ["u", "v", "jf"]
H       = q3.H_IRF
FILTERS = {"hp": ("hp", q3.BED_START), "ham": ("hamilton", "1995Q2")}
OUT_CSV = p6b.RESULTS_DIR / "estimation_targets.csv"


def moments(Theta, shocks):
    """Moment vector: for each shock, response, horizon (in that nesting order)."""
    return np.array([Theta[h, NAMES.index(r), NAMES.index(s)]
                     for s in shocks for r in RESP for h in range(H + 1)])


if __name__ == "__main__":
    rows = []
    for fkey, (method, start) in FILTERS.items():
        cyc = q3.build_cycles(method)
        frame = cyc.loc[cyc.index >= start, NAMES].dropna()
        Y = frame.to_numpy()
        p, _ = q3.bic_lag(Y)
        fit = q3.fit_var(Y, p)
        Theta = q3.irf(fit, NAMES, NAMES, H)
        draws = q3.bootstrap_irf(Y, fit, NAMES, NAMES, H, q3.N_BOOT)
        sd = draws.std(axis=0)
        for s in ("z", "s", "delta"):
            block = "test" if s == "delta" else "estimation"
            for r in RESP:
                for h in range(H + 1):
                    i, j = NAMES.index(r), NAMES.index(s)
                    rows.append({"filter": fkey, "p": p, "T": len(frame),
                                 "start": str(frame.index[0]), "end": str(frame.index[-1]),
                                 "shock": s, "response": r, "h": h,
                                 "data": Theta[h, i, j], "sd": sd[h, i, j], "block": block})
            # each shock's response of its own series (validation of D20; not estimated or tested)
            i = j = NAMES.index(s)
            for h in range(H + 1):
                rows.append({"filter": fkey, "p": p, "T": len(frame),
                             "start": str(frame.index[0]), "end": str(frame.index[-1]),
                             "shock": s, "response": s, "h": h,
                             "data": Theta[h, i, j], "sd": sd[h, i, j], "block": "own"})
        est = np.array([moments(D, ("z", "s")) for D in draws])
        tst = np.array([moments(D, ("delta",)) for D in draws])
        np.savez(p6b.RESULTS_DIR / f"estimation_draws_{fkey}.npz", estimation=est, test=tst,
                 point_estimation=moments(Theta, ("z", "s")), point_test=moments(Theta, ("delta",)))
        cond = np.linalg.cond(fit["Sigma"])
        print(f"{fkey}: {frame.index[0]}-{frame.index[-1]}, T = {len(frame)}, BIC p = {p}, "
              f"cond(Sigma) = {cond:.0f}, max root {q3.max_root(fit['A']):.3f}")
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"Saved: {OUT_CSV}")
