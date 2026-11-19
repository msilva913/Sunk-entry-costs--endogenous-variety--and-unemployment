"""
q3_s_jf_residuals.py -- is there a shared measurement component in s and jf?

Both series are built from the same CPS short-term unemployment u^s_{t+1} (observables.py,
shimer_flows): jf_t = 1 - (u_{t+1} - u^s_{t+1})/u_t, and the separation probability rises with
u^s_{t+1} given jf. An error in u^s therefore moves measured s and jf in the same direction in the
same month. If that component matters, the reduced-form innovations of s and jf are more
positively correlated in the data than in the model, which has no such measurement.

Reduced-form VAR(1) (z, s, delta, u, v, jf), as the estimation system (D18):
    data   HP (1992Q3-2019Q4) and Hamilton (1995Q2-2019Q4) cycles;
    model  the simulated windows of ags_svar_filters.jl (default calibration, reactivation rates
           1 to 0.03), same filter, same VAR; correlations averaged over 200 windows, with the
           5th-95th percentile across windows (the data correlation's sampling range if the model
           were true).
Reported: the correlations of the innovations of s with jf, u and v, and of jf with u.

Output: data/results/q3_s_jf_residuals.txt
"""

from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

NAMES = ["z", "s", "delta", "u", "v", "jf"]
PAIRS = [("s", "jf"), ("s", "u"), ("s", "v"), ("jf", "u")]
MODEL_DIR = Path("../Programs AGS")
CASES = pd.read_csv(MODEL_DIR / "ags_svar_filters_cases.csv")
FILTERS = {"hp": ("hp", q3.BED_START, "HP"), "ham": ("hamilton", "1995Q2", "Hamilton")}
OUT_TXT = p6b.RESULTS_DIR / "q3_s_jf_residuals.txt"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def resid_corr(Y):
    U = q3.fit_var(Y, 1)["U"]
    return np.corrcoef(U, rowvar=False)


def pick(C, a, b):
    return C[NAMES.index(a), NAMES.index(b)]


if __name__ == "__main__":
    say("=" * 78)
    say("Shared measurement in s and jf? Reduced-form innovation correlations (q3_s_jf_residuals.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])
    for fkey, (method, start, flabel) in FILTERS.items():
        cyc = q3.build_cycles(method)
        frame = cyc.loc[cyc.index >= start, NAMES].dropna()
        Cd = resid_corr(frame.to_numpy())
        tab = [["data", f"{frame.index[0]}-{frame.index[-1]}"] + [f"{pick(Cd, a, b):+.2f}" for a, b in PAIRS]]
        for _, c in CASES[CASES["filter"] == fkey].iterrows():
            reps = pd.read_csv(MODEL_DIR / f"ags_svar_reps_{fkey}_{c.case}.csv")
            vals = np.array([[pick(resid_corr(g[NAMES].dropna().to_numpy()), a, b) for a, b in PAIRS]
                             for _, g in reps.groupby("rep")])
            lo, hi = np.percentile(vals, [5, 95], axis=0)
            tab.append([f"model, Λ̄_r {c.lambda_r:.2f}", ""] +
                       [f"{m:+.2f} [{l:+.2f}, {h:+.2f}]" for m, l, h in zip(vals.mean(axis=0), lo, hi)])
        say(f"\n{'-' * 78}\n{flabel}\n{'-' * 78}")
        say(tabulate(tab, ["", "window"] + [f"corr({a},{b})" for a, b in PAIRS],
                     tablefmt="simple", disable_numparse=True))
    say("\n  Model: mean over 200 windows [5th, 95th percentile across windows].")
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}")
