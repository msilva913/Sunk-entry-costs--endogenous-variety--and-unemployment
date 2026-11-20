"""
q3_stylized_moments.py -- business-cycle moments of the labor-market series, for the Part 1 draft.

HP(1,600) log cycles from q3_data_svar.py (q3_data_svar_cycles_hp.csv): every series cut at 2019Q4
before filtering and filtered over its own window (u, v, theta, jf, z, s: 1951Q1-2019Q4; delta:
1992Q3-2019Q4). Reported for each series: SD of the cycle (x 100), first-order autocorrelation, and
the correlations with u and v over the common sample of the pair; and the lead-lag correlations
cor(x_t, u_{t+k}) for x in {z, s, delta}, k = -4..8.

Output: data/results/q3_stylized_moments.txt
"""

import pandas as pd
from tabulate import tabulate

import part6b_shock_persistence_cyclical as p6b

SERIES = ("u", "v", "theta", "jf", "z", "s", "delta")
OUT = p6b.RESULTS_DIR / "q3_stylized_moments.txt"

if __name__ == "__main__":
    cyc = pd.read_csv(p6b.RESULTS_DIR / "q3_data_svar_cycles_hp.csv", index_col="quarter")
    rows = []
    for x in SERIES:
        s = cyc[x].dropna()
        rows.append([x, f"{s.index[0]}-{s.index[-1]}", len(s), f"{100 * s.std():.2f}",
                     f"{s.autocorr(1):.2f}", f"{cyc[[x, 'u']].dropna().corr().iloc[0, 1]:+.3f}",
                     f"{cyc[[x, 'v']].dropna().corr().iloc[0, 1]:+.3f}"])
    text = ("Business-cycle moments, HP(1,600) log cycles, cut at 2019Q4 before filtering "
            "(q3_stylized_moments.py)\n\n"
            + tabulate(rows, ["series", "window", "T", "SD x 100", "rho_1", "cor with u", "cor with v"],
                       tablefmt="simple", disable_numparse=True) + "\n")
    lags = range(-4, 9)
    ll = []
    for x in ("z", "s", "delta"):
        ll.append([x] + [f"{cyc[x].corr(cyc['u'].shift(-k)):+.2f}" for k in lags])
    text += ("\nLead-lag correlations cor(x_t, u_{t+k}) (k > 0: x leads u), common sample of each pair\n\n"
             + tabulate(ll, ["x"] + [f"k={k}" for k in lags], tablefmt="simple", disable_numparse=True) + "\n")
    OUT.write_text(text, encoding="utf-8")
    print(text)
