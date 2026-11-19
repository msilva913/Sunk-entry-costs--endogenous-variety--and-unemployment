"""
q3_model_svar_components.py -- SVAR A on model data when measured delta is a persistent
component plus measurement noise (interpretation (a); findings.md, Nov 19).

Model side: ../Programs AGS/ags_svar_components.jl, one case per noise share w on part6d's
profile (shock_components_profile.csv): the economy's delta shock is the persistent component,
and the measured delta the SVAR sees carries the noise. z and s keep Block B''s processes and
the calibration is the default, so the comparison with q3_model_svar.py's default case isolates
the delta process. Same SVAR code, lag (data BIC) and filter as the data.

Reported for the delta shock: delta's own response, u and v per one SD and per 1 % of own
impact, peak and quarter, cumulative over h = 0-12, the distance J_delta to the data (u and v,
h = 0-12), and the true (unfiltered, noise-free) response of u to the persistent shock.

Outputs:
    data/results/q3_model_svar_components.txt
    data/results/q3_model_svar_components.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

MODEL_DIR = Path("../Programs AGS")
CASES     = pd.read_csv(MODEL_DIR / "ags_svar_components_cases.csv").set_index("case")
BASELINE  = pd.read_csv(p6b.RESULTS_DIR / "q3_model_svar.csv")
NAMES     = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
H         = q3.H_IRF
OUT_TXT   = p6b.RESULTS_DIR / "q3_model_svar_components.txt"
OUT_CSV   = p6b.RESULTS_DIR / "q3_model_svar_components.csv"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def svar_irf(frame, p):
    Y = frame[NAMES].to_numpy()
    return q3.irf(q3.fit_var(Y, p), NAMES, NAMES, H)


def peak(path):
    k = int(np.abs(path).argmax())
    return f"{path[k]:+.2f} (h {k})"


if __name__ == "__main__":
    cyc = q3.build_cycles("hp")
    svA = q3.prepare_svar("A", q3.SVARS["A"], cyc)
    p = svA["p_bic"]
    fitA = q3.fit_var(svA["Y"], p)
    Theta_d = q3.irf(fitA, NAMES, NAMES, H)
    draws = q3.bootstrap_irf(svA["Y"], fitA, NAMES, NAMES, H, q3.N_BOOT)
    sd_d = draws.std(axis=0)
    jd, iu, iv = NAMES.index("delta"), NAMES.index("u"), NAMES.index("v")

    say("=" * 78)
    say("SVAR A with measured delta = persistent component + noise (q3_model_svar_components.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])
    say(f"\nData: SVAR A, T = {len(svA['frame'])}, p = {p}; {q3.N_BOOT} bootstrap draws for the data SDs.")

    paths = {"data": {r: 100 * Theta_d[:, NAMES.index(r), jd] for r in ("delta", "u", "v")}}
    base = BASELINE[(BASELINE.case == "default") & (BASELINE.shock == "delta")]
    paths["baseline (one AR(1))"] = {r: 100 * base[base.response == r].sort_values("h").model_mean.to_numpy()
                                     for r in ("delta", "u", "v")}
    true_u = {}
    for case in CASES.index:
        reps = pd.read_csv(MODEL_DIR / f"ags_svar_reps_{case}.csv")
        rep_irf = np.array([svar_irf(g, p) for _, g in reps.groupby("rep")])
        m = rep_irf.mean(axis=0)
        paths[case] = {r: 100 * m[:, NAMES.index(r), jd] for r in ("delta", "u", "v")}
        tru = pd.read_csv(MODEL_DIR / f"ags_svar_true_{case}.csv").sort_values("h")
        true_u[case] = 100 * tru.u.to_numpy()[:H + 1] / tru.delta.to_numpy()[0]

    say(f"\n{'-' * 78}\nDelta shock: SVAR responses (one SD; log points x 100), h = 0-12\n{'-' * 78}")
    rows, tab = [], []
    for lab, P in paths.items():
        own = P["delta"][0]
        J_text = ""
        if lab != "data":
            J = float((((P["u"] - paths["data"]["u"]) / (100 * sd_d[:, iu, jd])) ** 2).sum()
                      + (((P["v"] - paths["data"]["v"]) / (100 * sd_d[:, iv, jd])) ** 2).sum())
            J_text = f"{J:.1f}"
        tab.append([lab, f"{own:.2f}", f"{P['delta'][2] / own:+.2f}", f"{P['delta'][4] / own:+.2f}",
                    peak(P["u"]), f"{P['u'].sum():+.1f}", peak(P["u"] / own), peak(P["v"]), J_text])
        for h in range(H + 1):
            rows.append({"case": lab, "h": h, "delta": P["delta"][h], "u": P["u"][h], "v": P["v"][h]})
    say(tabulate(tab, ["case", "own impact", "own h2/h0", "own h4/h0", "u peak", "u cum 0-12",
                       "u peak per 1 %", "v peak", "J δ"], tablefmt="simple", disable_numparse=True))
    say("\n  J δ: Σ ((model − data)/sd_data)² over u and v, h = 0–12 (26 moments).")

    say(f"\n  u by quarter (one SD):")
    say("    h            " + " ".join(f"{h:6d}" for h in range(H + 1)))
    for lab, P in paths.items():
        say(f"    {lab[:12]:<12} " + " ".join(f"{x:+6.2f}" for x in P["u"]))

    say(f"\n{'-' * 78}\nTrue response of u to the persistent delta shock (not filtered, no noise),\n"
        f"per 1 % of delta on impact\n{'-' * 78}")
    say("    h            " + " ".join(f"{h:6d}" for h in range(H + 1)))
    for case, x in true_u.items():
        say(f"    {case:<12} " + " ".join(f"{v:+6.2f}" for v in x))
    say()
    for case in CASES.index:
        c = CASES.loc[case]
        say(f"  {case}: persistent δ ρ_m {c.rho_delta:.4f} (half-life {c.half_life_m:.0f} months), "
            f"σ {c.sigma_delta:.4f}; noise σ {c.sigma_noise:.4f}")

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}")
