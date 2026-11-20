"""
q3_model_svar_filters.py -- the overidentifying test through SVAR A, under HP and Hamilton.

Model side: ../Programs AGS/ags_svar_filters.jl, reactivation rates 1 to 0.03, default calibration
otherwise, shock processes matched per filter on SVAR A's window (part6e). Data side: SVAR A on
the same filter's cycles (q3_data_svar.py), lag by BIC on the data, used on both sides.

The test: the separation responses identify the reactivation rate (the case with the smallest
J_s); the destruction response at that rate is the test (J_delta, u's peak and cumulative
response). Every case's J for all three shocks is reported, so the identification and the test
can be read together. A validation block reports how well each filter's SVAR recovers the
model's true responses (long simulated sample, default reactivation rate).

J: sum over u and v and h = 0-12 of ((model - data)/sd_data)^2, 26 moments per shock; data SDs
from a 1,000-draw residual bootstrap under the same filter.

Outputs:
    data/results/q3_model_svar_filters.txt
    data/results/q3_model_svar_filters.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

MODEL_DIR = Path("../Programs AGS")
CASES     = pd.read_csv(MODEL_DIR / "ags_svar_filters_cases.csv")
NAMES     = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
SHOCKS    = ("z", "delta", "s")
H         = q3.H_IRF
FILTERS   = {"hp": ("HP-1600", "hp"), "ham": ("Hamilton (h 8, p 4)", "hamilton")}
OUT_TXT   = p6b.RESULTS_DIR / "q3_model_svar_filters.txt"
OUT_CSV   = p6b.RESULTS_DIR / "q3_model_svar_filters.csv"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def svar_irf(frame, p):
    Y = frame[NAMES].dropna().to_numpy()
    return q3.irf(q3.fit_var(Y, p), NAMES, NAMES, H)


def J(model, data, sd, shock):
    j = NAMES.index(shock)
    out = 0.0
    for resp in ("u", "v"):
        i = NAMES.index(resp)
        out += float((((model[:, i, j] - data[:, i, j]) / sd[:, i, j]) ** 2).sum())
    return out


def u_peak(Theta, shock):
    path = 100 * Theta[:, NAMES.index("u"), NAMES.index(shock)]
    k = int(np.abs(path).argmax())
    return path[k], k, path.sum()


if __name__ == "__main__":
    say("=" * 78)
    say("The overidentifying test through SVAR A under HP and Hamilton (q3_model_svar_filters.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])
    rows = []
    for fkey, (flabel, method) in FILTERS.items():
        cyc = q3.build_cycles(method)
        svA = q3.prepare_svar("A", q3.SVARS["A"], cyc)
        p = svA["p_bic"]
        Yd = svA["Y"]
        fitA = q3.fit_var(Yd, p)
        Theta_d = q3.irf(fitA, NAMES, NAMES, H)
        sd = q3.bootstrap_irf(Yd, fitA, NAMES, NAMES, H, q3.N_BOOT).std(axis=0)
        say(f"\n{'=' * 78}\n{flabel}: data {svA['frame'].index[0]}-{svA['frame'].index[-1]}, "
            f"T = {len(svA['frame'])}, p = {p}\n{'=' * 78}")
        tab = []
        d_rows = {s: u_peak(Theta_d, s) for s in SHOCKS}
        tab.append(["data", "", "", "", "",
                    f"{d_rows['z'][0]:+.2f} (h {d_rows['z'][1]})",
                    f"{d_rows['delta'][0]:+.2f} (h {d_rows['delta'][1]}), {d_rows['delta'][2]:+.1f}",
                    f"{d_rows['s'][0]:+.2f} (h {d_rows['s'][1]}), {d_rows['s'][2]:+.1f}",
                    f"{100 * Theta_d[4, NAMES.index('v'), NAMES.index('s')]:+.2f}"])
        cases = CASES[CASES["filter"] == fkey]
        for _, c in cases.iterrows():
            reps = pd.read_csv(MODEL_DIR / f"ags_svar_reps_{fkey}_{c.case}.csv")
            m = np.mean([svar_irf(g, p) for _, g in reps.groupby("rep")], axis=0)
            Js = {s: J(m, Theta_d, sd, s) for s in SHOCKS}
            pk = {s: u_peak(m, s) for s in SHOCKS}
            sv4 = 100 * m[4, NAMES.index("v"), NAMES.index("s")]
            tab.append([c.case, f"{c.lambda_r:.2f}", f"{Js['z']:.0f}", f"{Js['delta']:.0f}", f"{Js['s']:.0f}",
                        f"{pk['z'][0]:+.2f} (h {pk['z'][1]})",
                        f"{pk['delta'][0]:+.2f} (h {pk['delta'][1]}), {pk['delta'][2]:+.1f}",
                        f"{pk['s'][0]:+.2f} (h {pk['s'][1]}), {pk['s'][2]:+.1f}", f"{sv4:+.2f}"])
            rows.append({"filter": fkey, "case": c.case, "lambda_r": c.lambda_r, "phi": c.phi,
                         "J_z": Js["z"], "J_delta": Js["delta"], "J_s": Js["s"],
                         "delta_u_peak": pk["delta"][0], "delta_u_hpeak": pk["delta"][1],
                         "delta_u_cum": pk["delta"][2], "s_u_peak": pk["s"][0], "s_u_hpeak": pk["s"][1],
                         "s_u_cum": pk["s"][2], "z_u_peak": pk["z"][0], "s_v_h4": sv4,
                         "data_delta_u_cum": d_rows["delta"][2], "data_s_u_cum": d_rows["s"][2]})
        say(tabulate(tab, ["case", "Λ̄_r", "J z", "J δ", "J s", "z→u peak", "δ→u peak, cum",
                           "s→u peak, cum", "s→v h4"], tablefmt="simple", disable_numparse=True))
        R = pd.DataFrame([r for r in rows if r["filter"] == fkey])
        best = R.loc[R.J_s.idxmin()]
        say(f"\n  Identified by s: Λ̄_r = {best.lambda_r:.2f} (J_s {best.J_s:.0f}). At that rate: J_δ "
            f"{best.J_delta:.0f}; δ→u cum {best.delta_u_cum:+.1f} against the data's "
            f"{best.data_delta_u_cum:+.1f}; s→u cum {best.s_u_cum:+.1f} against {best.data_s_u_cum:+.1f}.")
        say(f"  Smallest J_δ over the grid: {R.J_delta.min():.0f} (Λ̄_r = {R.loc[R.J_delta.idxmin()].lambda_r:.2f}).")
        say("  26 moments per shock: J ≈ 26 is a fit within one data SD on average.")

        # validation: long sample SVAR against the true responses, Λ̄_r = 0.25
        long = pd.read_csv(MODEL_DIR / f"ags_svar_long_{fkey}_L0.25.csv")
        Th_long = svar_irf(long, p)
        tru = pd.read_csv(MODEL_DIR / f"ags_svar_true_{fkey}_L0.25.csv")
        vt = []
        for s in SHOCKS:
            g = tru[tru.shock == s].sort_values("h")
            t_u = g.u.to_numpy()[:H + 1] / g[s].to_numpy()[0]
            j = NAMES.index(s)
            l_u = Th_long[:, NAMES.index("u"), j] / Th_long[0, j, j]
            gap = np.abs(l_u - t_u).max() / np.abs(t_u).max()
            vt.append([s, f"{t_u[np.abs(t_u).argmax()]:+.2f} (h {np.abs(t_u).argmax()})",
                       f"{l_u[np.abs(l_u).argmax()]:+.2f} (h {np.abs(l_u).argmax()})", f"{100 * gap:.0f} %"])
        say(f"\n  Validation (Λ̄_r = 0.25, long sample, p = {p}): u per unit of own impact")
        say(tabulate(vt, ["shock", "true peak", "SVAR peak", "gap"], tablefmt="simple", disable_numparse=True))

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}")
