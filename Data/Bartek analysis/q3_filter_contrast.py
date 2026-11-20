"""
q3_filter_contrast.py -- does the destruction-separation contrast survive the Hamilton filter?

The draft's empirical core: after a destruction shock u builds to a peak at quarter 4, mostly
through a lower job-finding rate, whereas a separation shock acts sooner and fades sooner. Under
HP (SVAR A) s -> u peaks at h 2 and delta -> u at h 4; under the Hamilton filter both peak at
h 5 (q3_data_svar.txt section 7). This program asks, on the data only, whether the contrast
holds under each filter, and separates the filter from the sample.

Three specifications, SVAR A (z, s, delta, u, v) at p = 1, one-SD shocks:
    HP          HP-1600 cycles, 1992Q3-2019Q4 (the baseline)
    HP, 1995Q2  the same cycles on the Hamilton window (sample effect only)
    Hamilton    Hamilton (h = 8, p = 4) cycles, 1995Q2-2019Q4 (filter + sample)

For delta and s, per specification: u's peak (value, quarter), u(h 8)/peak, cumulative u over
h = 0-12, and the inflow/outflow split of u from SVAR (z, s, delta, u, v, f) with the steady-
state flow identity (as q3_persistence_tests.py test 2): the outflow share of the cumulative
response. The contrast statistics are computed draw by draw in one residual bootstrap
(1,000 draws), so their bands account for the correlation between the two responses:
    peak quarter (delta) - peak quarter (s)
    h8/peak (delta) - h8/peak (s)                       persistence after the peak
    outflow share (delta) - outflow share (s)           composition
Lag robustness: the same statistics at p = 2 (point estimates).

Outputs:
    data/results/q3_filter_contrast.txt
    data/results/q3_filter_contrast.csv
"""

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

NAMES   = list(q3.SVARS["A"]["vars"])           # z, s, delta, u, v
NAMES_F = NAMES + ["jf"]
H       = q3.H_IRF
W_DELTA, U_BAR = 0.087, 0.06                    # as q3_persistence_tests.py
HAM_START = "1995Q2"
OUT_TXT = p6b.RESULTS_DIR / "q3_filter_contrast.txt"
OUT_CSV = p6b.RESULTS_DIR / "q3_filter_contrast.csv"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def stats(Theta, names, shock):
    """u's peak and quarter, h8/peak, cumulative, and the outflow share (needs jf in names)."""
    j, iu = names.index(shock), names.index("u")
    u = Theta[:, iu, j]
    k = int(np.abs(u).argmax())
    out = {"peak": 100 * u[k], "h_peak": k, "ratio8": u[8] / u[k], "cum": 100 * u.sum()}
    if "jf" in names:
        tau = W_DELTA * Theta[:, names.index("delta"), j] + (1 - W_DELTA) * Theta[:, names.index("s"), j]
        inflow = (1 - U_BAR) * tau.sum()
        outflow = -(1 - U_BAR) * Theta[:, names.index("jf"), j].sum()
        out["outflow_share"] = outflow / (inflow + outflow)
    return out


def contrast(sd, ss):
    """delta minus s."""
    return {"d_hpeak": sd["h_peak"] - ss["h_peak"], "d_ratio8": sd["ratio8"] - ss["ratio8"],
            "d_outflow": sd["outflow_share"] - ss["outflow_share"]}


def spec_frame(cyc, start):
    return cyc.loc[cyc.index >= start, NAMES_F].dropna()


if __name__ == "__main__":
    cyc_hp, cyc_ham = q3.build_cycles("hp"), q3.build_cycles("hamilton")
    specs = {"HP": spec_frame(cyc_hp, q3.BED_START),
             "HP, 1995Q2": spec_frame(cyc_hp, HAM_START),
             "Hamilton": spec_frame(cyc_ham, HAM_START)}

    say("=" * 78)
    say("The destruction-separation contrast under HP and Hamilton (q3_filter_contrast.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])

    rows, point_tab, con_tab = [], [], []
    for lab, fr in specs.items():
        Y5, Y6 = fr[NAMES].to_numpy(), fr[NAMES_F].to_numpy()
        fit5, fit6 = q3.fit_var(Y5, 1), q3.fit_var(Y6, 1)
        Th5 = q3.irf(fit5, NAMES, NAMES, H)
        Th6 = q3.irf(fit6, NAMES_F, NAMES_F, H)
        # u statistics from SVAR A; the outflow share from the six-variable system
        pt = {}
        for s in ("delta", "s"):
            st = stats(Th5, NAMES, s)
            st["outflow_share"] = stats(Th6, NAMES_F, s)["outflow_share"]
            pt[s] = st
        c_pt = contrast(pt["delta"], pt["s"])
        D5 = q3.bootstrap_irf(Y5, fit5, NAMES, NAMES, H, q3.N_BOOT)
        D6 = q3.bootstrap_irf(Y6, fit6, NAMES_F, NAMES_F, H, q3.N_BOOT)
        draws = []
        for b in range(q3.N_BOOT):
            sd, ss = stats(D5[b], NAMES, "delta"), stats(D5[b], NAMES, "s")
            sd["outflow_share"] = stats(D6[b], NAMES_F, "delta")["outflow_share"]
            ss["outflow_share"] = stats(D6[b], NAMES_F, "s")["outflow_share"]
            draws.append(contrast(sd, ss))
        draws = pd.DataFrame(draws)
        # p = 2 point estimates
        Th5_2 = q3.irf(q3.fit_var(Y5, 2), NAMES, NAMES, H)
        Th6_2 = q3.irf(q3.fit_var(Y6, 2), NAMES_F, NAMES_F, H)
        p2 = {}
        for s in ("delta", "s"):
            st = stats(Th5_2, NAMES, s)
            st["outflow_share"] = stats(Th6_2, NAMES_F, s)["outflow_share"]
            p2[s] = st
        c_p2 = contrast(p2["delta"], p2["s"])

        for s in ("delta", "s"):
            r = pt[s]
            point_tab.append([lab, f"{fr.index[0]}-{fr.index[-1]} (T {len(fr)})", q3.LABEL[s],
                              f"{r['peak']:+.2f} (h {r['h_peak']})", f"{r['ratio8']:+.2f}",
                              f"{r['cum']:+.1f}", f"{100 * r['outflow_share']:.0f} %",
                              f"h {p2[s]['h_peak']}, {100 * p2[s]['outflow_share']:.0f} %"])
            rows.append({"spec": lab, "shock": s, **r, "p2_h_peak": p2[s]["h_peak"],
                         "p2_outflow_share": p2[s]["outflow_share"]})
        for key, name, fmt in (("d_hpeak", "peak quarter, δ − s", "{:+.0f}"),
                               ("d_ratio8", "h8/peak, δ − s", "{:+.2f}"),
                               ("d_outflow", "outflow share, δ − s", "{:+.2f}")):
            lo68, hi68 = np.percentile(draws[key], [16, 84])
            lo90, hi90 = np.percentile(draws[key], [5, 95])
            share_pos = float((draws[key] > 0).mean())
            con_tab.append([lab, name, fmt.format(c_pt[key]), f"[{fmt.format(lo68)}, {fmt.format(hi68)}]",
                            f"[{fmt.format(lo90)}, {fmt.format(hi90)}]", f"{share_pos:.2f}",
                            fmt.format(c_p2[key])])
            rows.append({"spec": lab, "shock": "contrast", "stat": key, "point": c_pt[key],
                         "lo68": lo68, "hi68": hi68, "lo90": lo90, "hi90": hi90,
                         "share_pos": share_pos, "p2": c_p2[key]})

    say(f"\n{'-' * 78}\n1. u's response to δ and s (one SD, p = 1; log points × 100)\n{'-' * 78}")
    say(tabulate(point_tab, ["spec", "window", "shock", "u peak", "h8/peak", "cum 0–12",
                             "outflow share", "p = 2: peak, outflow"],
                 tablefmt="simple", disable_numparse=True))
    say("\n  Outflow share: outflow part of the cumulative u response (h = 0–12) over inflow + outflow,")
    say("  from SVAR (z, s, δ, u, v, f); above 100 % when the inflow part is negative.")

    say(f"\n{'-' * 78}\n2. The contrast, δ minus s: point, 68 % and 90 % bootstrap bands, share of draws > 0\n{'-' * 78}")
    say(tabulate(con_tab, ["spec", "statistic", "point", "68 %", "90 %", "P(>0)", "p = 2"],
                 tablefmt="simple", disable_numparse=True))

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}")
