"""
q3_delta_joint.py -- the delta shock net of demand and financial conditions together.

SVAR (z, CU, EBP, s, delta, u, v), recursive, 1992Q3-2019Q4: total capacity utilization
(q3_delta_cu.py) and the excess bond premium (q3_finance_delta.py) both ordered before delta. The
delta shock is then exit net of same-quarter and lagged utilization and credit conditions. The
order of CU and EBP between themselves does not affect the delta shock (both precede it in the
Cholesky ordering). An upper bound on contamination by demand and financial disturbances: the
two proxies take all same-quarter comovement with delta.
Reported against SVAR A on the same window and lag: delta -> u and v (peak, path), the share of the
cumulative response (h = 0-12) absorbed, with a joint residual-bootstrap band (SVAR A re-fitted on
each bootstrap sample), and delta's share of the forecast-error variance of u.
Specifications: HP, BIC lag and p = 4; Hamilton, BIC lag.

Output:
    data/results/q3_delta_joint.txt
"""

import numpy as np
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b
import q3_finance_delta as rev
from q3_delta_cu import cycles

NAMES   = ["z", "cu", "ebp", "s", "delta", "u", "v"]
H       = q3.H_IRF
OUT_TXT = p6b.RESULTS_DIR / "q3_delta_joint.txt"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def run(cyc, p=None):
    frame = cyc.loc[cyc.index >= q3.BED_START, NAMES].dropna()
    Y = frame.to_numpy()
    if p is None:
        p = q3.bic_lag(Y)[0]
    idxA = [NAMES.index(v) for v in rev.NAMES_A]
    fit, fitA = q3.fit_var(Y, p), q3.fit_var(Y[:, idxA], p)
    Th = q3.irf(fit, NAMES, NAMES, H)
    ThA = q3.irf(fitA, rev.NAMES_A, rev.NAMES_A, H)
    fev, fevA = q3.fevd(fit, NAMES, NAMES, (8,)), q3.fevd(fitA, rev.NAMES_A, rev.NAMES_A, (8,))
    return {"p": p, "frame": frame, "Th": Th, "ThA": ThA,
            "pt": rev.stats(Th, ThA, fit, fitA, NAMES),
            "draws": rev.bootstrap(Y, fit, NAMES, idxA, q3.N_BOOT),
            "fevd": (100 * fev[8][NAMES.index("u"), NAMES.index("delta")],
                     100 * fevA[8][rev.NAMES_A.index("u"), rev.NAMES_A.index("delta")])}


if __name__ == "__main__":
    say("=" * 78)
    say("The delta shock net of utilization and the excess bond premium (q3_delta_joint.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])
    cyc_hp, cyc_ham = cycles("hp"), cycles("hamilton")
    runs = []
    r = run(cyc_hp)
    runs.append((f"HP, p = {r['p']} (BIC)", r))
    if r["p"] != 4:
        runs.append(("HP, p = 4", run(cyc_hp, p=4)))
    r = run(cyc_ham)
    runs.append((f"Hamilton, p = {r['p']} (BIC)", r))
    fr, fh = runs[0][1]["frame"], runs[-1][1]["frame"]
    say(f"\nHP sample {fr.index[0]}-{fr.index[-1]}, T = {len(fr)}; Hamilton {fh.index[0]}-{fh.index[-1]}, "
        f"T = {len(fh)}. Bands: residual bootstrap, {q3.N_BOOT} draws, [68 %].\n")
    tab = []
    jd, iu, iv = NAMES.index("delta"), NAMES.index("u"), NAMES.index("v")
    jA, iuA, ivA = (rev.NAMES_A.index(k) for k in ("delta", "u", "v"))
    for lab, r in runs:
        row = [lab]
        for i, iA in ((iu, iuA), (iv, ivA)):
            pr, pa = r["Th"][:, i, jd], r["ThA"][:, iA, jA]
            kr, ka = int(np.abs(pr).argmax()), int(np.abs(pa).argmax())
            row += [f"{100 * pr[kr]:+.2f} (h {kr})", f"{100 * pa[ka]:+.2f} (h {ka})"]
        bu, bv = rev.band(r["draws"], "absorbed_u"), rev.band(r["draws"], "absorbed_v")
        row += [f"{r['pt']['absorbed_u']:.0f} [{bu[0]:.0f}, {bu[1]:.0f}]",
                f"{r['pt']['absorbed_v']:.0f} [{bv[0]:.0f}, {bv[1]:.0f}]",
                f"{r['fevd'][0]:.0f} / {r['fevd'][1]:.0f}"]
        tab.append(row)
    say(tabulate(tab, ["spec", "δ→u peak", "A peak", "δ→v peak", "A peak", "absorbed u", "absorbed v",
                       "FEVD δ→u 8q"], tablefmt="simple", disable_numparse=True))
    say("\n  absorbed = 1 − cum/cum_A over h = 0–12; FEVD: δ's share of var(u) at 8 quarters, this SVAR / A.")
    for lab, r in runs:
        say(f"\n  δ → u by quarter, {lab} (log points × 100):")
        pr = 100 * r["Th"][:, iu, jd]
        pa = 100 * r["ThA"][:, iuA, jA]
        say("    h        " + " ".join(f"{h:6d}" for h in range(H + 1)))
        say("    A        " + " ".join(f"{x:+6.2f}" for x in pa))
        say("    net both " + " ".join(f"{x:+6.2f}" for x in pr))
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}")
