"""
q3_z_purged.py -- how much of the z shock is demand? z purged of capacity utilization.

SVAR A orders z (labor productivity) first and reads its shock as technology. In goods-market
search models, measured productivity also moves with demand: in Silva and Urias (2025, Economic
Modelling 153, 107317, Table 4) shopping-effort shocks explain 46 % of the Solow residual's
forecast-error variance, and the Solow residual's growth equals utilization growth plus
technology plus input-share terms (their eq. 26). Ordering total capacity utilization (CU,
q3_delta_cu.py) before z purges the z shock of same-quarter utilization:
    SVAR P: (CU, z, s, delta, u, v), recursive, 1992Q3-2019Q4,
so the z shock is productivity net of utilization (closer to technology) and the CU shock is
demand. The model has no utilization, so its z is technology by construction: the purged data
response is the counterpart of the model's z shock (as A_F is a bound for delta).
Reported, against SVAR A on the same window and lag:
  1. z -> u and v: peak, cumulative (h = 0-12), per one-SD shock and per 1 % impact on z; the
     ratio cum(P)/cum(A) with a joint residual-bootstrap band (A re-fitted on each sample);
  2. the variance shares of u at 4 and 8 quarters for z, CU and delta;
  3. the AGS model's z -> u through SVAR A (q3_model_svar.csv, cases default and compromise), for
     the size of the gap against the original and the purged data responses;
  4. delta -> u in SVAR P (the delta shock is that of q3_delta_cu.py's first ordering).
Specifications: HP, BIC lag and p = 4; Hamilton, BIC lag.

Output:
    data/results/q3_z_purged.txt
"""

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b
from q3_delta_cu import cycles

NAMES   = ["cu", "z", "s", "delta", "u", "v"]
NAMES_A = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
H       = q3.H_IRF
MODEL   = p6b.RESULTS_DIR / "q3_model_svar.csv"
OUT_TXT = p6b.RESULTS_DIR / "q3_z_purged.txt"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def zstats(Th, ThA):
    """z -> u, v in P and A: paths per one SD and per 1 % impact on z."""
    jz, iz, iu, iv = NAMES.index("z"), NAMES.index("z"), NAMES.index("u"), NAMES.index("v")
    jzA, izA, iuA, ivA = (NAMES_A.index(k) for k in ("z", "z", "u", "v"))
    out = {"impact": (Th[0, iz, jz], ThA[0, izA, jzA])}
    for r, i, iA in (("u", iu, iuA), ("v", iv, ivA)):
        out[r] = (Th[:, i, jz], ThA[:, iA, jzA])
        out[f"ratio_{r}"] = Th[:, i, jz].sum() / ThA[:, iA, jzA].sum()
        out[f"ratio_{r}_unit"] = (Th[:, i, jz].sum() / Th[0, iz, jz]) / (ThA[:, iA, jzA].sum() / ThA[0, izA, jzA])
    return out


def run(cyc, p=None, n_boot=q3.N_BOOT, seed=q3.SEED):
    frame = cyc.loc[cyc.index >= q3.BED_START, NAMES].dropna()
    Y = frame.to_numpy()
    if p is None:
        p = q3.bic_lag(Y)[0]
    idxA = [NAMES.index(v) for v in NAMES_A]
    fit, fitA = q3.fit_var(Y, p), q3.fit_var(Y[:, idxA], p)
    Th, ThA = q3.irf(fit, NAMES, NAMES, H), q3.irf(fitA, NAMES_A, NAMES_A, H)
    fev = q3.fevd(fit, NAMES, NAMES, (4, 8))
    fevA = q3.fevd(fitA, NAMES_A, NAMES_A, (4, 8))
    # joint bootstrap
    rng = np.random.default_rng(seed)
    U = fit["U"] - fit["U"].mean(axis=0)
    T = U.shape[0]
    draws = []
    for _ in range(n_boot):
        e = U[rng.integers(0, T, T)]
        Ys = np.empty_like(Y)
        Ys[:p] = Y[:p]
        for t in range(p, p + T):
            Ys[t] = fit["c"] + sum(fit["A"][j] @ Ys[t - 1 - j] for j in range(p)) + e[t - p]
        draws.append(zstats(q3.irf(q3.fit_var(Ys, p), NAMES, NAMES, H),
                            q3.irf(q3.fit_var(Ys[:, idxA], p), NAMES_A, NAMES_A, H)))
    return {"p": p, "frame": frame, "Th": Th, "ThA": ThA, "pt": zstats(Th, ThA), "draws": draws,
            "fev": fev, "fevA": fevA}


def pk(path):
    k = int(np.abs(path).argmax())
    return f"{100 * path[k]:+.2f} (h {k})"


if __name__ == "__main__":
    say("=" * 78)
    say("How much of the z shock is demand? z purged of capacity utilization (q3_z_purged.py)")
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
        f"T = {len(fh)}. Bands: joint residual bootstrap, {q3.N_BOOT} draws, [68 %].")

    say(f"\n{'-'*78}\n1. z -> u and v: purged (P) against SVAR A\n{'-'*78}")
    say("  impact: z's own response on impact (one-SD shock); peak per one SD; ratio = cum_P/cum_A over")
    say("  h = 0–12, per one SD and per 1 % impact on z [68 %].\n")
    tab = []
    for lab, r in runs:
        pt = r["pt"]
        for resp in ("u", "v"):
            b1 = np.percentile([d[f"ratio_{resp}"] for d in r["draws"]], [16, 84])
            b2 = np.percentile([d[f"ratio_{resp}_unit"] for d in r["draws"]], [16, 84])
            tab.append([lab, resp, f"{100 * pt['impact'][0]:.2f} / {100 * pt['impact'][1]:.2f}",
                        pk(pt[resp][0]), pk(pt[resp][1]),
                        f"{100 * pt[resp][0].sum():+.1f} / {100 * pt[resp][1].sum():+.1f}",
                        f"{pt[f'ratio_{resp}']:.2f} [{b1[0]:.2f}, {b1[1]:.2f}]",
                        f"{pt[f'ratio_{resp}_unit']:.2f} [{b2[0]:.2f}, {b2[1]:.2f}]"])
    say(tabulate(tab, ["spec", "", "z impact P / A", "peak P", "peak A", "cum P / A", "ratio, one SD",
                       "ratio, per 1 %"], tablefmt="simple", disable_numparse=True))

    say(f"\n{'-'*78}\n2. Variance shares of u (FEVD, %), 4q / 8q\n{'-'*78}\n")
    tab = []
    iu, iuA = NAMES.index("u"), NAMES_A.index("u")
    for lab, r in runs:
        f, fA = r["fev"], r["fevA"]
        cell = lambda F, i, j: f"{100 * F[4][i, j]:.0f} / {100 * F[8][i, j]:.0f}"
        tab.append([lab, cell(f, iu, NAMES.index("z")), cell(fA, iuA, NAMES_A.index("z")),
                    cell(f, iu, NAMES.index("cu")), cell(f, iu, NAMES.index("delta")),
                    cell(fA, iuA, NAMES_A.index("delta"))])
    say(tabulate(tab, ["spec", "z, P", "z, A", "CU, P", "δ, P", "δ, A"], tablefmt="simple",
                 disable_numparse=True))

    say(f"\n{'-'*78}\n3. Against the AGS model (model side: SVAR A on simulated windows)\n{'-'*78}")
    say("  z -> u by quarter, one-SD shock, log points × 100; model: mean over windows (q3_model_svar.csv);")
    say("  data purged rescaled to A's impact on z (per 1 % comparison) in the last row.\n")
    r = runs[0][1]
    pt = r["pt"]
    md = pd.read_csv(MODEL)
    rows = {"data, A": 100 * pt["u"][1], "data, purged (P)": 100 * pt["u"][0],
            "data, purged, rescaled": 100 * pt["u"][0] * pt["impact"][1] / pt["impact"][0]}
    for case in ("default", "compromise"):
        g = md[(md.case == case) & (md.shock == "z") & (md.response == "u")].sort_values("h")
        rows[f"model, {case}"] = 100 * g.model_mean.to_numpy()[:H + 1]
    say("    " + "".ljust(24) + " ".join(f"{h:6d}" for h in range(H + 1)) + "    cum")
    for k, v in rows.items():
        say("    " + k.ljust(24) + " ".join(f"{x:+6.2f}" for x in v) + f"  {v.sum():+6.1f}")
    say(f"\n  Model impact on z (default): {100 * md[(md.case == 'default') & (md.shock == 'z') & (md.response == 'z') & (md.h == 0)].model_mean.iloc[0]:.2f}; "
        f"data A {100 * pt['impact'][1]:.2f}; purged {100 * pt['impact'][0]:.2f} ({runs[0][0]}).")

    say(f"\n{'-'*78}\n4. delta -> u in SVAR P\n{'-'*78}")
    for lab, r in runs:
        jd = NAMES.index("delta")
        pP = r["Th"][:, NAMES.index("u"), jd]
        pA = r["ThA"][:, NAMES_A.index("u"), NAMES_A.index("delta")]
        say(f"  {lab}: peak {pk(pP)} (A {pk(pA)}); cum {100 * pP.sum():+.1f} (A {100 * pA.sum():+.1f})")
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}")
