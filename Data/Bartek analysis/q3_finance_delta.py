"""
q3_finance_delta.py -- do financial shocks move delta, and how much of delta's effect survives them?

The reverse of q3_delta_finance.py. A financial variable F is ordered after z and before s and
delta:
    (z, F, s, delta, u, v), recursive, 1992Q3-2019Q4,
so F can move delta within the quarter, and the delta shock is exit net of same-quarter and
lagged financial conditions. This is A_F of q3_data_svar.py, extended to the excess bond premium.
    F = ebp        Gilchrist-Zakrajsek excess bond premium (headline; part7d_ebp.py)
        gz_spread  the GZ credit spread (EBP plus expected default)
        nfci       Chicago Fed NFCI (as in A_F)
        anfci      NFCI adjusted for macro conditions (as in A_Fa)
Cycles as in q3_delta_finance.py (F HP-filtered in levels over its own window).

Questions (Notes/ags_equilibrium.md section 10; the next steps of Nov 15):
  1. Do financial shocks raise delta, on impact or with a lag? The impact response is the
     correlation of the F and delta innovations, so it does not depend on the ordering's
     direction; the lagged response does not either, given the VAR.
  2. How much of delta's effect on u and v is absorbed once F is ordered first? Cumulative
     delta -> u, v (h = 0-12) against SVAR A on the same window and lag; share absorbed =
     1 - cum(this SVAR)/cum(A); and delta's share of the forecast-error variance of u.
  3. How much of a financial shock's effect on u runs through delta? Mediation as in
     q3_data_svar.py (delta's response held at zero by offsetting delta shocks).
Bands: residual bootstrap of the six-variable VAR; SVAR A is re-fitted on the same bootstrap
sample (its five variables), so the share absorbed has a joint band.
Specifications: BIC lag and p = 4 on HP cycles; BIC lag on Hamilton cycles (EBP).

Outputs:
    data/results/q3_finance_delta.txt
    data/results/q3_finance_delta.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b
from q3_delta_finance import cycles, FIN_VARS

H        = q3.H_IRF
N_BOOT   = q3.N_BOOT
NAMES_A  = list(q3.SVARS["A"]["vars"])        # z, s, delta, u, v
FFFF_H   = (4, 8)                             # FEVD horizons
OUT_TXT  = p6b.RESULTS_DIR / "q3_finance_delta.txt"
OUT_PNG  = p6b.RESULTS_DIR / "q3_finance_delta.png"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def names_for(F):
    return ["z", F, "s", "delta", "u", "v"]


def stats(Th, ThA, fit, fitA, names):
    """Everything reported for one draw (or the point estimate)."""
    F = names[1]
    jF, jd, iu, iv = names.index(F), names.index("delta"), names.index("u"), names.index("v")
    jdA, iuA, ivA = NAMES_A.index("delta"), NAMES_A.index("u"), NAMES_A.index("v")
    out = {"F_delta": Th[:, jd, jF]}
    for r, i, iA in (("u", iu, iuA), ("v", iv, ivA)):
        cum, cumA = Th[:, i, jd].sum(), ThA[:, iA, jdA].sum()
        out[f"absorbed_{r}"] = 100 * (1 - cum / cumA)
    full = Th[:, iu, jF]
    cf = q3.mediated(Th, jF, jd)[:, iu]
    out["F_u"] = full
    out["F_u_cf"] = cf
    out["via_delta_u"] = 100 * (1 - cf.sum() / full.sum())
    return out


def bootstrap(Y, fit, names, idxA, n_boot, seed=q3.SEED):
    """Residual bootstrap of the six-variable VAR; SVAR A re-fitted on each sample."""
    rng = np.random.default_rng(seed)
    p, c, A = fit["p"], fit["c"], fit["A"]
    U = fit["U"] - fit["U"].mean(axis=0)
    T = U.shape[0]
    draws = []
    for _ in range(n_boot):
        e = U[rng.integers(0, T, T)]
        Ys = np.empty_like(Y)
        Ys[:p] = Y[:p]
        for t in range(p, p + T):
            Ys[t] = c + sum(A[j] @ Ys[t - 1 - j] for j in range(p)) + e[t - p]
        f6, f5 = q3.fit_var(Ys, p), q3.fit_var(Ys[:, idxA], p)
        draws.append(stats(q3.irf(f6, names, names, H), q3.irf(f5, NAMES_A, NAMES_A, H), f6, f5, names))
    return draws


def run(cyc, F, p=None):
    names = names_for(F)
    frame = cyc.loc[cyc.index >= q3.BED_START, names].dropna()
    Y = frame.to_numpy()
    if p is None:
        p = q3.bic_lag(Y)[0]
    idxA = [names.index(v) for v in NAMES_A]
    fit, fitA = q3.fit_var(Y, p), q3.fit_var(Y[:, idxA], p)
    Th, ThA = q3.irf(fit, names, names, H), q3.irf(fitA, NAMES_A, NAMES_A, H)
    U = fit["U"]
    rho = np.corrcoef(U[:, names.index(F)], U[:, names.index("delta")])[0, 1]
    fev = q3.fevd(fit, names, names, FFFF_H)
    fevA = q3.fevd(fitA, NAMES_A, NAMES_A, FFFF_H)
    share = {h: (100 * fev[h][names.index("u"), names.index("delta")],
                 100 * fevA[h][NAMES_A.index("u"), NAMES_A.index("delta")],
                 100 * fev[h][names.index("u"), names.index(F)]) for h in FFFF_H}
    return {"F": F, "p": p, "frame": frame, "names": names, "rho": rho, "share": share,
            "Th": Th, "ThA": ThA, "pt": stats(Th, ThA, fit, fitA, names),
            "draws": bootstrap(Y, fit, names, idxA, N_BOOT)}


def band(draws, key, h=None, q=(16, 84)):
    if h is None:
        x = [d[key] for d in draws]
    else:
        x = [d[key][h] for d in draws]
    return np.percentile(x, q)


if __name__ == "__main__":
    say("=" * 78)
    say("Do financial shocks move delta? (q3_finance_delta.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])

    cyc_hp, cyc_ham = cycles("hp"), cycles("hamilton")
    runs = []
    for F in FIN_VARS:
        r = run(cyc_hp, F)
        runs.append((f"HP, p = {r['p']} (BIC)", r))
        if r["p"] != 4:
            runs.append(("HP, p = 4", run(cyc_hp, F, p=4)))
    r = run(cyc_ham, "ebp")
    runs.append((f"Hamilton, p = {r['p']} (BIC)", r))
    fr, fh = runs[0][1]["frame"], runs[-1][1]["frame"]
    say(f"\nHP sample {fr.index[0]}-{fr.index[-1]}, T = {len(fr)}; Hamilton {fh.index[0]}-{fh.index[-1]}, "
        f"T = {len(fh)}. Bands: residual bootstrap, {N_BOOT} draws; [68 %].")

    # 1. F -> delta
    say(f"\n{'-'*78}\n1. Do financial shocks move delta? One-SD F shock -> delta (log points × 100)\n{'-'*78}")
    say("  corr: correlation of the F and delta innovations (the impact link, independent of the")
    say("  ordering's direction). peak: largest |response| and its quarter [68 %, 90 % at that quarter].\n")
    tab = []
    for lab, r in runs:
        path = r["pt"]["F_delta"]
        k = int(np.abs(path).argmax())
        b68 = 100 * band(r["draws"], "F_delta", k)
        b90 = 100 * band(r["draws"], "F_delta", k, (5, 95))
        tab.append([lab, r["F"], f"{r['rho']:+.2f}", f"{100 * path[0]:+.2f}",
                    f"{100 * path[k]:+.2f} (h {k})", f"[{b68[0]:+.2f}, {b68[1]:+.2f}]",
                    f"[{b90[0]:+.2f}, {b90[1]:+.2f}]", f"{100 * path.sum():+.1f}"])
    say(tabulate(tab, ["spec", "F", "corr", "impact", "peak", "68 %", "90 %", "cum 0–12"],
                 tablefmt="simple", disable_numparse=True))
    say("\n  For scale: a one-SD delta shock raises delta by about 5.2 on impact (SVAR A).")

    # 2. delta's effect net of F
    say(f"\n{'-'*78}\n2. delta's effect on u and v, net of financial conditions\n{'-'*78}")
    say("  peak and cum (h = 0–12) of delta -> u in this SVAR and in SVAR A (same window and lag);")
    say("  absorbed = 1 − cum/cum_A, % [68 %]; FEVD: delta's share of var(u) at 4q and 8q, this SVAR /")
    say("  SVAR A, and F's own share.\n")
    tab = []
    for lab, r in runs:
        names = r["names"]
        jd, iu = names.index("delta"), names.index("u")
        pr, pa = r["Th"][:, iu, jd], r["ThA"][:, NAMES_A.index("u"), NAMES_A.index("delta")]
        kr, ka = int(np.abs(pr).argmax()), int(np.abs(pa).argmax())
        bu, bv = band(r["draws"], "absorbed_u"), band(r["draws"], "absorbed_v")
        sh = r["share"]
        tab.append([lab, r["F"], f"{100 * pr[kr]:+.2f} (h {kr})", f"{100 * pa[ka]:+.2f} (h {ka})",
                    f"{100 * pr.sum():+.1f} / {100 * pa.sum():+.1f}",
                    f"{r['pt']['absorbed_u']:.0f} [{bu[0]:.0f}, {bu[1]:.0f}]",
                    f"{r['pt']['absorbed_v']:.0f} [{bv[0]:.0f}, {bv[1]:.0f}]",
                    " ".join(f"{sh[h][0]:.0f}/{sh[h][1]:.0f}" for h in FFFF_H),
                    " ".join(f"{sh[h][2]:.0f}" for h in FFFF_H)])
    say(tabulate(tab, ["spec", "F", "δ→u peak", "A peak", "cum / cum_A", "absorbed u",
                       "absorbed v", "FEVD δ 4q 8q", "FEVD F"], tablefmt="simple", disable_numparse=True))

    # 3. F -> u through delta
    say(f"\n{'-'*78}\n3. How much of a financial shock's effect on u runs through delta?\n{'-'*78}")
    say("  F -> u: peak and cum (h = 0–12); with delta held at zero; share via delta = 1 − cf/full, % [68 %].\n")
    tab = []
    for lab, r in runs:
        full, cf = r["pt"]["F_u"], r["pt"]["F_u_cf"]
        k = int(np.abs(full).argmax())
        b = band(r["draws"], "via_delta_u")
        tab.append([lab, r["F"], f"{100 * full[k]:+.2f} (h {k})", f"{100 * full.sum():+.1f}",
                    f"{100 * cf.sum():+.1f}", f"{r['pt']['via_delta_u']:.0f} [{b[0]:.0f}, {b[1]:.0f}]"])
    say(tabulate(tab, ["spec", "F", "F→u peak", "cum full", "cum δ fixed", "share via δ"],
                 tablefmt="simple", disable_numparse=True))

    # figure (HP, BIC lag, EBP)
    r = runs[0][1]
    names = r["names"]
    hh = np.arange(H + 1)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    lo, hi = 100 * np.percentile([d["F_delta"] for d in r["draws"]], [16, 84], axis=0)
    axes[0].fill_between(hh, lo, hi, color="#2a78d6", alpha=0.15, lw=0)
    axes[0].plot(hh, 100 * r["pt"]["F_delta"], color="#2a78d6", lw=2, label="EBP shock → δ")
    axes[0].set_title("EBP shock → δ (log points × 100), 68 % band", fontsize=10)
    jd, iu = names.index("delta"), names.index("u")
    axes[1].plot(hh, 100 * r["ThA"][:, NAMES_A.index("u"), NAMES_A.index("delta")], color="#52514e",
                 lw=2, label="SVAR A")
    axes[1].plot(hh, 100 * r["Th"][:, iu, jd], color="#2a78d6", lw=2, ls="--", label="EBP ordered first")
    axes[1].set_title("δ shock → u (log points × 100)", fontsize=10)
    axes[2].plot(hh, 100 * r["pt"]["F_u"], color="#52514e", lw=2, label="full")
    axes[2].plot(hh, 100 * r["pt"]["F_u_cf"], color="#2a78d6", lw=2, ls="--", label="δ held at zero")
    axes[2].set_title("EBP shock → u (log points × 100)", fontsize=10)
    for ax in axes:
        ax.axhline(0, color="#8a8984", lw=1)
        ax.set_xlabel("quarters")
        ax.grid(alpha=0.15)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"SVAR (z, EBP, s, δ, u, v), {runs[0][0]}", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_PNG}")
