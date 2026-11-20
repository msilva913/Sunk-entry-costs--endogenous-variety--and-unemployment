"""
q3_delta_cu.py -- is the delta shock contaminated by aggregate demand? Capacity utilization as the proxy.

Demand moves the share of existing capacity in use (Silva and Urias 2025, Economic Modelling 153,
107317; Bai, Rios-Rull and Storesletten 2025, Review of Economic Studies). Total capacity
utilization (CU; Federal Reserve G.17, total index, part7e_tcu.py) is the proxy, log HP cycle over
its own window (1967Q1-2019Q4, cut before filtering). It covers manufacturing, mining and
utilities only; BED deaths are mostly services. Monetary policy works through it too (lower real
rates raise spending and utilization), so it proxies demand broadly, not search effort alone.

Two orderings, as in the finance tests (q3_finance_delta.py, q3_delta_finance.py):
  1. (z, CU, s, delta, u, v): the delta shock net of same-quarter and lagged utilization. An upper
     bound on demand contamination: CU first takes all same-quarter comovement with delta,
     including any mechanical effect of closures on utilization (closed plants leave capacity).
     Reported: the CU-delta innovation correlation, CU -> delta, delta -> u and v against SVAR A
     (share absorbed, joint bootstrap), delta's FEVD share of u, and the share of a CU shock's
     effect on u that runs through delta.
  2. (z, s, delta, CU, u, v): the signature. A delta shock that is demand in disguise lowers
     utilization; destruction of under-used plants need not. Reported: delta -> CU, and the share
     of delta -> u, v running through CU.
Specifications: HP, BIC lag and p = 4; Hamilton, BIC lag.

Outputs:
    data/results/q3_delta_cu.txt
    data/results/q3_delta_cu.png
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b
import q3_delta_finance as fwd                  # (z, s, delta, F, u, v)
import q3_finance_delta as rev                  # (z, F, s, delta, u, v)

CU_CACHE = Path("data/cache/tcu_quarterly.parquet")
H        = q3.H_IRF
OUT_TXT  = p6b.RESULTS_DIR / "q3_delta_cu.txt"
OUT_PNG  = p6b.RESULTS_DIR / "q3_delta_cu.png"
LINES = []
fwd.UNIT["cu"] = ("log pts × 100", 100)


def say(text=""):
    LINES.append(text)
    print(text)


def cycles(method):
    """The finance tests' cycles plus CU (log cycle over its own window)."""
    cyc = fwd.cycles(method)
    tcu = pd.read_parquet(CU_CACHE)
    tcu.index = pd.PeriodIndex(tcu["quarter_label"], freq="Q")
    x = np.log(p6b.window(tcu["tcu"], start=tcu.index[0]))
    if method == "hp":
        cyc["cu"] = p6b.hpfilter(x, lamb=p6b.HP_LAMBDA)[0]
    else:
        cyc["cu"] = p6b.hamilton_cycle(x)
    return cyc


def specs(cyc_hp, cyc_ham, fit):
    """[(label, result)] for the BIC lag, p = 4 (if BIC differs), and Hamilton at its BIC lag."""
    out = []
    r = fit(cyc_hp, "cu")
    out.append((f"HP, p = {r['p']} (BIC)", r, cyc_hp))
    if r["p"] != 4:
        out.append(("HP, p = 4", fit(cyc_hp, "cu", p=4), cyc_hp))
    r = fit(cyc_ham, "cu")
    out.append((f"Hamilton, p = {r['p']} (BIC)", r, cyc_ham))
    return out


if __name__ == "__main__":
    say("=" * 78)
    say("Is the delta shock demand in disguise? Capacity utilization (q3_delta_cu.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])
    cyc_hp, cyc_ham = cycles("hp"), cycles("hamilton")

    # ── 1. CU first ─────────────────────────────────────────────────────────
    runs = specs(cyc_hp, cyc_ham, rev.run)
    fr, fh = runs[0][1]["frame"], runs[-1][1]["frame"]
    say(f"\nHP sample {fr.index[0]}-{fr.index[-1]}, T = {len(fr)}; Hamilton {fh.index[0]}-{fh.index[-1]}, "
        f"T = {len(fh)}. Bands: residual bootstrap, {q3.N_BOOT} draws.")
    say(f"\n{'-'*78}\n1. CU ordered before delta: (z, CU, s, delta, u, v)\n{'-'*78}")
    say("  corr: CU and delta innovations. CU→δ: one-SD CU shock, peak (quarter) [68 %].")
    say("  δ→u: peak in this SVAR and in SVAR A (same window and lag); absorbed = 1 − cum/cum_A over")
    say("  h = 0–12, % [68 %]; FEVD: δ's share of var(u) at 8q, this SVAR / A; via δ: share of a CU")
    say("  shock's cumulative effect on u running through δ, % [68 %].\n")
    tab = []
    for lab, r, _ in runs:
        names = r["names"]
        jd, iu, iv = names.index("delta"), names.index("u"), names.index("v")
        pF = r["pt"]["F_delta"]
        kF = int(np.abs(pF).argmax())
        bF = 100 * rev.band(r["draws"], "F_delta", kF)
        pr = r["Th"][:, iu, jd]
        pa = r["ThA"][:, rev.NAMES_A.index("u"), rev.NAMES_A.index("delta")]
        kr, ka = int(np.abs(pr).argmax()), int(np.abs(pa).argmax())
        bu, bv, bm = (rev.band(r["draws"], k) for k in ("absorbed_u", "absorbed_v", "via_delta_u"))
        tab.append([lab, f"{r['rho']:+.2f}", f"{100 * pF[kF]:+.2f} (h {kF}) [{bF[0]:+.2f}, {bF[1]:+.2f}]",
                    f"{100 * pr[kr]:+.2f} (h {kr})", f"{100 * pa[ka]:+.2f} (h {ka})",
                    f"{r['pt']['absorbed_u']:.0f} [{bu[0]:.0f}, {bu[1]:.0f}]",
                    f"{r['pt']['absorbed_v']:.0f} [{bv[0]:.0f}, {bv[1]:.0f}]",
                    f"{r['share'][8][0]:.0f} / {r['share'][8][1]:.0f}",
                    f"{r['pt']['via_delta_u']:.0f} [{bm[0]:.0f}, {bm[1]:.0f}]"])
    say(tabulate(tab, ["spec", "corr", "CU→δ peak", "δ→u peak", "A peak", "absorbed u", "absorbed v",
                       "FEVD δ 8q", "via δ"], tablefmt="simple", disable_numparse=True))
    say("\n  δ → u by quarter, HP at the BIC lag: SVAR A and net of CU (log points × 100):")
    r = runs[0][1]
    names = r["names"]
    pr = 100 * r["Th"][:, names.index("u"), names.index("delta")]
    pa = 100 * r["ThA"][:, rev.NAMES_A.index("u"), rev.NAMES_A.index("delta")]
    say("    h      " + " ".join(f"{h:6d}" for h in range(H + 1)))
    say("    A      " + " ".join(f"{x:+6.2f}" for x in pa))
    say("    net CU " + " ".join(f"{x:+6.2f}" for x in pr))

    # ── 2. CU after delta ──────────────────────────────────────────────────
    say(f"\n{'-'*78}\n2. CU ordered after delta: (z, s, delta, CU, u, v), the signature\n{'-'*78}")
    say("  δ→CU: one-SD δ shock, peak (quarter), 68 % and 90 % bands at that quarter, cum 0–12;")
    say("  δ→u, v peak; share of cum δ→u, v running through CU, % [68 %].\n")
    tab = []
    fwd_specs = specs(cyc_hp, cyc_ham, fwd.fit_spec)
    for lab, s, cyc in fwd_specs:
        f_row, lm_rows = fwd.rows_for(s, lab, cyc)
        tab.append([lab] + f_row[2:6] + [lm_rows[0][3], lm_rows[0][7], lm_rows[1][3], lm_rows[1][7]])
    say(tabulate(tab, ["spec", "δ→CU peak", "68 %", "90 %", "cum 0–12", "δ→u peak", "via CU u",
                       "δ→v peak", "via CU v"], tablefmt="simple", disable_numparse=True))
    s = fwd_specs[0][1]
    jd, iF = s["names"].index("delta"), s["names"].index("cu")
    lo, hi = 100 * np.percentile(s["draws"][:, :, iF, jd], [16, 84], axis=0)
    say(f"\n  δ → CU by quarter ({fwd_specs[0][0]}), log points × 100, point [68 %]:")
    for h in range(H + 1):
        say(f"    h {h:2d}: {100 * s['Theta'][h, iF, jd]:+6.2f}  [{lo[h]:+6.2f}, {hi[h]:+6.2f}]")

    # s -> CU: an unemployment-risk (demand) channel is blind to whether the job loss destroys the
    # position, so it predicts utilization falls after s at least as much as after delta
    # (Broer et al. 2025, IER: front-loaded job-loss risk drives precautionary saving).
    js = s["names"].index("s")
    lo_s, hi_s = 100 * np.percentile(s["draws"][:, :, iF, js], [16, 84], axis=0)
    say(f"\n  s → CU by quarter ({fwd_specs[0][0]}), log points × 100, point [68 %]:")
    for h in range(H + 1):
        say(f"    h {h:2d}: {100 * s['Theta'][h, iF, js]:+6.2f}  [{lo_s[h]:+6.2f}, {hi_s[h]:+6.2f}]")
    say(f"  cum 0–12: s → CU {100 * s['Theta'][:, iF, js].sum():+.2f}, "
        f"δ → CU {100 * s['Theta'][:, iF, jd].sum():+.2f}")

    # ── figure ──────────────────────────────────────────────────────────────
    hh = np.arange(H + 1)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(hh, pa, color="#52514e", lw=2, label="SVAR A")
    axes[0].plot(hh, pr, color="#2a78d6", lw=2, ls="--", label="CU ordered first")
    axes[0].set_title("δ shock → u (log points × 100)", fontsize=10)
    axes[1].fill_between(hh, lo, hi, color="#2a78d6", alpha=0.15, lw=0)
    axes[1].plot(hh, 100 * s["Theta"][:, iF, jd], color="#2a78d6", lw=2, label="δ shock → CU, 68 % band")
    axes[1].set_title("δ shock → CU (CU after δ)", fontsize=10)
    r = runs[0][1]
    axes[2].plot(hh, 100 * r["pt"]["F_u"], color="#52514e", lw=2, label="CU shock → u")
    axes[2].plot(hh, 100 * r["pt"]["F_u_cf"], color="#2a78d6", lw=2, ls="--", label="δ held at zero")
    axes[2].set_title("CU shock → u (CU first)", fontsize=10)
    for ax in axes:
        ax.axhline(0, color="#8a8984", lw=1)
        ax.set_xlabel("quarters")
        ax.grid(alpha=0.15)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"Capacity utilization and the δ shock, {runs[0][0]}", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_PNG}")
