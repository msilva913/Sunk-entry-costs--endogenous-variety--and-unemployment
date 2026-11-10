"""
q3_delta_finance.py -- does a delta shock tighten financial conditions?

The financial friction proposed for AGS (Notes/ags_equilibrium.md section 10) makes the
death of product lines destroy the lines' internal funds. Leverage and the external finance
premium then rise, and stay high while funds are rebuilt. Two predictions:
  (a) after a delta shock, credit spreads rise, with a lag (funds rebuild slowly when moving
      funds between households and lines is costly, kappa_d > 0);
  (b) a sizable share of delta's effect on u and v runs through financial conditions.
If spreads do not respond, the friction explains at most the shape of delta's effect, not
its size, and the data point to kappa_d near zero (section 10.3: a level wedge only).

SVAR D_F: SVAR A with a financial variable F ordered after delta and before (u, v):
    (z, s, delta, F, u, v), recursive, 1992Q3-2019Q4.
delta is then not allowed to respond to F within the quarter (it responds with a lag), and
F responds to delta within the quarter. This is the reverse of A_F in q3_data_svar.py
(F before delta), which asks how much of delta is a financial shock. Here the question is
how much financial conditions respond to delta.
    F = ebp        Gilchrist-Zakrajsek excess bond premium (the headline; part7d_ebp.py)
        gz_spread  the GZ credit spread (EBP plus expected default)
        nfci       Chicago Fed NFCI (part7_nfci.py)
        anfci      NFCI adjusted for macro conditions
F is in percentage points (EBP, GZ) or index units (NFCI), HP-filtered in levels over its own
window (cut at 2019Q4 first), like the NFCI in q3_data_svar.py.

Reported, for the BIC lag, p = 4, and Hamilton cycles (EBP only):
  1. delta -> F: the response to a one-SD delta shock (basis points for EBP and GZ), its peak
     and quarter, 68 % and 90 % residual-bootstrap bands at the peak.
  2. delta -> u, v: peak and cumulative response over h = 0-12, against SVAR A (no F).
  3. Mediation: the share of the cumulative delta -> u and delta -> v response (h = 0-12) that
     runs through F, holding F's response at zero with offsetting F shocks (q3_data_svar.mediated;
     a reduced-form counterfactual), with a 68 % band.

Outputs:
    data/results/q3_delta_finance.txt
    data/results/q3_delta_finance.png
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

EBP_CACHE = Path("data/cache/ebp_quarterly.parquet")
FIN_VARS  = ("ebp", "gz_spread", "nfci", "anfci")
UNIT      = {"ebp": ("bp", 100), "gz_spread": ("bp", 100), "nfci": ("index × 100", 100),
             "anfci": ("index × 100", 100)}
H         = q3.H_IRF
N_BOOT    = q3.N_BOOT
OUT_TXT   = p6b.RESULTS_DIR / "q3_delta_finance.txt"
OUT_PNG   = p6b.RESULTS_DIR / "q3_delta_finance.png"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def cycles(method):
    """q3_data_svar's cycles plus the EBP and GZ spread (levels, filtered over their window)."""
    cyc = q3.build_cycles(method)
    ebp = pd.read_parquet(EBP_CACHE)
    ebp.index = pd.PeriodIndex(ebp["quarter_label"], freq="Q")
    for name in ("ebp", "gz_spread"):
        x = p6b.window(ebp[name], start=ebp.index[0])
        if method == "hp":
            cyc[name] = p6b.hpfilter(x, lamb=p6b.HP_LAMBDA)[0]
        else:
            cyc[name] = p6b.hamilton_cycle(x)
    return cyc


def names_for(F):
    return ["z", "s", "delta", F, "u", "v"]


def fit_spec(cyc, F, p=None):
    """SVAR D_F on `cyc`: frame, lag, point responses and bootstrap draws."""
    names = names_for(F)
    frame = cyc.loc[cyc.index >= q3.BED_START, names].dropna()
    Y = frame.to_numpy()
    if p is None:
        p = q3.bic_lag(Y)[0]
    fit = q3.fit_var(Y, p)
    Theta = q3.irf(fit, names, names, H)
    draws = q3.bootstrap_irf(Y, fit, names, names, H, N_BOOT)
    return {"F": F, "names": names, "frame": frame, "p": p, "Theta": Theta, "draws": draws}


def svar_a(cyc, frame_index, p):
    """SVAR A (no F) on the same window and lag, for comparison."""
    names = list(q3.SVARS["A"]["vars"])
    frame = cyc.loc[frame_index, names]
    return q3.irf(q3.fit_var(frame.to_numpy(), p), names, names, H), names


def share_through(Theta, names, resp, F):
    """Share (%) of the cumulative delta -> resp response (h = 0-12) that runs through F."""
    j, m, i = names.index("delta"), names.index(F), names.index(resp)
    full = Theta[:, i, j]
    cf = q3.mediated(Theta, j, m)[:, i]
    return 100 * (1 - cf.sum() / full.sum()), full, cf


def peak(path):
    k = int(np.abs(path).argmax())
    return k, path[k]


def rows_for(spec, label, cyc):
    """Table rows: delta -> F, delta -> u and v (with SVAR A), mediation shares."""
    F, names, Th, D = spec["F"], spec["names"], spec["Theta"], spec["draws"]
    unit, scale = UNIT[F]
    j, iF = names.index("delta"), names.index(F)
    k, pk = peak(Th[:, iF, j])
    band68 = np.percentile(D[:, k, iF, j], [16, 84]) * scale
    band90 = np.percentile(D[:, k, iF, j], [5, 95]) * scale
    f_row = [label, F, f"{scale * pk:+.1f} (h {k})", f"[{band68[0]:+.1f}, {band68[1]:+.1f}]",
             f"[{band90[0]:+.1f}, {band90[1]:+.1f}]", f"{scale * Th[:, iF, j].sum():+.1f}", unit]
    ThA, namesA = svar_a(cyc, spec["frame"].index, spec["p"])
    lm_rows = []
    for resp in ("u", "v"):
        sh, full, cf = share_through(Th, names, resp, F)
        bs = [share_through(Dd, names, resp, F)[0] for Dd in D]
        kk, pv = peak(full)
        a = ThA[:, namesA.index(resp), namesA.index("delta")]
        ka, pa = peak(a)
        lm_rows.append([label, F, resp, f"{100 * pv:+.2f} (h {kk})", f"{100 * pa:+.2f} (h {ka})",
                        f"{100 * full.sum():+.1f}", f"{100 * cf.sum():+.1f}",
                        f"{sh:.0f} [{np.percentile(bs, 16):.0f}, {np.percentile(bs, 84):.0f}]"])
    return f_row, lm_rows


if __name__ == "__main__":
    say("=" * 78)
    say("Does a delta shock tighten financial conditions? (q3_delta_finance.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])

    cyc_hp = cycles("hp")
    cyc_ham = cycles("hamilton")
    specs = []
    for F in FIN_VARS:
        s_bic = fit_spec(cyc_hp, F)
        specs.append((f"HP, p = {s_bic['p']} (BIC)", s_bic, cyc_hp))
        if s_bic["p"] != 4:
            specs.append(("HP, p = 4", fit_spec(cyc_hp, F, p=4), cyc_hp))
    s_ham = fit_spec(cyc_ham, "ebp")
    specs.append((f"Hamilton, p = {s_ham['p']} (BIC)", s_ham, cyc_ham))

    f_tab, lm_tab = [], []
    for label, spec, cyc in specs:
        f_row, lm_rows = rows_for(spec, label, cyc)
        f_tab.append(f_row)
        lm_tab += lm_rows
    fr = specs[0][1]["frame"]
    fh = specs[-1][1]["frame"]
    say(f"\nHP sample {fr.index[0]}-{fr.index[-1]}, T = {len(fr)}; Hamilton {fh.index[0]}-{fh.index[-1]}, "
        f"T = {len(fh)}. Bands: residual bootstrap, {N_BOOT} draws.")

    say(f"\n{'-'*78}\n1. delta -> F: response to a one-SD delta shock\n{'-'*78}")
    say("  peak: largest |response| and its quarter; bands at the peak's quarter; cum: sum over h = 0–12.\n")
    say(tabulate(f_tab, ["spec", "F", "peak", "68 %", "90 %", "cum 0–12", "units"], tablefmt="simple",
                 disable_numparse=True))

    say(f"\n{'-'*78}\n2–3. delta -> u, v, and the share running through F\n{'-'*78}")
    say("  peak (log points × 100) in D_F and in SVAR A on the same window and lag; cum: sum over")
    say("  h = 0–12, full and with F held at zero; share through F = 1 − cf/full, % [68 % band].\n")
    say(tabulate(lm_tab, ["spec", "F", "", "peak D_F", "peak A", "cum full", "cum F fixed",
                          "share via F"], tablefmt="simple", disable_numparse=True))

    # path of delta -> EBP, HP at the BIC lag, for the record
    head = specs[0][1]
    j, iF = head["names"].index("delta"), head["names"].index("ebp")
    lo, hi = np.percentile(head["draws"][:, :, iF, j], [16, 84], axis=0)
    say(f"\n  delta -> EBP by quarter (bp; {specs[0][0]}), point [68 %]:")
    for h in range(H + 1):
        say(f"    h {h:2d}: {100 * head['Theta'][h, iF, j]:+6.1f}  [{100 * lo[h]:+6.1f}, {100 * hi[h]:+6.1f}]")

    # figure: delta -> F (HP, BIC lag) for each F; delta -> u and v, full and with EBP fixed
    hh = np.arange(H + 1)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    pal = {"ebp": "#2a78d6", "gz_spread": "#eb6834", "nfci": "#1baf7a", "anfci": "#eda100"}
    for label, spec, _ in specs:
        if not label.startswith("HP") or "BIC" not in label:
            continue
        F = spec["F"]
        jj, ii = spec["names"].index("delta"), spec["names"].index(F)
        path = 100 * spec["Theta"][:, ii, jj]
        axes[0].plot(hh, path, color=pal[F], lw=2, label=F)
        if F == "ebp":
            lo, hi = np.percentile(100 * spec["draws"][:, :, ii, jj], [16, 84], axis=0)
            axes[0].fill_between(hh, lo, hi, color=pal[F], alpha=0.15, lw=0)
    axes[0].set_title("δ shock → F (bp; NFCI index × 100), 68 % band for EBP", fontsize=10)
    for a, resp in enumerate(("u", "v"), start=1):
        _, full, cf = share_through(head["Theta"], head["names"], resp, "ebp")
        axes[a].plot(hh, 100 * full, color="#52514e", lw=2, label="full")
        axes[a].plot(hh, 100 * cf, color="#2a78d6", lw=2, ls="--", label="EBP held at zero")
        axes[a].set_title(f"δ shock → {resp} (log points × 100)", fontsize=10)
    for ax in axes:
        ax.axhline(0, color="#8a8984", lw=1)
        ax.set_xlabel("quarters")
        ax.grid(alpha=0.15)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"SVAR (z, s, δ, F, u, v), {specs[0][0]}", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_PNG}")
