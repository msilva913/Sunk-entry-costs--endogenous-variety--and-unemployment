"""
q3_persistence_tests.py -- two checks on the persistence gap between the data and AGS.

Through SVAR A on both sides (q3_data_svar.py's functions; the model's 200 simulated windows of
110 quarters from Programs AGS/ags_svar_sim.jl, headline cases "default" and "compromise").

Test 1. Lag length. In the data, every shock's effect on u is persistent (rho_1 about 0.95). With
p = 1, all impulse responses are combinations of the powers of one matrix, so they decay at the
same roots; those roots reflect the persistence of the data as a whole, including the shocks the
model lacks. Re-estimate SVAR A with p = 1, 2, 4 on data and model and compare, per shock: the
peak of u's response and its quarter, u at h = 8 relative to the peak, the cumulative response
over h = 0-12, and the conditional rho_1(u) and cor(u, v) (sums over h = 0-200, as in
q3_model_svar.py section 5). If more lags shorten the data's responses, part of the common
persistence is the VAR(1).

Test 2. Inflows against outflows. Add the job-finding rate f (the data's jf; the model's gross
f = (1 - delta) A theta^(1-eta)) last: SVAR (z, s, delta, u, v, f), p = 1. The z, s and delta
shocks keep their identification (they are ordered first). The steady-state flow identity
u = tau / (tau + f) gives, in log deviations,
    d ln u ~ (1 - u_bar) (d ln tau - d ln f),   d ln tau ~ w_d d ln delta + (1 - w_d) d ln s,
with w_d = delta/tau = 0.087 (D1) and u_bar = 0.06. Each shock's u response is split into an
inflow part (s and delta) and an outflow part (f); their sum is compared with u's response (the
identity holds approximately: u adjusts within a quarter or two, Shimer 2012, citation
unverified). If u stays high after delta and s because f stays low, the missing mechanism is on
the job-creation side.

Outputs:
    data/results/q3_persistence_tests.txt   the report
    data/results/q3_persistence_tests.png   u's response to each shock at p = 1 and p = 4, data
                                            and model; and test 2's decomposition
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

MODEL_DIR = Path("../Programs AGS")
CASES     = ("default", "compromise")
NAMES     = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
NAMES_F   = NAMES + ["jf"]                       # test 2
SHOCKS    = q3.SHOCKS                            # z, delta, s
LAGS      = (1, 2, 4)
H, H_MOM  = q3.H_IRF, 200
W_DELTA, U_BAR = 0.087, 0.06
OUT_TXT   = p6b.RESULTS_DIR / "q3_persistence_tests.txt"
OUT_PNG   = p6b.RESULTS_DIR / "q3_persistence_tests.png"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def irfs(Y, p, names, Hh):
    """Recursive SVAR in `names` order on the data matrix Y; (Hh+1, var, shock)."""
    return q3.irf(q3.fit_var(Y, p), names, names, Hh)


def persistence_stats(Theta_short, Theta_long, names, shock):
    """u's response to `shock`: peak (value, quarter), u(h=8)/peak, cumulative 0-12, and the
    conditional rho_1(u) and cor(u, v) from the long-horizon responses."""
    j, iu, iv = names.index(shock), names.index("u"), names.index("v")
    path = Theta_short[:, iu, j]
    k = int(np.abs(path).argmax())
    U, V = Theta_long[:, iu, j], Theta_long[:, iv, j]
    return {"peak": 100 * path[k], "h_peak": k, "ratio8": path[8] / path[k],
            "cum": 100 * path.sum(), "rho1": (U[:-1] * U[1:]).sum() / (U ** 2).sum(),
            "cor_uv": (U * V).sum() / np.sqrt((U ** 2).sum() * (V ** 2).sum())}


def decomposition(Theta, names, shock):
    """Inflow and outflow parts of u's response to `shock` (steady-state flow identity)."""
    j = names.index(shock)
    tau = W_DELTA * Theta[:, names.index("delta"), j] + (1 - W_DELTA) * Theta[:, names.index("s"), j]
    inflow = (1 - U_BAR) * tau
    outflow = -(1 - U_BAR) * Theta[:, names.index("jf"), j]
    return inflow, outflow, Theta[:, names.index("u"), j]


if __name__ == "__main__":
    cyc = q3.build_cycles("hp")
    svA = q3.prepare_svar("A", q3.SVARS["A"], cyc)
    frame = svA["frame"]
    Yd = frame[NAMES].to_numpy()
    reps = {c: pd.read_csv(MODEL_DIR / f"ags_svar_reps_{c}.csv") for c in CASES}

    say("=" * 78)
    say("Persistence tests: lag length, and inflows against outflows (q3_persistence_tests.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])
    say(f"\nData: SVAR A, {frame.index[0]}-{frame.index[-1]}, T = {len(frame)}. Model: mean over "
        f"{reps['default'].rep.nunique()} windows of 110 quarters; σ = 1, Block B′ shock processes.")

    # ── Test 1: lag length ──────────────────────────────────────────────────
    say(f"\n{'-'*78}\nTest 1. Lag length: u's response to each shock, p = {', '.join(map(str, LAGS))}\n{'-'*78}")
    say("  peak: largest |response| (one SD, log points × 100) and its quarter; h8/peak: u at h = 8")
    say("  relative to the peak; cum: sum over h = 0–12; ρ1, cor(u,v): conditional on the shock.\n")
    stats = {}
    for p in LAGS:
        stats[("data", p)] = {s: persistence_stats(irfs(Yd, p, NAMES, H), irfs(Yd, p, NAMES, H_MOM), NAMES, s)
                              for s in SHOCKS}
        for c in CASES:
            per = {s: [] for s in SHOCKS}
            for _, g in reps[c].groupby("rep"):
                Y = g[NAMES].to_numpy()
                Ts, Tl = irfs(Y, p, NAMES, H), irfs(Y, p, NAMES, H_MOM)
                for s in SHOCKS:
                    per[s].append(persistence_stats(Ts, Tl, NAMES, s))
            stats[(c, p)] = {s: {k: np.mean([x[k] for x in per[s]]) for k in per[s][0]} for s in SHOCKS}
    for s in SHOCKS:
        say(f"\n  {q3.LABEL[s]} shock")
        tab = []
        for who in ("data",) + CASES:
            for p in LAGS:
                r = stats[(who, p)][s]
                tab.append([who, p, f"{r['peak']:+.2f} (h {r['h_peak']:.0f})", f"{r['ratio8']:+.2f}",
                            f"{r['cum']:+.1f}", f"{r['rho1']:.3f}", f"{r['cor_uv']:+.3f}"])
        say(tabulate(tab, ["", "p", "u peak", "h8/peak", "cum 0–12", "ρ1(u)", "cor(u,v)"],
                     tablefmt="simple", disable_numparse=True))
    say("\n  Model: the peak's quarter is the mean over windows.")

    # ── Test 2: inflows against outflows ────────────────────────────────────
    say(f"\n{'-'*78}\nTest 2. Inflows against outflows: SVAR (z, s, δ, u, v, f), p = 1\n{'-'*78}")
    say(f"  u ≈ (1 − ū)(τ − f) in log deviations, τ = {W_DELTA}·δ + {1 - W_DELTA:.3f}·s, ū = {U_BAR}.")
    say("  Inflow = (1 − ū)τ, outflow = −(1 − ū)f; each shock's u response, log points × 100.\n")
    framef = cyc.loc[cyc.index >= q3.BED_START, NAMES_F].dropna()
    Th_f = {"data": irfs(framef[NAMES_F].to_numpy(), 1, NAMES_F, H)}
    for c in CASES:
        Th_f[c] = np.mean([irfs(g[NAMES_F].to_numpy(), 1, NAMES_F, H) for _, g in reps[c].groupby("rep")], axis=0)
    dec = {}
    for s in SHOCKS:
        say(f"  {q3.LABEL[s]} shock")
        tab = []
        for who in ("data",) + CASES:
            inflow, outflow, u = decomposition(Th_f[who], NAMES_F, s)
            dec[(who, s)] = (inflow, outflow, u)
            for h in (0, 2, 4, 8, 12):
                tab.append([who, h, f"{100 * u[h]:+.2f}", f"{100 * inflow[h]:+.2f}", f"{100 * outflow[h]:+.2f}",
                            f"{100 * (inflow[h] + outflow[h]):+.2f}"])
            tab.append([who, "cum", f"{100 * u.sum():+.1f}", f"{100 * inflow.sum():+.1f}",
                        f"{100 * outflow.sum():+.1f}", f"{100 * (inflow + outflow).sum():+.1f}"])
        say(tabulate(tab, ["", "h", "u", "inflow", "outflow", "sum"], tablefmt="simple",
                     disable_numparse=True))
        say("")
    say(f"  Data sample for test 2: {framef.index[0]}-{framef.index[-1]}, T = {len(framef)}.")

    # ── figure ──────────────────────────────────────────────────────────────
    hh = np.arange(H + 1)
    pal = {"data": "#52514e", "default": "#2a78d6", "compromise": "#eb6834"}
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    for b, s in enumerate(SHOCKS):
        ax = axes[0, b]
        j, iu = NAMES.index(s), NAMES.index("u")
        for who in ("data",) + CASES:
            if who == "data":
                Y_list = [Yd]
            else:
                Y_list = [g[NAMES].to_numpy() for _, g in reps[who].groupby("rep")]
            for p, ls in ((1, "-"), (4, "--")):
                path = np.mean([irfs(Y, p, NAMES, H)[:, iu, j] for Y in Y_list], axis=0)
                ax.plot(hh, 100 * path, color=pal[who], lw=2, ls=ls, label=f"{who}, p = {p}")
        ax.axhline(0, color="#8a8984", lw=1, ls=":")
        ax.set_title(f"{q3.LABEL[s]} shock → u, p = 1 (solid) and 4 (dashed)", fontsize=10)
        ax.grid(alpha=0.15)
        ax = axes[1, b]
        for who in ("data", "default"):
            inflow, outflow, u = dec[(who, s)]
            ax.plot(hh, 100 * u, color=pal[who], lw=2, label=f"{who}: u")
            ax.plot(hh, 100 * inflow, color=pal[who], lw=1.5, ls="--", label=f"{who}: inflow")
            ax.plot(hh, 100 * outflow, color=pal[who], lw=1.5, ls=":", label=f"{who}: outflow")
        ax.axhline(0, color="#8a8984", lw=1, ls=":")
        ax.set_title(f"{q3.LABEL[s]} shock: u = inflow + outflow (p = 1)", fontsize=10)
        ax.set_xlabel("quarters")
        ax.grid(alpha=0.15)
    axes[0, 0].set_ylabel("log points × 100")
    axes[1, 0].set_ylabel("log points × 100")
    axes[0, 0].legend(fontsize=7, frameon=False)
    axes[1, 0].legend(fontsize=7, frameon=False)
    fig.suptitle("Persistence tests: lag length (top) and inflow/outflow decomposition (bottom)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_PNG}")
