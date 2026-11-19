"""
q3_estimation_figure.py -- the estimated AGS model against the data through the SVAR (D18, D20).

For each filter (HP, Hamilton): a 3 x 3 grid, rows = shocks (z, s, delta), columns = responses
(u, v, f), horizons 0-12 quarters, one-SD shocks, log points x 100.
    data   SVAR (z, s, delta, u, v, jf) point estimate (estimation_targets.csv), with 68 and 90 %
           percentile bands from the 1,000-draw residual bootstrap (estimation_draws_<f>.npz)
    model  mean over 200 simulated windows at the final estimates
           (Programs AGS/ags_estimate_irf_<f>.csv, from ags_estimate.jl AGS_MODE = final)
The z and s rows are estimated; the delta row is held out (the test).

Outputs: data/results/ags_estimate_irf_hp.png, data/results/ags_estimate_irf_ham.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import part6b_shock_persistence_cyclical as p6b

AGS_DIR  = Path("../Programs AGS")
TARGETS  = p6b.RESULTS_DIR / "estimation_targets.csv"
SHOCKS   = ["z", "s", "delta"]
RESP     = ["u", "v", "jf"]
H        = 13
LAB_SHK  = {"z": r"Productivity $z$", "s": r"Separation $s$", "delta": r"Destruction $\delta$ (test)"}
LAB_RSP  = {"u": r"Unemployment $u$", "v": r"Vacancies $v$", "jf": r"Job finding $f$"}
INK, INK2, MODEL = "#0b0b0b", "#52514e", "#2a78d6"
BAND68, BAND90 = "#b9b8b2", "#deddd8"


def bands(fkey):
    """Percentile bands (16/84, 5/95) per shock, response, horizon, from the bootstrap draws."""
    d = np.load(p6b.RESULTS_DIR / f"estimation_draws_{fkey}.npz")
    draws = np.hstack([d["estimation"], d["test"]])            # moment order: shock, response, h
    pct = np.percentile(draws, [5, 16, 84, 95], axis=0)
    point = np.concatenate([d["point_estimation"], d["point_test"]])
    return pct, point


if __name__ == "__main__":
    tg = pd.read_csv(TARGETS)
    for fkey in ("hp", "ham"):
        mod = pd.read_csv(AGS_DIR / f"ags_estimate_irf_{fkey}.csv")
        pct, point = bands(fkey)
        data_pt = tg[(tg["filter"] == fkey) & (tg["block"] != "own")]
        gap = np.abs(point - data_pt["data"].to_numpy()).max()
        if gap > 1e-12:
            raise ValueError(f"bootstrap point and targets differ ({gap})")
        fig, axes = plt.subplots(3, 3, figsize=(9.0, 7.6), sharex=True)
        hs = np.arange(H)
        for a, s in enumerate(SHOCKS):
            for b, r in enumerate(RESP):
                ax = axes[a, b]
                k = (a * len(RESP) + b) * H
                sl = slice(k, k + H)
                ax.fill_between(hs, 100 * pct[0, sl], 100 * pct[3, sl], color=BAND90, lw=0, label="Data, 90 % band")
                ax.fill_between(hs, 100 * pct[1, sl], 100 * pct[2, sl], color=BAND68, lw=0, label="Data, 68 % band")
                ax.plot(hs, 100 * point[sl], color=INK, lw=2, label="Data")
                m = mod[(mod["shock"] == s) & (mod["response"] == r)].sort_values("h")
                ax.plot(hs, 100 * m["model"].to_numpy(), color=MODEL, lw=2, label="Model")
                ax.axhline(0, color=INK2, lw=0.6)
                ax.grid(axis="y", color="#e8e7e3", lw=0.5)
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
                for side in ("left", "bottom"):
                    ax.spines[side].set_color(INK2)
                ax.tick_params(colors=INK2, labelsize=8)
                ax.set_xticks([0, 4, 8, 12])
                if a == 0:
                    ax.set_title(LAB_RSP[r], fontsize=10, color=INK)
                if b == 0:
                    ax.set_ylabel(LAB_SHK[s], fontsize=10, color=INK)
                if a == 2:
                    ax.set_xlabel("Quarters", fontsize=9, color=INK2)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        order = [2, 3, 1, 0]                                   # Data, Model, 68 %, 90 %
        fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper center",
                   ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout(rect=(0, 0, 1, 0.955))
        out = p6b.RESULTS_DIR / f"ags_estimate_irf_{fkey}.png"
        fig.savefig(out, dpi=300, facecolor="white")
        plt.close(fig)
        print(f"Saved: {out}")
