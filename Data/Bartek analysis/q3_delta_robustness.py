"""
q3_delta_robustness.py -- is the data's persistent delta response a VAR artifact?

The delta shock's effect on u in SVAR A (HP, p = 1) builds for four quarters and runs through
depressed job finding for about two years (q3_persistence_tests.txt). Every estimate so far
extrapolates with the VAR's dynamics. Two checks that do not:

1. Local projections on the identified delta shock. The structural shocks of SVAR A (recursive,
   z, s, delta, u, v; the data's BIC lag) are recovered from the reduced-form residuals,
   eps_t = P^-1 u_t. For each horizon h = 0..12,
       y_{t+h} = a_h + b_h eps^delta_t + sum_{l=1..4} c_{h,l}' X_{t-l} + e_{t+h},
   with X = (z, s, delta, u, v) and y in {u, v, jf}; b_h is the response to a one-SD delta shock.
   HAC (Newey-West) standard errors with h + 1 lags. In population, a local projection on the
   identified shock with these controls estimates the same response as the VAR (Plagborg-Moller
   and Wolf 2021, citation unverified); in a finite sample it does not impose the VAR's decay.
2. The Hamilton filter (h = 8, p = 4; one-sided) in place of HP: SVAR A and the local projection
   on Hamilton cycles (q3_data_svar.build_cycles("hamilton")). delta's Hamilton cycle starts in
   1995Q2, so the sample is shorter.

Reported for each specification: the peak of u's response and its quarter, u at h = 8 relative to
the peak, the cumulative response over h = 0-12, and the same for v and job finding.

Outputs:
    data/results/q3_delta_robustness.txt   the report
    data/results/q3_delta_robustness.png   delta -> u, v, jf: SVAR and local projections, HP and Hamilton
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

NAMES    = list(q3.SVARS["A"]["vars"])           # z, s, delta, u, v
OUTS     = ("u", "v", "jf")
H        = q3.H_IRF
N_CTRL   = 4
OUT_TXT  = p6b.RESULTS_DIR / "q3_delta_robustness.txt"
OUT_PNG  = p6b.RESULTS_DIR / "q3_delta_robustness.png"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def svar_delta(cyc, p=None):
    """SVAR A on `cyc`: delta-shock responses of u, v (and jf with jf appended last), the
    structural delta shock series, and the frame. p = None uses the BIC lag."""
    sv = q3.prepare_svar("A", q3.SVARS["A"], cyc)
    if p is None:
        p = sv["p_bic"]
    fit = q3.fit_var(sv["Y"], p)
    P = q3.structural(fit, NAMES, NAMES)
    eps = fit["U"] @ np.linalg.inv(P).T                     # structural shocks, NAMES order
    shock = pd.Series(eps[:, NAMES.index("delta")], index=sv["frame"].index[p:])
    shock = shock / shock.std()
    Th = q3.irf(fit, NAMES, NAMES, H)
    resp = {r: Th[:, NAMES.index(r), NAMES.index("delta")] for r in ("u", "v")}
    # jf: the same SVAR with jf appended last (the delta shock's identification is unchanged)
    names_f = NAMES + ["jf"]
    frame_f = cyc.loc[cyc.index >= sv["frame"].index[0], names_f].dropna()
    Th_f = q3.irf(q3.fit_var(frame_f.to_numpy(), p), names_f, names_f, H)
    resp["jf"] = Th_f[:, names_f.index("jf"), names_f.index("delta")]
    return resp, shock, sv["frame"], p


def local_projection(cyc, shock, frame):
    """b_h and HAC standard errors for each outcome and h = 0..H."""
    data = cyc.loc[frame.index, NAMES + ["jf"]].copy()
    data["eps"] = shock
    ctrl = pd.concat({f"{v}_l{l}": data[v].shift(l) for v in NAMES for l in range(1, N_CTRL + 1)}, axis=1)
    out = {}
    for y in OUTS:
        b, se = np.full(H + 1, np.nan), np.full(H + 1, np.nan)
        for h in range(H + 1):
            df = pd.concat([data[y].shift(-h).rename("y"), data["eps"], ctrl], axis=1).dropna()
            fit = OLS(df["y"], add_constant(df.drop(columns="y"))).fit(cov_type="HAC",
                                                                      cov_kwds={"maxlags": h + 1})
            b[h], se[h] = fit.params["eps"], fit.bse["eps"]
        out[y] = (b, se)
    return out


def summary(path):
    k = int(np.abs(path).argmax())
    return f"{100 * path[k]:+.2f} (h {k})", f"{path[8] / path[k]:+.2f}", f"{100 * np.nansum(path):+.1f}"


if __name__ == "__main__":
    say("=" * 78)
    say("Is the data's persistent delta response a VAR artifact? (q3_delta_robustness.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])

    specs, lps = {}, {}
    for filt in ("hp", "hamilton"):
        cyc = q3.build_cycles(filt)
        for p in (None, 4):
            resp, shock, frame, p_used = svar_delta(cyc, p)
            specs[(filt, f"SVAR p = {p_used}")] = resp
            if p is None:
                lps[filt] = (local_projection(cyc, shock, frame), frame, p_used)
                say(f"\n  {filt.upper()}: SVAR A {frame.index[0]}-{frame.index[-1]}, T = {len(frame)}, "
                    f"BIC p = {p_used}; local projection on the p = {p_used} shock, {N_CTRL} lags of controls.")

    say(f"\n{'-'*78}\nResponses to a one-SD delta shock (log points × 100)\n{'-'*78}")
    say("  peak: largest |response| and its quarter; h8/peak; cum: sum over h = 0–12.\n")
    for y in OUTS:
        tab = []
        for filt in ("hp", "hamilton"):
            for (f, lab), resp in specs.items():
                if f == filt:
                    tab.append([filt, lab] + list(summary(resp[y])))
            b, se = lps[filt][0][y]
            k = int(np.abs(b).argmax())
            pk, r8, cum = summary(b)
            tab.append([filt, "local projection", f"{pk} [± {100 * se[k]:.2f}]", r8, cum])
        say(f"  {y}")
        say(tabulate(tab, ["filter", "estimator", "peak", "h8/peak", "cum 0–12"], tablefmt="simple",
                     disable_numparse=True))
        say("")
    say("  Local projection: [± one HAC standard error at the peak].")

    # figure
    hh = np.arange(H + 1)
    col = {"hp": "#2a78d6", "hamilton": "#eb6834"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for a, y in enumerate(OUTS):
        ax = axes[a]
        for filt in ("hp", "hamilton"):
            b, se = lps[filt][0][y]
            ax.fill_between(hh, 100 * (b - se), 100 * (b + se), color=col[filt], alpha=0.15, lw=0)
            ax.plot(hh, 100 * b, color=col[filt], lw=2, label=f"{filt}: local projection ± 1 SE")
            for (f, lab), resp in specs.items():
                if f == filt:
                    ls = "--" if lab.endswith("= 4") else ":"
                    ax.plot(hh, 100 * resp[y], color=col[filt], lw=1.5, ls=ls, label=f"{filt}: {lab}")
        ax.axhline(0, color="#8a8984", lw=1)
        ax.set_title(f"δ shock → {y}", fontsize=10)
        ax.set_xlabel("quarters")
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("log points × 100")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("The data's δ response: SVAR A against local projections, HP and Hamilton", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_PNG}")
