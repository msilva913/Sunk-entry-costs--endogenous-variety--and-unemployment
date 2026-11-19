"""
q3_gross_flows_compare.py -- duration-based (Shimer) against gross-flow (CPS E->U, U->E) separation
and job-finding rates: the series, and the empirical impulse responses.

Series. Old: tau and jf from observables.py (shimer_flows); s = (tau - delta)/(1 - delta) with dbar
before BED (part6b build_levels). New: tau_gf = E->U and jf_gf = U->E probabilities (part6f), with
s built by the same identity and splice. Old series are filtered over 1951-2019 (as in every SVAR
so far); new ones over their own span, 1990Q2-2019Q4 (cut before filtering, S13). A third column
re-filters the old series over 1990Q2-2019Q4, to separate the measure from the filter span.

Impulse responses. SVAR (z, s, delta, u, v, jf), p = 1, one-SD shocks, HP and Hamilton, with old and
new (s, jf); 90 % bands from 1,000-draw residual bootstraps. The impact matching-efficiency residual
of a shock, jf - (1 - eta)(v - u) at h = 0 with eta = 0.6, is what a stable Cobb-Douglas matching
function requires to be zero (shared measurement, recall or non-participation inflows show up
here).

Outputs:
    data/results/q3_gross_flows_compare.txt
    data/results/q3_gross_flows_compare.csv
"""

import os
from pathlib import Path

# The "old" columns are the duration-based baseline: pin the pipeline switch (part6b, D19) to it.
os.environ["FLOW_MEASURE"] = "duration"

import numpy as np
import pandas as pd
from tabulate import tabulate

import q3_data_svar as q3
import part6b_shock_persistence_cyclical as p6b

NAMES   = ["z", "s", "delta", "u", "v", "jf"]
H       = q3.H_IRF
ETA     = 0.6
GF_Q    = Path("data/cache/gross_flows_quarterly.parquet")
GF_START = "1990Q2"
OUT_TXT = Path("data/results") / "q3_gross_flows_compare.txt"
OUT_CSV = Path("data/results") / "q3_gross_flows_compare.csv"
FILTERS = {"hp": ("hp", q3.BED_START, "HP"), "hamilton": ("hamilton", "1995Q2", "Hamilton")}
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def raw_levels():
    dat = pd.read_pickle(p6b.RAW_DATA)
    if not isinstance(dat.index, pd.PeriodIndex):
        dat.index = dat.index.to_period("Q")
    lev = p6b.build_levels(dat)                         # tau, delta, delta_used, s (old)
    gf = pd.read_parquet(GF_Q)
    gf.index = pd.PeriodIndex(gf.index, freq="Q")
    gf = gf[(gf.index >= GF_START) & (gf.index <= p6b.SAMPLE_END)]
    du = lev["delta_used"].reindex(gf.index)
    new = pd.DataFrame({"tau": gf["tau_gf"], "jf": gf["jf_gf"],
                        "s": (gf["tau_gf"] - du) / (1.0 - du)})
    old = pd.DataFrame({"tau": lev["tau"], "jf": dat["jf"], "s": lev["s"]})
    return old, new, lev


def cycles(method, version, old, new):
    """q3_data_svar cycles with s and jf replaced by the requested version."""
    cyc = q3.build_cycles(method)
    if version == "new":
        for c in ("s", "jf"):
            cyc[c] = p6b.log_cycle(p6b.window(new[c], start=GF_START), method).reindex(cyc.index)
    elif version == "old, 1990 span":
        for c in ("s", "jf"):
            cyc[c] = p6b.log_cycle(p6b.window(old[c], start=GF_START), method).reindex(cyc.index)
    return cyc


def me_residual(Theta, shock):
    j = NAMES.index(shock)
    jf, v, u = (Theta[0, NAMES.index(x), j] for x in ("jf", "v", "u"))
    return jf - (1 - ETA) * (v - u)


if __name__ == "__main__":
    old, new, lev = raw_levels()
    win = lambda x: x[(x.index >= q3.BED_START) & (x.index <= p6b.SAMPLE_END)]
    say("=" * 78)
    say("Duration-based against gross-flow s and jf (q3_gross_flows_compare.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])

    # ── 1. levels ──────────────────────────────────────────────────────────
    d_mean = win(lev["delta"]).mean()
    tab = []
    for c, lab in (("tau", "separation prob. tau (monthly)"), ("jf", "job-finding prob. (monthly)"),
                   ("s", "s = (tau - delta)/(1 - delta)")):
        tab.append([lab, f"{win(old[c]).mean():.4f}", f"{win(new[c]).mean():.4f}"])
    tab.append(["delta / tau", f"{d_mean / win(old['tau']).mean():.3f}", f"{d_mean / win(new['tau']).mean():.3f}"])
    say(f"\n{'-' * 78}\n1. Levels, means over 1992Q3-2019Q4\n{'-' * 78}")
    say(tabulate(tab, ["", "duration (old)", "gross flow (new)"], tablefmt="simple", disable_numparse=True))

    # ── 2. cycles ──────────────────────────────────────────────────────────
    say(f"\n{'-' * 78}\n2. HP log cycles on 1992Q3-2019Q4: SD (× 100), rho_1, correlations\n{'-' * 78}")
    cy = {v: cycles("hp", v, old, new) for v in ("old", "old, 1990 span", "new")}
    tab = []
    for c in ("s", "jf"):
        for v, C in cy.items():
            x = win(C[c].dropna())
            u = win(C["u"]).reindex(x.index)
            other = "jf" if c == "s" else "s"
            tab.append([c, v, f"{100 * x.std():.2f}", f"{x.autocorr():.2f}", f"{x.corr(u):+.2f}",
                        f"{x.corr(win(C[other]).reindex(x.index)):+.2f}",
                        f"{x.corr(win(cy['new'][c]).reindex(x.index)):+.2f}"])
    say(tabulate(tab, ["series", "version", "SD", "ρ1", "corr with u", "corr(s, jf)", "corr with new"],
                 tablefmt="simple", disable_numparse=True))

    # ── 3. impulse responses ───────────────────────────────────────────────
    rows = []
    for fkey, (method, start, flabel) in FILTERS.items():
        say(f"\n{'-' * 78}\n3. SVAR (z, s, δ, u, v, jf), p = 1, {flabel}: one-SD responses × 100 [90 %]\n{'-' * 78}")
        tab = []
        for v in ("old", "new"):
            C = cy[v] if fkey == "hp" else cycles(method, v, old, new)
            frame = C.loc[C.index >= start, NAMES].dropna()
            Y = frame.to_numpy()
            fit = q3.fit_var(Y, 1)
            Th = q3.irf(fit, NAMES, NAMES, H)
            D = q3.bootstrap_irf(Y, fit, NAMES, NAMES, H, q3.N_BOOT)
            lo, hi = np.percentile(D, 5, axis=0), np.percentile(D, 95, axis=0)
            for shock, resp, hs in (("s", "u", (0, 1, 2, 4, 8)), ("s", "v", (0, 2, 4)), ("s", "jf", (0, 2, 4)),
                                    ("delta", "u", (0, 2, 4, 8)), ("delta", "jf", (0, 4)), ("z", "u", (0, 4))):
                i, j = NAMES.index(resp), NAMES.index(shock)
                path = 100 * Th[:, i, j]
                k = int(np.abs(path).argmax())
                cells = [f"h{h} {100 * Th[h, i, j]:+.2f} [{100 * lo[h, i, j]:+.2f}, {100 * hi[h, i, j]:+.2f}]" for h in hs]
                tab.append([v, f"{q3.LABEL[shock]}→{resp}", f"{path[k]:+.2f} (h {k})", f"{path.sum():+.1f}", "; ".join(cells)])
                for h in range(H + 1):
                    rows.append({"filter": fkey, "version": v, "shock": shock, "response": resp, "h": h,
                                 "irf": Th[h, i, j], "lo90": lo[h, i, j], "hi90": hi[h, i, j]})
            me = me_residual(Th, "s")
            me_d = np.array([me_residual(Db, "s") for Db in D])
            tab.append([v, "s: matching-eff. residual h0", f"{100 * me:+.2f}", "",
                        f"[{100 * np.percentile(me_d, 5):+.2f}, {100 * np.percentile(me_d, 95):+.2f}]"])
            tab.append(["", "", "", "", ""])
            say(f"  {v}: {frame.index[0]}-{frame.index[-1]}, T = {len(frame)}")
        say(tabulate(tab, ["version", "response", "peak", "cum 0–12", "by horizon [90 %]"],
                     tablefmt="simple", disable_numparse=True))
    say(f"\n  Matching-efficiency residual: jf − (1 − η)(v − u) on impact, η = {ETA}; zero under a stable")
    say("  Cobb-Douglas matching function with no measurement error.")
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}")
