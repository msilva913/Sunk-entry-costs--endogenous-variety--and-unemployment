"""
q3_model_svar.py -- Block B' model side: SVAR A on simulated AGS data, against the data.

The same estimator as q3_data_svar.py (imported, not copied): SVAR A on HP-1600 cycles of
(z, s, delta, u, v), recursive in that order, a constant, the data's BIC lag. Only the responses
to the z, s and delta shocks are used. The model has three shocks and the SVAR five variables, so
the model's u and v residuals are close to combinations of the first three; with u and v ordered
last, the Cholesky columns of z, s and delta stay well defined. Variance shares are not compared
(they need measurement error or residual shocks: the Block M decision).

Simulated data: Programs AGS/ags_svar_sim.jl, two cases:
    default      AGS_TARGETS (x_v 0.5, xi_inv 1)
    compromise   x_v 0.25, xi_inv 0.5 (Notes/ags_equilibrium.md section 7.6)
For each case:
    small sample  N_REP windows of T = 110 quarters (SVAR A's length), filtered as the data are:
                  z, s, u, v over the preceding 276 quarters (1951Q1-2019Q4 in the data), delta
                  over the 110 alone (window_filter in time_series_fun.jl). The model's IRF is the
                  mean over windows (indirect inference: model and data go through the same
                  small-sample estimator);
    long sample   one HP-filtered run of 30,000 quarters (the population SVAR);
    true          the model's quarterly-averaged IRFs to one-SD monthly shocks (validation).

Validation. In the model the shocks are close to observed (z and s exactly, delta with a one-month
lag), so a long-sample SVAR should recover the true responses. The SVAR shock is a one-SD
quarterly innovation and the true one a one-SD monthly shock, so both are normalized by the shock
variable's own response at h = 0 before they are compared. The true responses are not filtered,
so the clean test is the SVAR on the UNFILTERED long sample; it is run at p = 1 (the data's lag)
and at the long sample's BIC lag, which separates lag truncation from identification. The SVAR on
the HP-filtered long sample is shown too: its gap to the truth is the HP filter's effect, which
the data and the model share and which therefore does not bias the comparison.

Distance. J = sum over responses (u, v), shocks (z, s, delta) and h = 0..12 of
((model - data) / sd_data)^2, with sd_data the data bootstrap standard deviation of each
response (a diagonal weighting matrix, as in Christiano, Eichenbaum and Evans 2005; citation
unverified). Model = the small-sample mean. Reported in total and by shock, twice:
    one SD     responses to one-SD shocks, as estimated;
    scaled     responses divided by the shock variable's own impact response (per 1 % of z, of
               delta, of s on impact), in the data draw by draw. The one-SD distance also
               measures whether the SVAR sees shocks of the same size in model and data; the
               scaled one measures propagation alone.

A second target for the delta shock (section 4). In the data, delta is also estimated net of
financial conditions: A_F = (z, NFCI, s, delta, u, v), q3_data_svar.py section 4b. The model has
no NFCI, so A_F cannot be run on simulated data; the model's SVAR A delta responses are compared
with the data's A_F delta responses instead. This breaks the same-estimator rule of Block B'
(the data side conditions on NFCI, the model side does not), so A_F is a bound, not the target:
NFCI also tightens when the economy weakens, so conditioning on it removes some true delta
effect as well as financial contamination. A and A_F bracket the data's delta responses.

Outputs:
    data/results/q3_model_svar.txt   the report
    data/results/q3_model_svar.csv   IRFs: data (point, 16/84 %), model small-sample mean and
                                     16/84 %, long sample, true (normalized), by case; for the
                                     delta shock also data_AF (point, sd) from A_F
    data/results/q3_model_svar.png   u and v responses to z, s, delta: data band, model cases
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
CASE_LABEL = {"default": "default (x_v 0.5, ξ_inv 1)", "compromise": "compromise (x_v 0.25, ξ_inv 0.5)"}
NAMES     = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
SHOCKS    = q3.SHOCKS                            # z, delta, s
OUTCOMES  = q3.OUTCOMES                          # u, v
H         = q3.H_IRF
OUT_TXT   = p6b.RESULTS_DIR / "q3_model_svar.txt"
OUT_CSV   = p6b.RESULTS_DIR / "q3_model_svar.csv"
OUT_PNG   = p6b.RESULTS_DIR / "q3_model_svar.png"
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def svar_irf(frame, p):
    """SVAR A's IRFs (H+1, var, shock), recursive in NAMES order, on a frame of cycles."""
    Y = frame[NAMES].to_numpy()
    return q3.irf(q3.fit_var(Y, p), NAMES, NAMES, H)


def normalized(Theta):
    """Each shock's responses divided by the shock variable's own response at h = 0."""
    out = Theta.copy()
    for shk in SHOCKS:
        j = NAMES.index(shk)
        out[:, :, j] = Theta[:, :, j] / Theta[0, j, j]
    return out


if __name__ == "__main__":
    # ── data ────────────────────────────────────────────────────────────────
    cyc = q3.build_cycles("hp")
    svA = q3.prepare_svar("A", q3.SVARS["A"], cyc)
    p = svA["p_bic"]
    fitA = q3.fit_var(svA["Y"], p)
    Theta_d = q3.irf(fitA, NAMES, NAMES, H)
    draws = q3.bootstrap_irf(svA["Y"], fitA, NAMES, NAMES, H, q3.N_BOOT)
    sd_d = draws.std(axis=0)
    lo_d, hi_d = np.percentile(draws, 16, axis=0), np.percentile(draws, 84, axis=0)
    frame = svA["frame"]

    say("=" * 78)
    say("Q3 model SVAR — Block B' model side (q3_model_svar.py)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[0])
    say(f"\nData: SVAR A, {frame.index[0]}-{frame.index[-1]}, T = {len(frame)}, p = {p} (BIC), "
        f"HP-1600, ordering {', '.join(NAMES)}; {q3.N_BOOT} bootstrap draws.")

    # ── model ───────────────────────────────────────────────────────────────
    rows, res = [], {}
    for case in CASES:
        reps = pd.read_csv(MODEL_DIR / f"ags_svar_reps_{case}.csv")
        long = pd.read_csv(MODEL_DIR / f"ags_svar_long_{case}.csv")
        true = pd.read_csv(MODEL_DIR / f"ags_svar_true_{case}.csv")
        rep_irf = np.array([svar_irf(g, p) for _, g in reps.groupby("rep")])
        m_mean = rep_irf.mean(axis=0)
        m_lo, m_hi = np.percentile(rep_irf, 16, axis=0), np.percentile(rep_irf, 84, axis=0)
        raw = pd.read_csv(MODEL_DIR / f"ags_svar_raw_{case}.csv")
        p_long, _ = q3.bic_lag(long[NAMES].to_numpy())
        p_raw, _ = q3.bic_lag(raw[NAMES].to_numpy())
        Theta_long = svar_irf(long, p)
        Theta_raw1 = svar_irf(raw, p)
        Theta_rawb = svar_irf(raw, p_raw)
        # true IRFs as (H+1, var, shock)
        Theta_true = np.zeros_like(Theta_long)
        for shk in SHOCKS:
            g = true[true.shock == shk].sort_values("h")
            Theta_true[:, :, NAMES.index(shk)] = g[NAMES].to_numpy()[:H + 1]
        res[case] = dict(mean=m_mean, lo=m_lo, hi=m_hi, long=Theta_long, true=Theta_true,
                         raw1=Theta_raw1, rawb=Theta_rawb, p_long=p_long, p_raw=p_raw,
                         n=len(rep_irf), reps=rep_irf)
        for shk in SHOCKS:
            j = NAMES.index(shk)
            for resp in NAMES:
                i = NAMES.index(resp)
                for h in range(H + 1):
                    rows.append({"case": case, "shock": shk, "response": resp, "h": h,
                                 "data": Theta_d[h, i, j], "data_lo": lo_d[h, i, j],
                                 "data_hi": hi_d[h, i, j], "data_sd": sd_d[h, i, j],
                                 "model_mean": m_mean[h, i, j], "model_lo": m_lo[h, i, j],
                                 "model_hi": m_hi[h, i, j], "model_long": Theta_long[h, i, j],
                                 "long_norm": normalized(Theta_long)[h, i, j],
                                 "true_norm": normalized(Theta_true)[h, i, j]})
    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)

    # ── 1. validation: long-sample SVAR against the true responses ─────────
    say(f"\n{'-'*78}\n1. Validation: long-sample SVARs against the model's true responses\n{'-'*78}")
    say("  Responses of u and v normalized by the shock variable's own response at h = 0; peak")
    say("  (largest |response|) and its quarter. Gap = max over h of |SVAR − true| / max |true|.")
    say("  raw: unfiltered long sample (the clean test); HP: filtered long sample (filter effect).\n")
    tab = []
    for case in CASES:
        R = res[case]
        T = normalized(R["true"])
        for shk in SHOCKS:
            j = NAMES.index(shk)
            for resp in OUTCOMES:
                i = NAMES.index(resp)
                def cell(Theta):
                    N = normalized(Theta)
                    k = int(np.abs(N[:, i, j]).argmax())
                    gap = np.abs(N[:, i, j] - T[:, i, j]).max() / np.abs(T[:, i, j]).max()
                    return f"{N[k, i, j]:+.3f} (h {k}) {100 * gap:3.0f} %"
                kt = int(np.abs(T[:, i, j]).argmax())
                tab.append([case, q3.LABEL[shk], resp, f"{T[kt, i, j]:+.3f} (h {kt})",
                            cell(R["raw1"]), cell(R["rawb"]), cell(R["long"])])
    say(tabulate(tab, ["case", "shock", "resp", "true", f"raw, p = {p}", "raw, BIC p", f"HP, p = {p}"],
                 tablefmt="simple", disable_numparse=True))
    say("\n  BIC lag on the long samples: " +
        ", ".join(f"{c}: raw p = {res[c]['p_raw']}, HP p = {res[c]['p_long']}" for c in CASES) +
        f" (data p = {p}, used for all small-sample comparisons).")

    # ── 2. data against model, small sample ────────────────────────────────
    say(f"\n{'-'*78}\n2. Impulse responses: data against the model (small-sample mean, T = {len(frame)})\n{'-'*78}")
    say("  One-SD shocks, log points x 100; peak (largest |response|) and its quarter.")
    say("  Data: point [16, 84 % bootstrap]. Model: mean over windows [16, 84 % across windows].\n")
    tab = []
    for shk in SHOCKS:
        j = NAMES.index(shk)
        for resp in NAMES:
            i = NAMES.index(resp)
            if resp not in OUTCOMES and resp != shk:
                continue
            kd = int(np.abs(Theta_d[:, i, j]).argmax())
            cells = [q3.LABEL[shk], q3.LABEL.get(resp, resp),
                     f"{100*Theta_d[kd,i,j]:+.2f} (h {kd}) [{100*lo_d[kd,i,j]:+.2f}, {100*hi_d[kd,i,j]:+.2f}]"]
            for case in CASES:
                M = res[case]
                km = int(np.abs(M["mean"][:, i, j]).argmax())
                cells.append(f"{100*M['mean'][km,i,j]:+.2f} (h {km}) [{100*M['lo'][km,i,j]:+.2f}, {100*M['hi'][km,i,j]:+.2f}]")
            tab.append(cells)
    say(tabulate(tab, ["shock", "resp", "data"] + list(CASES), tablefmt="simple", disable_numparse=True))
    say("  Own responses (z → z, δ → δ, s → s) check the shock processes (part6b), not the propagation.")

    # ── 3. distance ─────────────────────────────────────────────────────────
    say(f"\n{'-'*78}\n3. Distance: J = Σ ((model − data)/sd_data)², u and v, h = 0–{H}\n{'-'*78}")
    draws_scaled = np.array([normalized(D) for D in draws])
    sd_scaled = draws_scaled.std(axis=0)
    Theta_d_scaled = normalized(Theta_d)

    def distance(M, D, sd):
        parts = {}
        for shk in SHOCKS:
            j = NAMES.index(shk)
            parts[shk] = sum((((M[:, NAMES.index(r), j] - D[:, NAMES.index(r), j])
                               / sd[:, NAMES.index(r), j]) ** 2).sum() for r in OUTCOMES)
        return parts

    tab = []
    for case in CASES:
        for lab, M, D, sd in (("one SD", res[case]["mean"], Theta_d, sd_d),
                              ("scaled", np.array([normalized(X) for X in res[case]["reps"]]).mean(axis=0),
                               Theta_d_scaled, sd_scaled)):
            parts = distance(M, D, sd)
            tab.append([CASE_LABEL[case], lab] + [f"{parts[s]:.1f}" for s in SHOCKS] +
                       [f"{sum(parts.values()):.1f}"])
    say(tabulate(tab, ["case", "responses", "J z", "J δ", "J s", "J total"], tablefmt="simple",
                 disable_numparse=True))
    n_mom = len(OUTCOMES) * len(SHOCKS) * (H + 1)
    say(f"  {n_mom} moments; J ≈ {n_mom} would mean the model is within about one data SD on average.")
    say("\n  Scaled peaks (response per 1 % own impact; data, default, compromise):")
    tab = []
    for shk in SHOCKS:
        j = NAMES.index(shk)
        for resp in OUTCOMES:
            i = NAMES.index(resp)
            cells = [q3.LABEL[shk], resp]
            for X in [Theta_d_scaled] + [np.array([normalized(Y) for Y in res[c]["reps"]]).mean(axis=0) for c in CASES]:
                k = int(np.abs(X[:, i, j]).argmax())
                cells.append(f"{X[k, i, j]:+.3f} (h {k})")
            tab.append(cells)
    say(tabulate(tab, ["shock", "resp", "data", "default", "compromise"], tablefmt="simple",
                 disable_numparse=True))

    # ── 4. delta: the A_F target ────────────────────────────────────────────
    svF = q3.prepare_svar("A_F", q3.SVARS_F["A_F"], cyc)
    pF, namesF = svF["p_bic"], svF["names"]
    fitF = q3.fit_var(svF["Y"], pF)
    ThF = q3.irf(fitF, namesF, namesF, H)
    drF = q3.bootstrap_irf(svF["Y"], fitF, namesF, namesF, H, q3.N_BOOT)
    jF = namesF.index("delta")                   # A_F is in its own variable order
    iFd = namesF.index("delta")
    dF   = {r: ThF[:, namesF.index(r), jF] for r in ("delta",) + OUTCOMES}
    sdF  = {r: drF[:, :, namesF.index(r), jF].std(axis=0) for r in OUTCOMES}
    drF_s = drF[:, :, :, jF] / drF[:, 0, iFd, jF][:, None, None]          # scaled, draw by draw
    dF_s  = {r: dF[r] / dF["delta"][0] for r in OUTCOMES}
    sdF_s = {r: drF_s[:, :, namesF.index(r)].std(axis=0) for r in OUTCOMES}
    loF = {r: np.percentile(drF[:, :, namesF.index(r), jF], 16, axis=0) for r in OUTCOMES}
    hiF = {r: np.percentile(drF[:, :, namesF.index(r), jF], 84, axis=0) for r in OUTCOMES}
    frF = svF["frame"]

    say(f"\n{'-'*78}\n4. The delta shock against two data targets: A and A_F (net of NFCI)\n{'-'*78}")
    say(f"  A_F: {frF.index[0]}-{frF.index[-1]}, T = {len(frF)}, p = {pF} (BIC), ordering {', '.join(namesF)}.")
    say("  The model side is SVAR A in both comparisons (the model has no NFCI): A_F is a bound.\n")
    jd, idd = NAMES.index("delta"), NAMES.index("delta")
    scaled_mean = {c: np.array([normalized(X) for X in res[c]["reps"]]).mean(axis=0) for c in CASES}

    def pk(x):
        k = int(np.abs(x).argmax())
        return f"{100 * x[k]:+.2f} (h {k})"

    def pk_s(x):
        k = int(np.abs(x).argmax())
        return f"{x[k]:+.3f} (h {k})"

    tab = [["own impact", f"{100 * Theta_d[0, idd, jd]:.2f}", f"{100 * dF['delta'][0]:.2f}"] +
           [f"{100 * res[c]['mean'][0, idd, jd]:.2f}" for c in CASES]]
    for r in OUTCOMES:
        i = NAMES.index(r)
        tab.append([f"δ → {r}, one SD", pk(Theta_d[:, i, jd]), pk(dF[r])] +
                   [pk(res[c]["mean"][:, i, jd]) for c in CASES])
        tab.append([f"δ → {r}, per 1 %", pk_s(Theta_d_scaled[:, i, jd]), pk_s(dF_s[r])] +
                   [pk_s(scaled_mean[c][:, i, jd]) for c in CASES])
    say(tabulate(tab, ["", "data A", "data A_F", "model default", "model compromise"],
                 tablefmt="simple", disable_numparse=True))

    tab = []
    for c in CASES:
        cells = [CASE_LABEL[c]]
        for scale in ("one SD", "scaled"):
            for target in ("A", "A_F"):
                J = 0.0
                for r in OUTCOMES:
                    i = NAMES.index(r)
                    if scale == "one SD":
                        M = res[c]["mean"][:, i, jd]
                        if target == "A":
                            D, sd = Theta_d[:, i, jd], sd_d[:, i, jd]
                        else:
                            D, sd = dF[r], sdF[r]
                    else:
                        M = scaled_mean[c][:, i, jd]
                        if target == "A":
                            D, sd = Theta_d_scaled[:, i, jd], sd_scaled[:, i, jd]
                        else:
                            D, sd = dF_s[r], sdF_s[r]
                    J += (((M - D) / sd) ** 2).sum()
                cells.append(f"{J:.1f}")
        tab.append(cells)
    say("\n  J δ (u and v, h = 0–12, 26 moments), against each data target:")
    say(tabulate(tab, ["case", "one SD, A", "one SD, A_F", "scaled, A", "scaled, A_F"],
                 tablefmt="simple", disable_numparse=True))
    for r in OUTCOMES:
        sel = (out.shock == "delta") & (out.response == r)
        for c in CASES:
            out.loc[sel & (out.case == c), "data_AF"] = dF[r]
            out.loc[sel & (out.case == c), "data_AF_sd"] = sdF[r]
    out.to_csv(OUT_CSV, index=False)

    # ── figure ──────────────────────────────────────────────────────────────
    pal = {"default": "#2a78d6", "compromise": "#eb6834"}
    hh = np.arange(H + 1)
    fig, axes = plt.subplots(len(SHOCKS), len(OUTCOMES), figsize=(10, 9), sharex=True)
    for a, shk in enumerate(SHOCKS):
        j = NAMES.index(shk)
        for b, resp in enumerate(OUTCOMES):
            i, ax = NAMES.index(resp), axes[a, b]
            ax.fill_between(hh, 100 * lo_d[:, i, j], 100 * hi_d[:, i, j], color="#c3c2b7",
                            alpha=0.6, lw=0, label="data, 68 % band")
            ax.plot(hh, 100 * Theta_d[:, i, j], color="#52514e", lw=2, label="data")
            if shk == "delta":
                ax.plot(hh, 100 * dF[resp], color="#52514e", lw=1.5, ls=":",
                        label="data, A_F (net of NFCI)")
            for case in CASES:
                ax.plot(hh, 100 * res[case]["mean"][:, i, j], color=pal[case], lw=2,
                        label=f"model, {CASE_LABEL[case]}")
            ax.axhline(0, color="#8a8984", lw=1, ls="--")
            ax.set_title(f"{q3.LABEL[shk]} shock → {resp}", fontsize=10)
            ax.grid(alpha=0.15)
            if b == 0:
                ax.set_ylabel("log points × 100")
            if a == len(SHOCKS) - 1:
                ax.set_xlabel("quarters")
    axes[0, 0].legend(fontsize=7, frameon=False)
    axes[1, 0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"SVAR A: data ({frame.index[0]}–{frame.index[-1]}) against the AGS model, "
                 f"small-sample mean (T = {len(frame)})", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}\nSaved: {OUT_PNG}")
