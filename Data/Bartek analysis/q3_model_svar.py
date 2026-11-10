"""
q3_model_svar.py -- Block B' model side: SVAR A on simulated AGS data, against the data.

The same estimator as q3_data_svar.py (imported, not copied): SVAR A on HP-1600 cycles of
(z, s, delta, u, v), recursive in that order, a constant, the data's BIC lag. Only the responses
to the z, s and delta shocks are used. The model has three shocks and the SVAR five variables, so
the model's u and v residuals are close to combinations of the first three; with u and v ordered
last, the Cholesky columns of z, s and delta stay well defined. Variance shares are not compared
(they need measurement error or residual shocks: the Block M decision).

Simulated data: Programs AGS/ags_svar_sim.jl, with Block B''s shock processes (part6b section 9:
z and s re-matched on SVAR A's window, delta as in S14; decisions.md D13), two cases:
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

Conditional second moments (section 5). Reported, not targeted (decisions.md S16). The second
moments of u, v, theta (= v - u in logs) and z implied by the identified z, delta and s shocks,
each alone and together: with orthonormal shocks, var(x | J) = sum_h sum_{j in J} Theta_h(x, j)^2,
cov(x, y | J) = sum_h sum_j Theta_h(x, j) Theta_h(y, j), and the first autocorrelation from
sum_h Theta_h Theta_{h+1}, summed to h = 200. They are functions of the impulse responses (all
horizons), so they add no information to them; they summarize volatility and comovement, and flag
joint misfits (e.g. cor(u, v | s)) that a diagonal distance underweights. Both sides go through the
SVAR: the data's SVAR A, and the same SVAR on each simulated window (mean over windows). The model's
true-shock statistics (the sweeps) are not comparable with these: the SVAR recovers only about half
of delta's true effect (section 1).

Outputs:
    data/results/q3_model_svar.txt   the report
    data/results/q3_model_svar.csv   IRFs: data (point, 16/84 %), model small-sample mean and
                                     16/84 %, long sample, true (normalized), by case; for the
                                     delta shock also data_AF (point, sd) from A_F
    data/results/q3_model_svar.png   u and v responses to z, s, delta: data band, model cases
    data/results/q3_model_svar_moments.csv  conditional second moments: data (point, 16/84 %),
                                     model cases (mean, 16/84 % across windows)
    data/results/q3_model_svar_grid.png  z -> u, v by per-hire share, at each xi_inv (section 6)

Case grid (section 6; Nov 10, 2026). Can a large per-hire cost kappa amplify z at moderate
inelasticity of vacancy creation (xi_inv 0.5 and 1), where the comovement moments survive? The
grid crosses xi_inv in {1, 0.5} with x_v in {0.5, 0.25, 0.1, 0.05} (kappa's share of a hire
50-95 %), sigma = 1, Block B''s processes; one row per case with the conditional moments and J.
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
# Cases: the list ags_svar_sim.jl simulated (Programs AGS/ags_svar_cases.csv). The detailed
# sections (1, 2, 4, 5, figure) use the two headline cases; section 6 and the distance cover all.
MANIFEST  = pd.read_csv(MODEL_DIR / "ags_svar_cases.csv").set_index("case")
CASES     = tuple(MANIFEST.index)
HEADLINE  = ("default", "compromise")


def case_label(c):
    r = MANIFEST.loc[c]
    lab = f"x_v {r.x_v:g}, ξ_inv {r.xi_inv:g}"
    if c in HEADLINE:
        lab = f"{c} ({lab})"
    return lab


CASE_LABEL = {c: case_label(c) for c in CASES}
OUT_GRID  = p6b.RESULTS_DIR / "q3_model_svar_grid.png"
NAMES     = list(q3.SVARS["A"]["vars"])          # z, s, delta, u, v
SHOCKS    = q3.SHOCKS                            # z, delta, s
OUTCOMES  = q3.OUTCOMES                          # u, v
H         = q3.H_IRF
OUT_TXT   = p6b.RESULTS_DIR / "q3_model_svar.txt"
OUT_CSV   = p6b.RESULTS_DIR / "q3_model_svar.csv"
OUT_PNG   = p6b.RESULTS_DIR / "q3_model_svar.png"
OUT_MOM   = p6b.RESULTS_DIR / "q3_model_svar_moments.csv"
H_MOM     = 200                                  # horizons summed for the conditional moments
COND_SETS = (("z", ("z",)), ("δ", ("delta",)), ("s", ("s",)), ("z, δ, s", ("z", "delta", "s")))
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def svar_irf(frame, p):
    """SVAR A's IRFs (H+1, var, shock), recursive in NAMES order, on a frame of cycles."""
    Y = frame[NAMES].to_numpy()
    return q3.irf(q3.fit_var(Y, p), NAMES, NAMES, H)


def cond_moments(Theta):
    """
    Second moments implied by the identified shocks (orthonormal), for each set in COND_SETS:
    σ(u), σ(v), σ(θ) in %, cor(u, v), cor(u, z), ρ1(u), ρ1(v). Theta: (H+1, var, shock) in
    NAMES order; θ = v − u. cor(u, z) is reported only when z moves (sets containing z).
    """
    iu, iv, iz = NAMES.index("u"), NAMES.index("v"), NAMES.index("z")
    out = {}
    for lab, shocks in COND_SETS:
        J = [NAMES.index(s) for s in shocks]
        U, V, Z = Theta[:, iu, J], Theta[:, iv, J], Theta[:, iz, J]
        TH = V - U
        var = lambda X: (X ** 2).sum()
        cov = lambda X, Y: (X * Y).sum()
        ac1 = lambda X: (X[:-1] * X[1:]).sum() / var(X)
        row = {"σ(u)": 100 * np.sqrt(var(U)), "σ(v)": 100 * np.sqrt(var(V)),
               "σ(θ)": 100 * np.sqrt(var(TH)), "cor(u,v)": cov(U, V) / np.sqrt(var(U) * var(V)),
               "ρ1(u)": ac1(U), "ρ1(v)": ac1(V)}
        if "z" in shocks:
            row["cor(u,z)"] = cov(U, Z) / np.sqrt(var(U) * var(Z))
            row["σ(θ)/σ(z)"] = np.sqrt(var(TH) / var(Z))
        out[lab] = row
    return out


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
    for case in HEADLINE:
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
        ", ".join(f"{c}: raw p = {res[c]['p_raw']}, HP p = {res[c]['p_long']}" for c in HEADLINE) +
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
            for case in HEADLINE:
                M = res[case]
                km = int(np.abs(M["mean"][:, i, j]).argmax())
                cells.append(f"{100*M['mean'][km,i,j]:+.2f} (h {km}) [{100*M['lo'][km,i,j]:+.2f}, {100*M['hi'][km,i,j]:+.2f}]")
            tab.append(cells)
    say(tabulate(tab, ["shock", "resp", "data"] + list(HEADLINE), tablefmt="simple", disable_numparse=True))
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
            for X in [Theta_d_scaled] + [np.array([normalized(Y) for Y in res[c]["reps"]]).mean(axis=0) for c in HEADLINE]:
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
           [f"{100 * res[c]['mean'][0, idd, jd]:.2f}" for c in HEADLINE]]
    for r in OUTCOMES:
        i = NAMES.index(r)
        tab.append([f"δ → {r}, one SD", pk(Theta_d[:, i, jd]), pk(dF[r])] +
                   [pk(res[c]["mean"][:, i, jd]) for c in HEADLINE])
        tab.append([f"δ → {r}, per 1 %", pk_s(Theta_d_scaled[:, i, jd]), pk_s(dF_s[r])] +
                   [pk_s(scaled_mean[c][:, i, jd]) for c in HEADLINE])
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

    # ── 5. conditional second moments, through the SVAR on both sides ───────
    say(f"\n{'-'*78}\n5. Conditional second moments implied by the identified shocks (reported, not targeted)\n{'-'*78}")
    say(f"  Sums over h = 0–{H_MOM} of the SVAR's impulse responses; data: point [16, 84 % bootstrap];")
    say("  model: mean over windows [16, 84 % across windows]. Both sides through SVAR A.\n")
    Y = svA["Y"]
    m_data = cond_moments(q3.irf(fitA, NAMES, NAMES, H_MOM))
    drawsL = q3.bootstrap_irf(Y, fitA, NAMES, NAMES, H_MOM, q3.N_BOOT)
    m_draws = [cond_moments(D) for D in drawsL]
    m_model = {}
    for case in CASES:
        reps = pd.read_csv(MODEL_DIR / f"ags_svar_reps_{case}.csv")
        m_model[case] = [cond_moments(q3.irf(q3.fit_var(g[NAMES].to_numpy(), p), NAMES, NAMES, H_MOM))
                         for _, g in reps.groupby("rep")]
    # total VAR-implied σ(u) (all five shocks), for the three shocks' share
    tot_data = 100 * np.sqrt((q3.irf(fitA, NAMES, NAMES, H_MOM)[:, NAMES.index("u"), :] ** 2).sum())
    rows_csv, tab = [], []
    for lab, _ in COND_SETS:
        for stat in m_data[lab]:
            d = m_data[lab][stat]
            bs = np.array([x[lab][stat] for x in m_draws])
            cells = [lab, stat, f"{d:+.2f} [{np.percentile(bs, 16):+.2f}, {np.percentile(bs, 84):+.2f}]"]
            rec = {"set": lab, "stat": stat, "data": d, "data_lo": np.percentile(bs, 16),
                   "data_hi": np.percentile(bs, 84)}
            for case in CASES:
                ms = np.array([x[lab][stat] for x in m_model[case]])
                if case in HEADLINE:
                    cells.append(f"{ms.mean():+.2f} [{np.percentile(ms, 16):+.2f}, {np.percentile(ms, 84):+.2f}]")
                rec |= {f"{case}_mean": ms.mean(), f"{case}_lo": np.percentile(ms, 16),
                        f"{case}_hi": np.percentile(ms, 84)}
            tab.append(cells)
            rows_csv.append(rec)
    say(tabulate(tab, ["shocks", "statistic", "data", "default", "compromise"], tablefmt="simple",
                 disable_numparse=True))
    say(f"\n  Data: σ(u) from all five SVAR shocks {tot_data:.2f}; the three model shocks give "
        f"{m_data['z, δ, s']['σ(u)']:.2f}, {100 * (m_data['z, δ, s']['σ(u)'] / tot_data) ** 2:.0f} % of the variance.")
    say("  σ in % (log points × 100) of the HP cycle; θ = v − u. cor(u,z) and σ(θ)/σ(z) only where z moves.")
    pd.DataFrame(rows_csv).to_csv(OUT_MOM, index=False)

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
            for case in HEADLINE:
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

    # ── 6. the case grid: z amplification through κ at moderate inelasticity ─
    say(f"\n{'-'*78}\n6. Case grid: ξ_inv × x_v (κ's share of a hire = 1 − x_v), through the SVAR\n{'-'*78}")
    say("  Conditional moments (section 5) and the one-SD distance (section 3), every case. Model: mean")
    say("  over windows. Data: SVAR A. 'cor uz|z' = cor(u, z) under z alone; 'θ/z' = σ(θ)/σ(z) under z.\n")
    head = ["case", "ϕ", "σu z", "σu δ", "σu s", "σu all", "cor uv z", "cor uv δ", "cor uv s",
            "cor uz|z", "θ/z", "ρ1 u|z", "J z", "J δ", "J s", "J"]

    def grid_row(label, mom, J=None, phi=None):
        cells = [label, "" if phi is None else f"{phi:.2f}"]
        cells += [f"{mom[s]['σ(u)']:.2f}" for s in ("z", "δ", "s", "z, δ, s")]
        cells += [f"{mom[s]['cor(u,v)']:+.2f}" for s in ("z", "δ", "s")]
        cells += [f"{mom['z']['cor(u,z)']:+.2f}", f"{mom['z']['σ(θ)/σ(z)']:.1f}", f"{mom['z']['ρ1(u)']:.2f}"]
        if J is None:
            cells += ["", "", "", ""]
        else:
            cells += [f"{J[s]:.0f}" for s in SHOCKS] + [f"{sum(J.values()):.0f}"]
        return cells

    tab = [grid_row("data", m_data)]
    for case in CASES:
        mean_mom = {lab: {st: np.mean([x[lab][st] for x in m_model[case]]) for st in m_model[case][0][lab]}
                    for lab, _ in COND_SETS}
        J = distance(res[case]["mean"], Theta_d, sd_d)
        tab.append(grid_row(CASE_LABEL[case], mean_mom, J, MANIFEST.loc[case].phi))
    say(tabulate(tab, head, tablefmt="simple", disable_numparse=True))
    say("  J: one-SD responses of u and v, h = 0–12 (26 moments per shock, 78 in total).")

    # figure: z → u and z → v by x_v, one row per ξ_inv
    grid_pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]     # dataviz reference palette, fixed order
    xis = sorted(MANIFEST.xi_inv.unique(), reverse=True)
    fig2, axes2 = plt.subplots(len(xis), len(OUTCOMES), figsize=(10, 3.6 * len(xis)), sharex=True)
    jz = NAMES.index("z")
    for a, xi in enumerate(xis):
        cs = MANIFEST[MANIFEST.xi_inv == xi].sort_values("x_v", ascending=False)
        for b, resp in enumerate(OUTCOMES):
            i, ax = NAMES.index(resp), axes2[a, b]
            ax.fill_between(hh, 100 * lo_d[:, i, jz], 100 * hi_d[:, i, jz], color="#c3c2b7", alpha=0.6,
                            lw=0, label="data, 68 % band")
            ax.plot(hh, 100 * Theta_d[:, i, jz], color="#52514e", lw=2, label="data")
            for k, (case, r) in enumerate(cs.iterrows()):
                ax.plot(hh, 100 * res[case]["mean"][:, i, jz], color=grid_pal[k], lw=2,
                        label=f"x_v {r.x_v:g} (κ share {100 * (1 - r.x_v):.0f} %)")
            ax.axhline(0, color="#8a8984", lw=1, ls="--")
            ax.set_title(f"z shock → {resp}, ξ_inv = {xi:g}", fontsize=10)
            ax.grid(alpha=0.15)
            if b == 0:
                ax.set_ylabel("log points × 100")
            if a == len(xis) - 1:
                ax.set_xlabel("quarters")
    axes2[0, 0].legend(fontsize=7, frameon=False)
    fig2.suptitle("SVAR A, z shock: data against AGS by per-hire share (small-sample mean)", fontsize=11)
    fig2.tight_layout()
    fig2.savefig(OUT_GRID, dpi=150)

    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}\nSaved: {OUT_PNG}\nSaved: {OUT_MOM}\nSaved: {OUT_GRID}")
