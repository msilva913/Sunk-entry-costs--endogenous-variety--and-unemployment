"""
part6d_shock_components.py -- is each shock a persistent AR(1) plus short-lived noise?

part6b matches a single monthly AR(1) to each shock's filtered quarterly data, under HP and
under the Hamilton filter separately (shock_calibration.txt, sections 7-8). Because the
simulated series is filtered exactly as the data are, the two filters should recover the
same monthly process if the AR(1) is right. They do not: the half-life of delta is about
one month under HP and six under Hamilton, and similar gaps hold for s and z. The AR(1)
therefore misses variation at the frequencies where the two filters differ.

The candidate tested here: the observed monthly series is a persistent AR(1) plus
independent monthly noise,
    observed_t = sqrt(1 - w) x_t + sqrt(w) n_t,   x_t = rho x_{t-1} + sqrt(1 - rho^2) e_t,
with x and n of unit variance and w the noise share of the monthly variance. The noise can
be measurement error (BED deaths) or a genuinely short-lived component; the statistics
cannot tell the two apart, and the interpretation is a modeling decision.

Targets (part6b, shock_calibration.csv): the quarterly AR(1) coefficient and residual SD
under HP(1,600) and under Hamilton (h = 8, p = 4) -- four moments. Parameters: rho, w and
an overall scale -- three. Fit: the three scale-free moments (rho_Q under each filter and
the ratio of the two residual SDs) over a (rho, w) grid, then the scale from the two
residual SDs. One overidentifying moment remains, reported as the fit's residual.

Simulation as in part6b Step 4b: monthly draws, averaged to quarters, filtered over a
data-length sample (Hamilton drops its first 11 quarters, as in the data), AR(1) by OLS
with a constant; means over replications with common random numbers.

delta's targets come from the (z, delta) VAR (own lag, innovation purged of z); fitting the
components to them treats delta as univariate, an approximation noted in the output.

Diagnostic only: SHOCKS_BASE (run_solution_core.jl) is not changed.

Outputs:
    data/results/shock_components.txt
    data/results/shock_components.csv
"""

import numpy as np
import pandas as pd

import part6b_shock_persistence_cyclical as p6b

OUT_TXT = p6b.RESULTS_DIR / "shock_components.txt"
OUT_CSV = p6b.RESULTS_DIR / "shock_components.csv"
SERIES  = ["z", "delta", "s"]
RHO_GRID = np.concatenate([np.linspace(0.0, 0.95, 39), np.linspace(0.955, 0.999, 23)])
W_GRID   = np.concatenate([[0.0], np.linspace(0.02, 0.98, 49)])
LINES = []


def say(text=""):
    LINES.append(text)
    print(text)


def ar1_stats(c):
    """part6b's AR(1) estimator, replication by replication: mean rho_Q and residual SD."""
    y, z = c[:, 1:], c[:, :-1]
    yd = y - y.mean(axis=1, keepdims=True)
    zd = z - z.mean(axis=1, keepdims=True)
    b = (yd * zd).sum(axis=1) / (zd * zd).sum(axis=1)
    resid = yd - b[:, None] * zd
    sig = np.sqrt((resid ** 2).sum(axis=1) / (resid.shape[1] - 2))
    return float(b.mean()), float(sig.mean())


def quarterly(monthly, nq):
    """Drop the burn-in and average months to quarters."""
    R = monthly.shape[0]
    return monthly[:, p6b.BURN_M:p6b.BURN_M + 3 * nq].reshape(R, nq, 3).mean(axis=2)


def persistent_paths(rho, e):
    """Unit-variance monthly AR(1) on the draws e (replications x months)."""
    x = np.empty_like(e)
    x[:, 0] = e[:, 0]
    scale = np.sqrt(1 - rho ** 2)
    for k in range(1, e.shape[1]):
        x[:, k] = rho * x[:, k - 1] + scale * e[:, k]
    return x


def moment_surface(nq, rng):
    """
    For every (rho, w) on the grid: rho_Q and the residual SD under HP and Hamilton,
    for an observed series of unit monthly variance on nq quarters.
    """
    months = p6b.BURN_M + 3 * nq
    e = rng.standard_normal((p6b.N_MATCH, months))
    n = rng.standard_normal((p6b.N_MATCH, months))
    C = p6b.hp_cycle_matrix(nq)
    noise_q = quarterly(n, nq)
    rows = []
    for rho in RHO_GRID:
        x_q = quarterly(persistent_paths(rho, e), nq)
        for w in W_GRID:
            q = np.sqrt(1 - w) * x_q + np.sqrt(w) * noise_q
            r_hp, s_hp = ar1_stats(q @ C.T)
            r_ham, s_ham = ar1_stats(p6b.hamilton_rows(q))
            rows.append({"rho_m": rho, "w": w, "rho_hp": r_hp, "sig_hp": s_hp,
                         "rho_ham": r_ham, "sig_ham": s_ham})
    return pd.DataFrame(rows)


def fit(surface, target):
    """Best (rho, w) on the three scale-free moments; the scale from both residual SDs."""
    d = surface.copy()
    log_ratio_data = np.log(target["sigma_hp"] / target["sigma_ham"])
    d["dist"] = ((d.rho_hp - target["rho_hp"]) ** 2 + (d.rho_ham - target["rho_ham"]) ** 2
                 + (np.log(d.sig_hp / d.sig_ham) - log_ratio_data) ** 2)
    best = d.loc[d.dist.idxmin()].copy()
    # One scale for both filters: the geometric mean of the two implied scales.
    best["scale"] = np.sqrt(target["sigma_hp"] / best.sig_hp * target["sigma_ham"] / best.sig_ham)
    return best


def half_life(rho):
    """Half-life in months of a monthly AR(1) with coefficient rho."""
    if rho <= 0:
        return 0.0
    return float(np.log(0.5) / np.log(rho))


if __name__ == "__main__":
    calib = pd.read_csv(p6b.RESULTS_DIR / "shock_calibration.csv").set_index("series")
    say("=" * 78)
    say("Persistent AR(1) plus noise: one monthly process for both filters (part6d)")
    say("=" * 78)
    say(__doc__.strip().split("\n\n")[1])
    say()
    say("Targets (shock_calibration.csv): quarterly AR(1) rho_Q and residual SD x 100.")
    say("Single AR(1) = part6b's matched monthly process under each filter (w = 0 by construction).")
    say("Components = this program's fit of (rho_m, w) to both filters at once.\n")

    rng = np.random.default_rng(p6b.MATCH_SEED)
    out = []
    for name in SERIES:
        c = calib.loc[name]
        nq = int(c["nobs"]) + 1                    # quarters before filtering (AR(1) loses one)
        target = {"rho_hp": c["rho_Q"], "sigma_hp": c["sigma_Q"],
                  "rho_ham": c["hamilton_rho_Q"], "sigma_ham": c["hamilton_sigma_Q"]}
        surface = moment_surface(nq, rng)
        best = fit(surface, target)
        ar1 = fit(surface[surface.w == 0.0], target)
        for lab, b in (("single AR(1), both filters", ar1), ("AR(1) + noise", best)):
            out.append({"series": name, "spec": lab, "nq": nq,
                        "rho_m": b.rho_m, "half_life_m": half_life(b.rho_m), "w": b.w,
                        "rho_hp_fit": b.rho_hp, "rho_ham_fit": b.rho_ham,
                        "sig_hp_fit": 100 * b.sig_hp * b.scale, "sig_ham_fit": 100 * b.sig_ham * b.scale,
                        "rho_hp_data": target["rho_hp"], "rho_ham_data": target["rho_ham"],
                        "sig_hp_data": 100 * target["sigma_hp"], "sig_ham_data": 100 * target["sigma_ham"],
                        "dist": b.dist,
                        "sd_persistent_m": b.scale * np.sqrt(1 - b.w),
                        "sd_noise_m": b.scale * np.sqrt(b.w)})
        # How sharply the data pin w down: best distance at each w.
        prof = surface.copy()
        log_ratio_data = np.log(target["sigma_hp"] / target["sigma_ham"])
        prof["dist"] = ((prof.rho_hp - target["rho_hp"]) ** 2 + (prof.rho_ham - target["rho_ham"]) ** 2
                        + (np.log(prof.sig_hp / prof.sig_ham) - log_ratio_data) ** 2)
        profile = prof.loc[prof.groupby("w").dist.idxmin()]
        say(f"{'-' * 78}\n{name}: {nq} quarters before filtering\n{'-' * 78}")
        say(f"  data: rho_Q HP {target['rho_hp']:.3f}, Hamilton {target['rho_ham']:.3f}; "
            f"sigma_Q HP {100 * target['sigma_hp']:.2f}, Hamilton {100 * target['sigma_ham']:.2f}")
        say(f"  part6b, one AR(1) per filter: rho_m HP {c['matched_rho_m']:.4f} "
            f"(half-life {half_life(c['matched_rho_m']):.1f} months), Hamilton "
            f"{c['hamilton_matched_rho_m']:.4f} ({half_life(c['hamilton_matched_rho_m']):.1f} months)")
        for row in out[-2:]:
            say(f"  {row['spec']:<27} rho_m {row['rho_m']:.4f} (half-life {row['half_life_m']:6.1f} m), "
                f"noise share w {row['w']:.2f}")
            say(f"  {'':<27} fit: rho_Q HP {row['rho_hp_fit']:.3f}, Ham {row['rho_ham_fit']:.3f}; "
                f"sigma_Q HP {row['sig_hp_fit']:.2f}, Ham {row['sig_ham_fit']:.2f}; distance {row['dist']:.4f}")
        say("  profile (best rho_m at each noise share w; distance):")
        for _, p in profile.iloc[::5].iterrows():
            say(f"    w {p.w:.2f}: rho_m {p.rho_m:.4f} (half-life {half_life(p.rho_m):6.1f} m), "
                f"distance {p.dist:.4f}")
        say()
    say("delta's targets are the (z, delta) VAR's own lag and z-purged innovation; the")
    say("components treat delta as univariate.")
    pd.DataFrame(out).to_csv(OUT_CSV, index=False)
    OUT_TXT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
    print(f"\nSaved: {OUT_TXT}\nSaved: {OUT_CSV}")
