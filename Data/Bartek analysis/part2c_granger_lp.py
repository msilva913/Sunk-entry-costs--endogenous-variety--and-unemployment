"""
part2c_granger_lp.py  --  Panel LP Granger causality between shock series
==========================================================================
Tests whether δ (establishment closings) and s (match separations) are
causally ordered or reflect independent structural disturbances, conditional
on common and industry-specific productivity variation.

Economic motivation
-------------------
The model has three primitive shocks: productivity (z), product-line
destruction (δ), and match separation (s).  A concern is that δ might be
largely endogenous — rising as a lagged response when s was high and firms
were losing workers — in which case it is redundant to include both as
separate structural shocks.

We test causal ordering using a panel LP Granger causality design.  The key
questions are:

  (1) Does LD Granger-cause δ?
      If yes: closings are partly an endogenous response to prior separations.
      If no: δ is not forecastable from the history of s — consistent with
             it being a primitive disturbance.

  (2) Does δ Granger-cause LD?
      If yes and (1) is rejected: δ is the primitive shock, and separations
             follow as an endogenous consequence — the model's endogenous-
             exit mechanism.

Specification
-------------
For each pair (Y, X) and each horizon h = 0, ..., H:

    Y_{j,t+h} = α_j + α_t + Σ_{l=1}^{L} β_l Y_{j,t-l}
                            + Σ_{l=1}^{L} γ_l X_{j,t-l}  +  ε_{j,t+h}

Industry FEs (α_j) absorb industry-level means; time FEs (α_t) absorb
aggregate shocks including the aggregate productivity cycle.  Estimated by
two-way within transformation.  SEs are HC1 (White) robust.

We report γ_l coefficients and a joint Wald F-test on {γ_1,...,γ_L}.

Productivity conditioning
-------------------------
Both δ and s respond to productivity shocks.  Time FEs absorb aggregate
productivity variation, but industry-specific productivity is not absorbed.
We therefore run two versions of every test:

  Baseline  : raw log shock rates, full sample 2001Q1–2021Q4
  Preferred : residualized shock series ν^k from part2b (industry-specific
              productivity and tightness partialled out), sample 2005Q3–2021Q4

The comparison isolates whether Granger patterns survive conditioning on
industry-level productivity and demand variation — the same controls used
to construct the Bartik instruments.

The residualization is performed internally from cached parquets to avoid
relying on NTFS-mount parquet files (which can be corrupted by Windows I/O).

Run
---
    python part2c_granger_lp.py

Prerequisites
-------------
    python part2_shock_rates.py
    python part2_shock_rates_s.py
    python part2b_shock_comovement.py  (populates cache parquets)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── working directory ─────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import (
    DEFAULT_CACHE_DIR, SHOCK_RATES_PATH, INDUSTRY_LABELS,
)
from construct_s_instrument import (
    SHOCK_RATES_S_PATH, SHOCK_RATES_LD_PATH, SHOCK_RATES_QU_PATH,
)

RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── parameters ────────────────────────────────────────────────────────────────
N_LAGS   = 4
HORIZONS = list(range(0, 13))
CI_LEVEL = 0.90

# ── helpers ───────────────────────────────────────────────────────────────────
def _ql_to_dt(ql):
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)

def _dt_to_ql(dt):
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"

def _open_file(path):
    try:
        os.startfile(path)
    except AttributeError:
        pass

REC_SPANS = [("2001-03-01","2001-11-01"),
             ("2007-12-01","2009-06-01"),
             ("2020-01-01","2020-07-01")]


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 1: BUILD PANELS
#  Two panels are constructed:
#    nat_raw  -- raw log shock rates, 2001Q1–2021Q4 (12 industries × ~84 qtrs)
#    nat_resid -- residualized ν^k series, 2005Q3–2021Q4 (controls partialled)
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 65)
print("Part 2c  --  Panel LP Granger causality between shock series")
print("=" * 65)

missing = [p for p in [SHOCK_RATES_PATH, SHOCK_RATES_S_PATH,
                        SHOCK_RATES_LD_PATH, SHOCK_RATES_QU_PATH]
           if not p.exists()]
if missing:
    sys.exit("Missing shock rate files:\n" + "\n".join(f"  {p}" for p in missing))

def _nat_avg(path, col, new_col):
    """Average LOO rates across states → national industry-quarter series."""
    df = pd.read_parquet(path)
    return (df.groupby(["industry_code","quarter_label"])[col]
              .mean().reset_index().rename(columns={col: new_col}))

# ── 1a. raw panel ─────────────────────────────────────────────────────────────
print("\n--- Building raw shock panel ---")
nat_d  = _nat_avg(SHOCK_RATES_PATH,    "g_delta_loo", "g_delta")
nat_s  = _nat_avg(SHOCK_RATES_S_PATH,  "g_s_loo",     "g_s")
nat_ld = _nat_avg(SHOCK_RATES_LD_PATH, "g_ld_loo",    "g_ld")
nat_qu = _nat_avg(SHOCK_RATES_QU_PATH, "g_qu_loo",    "g_qu")

nat_raw = (nat_d.merge(nat_s,  on=["industry_code","quarter_label"], how="inner")
                .merge(nat_ld, on=["industry_code","quarter_label"], how="inner")
                .merge(nat_qu, on=["industry_code","quarter_label"], how="inner"))
nat_raw = nat_raw[(nat_raw.g_delta>0)&(nat_raw.g_s>0)&
                  (nat_raw.g_ld>0)&(nat_raw.g_qu>0)].copy()
for raw, log in [("g_delta","log_g_delta"),("g_s","log_g_s"),
                 ("g_ld","log_g_ld"),("g_qu","log_g_qu")]:
    nat_raw[log] = np.log(nat_raw[raw])
nat_raw["industry_label"] = nat_raw["industry_code"].map(INDUSTRY_LABELS)
nat_raw = nat_raw.sort_values(["industry_code","quarter_label"]).reset_index(drop=True)
print(f"  Raw panel: {nat_raw['quarter_label'].min()}–{nat_raw['quarter_label'].max()}, "
      f"{nat_raw['industry_code'].nunique()} industries, {len(nat_raw):,} obs")

# ── 1b. residualized panel: rebuild ν^k in-process from cache ─────────────────
# Avoids reading NTFS-mount parquet files which can be corrupted by Windows I/O.
# Replicates the part2b v2 residualization: project log g^k on industry FEs,
# Δlog p_t, Δlog VA_{j,t-1}, and log θ_{t-1}; retain residuals.
print("\n--- Building residualized shock panel (v2 spec) ---")

required_caches = [
    DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet",
    DEFAULT_CACHE_DIR / "national_tightness_quarterly.parquet",
    DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet",
]
missing_caches = [p for p in required_caches if not p.exists()]
if missing_caches:
    print("  WARNING: missing cache files — residualized panel unavailable.")
    print("  Run part2b_shock_comovement.py first with FRED_API_KEY set.")
    for p in missing_caches:
        print(f"    {p}")
    nat_resid = None
else:
    # Attach controls to raw panel
    nat_resid = nat_raw[["industry_code","quarter_label","industry_label",
                          "log_g_delta","log_g_ld","log_g_qu"]].copy()
    nat_resid["industry_code"] = nat_resid["industry_code"].astype(int)

    # Aggregate productivity
    prod = pd.read_parquet(DEFAULT_CACHE_DIR / "ophnfb_quarterly.parquet")
    prod = prod.sort_values("quarter_label").reset_index(drop=True)
    prod["log_p"]  = np.log(prod["productivity"])
    prod["dlog_p"] = prod["log_p"].diff()
    nat_resid = nat_resid.merge(prod[["quarter_label","dlog_p"]].dropna(),
                                on="quarter_label", how="inner")

    # Lagged market tightness
    tight = pd.read_parquet(DEFAULT_CACHE_DIR / "national_tightness_quarterly.parquet")
    tight = tight.sort_values("quarter_label").reset_index(drop=True)
    tight["log_theta_lag"] = tight["log_theta"].shift(1)
    nat_resid = nat_resid.merge(
        tight[["quarter_label","log_theta_lag"]].dropna(),
        on="quarter_label", how="inner"
    )

    # Lagged industry VA growth
    bea = pd.read_parquet(DEFAULT_CACHE_DIR / "bea_va_quarterly_12ind.parquet")
    bea = bea.sort_values(["industry_code","quarter_label"]).reset_index(drop=True)
    bea["dlog_va_lag"] = bea.groupby("industry_code")["dlog_va"].shift(1)
    bea_lag = bea[["industry_code","quarter_label","dlog_va_lag"]].dropna().copy()
    bea_lag["industry_code"] = bea_lag["industry_code"].astype(int)
    nat_resid = nat_resid.merge(bea_lag, on=["industry_code","quarter_label"],
                                how="inner")
    nat_resid = nat_resid.dropna().sort_values(
        ["industry_code","quarter_label"]).reset_index(drop=True)

    # Within-industry residualization (same estimator as part2b)
    def _residualize(df, dep_col, resid_col, reg_cols):
        work = df[["industry_code", dep_col, *reg_cols]].copy().reset_index(drop=True)
        work["y_dm"] = (work[dep_col]
                        - work.groupby("industry_code")[dep_col].transform("mean"))
        dm_cols = []
        for rc in reg_cols:
            dc = f"_dm_{rc}"
            work[dc] = work[rc] - work.groupby("industry_code")[rc].transform("mean")
            dm_cols.append(dc)
        X      = work[dm_cols].to_numpy()
        y      = work["y_dm"].to_numpy()
        gammas = np.linalg.lstsq(X, y, rcond=None)[0]
        resid  = y - X @ gammas
        r2 = 1 - (resid**2).sum() / (y**2).sum()
        coef_str = "  ".join(f"{rc}={g:+.4f}" for rc, g in zip(reg_cols, gammas))
        print(f"    {resid_col:<12}  {coef_str}   R²={r2:.4f}")
        out = df[["industry_code","quarter_label"]].copy().reset_index(drop=True)
        out[resid_col] = resid
        return out

    reg_base  = ("dlog_p", "dlog_va_lag")
    reg_tight = ("dlog_p", "dlog_va_lag", "log_theta_lag")

    print("  Residualizing:")
    rd  = _residualize(nat_resid, "log_g_delta", "nu_delta", reg_base)
    rld = _residualize(nat_resid, "log_g_ld",    "nu_ld",    reg_tight)
    rqu = _residualize(nat_resid, "log_g_qu",    "nu_qu",    reg_tight)

    nat_resid["nu_delta"] = rd["nu_delta"].values
    nat_resid["nu_ld"]    = rld["nu_ld"].values
    nat_resid["nu_qu"]    = rqu["nu_qu"].values

    print(f"  Residualized panel: {nat_resid['quarter_label'].min()}–"
          f"{nat_resid['quarter_label'].max()}, "
          f"{nat_resid['industry_code'].nunique()} industries, "
          f"{len(nat_resid):,} obs")

    # Correlation reduction summary
    print("\n  Pooled correlations  (raw → residualized):")
    for (c1r, c2r, c1n, c2n, lbl) in [
        ("log_g_delta","log_g_ld","nu_delta","nu_ld","r(δ,LD)"),
        ("log_g_delta","log_g_qu","nu_delta","nu_qu","r(δ,QU)"),
        ("log_g_ld",   "log_g_qu","nu_ld",  "nu_qu","r(LD,QU)"),
    ]:
        sub = nat_resid.dropna(subset=[c1r,c2r,c1n,c2n])
        rr = sub[c1r].corr(sub[c2r])
        rn = sub[c1n].corr(sub[c2n])
        print(f"    {lbl}: {rr:+.3f} → {rn:+.3f}  (Δ={rn-rr:+.3f})")


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 2: PANEL LP GRANGER ESTIMATOR
# ═══════════════════════════════════════════════════════════════════════════════

def _build_lag_matrix(df, y_col, x_col, n_lags, h,
                      industry_col="industry_code", time_col="quarter_label"):
    industries = sorted(df[industry_col].unique())
    quarters   = sorted(df[time_col].unique())
    q_index    = {q: i for i, q in enumerate(quarters)}
    ind_index  = {ind: i for i, ind in enumerate(industries)}

    rows = []
    for ind in industries:
        sub = df[df[industry_col]==ind].sort_values(time_col).reset_index(drop=True)
        T = len(sub)
        for t_idx in range(n_lags, T - h):
            y_val = sub.loc[t_idx + h, y_col]
            own   = [sub.loc[t_idx - l, y_col] for l in range(1, n_lags+1)]
            cross = [sub.loc[t_idx - l, x_col] for l in range(1, n_lags+1)]
            qt    = sub.loc[t_idx, time_col]
            rows.append({
                "y": y_val,
                **{f"own_lag{l}":   own[l-1]   for l in range(1, n_lags+1)},
                **{f"cross_lag{l}": cross[l-1] for l in range(1, n_lags+1)},
                "industry_id": ind_index[ind],
                "time_id":     q_index[qt],
                "industry":    ind,
                "quarter":     qt,
            })
    return pd.DataFrame(rows).dropna()


def _within_transform(data, group_col, cols):
    out = data.copy()
    for col in cols:
        out[col] = out[col] - out.groupby(group_col)[col].transform("mean")
    return out


def _norm_cdf(x):
    from math import erf, sqrt
    return np.array([0.5*(1+erf(xi/sqrt(2))) for xi in np.atleast_1d(x)])


def _ols_with_f(X, y, n_cross_lags):
    n, k = X.shape
    if n <= k:
        return None
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta    = XtX_inv @ X.T @ y
    resid   = y - X @ beta
    df_res  = n - k
    scale   = n / df_res
    meat    = (X * resid[:,None]).T @ (X * resid[:,None])
    V_hc1   = scale * XtX_inv @ meat @ XtX_inv
    se      = np.sqrt(np.diag(V_hc1))
    tstat   = beta / np.where(se>0, se, np.nan)
    pval    = 2*(1 - _norm_cdf(np.abs(tstat)))

    idx_cross = list(range(k - n_cross_lags, k))
    R = np.zeros((n_cross_lags, k))
    for i, j in enumerate(idx_cross):
        R[i, j] = 1.0
    Rb  = R @ beta
    RVR = R @ V_hc1 @ R.T
    try:
        f_stat = (Rb @ np.linalg.solve(RVR, Rb)) / n_cross_lags
    except np.linalg.LinAlgError:
        f_stat = np.nan

    return {"beta": beta, "se": se, "tstat": tstat, "pval": pval,
            "f_stat": f_stat, "nobs": n, "df_res": df_res}


def run_granger_lp(df, y_col, x_col, y_label, x_label,
                   n_lags=N_LAGS, horizons=HORIZONS):
    """
    Panel LP Granger test: does X Granger-cause Y?
    Returns DataFrame with h, gamma coefficients, SEs, p-values, joint F.
    """
    all_regs = (["y"] + [f"own_lag{l}" for l in range(1,n_lags+1)]
                       + [f"cross_lag{l}" for l in range(1,n_lags+1)])
    results = []
    for h in horizons:
        data = _build_lag_matrix(df, y_col, x_col, n_lags, h)
        if len(data) < 30:
            continue
        dm = _within_transform(data, "industry", all_regs)
        dm = _within_transform(dm,   "quarter",  all_regs)
        X_mat = dm[[f"own_lag{l}"   for l in range(1,n_lags+1)] +
                   [f"cross_lag{l}" for l in range(1,n_lags+1)]].to_numpy()
        y_vec = dm["y"].to_numpy()
        res   = _ols_with_f(X_mat, y_vec, n_cross_lags=n_lags)
        if res is None:
            continue
        gammas     = res["beta"][-n_lags:]
        gamma_ses  = res["se"][-n_lags:]
        gamma_pval = res["pval"][-n_lags:]
        betas      = res["beta"][:n_lags]
        row = {"h": h, "f_stat": res["f_stat"], "nobs": res["nobs"]}
        for l in range(1, n_lags+1):
            row[f"gamma_{l}"]    = gammas[l-1]
            row[f"gamma_{l}_se"] = gamma_ses[l-1]
            row[f"gamma_{l}_p"]  = gamma_pval[l-1]
            row[f"beta_{l}"]     = betas[l-1]
        results.append(row)
    return pd.DataFrame(results)


def _irf_lp(df, y_col, x_col, n_lags=N_LAGS, H=max(HORIZONS)):
    """LP IRF: response of Y to a 1-SD shock in X_{t-1} at each horizon h."""
    sd_x = df[x_col].std()
    all_regs = (["y"] + [f"own_lag{l}" for l in range(1,n_lags+1)]
                       + [f"cross_lag{l}" for l in range(1,n_lags+1)])
    rows = []
    for h in range(0, H+1):
        data = _build_lag_matrix(df, y_col, x_col, n_lags, h)
        if len(data) < 30:
            continue
        dm = _within_transform(data, "industry", all_regs)
        dm = _within_transform(dm,   "quarter",  all_regs)
        X_mat = dm[[f"own_lag{l}"   for l in range(1,n_lags+1)] +
                   [f"cross_lag{l}" for l in range(1,n_lags+1)]].to_numpy()
        y_vec = dm["y"].to_numpy()
        res   = _ols_with_f(X_mat, y_vec, n_cross_lags=n_lags)
        if res is None:
            continue
        gamma_1    = res["beta"][-n_lags]
        gamma_1_se = res["se"][-n_lags]
        rows.append({"h": h,
                     "irf":     gamma_1 * sd_x,
                     "irf_se":  gamma_1_se * sd_x,
                     "gamma_1": gamma_1,
                     "pval":    res["pval"][-n_lags]})
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 3: RUN TESTS — BASELINE (RAW) AND PREFERRED (RESIDUALIZED)
# ═══════════════════════════════════════════════════════════════════════════════

# The four economically meaningful pairs:
#   δ → LD: does destruction cause subsequent separations?    (model prediction)
#   LD → δ: do separations predict future closings?           (endogeneity concern)
#   δ → QU: does destruction raise quits?                     (placebo-type check)
#   QU → δ: do quits predict closings?                        (tightness channel)
PAIRS = [
    ("log_g_ld",    "log_g_delta", "nu_ld",    "nu_delta", "LD",  "δ",
     "Does δ Granger-cause LD?\n(destruction → separations)"),
    ("log_g_delta", "log_g_ld",    "nu_delta", "nu_ld",    "δ",   "LD",
     "Does LD Granger-cause δ?\n(separations → closings)"),
    ("log_g_qu",    "log_g_delta", "nu_qu",    "nu_delta", "QU",  "δ",
     "Does δ Granger-cause QU?\n(destruction → quits)"),
    ("log_g_delta", "log_g_qu",    "nu_delta", "nu_qu",    "δ",   "QU",
     "Does QU Granger-cause δ?\n(quits → closings)"),
]

raw_gr, raw_irf, res_gr, res_irf = {}, {}, {}, {}

for (y_raw, x_raw, y_res, x_res, y_lbl, x_lbl, title) in PAIRS:
    key = (y_lbl, x_lbl)

    print(f"\n{'='*65}")
    print(f"BASELINE (raw log rates, 2001Q1–2021Q4)")
    print(f"{title.replace(chr(10),' | ')}")
    raw_gr[key]  = run_granger_lp(nat_raw,  y_raw, x_raw, y_lbl, x_lbl)
    raw_irf[key] = _irf_lp(nat_raw, y_raw, x_raw)

    if nat_resid is not None:
        print(f"\n  PREFERRED (residualized ν^k, 2005Q3–2021Q4)")
        res_gr[key]  = run_granger_lp(nat_resid, y_res, x_res, y_lbl, x_lbl)
        res_irf[key] = _irf_lp(nat_resid, y_res, x_res)

# ── summary: joint F-stats side by side ──────────────────────────────────────
print(f"\n{'='*65}")
print("JOINT F-STATS  (H0: all cross-series lags γ_l = 0 jointly)")
print(f"  HC1 robust SEs  |  L={N_LAGS} lags  |  Industry + time FEs")
print(f"\n  {'Direction':<38}  {'h=0':>6}  {'h=2':>6}  {'h=4':>6}  {'h=8':>6}  {'h=12':>6}")
print(f"  {'-'*70}")

for (y_raw, x_raw, y_res, x_res, y_lbl, x_lbl, title) in PAIRS:
    key   = (y_lbl, x_lbl)
    label = title.replace("\n"," ")
    for tag, grdict in [("raw", raw_gr), ("resid", res_gr)]:
        gr = grdict.get(key, pd.DataFrame())
        if gr.empty:
            continue
        vals = []
        for h in [0, 2, 4, 8, 12]:
            row = gr[gr["h"]==h]
            vals.append(f"{row['f_stat'].values[0]:6.2f}" if len(row) else "   N/A")
        spec_tag = f"[{tag}]"
        print(f"  {label[:33]:<33} {spec_tag:<7}  {'  '.join(vals)}")
    print()

print(f"  Approx critical values (large N): F(4,∞) 10%≈2.0  5%≈2.4  1%≈3.3")


# ═══════════════════════════════════════════════════════════════════════════════
#  SECTION 4: PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

z_ci = 1.645 if CI_LEVEL == 0.90 else 1.960

# ── Plot 1: IRF grid — raw vs residualized, 4 pairs ─────────────────────────
# Two columns per pair: left=raw, right=residualized

n_pairs = len(PAIRS)
has_resid = nat_resid is not None
n_cols = 2 if has_resid else 1
fig, axes = plt.subplots(n_pairs, n_cols, figsize=(7*n_cols, 4*n_pairs),
                         squeeze=False)

colors_raw   = ["#1f77b4","#d62728","#2ca02c","#ff7f0e"]
colors_resid = ["#1f77b4","#d62728","#2ca02c","#ff7f0e"]

for row_i, (y_raw, x_raw, y_res, x_res, y_lbl, x_lbl, title) in enumerate(PAIRS):
    key = (y_lbl, x_lbl)
    for col_i, (tag, irf_dict, color) in enumerate(
        [("Baseline (raw)", raw_irf, colors_raw[row_i]),
         ("Preferred (residualized ν)", res_irf, colors_resid[row_i])]
        if has_resid else
        [("Baseline (raw)", raw_irf, colors_raw[row_i])]
    ):
        ax  = axes[row_i, col_i]
        irf = irf_dict.get(key, pd.DataFrame())
        if irf.empty:
            ax.set_visible(False)
            continue
        hs    = irf["h"].values
        vals  = irf["irf"].values
        ci_lo = vals - z_ci * irf["irf_se"].values
        ci_hi = vals + z_ci * irf["irf_se"].values

        ax.axhline(0, color="black", linewidth=0.8)
        ax.fill_between(hs, ci_lo, ci_hi, alpha=0.18, color=color)
        ax.plot(hs, vals, color=color, linewidth=2.0, marker="o", markersize=4)
        sig_mask = irf["pval"] < 0.10
        if sig_mask.any():
            ax.scatter(hs[sig_mask.values], vals[sig_mask.values],
                       color=color, s=55, zorder=5)

        ax.set_title(f"{title}\n[{tag}]", fontsize=9, fontweight="bold")
        ax.set_xlabel("Horizon h (quarters)", fontsize=8)
        ax.set_ylabel("Response in log shock rate\n(1-SD impulse in X_{t-1})",
                      fontsize=7)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        ax.set_xticks(range(0, max(HORIZONS)+1, 2))

fig.suptitle(
    f"Panel LP Granger causality  |  {CI_LEVEL:.0%} CI  "
    f"|  L={N_LAGS} own+cross lags  |  Industry + time FEs\n"
    "Filled markers = significant at 10%  |  "
    "Preferred spec conditions on Δlog p, Δlog VA lag, log θ lag",
    fontsize=10, y=1.01,
)
fig.tight_layout()
p = RESULTS_DIR / "granger_lp_irf_grid.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"\nPlot 1 saved: {p}")

# ── Plot 2: F-stat bars — δ↔LD, raw vs residualized ─────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 8))

plot_pairs = [
    ("LD", "δ", "δ → LD\n(destruction causes separations?)", "#1f77b4"),
    ("δ",  "LD","LD → δ\n(separations cause closings?)",     "#d62728"),
]
for col_i, (y_lbl, x_lbl, title, color) in enumerate(plot_pairs):
    key = (y_lbl, x_lbl)
    for row_i, (tag, grdict, alpha) in enumerate(
        [("Baseline (raw)",           raw_gr, 0.85),
         ("Preferred (residualized)", res_gr, 0.85)]
    ):
        ax = axes[row_i, col_i]
        gr = grdict.get(key, pd.DataFrame())
        if gr.empty:
            ax.set_visible(False)
            continue
        ax.bar(gr["h"], gr["f_stat"], color=color, alpha=alpha, width=0.7)
        ax.axhline(2.0, color="black", linewidth=1.1, linestyle="--",
                   label="F crit ≈ 2.0 (10%)")
        ax.axhline(3.3, color="black", linewidth=1.1, linestyle=":",
                   label="F crit ≈ 3.3 (1%)")
        ax.set_title(f"{title}\n[{tag}]", fontsize=10, fontweight="bold")
        ax.set_xlabel("Horizon h", fontsize=12)
        ax.set_ylabel(f"Joint F  ({N_LAGS} cross-series lags)", fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        ax.set_xticks(range(0, max(HORIZONS)+1, 2))

fig.suptitle(
    f"Joint Granger F-statistics: δ vs LD\n"
    f"H0: all {N_LAGS} cross-series lags jointly = 0  |  HC1 robust SEs  "
    "|  Industry + time FEs",
    fontsize=11,
)
fig.tight_layout()
p = RESULTS_DIR / "granger_lp_fstat_delta_ld.png"
fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); _open_file(p)
print(f"Plot 2 saved: {p}")

print("\nDone.")
