"""
part5_joint_lp.py -- Joint Panel Local Projections (delta + LD simultaneously)
===============================================================================
Replication check: runs delta and LD instruments jointly in a single regression
at each horizon h, on v2-residualized instruments.  Addresses the OVB concern
that separate LPs blend the two structural effects when corr(B^delta, B^LD) != 0.

Specification (joint, horizon h):
    y_{s,t+h} - y_{s,t-1} = a_s + a_t
                             + beta^delta_h * B^delta_{s,t}
                             + beta^LD_h    * B^LD_{s,t}
                             + gamma_1 * y_{s,t-1} + gamma_2 * log(LF_{s,t-1})
                             + eps_{s,t,h}

Both instruments scaled to 1-SD units *before* regression so coefficients are
directly comparable to the separate LPs in part5_lp.py (which post-multiply).
Here we pre-standardize so that beta^delta and beta^LD are each interpretable as
"effect of a 1-SD move in that instrument, holding the other fixed."

Outputs (data/results/)
-----------------------
  lp_joint_delta_ld_unemp.csv     -- h x (beta^delta, beta^LD, SEs, p-vals)
  lp_joint_delta_ld_vacancy.csv   -- same for vacancy outcome
  lp_joint_delta_ld_unemp.png     -- overlay: joint vs separate beta^delta & beta^LD
  lp_joint_delta_ld_vacancy.png   -- same for vacancy

Diagnostics printed to console
-------------------------------
  - Post-v2 corr(B^delta, B^LD) in the merged sample
  - VIF for both instruments at each horizon
  - Side-by-side comparison table: joint vs separate coefficients

Prerequisites
-------------
  part3_resid_instruments.py (v2 run)
  part4_outcomes.py
  part5_lp.py results CSVs (for comparison loading)
"""

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
    pass

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HORIZONS  = list(range(21))   # h = 0 .. 20 (match part5_lp.py)
BASE_YEAR = 2006
Z90       = 1.645
Z95       = 1.960
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DELTA_FILE = INSTR_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
LD_FILE    = INSTR_DIR / f"ld_instrument_resid_base{BASE_YEAR}.csv"
LAUS_FILE  = INSTR_DIR / "laus_quarterly.parquet"

# Separate LP CSVs from part5_lp.py -- loaded for comparison only
SEP_DELTA_UNEMP  = RESULTS_DIR / "lp_irf_delta_resid.csv"
SEP_LD_UNEMP     = RESULTS_DIR / "lp_irf_ld_resid.csv"
SEP_DELTA_VAC    = RESULTS_DIR / "lp_irf_delta_vacancy.csv"
SEP_LD_VAC       = RESULTS_DIR / "lp_irf_ld_vacancy.csv"

# ---------------------------------------------------------------------------
# Quarter arithmetic
# ---------------------------------------------------------------------------
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


# ---------------------------------------------------------------------------
# Data loading and panel construction
# ---------------------------------------------------------------------------
def load_data():
    delta = pd.read_csv(DELTA_FILE, dtype={"state_fips": str})
    ld    = pd.read_csv(LD_FILE,    dtype={"state_fips": str})
    delta["state_fips"] = delta["state_fips"].str.zfill(2)
    ld["state_fips"]    = ld["state_fips"].str.zfill(2)

    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)

    return delta, ld, laus


def build_joint_panel(delta: pd.DataFrame, ld: pd.DataFrame,
                      laus: pd.DataFrame) -> pd.DataFrame:
    """
    Merge both instruments with outcomes on (state_fips, quarter_label).
    Pre-standardize both instrument columns to unit SD so joint betas are
    directly on a 1-SD scale (matching the post-multiply scaling in part5_lp.py).
    """
    keep_laus = ["state_fips", "quarter_label", "unemp_rate", "labor_force"]
    has_vac   = "vacancies" in laus.columns
    if has_vac:
        keep_laus.append("vacancies")

    panel = (delta[["state_fips", "state", "quarter_label", "bartik_delta"]]
             .merge(ld[["state_fips", "quarter_label", "bartik_ld"]],
                    on=["state_fips", "quarter_label"], how="inner")
             .merge(laus[keep_laus],
                    on=["state_fips", "quarter_label"], how="inner"))

    # Convert to pp units then standardize
    panel["bartik_delta"] = panel["bartik_delta"] * 100
    panel["bartik_ld"]    = panel["bartik_ld"]    * 100

    sd_delta = panel["bartik_delta"].std()
    sd_ld    = panel["bartik_ld"].std()
    panel["B_delta"] = panel["bartik_delta"] / sd_delta
    panel["B_ld"]    = panel["bartik_ld"]    / sd_ld

    print(f"\n  Merged panel: {panel['state_fips'].nunique()} states x "
          f"{panel['quarter_label'].nunique()} quarters = {len(panel):,} rows")
    print(f"  Quarter range: {panel['quarter_label'].min()} - "
          f"{panel['quarter_label'].max()}")
    print(f"  SD(B^delta) before std: {sd_delta:.4f}   "
          f"SD(B^LD) before std: {sd_ld:.4f}")
    print(f"  corr(B^delta, B^LD): "
          f"{panel['bartik_delta'].corr(panel['bartik_ld']):.4f}")

    # Lagged controls
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)
    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )
    if has_vac:
        panel["vac_rate"]     = panel["vacancies"] * 1000 / panel["labor_force"] * 100
        panel["vac_rate_lag1"] = panel.groupby("state_fips")["vac_rate"].shift(1)

    # Zero out lags across quarter gaps
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    bad = panel["prev_ql"].notna() & (
        panel["prev_ql"] != quarter_shift(panel["quarter_label"], -1)
    )
    lag_cols = ["unemp_lag1", "lf_log_lag1"] + (["vac_rate_lag1"] if has_vac else [])
    panel.loc[bad, lag_cols] = np.nan
    panel.drop(columns="prev_ql", inplace=True)

    return panel, sd_delta, sd_ld


# ---------------------------------------------------------------------------
# VIF helper
# ---------------------------------------------------------------------------
def vif_pair(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float]:
    """Return VIF for x1 and x2 in a bivariate regression."""
    r = np.corrcoef(x1, x2)[0, 1]
    vif = 1.0 / (1.0 - r**2) if abs(r) < 1.0 else np.inf
    return round(vif, 3), round(vif, 3)


# ---------------------------------------------------------------------------
# Joint LP -- single horizon
# ---------------------------------------------------------------------------
def run_joint_horizon(panel: pd.DataFrame, h: int,
                      outcome: str = "unemp") -> dict | None:
    """
    Estimate joint LP at horizon h.  Returns dict with coefficients for
    both B_delta and B_ld, plus diagnostics.
    """
    df = panel.copy()
    df["future_ql"] = quarter_shift(df["quarter_label"], h)

    if outcome == "unemp":
        lookup = (df[["state_fips", "quarter_label", "unemp_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "unemp_rate": "y_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"]   = df["y_future"] - df["unemp_lag1"]
        lag_control     = "unemp_lag1"

    elif outcome == "vacancy":
        if "vac_rate" not in df.columns:
            return None
        lookup = (df[["state_fips", "quarter_label", "vac_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "vac_rate": "vac_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["vac_future"] - df["vac_rate_lag1"]
        lag_control   = "vac_rate_lag1"
    else:
        raise ValueError(f"outcome must be 'unemp' or 'vacancy', got {outcome!r}")

    required = ["dep_var", "B_delta", "B_ld", lag_control, "lf_log_lag1"]
    df = df.dropna(subset=required).copy()
    if len(df) < 100:
        return None

    state_dummies = pd.get_dummies(df["state_fips"],    prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)
    core = ["B_delta", "B_ld", lag_control, "lf_log_lag1"]
    X = sm.add_constant(
        pd.concat([df[core], state_dummies, time_dummies], axis=1),
        has_constant="add"
    )
    model = sm.OLS(df["dep_var"], X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values}
    )

    def _extract(col):
        b = float(model.params[col])
        s = float(model.bse[col])
        t = float(model.tvalues[col])
        p = float(model.pvalues[col])
        return b, s, t, p

    b_d, se_d, t_d, p_d = _extract("B_delta")
    b_l, se_l, t_l, p_l = _extract("B_ld")
    vif_d, vif_l = vif_pair(df["B_delta"].values, df["B_ld"].values)

    return {
        "h": h,
        # delta
        "beta_delta":   b_d,   "se_delta":  se_d,
        "tstat_delta":  t_d,   "pval_delta": p_d,
        "ci90_lo_delta": b_d - Z90*se_d, "ci90_hi_delta": b_d + Z90*se_d,
        "ci95_lo_delta": b_d - Z95*se_d, "ci95_hi_delta": b_d + Z95*se_d,
        # LD
        "beta_ld":      b_l,   "se_ld":     se_l,
        "tstat_ld":     t_l,   "pval_ld":   p_l,
        "ci90_lo_ld":   b_l - Z90*se_l, "ci90_hi_ld":   b_l + Z90*se_l,
        "ci95_lo_ld":   b_l - Z95*se_l, "ci95_hi_ld":   b_l + Z95*se_l,
        # diagnostics
        "vif_delta": vif_d, "vif_ld": vif_l,
        "nobs": int(model.nobs),
        "r2": float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range": f"{df['quarter_label'].min()}-{df['quarter_label'].max()}",
        "outcome": outcome,
    }


# ---------------------------------------------------------------------------
# Joint LP -- all horizons
# ---------------------------------------------------------------------------
def run_joint_lp(panel: pd.DataFrame, outcome: str = "unemp") -> pd.DataFrame:
    out_tag = "vacancy rate" if outcome == "vacancy" else "unemp rate"
    print(f"\n  Joint LP -> {out_tag}   h=0..{max(HORIZONS)}")
    print(f"  {'h':>3}  {'b_d':>8}  {'p_d':>6}  {'b_ld':>8}  {'p_ld':>6}  "
          f"{'VIF':>5}  {'N':>6}")
    print(f"  {'-'*60}")

    rows = []
    for h in HORIZONS:
        res = run_joint_horizon(panel, h, outcome=outcome)
        if res is None:
            continue
        rows.append(res)
        sig_d = "***" if res["pval_delta"] < .01 else "**" if res["pval_delta"] < .05 \
                else "*" if res["pval_delta"] < .10 else ""
        sig_l = "***" if res["pval_ld"] < .01 else "**" if res["pval_ld"] < .05 \
                else "*" if res["pval_ld"] < .10 else ""
        print(f"  {h:3d}  {res['beta_delta']:8.4f}{sig_d:<3}  {res['pval_delta']:6.3f}  "
              f"{res['beta_ld']:8.4f}{sig_l:<3}  {res['pval_ld']:6.3f}  "
              f"{res['vif_delta']:5.2f}  {res['nobs']:6d}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Comparison table: joint vs separate
# ---------------------------------------------------------------------------
def print_comparison(joint: pd.DataFrame, sep_csv: Path,
                     shock: str, outcome: str) -> None:
    """Print side-by-side joint vs separate coefficients at h=0..8."""
    if not sep_csv.exists():
        print(f"  [skip comparison] {sep_csv.name} not found")
        return
    sep = pd.read_csv(sep_csv).set_index("h")
    beta_col = "beta_delta" if shock == "delta" else "beta_ld"
    print(f"\n  {'h':>3}  {'joint':>9}  {'p_jt':>6}  {'separate':>9}  {'p_sep':>6}  "
          f"  diff = joint - separate")
    print(f"  {'-'*60}")
    for h in range(min(17, max(joint["h"]) + 1)):
        if h not in joint["h"].values:
            continue
        r    = joint[joint["h"] == h].iloc[0]
        b_jt = r[beta_col]
        p_jt = r[f"pval_{shock}"]
        if h in sep.index:
            b_sep = float(sep.loc[h, "beta"])
            p_sep = float(sep.loc[h, "pval"])
        else:
            b_sep, p_sep = np.nan, np.nan
        diff = b_jt - b_sep if not np.isnan(b_sep) else np.nan
        print(f"  {h:3d}  {b_jt:9.4f}  {p_jt:6.3f}  {b_sep:9.4f}  {p_sep:6.3f}  "
              f"  {diff:+.4f}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
_PALETTE = {
    "delta_joint":    ("#d62728", "-"),   # red  -- joint
    "delta_separate": ("#1f77b4", "-"),   # blue -- separate
    "ld_joint":       ("#d62728", "-"),   # red  -- joint
    "ld_separate":    ("#1f77b4", "-"),   # blue -- separate
}


def plot_joint_vs_separate(joint: pd.DataFrame,
                            sep_delta_csv: Path, sep_ld_csv: Path,
                            out_path: Path, outcome: str = "unemp") -> None:
    """
    Two-panel figure:
      Left:  delta -- joint (solid) vs separate (dashed)
      Right: LD    -- joint (solid) vs separate (dashed)
    """
    ylabel = ("pp change in vacancy rate\nper 1-SD shock"
              if outcome == "vacancy"
              else "pp change in unemp. rate\nper 1-SD shock")

    sep_d = pd.read_csv(sep_delta_csv).set_index("h") if sep_delta_csv.exists() else None
    sep_l = pd.read_csv(sep_ld_csv).set_index("h")   if sep_ld_csv.exists()    else None

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    COLOR_JOINT    = "#d62728"   # red  -- joint LP (controlling for other shock)
    COLOR_SEPARATE = "#1f77b4"   # blue -- separate LP (bivariate)

    for ax, shock, beta_col, pval_col, sep_df, title_str in [
        (axes[0], "delta", "beta_delta", "pval_delta", sep_d,
         r"$\delta$ shock: joint vs separate LP"),
        (axes[1], "ld",    "beta_ld",    "pval_ld",    sep_l,
         "LD shock: joint vs separate LP"),
    ]:
        h_vals = joint["h"].values

        # Joint -- red
        b_jt = joint[beta_col].values
        lo90 = joint[f"ci90_lo_{shock}"].values
        hi90 = joint[f"ci90_hi_{shock}"].values
        ax.fill_between(h_vals, lo90, hi90, color=COLOR_JOINT, alpha=0.15)
        ax.plot(h_vals, b_jt, color=COLOR_JOINT, lw=2.0, ls="-",
                marker="o", ms=3.5, label="joint LP")

        # Separate -- blue
        if sep_df is not None:
            hs   = [h for h in h_vals if h in sep_df.index]
            b_sp = [float(sep_df.loc[h, "beta"]) for h in hs]
            lo_s = [float(sep_df.loc[h, "ci90_lo"]) for h in hs]
            hi_s = [float(sep_df.loc[h, "ci90_hi"]) for h in hs]
            ax.fill_between(hs, lo_s, hi_s, color=COLOR_SEPARATE, alpha=0.15)
            ax.plot(hs, b_sp, color=COLOR_SEPARATE, lw=1.5, ls="-",
                    marker="s", ms=3.0, label="separate LP")

        ax.axhline(0, color="black", lw=0.8)
        for hh in range(4, int(h_vals.max()) + 1, 4):
            ax.axvline(hh, color="grey", lw=0.4, ls=":", alpha=0.5)
        ax.set_title(title_str, fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xticks(h_vals[::2] if h_vals.max() > 16 else h_vals)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: str(int(x)))
        )
        ax.legend(fontsize=9, framealpha=0.85)
        ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle(
        "Joint vs separate LP: v2-residualized instruments\n"
        "(solid = joint; dashed = separate; shading = 90% CI)",
        fontsize=12
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5 -- Joint LP (delta + LD, v2-residualized)")
print("=" * 60)

# [1] Load
print("\n[1] Loading data")
delta_resid, ld_resid, laus = load_data()

# [2] Build joint panel
print("\n[2] Building joint panel")
panel, sd_delta, sd_ld = build_joint_panel(delta_resid, ld_resid, laus)
has_vac = "vac_rate" in panel.columns

#Inspect 
panel.head(5)
panel.columns

# Variance inflation factor: collinearity diagnostics for joint LP
r_instr    = panel["bartik_delta"].corr(panel["bartik_ld"])
vif        = 1 / (1 - r_instr**2)
eigenvals  = np.linalg.eigvalsh(
    np.array([[1, r_instr],[r_instr, 1]])
)
kappa      = np.sqrt(eigenvals.max() / eigenvals.min())
se_inflate = np.sqrt(vif)

print(f"  r(δ,LD) = {r_instr:.4f}")
print(f"  VIF     = {vif:.4f}  (SE inflation = {se_inflate:.4f}×)")
print(f"  κ       = {kappa:.4f}  (condition number)")

# [3] Joint LP -- unemployment
print("\n[3] Joint LP -- unemployment outcome")
joint_unemp = run_joint_lp(panel, outcome="unemp")

# [4] Joint LP -- vacancy
if has_vac:
    print("\n[4] Joint LP -- vacancy outcome")
    joint_vac = run_joint_lp(panel, outcome="vacancy")
else:
    print("\n[4] Vacancy data not available -- skipping vacancy LP")
    joint_vac = pd.DataFrame()

# [5] Save results
print("\n[5] Saving results")
if not joint_unemp.empty:
    p = RESULTS_DIR / "lp_joint_delta_ld_unemp.csv"
    joint_unemp.to_csv(p, index=False)
    print(f"  Saved: {p}")

if not joint_vac.empty:
    p = RESULTS_DIR / "lp_joint_delta_ld_vacancy.csv"
    joint_vac.to_csv(p, index=False)
    print(f"  Saved: {p}")

# [6] Comparison tables
print("\n[6] Comparison: joint vs separate")

print(f"\n  --- delta -> unemployment ---")
print_comparison(joint_unemp, SEP_DELTA_UNEMP, "delta", "unemp")

print(f"\n  --- LD -> unemployment ---")
print_comparison(joint_unemp, SEP_LD_UNEMP, "ld", "unemp")

if not joint_vac.empty:
    print(f"\n  --- delta -> vacancy ---")
    print_comparison(joint_vac, SEP_DELTA_VAC, "delta", "vacancy")

    print(f"\n  --- LD -> vacancy ---")
    print_comparison(joint_vac, SEP_LD_VAC, "ld", "vacancy")

# [7] Plots
print("\n[7] Plotting")
plot_joint_vs_separate(
    joint_unemp,
    SEP_DELTA_UNEMP, SEP_LD_UNEMP,
    RESULTS_DIR / "lp_joint_delta_ld_unemp.png",
    outcome="unemp"
)

if not joint_vac.empty:
    plot_joint_vs_separate(
        joint_vac,
        SEP_DELTA_VAC, SEP_LD_VAC,
        RESULTS_DIR / "lp_joint_delta_ld_vacancy.png",
        outcome="vacancy"
    )

# [8] Summary: peak responses
print("\n[8] Peak responses -- joint LP (1-SD scale)")
print(f"  {'Shock':<12}  {'outcome':>8}  {'h_peak':>6}  {'beta_peak':>9}  "
      f"{'SE':>7}  {'p':>6}")
print(f"  {'-'*58}")

for shock, bcol, secol, pcol, irf in [
    ("delta", "beta_delta", "se_delta", "pval_delta", joint_unemp),
    ("LD",    "beta_ld",    "se_ld",    "pval_ld",    joint_unemp),
    ("delta", "beta_delta", "se_delta", "pval_delta", joint_vac),
    ("LD",    "beta_ld",    "se_ld",    "pval_ld",    joint_vac),
]:
    if irf is None or irf.empty or bcol not in irf.columns:
        continue
    out = irf["outcome"].iloc[0] if "outcome" in irf.columns else "?"
    idx = irf[bcol].abs().idxmax()
    row = irf.loc[idx]
    print(f"  {shock:<12}  {out:>8}  {int(row['h']):6d}  "
          f"{row[bcol]:9.4f}  {row[secol]:7.4f}  {row[pcol]:6.3f}")

print("\nDone.")
