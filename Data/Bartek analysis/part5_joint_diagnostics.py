"""
part5_joint_diagnostics.py -- Diagnostics for Joint LP Identification
=====================================================================
Distinguishes two competing interpretations of the zero beta^LD vacancy
coefficient in the joint LP (part5_joint_lp.py):

  Story 1: LD truly has no vacancy effect (reposting works; model confirmed)
  Story 2: LD has a negative vacancy effect but the joint LP cannot identify
           it because delta absorbs the shared variation (identification failure)

Two tests, each run at every LP horizon h=0..20:

Test 1 — Partial R² and Partial F
----------------------------------
Compares R² from three regressions at each h:
  (a) vacancy outcome on delta alone + controls + FEs
  (b) vacancy outcome on LD alone + controls + FEs
  (c) vacancy outcome on both delta + LD (joint) + controls + FEs

Partial F for LD = (R²_joint - R²_delta_only) / (1 - R²_joint) * (N-K) / 1
If partial F > 10: LD has sufficient independent variation to identify its
coefficient; the zero in the joint LP is informative (supports Story 1).
If partial F < ~3: LD's independent variation is too weak; zero is uninformative.

Test 2 — Artificial Orthogonalization
--------------------------------------
Construct B^LD_orth = residual from pooled panel regression of B_ld on B_delta
(controlling for state + time FEs, matching the LP specification). This purges
ALL delta-correlated variation from the LD instrument by construction.

Run the separate LP using B^LD_orth as the sole instrument. If the vacancy IRF
is still flat/zero: the zero is genuine, not a collinearity artifact (Story 1).
If the vacancy IRF is significantly negative: the zero in the joint LP was an
artifact of delta absorbing shared variation (Story 2).

Outputs (data/results/)
-----------------------
  joint_diagnostics_partial_f.csv    -- partial F, partial R² at each h
  joint_diagnostics_orth_vacancy.csv -- orthogonalized LD vacancy IRF
  joint_diagnostics_orth_unemp.csv   -- orthogonalized LD unemployment IRF
  joint_diagnostics.png              -- summary figure

Prerequisites
-------------
  part5_lp.py and part5_joint_lp.py results
  v2-residualized instruments in data/instruments/
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
HORIZONS  = list(range(21))
Z90       = 1.645
Z95       = 1.960
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DELTA_FILE = INSTR_DIR / "delta_instrument_resid_base2006.csv"
LD_FILE    = INSTR_DIR / "ld_instrument_resid_base2006.csv"
LAUS_FILE  = INSTR_DIR / "laus_quarterly.parquet"

# Separate LP results for comparison
SEP_LD_VAC   = RESULTS_DIR / "lp_irf_ld_vacancy.csv"
SEP_LD_UNEMP = RESULTS_DIR / "lp_irf_ld_resid.csv"

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
def load_and_build_panel():
    delta = pd.read_csv(DELTA_FILE, dtype={"state_fips": str})
    ld    = pd.read_csv(LD_FILE,    dtype={"state_fips": str})
    delta["state_fips"] = delta["state_fips"].str.zfill(2)
    ld["state_fips"]    = ld["state_fips"].str.zfill(2)

    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)

    keep_laus = ["state_fips", "quarter_label", "unemp_rate", "labor_force",
                 "vacancies"]
    panel = (delta[["state_fips", "state", "quarter_label", "bartik_delta"]]
             .merge(ld[["state_fips", "quarter_label", "bartik_ld"]],
                    on=["state_fips", "quarter_label"], how="inner")
             .merge(laus[keep_laus],
                    on=["state_fips", "quarter_label"], how="inner"))

    # pp units, then standardize to 1-SD
    panel["bartik_delta"] = panel["bartik_delta"] * 100
    panel["bartik_ld"]    = panel["bartik_ld"]    * 100
    sd_delta = panel["bartik_delta"].std()
    sd_ld    = panel["bartik_ld"].std()
    panel["B_delta"] = panel["bartik_delta"] / sd_delta
    panel["B_ld"]    = panel["bartik_ld"]    / sd_ld

    # Lagged controls
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)
    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )
    panel["vac_rate"]      = panel["vacancies"] * 1000 / panel["labor_force"] * 100
    panel["vac_rate_lag1"] = panel.groupby("state_fips")["vac_rate"].shift(1)

    # Zero out lags across quarter gaps
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    bad = panel["prev_ql"].notna() & (
        panel["prev_ql"] != quarter_shift(panel["quarter_label"], -1)
    )
    panel.loc[bad, ["unemp_lag1", "lf_log_lag1", "vac_rate_lag1"]] = np.nan
    panel.drop(columns="prev_ql", inplace=True)

    corr = panel["B_delta"].corr(panel["B_ld"])
    print(f"  Panel: {panel['state_fips'].nunique()} states x "
          f"{panel['quarter_label'].nunique()} quarters = {len(panel):,} rows")
    print(f"  corr(B^delta, B^LD) = {corr:.4f}")
    print(f"  SD(B^delta) raw pp = {sd_delta:.4f},  SD(B^LD) raw pp = {sd_ld:.4f}")

    return panel, sd_delta, sd_ld


# ---------------------------------------------------------------------------
# Orthogonalize LD on delta (pooled panel, with state + time FEs)
# ---------------------------------------------------------------------------
def orthogonalize_ld(panel: pd.DataFrame) -> pd.DataFrame:
    """
    Regress B_ld on B_delta + state FEs + time FEs (pooled across all
    state-quarters). Residual = B_ld_orth, the component of LD variation
    that is orthogonal to delta after absorbing fixed effects.
    """
    df = panel.dropna(subset=["B_delta", "B_ld"]).copy()
    state_d = pd.get_dummies(df["state_fips"],    prefix="st", drop_first=True, dtype=float)
    time_d  = pd.get_dummies(df["quarter_label"], prefix="qt", drop_first=True, dtype=float)
    X = sm.add_constant(
        pd.concat([df[["B_delta"]], state_d, time_d], axis=1),
        has_constant="add"
    )
    model = sm.OLS(df["B_ld"], X).fit()
    df["B_ld_orth"] = model.resid

    corr_check = df["B_delta"].corr(df["B_ld_orth"])
    sd_orth    = df["B_ld_orth"].std()
    print(f"  Orthogonalization: corr(B^delta, B^LD_orth) = {corr_check:.6f}")
    print(f"  SD(B^LD_orth) = {sd_orth:.4f}  "
          f"(vs SD(B^LD) = {df['B_ld'].std():.4f}, "
          f"ratio = {sd_orth / df['B_ld'].std():.3f})")

    panel = panel.merge(
        df[["state_fips", "quarter_label", "B_ld_orth"]],
        on=["state_fips", "quarter_label"], how="left"
    )
    return panel


# ---------------------------------------------------------------------------
# Build horizon-h estimation sample
# ---------------------------------------------------------------------------
def _build_h_sample(panel: pd.DataFrame, h: int,
                    outcome: str = "vacancy") -> pd.DataFrame | None:
    df = panel.copy()
    df["future_ql"] = quarter_shift(df["quarter_label"], h)

    if outcome == "vacancy":
        lookup = (df[["state_fips", "quarter_label", "vac_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "vac_rate": "vac_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["vac_future"] - df["vac_rate_lag1"]
        lag_control   = "vac_rate_lag1"
    else:  # unemp
        lookup = (df[["state_fips", "quarter_label", "unemp_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "unemp_rate": "y_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["y_future"] - df["unemp_lag1"]
        lag_control   = "unemp_lag1"

    required = ["dep_var", "B_delta", "B_ld", lag_control, "lf_log_lag1"]
    if "B_ld_orth" in df.columns:
        required.append("B_ld_orth")
    df = df.dropna(subset=required).copy()
    if len(df) < 100:
        return None

    df["_lag_control"] = df[lag_control]
    return df


# ---------------------------------------------------------------------------
# Test 1: Partial R² and Partial F
# ---------------------------------------------------------------------------
def test1_partial_f(panel: pd.DataFrame,
                    outcome: str = "vacancy") -> pd.DataFrame:
    """
    At each horizon h, run three regressions on the chosen outcome:
      (a) delta only   (b) LD only   (c) joint (delta + LD)
    Compute partial F for LD = marginal contribution of LD over delta-only,
    and partial F for delta = marginal contribution of delta over LD-only.
    """
    out_tag = "vacancy" if outcome == "vacancy" else "unemployment"
    print(f"\n  Test 1: Partial R² and Partial F ({out_tag} outcome)")
    print(f"  {'h':>3}  {'R2_d':>7}  {'R2_ld':>7}  {'R2_jt':>7}  "
          f"{'pF_ld':>7}  {'pR2_ld':>7}  {'pF_d':>7}  {'pR2_d':>7}  {'N':>6}")
    print(f"  {'-'*72}")

    rows = []
    for h in HORIZONS:
        df = _build_h_sample(panel, h, outcome=outcome)
        if df is None:
            continue

        state_d = pd.get_dummies(df["state_fips"],    prefix="st", drop_first=True, dtype=float)
        time_d  = pd.get_dummies(df["quarter_label"], prefix="qt", drop_first=True, dtype=float)
        controls = pd.concat([df[["_lag_control", "lf_log_lag1"]], state_d, time_d], axis=1)

        def _run(instr_cols):
            X = sm.add_constant(
                pd.concat([df[instr_cols], controls], axis=1),
                has_constant="add"
            )
            return sm.OLS(df["dep_var"], X).fit()

        mod_delta = _run(["B_delta"])
        mod_ld    = _run(["B_ld"])
        mod_joint = _run(["B_delta", "B_ld"])

        r2_d  = mod_delta.rsquared
        r2_l  = mod_ld.rsquared
        r2_j  = mod_joint.rsquared
        N     = int(mod_joint.nobs)
        K     = mod_joint.df_model + 1  # total params including constant

        # Partial R² and F for LD (marginal over delta-only)
        partial_r2_ld = (r2_j - r2_d) / (1.0 - r2_d) if r2_d < 1 else 0
        partial_f_ld  = (r2_j - r2_d) / (1.0 - r2_j) * (N - K) if r2_j < 1 else 0

        # Partial R² and F for delta (marginal over LD-only)
        partial_r2_d  = (r2_j - r2_l) / (1.0 - r2_l) if r2_l < 1 else 0
        partial_f_d   = (r2_j - r2_l) / (1.0 - r2_j) * (N - K) if r2_j < 1 else 0

        rows.append({
            "h": h, "outcome": outcome,
            "r2_delta_only": r2_d, "r2_ld_only": r2_l,
            "r2_joint": r2_j, "nobs": N,
            "partial_f_ld": partial_f_ld, "partial_r2_ld": partial_r2_ld,
            "partial_f_delta": partial_f_d, "partial_r2_delta": partial_r2_d,
        })

        print(f"  {h:3d}  {r2_d:7.4f}  {r2_l:7.4f}  {r2_j:7.4f}  "
              f"{partial_f_ld:7.2f}  {partial_r2_ld:7.4f}  "
              f"{partial_f_d:7.2f}  {partial_r2_d:7.4f}  {N:6d}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Test 2: Orthogonalized LD LP -- single horizon
# ---------------------------------------------------------------------------
def run_orth_lp_horizon(df: pd.DataFrame, h: int,
                        outcome: str = "vacancy") -> dict | None:
    """Run separate LP using B_ld_orth as the sole instrument."""
    sample = _build_h_sample(df, h, outcome=outcome)
    if sample is None or "B_ld_orth" not in sample.columns:
        return None

    state_d = pd.get_dummies(sample["state_fips"],    prefix="st", drop_first=True, dtype=float)
    time_d  = pd.get_dummies(sample["quarter_label"], prefix="qt", drop_first=True, dtype=float)
    core    = ["B_ld_orth", "_lag_control", "lf_log_lag1"]
    X = sm.add_constant(
        pd.concat([sample[core], state_d, time_d], axis=1),
        has_constant="add"
    )
    model = sm.OLS(sample["dep_var"], X).fit(
        cov_type="cluster",
        cov_kwds={"groups": sample["state_fips"].values}
    )

    b  = float(model.params["B_ld_orth"])
    se = float(model.bse["B_ld_orth"])
    t  = float(model.tvalues["B_ld_orth"])
    p  = float(model.pvalues["B_ld_orth"])

    return {
        "h": h, "beta": b, "se": se, "tstat": t, "pval": p,
        "ci90_lo": b - Z90*se, "ci90_hi": b + Z90*se,
        "ci95_lo": b - Z95*se, "ci95_hi": b + Z95*se,
        "nobs": int(model.nobs), "r2": float(model.rsquared),
        "outcome": outcome,
    }


def run_orth_lp(panel: pd.DataFrame, outcome: str = "vacancy") -> pd.DataFrame:
    out_tag = "vacancy" if outcome == "vacancy" else "unemp"
    print(f"\n  Test 2: Orthogonalized LD LP -> {out_tag}   h=0..{max(HORIZONS)}")
    print(f"  {'h':>3}  {'beta':>9}  {'SE':>7}  {'t':>7}  {'p':>6}  {'N':>6}")
    print(f"  {'-'*50}")

    rows = []
    for h in HORIZONS:
        res = run_orth_lp_horizon(panel, h, outcome=outcome)
        if res is None:
            continue
        rows.append(res)
        sig = "***" if res["pval"] < .01 else "**" if res["pval"] < .05 \
              else "*" if res["pval"] < .10 else ""
        print(f"  {h:3d}  {res['beta']:9.4f}  {res['se']:7.4f}  "
              f"{res['tstat']:7.3f}  {res['pval']:6.3f}  {res['nobs']:6d}  {sig}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Summary plot — 2 × 3 grid (rows = vacancy / unemployment)
# ---------------------------------------------------------------------------
def _plot_row(axes_row, pf_df: pd.DataFrame,
              orth_irf: pd.DataFrame, sep_csv: Path,
              outcome_label: str, ylabel_irf: str) -> None:
    """
    Fill one row of the 2×3 figure:
      Col 0: Partial F (delta vs LD)
      Col 1: Orth LD IRF vs separate LD IRF
      Col 2: Partial R²
    """
    sep = pd.read_csv(sep_csv).set_index("h") if sep_csv.exists() else None
    h_vals = pf_df["h"].values
    COLOR_D  = "#1f77b4"
    COLOR_LD = "#d62728"

    # Col 0 — Partial F
    ax = axes_row[0]
    ax.plot(h_vals, pf_df["partial_f_delta"], color=COLOR_D, lw=2,
            marker="o", ms=3.5, label=r"Partial F: $\delta$ (marginal over LD)")
    ax.plot(h_vals, pf_df["partial_f_ld"], color=COLOR_LD, lw=2,
            marker="s", ms=3.5, label="Partial F: LD (marginal over δ)")
    ax.axhline(10, color="grey", lw=1.2, ls="--", alpha=0.7, label="F = 10")
    ax.set_title(f"Partial F — {outcome_label}", fontsize=10)
    ax.set_xlabel("Horizon h (quarters)", fontsize=9)
    ax.set_ylabel("Partial F-statistic", fontsize=9)
    ax.legend(fontsize=7.5, framealpha=0.85)
    ax.grid(axis="y", lw=0.4, alpha=0.4)
    ax.set_xticks(h_vals[::2])

    # Col 1 — Orth LD IRF vs separate
    ax = axes_row[1]
    if not orth_irf.empty:
        h_v = orth_irf["h"].values
        ax.fill_between(h_v, orth_irf["ci90_lo"], orth_irf["ci90_hi"],
                        color=COLOR_LD, alpha=0.15)
        ax.plot(h_v, orth_irf["beta"], color=COLOR_LD, lw=2,
                marker="o", ms=3.5, label="LD orth (δ-purged)")
        ax.set_xticks(h_v[::2])
    if sep is not None:
        hs   = [h for h in HORIZONS if h in sep.index]
        betas = [float(sep.loc[h, "beta"]) for h in hs]
        ax.plot(hs, betas, color=COLOR_D, lw=1.5, ls="--",
                marker="s", ms=3, label="LD separate", alpha=0.8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_title(f"Orth LD IRF vs separate — {outcome_label}", fontsize=10)
    ax.set_xlabel("Horizon h (quarters)", fontsize=9)
    ax.set_ylabel(ylabel_irf, fontsize=9)
    ax.legend(fontsize=7.5, framealpha=0.85)
    ax.grid(axis="y", lw=0.4, alpha=0.4)

    # Col 2 — Partial R²
    ax = axes_row[2]
    ax.plot(h_vals, pf_df["partial_r2_delta"], color=COLOR_D, lw=2,
            marker="o", ms=3.5, label=r"Partial R²: $\delta$")
    ax.plot(h_vals, pf_df["partial_r2_ld"], color=COLOR_LD, lw=2,
            marker="s", ms=3.5, label="Partial R²: LD")
    ax.set_title(f"Partial R² — {outcome_label}", fontsize=10)
    ax.set_xlabel("Horizon h (quarters)", fontsize=9)
    ax.set_ylabel("Partial R²", fontsize=9)
    ax.legend(fontsize=7.5, framealpha=0.85)
    ax.grid(axis="y", lw=0.4, alpha=0.4)
    ax.set_xticks(h_vals[::2])


def plot_diagnostics(pf_vac: pd.DataFrame, pf_unemp: pd.DataFrame,
                     orth_vac: pd.DataFrame, orth_unemp: pd.DataFrame,
                     sep_vac_csv: Path, sep_unemp_csv: Path,
                     out_path: Path) -> None:
    """
    2 × 3 figure:
      Row 0 (top):    vacancy outcome
      Row 1 (bottom): unemployment outcome
      Cols: Partial F | Orth LD IRF vs separate | Partial R²
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    _plot_row(axes[0], pf_vac, orth_vac, sep_vac_csv,
              outcome_label="vacancy",
              ylabel_irf="pp change in vacancy rate\nper 1-SD shock")

    _plot_row(axes[1], pf_unemp, orth_unemp, sep_unemp_csv,
              outcome_label="unemployment",
              ylabel_irf="pp change in unemp. rate\nper 1-SD shock")

    fig.suptitle(
        "Joint LP identification diagnostics: Story 1 (true zero) vs Story 2 (underpowered)\n"
        "Top row = vacancy outcome | Bottom row = unemployment outcome",
        fontsize=12
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Joint LP Identification Diagnostics")
print("=" * 60)

# [1] Load and build panel
print("\n[1] Loading data and building panel")
panel, sd_delta, sd_ld = load_and_build_panel()

# [2] Orthogonalize LD on delta (for Test 2)
print("\n[2] Orthogonalizing B^LD on B^delta (pooled, with FEs)")
panel = orthogonalize_ld(panel)

# [3] Test 1: Partial R² and Partial F — both outcomes
print("\n[3] Test 1: Partial R² and Partial F")
partial_f_vac   = test1_partial_f(panel, outcome="vacancy")
partial_f_unemp = test1_partial_f(panel, outcome="unemp")

# [4] Test 2: Orthogonalized LD LP — both outcomes
print("\n[4] Test 2: Orthogonalized LD LP")
orth_vac   = run_orth_lp(panel, outcome="vacancy")
orth_unemp = run_orth_lp(panel, outcome="unemp")

# [5] Save results
print("\n[5] Saving results")
for df, fname in [
    (partial_f_vac,   "joint_diagnostics_partial_f_vac.csv"),
    (partial_f_unemp, "joint_diagnostics_partial_f_unemp.csv"),
    (orth_vac,        "joint_diagnostics_orth_vacancy.csv"),
    (orth_unemp,      "joint_diagnostics_orth_unemp.csv"),
]:
    if not df.empty:
        p = RESULTS_DIR / fname
        df.to_csv(p, index=False)
        print(f"  Saved: {p}")

# [6] Plot — 2×3 grid: top row = vacancy, bottom row = unemployment
print("\n[6] Plotting")
plot_diagnostics(
    partial_f_vac, partial_f_unemp,
    orth_vac, orth_unemp,
    SEP_LD_VAC, SEP_LD_UNEMP,
    RESULTS_DIR / "joint_diagnostics.png"
)

# [7] Summary interpretation
print("\n[7] Summary")
for label, pf, orth in [
    ("VACANCY",      partial_f_vac,   orth_vac),
    ("UNEMPLOYMENT", partial_f_unemp, orth_unemp),
]:
    print(f"\n  --- {label} ---")
    if not pf.empty:
        avg_pf_ld = pf.loc[pf["h"] <= 8, "partial_f_ld"].mean()
        avg_pf_d  = pf.loc[pf["h"] <= 8, "partial_f_delta"].mean()
        print(f"  Avg partial F (h=0..8):  delta = {avg_pf_d:.1f},  LD = {avg_pf_ld:.1f}")
        if avg_pf_ld > 10:
            print("  -> LD has sufficient independent variation: zero is INFORMATIVE")
            print("     (supports Story 1: genuine zero effect)")
        elif avg_pf_ld > 3:
            print("  -> LD has moderate independent variation: zero is AMBIGUOUS")
        else:
            print("  -> LD has weak independent variation: zero is UNINFORMATIVE")
            print("     (supports Story 2: identification failure)")

    if not orth.empty:
        n_sig  = (orth["pval"] < 0.10).sum()
        mean_b = orth.loc[orth["h"] <= 8, "beta"].mean()
        print(f"  Orth LD {label.lower()}: {n_sig}/21 horizons significant at 10%")
        print(f"  Mean beta (h=0..8) = {mean_b:.4f}")
        if n_sig <= 2 and abs(mean_b) < 0.05:
            print("  -> Orthogonalized LD shows no effect: supports Story 1")
        else:
            print("  -> Orthogonalized LD shows significant effect: supports Story 2")

print("\nDone.")
