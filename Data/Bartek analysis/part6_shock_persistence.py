"""
part6_shock_persistence.py — Aggregate AR(1) persistence of δ and s shocks
===========================================================================
Estimates the persistence of the aggregate national job-destruction (δ) and
separation (s) shock rates, which the structural model treats as AR(1)
processes. The estimated AR(1) coefficients discipline interpretation of the
local projection IRFs in part5: a highly persistent shock should produce a
slowly-reverting IRF, so the empirical IRF shape can be compared against the
model-predicted IRF given these parameters.

Aggregate shock construction
-----------------------------
The aggregate national shock rate in quarter t is the national employment-
weighted average of the industry-level LOO shock rates:

    ḡ_t^(k) = Σ_j ω_j · ḡ_{j,t}^(k)

where ω_j = emp_nat_j / Σ_j emp_nat_j  (national employment share of
industry j at base year), and ḡ_{j,t}^(k) is the national (non-LOO) shock
rate for industry j in quarter t.

Note on LOO vs. national rates
--------------------------------
The LOO rates g_{-s,j,t} used in the Bartik instrument are state-specific.
For the aggregate time series we want the pure national rate — the simple
employment-weighted average across industries with no LOO adjustment. We
recover this by averaging the LOO rates across states (the state-specific
adjustment is small and averages out), or equivalently by weighting the
industry-level shock rates directly.

AR(1) specification
--------------------
    ḡ_t^(k) = μ^(k) + ρ^(k) · ḡ_{t-1}^(k) + η_t^(k)

Estimated by OLS with HC3 robust standard errors (time series, n~80-120).
The intercept μ = (1 - ρ) · ḡ̄ pins the unconditional mean.

Impulse response
-----------------
Given ρ̂, the theoretical IRF of the shock rate itself to a unit innovation
at h=0 is simply:

    IRF_h^(k) = ρ̂^h,  h = 0, 1, ..., H

This is plotted alongside the OLS fit and the raw time series.

Outputs
-------
    data/results/shock_persistence.csv    — ρ̂, SE, CI, unconditional mean
    data/results/shock_persistence.png    — time series + AR(1) fit + IRF

Run
---
    python part6_shock_persistence.py

Prerequisites
-------------
    python part1_shares.py
    python part2_shock_rates.py       (δ shock rates)
    python part2_shock_rates_s.py     (s shock rates)
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

def _open_file(path):
    os.startfile(path)

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents"
        / "GitHub"
        / "Sunk-entry-costs--endogenous-variety--and-unemployment"
        / "Data"
        / "Bartek analysis"
    )

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_YEAR  = 2006
IRF_H      = 20       # horizons for theoretical AR(1) IRF plot

INSTR_DIR  = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SHARES_FILE      = INSTR_DIR / f"shares_base{BASE_YEAR}.parquet"
DELTA_RATES_FILE = INSTR_DIR / "shock_rates_1992Q3_2023Q1.parquet"
S_RATES_FILE     = INSTR_DIR / "shock_rates_s_2001Q1_2023Q1.parquet"


# ---------------------------------------------------------------------------
# Step 1: Load shares — national employment weights ω_j
# ---------------------------------------------------------------------------
def load_national_weights(shares_path: Path) -> pd.Series:
    """
    Compute national industry employment weights from the shares parquet.
    Returns a Series indexed by industry_code, summing to 1.
    """
    shares = pd.read_parquet(shares_path)
    # National employment by industry: sum emp_state_ind across all states
    nat_emp = (
        shares.groupby("industry_code")["emp_state_ind"]
        .sum()
    )
    weights = nat_emp / nat_emp.sum()
    print(f"  National weights: {len(weights)} industries, sum={weights.sum():.6f}")
    return weights


# ---------------------------------------------------------------------------
# Step 2: Construct aggregate time series ḡ_t^(k)
# ---------------------------------------------------------------------------
def build_aggregate_shock(rates: pd.DataFrame,
                          rate_col: str,
                          weights: pd.Series) -> pd.Series:
    """
    Construct the aggregate national shock rate as the employment-weighted
    average of industry-level LOO rates across states.

    For the aggregate time series we average LOO rates across states first
    (recovering approximately the national rate), then weight by industry.

    Returns a Series indexed by quarter_label, sorted chronologically.
    """
    # Average LOO rates across states within each (industry, quarter)
    # — the state-specific LOO adjustment is small and averages out
    nat_by_ind_qt = (
        rates.groupby(["industry_code", "quarter_label"])[rate_col]
        .mean()
        .reset_index()
    )

    # Merge national weights
    nat_by_ind_qt = nat_by_ind_qt.merge(
        weights.rename("weight").reset_index(),
        on="industry_code",
        how="inner",
    )

    # Weighted average across industries within each quarter
    agg = (
        nat_by_ind_qt.groupby("quarter_label")
        .apply(lambda g: np.average(g[rate_col], weights=g["weight"]))
        .rename("agg_rate")
        .sort_index()
    )

    print(f"  Aggregate series: {len(agg)} quarters  "
          f"({agg.index[0]} – {agg.index[-1]})")
    print(f"  Mean={agg.mean():.5f}  SD={agg.std():.5f}  "
          f"Min={agg.min():.5f}  Max={agg.max():.5f}")
    return agg


# ---------------------------------------------------------------------------
# Step 3: AR(1) OLS estimation
# ---------------------------------------------------------------------------
def estimate_ar1(series: pd.Series, label: str) -> dict:
    """
    Estimate ḡ_t = μ + ρ · ḡ_{t-1} + η_t by OLS with HC3 robust SEs.

    Returns dict with ρ̂, SE, 95% CI, unconditional mean, half-life.
    """
    df = pd.DataFrame({"y": series, "y_lag": series.shift(1)}).dropna()

    X = sm.add_constant(df["y_lag"], has_constant="add")
    y = df["y"]

    model = sm.OLS(y, X).fit(cov_type="HC3")

    mu  = float(model.params["const"])
    rho = float(model.params["y_lag"])
    se  = float(model.bse["y_lag"])
    ci  = model.conf_int(alpha=0.05).loc["y_lag"].values

    # Unconditional mean = μ / (1 - ρ)
    uncond_mean = mu / (1 - rho) if abs(1 - rho) > 1e-6 else np.nan

    # Half-life: ρ^h = 0.5 → h = log(0.5) / log(ρ)
    half_life = np.log(0.5) / np.log(rho) if 0 < rho < 1 else np.nan

    print(f"\n  {label} AR(1):")
    print(f"    ρ̂  = {rho:.4f}  SE={se:.4f}  "
          f"95% CI=[{ci[0]:.4f}, {ci[1]:.4f}]")
    print(f"    μ  = {mu:.6f}")
    print(f"    Unconditional mean = {uncond_mean:.5f}")
    print(f"    Half-life          = {half_life:.1f} quarters")
    print(f"    R²                 = {model.rsquared:.4f}")
    print(f"    N                  = {int(model.nobs)}")

    return {
        "label":        label,
        "rho":          rho,
        "se":           se,
        "ci95_lo":      ci[0],
        "ci95_hi":      ci[1],
        "mu":           mu,
        "uncond_mean":  uncond_mean,
        "half_life":    half_life,
        "r2":           float(model.rsquared),
        "nobs":         int(model.nobs),
        "fitted":       model.fittedvalues,
        "resid":        model.resid,
        "y":            y,
        "y_lag":        df["y_lag"],
    }


# ---------------------------------------------------------------------------
# Step 4: Theoretical AR(1) IRF
# ---------------------------------------------------------------------------
def ar1_irf(rho: float, H: int = IRF_H) -> np.ndarray:
    """
    Theoretical IRF of the shock rate to a unit innovation at h=0.
    IRF_h = ρ^h for h = 0, 1, ..., H.
    """
    return np.array([rho ** h for h in range(H + 1)])


# ---------------------------------------------------------------------------
# Step 5: Plot
# ---------------------------------------------------------------------------
def plot_persistence(results: dict[str, dict], out_path: Path) -> None:
    """
    Three-panel figure for each shock:
      Left:   Raw aggregate shock rate time series + AR(1) fitted values
      Middle: Scatter of g_t vs g_{t-1} with OLS line
      Right:  Theoretical AR(1) IRF given ρ̂

    Two rows — one per shock (δ top, s bottom).
    """
    shocks = list(results.keys())
    n = len(shocks)
    fig, axes = plt.subplots(n, 3, figsize=(15, 4.5 * n))
    if n == 1:
        axes = axes[np.newaxis, :]

    colors = {"δ": "#1f77b4", "s": "#d62728"}

    for row, label in enumerate(shocks):
        res   = results[label]
        rho   = res["rho"]
        color = colors.get(label, "#2ca02c")
        irf   = ar1_irf(rho, IRF_H)
        h_arr = np.arange(IRF_H + 1)

        # --- Panel 1: Time series ---
        ax1 = axes[row, 0]
        # Convert quarter labels to approximate dates for x-axis
        qs  = res["y"].index.tolist()
        dts = [pd.Period(q, freq="Q").to_timestamp() for q in qs]
        fit_dts = [pd.Period(q, freq="Q").to_timestamp()
                   for q in res["fitted"].index]

        ax1.plot(dts, res["y"].values * 100, color=color,
                 linewidth=1.2, alpha=0.8, label="Observed")
        ax1.plot(fit_dts, res["fitted"].values * 100,
                 color="black", linewidth=1.0, linestyle="--",
                 alpha=0.7, label="AR(1) fit")

        # NBER recession shading
        for rs, re in [("1990-07-01","1991-03-01"),
                       ("2001-03-01","2001-11-01"),
                       ("2007-12-01","2009-06-01")]:
            ax1.axvspan(pd.Timestamp(rs), pd.Timestamp(re),
                        alpha=0.10, color="grey")

        ax1.set_title(f"{label} shock — aggregate rate over time", fontsize=10)
        ax1.set_ylabel("Shock rate (×100, %)", fontsize=9)
        ax1.legend(fontsize=8)
        ax1.grid(axis="y", linewidth=0.4, alpha=0.4)

        # --- Panel 2: Scatter g_t vs g_{t-1} ---
        ax2 = axes[row, 1]
        x_sc = res["y_lag"].values * 100
        y_sc = res["y"].values * 100
        ax2.scatter(x_sc, y_sc, color=color, alpha=0.4, s=12)

        # OLS line
        x_line = np.linspace(x_sc.min(), x_sc.max(), 100)
        mu_pct  = res["mu"] * 100
        rho_val = res["rho"]
        ax2.plot(x_line, mu_pct + rho_val * x_line,
                 color="black", linewidth=1.2,
                 label=f"ρ̂ = {rho_val:.3f}")
        ax2.set_xlabel(r"$\bar{g}_{t-1}$ (×100)", fontsize=9)
        ax2.set_ylabel(r"$\bar{g}_t$ (×100)", fontsize=9)
        ax2.set_title(f"{label} shock — AR(1) scatter", fontsize=10)
        ax2.legend(fontsize=9)
        ax2.grid(linewidth=0.4, alpha=0.4)

        # --- Panel 3: Theoretical IRF ---
        ax3 = axes[row, 2]
        ax3.plot(h_arr, irf, color=color, linewidth=2.0,
                 marker="o", markersize=3.5,
                 label=f"ρ̂ = {rho:.3f}")
        ax3.axhline(0, color="black", linewidth=0.8)
        ax3.axhline(0.5, color="grey", linewidth=0.6,
                    linestyle=":", alpha=0.7,
                    label=f"Half-life ≈ {res['half_life']:.1f} qtrs")
        for hh in [4, 8, 12, 16, 20]:
            ax3.axvline(hh, color="grey", linewidth=0.4,
                        linestyle=":", alpha=0.5)
        ax3.set_title(f"{label} shock — theoretical AR(1) IRF\n"
                      f"(response of shock rate to unit innovation)",
                      fontsize=10)
        ax3.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax3.set_ylabel("IRF = ρ̂ʰ", fontsize=9)
        ax3.set_xticks(h_arr)
        ax3.legend(fontsize=9)
        ax3.grid(axis="y", linewidth=0.4, alpha=0.4)

    # Caption
    caption = (
        "Notes: Left panels show the aggregate national shock rate "
        r"$\bar{g}_t^{(k)} = \sum_j \omega_j \bar{g}_{j,t}^{(k)}$"
        " (employment-weighted average across industries, ×100) "
        "with the AR(1) fitted values overlaid. Grey shading marks NBER recessions.  "
        "Middle panels show the AR(1) scatter of $\\bar{g}_t$ on $\\bar{g}_{t-1}$ "
        "with the OLS regression line; $\\hat{\\rho}$ is estimated by OLS with HC3 "
        "robust standard errors.  "
        "Right panels show the theoretical impulse response of the shock rate "
        "to a unit innovation at $h=0$: $\\mathrm{IRF}_h = \\hat{\\rho}^h$. "
        "The dashed grey line marks the half-life (IRF = 0.5). "
        "These persistence estimates discipline interpretation of the LP IRFs "
        "in part5: an empirical unemployment IRF that does not revert within "
        "16 quarters is consistent with a structural shock whose own IRF "
        "(right panel) also remains well above zero at that horizon."
    )
    fig.text(0.5, -0.03, caption, ha="center", va="top", fontsize=7.5,
             wrap=True, transform=fig.transFigure,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                       edgecolor="#cccccc", linewidth=0.8),
             multialignment="left")

    fig.suptitle(
        "Aggregate AR(1) shock persistence — δ and s shock rates\n"
        f"(base year {BASE_YEAR} employment weights; OLS with HC3 SEs)",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 6 — Aggregate Shock Persistence (AR(1))")
print("=" * 60)

# -----------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------
print("\n[1] Loading national employment weights")
weights = load_national_weights(SHARES_FILE)

print("\n[2] Loading shock rate panels")
delta_rates = pd.read_parquet(DELTA_RATES_FILE)
s_rates     = pd.read_parquet(S_RATES_FILE)
print(f"  δ rates: {delta_rates.shape}")
print(f"  s rates: {s_rates.shape}")

# -----------------------------------------------------------------------
# Build aggregate time series
# -----------------------------------------------------------------------
print("\n[3] Constructing aggregate shock rates")
print("  δ:")
agg_delta = build_aggregate_shock(delta_rates, "g_delta_loo", weights)
print("  s:")
agg_s     = build_aggregate_shock(s_rates,     "g_s_loo",     weights)

# -----------------------------------------------------------------------
# AR(1) estimation
# -----------------------------------------------------------------------
print("\n[4] Estimating AR(1) persistence")
res_delta = estimate_ar1(agg_delta, "δ")
res_s     = estimate_ar1(agg_s,     "s")

# -----------------------------------------------------------------------
# Save results table
# -----------------------------------------------------------------------
summary = pd.DataFrame([
    {k: v for k, v in res_delta.items()
     if k not in ("fitted","resid","y","y_lag")},
    {k: v for k, v in res_s.items()
     if k not in ("fitted","resid","y","y_lag")},
])
summary.to_csv(RESULTS_DIR / "shock_persistence.csv", index=False)
print(f"\n  Saved: {RESULTS_DIR / 'shock_persistence.csv'}")

# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------
print("\n[5] Plotting")
plot_persistence(
    {"δ": res_delta, "s": res_s},
    out_path=RESULTS_DIR / "shock_persistence.png",
)
_open_file(RESULTS_DIR / "shock_persistence.png")

# -----------------------------------------------------------------------
# Summary for part5 use
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print("Summary for LP interpretation (part5)")
print("=" * 60)
for res in [res_delta, res_s]:
    irf = ar1_irf(res["rho"], IRF_H)
    print(f"\n  {res['label']} shock:")
    print(f"    ρ̂ = {res['rho']:.4f}  (half-life = {res['half_life']:.1f} qtrs)")
    print(f"    Shock IRF at h=4:  {irf[4]:.3f}")
    print(f"    Shock IRF at h=8:  {irf[8]:.3f}")
    print(f"    Shock IRF at h=16: {irf[16]:.3f}")
    print(f"    → An unemployment IRF that has not reverted by h=16 is "
          f"{'consistent' if irf[16] > 0.1 else 'harder to explain'} "
          f"with ρ̂={res['rho']:.3f}")
