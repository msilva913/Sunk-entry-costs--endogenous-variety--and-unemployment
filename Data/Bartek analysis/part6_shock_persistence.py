"""
part6_shock_persistence.py — Aggregate AR(1) persistence of δ and s shocks
===========================================================================
Estimates AR(1) persistence for the δ (establishment closings) and LD
(layoffs/discharges) shock processes at the aggregate national level, using
both raw and residualized shock rates.

Why residualized rates matter for calibration
---------------------------------------------
- The structural model treats δ and s as AR(1) processes whose innovations are
primitive disturbances. 
- The raw shock rates g^δ and g^{LD} conflate this
structural variation with demand and productivity shocks that are explicitly
purged in the Bartik instrument construction. 
- Estimating AR(1) on raw rates therefore overstates structural persistence: if a demand contraction is itself
persistent, g^δ will be persistent partly because demand is persistent, not
because the structural process is. 
- The residualized rates ν^δ and ν^{LD}
(from part2b) remove aggregate productivity growth, lagged industry VA growth,
and lagged market tightness before Bartik aggregation — they are the correct
object for calibrating the model's AR(1) shock processes.

We therefore report two sets of estimates:

  Raw:          AR(1) on ḡ_t^(k) = Σ_j ω_j · g^(k)_{j,t}
                (employment-weighted average of raw LOO shock rates)

  Residualized: AR(1) on ν̄_t^(k) = Σ_j ω_j · ν^(k)_{j,t}
                (employment-weighted average of residualized shock rates)

The residualized series uses a common 2001Q1+ window for both δ and LD,
matching the JOLTS start date for LD and avoiding mechanical AR(1)
inflation from the low-variance pre-2001 expansion years in the δ series.
The extended δ sample (1997Q1+) is used only in the part5 LP robustness check.

Innovation covariance across shocks
------------------------------------
Coles and Kelishomi (2018, AEJ Macro) allow the innovations to their
productivity and separation shocks to be correlated. The analogous object
here is corr(η^δ_t, η^{LD}_t) — the contemporaneous correlation of AR(1)
residuals across the two shock processes. Even if we do not model cross-
persistence dynamics (which have unclear structural interpretation, per the
discussion in part2d), allowing correlated innovations is both empirically
motivated and straightforward to incorporate in the model.

We estimate this innovation correlation in both the raw and residualized
specifications and report it as a calibration target alongside ρ̂.

Aggregate shock construction
-----------------------------
    ḡ_t^(k)  = Σ_j ω_j · ḡ_{j,t}^(k)   (raw)
    ν̄_t^(k)  = Σ_j ω_j · ν^(k)_{j,t}   (residualized)

where ω_j = national employment share of industry j at base year 2006.
Industry-level rates are obtained by averaging LOO rates across states (the
state-specific LOO adjustment is small and averages out in the aggregate).
Residualized series ν^(k)_{j,t} are built in-process from cache parquets
using the same estimator as part2b, to avoid NTFS-mount parquet corruption.

AR(1) specification
--------------------
    x_t^(k) = μ^(k) + ρ^(k) · x^(k)_{t-1} + η^(k)_t

Estimated by OLS with HC3 robust standard errors. #refined version of Huber-White

Output
------
    data/results/shock_persistence.csv      — ρ̂ table (raw + residualized)
    data/results/shock_persistence_raw.png  — raw rates: time series + scatter + IRF
    data/results/shock_persistence_resid.png — residualized rates: same layout

Run
---
    python part6_shock_persistence.py

Prerequisites
-------------
    python part2_shock_rates.py
    python part2_shock_rates_s.py
    python part2b_shock_comovement.py   (populates BEA VA cache)
"""

import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
     os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

from construct_delta_instrument import (
    DEFAULT_CACHE_DIR, DEFAULT_OUTPUT_DIR, SHOCK_RATES_PATH, INDUSTRY_LABELS,
)
from construct_s_instrument import SHOCK_RATES_LD_PATH

RESID_D_PATH  = DEFAULT_OUTPUT_DIR / "shock_rates_delta_resid.parquet"
RESID_LD_PATH = DEFAULT_OUTPUT_DIR / "shock_rates_ld_resid.parquet"

# ── parameters ────────────────────────────────────────────────────────────────
BASE_YEAR   = 2006
IRF_H       = 20
RESULTS_DIR = Path("data/results")
INSTR_DIR   = Path("data/instruments")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SHARES_FILE = INSTR_DIR / f"shares_base{BASE_YEAR}.parquet"


# ═══════════════════════════════════════════════════════════════════════════════
#  DATA HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _load_weights() -> pd.Series:
    """
    National employment weights ω_j from base-year shares parquet.
    Returns a Series indexed by industry_code, summing to 1.
    """
    shares  = pd.read_parquet(SHARES_FILE)
    nat_emp = shares.groupby("industry_code")["emp_state_ind"].sum()
    return nat_emp / nat_emp.sum()


def _nat_avg(path, col):
    """Average LOO rates across states → industry × quarter DataFrame."""
    df = pd.read_parquet(path)
    return (df.groupby(["industry_code", "quarter_label"])[col]
              .mean().reset_index())


def _weighted_agg(ind_qt: pd.DataFrame, rate_col: str,
                  weights: pd.Series) -> pd.Series:
    """
    Employment-weighted average across industries within each quarter.

    ḡ_t = Σ_j ω_j · g_{j,t}

    Returns a Series indexed by quarter_label, sorted chronologically.
    """
    w_df = weights.rename("weight").reset_index()
    w_df["industry_code"] = w_df["industry_code"].astype(
        ind_qt["industry_code"].dtype)
    merged = ind_qt.merge(w_df, on="industry_code", how="inner")
    agg = (merged.groupby("quarter_label")
                 .apply(lambda g: np.average(g[rate_col], weights=g["weight"]),
                        include_groups=False)
                 .rename("agg_rate")
                 .sort_index())
    return agg


def _load_resid_series(parquet_path: Path, col: str,
                       weights: pd.Series) -> pd.Series | None:
    """
    Load a residualized shock series from a part2b output parquet and return
    its employment-weighted aggregate time series.

    Each series is loaded independently so δ (1997Q1+) and LD (2001Q1+) keep
    their own sample ranges.  The old _build_residualized_panel() merged both
    into a single inner-joined DataFrame, clipping δ to LD's shorter JOLTS
    sample and then further to the 2005Q4 FRED-VA start — producing N=65 for
    both.  Loading separately gives the correct N for each shock.

    Returns None if the parquet does not exist.
    """
    if not parquet_path.exists():
        return None
    df = pd.read_parquet(parquet_path)
    df = df.rename(columns={col: "rate"})
    return _weighted_agg(df[["industry_code", "quarter_label", "rate"]],
                         "rate", weights)


# ═══════════════════════════════════════════════════════════════════════════════
#  AR(1) ESTIMATION
# ═══════════════════════════════════════════════════════════════════════════════

def estimate_ar1(series: pd.Series, label: str) -> dict:
    """
    Estimate x_t = μ + ρ · x_{t-1} + η_t by OLS with HC3 robust SEs.

    Returns dict with ρ̂, SE, 95% CI, half-life, residuals, and fit objects.
    """
    df    = pd.DataFrame({"y": series, "y_lag": series.shift(1)}).dropna()
    X     = sm.add_constant(df["y_lag"], has_constant="add")
    model = sm.OLS(df["y"], X).fit(cov_type="HC3")

    mu        = float(model.params["const"])
    rho       = float(model.params["y_lag"])
    se        = float(model.bse["y_lag"])
    ci        = model.conf_int(alpha=0.05).loc["y_lag"].values
    uncond    = mu / (1 - rho) if abs(1 - rho) > 1e-6 else np.nan
    # half life: (1/2)=phi^h
    half_life = np.log(0.5) / np.log(rho) if 0 < rho < 1 else np.nan

    return {
        "label":       label,
        "rho":         rho,
        "se":          se,
        "ci95_lo":     ci[0],
        "ci95_hi":     ci[1],
        "mu":          mu,
        "uncond_mean": uncond,
        "half_life":   half_life,
        "r2":          float(model.rsquared),
        "nobs":        int(model.nobs),
        "fitted":      model.fittedvalues,
        "resid":       model.resid,
        "y":           df["y"],
        "y_lag":       df["y_lag"],
    }


def innovation_correlation(res_a: dict, res_b: dict) -> float:
    """
    Pearson correlation between AR(1) residuals η^a and η^b,
    aligned on their common quarters.

    This is corr(η^δ_t, η^{LD}_t) — the contemporaneous innovation
    correlation across shock processes.  Even without cross-persistence
    dynamics in A_1, allowing this to be nonzero (as Coles and Kelishomi
    2018 do for their productivity-separation shock pair) is empirically
    motivated and easy to incorporate in the model's shock calibration.
    """
    df = pd.DataFrame({"a": res_a["resid"], "b": res_b["resid"]}).dropna()
    return float(df["a"].corr(df["b"]))


def ar1_irf(rho: float, H: int = IRF_H) -> np.ndarray:
    """Theoretical shock-rate IRF to a unit innovation: IRF_h = ρ^h."""
    return np.array([rho ** h for h in range(H + 1)])


# ═══════════════════════════════════════════════════════════════════════════════
#  PRINTING
# ═══════════════════════════════════════════════════════════════════════════════

def _print_ar1(res: dict) -> None:
    print(f"    ρ̂  = {res['rho']:.4f}  SE={res['se']:.4f}  "
          f"95% CI=[{res['ci95_lo']:.4f}, {res['ci95_hi']:.4f}]")
    print(f"    Half-life = {res['half_life']:.2f} qtrs  "
          f"R² = {res['r2']:.4f}  N = {res['nobs']}")


def _print_comparison(label_a: str, raw: dict, resid: dict) -> None:
    """Side-by-side comparison of raw vs. residualized AR(1) for one shock."""
    print(f"\n  {label_a}:")
    print(f"    {'':30s}  {'Raw':>10}  {'Residualized':>14}")
    print(f"    {'ρ̂':30s}  {raw['rho']:>10.4f}  {resid['rho']:>14.4f}")
    print(f"    {'SE':30s}  {raw['se']:>10.4f}  {resid['se']:>14.4f}")
    print(f"    {'Half-life (qtrs)':30s}  {raw['half_life']:>10.2f}  "
          f"{resid['half_life']:>14.2f}")
    print(f"    {'N':30s}  {raw['nobs']:>10d}  {resid['nobs']:>14d}")
    print(f"    {'Sample start':30s}  {raw['y'].index[0]:>10}  "
          f"{resid['y'].index[0]:>14}")
    bias = raw['rho'] - resid['rho']
    print(f"    {'Δρ (raw − resid)':30s}  {bias:>+10.4f}"
          f"  ← demand/productivity bias in raw ρ̂")


# ═══════════════════════════════════════════════════════════════════════════════
#  PLOTTING
# ═══════════════════════════════════════════════════════════════════════════════

REC_SPANS = [("1990-07-01", "1991-03-01"),
             ("2001-03-01", "2001-11-01"),
             ("2007-12-01", "2009-06-01")]


def _plot_persistence(results: dict, spec_label: str, out_path: Path) -> None:
    """
    Two-row figure (δ top, LD bottom), three panels each:
      Left:   Aggregate shock rate time series + AR(1) fitted values
      Middle: Scatter g_t vs g_{t-1} with OLS line
      Right:  Theoretical AR(1) IRF given ρ̂
    """
    shocks = list(results.keys())
    colors = {"δ": "#1f77b4", "LD": "#d62728"}
    fig, axes = plt.subplots(len(shocks), 3, figsize=(15, 4.5 * len(shocks)))
    if len(shocks) == 1:
        axes = axes[np.newaxis, :]

    for row, lbl in enumerate(shocks):
        res   = results[lbl]
        color = colors.get(lbl, "#2ca02c")
        irf   = ar1_irf(res["rho"], IRF_H)
        h_arr = np.arange(IRF_H + 1)

        # Panel 1: time series
        ax1  = axes[row, 0]
        qs   = res["y"].index.tolist()
        dts  = [pd.Period(q, freq="Q").to_timestamp() for q in qs]
        fdts = [pd.Period(q, freq="Q").to_timestamp()
                for q in res["fitted"].index]
        ax1.plot(dts, res["y"].values * 100, color=color,
                 lw=1.2, alpha=0.8, label="Observed")
        ax1.plot(fdts, res["fitted"].values * 100,
                 color="black", lw=1.0, ls="--", alpha=0.7, label="AR(1) fit")
        for rs, re in REC_SPANS:
            ax1.axvspan(pd.Timestamp(rs), pd.Timestamp(re),
                        alpha=0.10, color="grey")
        ax1.set_title(f"{lbl} — {spec_label} rate", fontsize=10)
        ax1.set_ylabel("Shock rate (×100)", fontsize=9)
        ax1.legend(fontsize=8)
        ax1.grid(axis="y", lw=0.4, alpha=0.4)

        # Panel 2: scatter
        ax2  = axes[row, 1]
        x_sc = res["y_lag"].values * 100
        y_sc = res["y"].values * 100
        ax2.scatter(x_sc, y_sc, color=color, alpha=0.4, s=12)
        x_line = np.linspace(x_sc.min(), x_sc.max(), 100)
        ax2.plot(x_line, res["mu"] * 100 + res["rho"] * x_line,
                 color="black", lw=1.2, label=f"ρ̂ = {res['rho']:.3f}")
        ax2.set_xlabel(r"$x_{t-1}$ (×100)", fontsize=9)
        ax2.set_ylabel(r"$x_t$ (×100)", fontsize=9)
        ax2.set_title(f"{lbl} — AR(1) scatter", fontsize=10)
        ax2.legend(fontsize=9)
        ax2.grid(lw=0.4, alpha=0.4)

        # Panel 3: IRF
        ax3 = axes[row, 2]
        ax3.plot(h_arr, irf, color=color, lw=2.0, marker="o", ms=3.5,
                 label=f"ρ̂ = {res['rho']:.3f}")
        ax3.axhline(0,   color="black", lw=0.8)
        ax3.axhline(0.5, color="grey",  lw=0.6, ls=":",
                    label=f"Half-life ≈ {res['half_life']:.1f} qtrs")
        ax3.set_title(f"{lbl} — theoretical IRF  (ρ̂ʰ)", fontsize=10)
        ax3.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax3.set_ylabel("IRF = ρ̂ʰ", fontsize=9)
        ax3.set_xticks(h_arr)
        ax3.legend(fontsize=9)
        ax3.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle(
        f"Aggregate AR(1) shock persistence — {spec_label} rates\n"
        f"(base-year {BASE_YEAR} employment weights; OLS with HC3 SEs)",
        fontsize=12, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


def _plot_irf_comparison(raw_results: dict, resid_results: dict,
                         out_path: Path) -> None:
    """
    Overlay raw vs. residualized theoretical IRFs for each shock.
    One panel per shock (δ left, LD right).
    """
    shocks = list(raw_results.keys())
    colors_raw   = {"δ": "#1f77b4", "LD": "#d62728"}
    colors_resid = {"δ": "#aec7e8", "LD": "#f7b6b2"}
    h_arr = np.arange(IRF_H + 1)

    fig, axes = plt.subplots(1, len(shocks), figsize=(7 * len(shocks), 4))
    if len(shocks) == 1:
        axes = [axes]

    for ax, lbl in zip(axes, shocks):
        rho_r  = raw_results[lbl]["rho"]
        rho_re = resid_results[lbl]["rho"]
        ax.plot(h_arr, ar1_irf(rho_r,  IRF_H), color=colors_raw[lbl],
                lw=2, label=f"Raw  ρ̂={rho_r:.3f}")
        ax.plot(h_arr, ar1_irf(rho_re, IRF_H), color=colors_resid[lbl],
                lw=2, ls="--", label=f"Residualized  ρ̂={rho_re:.3f}")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"{lbl} shock — IRF comparison", fontsize=10)
        ax.set_xlabel("Horizon h (quarters)", fontsize=9)
        ax.set_ylabel("IRF = ρ̂ʰ", fontsize=9)
        ax.legend(fontsize=9)
        ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle(
        "Raw vs. residualized AR(1) IRF comparison\n"
        "(gap = demand/productivity contamination in raw ρ̂)",
        fontsize=11
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Comparison plot saved: {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 65)
    print("Part 6 — Aggregate Shock Persistence (AR(1))")
    print("=" * 65)

    # ── 1. Employment weights ─────────────────────────────────────────────────
    print("\n[1] Loading employment weights")
    weights = _load_weights()
    print(f"  {len(weights)} industries, weight sum={weights.sum():.6f}")

    # ── 2. Raw aggregate shock rates ──────────────────────────────────────────
    print("\n[2] Raw aggregate shock rates")
    raw_d  = _weighted_agg(_nat_avg(SHOCK_RATES_PATH,    "g_delta_loo")
                            .rename(columns={"g_delta_loo": "rate"}),
                           "rate", weights)
    raw_ld = _weighted_agg(_nat_avg(SHOCK_RATES_LD_PATH, "g_ld_loo")
                            .rename(columns={"g_ld_loo": "rate"}),
                           "rate", weights)
    # _nat_avg returns a 3-column df; rename for _weighted_agg
    def _agg(path, col):
        df = _nat_avg(path, col)
        df = df.rename(columns={col: "rate"})
        return _weighted_agg(df, "rate", weights)

    agg_raw_d  = _agg(SHOCK_RATES_PATH,    "g_delta_loo")
    agg_raw_ld = _agg(SHOCK_RATES_LD_PATH, "g_ld_loo")
    print(f"  δ:  {agg_raw_d.index[0]}–{agg_raw_d.index[-1]},  "
          f"mean={agg_raw_d.mean():.5f},  SD={agg_raw_d.std():.5f}")
    print(f"  LD: {agg_raw_ld.index[0]}–{agg_raw_ld.index[-1]},  "
          f"mean={agg_raw_ld.mean():.5f},  SD={agg_raw_ld.std():.5f}")

    # ── 3. Residualized aggregate shock rates ─────────────────────────────────
    print("\n[3] Residualized aggregate shock rates (part2b outputs)")
    agg_resid_d  = _load_resid_series(RESID_D_PATH,  "nu_delta", weights)
    agg_resid_ld = _load_resid_series(RESID_LD_PATH, "nu_ld",    weights)

    if agg_resid_d is None or agg_resid_ld is None:
        missing = []
        if agg_resid_d  is None: missing.append(str(RESID_D_PATH))
        if agg_resid_ld is None: missing.append(str(RESID_LD_PATH))
        print(f"  WARNING: residualized parquets not found:\n"
              + "\n".join(f"    {p}" for p in missing))
        print("  Run part2b_residualize_shocks.py first.")
        has_resid = False
    else:
        # Restrict both residualized series to 2001Q1+ for a common window.
        # ν^δ covers 1997Q1+ after the part2b panel split, but the pre-2001
        # expansion years (low, stable exit rates) mechanically inflate the
        # AR(1) coefficient because consecutive quarters of near-identical
        # values look highly autocorrelated.  Using 2001Q1+ for both δ and LD
        # gives a comparable sample and avoids this bias.  The extended δ
        # series (1997Q1+) is used only in the part5 sample-extension LP.
        RESID_START = "2001Q1"
        agg_resid_d  = agg_resid_d[agg_resid_d.index  >= RESID_START]
        agg_resid_ld = agg_resid_ld[agg_resid_ld.index >= RESID_START]
        has_resid = True
        print(f"  Common window: {RESID_START}+ (δ extended series restricted "
              f"to match LD JOLTS start; avoids pre-2001 expansion bias)")
        print(f"  ν^δ:  {agg_resid_d.index[0]}–{agg_resid_d.index[-1]},  "
              f"mean={agg_resid_d.mean():.5f},  SD={agg_resid_d.std():.5f}")
        print(f"  ν^LD: {agg_resid_ld.index[0]}–{agg_resid_ld.index[-1]},  "
              f"mean={agg_resid_ld.mean():.5f},  SD={agg_resid_ld.std():.5f}")

    # ── 4. AR(1) estimation ───────────────────────────────────────────────────
    print("\n[4] AR(1) estimation")
    print("\n  --- Raw rates ---")
    res_raw_d  = estimate_ar1(agg_raw_d,  "δ  (raw)")
    res_raw_ld = estimate_ar1(agg_raw_ld, "LD (raw)")
    _print_ar1(res_raw_d)
    _print_ar1(res_raw_ld)

    corr_raw = innovation_correlation(res_raw_d, res_raw_ld)
    print(f"\n  Innovation correlation (raw):")
    print(f"    corr(η^δ, η^{{LD}}) = {corr_raw:+.4f}")
    print(f"    Model assumption corr=0 is "
          f"{'approximately satisfied' if abs(corr_raw) < 0.15 else 'violated'}.")

    if has_resid:
        print("\n  --- Residualized rates ---")
        res_re_d  = estimate_ar1(agg_resid_d,  "δ  (residualized)")
        res_re_ld = estimate_ar1(agg_resid_ld, "LD (residualized)")
        _print_ar1(res_re_d)
        _print_ar1(res_re_ld)

        corr_resid = innovation_correlation(res_re_d, res_re_ld)
        print(f"\n  Innovation correlation (residualized):")
        print(f"    corr(η^δ, η^{{LD}}) = {corr_resid:+.4f}")
        print(f"    Model assumption corr=0 is "
              f"{'approximately satisfied' if abs(corr_resid) < 0.15 else 'violated'}.")

    # ── 5. Comparison table ───────────────────────────────────────────────────
    print("\n[5] Raw vs. residualized comparison")
    if has_resid:
        _print_comparison("δ",  res_raw_d,  res_re_d)
        _print_comparison("LD", res_raw_ld, res_re_ld)

        print(f"\n  Innovation correlations:")
        print(f"    Raw:          corr(η^δ, η^{{LD}}) = {corr_raw:+.4f}")
        print(f"    Residualized: corr(η^δ, η^{{LD}}) = {corr_resid:+.4f}")
        print(f"\n  Interpretation:")
        print(f"    The gap Δρ^δ  = {res_raw_d['rho'] - res_re_d['rho']:+.4f} "
              f"is the upward bias in raw ρ̂^δ from demand/productivity.")
        print(f"    The gap Δρ^LD = {res_raw_ld['rho'] - res_re_ld['rho']:+.4f} "
              f"is the analogous bias for LD.")
        print(f"    The residualized estimates are the preferred calibration "
              f"targets for the model's AR(1) shock processes.")
        print(f"    The innovation correlation — even after residualization —")
        print(f"    motivates allowing corr(ε^δ, ε^s) ≠ 0 in the model,")
        print(f"    consistent with Coles and Kelishomi (2018).")

    # ── 6. Save results CSV ───────────────────────────────────────────────────
    print("\n[6] Saving results")
    rows = []
    for res, spec in [(res_raw_d, "raw"), (res_raw_ld, "raw")]:
        rows.append({k: v for k, v in res.items()
                     if k not in ("fitted", "resid", "y", "y_lag")})
        rows[-1]["spec"] = spec

    if has_resid:
        for res, spec in [(res_re_d, "residualized"), (res_re_ld, "residualized")]:
            rows.append({k: v for k, v in res.items()
                         if k not in ("fitted", "resid", "y", "y_lag")})
            rows[-1]["spec"] = spec

        rows.append({"label": "innovation_corr_raw",   "rho": corr_raw,   "spec": "raw"})
        rows.append({"label": "innovation_corr_resid", "rho": corr_resid, "spec": "residualized"})

    pd.DataFrame(rows).to_csv(RESULTS_DIR / "shock_persistence.csv", index=False)
    print(f"  Saved: {RESULTS_DIR / 'shock_persistence.csv'}")

    # ── 7. Plots ──────────────────────────────────────────────────────────────
    print("\n[7] Plotting")
    _plot_persistence(
        {"δ": res_raw_d, "LD": res_raw_ld},
        spec_label="raw",
        out_path=RESULTS_DIR / "shock_persistence_raw.png"
    )
    if has_resid:
        _plot_persistence(
            {"δ": res_re_d, "LD": res_re_ld},
            spec_label="residualized",
            out_path=RESULTS_DIR / "shock_persistence_resid.png"
        )
        _plot_irf_comparison(
            {"δ": res_raw_d,  "LD": res_raw_ld},
            {"δ": res_re_d,   "LD": res_re_ld},
            out_path=RESULTS_DIR / "shock_persistence_comparison.png"
        )

    # ── 8. Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("CALIBRATION SUMMARY")
    print("=" * 65)
    print(f"\n  {'Shock':<6}  {'Spec':>14}  {'ρ̂':>8}  {'Half-life':>10}  {'N':>6}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*8}  {'-'*10}  {'-'*6}")
    for res, spec in [(res_raw_d, "raw"), (res_raw_ld, "raw")]:
        print(f"  {res['label'].split()[0]:<6}  {spec:>14}  "
              f"{res['rho']:>8.4f}  {res['half_life']:>10.2f}  {res['nobs']:>6}")
    if has_resid:
        for res, spec in [(res_re_d, "residualized"), (res_re_ld, "residualized")]:
            print(f"  {res['label'].split()[0]:<6}  {spec:>14}  "
                  f"{res['rho']:>8.4f}  {res['half_life']:>10.2f}  {res['nobs']:>6}")
        print(f"\n  Innovation correlations:")
        print(f"    Raw:          {corr_raw:+.4f}")
        print(f"    Residualized: {corr_resid:+.4f}")
    print(f"\n  Preferred calibration targets (residualized):")
    if has_resid:
        print(f"    ρ_δ  = {res_re_d['rho']:.4f}  "
              f"(raw: {res_raw_d['rho']:.4f}, bias={res_raw_d['rho']-res_re_d['rho']:+.4f})")
        print(f"    ρ_LD = {res_re_ld['rho']:.4f}  "
              f"(raw: {res_raw_ld['rho']:.4f}, bias={res_raw_ld['rho']-res_re_ld['rho']:+.4f})")
        print(f"    corr(ε^δ, ε^LD) = {corr_resid:+.4f}")
    else:
        print(f"    ρ_δ  = {res_raw_d['rho']:.4f}  (raw only — residualized unavailable)")
        print(f"    ρ_LD = {res_raw_ld['rho']:.4f}  (raw only — residualized unavailable)")

    print("\nDone.")


if __name__ == "__main__":
    main()

