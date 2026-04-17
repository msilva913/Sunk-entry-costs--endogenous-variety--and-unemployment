"""
part7c_recession_placebo.py — Recession-severity placebo for δ, LD, QU
=======================================================================

Motivation
----------
The baseline LP in part5 estimates:

    y_{s,t+h} - y_{s,t-1} = α_s + α_t + β_h B^k_{s,t} + controls + ε

Time fixed effects α_t absorb the *level* of aggregate conditions in each
shock quarter t — they soak up the fact that all states are hit by the same
national recession in, say, 2008Q4.  But they cannot absorb a
*differential* effect: if states with high Bartik instrument values
systematically suffer more during bad aggregate quarters than during good
ones, β_h will pick up that composition-times-severity interaction rather
than the structural shock-transmission channel.

This is the central remaining identification threat after the GFC
interaction (part5) and SLOOS (part7b) checks.  The GFC interaction used a
binary dummy for 2008Q3–2009Q4, which tested one specific episode.  The
severity placebo uses a *continuous*, quarter-by-quarter measure of aggregate
conditions — the national unemployment rate change Δu^nat_t — that varies
across all three business cycles in the sample (2001 recession, GFC, and the
expansion years between them).

The intuition
-------------
Imagine two states, A and B.  State A has a high Manufacturing share (high
δ Bartik value); State B has a low Manufacturing share (low δ Bartik value).
During a recession, Manufacturing is hit hard and state A suffers more than
state B.  During expansions, the gap is smaller.  If this differential
cyclical sensitivity is what drives the unemployment gap between A and B —
not structural firm destruction — then β_h is a severity-composition
artifact, not a structural estimate.

The severity placebo adds B^k_{s,t} × Δu^nat_t as an additional regressor:

    y_{s,t+h} - y_{s,t-1} = α_s + α_t
                             + β_h  B^k_{s,t}                    ← structural channel
                             + δ_h  B^k_{s,t} × Δu^nat_t         ← severity interaction
                             + φ_h  Δu^nat_t                     ← severity level
                             + γ_1 y_{s,t-1} + γ_2 log LF_{s,t-1}
                             + ε_{s,t,h}

The level term φ_h Δu^nat_t is included so that δ_h identifies the
*differential* response of high-instrument states to a given unit of
aggregate severity, not any residual aggregate correlation that α_t might
miss at the LP horizon h (since α_t controls the shock quarter t, not the
outcome quarter t+h).

Reading the results
-------------------
β_h (conditional on average severity):
  The structural IRF estimate cleaned of severity-composition variation.
  This is the headline number.  Compare to the part5 baseline.

δ_h (severity interaction coefficient):
  The additional unemployment/vacancy response per unit of Δu^nat_t, for
  states with a one-unit-higher Bartik value.
  - δ_h ≈ 0, insignificant: the Bartik effect is constant across the
    business cycle → severity is NOT driving β_h → identification is clean.
  - δ_h > 0, significant: high-instrument states are hit harder during
    recessions → some severity confound in β_h → check whether β_h collapses.
  - β_h collapses to ≈ 0 when interaction included: severity was doing all
    the work → instrument validity is compromised.
  - β_h remains substantial alongside significant δ_h: genuine structural
    channel with amplification during downturns → consistent with
    state-dependent transmission, not a full instrument failure.

Why Δu^nat_t rather than NFCI
------------------------------
NFCI conflates credit conditions, financial-market volatility, and risk
premia — three channels with distinct structural interpretations.  A
significant NFCI interaction cannot be cleanly attributed to severity versus
a genuine credit amplification channel (consistent with the SLOOS results
suggesting a small credit effect for δ).

Δu^nat_t = national unemployment rate change (from FRED UNRATE, quarterly
average of monthly data) is a cleaner severity measure: there is no model in
which the structural effect of firm destruction on a state's unemployment
should depend on what the national unemployment rate happens to be doing.
A significant interaction with Δu^nat_t is difficult to rationalize as
anything other than a composition-severity confound.

We also report NFCI as a secondary severity measure for completeness and
comparability with prior results.

Interpretation for LD and QU
-----------------------------
For LD and QU the stakes are higher.  The model predicts vacancy *reposting*
following match separations — so a negative vacancy IRF for LD or QU
contradicts the theory.  If the severity placebo shows that the negative
vacancy response to LD and QU disappears once Δu^nat_t is controlled, this
is the key result: the apparent contradiction was a severity-composition
artifact, not a genuine refutation of the reposting mechanism.  Concretely:
high-layoff industries (Manufacturing, Construction) are also high-cyclicality
industries; states with large shares of those industries receive large LD
Bartik values exactly during aggregate contractions when vacancies are falling
everywhere for macro reasons.

If β_h for LD and QU vacancy becomes statistically insignificant or changes
sign after controlling for severity, the theoretical consistency of the
Beveridge asymmetry is restored: δ shocks suppress vacancies (structural),
LD/QU shocks do not (once severity confound is removed).

Outputs
-------
  data/results/placebo_sev_{delta,ld,qu}_{unemp,vac}.csv
      h, beta_h, se, t, p, delta_h (severity interaction), delta_h_p,
      baseline_beta (part5 baseline for comparison), N

  data/results/placebo_sev_delta_unemp_overlay.png
  data/results/placebo_sev_delta_vac_overlay.png
  data/results/placebo_sev_ld_vac_overlay.png
  data/results/placebo_sev_qu_vac_overlay.png
      Each: baseline IRF (blue) vs. severity-controlled IRF (orange dashed),
      with 90% CI bands.  The key visual test: does orange ≈ blue (pass) or
      collapse to zero (fail)?

  data/results/placebo_sev_interaction_summary.png
      δ_h interaction coefficients for all three instruments × both outcomes,
      with significance markers.  Flat at zero = clean identification.

Prerequisites
-------------
  part3_resid_instruments.py    (residualized Bartik CSVs)
  part4_outcomes.py             (LAUS quarterly parquet with vacancy data)
  part5_lp.py                   (baseline IRF CSVs for comparison)
  part7_nfci.py                 (NFCI cache; optional, secondary severity)

Run
---
  python part7c_recession_placebo.py
"""

import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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

# ── imports from pipeline ─────────────────────────────────────────────────────
from construct_delta_instrument import DEFAULT_CACHE_DIR

# ── FRED fetch (for national unemployment) ────────────────────────────────────
from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY", "")

def _open_file(path):
    try:
        os.startfile(path)
    except AttributeError:
        pass

# =============================================================================
# Configuration
# =============================================================================
HORIZONS            = list(range(21))      # h = 0 … 20 quarters
BASE_YEAR           = 2006
Z90                 = 1.645
Z95                 = 1.960
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DELTA_RESID_FILE = INSTR_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
LD_RESID_FILE    = INSTR_DIR / f"ld_instrument_resid_base{BASE_YEAR}.csv"
QU_RESID_FILE    = INSTR_DIR / f"qu_instrument_resid_base{BASE_YEAR}.csv"
LAUS_FILE        = INSTR_DIR / "laus_quarterly.parquet"

# Severity measure: national unemployment rate (FRED UNRATE)
UNRATE_CACHE = DEFAULT_CACHE_DIR / "unrate_quarterly.parquet"
# Secondary severity: NFCI (optional; used for comparison only)
NFCI_CACHE   = DEFAULT_CACHE_DIR / "nfci_quarterly.parquet"


# =============================================================================
# Quarter arithmetic
# =============================================================================
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


def _ql_to_dt(ql: str) -> pd.Timestamp:
    y, q = ql.split("Q")
    return pd.Timestamp(year=int(y), month=int(q) * 3 - 2, day=1)


# =============================================================================
# Severity series: Δu^nat_t
# =============================================================================
def load_severity() -> pd.DataFrame:
    """
    Load national unemployment rate from FRED (UNRATE), compute both the
    demeaned level and the first difference, and return both as severity
    measures.

    Two measures, two distinct confound hypotheses:

    u^nat_t  (level, demeaned):
        Captures the *regime* — how bad aggregate conditions are at the
        moment of the shock.  This is the primary severity measure.  A
        significant B × u^nat_t interaction means high-Bartik states are
        hit harder when the unemployment level is elevated, regardless of
        whether it is rising or falling.  This directly tests the hypothesis
        that the Bartik instrument is a severity-composition proxy.

    Δu^nat_t (first difference, demeaned):
        Captures *momentum* — how fast conditions are deteriorating.  A
        shock landing when unemployment is rising fast (early recession) is
        different from one landing when unemployment has stabilized at a high
        level (late recession).  Including this alongside the level rules out
        the additional hypothesis that Bartik effects scale with the speed of
        deterioration, not just the level of distress.

    Both are demeaned within the estimation sample (2001Q1–2019Q4) so that
    β_h is interpreted as the IRF at average aggregate conditions — identical
    interpretation convention to part5.

    Returns DataFrame: quarter_label, u_nat_dm (demeaned level),
                                      du_nat_dm (demeaned change).
    """
    if UNRATE_CACHE.exists():
        df = pd.read_parquet(UNRATE_CACHE)
        print("  National unemployment rate: loaded from cache")
    else:
        print("  Fetching UNRATE from FRED ...")
        if not FRED_API_KEY:
            sys.exit(
                "ERROR: FRED_API_KEY not set and no cached UNRATE found.\n"
                "  Set FRED_API_KEY or run part7_nfci.py first."
            )
        from fredapi import Fred
        fred = Fred(api_key=FRED_API_KEY)
        s    = fred.get_series("UNRATE")
        df   = (s.rename("unrate").reset_index()
                 .rename(columns={"index": "date"}))
        df["date"] = pd.to_datetime(df["date"])
        df["quarter_label"] = df["date"].apply(
            lambda d: f"{d.year}Q{(d.month - 1) // 3 + 1}"
        )
        df = df.groupby("quarter_label")["unrate"].mean().reset_index()
        df.to_parquet(UNRATE_CACHE, index=False)
        print(f"  Saved: {UNRATE_CACHE}")

    df = df.sort_values("quarter_label").reset_index(drop=True)
    df["du_nat"] = df["unrate"].diff()
    df = df.dropna(subset=["du_nat"])

    # Demean both in sample
    sample = (df["quarter_label"] >= "2001Q1") & (df["quarter_label"] <= "2019Q4")
    mean_u  = df.loc[sample, "unrate"].mean()
    mean_du = df.loc[sample, "du_nat"].mean()
    df["u_nat_dm"]  = df["unrate"] - mean_u    # primary: level
    df["du_nat_dm"] = df["du_nat"] - mean_du   # secondary: change

    print(f"  u^nat  (level): sample mean removed = {mean_u:.3f}%,  "
          f"SD(in-sample) = {df.loc[sample,'u_nat_dm'].std():.3f} pp")
    print(f"  Δu^nat (change): sample mean removed = {mean_du:.4f} pp,  "
          f"SD(in-sample) = {df.loc[sample,'du_nat_dm'].std():.3f} pp")
    print(f"  Corr(level, change) in sample: "
          f"{df.loc[sample,'u_nat_dm'].corr(df.loc[sample,'du_nat_dm']):.3f} "
          f"(low collinearity → both add distinct information)")

    return df[["quarter_label", "u_nat_dm", "du_nat_dm"]]


def load_nfci_severity() -> pd.DataFrame | None:
    """
    Secondary severity measure: NFCI demeaned in sample.
    Returns None if cache not available.
    """
    if not NFCI_CACHE.exists():
        print("  NFCI cache not found — skipping secondary severity measure.")
        return None
    nfci = pd.read_parquet(NFCI_CACHE)
    nfci.columns = [c.lower() for c in nfci.columns]
    if "quarter_label" not in nfci.columns:
        print("  NFCI cache missing quarter_label — skipping.")
        return None
    # Use headline NFCI (or risk subindex if available)
    col = "nfci_risk" if "nfci_risk" in nfci.columns else \
          "nfci"      if "nfci"      in nfci.columns else None
    if col is None:
        print("  No usable NFCI column found — skipping.")
        return None
    nfci = nfci[["quarter_label", col]].rename(columns={col: "nfci_dm"})
    nfci = nfci.sort_values("quarter_label")
    sample = (nfci["quarter_label"] >= "2001Q1") & (nfci["quarter_label"] <= "2019Q4")
    nfci["nfci_dm"] = nfci["nfci_dm"] - nfci.loc[sample, "nfci_dm"].mean()
    print(f"  NFCI ({col}): demeaned in sample")
    return nfci


# =============================================================================
# Data loading and panel construction
# =============================================================================
def load_resid(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"state_fips": str})
    df["state_fips"] = df["state_fips"].str.zfill(2)
    return df


def load_outcomes() -> pd.DataFrame:
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    keep = ["state_fips", "state", "quarter_label", "unemp_rate", "labor_force"]
    if "vac_rate" in laus.columns:
        keep.append("vac_rate")
    return laus[keep]


def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame,
                severity: pd.DataFrame) -> pd.DataFrame:
    """
    Merge instrument, outcomes, and severity into one panel.

    Creates:
      instr_col           : Bartik shock × 100 (pp)
      u_nat_dm            : u^nat_t demeaned (primary severity: level)
      du_nat_dm           : Δu^nat_t demeaned (secondary severity: change)
      instr_x_u_lev       : B^k_{s,t} × u^nat_t   (level interaction)
      instr_x_u_chg       : B^k_{s,t} × Δu^nat_t  (change interaction)
      unemp_lag1, lf_log_lag1, vac_rate_lag1 : predetermined controls
    """
    merge_cols = ["state_fips", "quarter_label", "unemp_rate", "labor_force"]
    _has_vac   = "vac_rate" in outcomes.columns
    if _has_vac:
        merge_cols.append("vac_rate")

    panel = (instr[["state_fips", "state", "quarter_label", instr_col]]
             .merge(outcomes[merge_cols], on=["state_fips", "quarter_label"],
                    how="inner")
             .merge(severity, on="quarter_label", how="inner"))

    panel[instr_col]       = panel[instr_col] * 100
    panel["instr_x_u_lev"] = panel[instr_col] * panel["u_nat_dm"]
    panel["instr_x_u_chg"] = panel[instr_col] * panel["du_nat_dm"]

    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)
    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )
    if _has_vac:
        panel["vac_rate_lag1"] = panel.groupby("state_fips")["vac_rate"].shift(1)

    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    bad = panel["prev_ql"].notna() & (
        panel["prev_ql"] != quarter_shift(panel["quarter_label"], -1)
    )
    lag_cols = ["unemp_lag1", "lf_log_lag1"] + (["vac_rate_lag1"] if _has_vac else [])
    panel.loc[bad, lag_cols] = np.nan
    panel.drop(columns="prev_ql", inplace=True)

    print(f"  {instr_col}: {panel['state_fips'].nunique()} states × "
          f"{panel['quarter_label'].nunique()} quarters = {len(panel):,} rows")
    return panel


# =============================================================================
# LP estimation — single horizon
# =============================================================================
def run_lp_horizon_sev(panel: pd.DataFrame, h: int,
                       shock_col: str, outcome: str = "unemp") -> dict | None:
    """
    Estimate the severity-placebo LP at horizon h.

    Specification:
      y_{s,t+h} - y_{s,t-1} = α_s + α_t
                               + β_h  B^k_{s,t}
                               + δ_h  B^k_{s,t} × Δu^nat_t
                               + φ_h  Δu^nat_t
                               + γ_1 y_{s,t-1}  +  γ_2 log LF_{s,t-1}
                               + ε_{s,t,h}

    β_h : IRF at average aggregate severity (Δu^nat_t = 0 after demeaning).
          Compare to part5 baseline β_h to assess severity bias.

    δ_h : severity interaction.
          ≈ 0   → identification clean; effect invariant to aggregate severity.
          > 0   → additional amplification during recessions; check β_h.
          β_h → 0 after adding interaction → instrument is a severity proxy.

    φ_h : level effect of national severity on the outcome, independent of
          the instrument.  Absorbs residual aggregate-to-local transmission
          at outcome horizon h (time FE α_t controls the shock quarter t
          but not the outcome quarter t+h).

    Returns dict with β_h, SE, CI, δ_h, φ_h and diagnostics, or None.
    """
    if outcome not in ("unemp", "vacancy"):
        raise ValueError(f"outcome must be 'unemp' or 'vacancy', got {outcome!r}")

    df = panel.copy()
    df["future_ql"] = quarter_shift(df["quarter_label"], h)

    if outcome == "unemp":
        lookup = (df[["state_fips", "quarter_label", "unemp_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "unemp_rate": "y_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["y_future"] - df["unemp_lag1"]
        lag_col = "unemp_lag1"

    else:
        if "vac_rate" not in df.columns or "vac_rate_lag1" not in df.columns:
            return None
        lookup = (df[["state_fips", "quarter_label", "vac_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "vac_rate": "vac_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["vac_future"] - df["vac_rate_lag1"]
        lag_col = "vac_rate_lag1"

    required = ["dep_var", shock_col,
                "instr_x_u_lev", "instr_x_u_chg",
                "u_nat_dm", "du_nat_dm",
                lag_col, "lf_log_lag1"]
    df = df.dropna(subset=required).copy()
    if len(df) < 100:
        return None

    state_dummies = pd.get_dummies(df["state_fips"],    prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)

    core = [shock_col,
            "instr_x_u_lev", "instr_x_u_chg",
            "u_nat_dm", "du_nat_dm",
            lag_col, "lf_log_lag1"]
    X = sm.add_constant(
        pd.concat([df[core], state_dummies, time_dummies], axis=1),
        has_constant="add"
    )
    model = sm.OLS(df["dep_var"], X).fit(
        cov_type="cluster",
        cov_kwds={"groups": df["state_fips"].values}
    )

    beta  = float(model.params[shock_col])
    se    = float(model.bse[shock_col])
    tstat = float(model.tvalues[shock_col])
    pval  = float(model.pvalues[shock_col])

    return {
        "h": h, "beta": beta, "se": se, "tstat": tstat, "pval": pval,
        "partial_f": tstat**2,
        "ci90_lo": beta - Z90*se, "ci90_hi": beta + Z90*se,
        "ci95_lo": beta - Z95*se, "ci95_hi": beta + Z95*se,
        # Level interaction δ^lev_h (primary)
        "delta_lev_h":  float(model.params["instr_x_u_lev"]),
        "delta_lev_se": float(model.bse["instr_x_u_lev"]),
        "delta_lev_t":  float(model.tvalues["instr_x_u_lev"]),
        "delta_lev_p":  float(model.pvalues["instr_x_u_lev"]),
        # Change interaction δ^chg_h (secondary)
        "delta_chg_h":  float(model.params["instr_x_u_chg"]),
        "delta_chg_se": float(model.bse["instr_x_u_chg"]),
        "delta_chg_t":  float(model.tvalues["instr_x_u_chg"]),
        "delta_chg_p":  float(model.pvalues["instr_x_u_chg"]),
        # Severity level effects
        "phi_lev_h": float(model.params["u_nat_dm"]),
        "phi_lev_p": float(model.pvalues["u_nat_dm"]),
        "phi_chg_h": float(model.params["du_nat_dm"]),
        "phi_chg_p": float(model.pvalues["du_nat_dm"]),
        "nobs":       int(model.nobs),
        "r2":         float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range":   f"{df['quarter_label'].min()}–{df['quarter_label'].max()}",
        "outcome":    outcome,
    }


# =============================================================================
# LP loop — all horizons, with scaling
# =============================================================================
def run_lp_sev(panel: pd.DataFrame, shock_col: str,
               label: str, outcome: str = "unemp") -> pd.DataFrame:
    """
    Run severity-placebo LP for all horizons.

    Prints both interaction coefficients (level and change) per horizon.
    Scales β_h, δ^lev_h, δ^chg_h to 1-SD instrument units (part5 convention).
    """
    tag = "→ vacancy rate" if outcome == "vacancy" else "→ unemp rate"
    print(f"\n  LP (severity placebo): {label}  {tag}  h=0..{max(HORIZONS)}")
    print(f"  {'h':>3}  {'beta_h':>9}  {'SE':>7}  {'p':>6}  "
          f"  {'dlev_h':>8}  {'lev_p':>6}  {'dchg_h':>8}  {'chg_p':>6}  {'N':>6}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon_sev(panel, h, shock_col, outcome)
        if res is None:
            continue
        rows.append(res)
        sig  = ("***" if res["pval"]        < .01 else "**" if res["pval"]        < .05
                else "*" if res["pval"]        < .10 else "")
        lsig = ("***" if res["delta_lev_p"] < .01 else "**" if res["delta_lev_p"] < .05
                else "*" if res["delta_lev_p"] < .10 else "")
        csig = ("***" if res["delta_chg_p"] < .01 else "**" if res["delta_chg_p"] < .05
                else "*" if res["delta_chg_p"] < .10 else "")
        print(f"  {h:3d}  {res['beta']:9.4f}  {res['se']:7.4f}  "
              f"{res['pval']:6.3f}{sig:<3}  "
              f"  {res['delta_lev_h']:8.4f}  {res['delta_lev_p']:6.3f}{lsig:<3}"
              f"  {res['delta_chg_h']:8.4f}  {res['delta_chg_p']:6.3f}{csig:<3}"
              f"  {res['nobs']:6d}")

    if not rows:
        return pd.DataFrame()

    irf = pd.DataFrame(rows)
    sd  = float(panel[shock_col].std())
    for col in ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]:
        irf[col] = irf[col] * sd
    # Scale both interaction coefficients
    for col in ["delta_lev_h", "delta_lev_se", "delta_chg_h", "delta_chg_se"]:
        irf[col] = irf[col] * sd
    irf["instr_sd"] = sd
    irf["label"]    = label
    irf["outcome"]  = outcome
    return irf


# =============================================================================
# Plotting
# =============================================================================
_PALETTE = {"delta": "#1f77b4", "ld": "#d62728", "qu": "#2ca02c"}


def _plot_overlay(irf_base: pd.DataFrame, irf_sev: pd.DataFrame,
                  label: str, shock: str, outcome: str,
                  out_path: Path) -> None:
    """
    Overlay plot: baseline IRF (blue solid) vs. severity-controlled IRF
    (orange dashed).

    The key visual test:
      - Lines overlap / close → severity is NOT driving the result (pass)
      - Orange collapses toward zero → severity WAS driving the result (fail)
      - Both remain significant, orange lower → genuine channel + amplification
    """
    plt.ioff()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    col  = _PALETTE.get(shock, "#555555")
    ytag = "pp change in vacancy rate" if outcome == "vacancy" \
           else "pp change in unemp. rate"

    # Left: point estimates
    ax = axes[0]
    ax.axhline(0, color="black", lw=0.8)
    for hh in range(4, 21, 4):
        ax.axvline(hh, color="grey", lw=0.4, ls=":", alpha=0.5)

    if not irf_base.empty:
        h, b = irf_base["h"].values, irf_base["beta"].values
        _ = ax.fill_between(h, irf_base["ci90_lo"], irf_base["ci90_hi"],
                            color=col, alpha=0.15)
        _ = ax.plot(h, b, color=col, lw=2.0, ls="-",
                    marker="o", ms=3.5, label="Baseline (part5)")

    if not irf_sev.empty:
        h, b = irf_sev["h"].values, irf_sev["beta"].values
        _ = ax.fill_between(h, irf_sev["ci90_lo"], irf_sev["ci90_hi"],
                            color="#ff7f0e", alpha=0.15)
        _ = ax.plot(h, b, color="#ff7f0e", lw=2.0, ls="--",
                    marker="o", ms=3.5,
                    label=r"+ $u^{nat}_t$ & $\Delta u^{nat}_t$ controls")

    _ = ax.set_title(f"{label} → {outcome}\nBaseline vs. severity-controlled",
                     fontsize=11)
    _ = ax.set_xlabel("Horizon h (quarters)", fontsize=10)
    _ = ax.set_ylabel(f"{ytag}\nper 1-SD shock", fontsize=10)
    _ = ax.set_xticks(range(0, 21, 2))
    _ = ax.legend(fontsize=9, framealpha=0.85)
    _ = ax.grid(axis="y", lw=0.4, alpha=0.4)

    # Right: both severity interaction coefficients δ^lev_h and δ^chg_h
    ax2 = axes[1]
    ax2.axhline(0, color="black", lw=0.8)
    for hh in range(4, 21, 4):
        ax2.axvline(hh, color="grey", lw=0.4, ls=":", alpha=0.5)

    if not irf_sev.empty:
        h = irf_sev["h"].values
        # Level interaction (primary) — purple solid
        dl  = irf_sev["delta_lev_h"].values
        dls = irf_sev["delta_lev_se"].values
        _ = ax2.fill_between(h, dl - Z90*dls, dl + Z90*dls,
                             color="#9467bd", alpha=0.13)
        _ = ax2.plot(h, dl, color="#9467bd", lw=2.0, ls="-",
                     marker="o", ms=3.5,
                     label=r"$\delta^{lev}_h$ (B × $u^{nat}_t$, primary)")
        # Change interaction (secondary) — brown dashed
        dc  = irf_sev["delta_chg_h"].values
        dcs = irf_sev["delta_chg_se"].values
        _ = ax2.fill_between(h, dc - Z90*dcs, dc + Z90*dcs,
                             color="#8c564b", alpha=0.10)
        _ = ax2.plot(h, dc, color="#8c564b", lw=1.8, ls="--",
                     marker="s", ms=3.0,
                     label=r"$\delta^{chg}_h$ (B × $\Delta u^{nat}_t$, secondary)")
        # Significance markers
        for _, row in irf_sev.iterrows():
            for val, pval, color, offset in [
                (row["delta_lev_h"], row["delta_lev_p"], "#9467bd",  4),
                (row["delta_chg_h"], row["delta_chg_p"], "#8c564b", -6),
            ]:
                if pval < 0.10:
                    mk = "***" if pval < 0.01 else "**" if pval < 0.05 else "*"
                    _ = ax2.annotate(mk, xy=(row["h"], val),
                                     fontsize=6, ha="center", va="bottom",
                                     color=color,
                                     xytext=(0, offset), textcoords="offset points")

    _ = ax2.set_title(
        r"Severity interactions: level $u^{nat}_t$ (solid) & change $\Delta u^{nat}_t$ (dashed)"
        "\nBoth ≈ 0 → identification robust; β_h not severity-driven",
        fontsize=10
    )
    _ = ax2.set_xlabel("Horizon h (quarters)", fontsize=10)
    _ = ax2.set_ylabel(r"$\delta_h$ per 1-SD shock per pp", fontsize=10)
    _ = ax2.set_xticks(range(0, 21, 2))
    _ = ax2.legend(fontsize=8, framealpha=0.85)
    _ = ax2.grid(axis="y", lw=0.4, alpha=0.4)

    fig.suptitle(
        f"Recession-severity placebo — {label} → {outcome}\n"
        r"Orange = β_h after controlling for $u^{nat}_t$ (level) & $\Delta u^{nat}_t$ (change). "
        "Overlap with baseline → severity not driving result.",
        fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


def _plot_interaction_summary(results: dict, out_path: Path) -> None:
    """
    Summary panel: both δ^lev_h (level) and δ^chg_h (change) interaction
    coefficients for all instrument × outcome combinations.
    Solid = level interaction (primary); dashed = change interaction (secondary).
    Both flat near zero = identification clean on both dimensions.
    """
    plt.ioff()
    items = [(k, df) for k, df in results.items() if not df.empty]
    n     = len(items)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5.5 * ncols, 4 * nrows),
                             squeeze=False)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    for idx, (key, irf) in enumerate(items):
        ax  = axes[idx // ncols][idx % ncols]
        h   = irf["h"].values
        col = colors[idx % len(colors)]

        ax.axhline(0, color="black", lw=0.8)
        for hh in range(4, 21, 4):
            ax.axvline(hh, color="grey", lw=0.4, ls=":", alpha=0.5)

        # Level interaction (primary) — solid
        dl  = irf["delta_lev_h"].values
        dls = irf["delta_lev_se"].values
        _ = ax.fill_between(h, dl - Z90*dls, dl + Z90*dls, color=col, alpha=0.13)
        _ = ax.plot(h, dl, color=col, lw=1.8, ls="-", marker="o", ms=3.0,
                    label=r"$u^{nat}$ (level)")

        # Change interaction (secondary) — dashed, lighter
        dc  = irf["delta_chg_h"].values
        dcs = irf["delta_chg_se"].values
        _ = ax.plot(h, dc, color=col, lw=1.4, ls="--", marker="s", ms=2.5,
                    alpha=0.6, label=r"$\Delta u^{nat}$ (change)")

        _ = ax.set_title(key, fontsize=9)
        _ = ax.set_xlabel("h", fontsize=8)
        _ = ax.set_ylabel(r"$\delta_h$", fontsize=8)
        _ = ax.set_xticks(range(0, 21, 4))
        _ = ax.legend(fontsize=7, framealpha=0.8)
        _ = ax.grid(axis="y", lw=0.4, alpha=0.4)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(
        r"Severity interactions summary — solid: $u^{nat}_t$ (level), "
        r"dashed: $\Delta u^{nat}_t$ (change)" + "\n"
        "Both flat near zero = identification robust to recession severity",
        fontsize=11, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


# =============================================================================
# Load part5 baseline for comparison
# =============================================================================
def _load_baseline(fname: str) -> pd.DataFrame:
    p = RESULTS_DIR / fname
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p)


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 65)
    print("Part 7c — Recession-severity placebo")
    print("=" * 65)

    # ── 1. Check required inputs ──────────────────────────────────────────────
    missing = [p for p in [DELTA_RESID_FILE, LD_RESID_FILE,
                            QU_RESID_FILE, LAUS_FILE] if not p.exists()]
    if missing:
        sys.exit(
            "Missing required files:\n" +
            "\n".join(f"  {p}" for p in missing) +
            "\nRun part3_resid_instruments.py and part4_outcomes.py first."
        )

    # ── 2. Load data ──────────────────────────────────────────────────────────
    print("\n[1] Loading data")
    outcomes = load_outcomes()
    vac_ok   = "vac_rate" in outcomes.columns and outcomes["vac_rate"].notna().any()
    print(f"  Vacancy data available: {vac_ok}")

    print("\n[2] Severity measure: Δu^nat_t (national unemployment change)")
    severity = load_severity()

    # ── 3. Build panels ───────────────────────────────────────────────────────
    print("\n[3] Building panels")
    delta_panel = build_panel(load_resid(DELTA_RESID_FILE),
                              "bartik_delta", outcomes, severity)
    ld_panel    = build_panel(load_resid(LD_RESID_FILE),
                              "bartik_ld",    outcomes, severity)
    qu_panel    = build_panel(load_resid(QU_RESID_FILE),
                              "bartik_qu",    outcomes, severity)

    # Restrict δ to 2001Q1+ for common window (matches LD/QU JOLTS start)
    delta_panel = delta_panel[delta_panel["quarter_label"] >= "2001Q1"].copy()

    # ── 4. Severity-placebo LPs ───────────────────────────────────────────────
    print("\n[4] Severity-placebo LPs — unemployment outcome")
    irf_d_u   = run_lp_sev(delta_panel, "bartik_delta", "δ  (sev+unemp)",  "unemp")
    irf_ld_u  = run_lp_sev(ld_panel,    "bartik_ld",    "LD (sev+unemp)",  "unemp")
    irf_qu_u  = run_lp_sev(qu_panel,    "bartik_qu",    "QU (sev+unemp)",  "unemp")

    if vac_ok:
        print("\n[5] Severity-placebo LPs — vacancy outcome")
        irf_d_v   = run_lp_sev(delta_panel, "bartik_delta", "δ  (sev+vac)",  "vacancy")
        irf_ld_v  = run_lp_sev(ld_panel,    "bartik_ld",    "LD (sev+vac)",  "vacancy")
        irf_qu_v  = run_lp_sev(qu_panel,    "bartik_qu",    "QU (sev+vac)",  "vacancy")
    else:
        irf_d_v = irf_ld_v = irf_qu_v = pd.DataFrame()

    # ── 5. Load part5 baselines for comparison ────────────────────────────────
    print("\n[6] Loading part5 baselines for comparison")
    base_d_u  = _load_baseline("lp_irf_delta_resid.csv")
    base_ld_u = _load_baseline("lp_irf_ld_resid.csv")
    base_qu_u = _load_baseline("lp_irf_qu_resid.csv")
    base_d_v  = _load_baseline("lp_irf_delta_vacancy.csv")
    base_ld_v = _load_baseline("lp_irf_ld_vacancy.csv")
    base_qu_v = _load_baseline("lp_irf_qu_vacancy.csv")

    for lbl, df in [("delta unemp", base_d_u), ("ld unemp", base_ld_u),
                    ("delta vac",   base_d_v)]:
        print(f"  {lbl}: {'loaded' if not df.empty else 'not found (run part5 first)'}")

    # ── 6. Save CSVs ──────────────────────────────────────────────────────────
    print("\n[7] Saving CSVs")
    def _save(irf, fname):
        if irf is not None and not irf.empty:
            p = RESULTS_DIR / fname
            irf.to_csv(p, index=False)
            print(f"  Saved: {p}")

    _save(irf_d_u,  "placebo_sev_delta_unemp.csv")
    _save(irf_ld_u, "placebo_sev_ld_unemp.csv")
    _save(irf_qu_u, "placebo_sev_qu_unemp.csv")
    _save(irf_d_v,  "placebo_sev_delta_vac.csv")
    _save(irf_ld_v, "placebo_sev_ld_vac.csv")
    _save(irf_qu_v, "placebo_sev_qu_vac.csv")

    # ── 7. Plots ──────────────────────────────────────────────────────────────
    print("\n[8] Plotting")

    _plot_overlay(base_d_u,  irf_d_u,  "δ",  "delta", "unemp",
                  RESULTS_DIR / "placebo_sev_delta_unemp_overlay.png")
    _plot_overlay(base_ld_u, irf_ld_u, "LD", "ld",    "unemp",
                  RESULTS_DIR / "placebo_sev_ld_unemp_overlay.png")
    _plot_overlay(base_qu_u, irf_qu_u, "QU", "qu",    "unemp",
                  RESULTS_DIR / "placebo_sev_qu_unemp_overlay.png")

    if vac_ok:
        _plot_overlay(base_d_v,  irf_d_v,  "δ",  "delta", "vacancy",
                      RESULTS_DIR / "placebo_sev_delta_vac_overlay.png")
        _plot_overlay(base_ld_v, irf_ld_v, "LD", "ld",    "vacancy",
                      RESULTS_DIR / "placebo_sev_ld_vac_overlay.png")
        _plot_overlay(base_qu_v, irf_qu_v, "QU", "qu",    "vacancy",
                      RESULTS_DIR / "placebo_sev_qu_vac_overlay.png")

    # Summary: all interaction δ_h coefficients
    interaction_results = {
        "δ → unemp":    irf_d_u,
        "LD → unemp":   irf_ld_u,
        "QU → unemp":   irf_qu_u,
    }
    if vac_ok:
        interaction_results.update({
            "δ → vacancy":  irf_d_v,
            "LD → vacancy": irf_ld_v,
            "QU → vacancy": irf_qu_v,
        })
    _plot_interaction_summary(
        interaction_results,
        RESULTS_DIR / "placebo_sev_interaction_summary.png"
    )

    # ── 8. Summary table ──────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("RECESSION-SEVERITY PLACEBO — SUMMARY")
    print("=" * 65)
    print(
        "\nKey question: does β_h survive severity control?\n"
        "  β_h ≈ baseline → severity NOT driving result (identification robust)\n"
        "  β_h → 0        → severity WAS driving result (instrument failure)\n"
        "  δ_h ≈ 0        → no amplification effect (cleanest result)\n"
        "  δ_h > 0 + β_h survives → genuine channel with recession amplification\n"
        "  Severity measures: u^nat_t = level (primary); Δu^nat_t = change (secondary)\n"
    )
    print(f"  {'Series':<20}  {'outcome':>8}  {'h_peak':>6}  "
          f"{'beta_peak':>9}  {'p':>6}  "
          f"{'dlev_peak':>10}  {'lev_p':>6}  "
          f"{'dchg_peak':>10}  {'chg_p':>6}")
    print(f"  {'-'*90}")
    for lbl, irf in [("δ (sev)",  irf_d_u),  ("LD (sev)", irf_ld_u),
                     ("QU (sev)", irf_qu_u), ("δ (sev)",  irf_d_v),
                     ("LD (sev)", irf_ld_v), ("QU (sev)", irf_qu_v)]:
        if irf is None or irf.empty:
            continue
        pk   = irf.loc[irf["beta"].abs().idxmax()]
        lpk  = irf.loc[irf["delta_lev_h"].abs().idxmax()]
        cpk  = irf.loc[irf["delta_chg_h"].abs().idxmax()]
        print(f"  {lbl:<20}  {pk['outcome']:>8}  {int(pk['h']):>6}  "
              f"{pk['beta']:>9.4f}  {pk['pval']:>6.3f}  "
              f"{lpk['delta_lev_h']:>10.4f}  {lpk['delta_lev_p']:>6.3f}  "
              f"{cpk['delta_chg_h']:>10.4f}  {cpk['delta_chg_p']:>6.3f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
