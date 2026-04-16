"""
part7b_sloos.py  --  SLOOS C&I Tightening as Alternative Interaction Variable
==============================================================================
Replicates the NFCI × δ interaction LP from part5_lp.py, replacing NFCI_risk
with the Senior Loan Officer Opinion Survey (SLOOS) C&I net tightening index.

Research question
-----------------
The NFCI interaction (part5) shows the δ unemployment effect concentrates
almost entirely in tight financial conditions (δ_h significant, β_h ≈ 0 at
average conditions).  Two competing explanations:

  (A) Credit supply channel: firm exit is amplified by restricted credit
      access.  SLOOS directly measures banks' *willingness* to lend (supply),
      so the δ × SLOOS interaction should survive if (A) is correct.

  (B) Recession-severity proxy: NFCI is correlated with recession depth.
      If NFCI is just picking up GDP contractions, the interaction should
      weaken once we use a pure credit-supply measure.

Comparing δ_h(SLOOS) vs δ_h(NFCI) directly tests (A) vs (B).

SLOOS series used
-----------------
  DRSDCILM  Net % banks tightening C&I standards, large/medium firms  (primary)
  DRSDCIS   Net % banks tightening C&I standards, small firms          (secondary)
  Both: quarterly, NSA, back to 1990Q4.  Source: Federal Reserve / FRED.

Interaction specification
-------------------------
  y_{s,t+h} - y_{s,t-1} = a_s + a_t
                           + beta_h  * B_{s,t}^delta
                           + delta_h * B_{s,t}^delta x SLOOS_dm_t
                           + gamma_1 * y_{s,t-1} + gamma_2 * log(LF_{s,t-1})
                           + eps

  SLOOS_dm_t demeaned within estimation sample so beta_h = IRF at average
  credit conditions.

Outputs (data/results/)
-----------------------
  sloos_quarterly.parquet             -- SLOOS quarterly series (cached)
  lp_irf_delta_sloos_lm.csv          -- δ × SLOOS_LM interaction IRF
  lp_irf_delta_sloos_sm.csv          -- δ × SLOOS_SM interaction IRF (robustness)
  lp_irf_delta_sloos_comparison.png  -- β_h and δ_h for SLOOS_LM vs NFCI

Prerequisites
-------------
  part3_resid_instruments.py  (produces delta_instrument_resid_base2006.csv)
  part4_outcomes.py           (produces laus_quarterly.parquet)
  part5_lp.py                 (produces lp_irf_delta_nfci.csv for comparison)
  FRED_API_KEY env var        (same key used in part2b)
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

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    pass

# ---------------------------------------------------------------------------
# Configuration  (mirror part5_lp.py exactly)
# ---------------------------------------------------------------------------
HORIZONS          = list(range(21))   # h = 0 .. 20
BASE_YEAR         = 2006
CI_LEVEL          = 0.90
Z90               = 1.645
Z95               = 1.960
MAX_OUTCOME_QUARTER = "2019Q4"        # COVID exclusion

INSTR_DIR   = Path("data/instruments")
CACHE_DIR   = Path("data/cache")
RESULTS_DIR = Path("data/results")
for d in [CACHE_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

DELTA_RESID_FILE = INSTR_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
LAUS_FILE        = INSTR_DIR / "laus_quarterly.parquet"
NFCI_FILE        = CACHE_DIR / "nfci_quarterly.parquet"
SLOOS_CACHE      = CACHE_DIR / "sloos_quarterly.parquet"

# FRED series IDs for SLOOS C&I net tightening
SLOOS_SERIES = {
    "sloos_lm": "DRSDCILM",   # large & middle-market firms  (primary)
    "sloos_sm": "DRSDCIS",    # small firms                  (robustness)
}

# ---------------------------------------------------------------------------
# FRED API key
# ---------------------------------------------------------------------------
from dotenv import load_dotenv
load_dotenv()
FRED_API_KEY = os.getenv("FRED_API_KEY")


# ---------------------------------------------------------------------------
# Quarter arithmetic  (identical to part5_lp.py)
# ---------------------------------------------------------------------------
def quarter_shift(ql_series: pd.Series, h: int) -> pd.Series:
    result = []
    for ql in ql_series:
        p = pd.Period(ql, freq="Q") + h
        result.append(f"{p.year}Q{p.quarter}")
    return pd.Series(result, index=ql_series.index)


def _dt_to_ql(dt: pd.Timestamp) -> str:
    return f"{dt.year}Q{(dt.month - 1) // 3 + 1}"


# ---------------------------------------------------------------------------
# SLOOS fetch & cache
# ---------------------------------------------------------------------------
def _get_fred_client():
    if not FRED_API_KEY:
        raise EnvironmentError(
            "FRED_API_KEY is not set.\n"
            "  Register free at https://fred.stlouisfed.org/docs/api/api_key.html\n"
            "  Then: export FRED_API_KEY=<your_key>  (or add to .env)"
        )
    from fredapi import Fred
    return Fred(api_key=FRED_API_KEY)


def fetch_sloos() -> pd.DataFrame:
    """
    Fetch SLOOS C&I net-tightening series from FRED and return a quarterly
    DataFrame with columns [quarter_label, sloos_lm, sloos_sm].

    SLOOS is already quarterly (survey conducted ~4 weeks after quarter-end).
    FRED returns it with a timestamp at the start of the reference quarter.
    """
    if SLOOS_CACHE.exists():
        print(f"  SLOOS: loading from cache ({SLOOS_CACHE})")
        return pd.read_parquet(SLOOS_CACHE)

    print("  SLOOS: fetching from FRED ...")
    fred = _get_fred_client()

    frames = {}
    for col, sid in SLOOS_SERIES.items():
        try:
            s = fred.get_series(sid, observation_start="1990-01-01")
            # Resample to quarter-start to align with other series
            s = s.resample("QS").first().dropna()
            frames[col] = s
            print(f"    {sid}: {s.index.min().date()} – {s.index.max().date()}"
                  f"  ({len(s)} quarters)")
        except Exception as e:
            raise ValueError(
                f"Could not fetch FRED series '{sid}': {e}\n"
                f"  Verify at https://fred.stlouisfed.org/series/{sid}"
            ) from e

    df = pd.DataFrame(frames).reset_index().rename(columns={"index": "date"})
    df["quarter_label"] = df["date"].apply(_dt_to_ql)
    df = df[["quarter_label"] + list(SLOOS_SERIES.keys())].copy()

    df.to_parquet(SLOOS_CACHE, index=False)
    print(f"  SLOOS: saved to {SLOOS_CACHE}")
    return df


# ---------------------------------------------------------------------------
# Data loading  (mirror part5_lp.py)
# ---------------------------------------------------------------------------
def load_delta_resid() -> pd.DataFrame:
    if not DELTA_RESID_FILE.exists():
        sys.exit(
            f"[ERROR] {DELTA_RESID_FILE} not found.\n"
            "  Run part2b_shock_comovement.py then part3_resid_instruments.py first."
        )
    df = pd.read_csv(DELTA_RESID_FILE, dtype={"state_fips": str})
    df["state_fips"] = df["state_fips"].str.zfill(2)
    print(f"  δ resid instrument: {df['quarter_label'].min()} – "
          f"{df['quarter_label'].max()}  ({df['state_fips'].nunique()} states)")
    return df


def load_outcomes() -> pd.DataFrame:
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    keep = ["state_fips", "state", "quarter_label", "unemp_rate", "labor_force"]
    if "vacancies" in laus.columns:
        keep.append("vacancies")
    return laus[keep]


def load_nfci() -> pd.DataFrame:
    """Load existing NFCI results for comparison plot."""
    if not NFCI_FILE.exists():
        print(f"  [warn] NFCI cache not found at {NFCI_FILE} -- skipping comparison")
        return pd.DataFrame()
    return pd.read_parquet(NFCI_FILE)[["quarter_label", "nfci_risk"]]


# ---------------------------------------------------------------------------
# Panel construction  (mirror part5_lp.py)
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame) -> pd.DataFrame:
    merge_cols = ["state_fips", "quarter_label", "unemp_rate", "labor_force"]
    _has_vac   = "vacancies" in outcomes.columns
    if _has_vac:
        merge_cols.append("vacancies")

    panel = instr[["state_fips", "state", "quarter_label", instr_col]].merge(
        outcomes[merge_cols], on=["state_fips", "quarter_label"], how="inner"
    )
    panel[instr_col] = panel[instr_col] * 100   # decimal -> pp
    panel = panel.sort_values(["state_fips", "quarter_label"]).reset_index(drop=True)

    panel["unemp_lag1"]  = panel.groupby("state_fips")["unemp_rate"].shift(1)
    panel["lf_log_lag1"] = panel.groupby("state_fips")["labor_force"].transform(
        lambda x: np.log(x.shift(1))
    )
    if _has_vac:
        panel["vac_rate"]      = panel["vacancies"] * 1000 / panel["labor_force"] * 100
        panel["vac_rate_lag1"] = panel.groupby("state_fips")["vac_rate"].shift(1)

    # Null out lags across sequence gaps
    panel["prev_ql"] = panel.groupby("state_fips")["quarter_label"].shift(1)
    bad = panel["prev_ql"].notna() & (
        panel["prev_ql"] != quarter_shift(panel["quarter_label"], -1)
    )
    lag_cols = ["unemp_lag1", "lf_log_lag1"] + (["vac_rate_lag1"] if _has_vac else [])
    panel.loc[bad, lag_cols] = np.nan
    panel.drop(columns="prev_ql", inplace=True)
    return panel


def attach_interaction(panel: pd.DataFrame, interact_df: pd.DataFrame,
                       interact_col: str, instr_col: str) -> pd.DataFrame:
    """
    Merge a national interaction variable and add demeaned interaction term.
    interact_df must have columns [quarter_label, <interact_col>].
    Demeaning over the panel's estimation sample ensures beta_h = IRF at
    average credit conditions (same convention as NFCI in part5_lp.py).
    """
    panel = panel.merge(interact_df[["quarter_label", interact_col]],
                        on="quarter_label", how="left")
    n_miss = panel[interact_col].isna().sum()
    if n_miss:
        print(f"  [warn] {n_miss} rows missing {interact_col} -- "
              "check date coverage")
    dm_col = f"{interact_col}_dm"
    panel[dm_col] = panel[interact_col] - panel[interact_col].mean()
    panel[f"instr_x_{interact_col}"] = panel[instr_col] * panel[dm_col]
    return panel


# ---------------------------------------------------------------------------
# LP estimation  (parameterized interaction_col; mirrors part5_lp.py)
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel: pd.DataFrame, h: int,
                   shock_col: str, interaction_col: str,
                   outcome: str = "unemp") -> dict | None:
    """
    Estimate LP at horizon h with a generic interaction term.
    interaction_col: column name of the pre-built B × X_dm term.
    """
    df = base_panel.copy()
    df["future_ql"] = quarter_shift(df["quarter_label"], h)

    if outcome == "unemp":
        lookup = (df[["state_fips", "quarter_label", "unemp_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "unemp_rate": "y_future"}))
        df        = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df    = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["y_future"] - df["unemp_lag1"]
        lag_control   = "unemp_lag1"

    elif outcome == "vacancy":
        if "vac_rate" not in df.columns:
            return None
        lookup = (df[["state_fips", "quarter_label", "vac_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "vac_rate": "vac_future"}))
        df        = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df    = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["vac_future"] - df["vac_rate_lag1"]
        lag_control   = "vac_rate_lag1"
    else:
        raise ValueError(f"outcome must be 'unemp' or 'vacancy', got {outcome!r}")

    required = ["dep_var", shock_col, interaction_col, lag_control, "lf_log_lag1"]
    df = df.dropna(subset=required).copy()
    if len(df) < 100:
        return None

    state_dummies = pd.get_dummies(df["state_fips"],    prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)
    core = [shock_col, interaction_col, lag_control, "lf_log_lag1"]
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
        "delta_h":    float(model.params[interaction_col]),
        "delta_h_se": float(model.bse[interaction_col]),
        "delta_h_t":  float(model.tvalues[interaction_col]),
        "delta_h_p":  float(model.pvalues[interaction_col]),
        "nobs": int(model.nobs), "r2": float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range": f"{df['quarter_label'].min()}-{df['quarter_label'].max()}",
        "outcome": outcome,
    }


def run_lp(panel: pd.DataFrame, shock_col: str, interaction_col: str,
           label: str, outcome: str = "unemp") -> pd.DataFrame:
    """Run LP for all horizons, scale to 1-SD units, return DataFrame."""
    outcome_tag = "-> vacancy" if outcome == "vacancy" else "-> unemp"
    print(f"\n  LP: {label}  {outcome_tag}  h=0..{max(HORIZONS)}")
    print(f"  {'h':>3}  {'beta':>9}  {'SE':>7}  {'t':>7}  {'p':>6}  "
          f"{'delta_h':>9}  {'dh_p':>6}  {'N':>6}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(panel, h, shock_col, interaction_col, outcome)
        if res is None:
            continue
        rows.append(res)
        sig = ("***" if res["pval"] < .01 else "**" if res["pval"] < .05
               else "*" if res["pval"] < .10 else "")
        print(f"  {h:3d}  {res['beta']:9.4f}  {res['se']:7.4f}  "
              f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
              f"{res['delta_h']:9.4f}  {res['delta_h_p']:6.3f}  "
              f"{res['nobs']:6d}  {sig}")

    irf = pd.DataFrame(rows)
    if irf.empty:
        return irf
    sd = float(panel[shock_col].std())
    for col in ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]:
        irf[col] = irf[col] * sd
    irf["delta_h"]    = irf["delta_h"]    * sd
    irf["delta_h_se"] = irf["delta_h_se"] * sd
    irf["instr_sd"]   = sd
    irf["outcome"]    = outcome
    return irf


# ---------------------------------------------------------------------------
# Comparison plot: SLOOS_LM vs NFCI interaction coefficients
# ---------------------------------------------------------------------------
def plot_sloos_vs_nfci(irf_sloos: pd.DataFrame, irf_nfci: pd.DataFrame,
                       out_path: Path,
                       outcome: str = "unemp") -> None:
    """
    Two-panel figure.  Left: β_h (main effect at average conditions).
    Right: δ_h (interaction coefficient).
    Both panels show SLOOS_LM (primary) and NFCI (comparison) on same axes.
    outcome: 'unemp' or 'vacancy' -- controls ylabel and suptitle only.
    """
    if irf_nfci.empty:
        print("  [warn] NFCI IRF not available -- plotting SLOOS only")
        _plot_single(irf_sloos, out_path, label="SLOOS C&I (large/med)")
        return

    ylabel  = ("pp change in vacancy rate\nper 1-SD shock"
               if outcome == "vacancy"
               else "pp change in unemp. rate\nper 1-SD shock")
    suptitle = (r"$\delta$ shock: SLOOS vs NFCI interaction — "
                + ("vacancies" if outcome == "vacancy" else "unemployment"))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = {"SLOOS C&I (LM)": "#d62728", "NFCI risk": "#1f77b4"}

    for ax, col, title in zip(
        axes,
        ["beta", "delta_h"],
        [r"$\beta_h$ — main effect (avg. credit conditions)",
         r"$\delta_h$ — interaction coefficient"],
    ):
        for (label, irf, se_col) in [
            ("SLOOS C&I (LM)", irf_sloos,
             "se"        if col == "beta" else "delta_h_se"),
            ("NFCI risk",      irf_nfci,
             "se"        if col == "beta" else "delta_h_se"),
        ]:
            if irf is None or irf.empty or col not in irf.columns:
                continue
            c   = colors[label]
            h   = irf["h"].values
            b   = irf[col].values
            se  = irf[se_col].values
            ax.fill_between(h, b - Z90*se, b + Z90*se, color=c, alpha=0.15)
            ax.plot(h, b, color=c, lw=2.0, marker="o", ms=3.5, label=label)

        ax.axhline(0, color="black", lw=0.8)
        for hh in range(4, max(HORIZONS)+1, 4):
            ax.axvline(hh, color="grey", lw=0.5, ls=":", alpha=0.6)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.legend(fontsize=9)

    fig.suptitle(suptitle, fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def _plot_single(irf: pd.DataFrame, out_path: Path, label: str) -> None:
    """Fallback: plot SLOOS results alone when NFCI comparison is unavailable."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, col, title in zip(
        axes,
        ["beta", "delta_h"],
        [r"$\beta_h$ — main effect", r"$\delta_h$ — interaction"],
    ):
        se_col = "se" if col == "beta" else "delta_h_se"
        h, b, se = irf["h"].values, irf[col].values, irf[se_col].values
        ax.fill_between(h, b - Z90*se, b + Z90*se, alpha=0.20)
        ax.plot(h, b, lw=2.0, marker="o", ms=3.5, label=label)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel("pp per 1-SD shock", fontsize=10)
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 68)
    print("part7b_sloos.py  --  SLOOS C&I interaction LP")
    print("=" * 68)

    # [1] Fetch / load SLOOS
    print("\n[1] SLOOS data")
    sloos = fetch_sloos()
    print(f"  Coverage: {sloos['quarter_label'].min()} – "
          f"{sloos['quarter_label'].max()}")
    for col in SLOOS_SERIES:
        dm = sloos[col] - sloos[col].mean()
        print(f"  {col}: mean={sloos[col].mean():.2f}  "
              f"SD={sloos[col].std():.2f}  "
              f"[demeaned SD={dm.std():.2f}]")

    # [2] Load instruments and outcomes
    print("\n[2] Loading instruments and outcomes")
    delta_resid = load_delta_resid()
    outcomes    = load_outcomes()
    nfci_df     = load_nfci()

    # [3] Build base panel (post-2001 to align with LD/QU LP sample)
    print("\n[3] Building panel")
    panel = build_panel(delta_resid, "bartik_delta", outcomes)
    panel = panel[panel["quarter_label"] >= "2001Q1"].copy().reset_index(drop=True)
    print(f"  Post-2001 panel: {panel['state_fips'].nunique()} states × "
          f"{panel['quarter_label'].nunique()} quarters = {len(panel):,} rows")

    # [4] Attach interaction terms
    print("\n[4] Attaching interaction terms")

    # Primary: SLOOS large/medium-market C&I
    panel_lm = attach_interaction(panel.copy(), sloos, "sloos_lm", "bartik_delta")
    # Robustness: SLOOS small-firm C&I
    panel_sm = attach_interaction(panel.copy(), sloos, "sloos_sm", "bartik_delta")
    # NFCI (for comparison plot only -- replicate part5 logic)
    if not nfci_df.empty:
        panel_nfci = panel.copy().merge(
            nfci_df, on="quarter_label", how="left"
        )
        panel_nfci["nfci_risk_dm"] = (
            panel_nfci["nfci_risk"] - panel_nfci["nfci_risk"].mean()
        )
        panel_nfci["instr_x_nfci"] = (
            panel_nfci["bartik_delta"] * panel_nfci["nfci_risk_dm"]
        )
    else:
        panel_nfci = pd.DataFrame()

    vac_ok = "vac_rate" in panel.columns
    if not vac_ok:
        print("  [warn] Vacancy data not in panel -- vacancy LPs will be skipped")

    # [5] Unemployment LPs
    print("\n[5] LP: δ × SLOOS_LM → unemployment (primary)")
    irf_lm_u = run_lp(panel_lm, "bartik_delta",
                      "instr_x_sloos_lm", "delta (resid) × SLOOS LM",
                      outcome="unemp")

    print("\n[6] LP: δ × SLOOS_SM → unemployment (robustness)")
    irf_sm_u = run_lp(panel_sm, "bartik_delta",
                      "instr_x_sloos_sm", "delta (resid) × SLOOS SM",
                      outcome="unemp")

    irf_nfci_u = pd.DataFrame()
    if not panel_nfci.empty:
        print("\n[7] LP: δ × NFCI → unemployment (comparison replication)")
        irf_nfci_u = run_lp(panel_nfci, "bartik_delta",
                            "instr_x_nfci", "delta (resid) × NFCI risk",
                            outcome="unemp")

    # [8] Vacancy LPs
    irf_lm_v = irf_sm_v = irf_nfci_v = pd.DataFrame()
    if vac_ok:
        print("\n[8] LP: δ × SLOOS_LM → vacancy (primary)")
        irf_lm_v = run_lp(panel_lm, "bartik_delta",
                          "instr_x_sloos_lm", "delta (resid) × SLOOS LM",
                          outcome="vacancy")

        print("\n[9] LP: δ × SLOOS_SM → vacancy (robustness)")
        irf_sm_v = run_lp(panel_sm, "bartik_delta",
                          "instr_x_sloos_sm", "delta (resid) × SLOOS SM",
                          outcome="vacancy")

        if not panel_nfci.empty:
            print("\n[10] LP: δ × NFCI → vacancy (comparison replication)")
            irf_nfci_v = run_lp(panel_nfci, "bartik_delta",
                                "instr_x_nfci", "delta (resid) × NFCI risk",
                                outcome="vacancy")

    # [11] Save CSVs
    print("\n[11] Saving results")
    def _save(df, fname):
        if df is not None and not df.empty:
            p = RESULTS_DIR / fname
            df.to_csv(p, index=False)
            print(f"  Saved: {p}")

    # unemployment
    _save(irf_lm_u,   "lp_irf_delta_sloos_lm_unemp.csv")
    _save(irf_sm_u,   "lp_irf_delta_sloos_sm_unemp.csv")
    _save(irf_nfci_u, "lp_irf_delta_nfci_rep_unemp.csv")
    # vacancy
    _save(irf_lm_v,   "lp_irf_delta_sloos_lm_vacancy.csv")
    _save(irf_sm_v,   "lp_irf_delta_sloos_sm_vacancy.csv")
    _save(irf_nfci_v, "lp_irf_delta_nfci_rep_vacancy.csv")

    # [12] Comparison plots: SLOOS_LM vs NFCI, unemployment and vacancy
    print("\n[12] Comparison plots")
    plot_sloos_vs_nfci(
        irf_lm_u, irf_nfci_u,
        RESULTS_DIR / "lp_irf_delta_sloos_comparison_unemp.png",
        outcome="unemp"
    )
    if vac_ok:
        plot_sloos_vs_nfci(
            irf_lm_v, irf_nfci_v,
            RESULTS_DIR / "lp_irf_delta_sloos_comparison_vacancy.png",
            outcome="vacancy"
        )

    # [13] Summary table -- both outcomes
    print("\n[13] Summary -- interaction coefficient δ_h at selected horizons")
    for outcome_tag, specs in [
        ("UNEMPLOYMENT", [
            ("SLOOS LM (large/med C&I)", irf_lm_u),
            ("SLOOS SM (small C&I)",     irf_sm_u),
            ("NFCI risk (replication)",  irf_nfci_u),
        ]),
        ("VACANCY", [
            ("SLOOS LM (large/med C&I)", irf_lm_v),
            ("SLOOS SM (small C&I)",     irf_sm_v),
            ("NFCI risk (replication)",  irf_nfci_v),
        ]),
    ]:
        print(f"\n  --- {outcome_tag} ---")
        print(f"  {'Spec':<28}  {'h':>3}  {'beta_h':>8}  {'delta_h':>9}  "
              f"{'dh_SE':>7}  {'dh_p':>6}")
        print(f"  {'-'*65}")
        for lbl, irf in specs:
            if irf is None or irf.empty or "delta_h" not in irf.columns:
                continue
            for h in [0, 4, 8, 12, 16, 20]:
                row = irf[irf["h"] == h]
                if row.empty:
                    continue
                r   = row.iloc[0]
                sig = ("***" if r["delta_h_p"] < .01
                       else "**" if r["delta_h_p"] < .05
                       else "*"  if r["delta_h_p"] < .10 else "")
                print(f"  {lbl:<28}  {h:3d}  {r['beta']:8.4f}  "
                      f"{r['delta_h']:9.4f}  {r['delta_h_se']:7.4f}  "
                      f"{r['delta_h_p']:6.3f}  {sig}")

    print("\nDone.")


if __name__ == "__main__":
    main()