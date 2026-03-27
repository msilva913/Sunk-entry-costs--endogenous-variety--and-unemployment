"""
part5_lp.py -- Panel Local Projections
======================================
Estimates IRFs of state unemployment and vacancy rates to the delta and s
Bartik shocks via Jorda (2005) local projections.

Specification (outcome y, shock k, horizon h):

    y_{s,t+h} - y_{s,t-1} = a_s + a_t + beta_h * B_{s,t}^k
                             + gamma_1 * y_{s,t-1} + gamma_2 * log(LF_{s,t-1})
                             + eps_{s,t,h}

a_s = state FE, a_t = time FE, B^k = leave-one-out Bartik instrument.
State-clustered SEs; outcomes capped at 2019Q4 to exclude COVID.

Instruments
-----------
  Raw:          delta_instrument_base2006.csv   (delta, BED destruction)
                s_instrument_base2006.csv        (s-TS, total separations -- baseline only)
  Residualized: *_resid_base2006.csv (delta, LD, QU) from part3_resid_instruments.py
                residualized on Dlog p_t [+ log theta_t for LD, QU]

Outputs
-------
  lp_irf_delta.csv                       raw delta IRF (full sample)
  lp_irf_s.csv                           raw s-TS IRF (baseline comparison)
  lp_irf_delta_post2001.csv              delta restricted to post-2001
  lp_irf_delta_nfci.csv                  delta with NFCI x B interaction
  lp_irf_s_nfci.csv                      s-TS with NFCI x B interaction
  lp_irf_{delta,ld,qu}_resid.csv         residualized unemployment IRFs
  lp_irf_{delta,ld,qu}_vacancy.csv       residualized vacancy-rate IRFs
  lp_irf_combined.png                    delta vs s-TS baseline
  lp_irf_delta_sample_check.png          sample robustness
  lp_irf_delta_resid_comparison.png      raw vs residualized delta
  lp_irf_delta_ld_qu_overlay.png         asymmetry -- unemployment outcome
  lp_irf_vacancy_decomp_overlay.png      asymmetry -- vacancy outcome
  lp_irf_beveridge_asymmetry.png         delta/LD u-and-v side-by-side panels
  lp_irf_beveridge_path.png              (u_h, v_h) trajectory in UV space

Prerequisites
-------------
  part1_shares.py, part2_shock_rates.py, part2_shock_rates_s.py,
  part3_instrument.py, part3_instrument_s.py, part4_outcomes.py,
  part7_nfci.py, part2b_shock_comovement.py, part3_resid_instruments.py
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
    try:
        os.startfile(path)
    except AttributeError:
        pass

try:
    os.chdir(Path(__file__).resolve().parent)
except (NameError, FileNotFoundError):
    pass

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HORIZONS  = list(range(21))   # h = 0, 1, ..., 20 quarters
BASE_YEAR = 2006
CI_LEVEL  = 0.90
Z90       = 1.645
Z95       = 1.960

# Cap outcome quarter at 2019Q4 to exclude COVID from the outcome window.
MAX_OUTCOME_QUARTER = "2019Q4"

INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DELTA_INSTR_FILE       = INSTR_DIR / f"delta_instrument_base{BASE_YEAR}.csv"
S_INSTR_FILE           = INSTR_DIR / f"s_instrument_base{BASE_YEAR}.csv"
DELTA_RESID_INSTR_FILE = INSTR_DIR / f"delta_instrument_resid_base{BASE_YEAR}.csv"
LD_RESID_INSTR_FILE    = INSTR_DIR / f"ld_instrument_resid_base{BASE_YEAR}.csv"
QU_RESID_INSTR_FILE    = INSTR_DIR / f"qu_instrument_resid_base{BASE_YEAR}.csv"
LAUS_FILE              = INSTR_DIR / "laus_quarterly.parquet"
NFCI_FILE              = Path("data/cache") / "nfci_quarterly.parquet"


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
# Data loading
# ---------------------------------------------------------------------------
def load_instruments():
    delta = pd.read_csv(DELTA_INSTR_FILE, dtype={"state_fips": str})
    s     = pd.read_csv(S_INSTR_FILE,     dtype={"state_fips": str})
    delta["state_fips"] = delta["state_fips"].str.zfill(2)
    s["state_fips"]     = s["state_fips"].str.zfill(2)
    print(f"  delta instrument: {delta['quarter_label'].min()} - "
          f"{delta['quarter_label'].max()}")
    print(f"  s instrument:     {s['quarter_label'].min()} - "
          f"{s['quarter_label'].max()}")
    return delta, s


def load_resid(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"state_fips": str})
    df["state_fips"] = df["state_fips"].str.zfill(2)
    return df


def load_outcomes() -> pd.DataFrame:
    laus = pd.read_parquet(LAUS_FILE)
    laus["state_fips"] = laus["state_fips"].str.zfill(2)
    keep = ["state_fips", "state", "quarter_label", "unemp_rate", "labor_force"]
    if "vacancies" in laus.columns:
        keep.append("vacancies")
        n = laus["vacancies"].notna().sum()
        print(f"  Vacancies: {n:,} non-NaN rows "
              f"({100*n/len(laus):.0f}% coverage)")
    else:
        print("  Vacancies: not present -- run part4_outcomes.py to add")
    return laus[keep]


def load_nfci() -> pd.DataFrame:
    nfci = pd.read_parquet(NFCI_FILE)
    print(f"  NFCI: {nfci['quarter_label'].min()} - "
          f"{nfci['quarter_label'].max()}")
    return nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]]


# ---------------------------------------------------------------------------
# Panel construction
# ---------------------------------------------------------------------------
def build_panel(instr: pd.DataFrame, instr_col: str,
                outcomes: pd.DataFrame) -> pd.DataFrame:
    """
    Merge instrument with outcomes and construct lagged controls.

    Controls are predetermined (lagged to t-1):
        y_{s,t-1}          unemployment rate or vacancy rate
        log(LF_{s,t-1})    log labor force
        vac_rate_{s,t-1}   vacancy rate (when vacancies present)

    vacancy rate = vacancies (thousands) * 1000 / labor_force (persons) * 100
    """
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
        panel["vac_rate"] = panel["vacancies"] * 1000 / panel["labor_force"] * 100
        panel["vac_rate_lag1"] = panel.groupby("state_fips")["vac_rate"].shift(1)

    # Zero out lags where the quarter sequence has a gap
    panel["prev_ql"]  = panel.groupby("state_fips")["quarter_label"].shift(1)
    bad = panel["prev_ql"].notna() & (panel["prev_ql"] != quarter_shift(panel["quarter_label"], -1))
    lag_cols = ["unemp_lag1", "lf_log_lag1"] + (["vac_rate_lag1"] if _has_vac else [])
    panel.loc[bad, lag_cols] = np.nan
    panel.drop(columns="prev_ql", inplace=True)

    print(f"  {instr_col}: {panel['state_fips'].nunique()} states x "
          f"{panel['quarter_label'].nunique()} quarters = {len(panel):,} rows")
    return panel


def attach_nfci_interaction(panel: pd.DataFrame, nfci: pd.DataFrame,
                            instr_col: str) -> pd.DataFrame:
    """
    Merge NFCI and add demeaned interaction B_{s,t} x nfci_risk_dm.
    Demeaning is over the panel's own estimation sample so beta_h is
    interpretable as the IRF at average sample financial conditions.
    """
    panel = panel.merge(
        nfci[["quarter_label", "nfci", "nfci_risk", "anfci"]],
        on="quarter_label", how="left"
    )
    n_miss = panel["nfci_risk"].isna().sum()
    if n_miss:
        print(f"  [warn] {n_miss} rows missing NFCI -- check NFCI date coverage")
    panel["nfci_risk_dm"] = panel["nfci_risk"] - panel["nfci_risk"].mean()
    panel["instr_x_nfci"] = panel[instr_col] * panel["nfci_risk_dm"]
    return panel


# ---------------------------------------------------------------------------
# LP estimation -- single horizon
# ---------------------------------------------------------------------------
def run_lp_horizon(base_panel, h, shock_col,
                   include_nfci=False, outcome="unemp"):
    """
    Estimate the LP at horizon h.

    outcome : "unemp"   -> dep_var = u_{s,t+h} - u_{s,t-1}  (pp)
              "vacancy" -> dep_var = vac_rate_{s,t+h} - vac_rate_{s,t-1}  (pp)
    """
    if outcome not in ("unemp", "vacancy"):
        raise ValueError(f"outcome must be 'unemp' or 'vacancy', got {outcome!r}")

    df = base_panel.copy()
    df["future_ql"] = quarter_shift(df["quarter_label"], h)

    if outcome == "unemp":
        lookup = (df[["state_fips", "quarter_label", "unemp_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "unemp_rate": "y_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["y_future"] - df["unemp_lag1"]
        lag_control   = "unemp_lag1"

    else:  # vacancy
        if "vac_rate" not in df.columns or "vac_rate_lag1" not in df.columns:
            return None
        lookup = (df[["state_fips", "quarter_label", "vac_rate"]]
                  .rename(columns={"quarter_label": "future_ql",
                                   "vac_rate": "vac_future"}))
        df = df.merge(lookup, on=["state_fips", "future_ql"], how="left")
        if MAX_OUTCOME_QUARTER:
            df = df[df["future_ql"] <= MAX_OUTCOME_QUARTER].copy()
        df["dep_var"] = df["vac_future"] - df["vac_rate_lag1"]
        lag_control   = "vac_rate_lag1"

    required = ["dep_var", shock_col, lag_control, "lf_log_lag1"]
    if include_nfci:
        required.append("instr_x_nfci")
    df = df.dropna(subset=required).copy()

    if len(df) < 100:
        return None

    state_dummies = pd.get_dummies(df["state_fips"],    prefix="st",
                                   drop_first=True, dtype=float)
    time_dummies  = pd.get_dummies(df["quarter_label"], prefix="qt",
                                   drop_first=True, dtype=float)
    core = [shock_col, lag_control, "lf_log_lag1"]
    if include_nfci:
        core.append("instr_x_nfci")
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
    row   = {
        "h": h, "beta": beta, "se": se, "tstat": tstat, "pval": pval,
        "partial_f": tstat**2,
        "ci90_lo": beta - Z90*se, "ci90_hi": beta + Z90*se,
        "ci95_lo": beta - Z95*se, "ci95_hi": beta + Z95*se,
        "nobs": int(model.nobs), "r2": float(model.rsquared),
        "n_clusters": df["state_fips"].nunique(),
        "qt_range": f"{df['quarter_label'].min()}-{df['quarter_label'].max()}",
        "outcome": outcome,
    }
    if include_nfci:
        row.update({
            "delta_h":    float(model.params["instr_x_nfci"]),
            "delta_h_se": float(model.bse["instr_x_nfci"]),
            "delta_h_t":  float(model.tvalues["instr_x_nfci"]),
            "delta_h_p":  float(model.pvalues["instr_x_nfci"]),
        })
    return row


# ---------------------------------------------------------------------------
# LP estimation -- all horizons
# ---------------------------------------------------------------------------
def run_lp(panel, shock_col, label, include_nfci=False,
           outcome="unemp") -> pd.DataFrame:
    """Run LP for all horizons, scale to 1-SD units, return DataFrame."""
    outcome_tag = "-> vacancy rate" if outcome == "vacancy" else "-> unemp rate"
    nfci_tag    = " [+NFCI]" if include_nfci else ""
    print(f"\n  LP: {label}{nfci_tag}  {outcome_tag}  h=0..{max(HORIZONS)}")
    if include_nfci:
        print(f"  {'h':>3}  {'beta':>9}  {'SE':>7}  {'t':>7}  {'p':>6}  "
              f"{'delta_h':>9}  {'dh_p':>6}  {'N':>6}")
    else:
        print(f"  {'h':>3}  {'beta':>9}  {'SE':>7}  {'t':>7}  {'p':>6}  "
              f"{'partial-F':>9}  {'N':>6}")

    rows = []
    for h in HORIZONS:
        res = run_lp_horizon(panel, h, shock_col,
                             include_nfci=include_nfci, outcome=outcome)
        if res is None:
            continue
        rows.append(res)
        sig = "***" if res["pval"] < .01 else "**" if res["pval"] < .05 \
              else "*" if res["pval"] < .10 else ""
        if include_nfci:
            print(f"  {h:3d}  {res['beta']:9.4f}  {res['se']:7.4f}  "
                  f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
                  f"{res['delta_h']:9.4f}  {res['delta_h_p']:6.3f}  "
                  f"{res['nobs']:6d}  {sig}")
        else:
            print(f"  {h:3d}  {res['beta']:9.4f}  {res['se']:7.4f}  "
                  f"{res['tstat']:7.3f}  {res['pval']:6.3f}  "
                  f"{res['partial_f']:9.2f}  {res['nobs']:6d}  {sig}")

    irf = pd.DataFrame(rows)
    if irf.empty:
        return irf
    sd = float(panel[shock_col].std())
    for col in ["beta", "se", "ci90_lo", "ci90_hi", "ci95_lo", "ci95_hi"]:
        irf[col] = irf[col] * sd
    if include_nfci:
        irf["delta_h"]    = irf["delta_h"]    * sd
        irf["delta_h_se"] = irf["delta_h_se"] * sd
    irf["instr_sd"] = sd
    irf["outcome"]  = outcome
    return irf


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
_PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
_LINESTYLES = ["-", "-", "-", "--", ":", "-."]


def plot_irf(results: dict, out_path: Path) -> None:
    """Side-by-side IRF panels. Each DataFrame already 1-SD scaled."""
    n    = len(results)
    fig, axes = plt.subplots(1, n, figsize=(7*n, 6), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, (label, df), color in zip(axes, results.items(), _PALETTE):
        if df.empty:
            ax.set_visible(False)
            continue
        h, beta = df["h"].values, df["beta"].values
        ax.fill_between(h, df["ci90_lo"], df["ci90_hi"], color=color, alpha=0.20)
        ax.plot(h, df["ci95_lo"], color=color, lw=0.7, ls="--", alpha=0.6)
        ax.plot(h, df["ci95_hi"], color=color, lw=0.7, ls="--", alpha=0.6)
        ax.plot(h, beta, color=color, lw=2.0, marker="o", ms=3.5, label=label)
        ax.axhline(0, color="black", lw=0.8)
        for hh in range(4, int(h.max())+1, 4):
            ax.axvline(hh, color="grey", lw=0.5, ls=":", alpha=0.6)
        _outcome = df["outcome"].iloc[0] if "outcome" in df.columns else "unemp"
        _ylabel  = ("pp change in vacancy rate\nper 1-SD shock"
                    if _outcome == "vacancy"
                    else "pp change in unemp. rate\nper 1-SD shock")
        ax.set_title(f"IRF -- {label}", fontsize=11)
        ax.set_xlabel("Horizon h (quarters)", fontsize=10)
        ax.set_ylabel(_ylabel, fontsize=10)
        ax.set_xticks(h[::2] if h.max() > 16 else h)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: str(int(x))))
        ax.legend(fontsize=9, framealpha=0.85)
        ax.grid(axis="y", lw=0.4, alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


def overlay_plot_irf(results: dict, out_path: Path,
                     title: str = "", ylabel: str = "") -> None:
    """Multiple IRFs overlaid on one panel for direct shape comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(0, color="black", lw=0.8)
    for hh in [4, 8, 12, 16, 20]:
        ax.axvline(hh, color="grey", lw=0.4, ls=":", alpha=0.5)

    for idx, (label, df) in enumerate(results.items()):
        if df is None or df.empty:
            continue
        color, ls = _PALETTE[idx], _LINESTYLES[idx]
        h = df["h"].values
        ax.fill_between(h, df["ci90_lo"], df["ci90_hi"], color=color, alpha=0.12)
        ax.plot(h, df["beta"], color=color, lw=2.0, ls=ls,
                marker="o", ms=3.0, label=label)

    ax.set_xlabel("Horizon h (quarters)", fontsize=11)
    ax.set_ylabel(ylabel or "pp change per 1-SD shock", fontsize=11)
    ax.set_title(title or "IRF comparison", fontsize=12)
    ax.legend(fontsize=10, framealpha=0.85)
    ax.grid(axis="y", lw=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


def beveridge_path_plot(shocks: dict, out_path: Path,
                        label_horizons: list | None = None) -> None:
    """
    Parametric (u_h, v_h) trajectory in unemployment-vacancy space.

    shocks maps label -> (unemp_irf, vac_irf) DataFrames from run_lp().
    x = beta_h^u (pp change in unemployment rate per 1-SD shock)
    y = beta_h^v (pp change in vacancy rate per 1-SD shock)
    Origin (0,0) is the pre-shock state at t-1; h increases along the path.

    Theory: delta -> path moves southeast (u up, v down; Beveridge curve shift)
            LD   -> path moves more due-east (u up, v roughly flat; along curve)
    """
    if label_horizons is None:
        label_horizons = [0, 4, 8, 12, 16, 20]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.axhline(0, color="black", lw=0.7, alpha=0.5)
    ax.axvline(0, color="black", lw=0.7, alpha=0.5)
    ax.scatter([0], [0], color="black", s=50, zorder=6, label="pre-shock (t-1)")

    for idx, (label, (u_irf, v_irf)) in enumerate(shocks.items()):
        if u_irf is None or u_irf.empty or v_irf is None or v_irf.empty:
            continue
        color = _PALETTE[idx]
        u = u_irf.set_index("h")["beta"]
        v = v_irf.set_index("h")["beta"]
        hs  = sorted(set(u.index) & set(v.index))
        uv  = [(u[h], v[h]) for h in hs]
        us  = np.array([x for x, _ in uv])
        vs  = np.array([y for _, y in uv])

        # Path line
        ax.plot(us, vs, color=color, lw=2.0, label=label, zorder=3)

        # Directional arrows every 4 steps
        for i in range(0, len(hs) - 1, 4):
            ax.annotate("",
                        xy=(us[i+1], vs[i+1]),
                        xytext=(us[i], vs[i]),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        lw=1.3, mutation_scale=12),
                        zorder=4)

        # Labelled markers at selected horizons
        for h_lab in label_horizons:
            if h_lab in hs:
                i = hs.index(h_lab)
                ax.scatter([us[i]], [vs[i]], color=color, s=30, zorder=5)
                ax.annotate(f"h={h_lab}", xy=(us[i], vs[i]),
                            fontsize=7.5, color=color,
                            xytext=(4, 3), textcoords="offset points")

    ax.set_xlabel("pp change in unemployment rate (per 1-SD shock)", fontsize=11)
    ax.set_ylabel("pp change in vacancy rate (per 1-SD shock)", fontsize=11)
    ax.set_title(
        "Beveridge path: joint (u, v) response by shock type\n"
        "southeast = Beveridge curve shift; due-east = movement along curve",
        fontsize=11
    )
    ax.legend(fontsize=10, framealpha=0.85)
    ax.grid(lw=0.4, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
    _open_file(out_path)


# ===========================================================================
# Main
# ===========================================================================
print("\n" + "=" * 60)
print("Part 5 -- Panel Local Projections")
print("=" * 60)

# ---------------------------------------------------------------------------
# [1] Load data
# ---------------------------------------------------------------------------
print("\n[1] Loading data")
delta_instr, s_instr = load_instruments()
outcomes = load_outcomes()
nfci     = load_nfci()

resid_ok = all(p.exists() for p in [DELTA_RESID_INSTR_FILE,
                                     LD_RESID_INSTR_FILE,
                                     QU_RESID_INSTR_FILE])
vac_ok   = ("vacancies" in outcomes.columns
            and outcomes["vacancies"].notna().any())
print(f"  Residualized instruments available: {resid_ok}")
print(f"  Vacancy data available:             {vac_ok}")

# ---------------------------------------------------------------------------
# [2] Build panels -- one pass; panels are reused for all LP outcomes
# ---------------------------------------------------------------------------
print("\n[2] Building panels")

# Raw panels
delta_panel    = build_panel(delta_instr, "bartik_delta", outcomes)
s_panel        = build_panel(s_instr,     "bartik_s",     outcomes)

# Subset delta to post-2001 BEFORE attaching NFCI so demeaning uses
# the post-2001 sample mean (keeps beta_h interpretation clean).
delta_post2001 = (delta_panel[delta_panel["quarter_label"] >= "2001Q1"]
                  .copy().reset_index(drop=True))

# Attach NFCI interaction to each panel (demean within own sample)
delta_panel    = attach_nfci_interaction(delta_panel,    nfci, "bartik_delta")
s_panel        = attach_nfci_interaction(s_panel,        nfci, "bartik_s")
delta_post2001 = attach_nfci_interaction(delta_post2001, nfci, "bartik_delta")

# Residualized panels (delta, LD, QU) -- used for both unemployment and
# vacancy LPs; no need to rebuild in the vacancy section.
if resid_ok:
    delta_resid_panel = build_panel(load_resid(DELTA_RESID_INSTR_FILE),
                                    "bartik_delta", outcomes)
    ld_resid_panel    = build_panel(load_resid(LD_RESID_INSTR_FILE),
                                    "bartik_ld",    outcomes)
    qu_resid_panel    = build_panel(load_resid(QU_RESID_INSTR_FILE),
                                    "bartik_qu",    outcomes)

    # delta restricted to post-2001 for comparability with LD/QU
    delta_resid_post  = (delta_resid_panel[
                             delta_resid_panel["quarter_label"] >= "2001Q1"]
                         .copy().reset_index(drop=True))

    delta_resid_post = attach_nfci_interaction(delta_resid_post, nfci, "bartik_delta")
    ld_resid_panel   = attach_nfci_interaction(ld_resid_panel,   nfci, "bartik_ld")
    qu_resid_panel   = attach_nfci_interaction(qu_resid_panel,   nfci, "bartik_qu")

# ---------------------------------------------------------------------------
# [3] Raw baseline LPs
# ---------------------------------------------------------------------------
print("\n[3] Raw baseline LPs")

irf_delta      = run_lp(delta_panel,    "bartik_delta", "delta (raw, full)")
irf_s          = run_lp(s_panel,        "bartik_s",     "s-TS (raw)")
irf_delta_post = run_lp(delta_post2001, "bartik_delta", "delta (raw, post-2001)")
irf_delta_nfci = run_lp(delta_panel,    "bartik_delta", "delta (raw)",
                         include_nfci=True)
irf_s_nfci     = run_lp(s_panel,        "bartik_s",     "s-TS (raw)",
                         include_nfci=True)

# ---------------------------------------------------------------------------
# [4] Residualized unemployment LPs
# ---------------------------------------------------------------------------
if resid_ok:
    print("\n[4] Residualized unemployment LPs")
    irf_delta_resid = run_lp(delta_resid_post, "bartik_delta",
                              "delta (resid, post-2001)")
    irf_ld_resid    = run_lp(ld_resid_panel,   "bartik_ld",
                              "LD (resid)")
    irf_qu_resid    = run_lp(qu_resid_panel,   "bartik_qu",
                              "QU (resid)")
else:
    print("\n[4] Residualized instruments not found -- skipping.")
    print("    Run part2b_shock_comovement.py then part3_resid_instruments.py")

# ---------------------------------------------------------------------------
# [5] Residualized vacancy LPs -- reuse panels built in [2]
# ---------------------------------------------------------------------------
if resid_ok and vac_ok:
    print("\n[5] Residualized vacancy LPs")
    irf_delta_vac = run_lp(delta_resid_post, "bartik_delta",
                            "delta -> vacancy", outcome="vacancy")
    irf_ld_vac    = run_lp(ld_resid_panel,   "bartik_ld",
                            "LD -> vacancy",    outcome="vacancy")
    irf_qu_vac    = run_lp(qu_resid_panel,   "bartik_qu",
                            "QU -> vacancy",    outcome="vacancy")
elif vac_ok and not resid_ok:
    print("\n[5] Vacancy LPs skipped -- residualized instruments required.")
else:
    print("\n[5] Vacancy LPs skipped -- vacancy data not available.")

# ---------------------------------------------------------------------------
# [6] Save results
# ---------------------------------------------------------------------------
print("\n[6] Saving results")

def _save(irf, fname):
    if irf is not None and not irf.empty:
        p = RESULTS_DIR / fname
        irf.to_csv(p, index=False)
        print(f"  Saved: {p}")

_save(irf_delta,      "lp_irf_delta.csv")
_save(irf_s,          "lp_irf_s.csv")
_save(irf_delta_post, "lp_irf_delta_post2001.csv")
_save(irf_delta_nfci, "lp_irf_delta_nfci.csv")
_save(irf_s_nfci,     "lp_irf_s_nfci.csv")

if resid_ok:
    _save(irf_delta_resid, "lp_irf_delta_resid.csv")
    _save(irf_ld_resid,    "lp_irf_ld_resid.csv")
    _save(irf_qu_resid,    "lp_irf_qu_resid.csv")

if resid_ok and vac_ok:
    _save(irf_delta_vac, "lp_irf_delta_vacancy.csv")
    _save(irf_ld_vac,    "lp_irf_ld_vacancy.csv")
    _save(irf_qu_vac,    "lp_irf_qu_vacancy.csv")

# ---------------------------------------------------------------------------
# [7] Plots
# ---------------------------------------------------------------------------
print("\n[7] Plotting")

# 7a. Baseline: delta (full) vs s-TS raw
plot_irf({"delta (raw)": irf_delta, "s-TS (raw)": irf_s},
         RESULTS_DIR / "lp_irf_combined.png")

# 7b. Sample robustness: delta full vs post-2001
plot_irf({"delta (full)": irf_delta, "delta (post-2001)": irf_delta_post},
         RESULTS_DIR / "lp_irf_delta_sample_check.png")

if resid_ok:
    # 7c. Raw vs residualized delta
    plot_irf({"delta (raw, post-2001)": irf_delta_post,
              "delta (resid)":          irf_delta_resid},
             RESULTS_DIR / "lp_irf_delta_resid_comparison.png")

    # 7d. Core asymmetry -- unemployment outcome
    overlay_plot_irf(
        {"delta (resid)": irf_delta_resid,
         "LD (resid)":    irf_ld_resid,
         "QU (resid)":    irf_qu_resid},
        out_path = RESULTS_DIR / "lp_irf_delta_ld_qu_overlay.png",
        title    = (r"IRF: $\delta$ vs Layoffs+Discharges vs Quits"
                    "\n(productivity-residualized, 1-SD scale, post-2001)"),
        ylabel   = "pp change in unemp. rate per 1-SD shock",
    )

if resid_ok and vac_ok:
    # 7e. Core asymmetry -- vacancy outcome
    overlay_plot_irf(
        {"delta (resid)": irf_delta_vac,
         "LD (resid)":    irf_ld_vac,
         "QU (resid)":    irf_qu_vac},
        out_path = RESULTS_DIR / "lp_irf_vacancy_decomp_overlay.png",
        title    = (r"Vacancy IRF: $\delta$ vs LD vs QU"
                    "\n(productivity-residualized, 1-SD scale, post-2001)"),
        ylabel   = "pp change in vacancy rate per 1-SD shock",
    )

    # 7f. Beveridge asymmetry -- u and v side-by-side panels for delta and LD
    plot_irf(
        {"delta -- unemployment": irf_delta_resid,
         "delta -- vacancies":   irf_delta_vac,
         "LD -- unemployment":   irf_ld_resid,
         "LD -- vacancies":      irf_ld_vac},
        RESULTS_DIR / "lp_irf_beveridge_asymmetry.png",
    )

    # 7g. Beveridge path -- (u_h, v_h) trajectory in UV space
    beveridge_path_plot(
        {"delta (resid)": (irf_delta_resid, irf_delta_vac),
         "LD (resid)":    (irf_ld_resid,    irf_ld_vac)},
        out_path = RESULTS_DIR / "lp_irf_beveridge_path.png",
    )

# ---------------------------------------------------------------------------
# [8] Summary table
# ---------------------------------------------------------------------------
print("\n[8] Summary -- peak responses (1-SD standardized)")
print(f"  {'Series':<35}  {'outcome':>8}  {'h_peak':>6}  "
      f"{'beta_peak':>9}  {'SE':>7}  {'p':>6}")
print(f"  {'-'*78}")

all_irfs = [
    ("delta (raw, full)",          "unemp",   irf_delta),
    ("delta (raw, post-2001)",     "unemp",   irf_delta_post),
    ("s-TS (raw)",                 "unemp",   irf_s),
]
if resid_ok:
    all_irfs += [
        ("delta (resid, post-2001)",   "unemp",   irf_delta_resid),
        ("LD (resid)",                 "unemp",   irf_ld_resid),
        ("QU (resid)",                 "unemp",   irf_qu_resid),
    ]
if resid_ok and vac_ok:
    all_irfs += [
        ("delta (resid)",              "vacancy", irf_delta_vac),
        ("LD (resid)",                 "vacancy", irf_ld_vac),
        ("QU (resid)",                 "vacancy", irf_qu_vac),
    ]

for lbl, outcome_tag, df in all_irfs:
    if df is None or df.empty:
        continue
    # Peak by absolute value (vacancies have negative peaks)
    pk = df.loc[df["beta"].abs().idxmax()]
    print(f"  {lbl:<35}  {outcome_tag:>8}  {int(pk['h']):>6}  "
          f"{pk['beta']:>9.4f}  {pk['se']:>7.4f}  {pk['pval']:>6.3f}")
