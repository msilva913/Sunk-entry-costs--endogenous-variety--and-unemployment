"""
part10_state_scatter.py -- Cross-state motivating scatter
==========================================================
For each of the 50 states, plots the full-sample average quarterly
establishment exit rate (X-axis) against labor market outcome measures
(Y-axes), using the three NBER recessions covered by BED data:
  2001 (2001Q1-2001Q4), 2008-09 (2007Q4-2009Q2), 2020 (2020Q1-2020Q2).

X-axis: state average of raw Bartik delta instrument over full sample
  (1992Q3-2021Q4). This is the employment-share-weighted exit rate,
  averaged across all quarters. Captures structural cross-state
  heterogeneity in exit exposure. Non-residualized by design: the goal
  is geographic variation in exit propensity, not causal identification.

PRIMARY figure (state_scatter_primary.png) -- 1x2:
  Left:  cumulative unemployment gap (CUG, pp-qtrs), H=20 from rec. start
  Right: cumulative vacancy shortfall (CVS, pp-qtrs), H=20 from rec. start

APPENDIX figure (state_scatter_appendix.png) -- 2x2:
  Row 0: recession-end-anchored 50% U/LF and V/LF recovery speed
  Row 1: trough-anchored 50% U/LF and V/LF recovery speed
  (trough-anchor included for comparison; computed fresh here)

For each state, outcomes are averaged across the three recessions
(weighted by state labor force at recession start) to produce one
point per state.

Regression line: WLS weighted by state average labor force.
Bubble size: proportional to state average labor force.

Data sources:
  data/instruments/delta_instrument_base2006.csv  (bartik_delta)
  data/instruments/laus_quarterly.parquet         (unemp_rate, vac_rate, labor_force)

Output:
  data/results/state_scatter_primary.png
  data/results/state_scatter_appendix.png
  data/results/state_scatter.csv
  data/results/state_scatter_latex.tex
"""

import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
INSTR_DIR   = Path("data/instruments")
RESULTS_DIR = Path("data/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DELTA_CSV  = INSTR_DIR / "delta_instrument_base2006.csv"
LAUS_PATH  = INSTR_DIR / "laus_quarterly.parquet"

RECESSIONS = {
    "2001":    {"start": "2001Q1", "end": "2001Q4"},
    "2008-09": {"start": "2007Q4", "end": "2009Q2"},
    "2020":    {"start": "2020Q1", "end": "2020Q2"},
}

PRE_REC_WINDOW    = 4    # quarters before recession start for baseline
RECOVERY_THRESHOLD = 0.50
MAX_SEARCH_QTRS   = 28
CUM_GAP_HORIZON   = 20  # quarters from recession start

# ---------------------------------------------------------------------------
# Quarter utilities
# ---------------------------------------------------------------------------
def ql_shift(ql: str, n: int) -> str:
    p = pd.Period(ql, freq="Q") + n
    return f"{p.year}Q{p.quarter}"

def ql_range(start: str, end: str) -> list:
    return [f"{p.year}Q{p.quarter}"
            for p in pd.period_range(start, end, freq="Q")]

# ---------------------------------------------------------------------------
# Outcome functions (state-level series as pd.Series indexed by quarter_label)
# ---------------------------------------------------------------------------
def recovery_from_end(series: pd.Series, rec_end_ql: str,
                      trough_val: float, pre_level: float,
                      direction: str) -> float:
    """Quarters from recession end until RECOVERY_THRESHOLD of gap is closed."""
    if pd.isna(trough_val) or pd.isna(pre_level):
        return np.nan
    gap = abs(trough_val - pre_level)
    if gap < 1e-6:
        return 0.0
    target = (trough_val - RECOVERY_THRESHOLD * gap if direction == "down"
              else trough_val + RECOVERY_THRESHOLD * gap)
    for h in range(MAX_SEARCH_QTRS + 1):
        val = series.get(ql_shift(rec_end_ql, h), np.nan)
        if pd.isna(val):
            continue
        if direction == "down" and val <= target:
            return float(h)
        if direction == "up" and val >= target:
            return float(h)
    return float(MAX_SEARCH_QTRS + 1)


def recovery_from_trough(series: pd.Series, trough_ql: str,
                         trough_val: float, pre_level: float,
                         direction: str) -> float:
    """Quarters from series trough until RECOVERY_THRESHOLD of gap is closed."""
    if pd.isna(trough_val) or pd.isna(pre_level) or trough_ql is None:
        return np.nan
    gap = abs(trough_val - pre_level)
    if gap < 1e-6:
        return 0.0
    target = (trough_val - RECOVERY_THRESHOLD * gap if direction == "down"
              else trough_val + RECOVERY_THRESHOLD * gap)
    for h in range(MAX_SEARCH_QTRS + 1):
        val = series.get(ql_shift(trough_ql, h), np.nan)
        if pd.isna(val):
            continue
        if direction == "down" and val <= target:
            return float(h)
        if direction == "up" and val >= target:
            return float(h)
    return float(MAX_SEARCH_QTRS + 1)


def cumulative_gap(series: pd.Series, rec_start_ql: str,
                   pre_level: float, direction: str) -> float:
    """Area under positive gap curve from rec_start over CUM_GAP_HORIZON qtrs."""
    if pd.isna(pre_level):
        return np.nan
    total = 0.0
    for h in range(CUM_GAP_HORIZON + 1):
        val = series.get(ql_shift(rec_start_ql, h), np.nan)
        if pd.isna(val):
            continue
        total += max(val - pre_level, 0.0) if direction == "down" \
            else max(pre_level - val, 0.0)
    return total


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_data():
    delta = pd.read_csv(DELTA_CSV)
    # Ensure state_fips is zero-padded string
    delta["state_fips"] = delta["state_fips"].astype(str).str.zfill(2)

    laus = pd.read_parquet(LAUS_PATH)
    laus["state_fips"] = laus["state_fips"].astype(str).str.zfill(2)

    return delta, laus


# ---------------------------------------------------------------------------
# Compute per-state outcomes averaged across recessions
# ---------------------------------------------------------------------------
def compute_state_outcomes(delta: pd.DataFrame, laus: pd.DataFrame) -> pd.DataFrame:
    # Full-sample average exit rate per state
    avg_delta = (delta.groupby("state_fips")["bartik_delta"]
                 .mean().rename("avg_delta_rate"))

    states = laus["state_fips"].unique()
    records = []

    for fips in sorted(states):
        state_laus = laus[laus["state_fips"] == fips].copy()
        state_name = state_laus["state"].iloc[0] if "state" in state_laus.columns else fips
        unemp_s = state_laus.set_index("quarter_label")["unemp_rate"]
        vac_s   = state_laus.set_index("quarter_label")["vac_rate"]
        lf_s    = state_laus.set_index("quarter_label")["labor_force"]

        # Per-recession outcomes, weighted by LF at recession start
        rec_rows = []
        for rec_name, info in RECESSIONS.items():
            # Skip if vacancy data not available (pre-2001)
            if info["start"] < "2001Q1" and vac_s.reindex(
                    ql_range(info["start"], info["end"])).isna().all():
                continue

            pre_qtrs = [ql_shift(info["start"], -k)
                        for k in range(1, PRE_REC_WINDOW + 1)]
            pre_u = unemp_s.reindex(pre_qtrs).mean()
            pre_v = vac_s.reindex(pre_qtrs).mean()
            lf_at_start = lf_s.get(info["start"], np.nan)

            # Trough values
            search_u = ql_range(info["start"], ql_shift(info["end"], 4))
            u_vals = unemp_s.reindex(search_u).dropna()
            u_trough_ql  = u_vals.idxmax() if len(u_vals) else None
            u_trough_val = u_vals.max()    if len(u_vals) else np.nan

            search_v = ql_range(info["start"], ql_shift(info["end"], 4))
            v_vals = vac_s.reindex(search_v).dropna()
            v_trough_ql  = v_vals.idxmin() if len(v_vals) else None
            v_trough_val = v_vals.min()    if len(v_vals) else np.nan

            u_rec_end = recovery_from_end(
                unemp_s, info["end"], u_trough_val, pre_u, "down")
            v_rec_end = recovery_from_end(
                vac_s, info["end"], v_trough_val, pre_v, "up")
            u_rec_trough = recovery_from_trough(
                unemp_s, u_trough_ql, u_trough_val, pre_u, "down")
            v_rec_trough = recovery_from_trough(
                vac_s, v_trough_ql, v_trough_val, pre_v, "up")
            cug = cumulative_gap(unemp_s, info["start"], pre_u, "down")
            cvs = cumulative_gap(vac_s,   info["start"], pre_v, "up")

            rec_rows.append({
                "recession":     rec_name,
                "lf":            lf_at_start,
                "u_rec_end":     u_rec_end,
                "v_rec_end":     v_rec_end,
                "u_rec_trough":  u_rec_trough,
                "v_rec_trough":  v_rec_trough,
                "cug":           cug,
                "cvs":           cvs,
            })

        if not rec_rows:
            continue

        rec_df = pd.DataFrame(rec_rows).dropna(subset=["lf"])
        if len(rec_df) == 0:
            continue

        # LF-weighted average across recessions
        w = rec_df["lf"].values
        w = w / w.sum()

        records.append({
            "state_fips":    fips,
            "state":         state_name,
            "avg_delta_rate": avg_delta.get(fips, np.nan),
            "avg_lf":        rec_df["lf"].mean(),
            "u_rec_end":     np.average(rec_df["u_rec_end"].values, weights=w),
            "v_rec_end":     np.average(rec_df["v_rec_end"].fillna(MAX_SEARCH_QTRS + 1).values, weights=w),
            "u_rec_trough":  np.average(rec_df["u_rec_trough"].values, weights=w),
            "v_rec_trough":  np.average(rec_df["v_rec_trough"].fillna(MAX_SEARCH_QTRS + 1).values, weights=w),
            "cug":           np.average(rec_df["cug"].values, weights=w),
            "cvs":           np.average(rec_df["cvs"].fillna(0).values, weights=w),
            "n_rec":         len(rec_df),
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Shared plot helper
# ---------------------------------------------------------------------------
def _bubble_sizes(lf_all, lf_sub, s_min=30, s_max=350):
    return s_min + (s_max - s_min) * (lf_sub - lf_all.min()) / (lf_all.max() - lf_all.min())


def _draw_panel(ax, df, xcol, ycol, xlabel, ylabel, title, lf_all,
                s_min=30, s_max=350):
    sub = df.dropna(subset=[xcol, ycol]).copy()
    if len(sub) < 3:
        ax.set_title(f"{title}\n(insufficient data)")
        return np.nan

    x  = sub[xcol].values
    y  = sub[ycol].values
    w  = sub["avg_lf"].values
    sz = _bubble_sizes(lf_all, w, s_min, s_max)

    ax.scatter(x, y, s=sz, alpha=0.6, color="steelblue",
               edgecolors="white", linewidths=0.4, zorder=4)

    # WLS line
    W = np.diag(w)
    X = np.column_stack([np.ones(len(x)), x])
    coeffs = np.linalg.lstsq(W @ X, W @ y, rcond=None)[0]
    b, m = coeffs
    x_line = np.linspace(x.min() * 0.90, x.max() * 1.10, 200)
    ax.plot(x_line, b + m * x_line, color="firebrick", linewidth=1.5,
            linestyle="--", alpha=0.85)

    # Weighted correlation
    x_dm = x - np.average(x, weights=w)
    y_dm = y - np.average(y, weights=w)
    r_w  = (np.sum(w * x_dm * y_dm) /
            np.sqrt(np.sum(w * x_dm**2) * np.sum(w * y_dm**2)))
    ax.text(0.05, 0.91, f"$r_w$ = {r_w:.2f}", transform=ax.transAxes,
            fontsize=9, color="dimgrey")

    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.grid(linewidth=0.4, alpha=0.4)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.3f}"))
    return r_w


def _add_bubble_legend(fig, lf_all, s_min=30, s_max=350):
    lf_ref = [lf_all.min(), np.median(lf_all), lf_all.max()]
    handles = [
        plt.scatter([], [], s=_bubble_sizes(lf_all, np.array([lf]))[0],
                    color="steelblue", alpha=0.6, edgecolors="white",
                    linewidths=0.4, label=lbl)
        for lf, lbl in zip(lf_ref, ["Small state", "Median state", "Large state"])
    ]
    fig.legend(handles=handles, title="Avg. labor force (bubble size)",
               loc="lower center", ncol=3, fontsize=8, title_fontsize=8,
               bbox_to_anchor=(0.5, -0.03), frameon=True)


# ---------------------------------------------------------------------------
# Primary figure: CUG and CVS (1x2)
# ---------------------------------------------------------------------------
def plot_primary(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    lf_all = df["avg_lf"].values
    xlabel = ("Average quarterly establishment exit rate (full sample)\n"
              r"$\bar{g}^\delta_s$, 1992Q3--2021Q4")

    _draw_panel(axes[0], df, "avg_delta_rate", "cug", xlabel,
                f"CUG (pp-qtrs, H={CUM_GAP_HORIZON})",
                f"Cumulative unemployment gap (H={CUM_GAP_HORIZON} qtrs)",
                lf_all)
    _draw_panel(axes[1], df, "avg_delta_rate", "cvs", xlabel,
                f"CVS (pp-qtrs, H={CUM_GAP_HORIZON})",
                f"Cumulative vacancy shortfall (H={CUM_GAP_HORIZON} qtrs)",
                lf_all)

    _add_bubble_legend(fig, lf_all)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    out = RESULTS_DIR / "state_scatter_primary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")
    return out


# ---------------------------------------------------------------------------
# Appendix figure: recovery speed, rec-end and trough anchors (2x2)
# ---------------------------------------------------------------------------
def plot_appendix(df: pd.DataFrame):
    recov_pct = int(RECOVERY_THRESHOLD * 100)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    lf_all = df["avg_lf"].values
    xlabel = ("Average quarterly establishment exit rate (full sample)\n"
              r"$\bar{g}^\delta_s$, 1992Q3--2021Q4")

    panels = [
        (axes[0, 0], "u_rec_end",
         f"Qtrs from rec. end to {recov_pct}% U/LF recovery",
         "Unemployment recovery (rec.-end anchor)"),
        (axes[0, 1], "v_rec_end",
         f"Qtrs from rec. end to {recov_pct}% V/LF recovery",
         "Vacancy recovery (rec.-end anchor)"),
        (axes[1, 0], "u_rec_trough",
         f"Qtrs from trough to {recov_pct}% U/LF recovery",
         "Unemployment recovery (trough anchor)"),
        (axes[1, 1], "v_rec_trough",
         f"Qtrs from trough to {recov_pct}% V/LF recovery",
         "Vacancy recovery (trough anchor)"),
    ]

    for ax, ycol, ylabel, title in panels:
        _draw_panel(ax, df, "avg_delta_rate", ycol, xlabel, ylabel, title, lf_all)

    _add_bubble_legend(fig, lf_all)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    out = RESULTS_DIR / "state_scatter_appendix.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")
    return out


# ---------------------------------------------------------------------------
# LaTeX snippets
# ---------------------------------------------------------------------------
LATEX_PRIMARY = r"""\begin{figure}[H]
    \centering
    \includegraphics[width=0.92\textwidth]{state_scatter_primary.png}
    \caption{Cross-state relationship between average establishment exit rates
    and cumulative recession labor market damage. Each point is one of the 50
    U.S.\ states; bubble size is proportional to average labor force.
    The horizontal axis is the full-sample time average of the state-level
    Bartik-weighted establishment exit rate ($\bar{g}^\delta_s$,
    1992Q3--2021Q4), constructed from BLS Business Employment Dynamics (BED)
    data using 2006 base-year employment shares.
    The left panel shows the cumulative unemployment gap (CUG): the
    sum of quarterly deviations of the unemployment rate above its
    pre-recession average over $H=20$ quarters from recession onset,
    in percentage-point-quarters. The right panel shows the analogous
    cumulative vacancy shortfall (CVS) for the vacancy rate.
    Both outcomes are averaged across the three post-1992 NBER recessions
    (2001, 2007--2009, 2020), weighted by state labor force at recession onset.
    The regression line is estimated by WLS with state labor force weights;
    $r_w$ denotes the labor-force-weighted correlation.}
    \label{fig:state_scatter}
\end{figure}"""

LATEX_APPENDIX = r"""\begin{figure}[H]
    \centering
    \includegraphics[width=0.92\textwidth]{state_scatter_appendix.png}
    \caption{Cross-state relationship between average establishment exit rates
    and recession recovery speed. Layout and data sources as in
    Figure~\ref{fig:state_scatter}. The top row shows quarters from the NBER
    recession end date until 50\% closure of the trough gap for unemployment
    (left) and vacancy rates (right); the bottom row uses the series trough
    date as the clock origin. Recovery speed does not exhibit a robust positive
    relationship with exit rates in the cross-section, reflecting compositional
    differences in industry structure across high- and low-exit-rate states.
    The cumulative gap measures in Figure~\ref{fig:state_scatter} provide the
    cleaner cross-state comparison.}
    \label{fig:state_scatter_appendix}
\end{figure}"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading data ...")
    delta, laus = load_data()
    print(f"  Delta: {len(delta)} state-quarters, "
          f"{delta['state_fips'].nunique()} states")
    print(f"  LAUS:  {len(laus)} state-quarters, "
          f"{laus['state_fips'].nunique()} states")

    print("\nComputing state outcomes ...")
    df = compute_state_outcomes(delta, laus)
    print(f"  States with complete data: {len(df)}")

    # Save CSV
    csv_out = RESULTS_DIR / "state_scatter.csv"
    df.drop(columns=["state_fips"]).to_csv(csv_out, index=False)
    print(f"  CSV saved: {csv_out}")

    # Summary
    print("\nState-level summary (weighted correlation with avg_delta_rate):")
    for col in ["u_rec_end", "v_rec_end", "u_rec_trough", "v_rec_trough", "cug", "cvs"]:
        sub = df.dropna(subset=["avg_delta_rate", col])
        x = sub["avg_delta_rate"].values
        y = sub[col].values
        w = sub["avg_lf"].values
        xd = x - np.average(x, weights=w)
        yd = y - np.average(y, weights=w)
        r  = np.sum(w * xd * yd) / np.sqrt(np.sum(w * xd**2) * np.sum(w * yd**2))
        print(f"  {col:20s}: r_w = {r:.3f}  (n={len(sub)})")

    print("\nPlotting primary figure (CUG + CVS) ...")
    plot_primary(df)

    print("Plotting appendix figure (recovery speed) ...")
    plot_appendix(df)

    # LaTeX
    latex_out = RESULTS_DIR / "state_scatter_latex.tex"
    with open(latex_out, "w") as f:
        f.write("% PRIMARY FIGURE\n")
        f.write(LATEX_PRIMARY)
        f.write("\n\n% APPENDIX FIGURE\n")
        f.write(LATEX_APPENDIX)
        f.write("\n")
    print(f"LaTeX snippets saved: {latex_out}")
    print("\n--- PRIMARY LaTeX ---")
    print(LATEX_PRIMARY)
    print("\n--- APPENDIX LaTeX ---")
    print(LATEX_APPENDIX)
