"""
part10b_state_scatter_sensitivity.py -- X-axis sensitivity for recovery speed
==============================================================================
Compares two definitions of the state-level establishment exit rate on the
X-axis for the recovery speed panels only (U and V, recession-end anchor):

  Full-sample average: time-average of bartik_delta over 1992Q3-2021Q4.
  Recession-quarter average: average of bartik_delta restricted to the 13
    quarters spanning the three NBER recessions (2001Q1-2001Q4,
    2007Q4-2009Q2, 2020Q1-2020Q2).

The two measures have a cross-state correlation of r=0.995 and a uniform
scaling ratio of ~1.18 (recession quarters have elevated exit rates). The
comparison isolates whether the level rescaling affects the WLS slope or
weighted correlation for recovery speed -- the outcomes where sensitivity
is most relevant.

Layout: 2x2 grid
  Row 0: U/LF recovery speed -- full-sample X (left) vs recession-quarter X (right)
  Row 1: V/LF recovery speed -- full-sample X (left) vs recession-quarter X (right)

Output:
  data/results/state_scatter_sensitivity.png
  data/results/state_scatter_sensitivity_latex.tex
"""

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

DELTA_CSV = INSTR_DIR / "delta_instrument_base2006.csv"
LAUS_PATH = INSTR_DIR / "laus_quarterly.parquet"

RECESSIONS = {
    "2001":    {"start": "2001Q1", "end": "2001Q4"},
    "2008-09": {"start": "2007Q4", "end": "2009Q2"},
    "2020":    {"start": "2020Q1", "end": "2020Q2"},
}

PRE_REC_WINDOW     = 4
RECOVERY_THRESHOLD = 0.50
MAX_SEARCH_QTRS    = 28

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
# Recovery speed (recession-end anchor)
# ---------------------------------------------------------------------------
def recovery_from_end(series: pd.Series, rec_end_ql: str,
                      trough_val: float, pre_level: float,
                      direction: str) -> float:
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

# ---------------------------------------------------------------------------
# Load and compute
# ---------------------------------------------------------------------------
def load_and_compute():
    delta = pd.read_csv(DELTA_CSV)
    delta["state_fips"] = delta["state_fips"].astype(str).str.zfill(2)

    laus = pd.read_parquet(LAUS_PATH)
    laus["state_fips"] = laus["state_fips"].astype(str).str.zfill(2)

    # Recession quarters for restricted average
    rec_qtrs = set()
    for info in RECESSIONS.values():
        rec_qtrs.update(ql_range(info["start"], info["end"]))

    avg_full = (delta.groupby("state_fips")["bartik_delta"]
                .mean().rename("avg_delta_full"))
    avg_rec  = (delta[delta["quarter_label"].isin(rec_qtrs)]
                .groupby("state_fips")["bartik_delta"]
                .mean().rename("avg_delta_rec"))

    states = sorted(laus["state_fips"].unique())
    records = []

    for fips in states:
        sl = laus[laus["state_fips"] == fips]
        state_name = sl["state"].iloc[0] if "state" in sl.columns else fips
        unemp_s = sl.set_index("quarter_label")["unemp_rate"]
        vac_s   = sl.set_index("quarter_label")["vac_rate"]
        lf_s    = sl.set_index("quarter_label")["labor_force"]

        rec_rows = []
        for rec_name, info in RECESSIONS.items():
            pre_qtrs = [ql_shift(info["start"], -k)
                        for k in range(1, PRE_REC_WINDOW + 1)]
            pre_u = unemp_s.reindex(pre_qtrs).mean()
            pre_v = vac_s.reindex(pre_qtrs).mean()
            lf_at_start = lf_s.get(info["start"], np.nan)

            search_u = ql_range(info["start"], ql_shift(info["end"], 4))
            u_vals = unemp_s.reindex(search_u).dropna()
            u_trough_val = u_vals.max() if len(u_vals) else np.nan

            search_v = ql_range(info["start"], ql_shift(info["end"], 4))
            v_vals = vac_s.reindex(search_v).dropna()
            v_trough_val = v_vals.min() if len(v_vals) else np.nan

            u_rec = recovery_from_end(
                unemp_s, info["end"], u_trough_val, pre_u, "down")
            v_rec = recovery_from_end(
                vac_s, info["end"], v_trough_val, pre_v, "up")

            rec_rows.append({
                "lf": lf_at_start,
                "u_rec_end": u_rec,
                "v_rec_end": v_rec,
            })

        rec_df = pd.DataFrame(rec_rows).dropna(subset=["lf"])
        if len(rec_df) == 0:
            continue

        w = rec_df["lf"].values / rec_df["lf"].values.sum()
        records.append({
            "state_fips":    fips,
            "state":         state_name,
            "avg_lf":        rec_df["lf"].mean(),
            "avg_delta_full": avg_full.get(fips, np.nan),
            "avg_delta_rec":  avg_rec.get(fips, np.nan),
            "u_rec_end":     np.average(rec_df["u_rec_end"].values, weights=w),
            "v_rec_end":     np.average(
                rec_df["v_rec_end"].fillna(MAX_SEARCH_QTRS + 1).values, weights=w),
        })

    return pd.DataFrame(records)

# ---------------------------------------------------------------------------
# Weighted correlation and WLS slope
# ---------------------------------------------------------------------------
def wls_fit(x, y, w):
    """Returns (intercept, slope) from WLS."""
    W  = np.diag(w)
    X  = np.column_stack([np.ones(len(x)), x])
    coeffs, *_ = np.linalg.lstsq(W @ X, W @ y, rcond=None)
    return coeffs[0], coeffs[1]   # intercept, slope

def wcorr(x, y, w):
    xd = x - np.average(x, weights=w)
    yd = y - np.average(y, weights=w)
    return np.sum(w * xd * yd) / np.sqrt(np.sum(w * xd**2) * np.sum(w * yd**2))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot(df: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    recov_pct = int(RECOVERY_THRESHOLD * 100)

    # Bubble sizes normalized to state average LF
    lf_vals = df["avg_lf"].values
    s_min, s_max = 30, 350
    def bubble_size(lf):
        return s_min + (s_max - s_min) * (lf - lf_vals.min()) / (lf_vals.max() - lf_vals.min())

    # Panel definitions: (ax, xcol, ycol, xlabel_short, ylabel, title)
    xlabel_full = ("Full-sample avg. exit rate\n"
                   r"$\bar{g}^\delta_s$, 1992Q3--2021Q4")
    xlabel_rec  = ("Recession-quarter avg. exit rate\n"
                   r"$\bar{g}^\delta_s|_{\mathrm{rec}}$, 13 recession quarters")

    panels = [
        (axes[0, 0], "avg_delta_full", "u_rec_end", xlabel_full,
         f"Qtrs from rec. end to {recov_pct}% U/LF recovery",
         "Unemployment recovery -- full-sample X"),
        (axes[0, 1], "avg_delta_rec",  "u_rec_end", xlabel_rec,
         f"Qtrs from rec. end to {recov_pct}% U/LF recovery",
         "Unemployment recovery -- recession-quarter X"),
        (axes[1, 0], "avg_delta_full", "v_rec_end", xlabel_full,
         f"Qtrs from rec. end to {recov_pct}% V/LF recovery",
         "Vacancy recovery -- full-sample X"),
        (axes[1, 1], "avg_delta_rec",  "v_rec_end", xlabel_rec,
         f"Qtrs from rec. end to {recov_pct}% V/LF recovery",
         "Vacancy recovery -- recession-quarter X"),
    ]

    for ax, xcol, ycol, xlabel, ylabel, title in panels:
        sub = df.dropna(subset=[xcol, ycol]).copy()
        x  = sub[xcol].values
        y  = sub[ycol].values
        w  = sub["avg_lf"].values
        sz = bubble_size(w)

        ax.scatter(x, y, s=sz, alpha=0.6, color="steelblue",
                   edgecolors="white", linewidths=0.4, zorder=4)

        # WLS line
        b, m = wls_fit(x, y, w)
        x_line = np.linspace(x.min() * 0.88, x.max() * 1.12, 200)
        ax.plot(x_line, b + m * x_line, color="firebrick", linewidth=1.5,
                linestyle="--", alpha=0.85)

        # Annotations
        r_w = wcorr(x, y, w)
        ax.text(0.05, 0.91, f"$r_w$ = {r_w:.2f}", transform=ax.transAxes,
                fontsize=9, color="dimgrey")
        ax.text(0.05, 0.83, f"slope = {m:.1f}", transform=ax.transAxes,
                fontsize=9, color="dimgrey")

        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.grid(linewidth=0.4, alpha=0.4)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{v:.3f}"))

    # Column labels
    axes[0, 0].set_title("Unemployment recovery\nFull-sample X",
                          fontsize=10, fontweight="bold")
    axes[0, 1].set_title("Unemployment recovery\nRecession-quarter X",
                          fontsize=10, fontweight="bold")
    axes[1, 0].set_title("Vacancy recovery\nFull-sample X",
                          fontsize=10, fontweight="bold")
    axes[1, 1].set_title("Vacancy recovery\nRecession-quarter X",
                          fontsize=10, fontweight="bold")

    # Bubble legend
    lf_ref = [lf_vals.min(), np.median(lf_vals), lf_vals.max()]
    handles = [
        plt.scatter([], [], s=bubble_size(lf), color="steelblue", alpha=0.6,
                    edgecolors="white", linewidths=0.4, label=lbl)
        for lf, lbl in zip(lf_ref, ["Small state", "Median state", "Large state"])
    ]
    fig.legend(handles=handles, title="Avg. labor force (bubble size)",
               loc="lower center", ncol=3, fontsize=8, title_fontsize=8,
               bbox_to_anchor=(0.5, -0.03), frameon=True)

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    out = RESULTS_DIR / "state_scatter_sensitivity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")
    return out

# ---------------------------------------------------------------------------
# LaTeX snippet
# ---------------------------------------------------------------------------
LATEX_SNIPPET = r"""\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\textwidth]{state_scatter_sensitivity.png}
    \caption{Sensitivity of cross-state recovery speed to the definition
    of the establishment exit rate. Each row shows the relationship between
    state-level exit rates and recession-end-anchored 50\% recovery times
    for unemployment (top) and vacancies (bottom). The left column uses
    the full-sample average exit rate ($\bar{g}^\delta_s$, 1992Q3--2021Q4)
    as in Figure~\ref{fig:state_scatter}; the right column restricts the
    average to the 13 recession quarters spanning the three NBER episodes
    (2001, 2007--2009, 2020). The two X-axis measures have a cross-state
    correlation of $r = 0.995$ and differ only by a uniform scaling factor
    of approximately 1.18. Regression lines are estimated by WLS with state
    labor force weights; $r_w$ denotes the labor-force-weighted correlation
    and slope is in quarters per unit exit rate.}
    \label{fig:state_scatter_sensitivity}
\end{figure}"""

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading and computing ...")
    df = load_and_compute()
    print(f"  States: {len(df)}")

    print("\nSensitivity comparison:")
    print(f"  Cross-state corr(avg_delta_full, avg_delta_rec) = "
          f"{df['avg_delta_full'].corr(df['avg_delta_rec']):.3f}")
    print(f"  Mean ratio rec/full = "
          f"{(df['avg_delta_rec'] / df['avg_delta_full']).mean():.3f}")

    for ycol, label in [("u_rec_end", "U recovery"), ("v_rec_end", "V recovery")]:
        sub = df.dropna(subset=["avg_delta_full", "avg_delta_rec", ycol])
        w   = sub["avg_lf"].values
        r_full = wcorr(sub["avg_delta_full"].values, sub[ycol].values, w)
        r_rec  = wcorr(sub["avg_delta_rec"].values,  sub[ycol].values, w)
        _, m_full = wls_fit(sub["avg_delta_full"].values, sub[ycol].values, w)
        _, m_rec  = wls_fit(sub["avg_delta_rec"].values,  sub[ycol].values, w)
        print(f"\n  {label}:")
        print(f"    Full-sample X:        r_w={r_full:.3f}  slope={m_full:.1f}")
        print(f"    Recession-quarter X:  r_w={r_rec:.3f}  slope={m_rec:.1f}")

    print("\nPlotting ...")
    plot(df)

    latex_out = RESULTS_DIR / "state_scatter_sensitivity_latex.tex"
    with open(latex_out, "w") as f:
        f.write(LATEX_SNIPPET)
        f.write("\n")
    print(f"LaTeX snippet saved: {latex_out}")
    print("\n--- LaTeX snippet ---")
    print(LATEX_SNIPPET)
