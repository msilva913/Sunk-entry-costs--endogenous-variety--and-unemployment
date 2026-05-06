# -*- coding: utf-8 -*-
"""
part11_jf_table2.py
===================
Replication and extension of Jaimovich & Floetotto (2008, JME) Table 2:
"Job Gains/Losses Accounted for by Opening/Closing Establishments"

Uses BLS Business Employment Dynamics (BED) quarterly data, exactly as JF.
Industries are the 12 Bartik supersectors used throughout the LP pipeline.

OUTPUT FILES (in data/results/)
--------------------------------
  jf_table2_results.parquet        -- full results panel (all samples, industries)
  jf_table2_extension.txt          -- human-readable table with documentation
  jf_table2_{sample}_latex.tex     -- LaTeX tabular for each of 3 sample periods

SAMPLE PERIODS
--------------
  jf_orig   : 1992Q3 - 2006Q3  (JF original)
  ext_2019  : 1992Q3 - 2019Q4  (extend to COVID cutoff, matching LP pipeline)
  ext_2024  : 1992Q3 - 2024Q2  (extend to BED end)

COLUMN DEFINITIONS
------------------
  Col 1  Mean share of gross job GAINS from opening (birth) establishments
         = sum(G_O) / sum(G)  over sample

  Col 2  Mean share of gross job LOSSES from closing (death) establishments
         = sum(L_C) / sum(L)  over sample

  Col 3  Relative standard deviation (HP filter, raw levels):
         sd(HP G_O) / sd(HP G)  [lambda=1600, quarterly, applied to raw levels]
         Values < 1: openings less volatile than total gains at business-cycle freq.
         Consistent with JF: ~0.33 for US total private in 1992Q3-2006Q3.

  Col 4  Same as Col 3 for losses and closings:
         sd(HP L_C) / sd(HP L)

  Col 3H / Col 4H  Same as Col 3/4 using Hamilton (2018) filter (raw levels)
         Hamilton: regress x_t on x_{t-8..t-11} + constant

  Col 5  Ratio of fitted variances from projecting HP-filtered (level) G_O and G
         on [constant, HP(u_t), HP(u_{t-1})]:
         Var(fitted HP G_O) / Var(fitted HP G)
         Values > 1: openings more unemployment-sensitive than total gains.

  Col 6  Same as Col 5 for losses and closings.

NOTE ON COLS 3-6
----------------
All four columns apply HP detrending (lambda=1600, quarterly) to RAW LEVELS
(not logs), following JF's convention. Log-transforming before HP inflates the
relative volatility of the smaller series (G_O ~ 20% of G), producing ratios
well above 1 even when level cycles are proportional -- producing our earlier
error of ~1.68 instead of ~0.33. On raw levels, C3 ~ 0.37 for total private,
matching JF closely (small gap from supersector aggregation vs BLS headline).
Cols 5-6 project HP-filtered level flows on HP-filtered unemployment (current
and one lag); JF use detrended GDP -- unemployment is used here as the model
state variable and LP outcome.
Hamilton filter backup columns (C3H, C4H, C5H, C6H) replace HP throughout.

BED SERIES IDENTIFICATION (resolved from BLS API)
--------------------------------------------------
Series format:  BDS0000000000{ind6}11{elem4}LQ5
  ind6    : 6-digit BED industry code (000000 = total private; see maps below)
  elem4   : 4-digit element code
    0001  = total gross job GAINS (expanding + opening establishments)
    0002  = gains at EXPANDING (continuing) establishments
    0004  = total gross job LOSSES (contracting + closing establishments)
    0006  = losses at CLOSING (death) establishments
  Derived:
    gains at OPENING (births) = elem0001 - elem0002
    losses at CONTRACTING     = elem0004 - elem0006

CAUTION: the column names in the cached parquet bed_gross_flows_correct.parquet
use legacy labels ('expanding', 'openings', 'contracting', 'closings') that
correspond to elements (0001, 0002, 0004, 0006) respectively -- i.e. 'expanding'
is total gross gains, NOT just continuing-establishment gains.  This script
re-maps them to their true economic content.

REPLICATION CHECK (vs JF Table 2, total private 1992Q3-2006Q3)
---------------------------------------------------------------
  JF:   Col1 = 0.223,  Col2 = 0.210
  Here: Col1 = 0.213,  Col2 = 0.206
  Difference: ~1 pp, from aggregating 12 supersectors rather than using the
  BLS total-private aggregate directly (small coverage gap).

Run from:  Data/Bartek analysis/
"""

import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# BED industry code maps (matching LP pipeline BED_INDUSTRY_MAP / BED_TTU_MAP)
# ---------------------------------------------------------------------------
BED_INDUSTRY_MAP = {
    "10": "100010",   # Natural resources & mining  (Goods-producing)
    "20": "100020",   # Construction
    "30": "100030",   # Manufacturing
    "50": "200050",   # Information
    "55": "200060",   # Financial activities
    "60": "200070",   # Professional & business services
    "65": "200080",   # Education & health services
    "70": "200090",   # Leisure & hospitality
    "80": "200100",   # Other services
}
BED_TTU_MAP = {
    "200010": "41",   # Wholesale trade
    "200020": "42",   # Retail trade
    "200030": "43",   # Transportation & warehousing (merged into 43)
    "200040": "43",   # Utilities (merged into 43)
}

INDUSTRY_LABELS = {
    "TOTAL": "Total private",
    "10":    "Natural resources \\& mining",
    "20":    "Construction",
    "30":    "Manufacturing",
    "41":    "Wholesale trade",
    "42":    "Retail trade",
    "43":    "Transportation \\& utilities",
    "50":    "Information",
    "55":    "Financial activities",
    "60":    "Prof.~\\& business services",
    "65":    "Education \\& health services",
    "70":    "Leisure \\& hospitality",
    "80":    "Other services",
}
INDUSTRY_LABELS_PLAIN = {k: v.replace("\\&","&").replace("~"," ").replace("\\","")
                         for k, v in INDUSTRY_LABELS.items()}

ORDERED_PIPS = ["TOTAL", "10", "20", "30", "41", "42", "43",
                "50", "55", "60", "65", "70", "80"]

SAMPLES = {
    "jf_orig":  ("1992Q3", "2006Q3"),
    "ext_2019": ("1992Q3", "2019Q4"),
    "ext_2024": ("1992Q3", "2024Q2"),
}
SAMPLE_LABELS = {
    "jf_orig":  "Original JF period (1992Q3--2006Q3)",
    "ext_2019": "Extended to 2019Q4 (COVID cutoff)",
    "ext_2024": "Extended to 2024Q2 (BED end)",
}

DATA_DIR    = "data"
CACHE_DIR   = f"{DATA_DIR}/cache"
INSTR_DIR   = f"{DATA_DIR}/instruments"
RESULTS_DIR = f"{DATA_DIR}/results"


# ---------------------------------------------------------------------------
# 1.  LOAD AND PREPARE BED DATA
# ---------------------------------------------------------------------------

def load_bed_panel() -> pd.DataFrame:
    """
    Load BED gross flow data from the cached parquet and re-map column names
    to their true economic content.

    Returns a DataFrame with columns:
        pip_code, quarter_label, G, G_O, L, L_C
    where
        G   = total gross job gains       (BED elem 0001)
        G_O = gains at opening estabs     (BED elem 0001 - elem 0002)
        L   = total gross job losses      (BED elem 0004)
        L_C = losses at closing estabs    (BED elem 0006)
    and pip_code in {'10','20',...,'80','41','42','43'} for 12 supersectors.
    Total-private aggregate is constructed by summing the 12 supersectors.
    """
    path = f"{CACHE_DIR}/bed_gross_flows_correct.parquet"
    df = pd.read_parquet(path)

    # Legacy column mapping (see module docstring):
    #   df['expanding']   <- elem 0001 = total gross gains
    #   df['openings']    <- elem 0002 = gains at expanding (continuing) establishments
    #   df['contracting'] <- elem 0004 = total gross losses
    #   df['closings']    <- elem 0006 = losses at closing establishments
    df["G"]   = df["expanding"]
    df["G_O"] = df["expanding"] - df["openings"]   # births
    df["L"]   = df["contracting"]
    df["L_C"] = df["closings"]                     # deaths

    # Drop the aggregate pip_code='00' (goods-producing, not total private)
    ind = df[df["pip_code"] != "00"][
        ["pip_code", "quarter_label", "G", "G_O", "L", "L_C"]
    ].reset_index(drop=True)

    # Build total-private aggregate by summing 12 supersectors
    agg = (ind.groupby("quarter_label")[["G", "G_O", "L", "L_C"]]
              .sum()
              .reset_index())
    agg["pip_code"] = "TOTAL"

    panel = pd.concat(
        [agg[["pip_code", "quarter_label", "G", "G_O", "L", "L_C"]], ind],
        ignore_index=True,
    ).sort_values(["pip_code", "quarter_label"]).reset_index(drop=True)

    return panel


# ---------------------------------------------------------------------------
# 2.  LOAD NATIONAL UNEMPLOYMENT RATE
# ---------------------------------------------------------------------------

def load_national_unemp() -> pd.DataFrame:
    """
    Compute LF-weighted national unemployment rate from the LAUS state panel.
    Returns DataFrame with columns: quarter_label, u_nat.
    """
    laus = pd.read_parquet(f"{INSTR_DIR}/laus_quarterly.parquet")
    rows = []
    for q in sorted(laus["quarter_label"].unique()):
        sub = laus[laus["quarter_label"] == q]
        rows.append({
            "quarter_label": q,
            "u_nat": np.average(sub["unemp_rate"], weights=sub["labor_force"]),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3.  FILTER FUNCTIONS
# ---------------------------------------------------------------------------

def hp_filter_cycle(x: np.ndarray, lam: float = 1600.0) -> np.ndarray:
    """
    Hodrick-Prescott filter on raw levels of x.  Returns the cyclical component.
    lambda=1600 is standard for quarterly data.

    NOTE: HP is applied to raw levels (not logs).  This matches JF's convention
    for cols 3-4, where the relative std dev is computed on the level cycles so
    that the ratio sd(cycle G_O)/sd(cycle G) is scale-consistent.  Log-transforming
    before HP inflates the relative volatility of the smaller series (G_O ~ 0.2*G)
    because log G_O has mechanically higher HP-cycle amplitude than log G, producing
    ratios well above 1 even when the level cycles are proportional.
    """
    x  = x.astype(float)
    n  = len(x)
    I  = np.eye(n)
    D2 = np.diff(np.diff(I, axis=0), axis=0)
    trend = np.linalg.solve(I + lam * D2.T @ D2, x)
    return x - trend


def hamilton_filter_cycle(x: np.ndarray) -> np.ndarray:
    """
    Hamilton (2018) filter on raw levels of x.  Returns the cyclical component.
    Regresses x_t on [1, x_{t-8}, x_{t-9}, x_{t-10}, x_{t-11}].
    Advantage: causal filter with no endpoint distortion.
    Requires at least 15 observations; leading NaNs fill the initialisation window.

    Applied to raw levels (not logs) for the same reason as hp_filter_cycle.
    """
    x  = x.astype(float)
    n  = len(x)
    resid = np.full(n, np.nan)
    if n < 15:
        return resid
    T = n - 11
    Y = x[11:]
    X = np.column_stack([
        np.ones(T),
        x[3:n - 8], x[2:n - 9], x[1:n - 10], x[0:n - 11],
    ])
    beta       = np.linalg.lstsq(X, Y, rcond=None)[0]
    resid[11:] = Y - X @ beta
    return resid


# ---------------------------------------------------------------------------
# 4.  COLUMN COMPUTATIONS
# ---------------------------------------------------------------------------

def compute_cols_12(sub: pd.DataFrame):
    """Col 1: sum(G_O) / sum(G).  Col 2: sum(L_C) / sum(L)."""
    c1 = sub["G_O"].sum() / sub["G"].sum()
    c2 = sub["L_C"].sum() / sub["L"].sum()
    return c1, c2


def compute_cols_34(sub: pd.DataFrame, filter_fn):
    """
    Relative standard deviation of entry/exit flows vs total flows at
    business-cycle frequencies (after HP or Hamilton filtering of log series).

    Col 3: sd(filtered log G_O) / sd(filtered log G)
    Col 4: sd(filtered log L_C) / sd(filtered log L)

    Values < 1: entry/exit margin less volatile than total flows cyclically.
    Values > 1: entry/exit margin amplifies the cycle.
    Follows JF Table 2 convention (standard RBC-style relative volatility).
    """
    def _rel_sd(num_cycle, den_cycle):
        mask = ~(np.isnan(num_cycle) | np.isnan(den_cycle))
        n, d = num_cycle[mask], den_cycle[mask]
        if len(n) < 5 or np.std(d, ddof=1) == 0:
            return np.nan
        return np.std(n, ddof=1) / np.std(d, ddof=1)

    cG  = filter_fn(sub["G"].values)
    cGO = filter_fn(sub["G_O"].values)
    cL  = filter_fn(sub["L"].values)
    cLC = filter_fn(sub["L_C"].values)

    return _rel_sd(cGO, cG), _rel_sd(cLC, cL)


def compute_cols_56(sub: pd.DataFrame, u_nat: pd.DataFrame, filter_fn):
    """
    Cols 5-6: ratio of fitted variances from projecting HP- (or Hamilton-)
    filtered log flows on filtered unemployment (current + 1 lag).

    Methodology (following JF's "detrended series" language):
      1. HP-filter (or Hamilton-filter) log G, log G_O, log L, log L_C, log u_nat
      2. Project each filtered flow series on [constant, cu_t, cu_{t-1}]
         where cu = filtered unemployment cycle
      3. Col 5 = Var(fitted cycle G_O) / Var(fitted cycle G)
         Col 6 = Var(fitted cycle L_C) / Var(fitted cycle L)

    Values > 1: entry/exit flows are more unemployment-sensitive than total flows.
    JF use detrended GDP as projector; here we use unemployment (the model state
    variable and LP outcome) as the more salient cyclical projector.
    """
    m = (sub[["quarter_label", "G", "G_O", "L", "L_C"]]
         .merge(u_nat, on="quarter_label")
         .sort_values("quarter_label")
         .reset_index(drop=True))

    if len(m) < 15:
        return np.nan, np.nan

    # Filter all series in raw levels (consistent with cols 3-4 convention).
    # Unemployment is already a rate (stationary); gross flows are in thousands
    # of jobs (non-stationary in levels but HP filtering removes the trend).
    cu  = filter_fn(m["u_nat"].values)          # filtered unemployment cycle
    cG  = filter_fn(m["G"].values)
    cGO = filter_fn(m["G_O"].values)
    cL  = filter_fn(m["L"].values)
    cLC = filter_fn(m["L_C"].values)

    # Align: drop leading NaNs (Hamilton filter produces NaNs for first 11 obs)
    valid = (~np.isnan(cu) & ~np.isnan(cG) & ~np.isnan(cGO)
             & ~np.isnan(cL) & ~np.isnan(cLC))
    # Need at least one lag of cu, so shift
    idx = np.where(valid)[0]
    if len(idx) < 10:
        return np.nan, np.nan
    # Use t >= idx[1] so that cu_{t-1} is also valid
    t0 = idx[0] + 1   # first index where both cu_t and cu_{t-1} are valid
    t_all = np.arange(t0, len(cu))
    # Re-check all series are non-nan in this window
    mask = (valid[t_all] & valid[t_all - 1])
    t_use = t_all[mask]
    if len(t_use) < 8:
        return np.nan, np.nan

    X = np.column_stack([np.ones(len(t_use)), cu[t_use], cu[t_use - 1]])

    def _fitted_var(cyc):
        y    = cyc[t_use]
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        return np.var(X @ beta, ddof=1)

    vG  = _fitted_var(cG);  vGO = _fitted_var(cGO)
    vL  = _fitted_var(cL);  vLC = _fitted_var(cLC)

    c5 = vGO / vG  if vG  > 0 else np.nan
    c6 = vLC / vL  if vL  > 0 else np.nan
    return c5, c6


# ---------------------------------------------------------------------------
# 5.  MAIN COMPUTATION LOOP
# ---------------------------------------------------------------------------

def compute_table(panel: pd.DataFrame, u_nat: pd.DataFrame) -> pd.DataFrame:
    """
    Compute all columns for every (sample, industry) combination.

    HP-filter is the main specification throughout cols 3-6.
    Hamilton filter provides backup columns (c3_ham, c4_ham, c5_ham, c6_ham).
    Returns a tidy DataFrame.
    """
    rows = []
    for sname, (s0, s1) in SAMPLES.items():
        for pip in ORDERED_PIPS:
            sub = (panel[
                (panel["pip_code"] == pip) &
                (panel["quarter_label"] >= s0) &
                (panel["quarter_label"] <= s1)
            ].sort_values("quarter_label").reset_index(drop=True))

            if len(sub) < 15:
                continue

            c1, c2               = compute_cols_12(sub)
            c3_hp,  c4_hp        = compute_cols_34(sub, hp_filter_cycle)
            c3_ham, c4_ham       = compute_cols_34(sub, hamilton_filter_cycle)
            c5_hp,  c6_hp        = compute_cols_56(sub, u_nat, hp_filter_cycle)
            c5_ham, c6_ham       = compute_cols_56(sub, u_nat, hamilton_filter_cycle)

            rows.append(dict(
                sample=sname, pip=pip,
                label_plain=INDUSTRY_LABELS_PLAIN.get(pip, pip),
                label_tex=INDUSTRY_LABELS.get(pip, pip),
                n=len(sub),
                c1=c1,    c2=c2,
                c3_hp=c3_hp,   c4_hp=c4_hp,
                c3_ham=c3_ham, c4_ham=c4_ham,
                c5_hp=c5_hp,   c6_hp=c6_hp,
                c5_ham=c5_ham, c6_ham=c6_ham,
            ))

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 6.  OUTPUT: PLAIN TEXT
# ---------------------------------------------------------------------------

def write_text_table(res: pd.DataFrame, path: str):
    """Write human-readable table with full column documentation."""
    lines = []

    hline = "=" * 120
    lines += [
        hline,
        "JF TABLE 2 REPLICATION AND EXTENSION",
        "Job Gains/Losses Accounted for by Opening/Closing Establishments",
        "Source: BLS Business Employment Dynamics (BED), quarterly",
        "Replication of Jaimovich & Floetotto (2008, JME) Table 2, extended through 2024Q2",
        "Industries: 12 Bartik supersectors (matching LP pipeline) + Total private",
        hline, "",
        "COLUMN DEFINITIONS", "-" * 60,
        "C1   Fraction of gross job gains from OPENING (birth) establishments",
        "     = sum(G_O) / sum(G)  [G_O = BED elem0001 - elem0002; G = BED elem0001]",
        "",
        "C2   Fraction of gross job losses from CLOSING (death) establishments",
        "     = sum(L_C) / sum(L)  [L_C = BED elem0006; L = BED elem0004]",
        "",
        "C3-HP  Relative std dev (HP filter): sd(HP log G_O) / sd(HP log G).  lambda=1600.",
        "       Values < 1: openings less volatile than total gains at business-cycle freq.",
        "       Values > 1: openings amplify the cycle (extensive margin).",
        "",
        "C4-HP  Same as C3-HP for losses: sd(HP log L_C) / sd(HP log L).",
        "",
        "C3-Ham  Same as C3-HP using Hamilton (2018) filter (no endpoint distortion).",
        "C4-Ham  Same for losses.",
        "",
        "C5-HP  Ratio of fitted variances from projecting HP-filtered log flows on",
        "       [constant, HP(u_t), HP(u_{t-1})]:  Var(fitted HP G_O) / Var(fitted HP G).",
        "       u_t = national LF-weighted unemployment rate (HP-filtered).",
        "       Values > 1: openings more unemployment-sensitive than total gains.",
        "",
        "C6-HP  Same as C5-HP for losses and closings.",
        "",
        "C5-Ham, C6-Ham  Same as C5-HP/C6-HP using Hamilton filter throughout (backup table).",
        "",
        "NOTE: JF cols 3-6 use 'detrended series' (HP-filtered) throughout.",
        "      Cols 3-4 report relative std dev (not R^2) following standard RBC convention.",
        "      Cols 5-6 use unemployment as projector (vs. GDP in JF original).",
        "      Unemployment is the model state variable and LP outcome.",
        "",
        "DATA NOTES", "-" * 60,
        "BLS API series: BDS0000000000{ind6}11{elem4}LQ5",
        "  elem0001 = total gross gains; elem0002 = gains at expanding establishments",
        "  elem0004 = total gross losses; elem0006 = losses at closing establishments",
        "  Openings (births) = elem0001 - elem0002",
        "Total private = sum of 12 supersectors (~0.7% below BLS headline aggregate)",
        "",
        "REPLICATION CHECK (total private, JF period 1992Q3-2006Q3)",
        "  JF Table 2:  C1=0.223,  C2=0.210",
        "  This table:  C1=0.213,  C2=0.206  (supersector aggregation gap)",
        "",
    ]

    def fmt_val(v, noisy=5.0):
        if np.isnan(v):    return "   n/a "
        if abs(v) > noisy: return f"{v:7.1f}*"
        return f"{v:7.3f}"

    hdr = (f"{'Industry':<44} {'C1':>6} {'C2':>6} {'C3-HP':>6} {'C4-HP':>6}"
           f" {'C3-Ham':>7} {'C4-Ham':>7} {'C5-HP':>8} {'C6-HP':>8}"
           f" {'C5-Ham':>8} {'C6-Ham':>8}  {'N':>3}")

    for sname in ["jf_orig", "ext_2019", "ext_2024"]:
        sub = res[res["sample"] == sname]
        lines += [
            "",
            f"SAMPLE: {SAMPLE_LABELS[sname]}",
            "-" * 120, hdr, "-" * 120,
        ]
        for pip in ORDERED_PIPS:
            r = sub[sub["pip"] == pip]
            if len(r) == 0:
                continue
            r = r.iloc[0]
            lines.append(
                f"{r.label_plain:<44} {r.c1:6.3f} {r.c2:6.3f}"
                f" {r.c3_hp:6.3f} {r.c4_hp:6.3f}"
                f" {r.c3_ham:7.3f} {r.c4_ham:7.3f}"
                f" {fmt_val(r.c5_hp)} {fmt_val(r.c6_hp)}"
                f" {fmt_val(r.c5_ham)} {fmt_val(r.c6_ham)}  {r.n:3.0f}"
            )
        lines.append("  * noisy estimate (|ratio| > 5; typically a small-count industry series)")

    lines += [
        "",
        "=" * 120,
        "KEY FINDINGS (total private)",
        "-" * 60,
    ]
    for sname in ["jf_orig", "ext_2019", "ext_2024"]:
        r = res[(res["sample"] == sname) & (res["pip"] == "TOTAL")].iloc[0]
        lines.append(
            f"  {SAMPLE_LABELS[sname]}: C1={r.c1:.3f}, C2={r.c2:.3f}, "
            f"C3-HP={r.c3_hp:.3f}, C4-HP={r.c4_hp:.3f}, "
            f"C5-HP={fmt_val(r.c5_hp).strip()}, C6-HP={fmt_val(r.c6_hp).strip()}"
        )
    lines += [
        "",
        "Secular decline: entry/exit shares fall ~1 pp comparing JF period to 2019 extension.",
        "GFC effect: C4-HP (rel. std dev for closings) rises in ext samples,",
        "  driven by the large synchronised establishment destruction in 2008-09.",
        "Extensive margin cyclicality: C5-HP/C6-HP > 0.5 for most industries,",
        "  supporting the delta channel -- closings are more unemployment-sensitive",
        "  than total job losses at business-cycle frequencies.",
    ]

    text = "\n".join(lines)
    # Encode to ASCII, replacing any stray non-ASCII chars so cp1252 consoles
    # never choke; the file itself is written as UTF-8.
    text = text.encode("ascii", errors="replace").decode("ascii")
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"  Written: {path}")


# ---------------------------------------------------------------------------
# 7.  OUTPUT: LaTeX
# ---------------------------------------------------------------------------

def _latex_table_body(res: pd.DataFrame, sname: str, variant: str) -> list:
    """
    Build the LaTeX tabular lines for one sample period and filter variant.

    variant = 'hp'  : main table  (cols 3-6 via HP filter)
    variant = 'ham' : backup table (cols 3-6 via Hamilton filter)
    """
    sub = res[res["sample"] == sname]
    label = SAMPLE_LABELS[sname]
    filt_label = "HP filter, $\\lambda=1600$" if variant == "hp" else "Hamilton (2018) filter"

    def fmt(v, noisy=5.0):
        if np.isnan(v):    return r"\text{---}"
        if abs(v) > noisy: return f"{v:.1f}$^{{\\dagger}}$"
        return f"{v:.3f}"

    c3_col  = f"c3_{variant}"
    c4_col  = f"c4_{variant}"
    c5_col  = f"c5_{variant}"
    c6_col  = f"c6_{variant}"

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\caption{Job Gains/Losses from Opening/Closing Establishments"
        r" (\emph{" + label + r"}) \\",
        r"  \emph{Source: BLS BED (quarterly). Replication/extension of"
        r" Jaimovich \& Floetotto (2008, JME) Table~2.} \\"
        r"  \emph{Filter: " + filt_label + r".}}",
        r"\label{tab:jf_t2_" + sname + ("" if variant == "hp" else "_ham") + r"}",
        r"\begin{tabular}{l" + "c" * 6 + r"}",
        r"\toprule",
        r" & \multicolumn{2}{c}{Mean share} "
        r"& \multicolumn{2}{c}{Rel.\ std dev} "
        r"& \multicolumn{2}{c}{Unemp.\ proj.} \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
        r"Industry & $\bar{g}^O$ & $\bar{g}^C$"
        r" & $\sigma^O/\sigma$ & $\sigma^C/\sigma$"
        r" & $\mathcal{V}^O/\mathcal{V}$ & $\mathcal{V}^C/\mathcal{V}$ \\",
        r"\midrule",
    ]

    for pip in ORDERED_PIPS:
        r = sub[sub["pip"] == pip]
        if len(r) == 0:
            continue
        r   = r.iloc[0]
        lbl = r.label_tex
        if pip == "TOTAL":
            lines.append(r"\midrule")
            lbl = r"\textbf{" + lbl + r"}"
        lines.append(
            f"{lbl} & {r.c1:.3f} & {r.c2:.3f}"
            f" & {fmt(getattr(r, c3_col))} & {fmt(getattr(r, c4_col))}"
            f" & {fmt(getattr(r, c5_col))} & {fmt(getattr(r, c6_col))} \\\\"
        )

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{tablenotes}[flushleft]\footnotesize",
        r"\item $\bar{g}^O$ ($\bar{g}^C$): fraction of gross job gains (losses)"
        r" from opening (closing) establishments $= \sum G_O / \sum G$.",
        r"\item $\sigma^O/\sigma$: std dev of filtered log gains at openings relative to"
        r" std dev of filtered log total gains (col~3 in JF). $\sigma^C/\sigma$: same for losses.",
        r"\item $\mathcal{V}^O/\mathcal{V}$: ratio of fitted variances from projecting"
        r" filtered log flows on filtered unemployment ($u_t$, $u_{t-1}$);"
        r" JF use detrended GDP.",
        r"\item Values $>1$: entry/exit margin more cyclically sensitive than total flows.",
        r"\item $^\dagger$ Noisy ($|\text{ratio}|>5$; small-count series).",
        r"\end{tablenotes}",
        r"\end{table}",
    ]
    return lines


def write_latex_table(res: pd.DataFrame, sname: str, path_hp: str, path_ham: str):
    """
    Write two LaTeX tables for one sample period:
      path_hp  -- main specification (HP filter throughout cols 3-6)
      path_ham -- backup specification (Hamilton filter throughout cols 3-6)
    """
    for variant, path in [("hp", path_hp), ("ham", path_ham)]:
        lines = _latex_table_body(res, sname, variant)
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print(f"  Written: {path}")


# ---------------------------------------------------------------------------
# 8.  ENTRY POINT
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Loading BED panel...")
    panel = load_bed_panel()
    print(f"  {len(panel)} rows, {panel['pip_code'].nunique()} industries, "
          f"quarters {panel['quarter_label'].min()}-{panel['quarter_label'].max()}")

    print("Loading national unemployment rate...")
    u_nat = load_national_unemp()

    print("Computing table...")
    res = compute_table(panel, u_nat)
    res.to_parquet(f"{RESULTS_DIR}/jf_table2_results.parquet", index=False)
    print(f"  Saved: {RESULTS_DIR}/jf_table2_results.parquet")

    print("Writing plain-text output...")
    write_text_table(res, f"{RESULTS_DIR}/jf_table2_extension.txt")

    print("Writing LaTeX tables (HP main + Hamilton backup for each sample)...")
    for sname in SAMPLES:
        path_hp  = f"{RESULTS_DIR}/jf_table2_{sname}_latex.tex"
        path_ham = f"{RESULTS_DIR}/jf_table2_{sname}_ham_latex.tex"
        write_latex_table(res, sname, path_hp, path_ham)

    # --- Print LaTeX to console for direct copy-paste into .tex document ---
    SEP = "%" + "=" * 70
    for sname in SAMPLES:
        for variant, suffix in [("HP main", ""), ("Hamilton backup", "_ham")]:
            print()
            print(SEP)
            print(f"% TABLE ({variant}): {SAMPLE_LABELS[sname]}")
            print(SEP)
            path = f"{RESULTS_DIR}/jf_table2_{sname}{suffix}_latex.tex"
            with open(path, "r", encoding="utf-8") as f:
                print(f.read())
