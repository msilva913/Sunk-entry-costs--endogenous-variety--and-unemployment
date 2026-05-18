import numpy as np
import pandas as pd
from tabulate import tabulate

def create_stats_table(data, caption="Statistical Summary", label="tab:stats"):
    """
    Create a LaTeX table with dynamic column handling.
    """
    columns = data.columns
    latex_str = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\caption{" + caption + "}",
        "\\label{" + label + "}",
        "\\begin{tabular}{l" + "".join(["r"] * len(columns)) + "}",
        "\\toprule"
    ]
    header_row = ["Variable"] + [
        col.replace("_", "\\_").replace("-", "$-$")
           .replace("(", "\\left(").replace(")", "\\right)")
        for col in columns
    ]
    latex_str.append(" & ".join(header_row) + " \\\\")
    latex_str.append("\\midrule")
    for idx, row in data.iterrows():
        formatted_row = [str(idx)] + [
            f"{val:.3f}" if isinstance(val, (int, float)) else str(val)
            for val in row
        ]
        latex_str.append(" & ".join(formatted_row) + " \\\\")
    latex_str.extend([
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}"
    ])
    return "\n".join(latex_str)


def generate_stacked_moments_latex_table(summ):
    """
    Generate LaTeX code for a table with subheadings for standard deviations,
    cross correlations, and autocorrelations using booktabs.
    """
    headers = [f"Value{i+1}" for i in range(summ.shape[1])]
    header_line = " & ".join(["\\textbf{" + h + "}" for h in headers])
    nc = str(len(headers) + 1)
    latex_table = (
        "\n    \\begin{table}[h]\n    \\centering\n"
        "    \\begin{tabular}{l" + "r" * len(headers) + "}\n"
        "    \\toprule\n"
        "    \\textbf{Moment Type} & " + header_line + " \\\\\n"
        "    \\midrule\n"
        "    \\multicolumn{" + nc + "}{c}{\\textbf{Standard Deviations}} \\\\\n"
        "    \\midrule\n"
    )
    stds = summ.loc[summ.index.str.startswith('std')]
    for idx in stds.index:
        values_line = " & ".join([f"{val:.3f}" for val in stds.loc[idx]])
        latex_table += f"    {idx} & {values_line} \\\\\n"

    latex_table += "    \\midrule \\multicolumn{" + nc + "}{c}{\\textbf{Cross Correlations}} \\\\ \\midrule\n"
    cross_corrs = summ.loc[summ.index.str.startswith('Cor(') & ~summ.index.str.contains('_{-1}')]
    for idx in cross_corrs.index:
        values_line = " & ".join([f"{val:.3f}" for val in cross_corrs.loc[idx]])
        latex_table += f"    {idx} & {values_line} \\\\\n"

    latex_table += "    \\midrule \\multicolumn{" + nc + "}{c}{\\textbf{Autocorrelations}} \\\\ \\midrule\n"
    autocorrs = summ.loc[summ.index.str.contains('_{-1}')]
    for idx in autocorrs.index:
        values_line = " & ".join([f"{val:.3f}" for val in autocorrs.loc[idx]])
        latex_table += f"    {idx} & {values_line} \\\\\n"

    latex_table += "    \\bottomrule \\end{tabular} \\caption{Stacked Moments Summary} \\end{table}"
    return latex_table


def generate_stacked_moments_table(summ):
    """
    Generate a plain text table using tabulate for moments, including subheadings.
    """
    headers = ["Moment Type"] + [f"Value{i+1}" for i in range(summ.shape[1])]
    data = []

    data.append(["Standard Deviations"] + [""] * summ.shape[1])
    stds = summ.loc[summ.index.str.startswith('std')]
    for idx in stds.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", stds.loc[idx])))

    data.append(["Cross Correlations"] + [""] * summ.shape[1])
    cross_corrs = summ.loc[summ.index.str.startswith('Cor(') & ~summ.index.str.contains('_{-1}')]
    for idx in cross_corrs.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", cross_corrs.loc[idx])))

    data.append(["Autocorrelations"] + [""] * summ.shape[1])
    autocorrs = summ.loc[summ.index.str.contains('_{-1}')]
    for idx in autocorrs.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", autocorrs.loc[idx])))

    return tabulate(data, headers=headers, tablefmt="plain")


# ===========================================================================
# SMM Moment Table -- three-block design
#   Block 1 : standard deviations  (sigma)
#   Block 2 : contemporaneous correlation matrix  (lower triangular)
#   Block 3 : first-order autocorrelations  (rho_1)
# ===========================================================================

# Display metadata for each variable key used in observables_moments.py.
# Format: key -> (short_label, long_display_name, LaTeX_symbol)
_VAR_META = {
    "u":     ("u",  "Unemployment",          "$u$"),
    "v":     ("v",  "Vacancies",             "$v$"),
    "s":     ("s",  "Separations",           "$s$"),
    "jf":    ("f",  "Job-finding rate",      "$f$"),
    "delta": ("d",  "Exit rate",             "$\\delta$"),
    "ba":    ("a",  "Business applications", "$a$"),
    "lp":    ("z",  "Labor productivity",    "$z$"),
}


def compute_moments_matrix(cycle, var_labels):
    """
    Compute the three-block SMM moment matrix.

    Parameters
    ----------
    cycle : pd.DataFrame
        HP-filtered log-level cycle components, columns superset of var_labels.
    var_labels : list[str]
        Ordered variable keys; must be columns of cycle.

    Returns
    -------
    mom : pd.DataFrame
        Rows indexed by var_labels.
        Columns: ['SD'] + var_labels + ['AC']
        Correlation block is lower-triangular; diagonal and upper = NaN.
    """
    sub = cycle[var_labels]
    n   = len(var_labels)

    sds  = sub.std()
    corr = sub.corr()
    acs  = pd.Series(
        [sub[v].autocorr() for v in var_labels],
        index=var_labels,
    )

    # Lower-triangular correlation array (upper triangle + diagonal -> NaN)
    arr = corr.values.astype(float)
    for i in range(n):
        for j in range(i, n):
            arr[i, j] = np.nan

    mom = pd.DataFrame(index=var_labels)
    mom["SD"] = sds.values
    for j, v in enumerate(var_labels):
        mom[v] = arr[:, j]
    mom["AC"] = acs.values

    return mom


def format_moments_text(mom, var_labels, derived=None):
    """
    Format the moment matrix as a clean, shareable console table.

    Layout per row:  long name  |  sigma(x)  |  lower-tri correlations  |  rho_1(x)
    Upper-triangle cells are blank; diagonal cells show '  1  '.

    Parameters
    ----------
    mom : pd.DataFrame
        Output of compute_moments_matrix().
    var_labels : list[str]
        Ordered variable keys.
    derived : list of (str, str) or None
        Optional derived statistics appended below the main rows.
        Each entry is (text_label, value_str), e.g. ("sigma(theta)/sigma(z)", "38.4").

    Returns
    -------
    str : Multi-line formatted table.
    """
    n     = len(var_labels)
    short = {v: _VAR_META[v][0] if v in _VAR_META else v for v in var_labels}
    long_ = {v: _VAR_META[v][1] if v in _VAR_META else v for v in var_labels}

    rows = []
    for i, vi in enumerate(var_labels):
        row = [long_[vi], f"{mom.loc[vi, 'SD']:.4f}"]
        for j, vj in enumerate(var_labels):
            val = mom.loc[vi, vj]
            if j > i:
                row.append("")
            elif j == i:
                row.append("  1  ")
            else:
                row.append(f"{val: .3f}")
        row.append(f"{mom.loc[vi, 'AC']:.3f}")
        rows.append(row)

    headers = (
        ["Variable", "sigma(x)"]
        + [short[v] for v in var_labels]
        + ["rho_1(x)"]
    )
    col_align = ("left",) + ("right",) * (n + 2)

    table_body = tabulate(rows, headers=headers, tablefmt="simple",
                          colalign=col_align)

    derived_lines = ""
    if derived:
        sep = "-" * len(table_body.splitlines()[1])
        derived_rows = [
            [label, val_str] + ["  --  "] * n + ["  --  "]
            for label, val_str in derived
        ]
        derived_body = tabulate(derived_rows, tablefmt="plain",
                                colalign=col_align)
        derived_lines = "\n" + sep + "\n Derived statistics\n" + derived_body

    title = "Business Cycle Moments  [HP-filter, lambda=100,000]"
    corr_w = max(6 * n, len("Contemporaneous correlations"))
    banner = (
        "\n" + "=" * len(title) + "\n"
        + title + "\n"
        + "=" * len(title) + "\n"
        + f"{'':28s}  {'':8s}  "
        + ("-- Contemporaneous correlations --").center(corr_w)
        + f"  {'':8s}"
    )
    return banner + "\n" + table_body + derived_lines


def format_moments_latex(mom, var_labels,
                          caption="Business Cycle Moments",
                          label="tab:smm_moments",
                          note=None,
                          derived=None):
    """
    Generate a publication-quality LaTeX table (booktabs) with three blocks:

        sigma(x)  |  lower-triangular contemporaneous correlations  |  rho_1(x)

    Requires: booktabs (always), threeparttable (if note is provided).

    Parameters
    ----------
    mom : pd.DataFrame
        Output of compute_moments_matrix().
    var_labels : list[str]
        Ordered variable keys.
    caption : str
    label : str
    note : str or None
        Table note in footnotesize below the rule.
    derived : list of (str, str) or None
        Derived statistics appended after main rows, separated by a midrule.
        Each entry is (latex_label, value_str),
        e.g. (r"$\\sigma(\\theta)/\\sigma(z)$", "38.4").
        Value appears in the SD column; all other cells are blank.

    Returns
    -------
    str : LaTeX source for the complete table environment.
    """
    n       = len(var_labels)
    sym     = {v: _VAR_META[v][2] if v in _VAR_META else ("$" + v + "$")
               for v in var_labels}
    long_nm = {v: _VAR_META[v][1] if v in _VAR_META else v
               for v in var_labels}

    # Column spec: name | SD | r*n corr | AC
    vline_sep = r" @{\hspace{6pt}\vline\hspace{6pt}}"
    col_spec  = "l r" + vline_sep + " r" * n + vline_sep + " r"

    corr_start = 3
    corr_end   = n + 2

    L = []
    L.append(r"\begin{table}[htbp]")
    L.append(r"\centering")
    if note:
        L.append(r"\begin{threeparttable}")
    L.append("\\caption{" + caption + "}")
    L.append("\\label{" + label + "}")
    L.append("\\begin{tabular}{" + col_spec + "}")
    L.append(r"\toprule")

    # Header row 1: block span labels
    mc_corr = ("\\multicolumn{" + str(n) + "}{c}"
               + "{Contemporaneous correlations}")
    L.append("  &  $\\sigma(x)$  &  " + mc_corr
             + "  &  $\\rho_1(x)$  \\\\")
    L.append("\\cmidrule(lr){"
             + str(corr_start) + "-" + str(corr_end) + "}")

    # Header row 2: variable symbols
    sym_row = ("  &  &  "
               + "  &  ".join([sym[v] for v in var_labels])
               + "  &  \\\\")
    L.append(sym_row)
    L.append(r"\midrule")

    # Data rows
    for i, vi in enumerate(var_labels):
        sd_str = f"{mom.loc[vi, 'SD']:.4f}"
        ac_str = f"{mom.loc[vi, 'AC']:.3f}"
        corr_cells = []
        for j in range(n):
            if j > i:
                corr_cells.append("")
            elif j == i:
                corr_cells.append("$1$")
            else:
                val = mom.iloc[i, 1 + j]
                corr_cells.append(f"{val:.3f}")
        corr_str = "  &  ".join(corr_cells)
        L.append("  " + long_nm[vi]
                 + "  &  " + sd_str
                 + "  &  " + corr_str
                 + "  &  " + ac_str + "  \\\\")

    # Derived statistics block
    if derived:
        L.append(r"\midrule")
        ncols_str = str(n + 3)
        L.append("  \\multicolumn{" + ncols_str + "}{l}"
                 + "{\\textit{Derived statistics}}  \\\\")
        blank_corr = "  &  ".join([""] * n)
        for latex_label, val_str in derived:
            L.append("  \\quad " + latex_label
                     + "  &  " + val_str
                     + "  &  " + blank_corr
                     + "  &    \\\\")

    L.append(r"\bottomrule")
    L.append(r"\end{tabular}")

    if note:
        L.append(r"\begin{tablenotes}[flushleft]")
        L.append(r"\footnotesize")
        L.append("\\item \\textit{Note:} " + note)
        L.append(r"\end{tablenotes}")
        L.append(r"\end{threeparttable}")

    L.append(r"\end{table}")
    return "\n".join(L)


def generate_latex_subtables(df):
    "Generate high-quality table comparing model and empirical moments."
    pat_lag = '_{-1}'
    std_devs_autocorr = df[
        df.index.str.startswith('std') | df.index.str.contains(pat_lag)
    ]
    correlations = df[
        df.index.str.startswith('Cor')
        & ~df.index.str.contains('theta')
        & ~df.index.str.contains(pat_lag)
    ]

    latex_code = (
        "\n\\begin{table}[h]\n\\centering\n\\caption{Model Moments}\n"
        "\\begin{minipage}{0.45\\linewidth}\n\\centering\n"
        "\\subcaption{Standard Deviations and Autocorrelations}\n"
        "\\begin{tabular}{lcc}\n\\hline\n"
        "\\textbf{Moment} & \\textbf{Data Value} & \\textbf{Model Value} \\\\\n"
        "\\hline\n"
    )
    for idx, row in std_devs_autocorr.iterrows():
        latex_code += f"{idx} & {row['Values']:.3f} & {row['Model_moments']:.3f} \\\\\n"
    latex_code += (
        "\\hline\n\\end{tabular}\n\\end{minipage}%\n\\hfill\n"
        "\\begin{minipage}{0.45\\linewidth}\n\\centering\n"
        "\\subcaption{Contemporaneous Correlations}\n"
        "\\begin{tabular}{lcc}\n\\hline\n"
        "\\textbf{Moment} & \\textbf{Data Value} & \\textbf{Model Value} \\\\\n"
        "\\hline\n"
    )
    for idx, row in correlations.iterrows():
        latex_code += f"{idx} & {row['Values']:.3f} & {row['Model_moments']:.3f} \\\\\n"
    latex_code += "\\hline\n\\end{tabular}\n\\end{minipage}\n\\end{table}\n"
    return latex_code
