from tabulate import tabulate

def create_stats_table(data, caption="Statistical Summary", label="tab:stats"):
     """
     Create a LaTeX table with dynamic column handling
     
     Parameters:
     -----------
     data : pandas.DataFrame
         DataFrame containing the statistical measures
     caption : str
         Table caption
     label : str
         Table reference label
     """
     # Get column names dynamically from the DataFrame
     columns = data.columns
     
     latex_str = [
         "\\begin{table}[htbp]",
         "\\centering",
         f"\\caption{{{caption}}}",
         f"\\label{{{label}}}",
     #    "\\begin{threeparttable}",
         # Create dynamic column format based on number of columns
         f"\\begin{{tabular}}{{l{''.join(['r'] * len(columns))}}}",
         "\\toprule"
     ]
     
     # Create header row dynamically
     # Replace potentially problematic characters and add LaTeX formatting
     header_row = ["Variable"] + [
         col.replace("_", "\\_")  # Escape underscores
            .replace("-", "$-$")  # Format minus signs
            .replace("(", "\\left(").replace(")", "\\right)")  # Format parentheses
         for col in columns
     ]
     latex_str.append(" & ".join(header_row) + " \\\\")
     
     latex_str.append("\\midrule")
     
     # Add data rows with proper formatting
     for idx, row in data.iterrows():
         formatted_row = [
             f"{idx}"  # Variable name
         ] + [
             f"{val:.3f}" if isinstance(val, (int, float)) else str(val)
             for val in row
         ]
         latex_str.append(" & ".join(formatted_row) + " \\\\")
     
     latex_str.extend([
         "\\bottomrule",
         "\\end{tabular}",
       #  "\\end{threeparttable}",
         "\\end{table}"
         ])
 
     return "\n".join(latex_str)
 
def generate_stacked_moments_latex_table(summ):
    """
    Generate LaTeX code for a table with subheadings for standard deviations,
    cross correlations, and autocorrelations using booktabs for publication quality.
    Handles any number of columns, dynamically generating column headers.

    Parameters:
    summ (pd.DataFrame): DataFrame containing the stacked moments with index labels indicating their type.

    Returns:
    str: A string containing the LaTeX code for the table.
    """
    # Generate headers for each column dynamically
    headers = [f"Value{i+1}" for i in range(summ.shape[1])]
    
    # Initialize LaTeX table string with booktabs formatting
    header_line = " & ".join([f"\\textbf{{{header}}}" for header in headers])
    latex_table = f"""
    \\begin{{table}}[h]
    \\centering
    \\begin{{tabular}}{{l{'r' * len(headers)}}}
    \\toprule
    \\textbf{{Moment Type}} & {header_line} \\\\
    \\midrule
    \\multicolumn{{{len(headers) + 1}}}{{c}}{{\\textbf{{Standard Deviations}}}} \\\\
    \\midrule
    """
    
    # Add standard deviations
    stds = summ.loc[summ.index.str.startswith('std')]
    for idx in stds.index:
        values_line = " & ".join([f"{val:.3f}" for val in stds.loc[idx]])
        latex_table += f"{idx} & {values_line} \\\\ \n"
    
    latex_table += f"\\midrule \\multicolumn{{{len(headers) + 1}}}{{c}}{{\\textbf{{Cross Correlations}}}} \\\\ \\midrule \n"
    
    # Add cross correlations
    cross_corrs = summ.loc[summ.index.str.startswith('Cor(') & ~summ.index.str.contains('_{-1}')]
    for idx in cross_corrs.index:
        values_line = " & ".join([f"{val:.3f}" for val in cross_corrs.loc[idx]])
        latex_table += f"{idx} & {values_line} \\\\ \n"
    
    latex_table += f"\\midrule \\multicolumn{{{len(headers) + 1}}}{{c}}{{\\textbf{{Autocorrelations}}}} \\\\ \\midrule \n"
    
    # Add autocorrelations
    autocorrs = summ.loc[summ.index.str.contains('_{-1}')]
    for idx in autocorrs.index:
        values_line = " & ".join([f"{val:.3f}" for val in autocorrs.loc[idx]])
        latex_table += f"{idx} & {values_line} \\\\ \n"
    
    latex_table += f"\\bottomrule \\end{{tabular}} \\caption{{Stacked Moments Summary}} \\end{{table}}"
    
    return latex_table

def generate_stacked_moments_table(summ):
    """
    Generate a plain text table using the tabulate package for moments, including subheadings.

    Parameters:
    summ (pd.DataFrame): DataFrame containing the stacked moments with index labels indicating their type.

    Returns:
    str: A string containing the formatted text table.
    """
    # Prepare headers
    headers = ["Moment Type"] + [f"Value{i+1}" for i in range(summ.shape[1])]
    
    # Initialize list for table data
    data = []
    
    # Add standard deviations subheading and data
    data.append(["Standard Deviations"] + [""] * summ.shape[1])
    stds = summ.loc[summ.index.str.startswith('std')]
    for idx in stds.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", stds.loc[idx])))
    
    # Add cross correlations subheading and data
    data.append(["Cross Correlations"] + [""] * summ.shape[1])
    cross_corrs = summ.loc[summ.index.str.startswith('Cor(') & ~summ.index.str.contains('_{-1}')]
    for idx in cross_corrs.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", cross_corrs.loc[idx])))
    
    # Add autocorrelations subheading and data
    data.append(["Autocorrelations"] + [""] * summ.shape[1])
    autocorrs = summ.loc[summ.index.str.contains('_{-1}')]
    for idx in autocorrs.index:
        data.append([idx] + list(map(lambda x: f"{x:.3f}", autocorrs.loc[idx])))
    
    # Use tabulate to format the table
    return tabulate(data, headers=headers, tablefmt="plain")


def generate_latex_subtables(df):
    " Generate high-quality table comparing model and empirical moments "
    std_devs_autocorr = df[df.index.str.startswith('std') | df.index.str.contains('_{-1}')]
    correlations = df[df.index.str.startswith('Cor') & ~df.index.str.contains('theta') & ~df.index.str.contains('_{-1}')]

    latex_code = r"""
\begin{table}[h]
\centering
\caption{Model Moments}
\begin{minipage}{0.45\linewidth}
\centering
\subcaption{Standard Deviations and Autocorrelations}
\begin{tabular}{lcc}
\hline
\textbf{Moment} & \textbf{Data Value} & \textbf{Model Value} \\
\hline
"""
    for idx, row in std_devs_autocorr.iterrows():
        latex_code += f"{idx} & {row['Values']:.3f} & {row['Model_moments']:.3f} \\\\\n"

    latex_code += r"""
\hline
\end{tabular}
\end{minipage}%
\hfill
\begin{minipage}{0.45\linewidth}
\centering
\subcaption{Contemporaneous Correlations}
\begin{tabular}{lcc}
\hline
\textbf{Moment} & \textbf{Data Value} & \textbf{Model Value} \\
\hline
"""
    for idx, row in correlations.iterrows():
        latex_code += f"{idx} & {row['Values']:.3f} & {row['Model_moments']:.3f} \\\\\n"

    latex_code += r"""
\hline
\end{tabular}
\end{minipage}
\end{table}
"""
    return latex_code