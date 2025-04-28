

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