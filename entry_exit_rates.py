import pandas as pd
pd.set_option('display.max_columns', 10) 
import numpy as np
from statsmodels.tsa.stattools import adfuller
import matplotlib.pyplot as plt
from IPython.display import display
from time_series_functions import crosscorr
from fredapi import Fred
from datetime import datetime
import matplotlib.dates as mdates

# Data

df = pd.read_csv('BDS_Extension.csv')
fred = Fred(api_key="8d302fb5f121b2be4f7d6e795194bcde")

unemployment = fred.get_series('UNRATE')
unemployment_annual = unemployment.groupby(unemployment.index.year).mean()
unemployment_df = pd.DataFrame({
    'year': unemployment_annual.index,
    'unemployment_rate': unemployment_annual.values
})

df = pd.merge(df, unemployment_df, on='year', how='left')
#df_adf = df[df['year'] >= 1979].reset_index(drop=True)
df['year'] = pd.to_datetime(df['year'], format='%Y')
df.set_index("year", inplace=True)
# Episodes

recessions = [
    (pd.Timestamp('1980-01-01'), pd.Timestamp('1980-07-31')),
    (pd.Timestamp('1981-07-01'), pd.Timestamp('1982-11-30')),
    (pd.Timestamp('1990-07-01'), pd.Timestamp('1991-03-31')),
    (pd.Timestamp('2001-03-01'), pd.Timestamp('2001-11-30')),
    (pd.Timestamp('2007-12-01'), pd.Timestamp('2009-06-30')),
    (pd.Timestamp('2020-02-01'), pd.Timestamp('2020-04-30')),
]


# Define

variables = ['job_creation_rate', 'job_destruction_rate', 'estabs_entry_rate', 'estabs_exit_rate', 
            'firms_entry_rate', 'firms_exit_rate', 'unemployment_rate']
variable_names = ['Job Creation Rate', 'Job Destruction Rate', 'Establishment Entry Rate', 
                 'Establishment Exit Rate', 'Firm Entry Rate', 'Firm Exit Rate', 'Unemployment Rate']
plot_variables = ['estabs_entry_rate', 'estabs_exit_rate', 'firms_entry_rate', 
                 'firms_exit_rate', 'unemployment_rate']
plot_variable_names = ['Establishment Entry Rate', 'Establishment Exit Rate', 'Firm Entry Rate', 
                      'Firm Exit Rate', 'Unemployment Rate']


# Plot

fig, ax = plt.subplots(figsize=(12, 8), nrows=2)
ax[0].plot(df.index, df['estabs_entry_rate'], label='Establishments', linewidth=2)
ax[0].plot(df.index, df['firms_entry_rate'], label='Firms', linewidth=2)
ax[0].plot(df.index, df['job_creation_rate'], label='Job creation', linewidth=2)
# ax[0].plot(df.year, df['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[0].set_title("Entry rates", fontsize=12, pad=10)

ax[1].plot(df.index, df['estabs_exit_rate'], label='Establishments', linewidth=2)
ax[1].plot(df.index, df['firms_exit_rate'], label='Firms', linewidth=2)
ax[1].plot(df.index, df['job_destruction_rate'], label='Job destruction', linewidth=2)
#ax[1].plot(df.year, df['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[1].set_title("Exit rates", fontsize=12, pad=10)

       
for j in range(2):
    ax[j].set_xlabel('Year', fontsize=10)
    ax[j].set_ylabel('Rate (%)', fontsize=10)
    ax[j].set_xticks(range(1978, 2022, 5))
    for start, end in recessions:
        ax[j].axvspan(start, end, color='gray', alpha=0.3)
    ax[j].legend(loc='upper right', fontsize=10)
    ax[j].grid(True, alpha=0.3)
    #ax[j].tick_params(axis='both', which='major', labelsize=9)

plt.tight_layout()
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.savefig("entry_exit_rates.pdf", bbox_inches='tight', dpi=300)
plt.show()

# Sum

# def create_summary_table(variable, var_name): 
#     summary = df[variable].describe()
#     summary.loc['mean'] = df[variable].mean()
#     return summary.to_frame(var_name)  

# summary_tables = pd.concat([create_summary_table(var, name) for var, name in zip(variables, variable_names)], axis=1)
# summary_tables.columns = variable_names
# print("\nSummary Statistics:")
# display(summary_tables)

df_red = df [['estabs_exit_rate', 'firms_exit_rate', 'job_creation_rate', 'job_destruction_rate']]
df_red.mean()
df_red.corr()

def create_latex_summary_table(df, caption="Summary Statistics", label="tab:summary"):
    """
    Create a LaTeX table with means and correlations of a dataframe.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with numerical columns
    caption : str
        Caption for the LaTeX table
    label : str
        Label for referencing the table in LaTeX
        
    Returns:
    --------
    str
        LaTeX formatted table string
    """
    # Calculate means
    means = df.mean().round(3)
    
    # Calculate correlation matrix
    corr_matrix = df.corr().round(3)
    
    # Create the LaTeX table string
    latex_str = "\\begin{table}[htbp]\n"
    latex_str += "\\centering\n"
    latex_str += "\\caption{" + caption + "}\n"
    latex_str += "\\label{" + label + "}\n"
    
    # Panel A: Means
    latex_str += "\\textbf{Panel A: Means}\\\\\n"
    latex_str += means.to_frame().to_latex(header=False)
    latex_str += "\\\\\n"
    
    # Panel B: Correlations
    latex_str += "\\textbf{Panel B: Correlation Matrix}\\\\\n"
    latex_str += corr_matrix.to_latex()
    
    latex_str += "\\end{table}"
    
    return latex_str

sum_table = create_latex_summary_table(df_red)
print(sum_table)
# Comovement

comovement_pairs = [
    ('job_creation_rate', 'estabs_entry_rate'),
    ('firms_entry_rate', 'estabs_entry_rate'),
    ('job_destruction_rate', 'estabs_exit_rate'),
    ('firms_exit_rate', 'estabs_exit_rate'),
    ('unemployment_rate', 'estabs_entry_rate'),
    ('unemployment_rate', 'firms_entry_rate'),
    ('unemployment_rate', 'estabs_exit_rate'),
    ('unemployment_rate', 'firms_exit_rate')
]

# for var1, var2 in comovement_pairs:
#     correlation = df[[var1, var2]].corr().iloc[0, 1]
#     print(f"\nCorrelation between {variable_names[variables.index(var1)]} and {variable_names[variables.index(var2)]}: {correlation:.2f}")
    
def dynamic_correlations(data, var1, var2, ylabel, nleads=12, nlags=12, title=None):
    fig, ax = plt.subplots(figsize=(14, 5))
    rs = []
    lags = range(-nlags, nleads+1)
    
    
    for lag in lags:
        rs.append(crosscorr(data[var1], data[var2], lag))
    rs = pd.Series(rs)
    
   
    max_corr_idx = np.argmax(abs(rs))
    max_corr_value = rs[max_corr_idx]
    max_corr_lag = lags[max_corr_idx]
    
    
    ax.axhline(y=0.0, color="black", linestyle="--")
    ax.plot(range(len(rs)), rs, '-', alpha=0.7, linewidth=2.0)
    ax.axvline(max_corr_idx, linestyle='--', color='red')
    
    
    ax.text(max_corr_idx, max(rs) + 0.1, 
            f'Max corr: {max_corr_value:.2f}\nLag: {max_corr_lag}', 
            horizontalalignment='center')
    
    ax.set_xlabel(r'$\Delta$ (years)', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(range(0, len(rs)))
    ax.set_xticklabels(lags)
    ax.grid(True, alpha=0.3)
    
    if title is not None:
        ax.set_title(title, fontsize=14)
    plt.tight_layout()
    plt.show()

# for var1, var2 in comovement_pairs:
#     name1 = variable_names[variables.index(var1)]
#     name2 = variable_names[variables.index(var2)]
#     ylabel = f'Corr({name1}, {name2})'
#     title = f'Dynamic Correlations between {name1} and {name2}'
#     dynamic_correlations(df, var1, var2, ylabel, nleads=5, nlags=5, title=title)

# # Stationary Test
# def adf_test(series, name):
#     result = adfuller(series)
#     print(f'ADF Statistic for {name}: {result[0]:.4f}')
#     print(f'p-value: {result[1]:.4f}')
#     print(f'Critical Values:')
#     for key, value in result[4].items():
#         print(f'   {key}: {value:.4f}')
#     if result[1] <= 0.05:
#         print(f'{name} is stationary.\n')
#     else:
#         print(f'{name} is non-stationary.\n')

# print("Stationarity Tests:")
# for var, name in zip(variables, variable_names):
#     adf_test(df_adf[var], name)
    
