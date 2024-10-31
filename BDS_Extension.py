import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import adfuller
import matplotlib.pyplot as plt
from IPython.display import display
from time_series_functions import crosscorr

# Data

df = pd.read_csv('BDS_Extension.csv')
df_adf = df[df['year'] >= 1979].reset_index(drop=True)

# Episodes

periods = {
    'Volcker Rate Rise': (1979, 1987),
    'Great Recession': (2007, 2009),
    'Pandemic Covid': (2020, 2023)
}

# Define

variables = ['job_creation_rate', 'job_destruction_rate', 'estabs_entry_rate', 'estabs_exit_rate', 'firms_entry_rate', 'firms_exit_rate']
variable_names = ['Job Creation Rate', 'Job Destruction Rate', 'Establishment Entry Rate', 'Establishment Exit Rate', 'Firm Entry Rate', 'Firm Exit Rate']
plot_variables = ['estabs_entry_rate', 'estabs_exit_rate', 'firms_entry_rate', 'firms_exit_rate']
plot_variable_names = ['Establishment Entry Rate', 'Establishment Exit Rate', 'Firm Entry Rate', 'Firm Exit Rate']


# Plot

fig, ax = plt.subplots(figsize=(12, 8), nrows=2)
ax[0].plot(df.year, df['estabs_entry_rate'], label='Establishments')
ax[0].plot(df.year, df['firms_entry_rate'], label='Firms')
ax[0].set_title("Entry rates")
ax[1].plot(df.year, df['estabs_exit_rate'], label='Establishments')
ax[1].plot(df.year, df['firms_exit_rate'], label='Firms')
ax[1].set_title("Exit rates")

for j in range(2):
    ax[j].set_xlabel('Year')
    ax[j].set_ylabel('Rate (%)')
    ax[j].set_xticks(range(1978, 2022, 5))
    for start, end in periods.values():
        ax[j].axvspan(start, end, color='gray', alpha=0.3)
    ax[j].legend(loc='upper right')
    ax[j].grid(True)
plt.tight_layout()
plt.savefig("entry_exit_rates.pdf")
plt.show()

# Sum

def create_summary_table(variable):
    summary = df[variable].describe()
    summary.loc['mean'] = df[variable].mean()
    return summary.to_frame(name)

summary_tables = pd.concat([create_summary_table(var) for var in variables], axis=1)
summary_tables.columns = variable_names
print("\nSummary Statistics:")
display(summary_tables)

# Comovement

comovement_pairs = [
    ('job_creation_rate', 'estabs_entry_rate'),
    ('firms_entry_rate', 'estabs_entry_rate'),
    ('job_destruction_rate', 'estabs_exit_rate'),
    ('firms_exit_rate', 'estabs_exit_rate')
]

for var1, var2 in comovement_pairs:
    correlation = df[[var1, var2]].corr().iloc[0, 1]
    print(f"\nCorrelation between {variable_names[variables.index(var1)]} and {variable_names[variables.index(var2)]}: {correlation:.2f}")
    
def dynamic_correlations(data, var1, var2, ylabel, nleads=12, nlags=12, title=None):
    fig, ax = plt.subplots(figsize=(14, 5))
    x = data[var1]
    y = data[var2]
    rs = pd.Series([crosscorr(data[var1], data[var2], -lag) for lag in range(-nlags, nleads)])
    ax.axhline(y=0.0, color="black", linestyle="--")
    ax.plot(rs, '-', alpha=0.7, linewidth=2.0)
    ax.axvline(np.argmax(abs(rs)), linestyle='--', color='red')
    ax.set_xlabel(r'$\Delta$ (years)', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(range(0, len(rs)))
    ax.set_xticklabels(range(-nlags, nleads))
    if title is not None:
        ax.set_title(title, fontsize=14)
    plt.tight_layout()
    plt.show()

for var1, var2 in comovement_pairs:
    name1 = variable_names[variables.index(var1)]
    name2 = variable_names[variables.index(var2)]
    ylabel = f'Corr({name1}, {name2})'
    title = f'Dynamic Correlations between {name1} and {name2}'
    dynamic_correlations(df, var1, var2, ylabel, title=title)    

# Staionary Test

def adf_test(series, name):
    result = adfuller(series)
    print(f'ADF Statistic for {name}: {result[0]:.4f}')
    print(f'p-value: {result[1]:.4f}')
    print(f'Critical Values:')
    for key, value in result[4].items():
        print(f'   {key}: {value:.4f}')
    if result[1] <= 0.05:
        print(f'{name} is stationary.\n')
    else:
        print(f'{name} is non-stationary.\n')

print("Stationarity Tests:")
for var, name in zip(variables, variable_names):
    adf_test(df_adf[var], name)
    
