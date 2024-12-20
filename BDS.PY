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
from BLSAPI import bls
fred = Fred(api_key="8d302fb5f121b2be4f7d6e795194bcde")

# Data

BDS = pd.read_csv('BDS.csv')

birth_rate = bls.get_series("BDS0000000000000000120007RQ5", freq="Q").resample("A").mean().dropna() * 4
death_rate = bls.get_series("BDS0000000000000000120008RQ5", freq="Q").resample("A").mean().dropna() * 4

birth_raw= bls.get_series("BDS0000000000000000120007RQ5", freq="Q").dropna()
death_raw = bls.get_series("BDS0000000000000000120008RQ5", freq="Q").dropna()

BED = pd.DataFrame({
    "year": birth_rate.index.year,
    "establishments_birth_rate": birth_rate.values, # establishments_birth_rate: BDS0000000000000000120007RQ5 (1) as a percent of total establishments in this sector
    "establishments_death_rate": death_rate.values, # establishments_death_rate: BDS0000000000000000120008RQ5 (1) as a percent of total establishments in this sector
})

unemployment = fred.get_series('UNRATE')
unemployment_annual = unemployment.groupby(unemployment.index.year).mean()
unemployment_df = pd.DataFrame({
    'year': unemployment_annual.index,
    'unemployment_rate': unemployment_annual.values
})

# df

df_temp = pd.merge(BDS, unemployment_df, on='year', how='left')
df = pd.merge(df_temp, BED, on='year', how='left')
df_adf = df[df['year'] >= 1979].reset_index(drop=True)
df_bed = df[df['year'] >= 1993].reset_index(drop=True)

# Episodes

periods = {
    'Volcker Rate Rise': (1979, 1987),
    'Great Recession': (2007, 2009),
    'Pandemic Covid': (2020, 2023)
}

variables = ['job_creation_rate', 'job_destruction_rate', 'estabs_entry_rate', 'estabs_exit_rate', 
             'firms_entry_rate', 'firms_exit_rate', 'unemployment_rate', 
             'establishments_birth_rate', 'establishments_death_rate']
variable_names = ['Job Creation Rate', 'Job Destruction Rate', 'Establishment Entry Rate', 
                  'Establishment Exit Rate', 'Firm Entry Rate', 'Firm Exit Rate', 
                  'Unemployment Rate', 'Establishment Birth Rate', 'Establishment Death Rate']
plot_variables = ['estabs_entry_rate', 'estabs_exit_rate', 'firms_entry_rate', 
                  'firms_exit_rate', 'unemployment_rate', 
                  'establishments_birth_rate', 'establishments_death_rate']
plot_variable_names = ['Establishment Entry Rate', 'Establishment Exit Rate', 'Firm Entry Rate', 
                       'Firm Exit Rate', 'Unemployment Rate', 
                       'Establishment Birth Rate', 'Establishment Death Rate']

# Plot in whole

fig, ax = plt.subplots(figsize=(12, 8), nrows=2)
ax[0].plot(df.year, df['estabs_entry_rate'], label='Establishments', linewidth=2)
ax[0].plot(df.year, df['firms_entry_rate'], label='Firms', linewidth=2)
ax[0].plot(df.year, df['job_creation_rate'], label='Job creation', linewidth=2)
ax[0].plot(df.year, df['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[0].plot(df.year, df['establishments_birth_rate'], label='Estab_Birth_Rate', linestyle='--', color='blue')
ax[0].set_title("Entry rates, Birth Rate, and Unemployment", fontsize=12, pad=10)

ax[1].plot(df.year, df['estabs_exit_rate'], label='Establishments', linewidth=2)
ax[1].plot(df.year, df['firms_exit_rate'], label='Firms', linewidth=2)
ax[1].plot(df.year, df['job_destruction_rate'], label='Job destruction', linewidth=2)
ax[1].plot(df.year, df['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[1].plot(df.year, df['establishments_death_rate'], label='Estab_Death_Rate', linestyle='--', color='green')
ax[1].set_title("Exit rates, Death Rate, and Unemployment", fontsize=12, pad=10)

for j in range(2):
    ax[j].set_xlabel('Year', fontsize=10)
    ax[j].set_ylabel('Rate (%)', fontsize=10)
    ax[j].set_xticks(range(1978, 2022, 5))
    ax[j].legend(loc='upper right', fontsize=10)
    ax[j].grid(True, alpha=0.3)
    ax[j].tick_params(axis='both', which='major', labelsize=9)
    
    for name, (start, end) in periods.items():
        ax[j].axvspan(start, end, color='gray', alpha=0.3)

plt.tight_layout()
plt.savefig("entry_exit_birth_death_unemployment_rates.pdf", bbox_inches='tight', dpi=300)
plt.show()

# Plot variables in df_bed

fig, ax = plt.subplots(figsize=(12, 8), nrows=2)
ax[0].plot(df_bed.year, df_bed['estabs_entry_rate'], label='Establishments', linewidth=2)
ax[0].plot(df_bed.year, df_bed['firms_entry_rate'], label='Firms', linewidth=2)
ax[0].plot(df_bed.year, df_bed['job_creation_rate'], label='Job creation', linewidth=2)
ax[0].plot(df_bed.year, df_bed['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[0].plot(df_bed.year, df_bed['establishments_birth_rate'], label='Estab_Birth_Rate', linestyle='--', color='blue')
ax[0].set_title("Entry rates, Birth Rate, and Unemployment (After 1993)", fontsize=12, pad=10)

ax[1].plot(df_bed.year, df_bed['estabs_exit_rate'], label='Establishments', linewidth=2)
ax[1].plot(df_bed.year, df_bed['firms_exit_rate'], label='Firms', linewidth=2)
ax[1].plot(df_bed.year, df_bed['job_destruction_rate'], label='Job destruction', linewidth=2)
ax[1].plot(df_bed.year, df_bed['unemployment_rate'], label='Unemployment', linestyle='--', color='red')
ax[1].plot(df_bed.year, df_bed['establishments_death_rate'], label='Estab_Death_Rate', linestyle='--', color='green')
ax[1].set_title("Exit rates, Death Rate, and Unemployment (After 1993)", fontsize=12, pad=10)

target_periods = list(periods.items())[1:]
for j in range(2):
    ax[j].set_xlabel('Year', fontsize=10)
    ax[j].set_ylabel('Rate (%)', fontsize=10)
    ax[j].set_xticks(range(1993, 2024, 5))
    ax[j].legend(loc='upper right', fontsize=10)
    ax[j].grid(True, alpha=0.3)
    ax[j].tick_params(axis='both', which='major', labelsize=9)
    
    for name, (start, end) in target_periods:
     ax[j].axvspan(start, end, color='gray', alpha=0.3)

plt.tight_layout()
plt.savefig("entry_exit_birth_death_unemployment_rates_after_1993.pdf", bbox_inches='tight', dpi=300)
plt.show()

# Statistical Summary

df_red = df[['estabs_exit_rate', 'firms_exit_rate', 'job_creation_rate', 'job_destruction_rate', 
             'establishments_birth_rate', 'establishments_death_rate']]
print("Means:\n", df_red.mean())
print("\nCorrelations:\n", df_red.corr())

# Stationarity Test

print("\nStationarity Tests for df_bed Variables:")
for var, name in zip(variables, variable_names):
    if var in df_bed.columns:
        series = df_bed[var].dropna()  # Drop NaNs in the current series
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

# Comovement Analysis

comovement_pairs = [
    ('job_creation_rate', 'estabs_entry_rate'),
    ('firms_entry_rate', 'estabs_entry_rate'),
    ('job_destruction_rate', 'estabs_exit_rate'),
    ('firms_exit_rate', 'estabs_exit_rate'),
    ('unemployment_rate', 'estabs_entry_rate'),
    ('unemployment_rate', 'firms_entry_rate'),
    ('unemployment_rate', 'estabs_exit_rate'),
    ('unemployment_rate', 'firms_exit_rate'),
    ('unemployment_rate', 'establishments_birth_rate'),
    ('unemployment_rate', 'establishments_death_rate'),
]

print("\nComovement Analysis for df_bed Variables:")
for var1, var2 in comovement_pairs:
    if var1 in df_bed.columns and var2 in df_bed.columns:
        correlation = df_bed[[var1, var2]].corr().iloc[0, 1]
        name1 = variable_names[variables.index(var1)]
        name2 = variable_names[variables.index(var2)]
        print(f"Correlation between {name1} and {name2}: {correlation:.2f}")