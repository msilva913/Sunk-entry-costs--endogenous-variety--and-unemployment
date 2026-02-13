
" Investigate vacancy-business dynamic facts "
" 1) Higher establishment number leads to faster vacancy creation"
" 2) More firm exit leads to slower vacancy creation "

""" BDS Data Tables 

Table 1. Private sector gross job gains and losses, seasonally adjusted
Table 2. Private sector gross job gains and losses as a percent of employment (1), seasonally adjusted
Table 3. Private sector gross job gains and losses by industry, seasonally adjusted
Table 4. Private sector gross job gains and losses by firm size, seasonally adjusted
Table 5. Components of private sector gross job gains and losses by firm size, seasonally adjusted
Table 6. Private sector gross job gains and losses by state, seasonally adjusted
Table 7. Private sector gross job gains and losses as a percent of total employment by state, seasonally adjusted
Table 8. Private sector establishment births and deaths, seasonally adjusted

"""
# Initial imports 
from fredapi import Fred
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
fred = Fred(api_key='9c70445138df124be4928605b7e08bd4')
import os
os.chdir(r"C:\Users\msilv\Documents\GitHub\Sunk-entry-costs--endogenous-variety--and-unemployment\Data")
load = True
from time_series_functions import filter_transform, moments

if load == False:
    # Dictionary of industries for the seasonally adjusted level of vacancies
    industries_level = {
        'JTSJOL': 'Total Nonfarm',
        'JTS1000JOL': 'Total Private',
        'JTS2300JOL': 'Construction',
        'JTS3000JOL': 'Manufacturing',
        'JTS3200JOL': 'Durable Goods Manufacturing',
        'JTS3400JOL': 'Nondurable Goods Manufacturing',
        'JTS4000JOL': 'Trade, Transportation, and Utilities',
        'JTS4400JOL': 'Retail Trade',
        'JTS540099JOL': 'Professional and Business Services',
        'JTS6000JOL': 'Private Education and Health Services',
        'JTS6200JOL': 'Health Care and Social Assistance',
        'JTS7000JOL': 'Leisure and Hospitality',
        'JTS7100JOL': 'Arts, Entertainment, and Recreation',
        'JTS7200JOL': 'Accommodation and Food Services',
        'JTS9000JOL': 'Government',  
        'JTS9200JOL': 'State and Local',
    }


    
    data_list = []
    for series_id, name in industries_level.items():
        series_data = fred.get_series(series_id)
        series_data.name = name
        data_list.append(series_data)
    
    # Combine into one DataFrame
    df = pd.concat(data_list, axis=1)

    #Save to CSV
    df.to_csv('jolts_vacancies_industry_adjusted_levels.csv')
else:
    df = pd.read_csv('jolts_vacancies_industry_adjusted_levels.csv')

print(df.shape, df.columns)
df.head()

# Set datetimeindex 
df.rename(columns={"Unnamed: 0": "date"}, inplace=True)
df.index = pd.to_datetime(df.date, format='%Y-%m-%d')
df.drop(["date"], axis=1, inplace=True)

# Select industry groups 
industry_gr = ["Nondurable Goods Manufacturing", 
            "Durable Goods Manufacturing",
            "Construction",
            "Trade, Transportation, and Utilities",
            "Retail Trade",
            "Professional and Business Services"]

# Reduced dataset corresponding to select industry groups
df_red = df[industry_gr]
ind_means = df_red.apply(np.mean, axis=0)
ind_means.sort_values(ascending=False, inplace=True)
print(ind_means)

fig, ax = plt.subplots()
sns.barplot(data=ind_means, ax=ax)
plt.xticks(rotation=45, ha="right") # rotate and align labels
plt.tight_layout()

# Time series plots 
# 1) Raw plots
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(14, 6), 
                       sharex=True, sharey=False)
axes_flat = axes.ravel()
for i, col in enumerate(df_red.columns):
    ax = axes_flat[i] # flatten the axes
    ax.plot(df_red.index, df_red[col], linestyle='-', alpha=0.9, linewidth=1.5, color='forestgreen')
    ax.set_title(col)
plt.tight_layout()
plt.show()

" Apply logarithm and construct cylical deviation"
init = df_red.index[0] 
final = df_red.index[-1]

cycle = pd.concat([filter_transform(df_red[x], init=init, final=final, transform_type='log',
                        filter_type="hp_filter", lamb=100_000) for x in industry_gr], axis=1)
cycle.columns = df_red.columns
# 2) Plots of cyclical deviations 
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(14, 6), 
                       sharex=True, sharey=False)
axes_flat = axes.ravel()

for i, col in enumerate(cycle.columns):
    ax = axes_flat[i] # flatten the axes
    ax.plot(cycle.index, cycle[col], linestyle='-', alpha=0.9, linewidth=1.5, color="forestgreen")
    ax.set_title(col)
plt.tight_layout()
plt.show()

# 3) Moments 
mom = moments(cycle, relative_std="Nondurable Goods Manufacturing", lab=["Nondurable Goods Manufacturing"])

