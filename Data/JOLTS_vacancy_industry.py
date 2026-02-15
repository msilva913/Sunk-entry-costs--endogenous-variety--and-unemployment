
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

##############################
import requests
import json
import matplotlib.ticker as mticker

# --- Configuration ---
# 1. Enter your BLS API key here.
API_KEY = 'YOUR_API_KEY_HERE' 

# 2. Define the industries and data types you want to fetch.
# I've included the exact series IDs for your requested industries and more.
series_map = {
    # --- Establishment Births (Number) ---
    'BTH_CON': 'BDU000000101200000050012', # Construction Births
    'BTH_DUR': 'BDU0000003200000050012', # Durable Goods Mfg Births
    'BTH_NDR': 'BDU0000003100000050012', # Nondurable Goods Mfg Births
    'BTH_RT':  'BDU0000004400000050012', # Retail Trade Births
    'BTH_PBS': 'BDU0000006000000050012', # Professional & Business Services Births
    
    # --- Establishment Deaths (Number) ---
    'DTH_CON': 'BDU000000101200000050013', # Construction Deaths
    'DTH_DUR': 'BDU0000003200000050013', # Durable Goods Mfg Deaths
    'DTH_NDR': 'BDU0000003100000050013', # Nondurable Goods Mfg Deaths
    'DTH_RT':  'BDU0000004400000050013', # Retail Trade Deaths
    'DTH_PBS': 'BDU0000006000000050013', # Professional & Business Services Deaths
}

# 3. Define the time range for the data.
start_year = '2000'
end_year = '2023'

# --- API Request and Data Processing ---

def fetch_bls_data(series_dict, start, end, api_key):
    """Fetches data from the BLS API and returns a clean DataFrame."""
    if api_key == 'YOUR_API_KEY_HERE':
        print("Warning: Please replace 'YOUR_API_KEY_HERE' with your actual BLS API key.")
        return None

    headers = {'Content-type': 'application/json'}
    
    # The BLS API can accept up to 50 series in a single request
    data = json.dumps({
        "seriesid": list(series_dict.values()),
        "startyear": start,
        "endyear": end,
        "registrationkey": api_key,
        "catalog": True # Good practice to get metadata
    })
    
    # Send the request
    p = requests.post('https://api.bls.gov/publicAPI/v2/timeseries/data/', data=data, headers=headers)
    
    if p.status_code != 200:
        print(f"Error fetching data: {p.status_code}")
        print(p.json())
        return None
        
    json_data = p.json()
    
    if json_data['status'] != 'REQUEST_SUCCEEDED':
        print(f"BLS API Error: {json_data['message']}")
        return None

    # Process the JSON response into a DataFrame
    all_series_data = []
    series_titles = {s_id: series['seriesTitle'] for s_id, series in zip(series_dict.values(), json_data['Results']['series'])}

    for series_id, series_name in series_dict.items():
        try:
            series_data = json_data['Results']['series']
            # Find the specific series data from the response
            for s in series_data:
                if s['seriesID'] == series_name:
                    df = pd.DataFrame(s['data'])
                    df['seriesID'] = series_name
                    df['seriesName'] = series_id
                    df['value'] = pd.to_numeric(df['value'])
                    # Create a proper date index
                    df['date'] = pd.to_datetime(df['year'] + '-' + df['period'].str.replace('Q0', '').astype(int) * 3)
                    df = df[['date', 'value', 'seriesName']]
                    all_series_data.append(df)
                    break
        except KeyError:
            print(f"Could not find data for series: {series_name}")

    if not all_series_data:
        print("No data processed.")
        return None
        
    # Combine all series into a single pivot table
    combined_df = pd.concat(all_series_data)
    final_df = combined_df.pivot(index='date', columns='seriesName', values='value')
    final_df = final_df.sort_index()
    
    return final_df

# Fetch the data
bed_data = fetch_bls_data(series_map, start_year, end_year, API_KEY)


# --- Plotting the Results ---

if bed_data is not None:
    print("Successfully fetched data. Columns available:")
    print(bed_data.columns)
    
    # Define your color palette
    colors = {
        'births': '#3776ab', # Deep Azure
        'deaths': '#C0C0C0', # Silver/Gray
    }

    # Plot for Construction
    fig, ax = plt.subplots(figsize=(12, 7))
    
    bed_data['BTH_CON'].plot(ax=ax, color=colors['births'], label='Establishment Births', lw=2.5)
    bed_data['DTH_CON'].plot(ax=ax, color=colors['deaths'], label='Establishment Deaths', lw=2.5, linestyle='--')
    
    ax.set_title('Establishment Births and Deaths: Construction', fontsize=16, pad=20)
    ax.set_ylabel('Number of Establishments')
    ax.set_xlabel('Year')
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    ax.legend()
    
    # Format y-axis to have commas for thousands
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    plt.tight_layout()
    plt.show()

    # Plot for Manufacturing (Durable vs. Nondurable)
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharex=True)

    # Durable Goods
    bed_data['BTH_DUR'].plot(ax=ax1, color='#00008B', label='Births', lw=2) # Dark Blue
    bed_data['DTH_DUR'].plot(ax=ax1, color='#A9A9A9', label='Deaths', lw=2, linestyle='--') # Dark Gray
    ax1.set_title('Durable Goods Manufacturing', fontsize=14)
    ax1.set_ylabel('Number of Establishments')
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.6)
    ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    # Nondurable Goods
    bed_data['BTH_NDR'].plot(ax=ax2, color='#00008B', label='Births', lw=2)
    bed_data['DTH_NDR'].plot(ax=ax2, color='#A9A9A9', label='Deaths', lw=2, linestyle='--')
    ax2.set_title('Nondurable Goods Manufacturing', fontsize=14)
    ax2.set_ylabel('Number of Establishments')
    ax2.set_xlabel('Year')
    ax2.legend()
    ax2.grid(axis='y', linestyle='--', alpha=0.6)
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    fig.suptitle('Manufacturing Establishment Dynamics', fontsize=18, y=0.96)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.show()

