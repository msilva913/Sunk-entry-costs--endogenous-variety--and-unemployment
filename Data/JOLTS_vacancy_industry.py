
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

"""
This script loads, cleans, and merges JOLTS vacancy data from FRED
and Business Employment Dynamics (BED) data from a local file.
The goal is to create a single panel DataFrame for analyzing the 
relationship between vacancies, establishment births, and deaths across key industries.

Analysis Hypotheses:
1) Higher establishment numbers (births) lead to faster vacancy creation.
2) More firm exits (deaths) lead to slower vacancy creation.
"""

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================
import pandas as pd
import numpy as np
from fredapi import Fred
import os
from time_series_functions import hp_filter
# --- Configuration ---
# Set your working directory if needed
# os.chdir(r"C:\Your\Project\Path\Data")

# API and File Paths
FRED_API_KEY = '9c70445138df124be4928605b7e08bd4'
JOLTS_DATA_CSV = 'jolts_vacancies_industry_adjusted_levels.csv'
BED_DATA_XLSX = "BED_establishment_industry.xlsx"
FORCE_JOLTS_RELOAD = False # Set to True to re-download JOLTS data from FRED

# Define the 6 industries of interest
SELECTED_INDUSTRIES = [
    "Nondurable Goods Manufacturing", "Durable Goods Manufacturing", "Construction",
    "Trade, Transportation, and Utilities", "Professional and Business Services", "Leisure and Hospitality"
]

# Mapping for BED Series IDs to internal codes
BED_ID_TO_NAME_MAP = {
    'BDS0000000000200020120007LQ5': 'BTH_NDR', 'BDS0000000000200010120007LQ5': 'BTH_DUR',
    'BDS0000000000100000120007LQ5': 'BTH_CON', 'BDS0000000000200070120007LQ5': 'BTH_TTU',
    #'BDS0000000000600000120007LQ5': 'BTH_PBS',
    'BDS0000000000200050120007LQ5': 'BTH_LHS',
    'BDS0000000000200020120008LQ5': 'DTH_NDR', 'BDS0000000000200010120008LQ5': 'DTH_DUR',
    'BDS0000000000100000120008LQ5': 'DTH_CON', 'BDS0000000000200070120008LQ5': 'DTH_TTU',
  #  'BDS0000000000600000120008LQ5': 'DTH_PBS', 
  'BDS0000000000200050120008LQ5': 'DTH_LHS',
}

# Mapping from internal codes to the full industry names used in the JOLTS data
BED_CODE_TO_INDUSTRY_MAP = {
    'NDR': 'Nondurable Goods Manufacturing', 'DUR': 'Durable Goods Manufacturing', 'CON': 'Construction',
    'TTU': 'Trade, Transportation, and Utilities', 'PBS': 'Professional and Business Services', 'LHS': 'Leisure and Hospitality',
}


# =============================================================================
# 2. DATA LOADING AND PREPARATION FUNCTIONS
# =============================================================================

def load_jolts_data(api_key, selected_industries, filepath, force_reload=False):
    """
    Loads JOLTS vacancy data from FRED API or local CSV, processes it into
    a clean quarterly panel format for the selected industries.
    """
    print("--- Loading JOLTS Vacancy Data ---")
    if force_reload or not os.path.exists(filepath):
        print("Fetching data from FRED API...")
        fred = Fred(api_key=api_key)
        
        industries_level_full = {
            'JTS3200JOL': 'Durable Goods Manufacturing', 'JTS3400JOL': 'Nondurable Goods Manufacturing',
            'JTS2300JOL': 'Construction', 'JTS4000JOL': 'Trade, Transportation, and Utilities',
            'JTS540099JOL': 'Professional and Business Services', 'JTS7000JOL': 'Leisure and Hospitality'
        }
        
        data_list = [fred.get_series(sid, name=name) for sid, name in industries_level_full.items() if name in selected_industries]
        df_wide = pd.concat(data_list, axis=1)
        df_wide.to_csv(filepath)
    else:
        print(f"Loading data from local file: {filepath}")
        df_wide = pd.read_csv(filepath, index_col=0, parse_dates=True)

    df_wide.index = pd.to_datetime(df_wide.index)
    df_quarterly = df_wide[selected_industries].resample("QE").mean()
    df_panel = df_quarterly.reset_index().melt(id_vars=["index"], var_name="industry", value_name="vacancies")
    df_panel.rename(columns={'index': 'date'}, inplace=True)
    
    print("JOLTS data processing complete.\n")
    return df_panel

# --- THIS IS THE CORRECTED FUNCTION ---
def load_bed_data(filepath, id_map, industry_map):
    """
    Loads and processes Business Employment Dynamics (BED) data from a wide-format
    Excel file into a clean quarterly panel format by creating the date index from scratch.
    """
    print("--- Loading BED Birth/Death Data ---")
    # 1. Read wide-format Excel file, using the first column (Series ID) as the index
    df_raw = pd.read_excel(filepath, skiprows=3, index_col=0)

    # 2. Filter for the series we care about
    df_filtered = df_raw.loc[id_map.keys()]

    # 3. Transpose so dates are rows and series are columns
    df_transposed = df_filtered.transpose()
    
    # 4. Clean data values (remove footnotes, etc.) and convert to numeric
    for col in df_transposed.columns:
        df_transposed[col] = pd.to_numeric(
            df_transposed[col].astype(str).str.replace(r'[^\d.]', '', regex=True),
            errors='coerce'
        )

    # 5. --- CREATE DATE INDEX FROM SCRATCH ---
    # Generate a new DatetimeIndex starting from Q3 1992, matching the BED data's start.
    # The number of periods is determined by the number of rows in the transposed data.
    num_periods = len(df_transposed)
    date_index = pd.period_range(start='1992Q3', periods=num_periods, freq='Q').to_timestamp('Q')
    
    # Assign the newly created index to the DataFrame
    df_transposed.index = date_index

    # 6. Rename columns from Series IDs to our short names (e.g., 'BTH_NDR')
    df_final = df_transposed.rename(columns=id_map)
    df_final = df_final.dropna(how='all').sort_index()

    # 7. Reshape the data into a long panel format
    df_long = df_final.reset_index().melt(id_vars='index', var_name='series_name', value_name='value')
    df_long.rename(columns={'index': 'date'}, inplace=True)
    df_long['measure'] = df_long['series_name'].str.split('_').str[0]
    df_long['industry_code'] = df_long['series_name'].str.split('_').str[1]
    
    # 8. Pivot to get 'births' and 'deaths' as separate columns
    bed_panel_wide = df_long.pivot_table(
        index=['date', 'industry_code'],
        columns='measure',
        values='value'
    ).reset_index()
    bed_panel_wide.rename(columns={'BTH': 'births', 'DTH': 'deaths'}, inplace=True)
    
    # 9. Map industry codes to full names for merging and final cleanup
    bed_panel = bed_panel_wide.copy()
    bed_panel['industry'] = bed_panel['industry_code'].map(industry_map)
    bed_panel.drop(columns=['industry_code'], inplace=True)
    
    print("BED data processing complete.\n")
    return bed_panel


# =============================================================================
# 3. EXECUTION: LOAD, MERGE, AND VERIFY
# =============================================================================

# Load both datasets into panel format
jolts_panel = load_jolts_data(FRED_API_KEY, SELECTED_INDUSTRIES, JOLTS_DATA_CSV, force_reload=FORCE_JOLTS_RELOAD)
bed_panel = load_bed_data(BED_DATA_XLSX, BED_ID_TO_NAME_MAP, BED_CODE_TO_INDUSTRY_MAP)

# Merge the JOLTS and BED panels on 'date' and 'industry'
print("--- Merging JOLTS and BED data ---")
# Use an inner merge to keep only the time periods and industries present in both datasets
merged_panel = pd.merge(jolts_panel, bed_panel, on=['date', 'industry'], how='inner')

# Set a multi-index for easy analysis and sort
merged_panel.set_index(['industry', 'date'], inplace=True)
merged_panel.sort_index(inplace=True)

# Final verification
print("\n=== FINAL MERGED DATAFRAME ===")
print("\nShape:", merged_panel.shape)
print("\nData Head:")
print(merged_panel.head(10))
print("\n...")
print("Data Tail:")
print(merged_panel.tail(10))
print("\nData Info:")
merged_panel.info()
##################

# =============================================================================
# 4. VISUALIZATION: PLOT TIME SERIES IN LOG LEVELS
# =============================================================================
print("\n--- Generating Time Series Plots ---")
import matplotlib.pyplot as plt
# =============================================================================
# --- Use .apply() for efficient log transformation ---
#  1. Select the columns to transform
cols_to_transform = ['vacancies', 'births', 'deaths']

# 2. Apply the np.log function and add a 'log_' prefix to the new columns
log_df = merged_panel[cols_to_transform].apply(np.log)
log_df = log_df[cols_to_transform].groupby('industry').transform(lambda x: x-x.iloc[0])
# Get the unique list of industries to iterate over
industries = log_df.index.get_level_values('industry').unique()

# --- Plot 1: Vacancies vs. Establishment Births ---

# Create a 2x3 grid of subplots
fig1, axes1 = plt.subplots(nrows=2, ncols=3, figsize=(18, 8), sharex=True)
axes1 = axes1.ravel() # Flatten the 2x3 grid for easy iteration

for i, industry in enumerate(industries):
    ax = axes1[i]
    data_for_industry = log_df.loc[industry]
    
    # Plot using the new 'log_' columns
    ax.plot(data_for_industry.index, data_for_industry['vacancies'], color='blue', label='Log Vacancies', linewidth=1.5)
    ax.plot(data_for_industry.index, data_for_industry['births'], color='forestgreen', label='Log Births', linewidth=1.5, alpha=0.9)
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_title(industry, fontsize=14)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

axes1[0].legend()
fig1.suptitle('Log Vacancies vs. Log Establishment Births by Industry', fontsize=16, y=1.02)
plt.tight_layout()
plt.show()


# --- Plot 2: Vacancies vs. Establishment Deaths ---

# Create another 2x3 grid of subplots
fig2, axes2 = plt.subplots(nrows=2, ncols=3, figsize=(18, 8), sharex=True)
axes2 = axes2.ravel() # Flatten the grid

for i, industry in enumerate(industries):
    ax = axes2[i]
    data_for_industry = log_df.loc[industry]
    
    # Plot using the new 'log_' columns
    ax.plot(data_for_industry.index, data_for_industry['vacancies'], color='blue', label='Log Vacancies', linewidth=1.5)
    ax.plot(data_for_industry.index, data_for_industry['deaths'], color='red', label='Log Deaths', linewidth=1.5, alpha=0.8)
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    
    ax.set_title(industry, fontsize=14)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

axes2[0].legend()
fig2.suptitle('Log Vacancies vs. Log Establishment Deaths by Industry', fontsize=16, y=1.02)
plt.tight_layout()
plt.show()

" Extract cycle "
import statsmodels.api as sm
hp_lambda=10**5
def get_hp_cycle(series):
    """Applies the HP filter to a single series and returns the cyclical component."""
    cycle, trend = sm.tsa.filters.hpfilter(series.dropna(), hp_lambda)
    return cycle

def apply_filter_to_group(df_chunk):
    return df_chunk.apply(get_hp_cycle)

cols_to_filter = ['vacancies', 'births', 'deaths']
log_df_cleaned = log_df.dropna(subset=cols_to_filter, how="all")
cycle = log_df_cleaned.groupby('industry').apply(apply_filter_to_group)
cycle.index = cycle.index.droplevel(1)

# Create another 2x3 grid of subplots
fig4, axes4 = plt.subplots(nrows=2, ncols=3, figsize=(18, 8), sharex=True)
axes4 = axes4.ravel() # Flatten the grid

for i, industry in enumerate(industries):
    ax = axes4[i]
    data_for_industry = cycle.loc[industry]
    
    # Plot using the new 'log_' columns
    ax.plot(data_for_industry.index, data_for_industry['vacancies'], color='blue', label='Cyclical Vacancies', linewidth=1.5)
    ax.plot(data_for_industry.index, data_for_industry['deaths'], color='red', label='Cyclical Deaths', linewidth=1.5, alpha=0.8)
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    
    ax.set_title(industry, fontsize=14)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

axes4[0].legend()
fig4.suptitle('Log Vacancies vs. Log Establishment Deaths by Industry', fontsize=16, y=1.02)
plt.tight_layout()
plt.show()

# 3) Moments 
#mom = moments(cycle, relative_std="Nondurable Goods Manufacturing", lab=["Nondurable Goods Manufacturing"])