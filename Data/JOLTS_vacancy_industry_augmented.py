"""
ANALYSIS: Vacancy-Business Dynamic Analysis
============================================

Investigates the relationship between vacancies, establishment births, deaths, and total
establishment counts across industries.

Core Hypotheses:
1. Higher establishment numbers (births) → faster vacancy creation
2. More firm exits (deaths) → slower vacancy creation
3. Total establishment count captures the extensive margin of the labor market

Data Sources:
- JOLTS (Job Openings and Labor Turnover Survey): Vacancy data from FRED
- BED (Business Employment Dynamics): Birth/death/total establishment data (local Excel file)

Outputs:
- Merged panel dataset with vacancies, births, deaths, and total establishments (QCEW-anchored) by industry/quarter
- Time series visualizations in log levels and cyclical components
- Comovement statistics (correlations, volatility ratios)
- LaTeX and markdown tables for publication
"""

# =============================================================================
# 1. IMPORTS AND CONFIGURATION
# =============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from fredapi import Fred
from tabulate import tabulate
import os

# Display settings
pd.set_option('display.max_columns', 8)


# --- Working Directory ---
os.chdir(r"C:\Users\msilv\Documents\GitHub\Sunk-entry-costs--endogenous-variety--and-unemployment\Data")
os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment")
# --- FRED API Configuration ---
FRED_API_KEY = '9c70445138df124be4928605b7e08bd4'
FORCE_JOLTS_RELOAD = False  # Set to True to re-download JOLTS data from FRED

# --- File Paths ---
JOLTS_DATA_CSV = 'jolts_vacancies_industry_adjusted_levels.csv'
BED_DATA_XLSX = "BED_establishment_industry.xlsx"
QCEW_ESTABS_CSV = 'qcew_establishments_industry.csv'
FORCE_QCEW_RELOAD = False  # Set to True to re-download QCEW data from BLS

# --- HP Filter Configuration for Cyclical Analysis ---
HP_LAMBDA = 10**5

# --- Industry Definitions ---
SELECTED_INDUSTRIES = [
    "Nondurable Goods Manufacturing",
    "Durable Goods Manufacturing",
    "Construction",
    "Trade, Transportation, and Utilities",
    "Professional and Business Services",
    "Leisure and Hospitality"
]

# --- BED Data Series ID to Internal Code Mapping ---
BED_ID_TO_NAME_MAP = {
    'BDS0000000000200020120007LQ5': 'BTH_NDR',  # Nondurable births
    'BDS0000000000200010120007LQ5': 'BTH_DUR',  # Durable births
    'BDS0000000000100000120007LQ5': 'BTH_CON',  # Construction births
    'BDS0000000000200070120007LQ5': 'BTH_TTU',  # Trade/Trans/Utils births
    'BDS0000000000200050120007LQ5': 'BTH_LHS',  # Leisure/Hospitality births
    'BDS0000000000200100120007LQ5': 'BTH_PBS',  # Professional and Business Services births
    'BDS0000000000200020120008LQ5': 'DTH_NDR',  # Nondurable deaths
    'BDS0000000000200010120008LQ5': 'DTH_DUR',  # Durable deaths
    'BDS0000000000100000120008LQ5': 'DTH_CON',  # Construction deaths
    'BDS0000000000200070120008LQ5': 'DTH_TTU',  # Trade/Trans/Utils deaths
    'BDS0000000000200050120008LQ5': 'DTH_LHS',  # Leisure/Hospitality deaths
    'BDS0000000000200100120008LQ5': 'DTH_PBS',  # Professional and Business Services deaths
}
# --- Industry Code to Full Name Mapping (for merging) ---
BED_CODE_TO_INDUSTRY_MAP = {
    'NDR': 'Nondurable Goods Manufacturing',
    'DUR': 'Durable Goods Manufacturing',
    'CON': 'Construction',
    'TTU': 'Trade, Transportation, and Utilities',
    'PBS': 'Professional and Business Services',
    'LHS': 'Leisure and Hospitality',
}


# --- QCEW Supersector/NAICS Code Map ---
# Maps QCEW industry codes to BED industry codes for merging.
# BED "Trade, Transportation, and Utilities" = QCEW supersector 1021
# BED "Professional and Business Services"   = QCEW supersector 1024
# BED "Leisure and Hospitality"              = QCEW supersector 1026
# BED "Construction"                         = QCEW supersector 1012
# BED "Nondurable Goods Manufacturing"       = NAICS 31-32 (sectors 31, 32)
# BED "Durable Goods Manufacturing"          = NAICS 33
QCEW_INDUSTRY_MAP = {
    'NDR': '31_32',   # Nondurable Goods Manufacturing (NAICS 31-32)
    'DUR': '33',      # Durable Goods Manufacturing (NAICS 33)
    'CON': '23',      # Construction (NAICS 23)
    'TTU': '1021',    # Trade, Transportation, and Utilities (supersector)
    'PBS': '1024',    # Professional and Business Services (supersector)
    'LHS': '1026',    # Leisure and Hospitality (supersector)
}

# =============================================================================
# 2. DATA LOADING AND PROCESSING FUNCTIONS
# =============================================================================

def load_jolts_data(api_key, selected_industries, filepath, force_reload=False):
    """
    Loads JOLTS vacancy data from FRED API or local CSV cache.
    
    Processes the data into a clean quarterly panel format for selected industries.
    
    Parameters
    ----------
    api_key : str
        FRED API key for data access
    selected_industries : list
        List of industry names to extract
    filepath : str
        Path to save/load cached JOLTS data
    force_reload : bool, optional
        If True, re-download from FRED instead of using cached file (default: False)
    
    Returns
    -------
    pd.DataFrame
        Long-format panel with columns: date, industry, vacancies
    """
    print("--- Loading JOLTS Vacancy Data ---")
    
    if force_reload or not os.path.exists(filepath):
        print("Fetching data from FRED API...")
        fred = Fred(api_key=api_key)
        
        # FRED series IDs for industry vacancy levels
        industries_level_full = {
            'JTS3200JOL': 'Durable Goods Manufacturing',
            'JTS3400JOL': 'Nondurable Goods Manufacturing',
            'JTS2300JOL': 'Construction',
            'JTS4000JOL': 'Trade, Transportation, and Utilities',
            'JTS540099JOL': 'Professional and Business Services',
            'JTS7000JOL': 'Leisure and Hospitality'
        }
        
        data_list = [
            fred.get_series(sid, name=name)
            for sid, name in industries_level_full.items()
            if name in selected_industries
        ]
        df_wide = pd.concat(data_list, axis=1)
        df_wide.to_csv(filepath)
    else:
        print(f"Loading data from local file: {filepath}")
        df_wide = pd.read_csv(filepath, index_col=0, parse_dates=True)

    # Clean data and convert to quarterly panel
    df_wide.index = pd.to_datetime(df_wide.index)
    df_quarterly = df_wide[selected_industries].resample("QE").mean()
    df_panel = df_quarterly.reset_index().melt(
        id_vars=["index"],
        var_name="industry",
        value_name="vacancies"
    )
    df_panel.rename(columns={'index': 'date'}, inplace=True)
    
    print("JOLTS data processing complete.\n")
    return df_panel



def load_qcew_data(industry_map, filepath, force_reload=False):
    """
    Loads QCEW quarterly establishment counts by industry from BLS open data API.
    
    Fetches the `qtrly_estabs` field (count of establishments) for each industry
    using the QCEW CSV slice API: https://data.bls.gov/cew/data/api/{year}/{qtr}/industry/{code}.csv
    
    Filters to: area_fips='US000' (national), own_code='5' (private sector),
    size_code='0' (all sizes). For multi-sector groupings (TTU, PBS, LHS), sums
    constituent supersector counts. For manufacturing sub-sectors (NDR, DUR),
    sums NAICS 31+32 or uses 33 respectively.
    
    Units: actual establishment counts (not thousands).
    
    Parameters
    ----------
    industry_map : dict
        Mapping from BED industry code (e.g. 'NDR') to QCEW industry/supersector code
    filepath : str
        Path to save/load cached QCEW data as CSV
    force_reload : bool
        If True, re-fetch from BLS API regardless of cache
    
    Returns
    -------
    pd.DataFrame
        Long-format panel with columns: date, industry_code, qcew_estabs
        Covers 1993Q1 onwards (QCEW NAICS data starts 1990, BED overlap starts ~1993).
    """
    print("--- Loading QCEW Establishment Counts ---")
    
    if not force_reload and os.path.exists(filepath):
        print(f"Loading QCEW data from cache: {filepath}")
        df = pd.read_csv(filepath, parse_dates=['date'])
        print("QCEW data loaded from cache.\n")
        return df
    
    print("Fetching QCEW data from BLS API...")
    import urllib.request
    import io
    
    # QCEW API parameters
    # own_code=5: private sector; agglvl_code=14 for supersectors, 13 for NAICS sectors
    # For national totals: area_fips='US000'
    # We need quarterly data from 1993 to present
    
    # Determine year range: BED starts Q2 1993, QCEW NAICS from 1990
    start_year = 1993
    current_year = 2025  # update as needed; script will handle missing future quarters
    
    records = []
    
    for ind_code, qcew_code in industry_map.items():
        print(f"  Fetching {ind_code} (QCEW code: {qcew_code})...")
        
        # For NAICS 31_32 (nondurable), we need to sum 31 and 32 separately
        # because QCEW doesn't have a combined 31-32 supersector for establishments
        if qcew_code == '31_32':
            codes_to_sum = ['31', '32']
        else:
            codes_to_sum = [qcew_code]
        
        for year in range(start_year, current_year + 1):
            for qtr in range(1, 5):
                total_estabs = 0
                any_data = False
                
                for code in codes_to_sum:
                    url = f"https://data.bls.gov/cew/data/api/{year}/{qtr}/industry/{code}.csv"
                    try:
                        with urllib.request.urlopen(url, timeout=30) as resp:
                            raw = resp.read().decode('utf-8')
                        df_slice = pd.read_csv(io.StringIO(raw))
                        # Filter: national, private sector, all sizes
                        mask = (
                            (df_slice['area_fips'] == 'US000') &
                            (df_slice['own_code'] == 5) &
                            (df_slice['size_code'] == 0)
                        )
                        row = df_slice[mask]
                        if len(row) == 1:
                            total_estabs += int(row['qtrly_estabs'].iloc[0])
                            any_data = True
                    except Exception:
                        pass  # missing quarter (future data); skip
                
                if any_data:
                    # Convert quarter to end-of-quarter timestamp matching BED convention
                    period = pd.Period(f"{year}Q{qtr}")
                    date = period.to_timestamp('Q')
                    records.append({
                        'date': date,
                        'industry_code': ind_code,
                        'qcew_estabs': total_estabs
                    })
    
    df_qcew = pd.DataFrame(records).sort_values(['industry_code', 'date']).reset_index(drop=True)
    df_qcew.to_csv(filepath, index=False)
    print(f"QCEW data saved to {filepath}\n")
    return df_qcew


def load_bed_data(filepath, id_map, industry_map, qcew_panel):
    """
    Loads and processes Business Employment Dynamics (BED) data.
    
    Converts wide-format Excel data into clean quarterly panel format.
    Recreates date index from scratch based on BED data structure (start: Q3 1992).
    Constructs total establishment count as a cumulative index (base = 100 at first
    observation) from births and deaths: N_t = N_{t-1} + births_t - deaths_t.
    
    Parameters
    ----------
    filepath : str
        Path to BED Excel file with wide format (series as rows, dates as columns)
    id_map : dict
        Mapping from FRED series IDs to short codes (e.g., 'BTH_NDR')
    industry_map : dict
        Mapping from short codes to full industry names
    
    Parameters (additional)
    -----------------------
    qcew_panel : pd.DataFrame
        QCEW establishment counts with columns: date, industry_code, qcew_estabs.
        Used to anchor the level of the establishments series. Units: actual counts.
    
    Returns
    -------
    pd.DataFrame
        Long-format panel with columns: date, industry, births, deaths, establishments
        where establishments is in actual counts (same units as QCEW qtrly_estabs).
    """
    print("--- Loading BED Birth/Death Data ---")
    
    # Read wide-format Excel file
    df_raw = pd.read_excel(filepath, skiprows=3, index_col=0)

    # Filter for selected series
    df_filtered = df_raw.loc[id_map.keys()]

    # Transpose to get dates as rows
    df_transposed = df_filtered.transpose()
    
    # Clean data values: remove footnotes and convert to numeric
    for col in df_transposed.columns:
        df_transposed[col] = pd.to_numeric(
            df_transposed[col].astype(str).str.replace(r'[^\d.]', '', regex=True),
            errors='coerce'
        )
    
    # Create date index from scratch (Q3 1992 start for BED data)
    num_periods = len(df_transposed)
    date_index = pd.period_range(
        start='1992Q3',
        periods=num_periods,
        freq='Q'
    ).to_timestamp('Q')
    
    df_transposed.index = date_index
    
    # Rename columns from series IDs to short codes
    df_final = df_transposed.rename(columns=id_map)
    df_final = df_final.dropna(how='all').sort_index()
    
    # Reshape to long panel format
    df_long = df_final.reset_index().melt(
        id_vars='index',
        var_name='series_name',
        value_name='value'
    )
    df_long.rename(columns={'index': 'date'}, inplace=True)
    
    # Extract measure type (births vs deaths) and industry code
    df_long['measure'] = df_long['series_name'].str.split('_').str[0]
    df_long['industry_code'] = df_long['series_name'].str.split('_').str[1]
    
    # Pivot to separate births and deaths columns
    bed_panel_wide = df_long.pivot_table(
        index=['date', 'industry_code'],
        columns='measure',
        values='value'
    ).reset_index()
    bed_panel_wide.columns.name = None
    bed_panel_wide.rename(columns={'BTH': 'births', 'DTH': 'deaths'}, inplace=True)
    bed_panel_wide = bed_panel_wide.sort_values(['industry_code', 'date']).reset_index(drop=True)
    
    # Construct total establishments anchored to QCEW quarterly establishment counts.
    #
    # Method (backwards recursion):
    #   1. BED births/deaths are in thousands; multiply by 1000 to match QCEW units.
    #   2. Find the most recent BED date with BOTH valid births AND deaths.
    #   3. Find the closest QCEW observation at or before that date to use as anchor level.
    #   4. Work backwards: N_{t-1} = N_t - (births_t - deaths_t)
    #   5. Work forwards: N_{t+1} = N_t + (births_{t+1} - deaths_{t+1})
    #   6. Stops when hitting NaN births or deaths.
    #
    # Units of output 'establishments': actual count (same as QCEW qtrly_estabs).
    
    # Convert BED births and deaths from thousands to actual counts
    bed_panel_wide['births_actual'] = bed_panel_wide['births'] * 1000
    bed_panel_wide['deaths_actual'] = bed_panel_wide['deaths'] * 1000
    
    results = []
    for code, grp in bed_panel_wide.groupby('industry_code'):
        grp = grp.sort_values('date').reset_index(drop=True).copy()
        
        # Get QCEW data for this industry code
        qcew_ind = qcew_panel[qcew_panel['industry_code'] == code].set_index('date')['qcew_estabs']
        
        # Find most recent BED date with BOTH valid births AND deaths
        both_valid = grp[grp['births_actual'].notna() & grp['deaths_actual'].notna()]
        if len(both_valid) == 0:
            grp['establishments'] = np.nan
            results.append(grp)
            continue
        
        last_valid_bed_date = both_valid['date'].max()
        print(f"DEBUG [{code}]: last_valid_bed_date = {last_valid_bed_date}, QCEW dates available: {len(qcew_ind) if len(qcew_ind) > 0 else 0}")
        
        # Find closest QCEW observation at or before that BED date
        if len(qcew_ind) == 0 or qcew_ind.isna().all():
            grp['establishments'] = np.nan
            results.append(grp)
            continue
        
        qcew_at_or_before = qcew_ind[qcew_ind.index <= last_valid_bed_date].dropna()
        if len(qcew_at_or_before) == 0:
            grp['establishments'] = np.nan
            results.append(grp)
            continue
        
        anchor_qcew_date = qcew_at_or_before.index.max()
        anchor_qcew_level = float(qcew_at_or_before[anchor_qcew_date])
        
        # Find position of last_valid_bed_date in grp
        anchor_pos = grp[grp['date'] == last_valid_bed_date].index[0]
        
        est = pd.Series(np.nan, index=grp.index)
        
        # Initialize at anchor
        est.iloc[anchor_pos] = anchor_qcew_level
        print(f"DEBUG [{code}]: anchor_pos = {anchor_pos}, anchor_qcew_level = {anchor_qcew_level}")
        
        # Work backwards from anchor to beginning
        back_count = 0
        for i in range(anchor_pos - 1, -1, -1):
            net_flow = grp.loc[i + 1, 'births_actual'] - grp.loc[i + 1, 'deaths_actual']
            if pd.notna(net_flow):
                est.iloc[i] = est.iloc[i + 1] - net_flow
                back_count += 1
            else:
                print(f"DEBUG [{code}]: backward loop stopped at i={i} (NaN births/deaths at i+1)")
                break  # Stop if we hit missing births/deaths
        
        print(f"DEBUG [{code}]: backward loop filled {back_count} rows")
        
        # Work forwards from anchor to end
        fwd_count = 0
        for i in range(anchor_pos + 1, len(grp)):
            net_flow = grp.loc[i, 'births_actual'] - grp.loc[i, 'deaths_actual']
            if pd.notna(net_flow):
                est.iloc[i] = est.iloc[i - 1] + net_flow
                fwd_count += 1
            else:
                print(f"DEBUG [{code}]: forward loop stopped at i={i} (NaN births/deaths)")
                break  # Stop if we hit missing births/deaths
        
        print(f"DEBUG [{code}]: forward loop filled {fwd_count} rows")
        print(f"DEBUG [{code}]: total non-null establishments = {est.notna().sum()}")
        print()
        
        grp['establishments'] = est
        results.append(grp)
    
    bed_panel_wide = pd.concat(results, ignore_index=True)
    # Drop intermediate columns
    bed_panel_wide.drop(columns=['births_actual', 'deaths_actual'], inplace=True)
    
    # Map industry codes to full names for merging
    bed_panel = bed_panel_wide.copy()
    bed_panel['industry'] = bed_panel['industry_code'].map(industry_map)
    bed_panel.drop(columns=['industry_code'], inplace=True)
    
    
    print("BED data processing complete.\n")
    return bed_panel


# =============================================================================
# 3. DATA LOADING, MERGING, AND VALIDATION
# =============================================================================

# Load datasets
jolts_panel = load_jolts_data(
    FRED_API_KEY,
    SELECTED_INDUSTRIES,
    JOLTS_DATA_CSV,
    force_reload=FORCE_JOLTS_RELOAD
)
qcew_panel = load_qcew_data(
    QCEW_INDUSTRY_MAP,
    QCEW_ESTABS_CSV,
    force_reload=FORCE_QCEW_RELOAD
)
bed_panel = load_bed_data(BED_DATA_XLSX, BED_ID_TO_NAME_MAP, BED_CODE_TO_INDUSTRY_MAP, qcew_panel)

# Merge JOLTS and BED panels
print("--- Merging JOLTS and BED data ---")
merged_panel = pd.merge(jolts_panel, bed_panel, on=['date', 'industry'], how='inner')
merged_panel.set_index(['industry', 'date'], inplace=True)
merged_panel.sort_index(inplace=True)

# Validation
print("\n=== MERGED DATAFRAME SUMMARY ===")
print(f"Shape: {merged_panel.shape}")
print(f"\nFirst 10 rows:\n{merged_panel.head(10)}")
print(f"\nLast 10 rows:\n{merged_panel.tail(10)}")
print(f"\nData Info:")
merged_panel.info()

# =============================================================================
# 4. LOG TRANSFORMATION AND DETRENDING
# =============================================================================

# Transform to log levels and demean by industry
print("\n--- Transforming to Log Levels ---")
cols_to_transform = ['vacancies', 'births', 'deaths', 'establishments']
log_df = merged_panel[cols_to_transform].apply(np.log)
log_df = log_df.groupby('industry').transform(lambda x: x - x.iloc[0])

industries = log_df.index.get_level_values('industry').unique()


# =============================================================================
# 5. CYCLICAL COMPONENT EXTRACTION (HP Filter)
# =============================================================================

def get_hp_cycle(series):
    """
    Applies HP filter to extract cyclical component.
    
    Parameters
    ----------
    series : pd.Series
        Time series to filter (should be in log levels or growth rates)
    
    Returns
    -------
    np.ndarray
        Cyclical component from HP filter
    """
    s = series.dropna()
    # HP filter requires at least 3 observations; use a higher threshold to avoid
    # numerical issues on very short series. If series is too short, return an
    # array of NaNs with the same index.
    if len(s) < 8:
        return pd.Series(data=np.nan, index=series.index)
    cycle, trend = sm.tsa.filters.hpfilter(s, HP_LAMBDA)
    # Reindex cycle to match original series index (fill missing with NaN)
    cycle = pd.Series(cycle, index=s.index).reindex(series.index)
    return cycle


def apply_filter_to_group(df_chunk):
    """Applies HP filter to each column within a grouped dataset."""
    return df_chunk.apply(get_hp_cycle)


print("--- Extracting Cyclical Components (HP Filter) ---")
log_df_cleaned = log_df.dropna(subset=cols_to_transform, how="all")
cycle = log_df_cleaned.groupby('industry').apply(apply_filter_to_group)
cycle.index = cycle.index.droplevel(1)


# =============================================================================
# 6. VISUALIZATION AND TIME SERIES PLOTTING
# =============================================================================

def create_industry_plots(data, industries, y_key, y_secondary, title_suffix):
    """
    Creates a 2x3 grid of time series plots (one per industry).
    
    REFACTORED to eliminate redundancy between vacancies vs births/deaths plotting.
    
    Parameters
    ----------
    data : pd.DataFrame
        Data with multi-index (industry, date)
    industries : pd.Index
        Unique industry names
    y_key : str
        Primary y-axis variable (e.g., 'vacancies')
    y_secondary : str
        Secondary variable to plot against (e.g., 'births' or 'deaths')
    title_suffix : str
        Suffix for the plot title
    
    Returns
    -------
    tuple of (fig, axes)
    """
    color_map = {'births': 'forestgreen', 'deaths': 'red', 'establishments': 'darkorange'}
    color_secondary = color_map.get(y_secondary, 'black')
    label_secondary = f'Log {y_secondary.capitalize()}'
    
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 8), sharex=True)
    axes = axes.ravel()
    
    for i, industry in enumerate(industries):
        ax = axes[i]
        industry_data = data.loc[industry]
        
        ax.plot(
            industry_data.index,
            industry_data[y_key],
            color='blue',
            label=f'Log {y_key.capitalize()}',
            linewidth=1.5
        )
        ax.plot(
            industry_data.index,
            industry_data[y_secondary],
            color=color_secondary,
            label=label_secondary,
            linewidth=1.5,
            alpha=0.9
        )
        ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        ax.set_title(industry, fontsize=14)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    axes[0].legend()
    fig.suptitle(f'Log Vacancies vs. Log {y_secondary.capitalize()} by Industry - {title_suffix}',
                 fontsize=16, y=1.02)
    plt.tight_layout()
    return fig, axes


# Plot 1: Log Levels - Vacancies vs. Births
print("\n--- Generating Time Series Plots (Log Levels) ---")
fig1, _ = create_industry_plots(log_df, industries, 'vacancies', 'births', 'Log Levels')
plt.show()

# Plot 2: Log Levels - Vacancies vs. Deaths
fig2, _ = create_industry_plots(log_df, industries, 'vacancies', 'deaths', 'Log Levels')
plt.show()

# Plot 3: Cyclical Components - Vacancies vs. Deaths
print("\n--- Generating Time Series Plots (Cyclical Components) ---")
fig3, _ = create_industry_plots(cycle, industries, 'vacancies', 'deaths', 'Cyclical Components')
plt.show()

# Plot 4: Log Levels - Vacancies vs. Total Establishments
print("\n--- Generating Time Series Plots (Log Levels, Establishments) ---")
fig4, _ = create_industry_plots(log_df, industries, 'vacancies', 'establishments', 'Log Levels')
plt.show()

# Plot 5: Cyclical Components - Vacancies vs. Total Establishments
fig5, _ = create_industry_plots(cycle, industries, 'vacancies', 'establishments', 'Cyclical Components')
plt.show()


# =============================================================================
# 7. COMOVEMENT STATISTICS AND CORRELATION ANALYSIS
# =============================================================================

def calculate_comovement_stats(group, series_name='deaths'):
    """
    Calculates key comovement statistics for vacancies and a secondary series (births or deaths).
    
    Computes contemporaneous correlations, lead/lag correlations, and relative volatilities
    to assess the dynamic relationship between vacancy openings and firm dynamics.
    
    Parameters
    ----------
    group : pd.DataFrame
        Cyclical components data for a single industry
    series_name : str, optional
        Name of the series to correlate with vacancies ('deaths' or 'births'), default 'deaths'
    
    Returns
    -------
    pd.Series
        Dictionary of statistics with keys:
        - corr(v_t, series_t): Contemporaneous correlation
        - corr(v_t, series_t-1): Lagged correlation (series leads)
        - corr(v_t, series_t+1): Leading correlation (series follows)
        - std(series) / std(vacancies): Relative volatility
        - std(vacancies): Vacancy volatility
        - std(series): Series volatility
    """
    if 'vacancies' not in group or series_name not in group:
        return pd.Series(dtype='float64')
    
    vacancies = group['vacancies']
    series = group[series_name]
    
    # Correlations
    corr_contemp = vacancies.corr(series)
    corr_lag1 = vacancies.corr(series.shift(1))
    corr_lead1 = vacancies.corr(series.shift(-1))
    
    # Volatility
    std_vac = vacancies.std()
    std_series = series.std()
    relative_vol = std_series / std_vac if std_vac > 0 else np.nan
    
    return pd.Series({
        f'corr(v_t, {series_name}_t)': corr_contemp,
        f'corr(v_t, {series_name}_t-1)': corr_lag1,
        f'corr(v_t, {series_name}_t+1)': corr_lead1,
        f'std({series_name}) / std(vacancies)': relative_vol,
        'std(vacancies)': std_vac,
        f'std({series_name})': std_series,
    })


# Apply function to each industry - Deaths, Births, and Total Establishments
print("\n--- Computing Industry Statistics ---")
comov_table_deaths = cycle.groupby('industry').apply(lambda x: calculate_comovement_stats(x, 'deaths'))
comov_table_births = cycle.groupby('industry').apply(lambda x: calculate_comovement_stats(x, 'births'))
comov_table_establishments = cycle.groupby('industry').apply(lambda x: calculate_comovement_stats(x, 'establishments'))

# Display in console
print("\nComovement Statistics: Cyclical Vacancies vs. Cyclical Deaths")
print(comov_table_deaths.to_string())

print("\n\nComovement Statistics: Cyclical Vacancies vs. Cyclical Births")
print(comov_table_births.to_string())

print("\n\nComovement Statistics: Cyclical Vacancies vs. Cyclical Total Establishments")
print(comov_table_establishments.to_string())



# =============================================================================
# 8. TABLE GENERATION AND EXPORT
# =============================================================================

# Helper function for formatting
def sig_fig_formatter(x):
    """Formats a number to 3 significant figures for publication."""
    return f"{x:.3g}"


# =============================================================================
# TABLE 8a: Vacancies vs. Deaths (Markdown)
# =============================================================================
print("\n--- TABLE 8a: Vacancies vs. Deaths (Markdown Format) ---")
comov_table_deaths_rounded = comov_table_deaths.round(3)
markdown_table_deaths = tabulate(
    comov_table_deaths_rounded,
    headers="keys",
    tablefmt="pipe",
    showindex=True
)
print(markdown_table_deaths)

# =============================================================================
# TABLE 8b: Vacancies vs. Births (Markdown)
# =============================================================================
print("\n--- TABLE 8b: Vacancies vs. Births (Markdown Format) ---")
comov_table_births_rounded = comov_table_births.round(3)
markdown_table_births = tabulate(
    comov_table_births_rounded,
    headers="keys",
    tablefmt="pipe",
    showindex=True
)
print(markdown_table_births)


# =============================================================================
# TABLE 8a: Vacancies vs. Deaths (LaTeX)
# =============================================================================
print("\n--- TABLE 8a: Vacancies vs. Deaths (LaTeX Format) ---")
latex_table_deaths = comov_table_deaths.to_latex(
    float_format=sig_fig_formatter,
    caption="Comovement Statistics for Cyclical Vacancies and Establishment Deaths",
    label="tab:comovement_deaths",
    header=True,
    index=True,
    column_format='l' + 'r' * len(comov_table_deaths.columns)
)
print(latex_table_deaths)

# =============================================================================
# TABLE 8b: Vacancies vs. Births (LaTeX)
# =============================================================================
print("\n--- TABLE 8b: Vacancies vs. Births (LaTeX Format) ---")
latex_table_births = comov_table_births.to_latex(
    float_format=sig_fig_formatter,
    caption="Comovement Statistics for Cyclical Vacancies and Establishment Births",
    label="tab:comovement_births",
    header=True,
    index=True,
    column_format='l' + 'r' * len(comov_table_births.columns)
)
print(latex_table_births)

# =============================================================================
# TABLE 8c: Vacancies vs. Total Establishments (Markdown)
# =============================================================================
print("\n--- TABLE 8c: Vacancies vs. Total Establishments (Markdown Format) ---")
comov_table_establishments_rounded = comov_table_establishments.round(3)
markdown_table_establishments = tabulate(
    comov_table_establishments_rounded,
    headers="keys",
    tablefmt="pipe",
    showindex=True
)
print(markdown_table_establishments)

# =============================================================================
# TABLE 8c: Vacancies vs. Total Establishments (LaTeX)
# =============================================================================
print("\n--- TABLE 8c: Vacancies vs. Total Establishments (LaTeX Format) ---")
latex_table_establishments = comov_table_establishments.to_latex(
    float_format=sig_fig_formatter,
    caption="Comovement Statistics for Cyclical Vacancies and Total Establishment Count",
    label="tab:comovement_establishments",
    header=True,
    index=True,
    column_format='l' + 'r' * len(comov_table_establishments.columns)
)
print(latex_table_establishments)


