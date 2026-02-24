"""
ANALYSIS: Vacancy-Business Dynamic Analysis
============================================

Investigates the relationship between vacancies, establishment births, and deaths across industries.

Core Hypotheses:
1. Higher establishment numbers (births) → faster vacancy creation
2. More firm exits (deaths) → slower vacancy creation

Data Sources:
- JOLTS (Job Openings and Labor Turnover Survey): Vacancy data from FRED
- BED (Business Employment Dynamics): Birth/death establishment data (local Excel file)

Outputs:
- Merged panel dataset with vacancies, births, and deaths by industry/quarter
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
import statsmodels.formula.api as smf
from fredapi import Fred
from tabulate import tabulate
import os

# Display settings
pd.set_option('display.max_columns', 8)


# --- Working Directory ---
os.chdir(r"C:\Users\msilva913\Documents\GitHub\Sunk_entry_costs_endogenous_variety_unemployment\Data")

# --- FRED API Configuration ---
FRED_API_KEY = '9c70445138df124be4928605b7e08bd4'
FORCE_JOLTS_RELOAD = False  # Set to True to re-download JOLTS data from FRED

# --- File Paths ---
JOLTS_DATA_CSV = 'jolts_vacancies_industry_adjusted_levels.csv'
BED_DATA_XLSX = "BED_establishment_industry.xlsx"

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
# --- QCEW Industry Code Map (corrected) ---
# For manufacturing, QCEW publishes NAICS sectors 31, 32, 33 separately.
# There is no combined nondurable (31+32) supersector code — must be summed manually.
QCEW_INDUSTRY_MAP = {
    'NDR': ['31', '32'],  # Nondurable: fetch NAICS 31 AND 32 separately, then sum
    'DUR': ['33'],         # Durable: NAICS sector 33
    'CON': ['23'],         # Construction: NAICS sector 23
    'TTU': ['1021'],       # Trade, Transportation, and Utilities (supersector)
    'PBS': ['1024'],       # Professional and Business Services (supersector)
    'LHS': ['1026'],       # Leisure and Hospitality (supersector)
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


def load_bed_data(filepath, id_map, industry_map):
    """
    Loads and processes Business Employment Dynamics (BED) data.
    
    Converts wide-format Excel data into clean quarterly panel format.
    Recreates date index from scratch based on BED data structure (start: Q3 1992).
    
    Parameters
    ----------
    filepath : str
        Path to BED Excel file with wide format (series as rows, dates as columns)
    id_map : dict
        Mapping from FRED series IDs to short codes (e.g., 'BTH_NDR')
    industry_map : dict
        Mapping from short codes to full industry names
    
    Returns
    -------
    pd.DataFrame
        Long-format panel with columns: date, industry, births, deaths
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
    bed_panel_wide.rename(columns={'BTH': 'births', 'DTH': 'deaths'}, inplace=True)
    
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
bed_panel = load_bed_data(BED_DATA_XLSX, BED_ID_TO_NAME_MAP, BED_CODE_TO_INDUSTRY_MAP)

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
cols_to_transform = ['vacancies', 'births', 'deaths']
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
    color_secondary = 'forestgreen' if y_secondary == 'births' else 'red'
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


# Apply function to each industry - Deaths and Births
print("\n--- Computing Industry Statistics ---")
comov_table_deaths = cycle.groupby('industry').apply(lambda x: calculate_comovement_stats(x, 'deaths'))
comov_table_births = cycle.groupby('industry').apply(lambda x: calculate_comovement_stats(x, 'births'))

# Display in console
print("\nComovement Statistics: Cyclical Vacancies vs. Cyclical Deaths")
print(comov_table_deaths.to_string())

print("\n\nComovement Statistics: Cyclical Vacancies vs. Cyclical Births")
print(comov_table_births.to_string())



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
# 9. PANEL REGRESSION: Vacancies on Lagged Deaths
# =============================================================================
"""
Panel Regression Analysis
==========================

Dependent variable : cyclical component of log(vacancies_{i,t})
Independent variable: cyclical component of log(deaths_{i,t-1})  [one-quarter lag]

Three specifications:
  (1) No controls           : vacancies_it = β·deaths_lag1_it + ε_it
  (2) Time fixed effects    : vacancies_it = β·deaths_lag1_it + γ_t + ε_it
  (3) Industry fixed effects: vacancies_it = β·deaths_lag1_it + α_i + ε_it

Standard errors clustered by industry in all specifications.

Note: Both industry and time FE together are not identified with only 6
industries and ~80 quarters (the within-group variation is exhausted). The
three separate specifications allow clear interpretation of each set of
controls.
"""

def _stars(p):
    """Return significance stars for a p-value."""
    if p < 0.01:  return '***'
    if p < 0.05:  return '**'
    if p < 0.10:  return '*'
    return ''


# --- Prepare regression data ---
print("\n--- Preparing Cyclical Panel for Regression ---")
cycle_panel = cycle.reset_index()

reg_df = cycle_panel.dropna(subset=['vacancies', 'deaths']).copy()
reg_df['deaths_lag1'] = reg_df.groupby('industry')['deaths'].shift(1)
reg_df = reg_df.dropna(subset=['deaths_lag1'])

# Categorical identifiers for fixed effects
reg_df['date_fe']     = reg_df['date'].astype(str)
reg_df['industry_fe'] = reg_df['industry'].astype(str)

print(f"Observations after lag: {len(reg_df)}")
print(f"Industries: {reg_df['industry'].nunique()}")
print(f"Quarters:   {reg_df['date'].nunique()}")

# --- Run three OLS models with clustered SEs ---
print("\n--- Running Panel OLS Regressions (3 specifications) ---")

specs = [
    ('(1) No Controls',          'vacancies ~ deaths_lag1'),
    ('(2) Time FE',               'vacancies ~ deaths_lag1 + C(date_fe)'),
    ('(3) Industry and Time FE',           'vacancies ~ deaths_lag1 + C(date_fe) + C(industry_fe)'),
]

results = {}
for label, formula in specs:
    ols_res  = smf.ols(formula, data=reg_df).fit()
    clustered = ols_res.get_robustcov_results(
        cov_type='cluster',
        groups=reg_df['industry']
    )

    coef_idx = ols_res.model.exog_names.index('deaths_lag1')

    results[label] = {
        'coef'  : float(np.asarray(clustered.params)[coef_idx]),
        'se'    : float(np.asarray(clustered.bse)[coef_idx]),
        'pval'  : float(np.asarray(clustered.pvalues)[coef_idx]),
        'nobs'  : int(ols_res.nobs),
        'r2'    : ols_res.rsquared,
    }

# --- Build summary table ---
print("\n" + "=" * 80)
print("TABLE 9: Effect of Lagged Establishment Deaths on Vacancies")
print("         Dependent variable: Cyclical log(vacancies)")
print("=" * 80)

summary_df = pd.DataFrame.from_dict(results, orient='index')
summary_df.index.name = 'Specification'
summary_df = summary_df.rename(columns={
    'coef': 'Coefficient',
    'se': 'Std. Error (clustered)',
    'pval': 'p-value',
    'r2': 'R²',
    'nobs': 'N'
})
summary_df['Coefficient'] = summary_df.apply(
    lambda r: f"{r['Coefficient']:.4f}{_stars(r['p-value'])}", axis=1
)
summary_df['Std. Error (clustered)'] = summary_df['Std. Error (clustered)'].map(lambda x: f"({x:.4f})")
summary_df['p-value'] = summary_df['p-value'].map(lambda x: f"{x:.4f}")
summary_df['R²'] = summary_df['R²'].map(lambda x: f"{x:.3f}")
summary_df['N'] = summary_df['N'].astype(int)
summary_df = summary_df[['Coefficient', 'Std. Error (clustered)', 'p-value', 'R²', 'N']]
print(summary_df.to_string())
print("\nSignificance: *** p<0.01  ** p<0.05  * p<0.10")
print("Standard errors clustered by industry.")
print("=" * 80)

# --- Markdown table ---
print("\n--- TABLE 9 (Markdown) ---")
print(tabulate(summary_df, headers='keys', tablefmt='pipe', showindex=True))

# --- LaTeX table ---
print("\n--- TABLE 9 (LaTeX) ---")
latex_reg = summary_df.to_latex(
    caption=(
        "Effect of Lagged Establishment Deaths on Vacancies. "
        "Dependent variable: cyclical component of log vacancies. "
        "Standard errors (in parentheses) clustered by industry. "
        "Significance: $^{***}$p$<$0.01, $^{**}$p$<$0.05, $^{*}$p$<$0.10."
    ),
    label='tab:panel_regression',
    column_format='l' + 'r' * len(summary_df.columns),
    escape=False,
)
print(latex_reg)

