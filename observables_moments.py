
import numpy as np
import pandas as pd
import scipy
from scipy.io import savemat
import statsmodels.api as sm
from fredapi import Fred
fred = Fred(api_key='d35aabd7dc07cd94481af3d1e2f0ecf3')


# Set display precision for pandas and numpy
pd.set_option('display.precision', 3)
np.set_printoptions(precision=3)

# Custom function imports
from time_series_functions import (
    moments,
    filter_transform,
    stacked_moments
)
from formatting_functions import (
    create_stats_table,
    generate_stacked_moments_table,
    generate_stacked_moments_latex_table,
    generate_latex_subtables
)

import matplotlib.dates as mdates

# FRED API key

# Set up date formatting for plots (if needed)
years = mdates.YearLocator(5, month=1)
years_fmt = mdates.DateFormatter('%Y')

# === 1. Load and Prepare Raw Data ===
lab = ['c', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'delta', 'w', 'bf', 'ba']
init = '1951-01-01'
final = '2025-12-30'

# Raw data series created by file observables.py
dat = pd.read_pickle("raw_data.pkl")

# === 2. Apply HP Filter for Business Cycle Detrending ===
filter_type="hp_filter"

if filter_type == "hp_filter":
    cycle_hp = pd.concat([
        filter_transform(dat[x], init=init, final=final, transform_type='log',
                         filter_type="hp_filter", lamb=100_000)
        for x in lab
    ], axis=1)
    cycle_hp.columns = lab
    cycle = cycle_hp
    #Optionally, consider Hamilton regression filer
elif filter_type == "hamilton":
    cycle_ham = pd.concat([
        filter_transform(dat[x], init=init, final=final, transform_type='log',
                         filter_type="hamilton")
        for x in lab
    ], axis=1)
    cycle_ham.columns = lab
    cycle = cycle_ham
    
  
# === 3. Compute Target Moments ===
mom_list = ["u", "v", "s", "jf", "delta", "bf", "lp"]


# Compute moments (relative std to 'lp', label unemployment and productivity)
mom = moments(cycle[mom_list], relative_std="lp", lab=["u", "lp"])
mom_stacked = stacked_moments(cycle, mom_list)
mom_stacked.columns = ["Values"]

# Summarize and print moments
mom_tex = create_stats_table(mom)
print(mom_tex)

# Export empirical moments to .mat file
savemat('moments_empirical.mat', mom_stacked.to_dict('list'))

# Generate and print LaTeX and text tables for stacked moments
tab_tex = generate_stacked_moments_latex_table(mom_stacked)
tab = generate_stacked_moments_table(mom_stacked)
print(tab)

# === 4. Compare Empirical Moments to Model Moments ===

compare_moments_model_data = False

if compare_moments_model_data:
    mom_model = scipy.io.loadmat("model_moments.mat")["model_moments"].flatten()
    mom_stacked["Model_moments"] = mom_model
    df = mom_stacked
    
    # Filter DataFrame for different types of moments
    std_devs = df[df.index.str.startswith('std')]
    correlations = df[df.index.str.startswith('Cor') & ~df.index.str.contains('theta')]
    autocorrelations = df[df.index.str.contains('_{-1}')]
    
    
    # === 5. Generate LaTeX Table Comparing Data and Model Moments ===
    
    latex_table_code = generate_latex_subtables(df)
    print(latex_table_code)
        
    # === 6. (Optional) Extended Moments with Delta ===
    # Uncomment if you want to include delta in the moments
    # mom_ext_delta = mom_list + ["delta"]
    # mom_stacked_delta = stacked_moments(cycle, mom_ext_delta)
    # mom_stacked_delta.columns = ["Values"]
    # tab_tex_delta = generate_stacked_moments_latex_table(mom_stacked_delta)
    # tab_delta = generate_stacked_moments_table(mom_stacked_delta)
    # print(tab_delta)
    
    # === 7. (Optional) Save additional observables for estimation ===
    # Uncomment and modify as needed
    # if save_observables:
    #     lab_obs = [x +'_obs' for x in lab]
    #     dic_data = dict(zip(lab_obs, [np.asarray(cycle_growth[x]) for x in cycle_growth.columns]))
    #     savemat('observables.mat', dic_data)


