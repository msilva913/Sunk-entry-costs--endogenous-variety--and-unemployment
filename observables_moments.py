
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
#import statsmodels.api as sm
from fredapi import Fred
import pickle
fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')
#from statsmodels.tsa.arima_model import ARMA
pd.set_option('display.precision', 3)
np.set_printoptions(precision=3)
from scipy.io import savemat
#pd.options.display.float_format = '{:5,.4g}'.format

from time_series_functions import (moments, stacked_moments, filter_transform)
from formatting_functions import (create_stats_table, generate_stacked_moments_table, 
                                  generate_stacked_moments_latex_table)

import matplotlib.dates as mdates
years = mdates.YearLocator(5, month=1)
years_fmt = mdates.DateFormatter('%Y')
#from statsmodels.tsa.x13 import x13_arima_analysis
#https://www.census.gov/srd/www/x13as/ (binaries necessary to do seasonal decomposition)
#import statsmodels
#arima =  statsmodels.tsa.x13.x13_arima_analysis

" Load raw data "
lab = ['c', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'delta', 'w', 'bf', 'ba']
init= '1951-01-01'
#final = '2024-10-30'
final='2020-01-01' # Just before pandemic shock
dat = pd.read_pickle("raw_data.pkl")

cycle_hp = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                    filter_type="hp_filter", lamb=100_000) for x in lab], axis=1)
cycle_hp.columns = lab

cycle_ham = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                    filter_type="hamilton") for x in lab], axis=1)
cycle_ham.columns = lab

#cycle_ham[["bf", "ba"]].corr()
#cycle = cycle_ham
cycle = cycle_hp

# Stacked moments 
mom_list = ["u", "v", "theta", "lp", "s", "bf"]
mom = moments(cycle, relative_std="lp", lab=["u", "lp"])
mom_stacked = stacked_moments(cycle, mom_list)
" Summarize moments in one column "

mom_stacked.columns = ["Values"]
mom_tex = create_stats_table(mom)
print(mom_tex)
# Export to mat file 
mom_stacked_dic = mom_stacked.to_dict('list')
#savemat('moments_empirical.mat', mom_stacked.to_dict('list'))

tab_tex = generate_stacked_moments_latex_table(mom_stacked)
tab = generate_stacked_moments_table(mom_stacked)
print(tab)

mom_ext_delta = mom_list + ["delta"]
mom_stacked_delta = stacked_moments(cycle, mom_ext_delta)
mom_stacked_delta.columns = ["Values"]
tab_tex_delta = generate_stacked_moments_latex_table(mom_stacked_delta)
tab_delta = generate_stacked_moments_table(mom_stacked_delta)
print(tab_delta)

#savemat('moments_bf_empirical.mat', mom_stacked_bf.to_dict('list'))
# if save_observables:
#     " Save relevant objects "
#     #save_object(cycle, 'cycle')
#     " Save output for estimation using growth filter"
#     lab_obs = [x +'_obs' for x in lab]
#     dic_data = dict(zip(lab_obs, [np.asarray(cycle_growth[x]) for x in cycle_growth.columns]))
#     sio.savemat('observables.mat', dic_data)



