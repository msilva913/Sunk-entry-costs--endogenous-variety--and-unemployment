import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio

from observables import construct_data
from time_series_functions import moments

init= '1951-01-01'
final='2019-12-31' # Just before pandemic shock

dat_simp = construct_data(init, final, 'ME')
dat_simp = pd.concat(dat_simp, axis=1)
lab = ["c", "cons_share", "u", 'v', 'theta', 'f', 'lab_prod', 'ls', 's', 'w', 'sbf4', 'bawba']
dat_simp.columns = lab
dat_simp = dat_simp.loc[init:final]
dat_simp = dat_simp[["u", "v", "s", "f", "lab_prod", "w", "cons_share"]]

dat_simp.to_pickle('labor_data_monthly.pkl')

dat_obs = dat_simp[["u", "v", "s"]]
rename_dict = {'u': 'u_obs', 'v':'v_obs', 's':'s_obs'}
dat_obs.rename(columns=rename_dict, inplace=True)
#dat_obs = np.log(dat_simp[["u", "v"]]).diff()

# Productivity to monthly
def monthly_to_quarterly_growth(x):
    x_growth = np.log(x).resample('Q').last().diff()
    x_growth = x_growth.resample('M').last()
    return x_growth

out = [monthly_to_quarterly_growth(x) for x in [dat_simp.lab_prod, dat_simp.w]]

x_obs, w_obs = out

dat_obs['x_obs'] = x_obs
dat_obs['w_obs'] = w_obs

# Align the consumption share with the quarterly end points: 3-31, 6-30, 9-30, 12-31
cons_share = dat_simp.cons_share.resample('Q-DEC').last().resample('M').last()
dat_obs['cons_share_obs'] = cons_share

# Calculate moments

mom = moments(100*dat_obs, lab=['u_obs','x_obs'], lags=[1])

# Demean the data
dat_obs_dem = dat_obs[:]
for x in ["u_obs", "v_obs", "s_obs", "x_obs", "w_obs"]:
    dat_obs_dem[x] = dat_obs[x] - dat_obs[x].mean()
    

" Save observables to MAT file "
lab_obs = dat_obs.columns
#dic_data = dict(zip(lab_obs, dat_obs))
dic_data = dict(zip(lab_obs, [np.asarray(dat_obs_dem[x]) for x in lab_obs]))
#sio.savemat('observables_u_growth.mat', dic_data)
sio.savemat('observables_u_level.mat', dic_data)

" In growth rates"
dat_obs['u_obs'] = np.log(dat_simp.u).diff()
dat_obs['v_obs'] = np.log(dat_simp.v).diff()

print(dat_obs.mean()*100)
mom_growth = moments(100*dat_obs, lab=['u_obs','x_obs'], lags=[1])

dat_obs_dem = dat_obs[:]
for x in ["u_obs", "v_obs", "s_obs", "x_obs", "w_obs"]:
    dat_obs_dem[x] = dat_obs[x] - dat_obs[x].mean()

dic_data = dict(zip(lab_obs, [np.asarray(dat_obs_dem[x]) for x in lab_obs]))
#sio.savemat('observables_u_growth.mat', dic_data)
sio.savemat('observables_u_growth.mat', dic_data)

