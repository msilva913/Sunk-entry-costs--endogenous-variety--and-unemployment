import pandas as pd
pd.set_option('display.max_columns', 10) 
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from fredapi import Fred


# Data
#download at https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html
#df = pd.read_csv('BDS_Extension.csv')
df = pd.read_csv('bds2022.csv')
fred = Fred(api_key="8d302fb5f121b2be4f7d6e795194bcde")
df['year'] = pd.to_datetime(df['year'], format='%Y')
df.set_index("year", inplace=True)

y = fred.get_series('GDPC1').resample('AS').mean()
pop = fred.get_series('CNP16OV').resample('AS').mean()
# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
# Population adjusted
y = y/pop
# Growth rate
y = np.log(y).diff()
y = y.to_frame()
y.rename(columns={0: "real_gdp_growth_rate"}, inplace=True)

df = pd.concat([df, y], axis=1)
df = df.loc["1978":"2022"]

#df.rename(columns={'firms': 'Firms', 'B': 'Column2'}, inplace=True)
# Episodes

df['firms_entry'] = df.firms.diff() + df.firmdeath_firms
df['firms_entry_rate'] = 100*df['firms_entry']/df['firms']
df['firms_exit_rate'] = 100*df['firmdeath_firms']/df['firms']

df.to_pickle("BDS_data_adj.pkl")
df = pd.read_pickle("BDS_data_adj.pkl")

" Estimate AR(1) process "
from statsmodels.tsa.filters.hp_filter import hpfilter
from statsmodels.tsa.arima.model import ARIMA


log_exit_rate = np.log(df["estabs_exit_rate"])
lam_ann = 100_000/256
exit_rate_cycle, exit_rate_trend = hpfilter(log_exit_rate, lam_ann)

fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
ax[0].plot(log_exit_rate, label='Logged firm exit rate')
ax[0].plot(exit_rate_trend, label='Trend', linestyle='--')
ax[1].plot(exit_rate_cycle, label='Cycle', linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()

# Fit AR(1)
model = ARIMA(exit_rate_cycle, order=(1, 0, 0)) # AR(1) process
fit = model.fit()
print(fit.summary())
rho_delta_ann = fit.params["ar.L1"]
sigma_sq_ann = fit.params["sigma2"]


# Convert to monthly
def AR1_conversion_upcast(rho, sigma_sq, n=3):
    rho_m = rho**(1/n) # conversion of AR1 persistence parameter to higher freq
    sigma_sq_n = (1-rho_m**(2*n))/(1-rho_m**2)*sigma_sq
    return rho_m, sigma_sq_n
    

rho_delta_mon, sigma_sq_mon = AR1_conversion_upcast(rho_delta_ann, sigma_sq_ann, n=12)
sigma_mon = np.sqrt(sigma_sq_mon)
print("rho_delta_mon=",rho_delta_mon)
print("sigma_mon=", sigma_mon)
# Test

labor_data_mon = pd.read_pickle('labor_data_monthly.pkl')

s = labor_data_mon.s.dropna()
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(s, label="Separation rate", color="red", alpha=0.6)
plt.tight_layout()
plt.show()
print("Mean separation rate= ", np.mean(s))

log_s = np.log(s) 
s_cycle, s_trend = hpfilter(log_s, 100_000)
#s_cycle = s_cycle.loc["1951":"2003"]
s_cycle = s_cycle - np.mean(s_cycle)
fig, ax = plt.subplots(nrows=2, figsize=(10, 10))
ax[0].plot(log_s, label='Log separation rate')
ax[0].plot(s_trend, label='Trend', color="red", linestyle='--', alpha=0.7)
ax[1].plot(s_cycle, label='Cycle', linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()

# Fit AR(1) to monthly cyclical separation rate data
model = ARIMA(s_cycle, order=(1, 0, 0)) # AR(1) process
fit = model.fit()
print(fit.summary())
rho_tau = fit.params["ar.L1"]
sigma_tau_sq = fit.params["sigma2"]

sigma_tau = np.sqrt(sigma_tau_sq)
print("rho_tau=",rho_tau)
print("sigma_tau", sigma_tau)

# Productivity shocks
lab_prod = labor_data_mon.lab_prod
lab_prod = lab_prod.resample("QE").mean()
lab_prod_log = np.log(lab_prod)
lp_cycle, lp_trend = hpfilter(lab_prod_log, 100_000)

model = ARIMA(lp_cycle, order=(1, 0, 0)) # AR(1) process
fit = model.fit()
print(fit.summary())
rho_lp = fit.params["ar.L1"]
sigma_lp_sq = fit.params["sigma2"]

# Convert to monthly
rho_lp_mon, sigma_lp_sq_mon = AR1_conversion_upcast(rho_lp, sigma_lp_sq, n=3)
sigma_lp_mon = np.sqrt(sigma_lp_sq_mon)

# Make simple table
from tabulate import tabulate
v1 = np.array((rho_lp_mon, sigma_lp_mon))
v2 = np.array((rho_delta_mon, sigma_mon))
v3 = np.array((rho_tau, sigma_tau))
comb = np.vstack([v1, v2, v3]).T
headers = ["Productivity", "Product destruction", "Aggregate separation"]
table = tabulate(comb, headers=headers, tablefmt = "outline")
print(table)



