import pandas as pd
pd.set_option('display.max_columns', 10) 
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.filters.hp_filter import hpfilter
from statsmodels.tsa.arima.model import ARIMA

# Load Business Dynamics Statistics data
df = pd.read_pickle("BDS_data_adj.pkl")


# Convert AR processes to monthly from higher frequencies
def AR1_conversion_upcast(rho, sigma_sq, n=3):
    rho_m = rho**(1/n) # conversion of AR1 persistence parameter to higher freq
    sigma_sq_n = (1-rho_m**(2*n))/(1-rho_m**2)*sigma_sq
    return rho_m, sigma_sq_n
    

" Destruction shocks: establishments and firms "
# Create lists
rho_delta_mon_seq = []
sigma_delta_mon_seq = []
exit_rate_cycle_seq = []
exit_rate_trend_seq = []

for x in ["estabs_exit_rate", "firms_exit_rate"]:
    log_exit_rate = np.log(df[x])
    lam_ann = 100_000/256 # adj to annual frequency
    exit_rate_cycle, exit_rate_trend = hpfilter(log_exit_rate, lam_ann)
    exit_rate_cycle_seq.append(exit_rate_cycle)
    exit_rate_trend_seq.append(exit_rate_trend)
    
    # Fit AR(1)
    model = ARIMA(exit_rate_cycle, order=(1, 0, 0)) # AR(1) process
    fit = model.fit()
    print(fit.summary())
    rho_delta_ann = fit.params["ar.L1"]
    sigma_sq_ann = fit.params["sigma2"]


    rho_delta_mon, sigma_sq_mon = AR1_conversion_upcast(rho_delta_ann, sigma_sq_ann, n=12)
    sigma_mon = np.sqrt(sigma_sq_mon)
    # Append results
    rho_delta_mon_seq.append(rho_delta_mon)
    sigma_delta_mon_seq.append(sigma_mon)
    
    print("rho_delta_mon=",rho_delta_mon)
    print("sigma_mon=", sigma_mon)

 # fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
 # ax[0].plot(log_exit_rate, label='Logged exit rate')
 # ax[0].plot(exit_rate_trend, label='Trend', linestyle='--')
 # ax[1].plot(exit_rate_cycle, label='Cycle', linestyle='--')
 # plt.legend()
 # plt.tight_layout()
 # plt.show()
 
# Test

" Aggregate separation shocks "
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

" Productivity shocks "
lab_prod = labor_data_mon.lab_prod
lab_prod = lab_prod.resample("QE").mean()
lab_prod_log = np.log(lab_prod)
lp_cycle, lp_trend = hpfilter(lab_prod_log, 100_000)

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(lp_cycle, label="Productivity shocks")
ax.plot(s_cycle, label="Separation shocks")
ax.legend()
plt.show()

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
v2 = np.array((rho_delta_mon_seq[0], sigma_delta_mon_seq[0]))
v3 = np.array((rho_delta_mon_seq[1], sigma_delta_mon_seq[1]))
v4 = np.array((rho_tau, sigma_tau))
comb = np.vstack([v1, v2, v3, v4]).T
headers = ["Productivity", "Establishment exit", "Firm exit", "Aggregate separation"]
table = tabulate(comb, headers=headers, tablefmt = "outline")
print(table)



