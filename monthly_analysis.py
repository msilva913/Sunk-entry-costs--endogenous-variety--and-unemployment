
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import statsmodels.api as sm
from fredapi import Fred
import pickle
fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')
#from statsmodels.tsa.arima_model import ARMA
pd.set_option('display.precision', 2)
#pd.options.display.float_format = '{:5,.4g}'.format
import matplotlib.dates as mdates
from time_series_functions import (moments, filter_transform, crosscorr)
#from observables import construct_data_monthly
# In levels 
# 1) Plot
years = mdates.YearLocator(5, month=1)
years_fmt = mdates.DateFormatter('%Y')

init= '1964-01-01'
final = '2019-12-30'
linestyle = ['-', ':', '-.']

def construct_data_monthly(init, final, freq='MS'):
    
    # Unemployment rates by sector
    #https://www.bls.gov/web/empsit/cpseea31.htm
    
  
    " Capacity utilization "
    #util = fred.get_series('TCU').resample(freq).mean().dropna()
    # Durable manufacturing
    util_D = fred.get_series('CAPUTLGMFDS').resample(freq).mean().dropna()
    # nondurable manufacturing
    util_ND = fred.get_series('CAPUTLGMFNS').resample(freq).mean().dropna()
    
    u = fred.get_series('UNRATE').resample(freq).mean().dropna()

    # Number of unemployed
    #u = df.iloc[:,1]
    unemp = fred.get_series("UNEMPLOY").resample(freq).mean().dropna()
    # Labor force
    #lf = df.iloc[:, 2]
    lf = fred.get_series("CLF16OV").resample(freq).mean().dropna()
    # Employed 
    e = lf - unemp
    # Vacancies
    df = pd.read_csv("https://raw.githubusercontent.com/letsgoexploring/economic-data/master/"
                      "dmp/csv/beveridge_curve_data.csv")
    dates = df.iloc[:,0]
    index = pd.DatetimeIndex(dates)
    v = df.iloc[:,3]
    v.index = index
    # Market tightness
    theta = v/u
    # Newly unemployed 
    u_new = fred.get_series('UEMPLT5').resample(freq).mean().dropna() #unemployed for less than 5 weeks
    f = 1-(unemp.shift(-1) - u_new.shift(-1))/unemp
  
    var_load_list = [util_D, util_ND, u, f] 
    return var_load_list


var_load_list = construct_data_monthly(init, final)

util_D, util_ND, u, f = var_load_list
e = 100 - u
dat_levels = pd.concat([util_D, util_ND, e, f], axis=1)
lab = ['util_D', 'util_ND', 'e', 'f']
dat_levels.columns = lab

fig, ax = plt.subplots(figsize=(11, 4))
ax.plot(dat_levels.util_D, linestyle[0], label= "Utilization: durables", lw=2, alpha=0.7)
ax.plot(dat_levels.util_ND, linestyle[1], label= "Utilization: nondurables", lw=2, alpha=0.7)
ax.plot(dat_levels.e, linestyle[2], label="Employment rate", lw=2, alpha=0.7)
ax.legend(loc="upper right")
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.set_ylabel("%")
ax.grid(True)
plt.tight_layout()
plt.savefig("util_e_levels.pdf")
plt.show()


# 1) Contemporaneous correlations
dat_levels.corr()

# 2) Dynamic correlations
def dynamic_correlations(cycle, var, ylabel, nleads=12, nlags=12, title=None,
                         savefig=None):
    fig, ax = plt.subplots(figsize=(14, 5))
    # Loop over variables
    x = cycle[var[0]]
    y = cycle[var[1]]
    rs = pd.Series([x.corr(y.shift(-lag))] for lag in range(-nlags, nleads))
    rs = pd.Series([crosscorr(cycle[var[0]], cycle[var[1]], -lag) for lag in range(-nlags, nleads)])
    ax.axhline(y=0.0, color="black", linestyle="--")
    ax.plot(rs, '-', alpha=0.7, linewidth=2.0)
    ax.axvline(np.argmax(abs(rs)),linestyle='--', color='red')
    ax.set_xlabel(r'$\Delta$ (monthly)', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(range(0, len(rs)))
    ax.set_xticklabels(range(-nlags, nleads))
    if title is not None:
        ax.set_title(title, fontsize=14)
    #ax.legend(loc = "upper left", fontsize=14)
    ax.axvline(nlags, color="black")
    plt.tight_layout()
    if savefig is not None:
        plt.savefig(savefig)

ylabel=r'$Corr(util_{D,t}, e_{t+\Delta})$'
labels = ["Utilization", "Employment"]
var = ["util_D", "e"]

dynamic_correlations(dat_levels, var, ylabel, title="Dynamic correlations between utilization and employment")
var = ["util_ND", "e"]
dynamic_correlations(dat_levels, var, ylabel, title="Dynamic correlations between utilization and employment")

labels = ["Utilization", "Employment"]
var = ["util_D", "f"]    
dynamic_correlations(dat_levels, var, ylabel, title="Dynamic correlations between utilization and job finding rate")

    # Granger causality tests
from statsmodels.tsa.stattools import grangercausalitytests
dat_levels_red = dat_levels[["util_D", "e"]].dropna()

def granger_causality_test(series1, series2, maxlag=12):
    # Tests wehther time series in second column Granger causes time series in the first column
    # Null hypothesis: x2 does not Granger cause x1
    data = pd.DataFrame({'Series1': series1, 'Series2': series2})
    data = data.dropna()
    res = grangercausalitytests(data, maxlag=maxlag)
    return res[3][0]['ssr_ftest'][1]  # Get p-value for lag 4

p_value_utilD_e = granger_causality_test(util_D, e)
p_value_e_utilD = granger_causality_test(e, util_D)

print("Granger causality from e to util_D:")
print("p-value:", p_value_utilD_e)

print("Granger causality from util_D to e:")
print("p-value:", p_value_e_utilD)

" Vector autoregression "
from statsmodels.tsa.api import VAR
# Slow moving variables should be first
data = pd.concat([e, util_D], axis=1)
labels = ["e", "util_D"]
data.columns = labels
data = data.dropna()


model = VAR(data)
res = model.fit(maxlags=12, trend="n")
res.sigma_u
irf = res.irf(36)
irf.plot(orth=False)
fevd = res.fevd(24)

" FEVD from utilization "
e_fevd = fevd.decomp[0]
fig, ax = plt.subplots()
ax.plot(e_fevd[:,1], label="Variance decomposition from utilization")
ax.legend()
plt.show()


def var_estimate(data, p):
    # Design matrix (nxk), k=(nlags*nvars)
    X = pd.concat([data.shift(periods=i) for i in range(1,p+1)], axis=1)
    # Lose p observations from missing values
    X = X.dropna()
    Y = data.iloc[p:] 
    # OLS regression: Y = Xbeta + eps
    A = np.linalg.solve(X.T@X, X.T@Y)
    # Want to rewrite as  A = [A_1, .... A_p], where A_i is nxn (2x2 in simple case)
    return A.T, X, Y


# def structural_irf(A, B0_inv, T):
#     # y_t = A_1 y_{t-1} + A_2 y_{t-2} + ... + A_p y_{t-p} + u_t
#     n = A.shape[0] #number of variables n
#     p = A.shape[1] // n #number of lags
#     # sirf[i, j, k]: response of variable j to one-unit shock to variable k after i periods
#     sirf = np.zeros((T + 1, n, n))
#     # Shocks are initialized to one (orthogonal
#     #sirf[0] = np.eye(n)
#     #eps_t = B0_inv u_t
#     sirf[0] = B0_inv
#     A_i = A.reshape((n, n, p), order='F')
#     for i in range(1, T + 1): # periods for which irf will be computed
#         sirf[i] = A_i[:,:,0]@sirf[i - 1]
#         for j in range(1, p):
#            # if i >= j+1:
#                 sirf[i] += A_i[:, :, j]@sirf[i-1-j]
#     return sirf

def structural_irf(A, B0_inv, T):
    # y_t = A_1 y_{t-1} + A_2 y_{t-2} + ... + A_p y_{t-p} + u_t
    n = A.shape[0] #number of variables n
    p = A.shape[1] // n #number of lags
    # sirf[i, j, k]: response of variable j to one-unit shock to variable k after i periods
    sirf = np.zeros((T + 1, n, n))
    # Shocks are initialized to one (orthogonal
    #sirf[0] = np.eye(n)
    #eps_t = B0_inv u_t
    sirf[0] = B0_inv
    A_i = A.reshape((n, n, p), order='F')
    for i in range(1, T + 1): # periods for which irf will be computed
        sirf[i] = sirf[i - 1].dot(A_i[:,:,0])
        for j in range(1, p):
           # if i >= j+1:
                sirf[i] += sirf[i-1-j].dot(A_i[:,:,j])
    return sirf


# Estimation to get A_i, reduced-form errors eps_t, and covariance matrix E(eps*eps')=Sigma_eps
A, X, Y = var_estimate(data, 12) #n x (nxp)
pred = X@(A.T)
pred.columns = labels
e_pred = pred.iloc[:,0]
fig, ax = plt.subplots()
ax.plot(e, label='Employment')
ax.plot(e_pred, label='Predicted employment', alpha=0.7)
ax.legend()

residuals = Y - pred
sigma = residuals.T.dot(residuals) / (residuals.shape[0] - 1)
B0_inv = np.linalg.cholesky(sigma)
assert (B0_inv@(B0_inv.T) - sigma).max().max() < 1e-4
#B = B_0*A
#B = np.linalg.inv(B0)@A
#B.T - A.T.dot(B0_inv)
# check
A - res.coefs.reshape((2, 24))
coefs_reshape = res.coefs.reshape((24, 2))

T = 36
sirf = structural_irf(A, B0_inv, T)
periods = range(T)
# Utilization shock on employment
imp = [sirf[i][0,1] for i in periods] # response of employment to utilization shock
fig, ax = plt.subplots()
ax.plot(periods, imp)
plt.show()


imp = [sirf[i][1, 0] for i in periods] # response of utilization to employment shock
fig, ax = plt.subplots()
ax.plot(periods, imp)
plt.show()
ax.plot(periods, imp, label="Employment shock on utilization")
plt.show()

imp = [sirf[i][1, 1] for i in periods] # response of utilization to utilization
fig, ax = plt.subplots()
ax.plot(periods, imp)
plt.show()
ax.plot(periods, imp, label="Utilization shock on utilization")
plt.show()