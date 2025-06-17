
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
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
from formatting_functions import create_stats_table
from statsmodels.stats.diagnostic import acorr_lm

import matplotlib.dates as mdates
years = mdates.YearLocator(5, month=1)
years_fmt = mdates.DateFormatter('%Y')
import matplotlib.ticker as mtick
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


#savemat('moments_bf_empirical.mat', mom_stacked_bf.to_dict('list'))
# if save_observables:
#     " Save relevant objects "
#     #save_object(cycle, 'cycle')
#     " Save output for estimation using growth filter"
#     lab_obs = [x +'_obs' for x in lab]
#     dic_data = dict(zip(lab_obs, [np.asarray(cycle_growth[x]) for x in cycle_growth.columns]))
#     sio.savemat('observables.mat', dic_data)



" Analysis on raw data "
"1) Unemployment, Vacancies, and ustar "
ustar = np.sqrt(dat.u*dat.v) # efficient unemployment rate (FERU) using approxiximation by Michaillat and Saez (2024)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(dat.u["1951":], alpha=0.7, label="Unemployment rate")
ax.plot(dat.v["1951":], alpha=0.7, label="Vacancy rate")
#ax.plot(ustar["1951":], alpha=0.7, label="FERU")
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.legend()
plt.savefig("u_vac_series")
plt.show()


"2) Beveridge curve (levels and filtered) "
fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
ax[0].scatter(dat.u, dat.v, alpha=0.7)
ax[0].set_xlabel("Unemployment rate")
ax[0].set_ylabel("Vacancy rate")
ax[0].set_title("Levels")
ax[1].scatter(cycle_hp.u, cycle_hp.v, alpha=0.7)
ax[1].set_xlabel("Unemployment rate")
ax[1].set_ylabel("Vacancy rate")
ax[1].set_title("HP-filtered")
plt.show()
dat[['u','v']].corr()


    
"3) Estimate elasticity of Beveridge curve "
#log_u = np.log(dat.u)
#log_v = np.log(dat.v)
u = cycle_hp.u
v = cycle_hp.v
#X = log_u
X = sm.add_constant(u)
Y = v
#Newey west standard errors to control for heteroskedasticity and autocorrelation
#model = sm.OLS(Y,X).fit(cov_type='HAC', cov_kwds={'maxlags':None})
model = sm.OLS(Y,X).fit()
print(model.summary())
# Plot
plt.figure(figsize=(10, 6))
plt.scatter(u, Y, alpha=0.5)
plt.plot(u, model.predict(X), color='red', linewidth=2, alpha=0.7, label=f"Slope={model.params.u:,.2f}")
plt.xlabel('Log unemployment rate')
plt.ylabel('Log vacancy rate')
plt.title('Beveridge curve with regression line')
plt.legend()
plt.savefig("Beveridge_logs.pdf")
plt.show()

"5) Estimate matching function "
# Impose m = Au^alpha v^(1-alpha)
# Implies f = Atheta^(1-alpha)
# In logs: log f = log A + (1-alpha)log theta
#log_theta = np.log(dat.v/dat.u)
#log_f = np.log(dat.f)
#X = sm.add_constant(log_theta)
#Y = log_f
#X = sm.add_constant(cycle_hp.theta)
reg = pd.concat([cycle_hp.theta, cycle_hp.jf], axis=1).dropna()
model = sm.OLS(reg.jf,reg.theta).fit(cov_type='HAC', cov_kwds={'maxlags':None})
print(model.summary())
alpha_hat = 1-model.params.theta
#A = np.exp(model.params.const)

# Test for heteroskedasticity and autocorrelation (shifts of Beveridge curve)
resid = model.resid
#white_test = sm.stats.diagnostic.het_white(resid, X)
#print(f"Test statistic: {white_test[0]}")
#print(f"P-value: {white_test[1]}")
# Autocorrelation

print("\nAutocorrelation Test (Breusch-Godfrey):")
bg_test = acorr_lm(resid, nlags=1, store=True)

# Estimate the covariance matrix for heteroskedasticity and autocorrelation
# from statsmodels.regression.linear_model import GLS
# cov_matrix = statsmodels.stats.sandwich_covariance.cov_hac(model, nlags=1)

# # Fit the GLS model
# gls_model = GLS("Y ~ X", sigma=cov_matrix)
# gls_results = gls_model.fit()



"6) Unemployment and labor productivity "

# fig = plt.figure(figsize=(12, 10))
# ax1 = fig.add_subplot(2, 1,1)
# ax1.plot(dat.u, label="Unemployment")
# ax1.plot(lab_prod_growth, label="Labor productivity growth")
# ax1.xaxis.set_major_locator(years)
# ax1.xaxis.set_major_formatter(years_fmt)
# ax1.legend()

fig, ax = plt.subplots(figsize=(8, 4))
#ax1 = fig.add_subplot(2,1,1)
ax.scatter(cycle_hp.u, cycle_hp.lp, alpha=0.5)
ax.set_xlabel("Unemployment rate")
ax.set_ylabel("Labor productivity")
ax.legend()


"""
Stylized facts
1) Sectoral comovement (idleness, unemployment)
2) Within-sector comovement of unemployment and idleness
"""

"7) Separation and job finding rates "
# Not in current data: but separation rates spike dramatically in March 2020
fig, ax = plt.subplots(figsize=(11, 5))
ax.plot(cycle_hp.s, label="Separation rate", alpha=0.7, color="red", linestyle="--")
ax.plot(cycle_hp.jf, label="Job finding rate", alpha=0.7, color="blue", linestyle=":")
#ax.plot(cycle_hp.lp, label="Labor productivity", alpha=0.7, color="black")
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.legend()
plt.show()

"8) Elasticity of wages to productivity"
X = sm.add_constant(cycle_hp.lp)
Y = cycle_hp.w
model = sm.OLS(Y,X).fit(cov_type='HAC', cov_kwds={'maxlags':None})
print(model.summary())
summ = model.summary()
elast_w_lp = model.params.lp



def first_moments(dat, alpha_hat):

    # Create dictionary of statistics
    stats_dict = {
        'Consumption share': dat.cons_share.mean(),
        'Unemployment rate': dat.u.mean(),
        'Vacancy rate': dat.v.mean(),
        'Market tightness': dat.theta.mean(),
        'Separation rate': dat.s.mean(),
        'Job finding rate': dat.jf.mean(),
        'Matching function elasticity': alpha_hat,
       # 'SBF4': dat.sbf4.mean(),
       # 'SBF8': dat.sbf8.mean(),
    }
    
    # Convert to DataFrame with descriptive index
    stats_df = pd.DataFrame({
        'Value': stats_dict
    })
    
    #stats_df['Value'] = np.around(stats_df['Value'], decimals=3)
    stats_df['Value'] = stats_df['Value'].apply(lambda x: f"{x:.3f}")
    return stats_df
    
data_means = first_moments(dat, alpha_hat)
print(data_means.to_latex())

choice_vars1 = ["u", "v", "lp"]
choice_vars2 = ["u", "s", "bf"] 
savefigs=["pairplot_labor.pdf","pairplot_new.pdf"]
###################

# Define a mapping of old labels to new labels
label_mapping = {
    'u': 'Unemployment',
    'v': 'Vacancies',
    'lp': 'Labor productivity'
}

for item, choice_vars in enumerate([choice_vars1, choice_vars2]):
    pairplot = sns.pairplot(cycle_hp[choice_vars], diag_kind="kde",
                            plot_kws={'alpha': 0.7},
                            diag_kws={'color': 'gold'})
    pairplot.fig.set_size_inches(12, 8)  # Width, Height
    
    # Add best-fit lines to scatter plots
    for ax in pairplot.diag_axes:
        ax.set_xlim(cycle_hp[choice_vars[:-1]].min().min(), cycle_hp[choice_vars[:-1]].max().max())
    for ax in pairplot.axes.flatten():
        if ax.get_xlabel() in choice_vars and ax.get_ylabel() in choice_vars:
            sns.regplot(x=ax.get_xlabel(), y=ax.get_ylabel(), data=cycle_hp, ax=ax,
                        scatter=False, color='gray')
    
    # Set custom axis labels using label_mapping
    for ax in pairplot.axes.flatten():
        if ax.get_xlabel() in choice_vars:
            ax.set_xlabel(label_mapping.get(ax.get_xlabel(), ax.get_xlabel()))
        if ax.get_ylabel() in choice_vars:
            ax.set_ylabel(label_mapping.get(ax.get_ylabel(), ax.get_ylabel()))
    
    for i, j in zip(*np.triu_indices_from(pairplot.axes, 1)):
        pairplot.axes[i, j].set_visible(False)
        x = choice_vars[j]
        y = choice_vars[i]
        corr = cycle_hp[x].corr(cycle_hp[y])
        pairplot.axes[j, i].annotate(f'Corr: {corr:.2f}', xy=(0.1, 0.9), 
                                     xycoords='axes fraction', fontsize=10, 
                                     bbox=dict(boxstyle="round, pad=0.3", edgecolor="gray", facecolor="white"))
    
    # Save the figure before showing it
    pairplot.savefig(savefigs[item])
    plt.show()
    
     
    
   
    
        
        
      
    
