import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio
import statsmodels.api as sm
#import statsmodels.api as sm
from fredapi import Fred
import pickle
fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')
#from statsmodels.tsa.arima_model import ARMA
pd.set_option('display.precision', 3)
np.set_printoptions(precision=3)
#pd.options.display.float_format = '{:5,.4g}'.format

from time_series_functions import (moments, stacked_moments, filter_transform, crosscorr, dynamic_correlations)
from statsmodels.tsa.seasonal import seasonal_decompose
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_lm

import matplotlib.dates as mdates
years = mdates.YearLocator(5, month=1)
years_fmt = mdates.DateFormatter('%Y')
import matplotlib.ticker as mtick
#from statsmodels.tsa.x13 import x13_arima_analysis
#https://www.census.gov/srd/www/x13as/ (binaries necessary to do seasonal decomposition)
#import statsmodels
#arima =  statsmodels.tsa.x13.x13_arima_analysis


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)  


def construct_data(init, final, freq):
    """
    Variables for estimation
    c: real per capita consumption 
    cons_share: consumption share of output
    u: unemployment rate (UNRATE)
    v: vacancy rate
    theta: market tightness 
    lp: labor productivity
    s: aggregate separation rate
    w: real hourly compensation for all workers
    SBF4: Business Formations Within FOUR Quarters
    SBF8: Business Formations Within Eight Quarters
    WBA: business applications with planned wages.
    High-Propensity Business Applications (HBA) that indicate a first wages-paid date on the IRS Form SS-4. 
    The indication of a wages-paid date is associated with a high likelihood of transitioning into a business with a payroll.
    """

    " GDP Deflator BEA code A191RD"
    #deflator =  fred.get_series('GDPDEF').resample(freq).mean()
    #Y = fred.get_series('GDPC1').resample(freq).mean().dropna() #real, quarterly
    
    C = fred.get_series('PCE').resample(freq).mean()
    #omega_SC = np.mean(C_S/(C))
    #print(rf'$\omega_{{SC}} =$ {omega_SC:.2f}')
    
    " Nominal investment: durables (PCDG), non-residential investment (PNFI), residential investment (PRFI) "
    #PCDG = fred.get_series("PCDG").resample(freq).mean().dropna()
    #PNFI = fred.get_series("PNFI").resample(freq).mean().dropna()
    #PRFI = fred.get_series("PRFI").resample(freq).mean().dropna()
    #I = PCDG + PNFI + PRFI
    I = fred.get_series('GPDI').resample(freq).mean()
    cons_share = C/(C+I)
    #np.mean(PCDG/I)
    " Price index of consumption goods "
    p_C = fred.get_series("PCEPI").resample(freq).mean().dropna()
    
    " Price index for investment goods "
    #p_I = fred.get_series("A006RD3Q086SBEA").resample(freq).mean().dropna()
    #p_I = fred.get_series("INVDEF").resample(freq).mean().dropna()
   
    #Equivalent to FPI (BEA A007RC) and change in business inventories (BEA CBI)
    #https://apps.bea.gov/iTable/?reqid=19&step=2&isuri=1&categories=underlying#eyJhcHBpZCI6MTksInN0ZXBzIjpbMSwyLDNdLCJkYXRhIjpbWyJjYXRlZ29yaWVzIiwiU3VydmV5Il0sWyJOSVBBX1RhYmxlX0xpc3QiLCIyMDcxIl1dfQ==
    #I = fred.get_series("GPDI").resample(freq).mean().dropna()
    
    " Non-institutional population "
    pop = fred.get_series('CNP16OV').resample(freq).mean()
    # Use HP-filtered trend for population to avoid discrete jumps around census dates
    pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
    
 
    
    c = C/(pop*p_C)
    #i = I/(pop*p_I)
    
    #Y = C + I
    #y = Y/(pop*deflator)
    #lab_prod = Y/(deflator*L)
    
    " Define labor productivity as real output per worker "
    #PRS85006163
    # seasonally adjusted real output per person in the non-farm business sector
    lp = fred.get_series("PRS85006163").resample(freq).mean().dropna()
    " Construct output from consumption and investment "
    #y = c + i
    " Relative price of investment: divide price indices"
    #p_I = p_I/p_C
    
    u = fred.get_series('UNRATE').resample(freq).mean().dropna()/100
    #fig, ax = plt.subplots()
    #ax.plot(u, label="UNRATE")
    #ax.plot(u_agg, label="Sectorally weighted unemployment rate")
    #ax.legend()
    " Construct idleness measures "
    
    
    " Vacancies "
    
    df = pd.read_csv("CompositeHWI.csv") # vacancy data based on Help-Wanted index
    date = pd.date_range(start="1951", end="2021-08", freq="MS")
    df.index = pd.DatetimeIndex(date)
    v = df["V/LF"]/100
    v = v.resample(freq).mean().dropna()
    # create dataframe
    theta = v/u
    
    """
     Separation rate and job finding rates: construct using Shimer 2005 approach
     u_{t+1}^0 = s_te_t(1-(1/2)f_t)
         u_t: number of workers unemployed with duration less than one month
         f_t: job finding rate
         e_t: number employed in month t
    Solve for s_t:
        s_t = u_{t+1}^0/(e_t(1-(1/2)f_t))
    
    We need job finding rate:
    u_{t+1} = u_t(1-f_t)+u_{t+1}^s =>
    f_t = 1 - (u_{t+1}-u_{t+1}^s)/u_t
    
    """
    
    u_level = fred.get_series('UNEMPLOY').loc[init:final] #number unemployed in thousands
    e_level = fred.get_series('CE16OV').loc[init:final] #number employed in thousands
    u_new = fred.get_series('UEMPLT5').loc[init:final] #unemployed for less than 5 weeks
    
    jf = 1-(u_level[1:len(u_level)-1]-u_new)/u_level# job finding rate series
    s = u_new/(e_level*(1-(1/2)*jf)) #separation rate series
    
    s = s.resample(freq).mean().dropna()
    f = jf.resample(freq).mean().dropna()
    
    " Labor share "
    ls = fred.get_series('PRS85006173').resample(freq).mean()
    
    " Wages "
    # Nonfarm Business Sector: Real Hourly Compensation for All Workers, index 2017=100
    #w = fred.get_series('COMPRNFB').resample(freq).mean()
    w = ls*lp
    
    " BFS "
    sbf4 = fred.get_series('BFBF4QTOTALSAUS').resample(freq).mean().dropna() / pop
    bawba = fred.get_series("BAWBATOTALSAUS").resample(freq).mean().dropna()/ pop
    #sbf8 = fred.get_series('BFBF8QTOTALSAUS').resample(freq).mean().dropna() / pop

    " Note: these series imply labor productivity in each sector "
    " List of data series "
    var_load_list = [c, cons_share, u, v, theta, f, lp, ls, s, w, sbf4, bawba] 
    return var_load_list
        
if __name__ == "__main__":       
        # Baseline
    init= '1951-01-01'
    #final = '2024-10-30'
    final='2020-01-01' # Just before pandemic shock
    # Comparison to earlier BRS
    #init = '1967-01-01'
    load = False
    #filter_type = 'hamilton'
    freq = 'QS'
    save_observables = False
    #load the BFS data, and convert the BFS data to quarterly frequency
    
    if load:
        var_load_list = pickle.load(open("var_load_list", "rb"))
    else:
        var_load_list = construct_data(init, final, freq)
        save_object(var_load_list, 'var_load_list')
    
    dat = pd.concat(var_load_list, axis=1)
    lab = ['c', 'cons_share', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'w', 'bf', 'ba']
    dat.columns = lab
    dat = dat.loc[init:final]
    
    " Data series in growth rates "
    #cycle_growth = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                      #  filter_type="growth", demean=False) for x in lab], axis=1)
    #cycle_growth.columns = lab
    #print(cycle_growth.mean())
    #cycle_growth = cycle_growth - cycle_growth.mean()
    cycle_hp = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                        filter_type="hp_filter", lamb=10_000) for x in lab], axis=1)
    cycle_hp.columns = lab
    cycle_hp.drop(['jf', 'cons_share'], axis=1, inplace=True)
    
    cycle_ham = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                        filter_type="hamilton") for x in lab], axis=1)
    cycle_ham.columns = lab
    cycle_ham.drop(['jf', 'cons_share'], axis=1, inplace=True)
    
    cycle_ham[["bf", "ba"]].corr()
    
    # Stacked moments 
    
    mom = moments(cycle_ham, relative_std="lp", lab=["u", "lp"])
    mom_stacked = stacked_moments(cycle_hp)
    " Summarize moments in one column "
 
    def create_stats_table(data, caption="Statistical Summary", label="tab:stats"):
        """
        Create a LaTeX table with dynamic column handling
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame containing the statistical measures
        caption : str
            Table caption
        label : str
            Table reference label
        """
        # Get column names dynamically from the DataFrame
        columns = data.columns
        
        latex_str = [
            "\\begin{table}[htbp]",
            "\\centering",
            f"\\caption{{{caption}}}",
            f"\\label{{{label}}}",
        #    "\\begin{threeparttable}",
            # Create dynamic column format based on number of columns
            f"\\begin{{tabular}}{{l{''.join(['r'] * len(columns))}}}",
            "\\toprule"
        ]
        
        # Create header row dynamically
        # Replace potentially problematic characters and add LaTeX formatting
        header_row = ["Variable"] + [
            col.replace("_", "\\_")  # Escape underscores
               .replace("-", "$-$")  # Format minus signs
               .replace("(", "\\left(").replace(")", "\\right)")  # Format parentheses
            for col in columns
        ]
        latex_str.append(" & ".join(header_row) + " \\\\")
        
        latex_str.append("\\midrule")
        
        # Add data rows with proper formatting
        for idx, row in data.iterrows():
            formatted_row = [
                f"{idx}"  # Variable name
            ] + [
                f"{val:.3f}" if isinstance(val, (int, float)) else str(val)
                for val in row
            ]
            latex_str.append(" & ".join(formatted_row) + " \\\\")
        
        latex_str.extend([
            "\\bottomrule",
            "\\end{tabular}",
          #  "\\end{threeparttable}",
            "\\end{table}"
            ])
    
        return "\n".join(latex_str)
    
    mom_tex = create_stats_table(mom)
    print(mom_tex)
    
    
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
    X = sm.add_constant(cycle_hp.theta)
    Y = cycle_hp.jf
    model = sm.OLS(Y,X).fit(cov_type='HAC', cov_kwds={'maxlags':None})
    print(model.summary())
    alpha_hat = 1-model.params.theta
    #A = np.exp(model.params.const)
    
    # Test for heteroskedasticity and autocorrelation (shifts of Beveridge curve)
    resid = model.resid
    white_test = sm.stats.diagnostic.het_white(resid, X)
    print(f"Test statistic: {white_test[0]}")
    print(f"P-value: {white_test[1]}")
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
    # Not in current data: but separation rates spike dramatically in March 2023
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
        
        stats_df['Value'] = np.around(stats_df['Value'], decimals=3)
        
        return stats_df
        
    data_means = first_moments(dat, alpha_hat)
    print(data_means.to_latex())
    # Comovement of unemployment, labor productivity, N
    # cycle_growth_labor = cycle_growth[['u', 'v','lab_prod', 'N']]
    # mom_growth_labor = moments(100*cycle_growth_labor, lab=['u','lab_prod'], lags=[1])
    
    # # Dependence on time period
    # corr_pre = cycle_growth_labor[init:"1983"].corr()
    # corr_post = cycle_growth_labor["1984":final].corr()
    # corrs_labor = pd.concat([corr_pre, corr_post])
    
    # Sectoral comovement
    
    
     
    
   
    
        
        
      
    
