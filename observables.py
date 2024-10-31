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

from time_series_functions import (moments, filter_transform, crosscorr, dynamic_correlations)
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
    Nominal GDP: BEA 191RC
		- Nominal C: BEA DNDGRC + DDURRC + DSERRC
		- Nominal I: BEA A006RC
		- Nominal H: BLS PRS85006023
		- Capacity utilization: Fed Board G.17, CAPUTL.B00004.S
		- Relative price of investment: BEA B006RG
		- Nonfarm Nominal Hourly Wage: BLS PRS85006103
		- GDP Deflator: A191RD
		- 
    
    Sectoral data Tables B6 and B7
    B6: number of nonsupervisory employees
    B7: average weekly hours
    https://www.bls.gov/ces/data/employment-situation-table-download.htm

    Table B6
    https://data.bls.gov/pdq/SurveyOutputServlet
    Employment
    CES2000000006
    CES3100000006
    CES3200000006
    CES0800000006
    
    Average weekly hours
    CES2000000007
    CES3100000007
    CES3200000007
    CES0800000007
    
     Unemployment rate by industry 
    https://www.bls.gov/webapps/legacy/cpsatab14.htm
    
    Only available by industry back to 2000 because of major change in industry classification system
    https://www.bls.gov/cps/cpsoccind.htm
    Documentation on changes: https://www.bls.gov/cps/rvcps03.pdf
    
    Construction of Beveridge curve data by Brian Jenkins following Petrosky-Nadeau and Zhang (2013)
    https://github.com/letsgoexploring/economic-data/blob/master/dmp/python/us_beveridge_curve_data.py
    Results saved in file "https://raw.githubusercontent.com/letsgoexploring/economic-data/master/"
                      "dmp/csv/beveridge_curve_data.csv"
    """
    
    sectoral = pd.read_csv("sectoral_labor.csv", sep= ",", header=0)
    date = pd.date_range(start='1/1947', periods=sectoral.shape[0], freq='MS')
    sectoral.index = date
    
    # Construct total hours
    sectoral["Durable_TH"] = sectoral["Durable Emp"]*sectoral["Durable Hours"]
    sectoral["Construction_TH"] = sectoral["Construction Emp"]*sectoral["Construction Hours"]
    sectoral["Nondurable_TH"] = sectoral["Nondurable goods Emp"]*sectoral["Nondurable Hours"]
    sectoral["Services_TH"] = sectoral["Services Emp"]*sectoral["Services Hours"]
    
    sectoral = sectoral[['Durable_TH', 'Construction_TH', 'Nondurable_TH', 'Services_TH']]
    
    " Consumption: nondurable plus services "
    sectoral["L_C"] = sectoral["Nondurable_TH"] + sectoral["Services_TH"]
    
    " Investment: construction plus durables "
    sectoral["L_I"] = sectoral["Construction_TH"] + sectoral["Durable_TH"]
    
    # Total hours
    " Aggregate to quarterly "
    sectoral = sectoral.resample(freq).mean().dropna()
    " Data limited to only after 1964 (services only available since then) "
    LC = sectoral.L_C
    LI = sectoral.L_I
    L = LC + LI
    
    " GDP Deflator BEA code A191RD"
    deflator =  fred.get_series('GDPDEF').resample(freq).mean()
    #Y = fred.get_series('GDPC1').resample(freq).mean().dropna() #real, quarterly
    
    " Investment goods deflator "
    #inv_deflator = fred.get_series("INVDEF").resample(freq).mean()
    #cons_deflator = fred.get_series("CONSDEF").resample(freq).mean()           
    
    
    " Nominal consumption (BEA codes DNDGRC + DDURRC + DSERRC) "
    # Personal consumption non-durables: BEA DNDGRC
    #C_ND = fred.get_series('PCND').resample(freq).mean().dropna()# monthly, nominal
    # Personal consumption expenditure services: BEA DSERRC
    #C_S = fred.get_series('PCESV').resample(freq).mean().dropna()
    #C = C_ND + C_S
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
    p_I = fred.get_series("INVDEF").resample(freq).mean().dropna()
   
    #Equivalent to FPI (BEA A007RC) and change in business inventories (BEA CBI)
    #https://apps.bea.gov/iTable/?reqid=19&step=2&isuri=1&categories=underlying#eyJhcHBpZCI6MTksInN0ZXBzIjpbMSwyLDNdLCJkYXRhIjpbWyJjYXRlZ29yaWVzIiwiU3VydmV5Il0sWyJOSVBBX1RhYmxlX0xpc3QiLCIyMDcxIl1dfQ==
    #I = fred.get_series("GPDI").resample(freq).mean().dropna()
    
    " Non-institutional population "
    pop = fred.get_series('CNP16OV').resample(freq).mean()
    # Use HP-filtered trend for population to avoid discrete jumps around census dates
    pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
    
    " Unemployment rate by sector "
    u_sectoral = pd.read_csv("unemployment_industry.csv", sep= ",", header=0)
    date = pd.date_range(start='1/2000', periods=u_sectoral.shape[0], freq='MS')
    u_sectoral.index = date
    u_sectoral.drop('Date', axis=1, inplace=True)
    
    # Remove white space
    u_sectoral.columns = [s.strip() for s in u_sectoral.columns]
    
    # Unemployment rates for durables and non-durables
    u_D = u_sectoral["Durable-U rate"]/100
    u_ND = u_sectoral["Nondurable-U rate"]/100
    
    # Remove seasonal effects
    u_D_seas = seasonal_decompose(u_D, model='additive', period=12).seasonal
    u_ND_seas = seasonal_decompose(u_ND, model='additive', period=12).seasonal
    u_D = u_D - u_D_seas
    u_ND = u_ND - u_ND_seas
    
    u_D = u_D.resample(freq).mean().dropna()
    u_ND = u_ND.resample(freq).mean().dropna()
    
    " Unemployment levels in consumption and investment "
    # Investment: construction + durables
    u_level_I = u_sectoral["Construction-U"] + u_sectoral["Durable-U"]
    # Consumption: non-durables + services
    service_labels = ["Retail-U", "Transportation-U", "Information-U", "Financial-U", 
                      "Professional-U", "Eductation-Health-U", "Leisure-U", "Other-U"]
    u_level_services = sum(u_sectoral[x] for x in service_labels)
    u_level_ND = u_sectoral["Nondurable-U"]
    u_level_C = u_level_services + u_level_ND 
    
    " Employment levels in consumption and investment, based on B1 table "
    e_sectoral = pd.read_csv("employment_industry.csv", sep= ",", header=0)
    date = pd.date_range(start='1/1964', periods=e_sectoral.shape[0], freq='MS')
    e_sectoral.index = date
    e_sectoral.drop('Date', axis=1, inplace=True)
    e_sectoral = e_sectoral.dropna()
    
    e_C = e_sectoral["Nondurable goods-E"] + e_sectoral["Services-E"]
    e_I = e_sectoral["Construction-E"] + e_sectoral["Durable goods-E"]
    
    
    # Construct unemployment rates in consumption and investment using levels
    u_C = u_level_C/(u_level_C + e_C)
    u_I = u_level_I/(u_level_I + e_I)
    u_C.dropna(inplace=True)
    u_I.dropna(inplace=True)
    
    u_C_seas = seasonal_decompose(u_C, model='additive', period=12).seasonal
    u_C_adj = u_C - u_C_seas
    u_I_seas = seasonal_decompose(u_I, model='additive', period=12).seasonal
    u_I_adj = u_I - u_I_seas
   
    
    "Implied aggregate unemployment rate "
    u = u_level_C + u_level_I
    e = e_C + e_I
    LF = u + e
    C_weight = (u_level_C + e_C)/LF
    I_weight = (u_level_I + e_I)/LF
    u_agg = C_weight*u_C_adj + I_weight*u_I_adj
    
    u_C = u_C_adj.resample(freq).mean()
    u_I = u_I_adj.resample(freq).mean()
    
    c = C/(pop*p_C)
    i = I/(pop*p_I)
    
    Y = C + I
    y = Y/(pop*deflator)
    #lab_prod = Y/(deflator*L)
    
    " Define labor productivity as real output per worker "
    #PRS85006163
    # seasonally adjusted real output per person in the non-farm business sector
    lab_prod = fred.get_series("PRS85006163").resample(freq).mean().dropna()
    " Construct output from consumption and investment "
    #y = c + i
    lc = LC/pop
    li = LI/pop
    l = L/pop
    " Relative price of investment: divide price indices"
    p_I = p_I/p_C
    
    " Capacity utilization "
    util = fred.get_series('TCU').resample(freq).mean().dropna()/100
    # Durable manufacturing
    util_D = fred.get_series('CAPUTLGMFDS').resample(freq).mean().dropna()/100
    # nondurable manufacturing
    util_ND = fred.get_series('CAPUTLGMFNS').resample(freq).mean().dropna()/100
    
    u = fred.get_series('UNRATE').resample(freq).mean().dropna()/100
    #fig, ax = plt.subplots()
    #ax.plot(u, label="UNRATE")
    #ax.plot(u_agg, label="Sectorally weighted unemployment rate")
    #ax.legend()
    " Construct idleness measures "
    
    idle = 1.0 - util
    idle_D = 1.0 - util_D
    idle_ND = 1.0 - util_ND
    
    " Vacancies "
    # Vacancies
    # df = pd.read_csv("https://raw.githubusercontent.com/letsgoexploring/economic-data/master/"
    #                   "dmp/csv/beveridge_curve_data.csv")
    # dates = df.iloc[:,0]
    # df.index = pd.DatetimeIndex(dates)
    # #u_bc = df.iloc[:,1]
    # u_level = df.iloc[:, 1]
    # lf = df.iloc[:,2]
    # v = df.iloc[:,3]
    # #u_bc = u_bc/lf
    # # Create vacancy rate series by dividing vacancies by labor force
    # v = v/lf
    # v = v.resample(freq).mean().dropna()
    
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
    
   #  u_short = fred.get_series('UEMPLT5').dropna()
   # #e = fred.get_series('PAYEMS').dropna()
   #  e_level = lf - u_level
   #  f = 1 - (u_level - u_short.shift(-1))/u_level
   #  print(f["1951":"2003"].mean())
   #  s = u_short.shift(-1)/(e_level*(1-(1/2)*f))
    
   #  s = s.resample(freq).mean().dropna()
   #  f = f.resample(freq).mean().dropna()
    
    " Wages "
    # Nonfarm Business Sector: Real Hourly Compensation for All Workers, index 2017=100
    w = fred.get_series('COMPRNFB').resample(freq).mean()
    

    " Note: these series imply labor productivity in each sector "
    " List of data series "
    var_load_list = [y, c, i, cons_share, lc, li, l, lab_prod, p_I, idle, idle_D, idle_ND,
                     u_D, u_ND, u_C, u_I, u, v, theta, s, f, w] 
    return var_load_list
        
if __name__ == "__main__":       
        # Baseline
    init= '1951-01-01'
    #final = '2024-05-30'
    final='2020-02-01' # Just before pandemic shock
    # Comparison to earlier BRS
    #init = '1967-01-01'
    load = True
    #filter_type = 'hamilton'
    freq = 'Q'
    save_observables = False
    
    if load:
        var_load_list = pickle.load(open("var_load_list", "rb"))
    else:
        var_load_list = construct_data(init, final, freq)
        save_object(var_load_list, 'var_load_list')
    
    dat = pd.concat(var_load_list, axis=1)
    lab = ['Y', 'C', 'I', 'cons_share', 'NC', 'NI', 'N', "lab_prod", 'p_I', 'idle', 'idle_D', 'idle_ND', 'u_D', 'u_ND', 'u_C', 'u_I', 'u', 'v', 'theta', 's','f', 'w']
    dat.columns = lab
    dat = dat.loc[init:final]
    
    " Data series in growth rates "
    cycle_growth = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                        filter_type="growth", demean=False) for x in lab], axis=1)
    cycle_growth.columns = lab
    print(cycle_growth.mean())
    cycle_growth = cycle_growth - cycle_growth.mean()
    
    if save_observables:
        " Save relevant objects "
        #save_object(cycle, 'cycle')
        " Save output for estimation using growth filter"
        lab_obs = [x +'_obs' for x in lab]
        dic_data = dict(zip(lab_obs, [np.asarray(cycle_growth[x]) for x in cycle_growth.columns]))
        sio.savemat('observables_un.mat', dic_data)
    
    " Analysis on raw data "
    
    lab_idle = ['idle_D', 'idle_ND', 'u_D', 'u_ND', 'u_C', 'u_I', 'u', 'v']
    dat_idle = dat[lab_idle]
    
    
    
    "1) Aggregate unemployment rate and idleness rate "
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_title("Rates of slack")
    ax.plot(dat.u, label='Unemployment rate')
    ax.plot(dat.idle, label='Idleness rate')
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.legend()
    plt.savefig("unemployment_idleness_aggregate.pdf")
    plt.show()
    dat[["u", "idle"]].corr()
    dat[["u", "idle"]].mean()
    
    "2) Unemployment rates and idleness (comovement)"
    fig, ax = plt.subplots(nrows=2, figsize=(12, 8))
    ax[0].set_title("Unemployment rates")
    ax[0].plot(dat_idle.u_C, label='Consumption')
    ax[0].plot(dat_idle.u_I, label='Investment')
    ax[0].plot(dat_idle.u, label="Aggregate")
    ax[1].set_title("Idleness rates")
    ax[1].plot(dat_idle.idle_ND, label='Nondurables')
    ax[1].plot(dat_idle.idle_D, label='Durables')
    ax[1].plot(dat.idle, label='Aggregate')
    for i in range(2):
        ax[i].xaxis.set_major_locator(years)
        ax[i].xaxis.set_major_formatter(years_fmt)
        ax[i].legend()
    plt.savefig("unemployment_idleness_comovement.pdf")
    plt.show()
    
    
    "3) Idleness and unemployment rate (within-sector)"
    fig, ax = plt.subplots(nrows=2, figsize=(12, 10))
    ax[0].plot(dat_idle.idle_D, lw=2, label="idle_D")
    ax[0].plot(dat_idle.u_D, lw=2, label="u_D")
    ax[1].plot(dat_idle.idle_ND, lw=2, label="idle_ND")
    ax[1].plot(dat_idle.u_ND, lw=2, label="u_ND")
    for i in range(2):
        ax[i].xaxis.set_major_locator(years)
        ax[i].xaxis.set_major_formatter(years_fmt)
        ax[i].legend()
        ax[1].legend()
    plt.savefig("idleness_u.pdf")
    plt.show()
    
    "4) Beveridge curve "
    fig, ax = plt.subplots()
    ax.scatter(dat.u, dat.v, alpha=0.5)
    ax.set_xlabel("Unemployment rate")
    ax.set_ylabel("Vacancy rate")
    plt.show()
    dat[['u','v']].corr()
    
    "5) Unemployment, Vacancies, and ustar "
    ustar = np.sqrt(dat.u*dat.v) # efficient unemployment rate (FERU) using approxiximation by Michaillat and Saez (2024)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(dat.u["1951":], alpha=0.7, label="Unemployment rate")
    ax.plot(dat.v["1951":], alpha=0.7, label="Vacancy rate")
    ax.plot(ustar["1951":], alpha=0.7, label="FERU")
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    ax.legend()
    plt.savefig("u_vac_series")
    plt.show()
    
    "6) Estimate elasticity of Beveridge curve "
    log_u = np.log(dat.u)
    log_v = np.log(dat.v)
    #X = log_u
    X = sm.add_constant(log_u)
    Y = log_v 
    #Newey west standard errors to control for heteroskedasticity and autocorrelation
    model = sm.OLS(Y,X).fit(cov_type='HAC', cov_kwds={'maxlags':None})
    print(model.summary())
    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(log_u, Y, alpha=0.5)
    plt.plot(log_u, model.predict(X), color='red', linewidth=2, alpha=0.7, label=f"Slope={model.params.u:,.2f}")
    plt.xlabel('Log unemployment rate')
    plt.ylabel('Log vacancy rate')
    plt.title('Beveridge curve with regression line')
    plt.legend()
    plt.savefig("Beveridge_logs.pdf")
    plt.show()
    
    "7) Estimate matching function "
    # Impose m = Au^alpha v^(1-alpha)
    # Implies f = Atheta^(1-alpha)
    # In logs: log f = log A + (1-alpha)log theta
    log_theta = np.log(dat.v/dat.u)
    log_f = np.log(dat.f)
    X = sm.add_constant(log_theta)
    Y = log_f
    model = sm.OLS(Y,X).fit(cov_type='HAC', cov_kwds={'maxlags':None})
    print(model.summary())
    alpha_hat = 1-model.params[0]
    A = np.exp(model.params.const)
    
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
    
    
    
    "8) Unemployment and labor productivity growth "
    lab_prod_growth = np.log(dat.lab_prod).diff()
    # Average annual growth
    print(lab_prod_growth.mean()*4)
    # Correlation
    print(lab_prod_growth.corr(dat.u))
    
    # fig = plt.figure(figsize=(12, 10))
    # ax1 = fig.add_subplot(2, 1,1)
    # ax1.plot(dat.u, label="Unemployment")
    # ax1.plot(lab_prod_growth, label="Labor productivity growth")
    # ax1.xaxis.set_major_locator(years)
    # ax1.xaxis.set_major_formatter(years_fmt)
    # ax1.legend()
    
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(2,1,1)
    ax1.scatter(dat.u, lab_prod_growth, alpha=0.5)
    ax1.set_xlabel("Unemployment rate")
    ax1.set_ylabel("Labor productivity growth")
    ax1.legend()
    
    ax2 = fig.add_subplot(2,1,2)
    ax2.scatter(np.log(dat.u).diff(), lab_prod_growth, alpha=0.5)
    ax2.set_xlabel("Unemployment rate growth")
    ax2.set_ylabel("Labor productivity growth")
    ax2.legend()
    plt.savefig("unemployment_labor_productivity.pdf")
    plt.show()
    
    
    
    "9) Dynamic correlation between unemployment and labor productivity growth"
    ylabel=r'$Corr(u_{,t}, lab_prod_{t+\Delta})$'
    labels = ["Unemployment rate", "Labor productivity growth"]
    
    
    dat_red = pd.concat([dat.u, np.log(dat.u).diff(), dat.v, lab_prod_growth], axis=1)
    index = ["u", "u_growth", "v", "lab_prod_growth"]
    
    dat_red.columns = index
    var = ["u", "lab_prod_growth"]
    dynamic_correlations(dat_red, var, ylabel, title="Dynamic correlations")
    
    var = ["u_growth", "lab_prod_growth"]
    dynamic_correlations(dat_red, var, ylabel, title="Dynamic correlations")
    
    rolling_corr = dat_red["u"].rolling(window=30).corr(dat_red.lab_prod_growth)
    plt.plot(rolling_corr)
    rolling_corr = dat_red["u_growth"].rolling(window=30).corr(dat_red.lab_prod_growth)
    
    var = ["v", "lab_prod_growth"]
    dynamic_correlations(dat_red, var, ylabel, title="Dynamic correlations")
    
    corr_mat = dat_idle[["idle_D", "idle_ND", "u_D", "u_ND"]].corr()
    print(corr_mat)
    """
    Stylized facts
    1) Sectoral comovement (idleness, unemployment)
    2) Within-sector comovement of unemployment and idleness
    """
    
    "7) Separation and job finding rates "
    # Not in current data: but separation rates spike dramatically in March 2023
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(dat.s, label="Separation rate")
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(dat.f, label="Job finding rate")
    for ax in [ax1, ax2]:
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(years_fmt)
        ax.legend()
    plt.show()
    
    " Summary statistics"
    print("Mean consumption share:", dat.cons_share.mean())
    print("Mean Unemployment rate:", dat.u.mean())
    print("Mean vacancy rate:", dat.v.mean())
    print("Mean annual labor productivity growth:", 4*lab_prod_growth.mean())
    print("Mean separation rate:", dat.s.mean())
    print("Mean job finding rate:", dat.f.mean())
    print("Estimated matching function elasticity:", alpha_hat)
    
    
    
    # Comovement of unemployment, labor productivity, N
    cycle_growth_labor = cycle_growth[['u', 'v','lab_prod', 'N']]
    mom_growth_labor = moments(100*cycle_growth_labor, lab=['u','lab_prod'], lags=[1])
    
    # Dependence on time period
    corr_pre = cycle_growth_labor[init:"1983"].corr()
    corr_post = cycle_growth_labor["1984":final].corr()
    corrs_labor = pd.concat([corr_pre, corr_post])
    
    # Sectoral comovement
    cycle_growth_com = cycle_growth[['C', 'I', 'u_C', 'u_I', 'NC', 'NI','idle_D','idle_ND']]
    mom_growth_com = moments(100*cycle_growth_com, lab=['I', 'NI'], lags=[1])
    
    
     
    
    
    " Summarize correlations"
    corrs = np.zeros((8, 1))
    corr_mat = cycle_growth.corr()
    # Sectoral comovement
    corrs[0] = corr_mat['C']['I']
    corrs[1] = corr_mat['NC']['NI']
    corrs[2] = corr_mat['u_C']['u_I']
    corrs[3] = corr_mat['idle_ND']['idle_D']
    corrs[4] = corr_mat['p_I']['I']
    # Labor
    corrs[5] = corr_mat['u']['N']
    corrs[6] = corr_mat['u']['lab_prod']
    corrs[7] = corr_mat['u']['v']
    
    cor_dat = pd.DataFrame(corrs)
    cor_dat.index = ["Cor(C, I)", "Cor(NC, NI)", "Corr(u_C, u_I)",  "Cor(idle_ND, idle_D)", "Cor(p_I, I)",
                      "Cor(u,N)" , "Corr(u, lab_prod)", "Cor(u,v)"]
    
    autocorr_dat = np.zeros((5, 1))
    autocorr_dat[0] = cycle_growth.NC.autocorr()
    autocorr_dat[1] = cycle_growth.NI.autocorr()
    autocorr_dat[2] = cycle_growth.u.autocorr()
    autocorr_dat[3] = cycle_growth.idle_ND.autocorr()
    autocorr_dat[4] = cycle_growth.idle_D.autocorr()
    autocorr_dat = pd.DataFrame(autocorr_dat)
    autocorr_dat.index = ["Cor(N_C, N_{C,-1})", "Cor(N_I, N_{I,-1})", "Corr(u, u_{u,-1})",
                          "Corr(idle_ND, idle_{ND,-1}", "Corr(idle_D, idle_{D,-1}"]
    summ = pd.concat([cor_dat, autocorr_dat])
    print(summ.style.format(precision=2).to_latex())
    
    # # Stacked moments
    # stds = cycle_growth.std(axis=0)
    # stds = 100*stds[["Y", "C", "I", "u_C", "u_I", "p_I", "idle_ND", "idle_D"]]
    # stds = pd.DataFrame(stds)
    # stds.index = ["std(Y)", "std(C)", "std(I)", "std(u_C)", "std(u_I)", "std(p_I)", 
    #               "std(idle_ND)", "std(idle_D)"]
    corrs = np.zeros((9, 1))
    corr_mat = cycle_growth.corr()
    # Outputs
    corrs[0] = corr_mat['C']['I']
    # Inputs
    corrs[4] = corr_mat['u_C']['u_I']
    corrs[5] = corr_mat['idle_ND']['idle_D']
    # Cross relations
    corrs[1] = corr_mat['C']['u_I']
    corrs[2] = corr_mat['I']['u_C']
    # u-v relation
    corrs[2] = corr_mat['u']['lab_prod']
              
    # Relative price of investment
    corrs[6] = corr_mat['p_I']['I']
    cor_dat = pd.DataFrame(corrs)
    cor_dat.index = ["Cor(C, I)", "Cor(C, N_I)",  "Cor(C, lab_prod)", "Cor(C, util)",
                      "Cor(N_C, N_I)" , "Corr(util_D, util_ND)", "Cor(p_I, I)", "Cor(u, util_D)", "Cor(u, util_ND)"]
    autocorr_dat = np.zeros((3, 1))
    autocorr_dat[0] = cycle_growth.NC.autocorr()
    autocorr_dat[1] = cycle_growth.NI.autocorr()
    autocorr_dat[2] = cycle_growth.u.autocorr()
    autocorr_dat = pd.DataFrame(autocorr_dat)
    autocorr_dat.index = ["Cor(N_C, N_{C,-1})", "Cor(N_I, N_{I,-1})", "Corr(u, u_{u,-1})"]
    summ = pd.concat([cor_dat, autocorr_dat])
    print(summ.style.format(precision=2).to_latex())
    
    " For description: Hamilton-filtered series "
    cycle_ham = pd.concat([filter_transform(dat[x], init=init, final=final, transform_type='log',
                                        filter_type="ham", demean=False) for x in lab], axis=1)
    cycle_ham.columns = lab
    print(cycle_ham.mean())
    cycle_growth = cycle_growth - cycle_growth.mean()
    
    
    
    
    # Utilization
    # fig, ax = plt.subplots(figsize=(11, 4))
    # ax.plot(100*cycle.util_D, linestyle[0], label= "Utilization: durables", lw=2, alpha=0.7)
    # ax.plot(100*cycle.util_ND, linestyle[1], label= "Utilization: nondurables", lw=2, alpha=0.7)
    # ax.plot(100*cycle.u, linestyle[2], label="Unemployment rate", lw=2, alpha=0.7)
    # #ax.plot(100*cycle.util, linestyle[2], label= "Utilization: total industry", lw=2, alpha=0.7)
    # ax.legend(loc="upper right")
    # ax.xaxis.set_major_locator(years)
    # ax.xaxis.set_major_formatter(years_fmt)
    # ax.set_ylabel("%")
    # ax.grid(True)
    # plt.tight_layout()
    # plt.savefig("util_u_comovement.pdf")
    # plt.show()
    
    # In levels 
    # 1) Plot
    
    
        
    
        
        
      
    
