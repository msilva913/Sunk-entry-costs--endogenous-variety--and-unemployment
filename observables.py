
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


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)  


def construct_data(init, final, freq):
    """
    Variables for estimation
    c: real per capita consumption 
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
    
    " BED "
    " Gross job gains and gross job losses from 1992 forward "
    BED_dat = pd.read_excel('BED_data.xlsx', sheet_name='Data Import')
    BED_dat['Series'] = pd.date_range(start="1992Q3", end="2024Q3", freq="QS")
    BED_dat.set_index("Series", inplace=True)
    
    #series_id = "BDS0000000000000000110008RQ5"
    series_id = "BDS0000000000000000120008RQ5"
    BED_dat.rename(columns={series_id: "estabs_exit_rate"}, inplace=True)
    delta = BED_dat["estabs_exit_rate"].resample(freq).mean().dropna()/(100*3)

    " Note: these series imply labor productivity in each sector "
    " List of data series "
    var_load_list = [c, u, v, theta, f, lp, ls, s, delta, w, sbf4, bawba] 
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
    lab = ['c', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'delta', 'w', 'bf', 'ba']
    dat.columns = lab
    dat = dat.loc[init:final]
    
    dat[["s", "delta"]].loc["1992":"2019"]
    
    s_ind = dat.s - dat.delta
    fig, ax = plt.subplots()
    ax.plot(dat.delta.index, s_ind, label="Idiosyncratic separation")
    ax.plot(dat.delta, label="Establishment exit rate")
    plt.legend()
    plt.show()
    dat.to_pickle("raw_data.pkl")
    #dat = pd.read_pickle("raw.pkl")
    
  
   
    
        
        
      
    
