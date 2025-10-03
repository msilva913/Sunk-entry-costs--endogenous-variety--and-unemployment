
#######################
# Imports and Settings #
#######################

import numpy as np
import pandas as pd
from fredapi import Fred
import pickle
import statsmodels.api as sm
import matplotlib.pyplot as plt


# Display and precision settings
pd.set_option('display.precision', 3)
np.set_printoptions(precision=3)

###########################
# FRED/API Configuration  #
###########################

fred = Fred(api_key='d35aabd7dc07cd94481af3d1e2f0ecf3')

####################
# Utility Function #
####################

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

##############################
# Economic Data Construction #
##############################

def construct_data(init, final, freq):
    """
    Returns a list of key economic series for estimation:
    c: Real per capita consumption
    u: Unemployment rate
    v: Vacancy rate
    theta: Labor market tightness
    lp: Labor productivity
    s: Aggregate separation rate
    w: Real hourly compensation for all workers
    SBF4: Business formations within four quarters
    WBA: Planned wage business applications
    delta: Establishment exit rate
    """
    # Consumption and Price levels
    C = fred.get_series('PCE').resample(freq).mean()
    p_C = fred.get_series("PCEPI").resample(freq).mean().dropna()

    # Non-institutional population (smoothed using HP filter)
    pop = fred.get_series('CNP16OV').resample(freq).mean()
    pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]

    # Real per capita consumption
    c = C / (pop * p_C)

    # Labor productivity
    lp = fred.get_series("PRS85006163").resample(freq).mean().dropna()

    # Unemployment rate
    u = fred.get_series('UNRATE').resample(freq).mean().dropna() / 100

    # Labor force
    lf = fred.get_series('CLF16OV')

    # Vacancy rate: Composite Help-Wanted Index & JOLTS openings
    df = pd.read_csv("CompositeHWI.csv")
    date = pd.date_range(start="1951", end="2021-08", freq="MS")
    df.index = pd.DatetimeIndex(date)
    v_hwi = df["V_hwi"]
    opens = fred.get_series("JTSJOL", start=init, stop=final)
    v = pd.concat([v_hwi[:"2000-11-01"], opens])
    v = v / lf
    v = v.resample(freq).mean().dropna()
    theta = v / u

    # Job finding/separation rates via Shimer method
    u_level = fred.get_series('UNEMPLOY').loc[init:final]
    e_level = fred.get_series('CE16OV').loc[init:final]
    u_short = fred.get_series('UEMPLT5').loc[init:final]
    u_next = u_level[1:len(u_level)]
    jf = 1 - (u_next - u_short) / u_level
    jfh = -np.log(1 - jf)
    denom = -np.expm1(-jfh)
    sh = (u_short * jfh) / (e_level * denom)
    s = 1.0 - np.exp(-sh)
    s = s.drop("2020-03-01")
    
    # Quarterly averages
    s = s.resample(freq).mean().dropna()
    f = jf.resample(freq).mean().dropna()

    # Labor share & wage
    ls = fred.get_series('PRS85006173').resample(freq).mean()
    w = ls * lp

    # Business formation series
    sbf4 = fred.get_series('BFBF4QTOTALSAUS').resample(freq).mean().dropna() / pop
    bawba = fred.get_series("BAWBATOTALSAUS").resample(freq).mean().dropna() / pop

    # BED: Establishment exit rate
    BED_dat = pd.read_excel('BED_data.xlsx', sheet_name='Data Import')
    BED_dat['Series'] = pd.date_range(start="1992Q3", end="2024Q3", freq="QS")
    BED_dat.set_index("Series", inplace=True)
    BED_dat.rename(
        columns={"BDS0000000000000000120008RQ5": "estabs_exit_rate"},
        inplace=True
    )
    delta = BED_dat["estabs_exit_rate"].resample(freq).mean().dropna() / (100 * 3)

    # Series output
    var_load_list = [c, u, v, theta, f, lp, ls, s, delta, w, sbf4, bawba]
    return var_load_list

######################
# Main Script/Plots  #
######################

if __name__ == "__main__":

    # Date range and configuration
    init = '1951-01-01'
    final = '2025-07-30'
    freq = 'QS'
    load = False

    # Load or construct data
    if load:
        var_load_list = pickle.load(open("var_load_list", "rb"))
    else:
        var_load_list = construct_data(init, final, freq)
        save_object(var_load_list, 'var_load_list')

    # Combine into DataFrame
    dat = pd.concat(var_load_list, axis=1)
    dat.columns = ['c', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'delta', 'w', 'bf', 'ba']
    dat = dat.loc[init:final]

    # Optional: View separation series since 1992
    separation_panel = dat[["s", "delta"]].loc["1992":final]
    dat.to_pickle("raw_data.pkl")

    # Plot idiosyncratic separation vs establishment exit rate
    s_ind = dat.s - dat.delta
    fig, ax = plt.subplots()
    ax.plot(dat.delta.index, s_ind, label="Idiosyncratic separation")
    ax.plot(dat.delta.index, dat.delta, label="Establishment exit rate")
    plt.legend()
    plt.show()  
   
    
        
        
      
    
