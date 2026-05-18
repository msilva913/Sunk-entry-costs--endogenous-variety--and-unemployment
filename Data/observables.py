
#######################
# Imports and Settings #
#######################

import time
import numpy as np
import pandas as pd
import requests
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


def fetch_bed_deaths_national():
    """
    Fetch BED national employment at dying establishments (dataclass=08) from
    the BLS public API v1.  Returns a quarterly pd.Series of deaths employment
    (thousands of workers), indexed by quarter-start dates.

    Series ID: BDS0000000000000000110008LQ5
      BD  = BED prefix
      S   = seasonally adjusted
      00000 = MSA 00000 (national)
      00    = state 00 (national)
      000   = county 000 (national)
      000000 = industry 000000 (all private)
      1     = unit analysis: establishments (employment reported)
      1     = data element: employment level
      00    = size class: all sizes
      08    = dataclass: Deaths (absent 4+ consecutive quarters)
      L     = rate/level: level (thousands of workers)
      Q     = periodicity: quarterly
      5     = ownership: private

    Deaths (dataclass=08) ⊆ Closings (dataclass=06): only establishments that
    do not reopen within four quarters are counted.  This is the conceptually
    clean measure of permanent exit, consistent with the Bartik instrument in
    construct_delta_instrument.py.
    """
    BLS_API = "https://api.bls.gov/publicAPI/v1/timeseries/data/"
    series_id = "BDS0000000000000000110008LQ5"
    year_chunks = [("1992", "2001"), ("2002", "2011"), ("2012", "2030")]
    rows = []
    for start_yr, end_yr in year_chunks:
        payload = {
            "seriesid" : [series_id],
            "startyear": start_yr,
            "endyear"  : end_yr,
        }
        resp   = requests.post(BLS_API, json=payload,
                               headers={"Content-Type": "application/json"},
                               timeout=60)
        resp.raise_for_status()
        result = resp.json()
        if result.get("status") != "REQUEST_SUCCEEDED":
            raise ValueError(f"BLS API error ({start_yr}-{end_yr}): "
                             f"{result.get('message', result)}")
        for s in result["Results"]["series"]:
            for obs in s["data"]:
                if obs.get("value", "-") == "-":
                    continue
                qnum  = int(obs["period"].lstrip("Q"))
                month = (qnum - 1) * 3 + 1
                rows.append({
                    "date"      : pd.Timestamp(int(obs["year"]), month, 1),
                    "deaths_emp": float(obs["value"].replace(",", "")),
                })
        time.sleep(1.5)   # respect BLS rate limit
    return (pd.DataFrame(rows)
              .set_index("date")["deaths_emp"]
              .sort_index())


##############################
# Economic Data Construction #
##############################

def construct_data(init, final, freq):
    """
    Returns a DataFrame of key economic series for estimation:
      c     : Real per capita consumption
      u     : Unemployment rate
      v     : Vacancy rate
      theta : Labor market tightness (v/u)
      jf    : Job-finding rate (Shimer continuous-time)
      lp    : Labor productivity index (2012=100)
      ls    : Labor share, nonfarm business
      s     : Aggregate separation rate (Shimer continuous-time)
      delta : Permanent establishment exit rate (BED Deaths, employment-weighted)
      w     : Real compensation index (ls * lp; index, not level)
      ba    : Per-capita high-propensity business applications (BFS)
    """
    # ------------------------------------------------------------------ #
    # Consumption and price level                                          #
    # ------------------------------------------------------------------ #
    C   = fred.get_series('PCE').resample(freq).mean()
    p_C = fred.get_series("PCEPI").resample(freq).mean().dropna()

    # ------------------------------------------------------------------ #
    # Non-institutional population (HP-filtered trend)                    #
    # HP filter (λ=10,000) smooths discrete Census benchmark-revision     #
    # jumps in CNP16OV; trend population is used for per-capita scaling.  #
    # ------------------------------------------------------------------ #
    pop = fred.get_series('CNP16OV').resample(freq).mean()
    pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]

    # Real per capita consumption
    c = C / (pop * p_C)

    # ------------------------------------------------------------------ #
    # Labor productivity                                                   #
    # PRS85006163: output per hour, nonfarm business sector, index        #
    # 2012=100.  Used as a growth-rate / second-moment target; the level  #
    # is not economically meaningful.                                      #
    # ------------------------------------------------------------------ #
    lp = fred.get_series("PRS85006163").resample(freq).mean().dropna()

    # ------------------------------------------------------------------ #
    # Unemployment rate                                                    #
    # ------------------------------------------------------------------ #
    u = fred.get_series('UNRATE').resample(freq).mean().dropna() / 100

    # ------------------------------------------------------------------ #
    # Labor force (for vacancy rate denominator)                          #
    # ------------------------------------------------------------------ #
    lf = fred.get_series('CLF16OV')

    # ------------------------------------------------------------------ #
    # Vacancy rate: spliced Help-Wanted Index + JOLTS openings            #
    # The HWI series (Barnichon 2010, AEJ Macro) provides consistent      #
    # composite vacancy data back to 1951.  The splice point (2000-11)    #
    # is a level discontinuity; for improved consistency consider         #
    # Barnichon's (2010) regression-based rescaling over the JOLTS        #
    # overlap period before concatenating.                                 #
    # ------------------------------------------------------------------ #
    df_hwi = pd.read_csv("CompositeHWI.csv")
    date   = pd.date_range(start="1951", end="2021-08", freq="MS")
    df_hwi.index = pd.DatetimeIndex(date)
    v_hwi  = df_hwi["V_hwi"]
    opens  = fred.get_series("JTSJOL", start=init, stop=final)
    v      = pd.concat([v_hwi[:"2000-11-01"], opens])
    v      = v / lf
    v      = v.resample(freq).mean().dropna()
    theta  = v / u

    # ------------------------------------------------------------------ #
    # Job-finding and separation rates via Shimer (2012) method           #
    # Continuous-time correction accounts for time-aggregation bias.      #
    #                                                                      #
    # 1994 CPS redesign correction: the January 1994 CPS questionnaire    #
    # redesign artificially increased measured short-term unemployment     #
    # (UEMPLT5) by approximately 16%.  We multiply UEMPLT5 by 1/1.16 for #
    # dates before 1994 (equivalently, multiply post-1994 implied s_h by  #
    # 1.16 if using pre-1994 as baseline — see Shimer 2012 footnote 8 and #
    # Elsby, Hobijn & Sahin 2013).  We apply the adjustment directly to   #
    # u_short so all subsequent formulas remain unchanged.                 #
    # ------------------------------------------------------------------ #
    u_level = fred.get_series('UNEMPLOY').loc[init:final]
    e_level = fred.get_series('CE16OV').loc[init:final]
    u_short = fred.get_series('UEMPLT5').loc[init:final]

    # Apply 1994 CPS redesign correction: scale up pre-1994 short-term
    # unemployment to be consistent with the post-redesign measurement.
    u_short.loc[u_short.index < '1994-01-01'] *= 1.16

    u_next  = u_level[1:len(u_level)]
    jf      = 1 - (u_next - u_short) / u_level
    jfh     = -np.log(1 - jf)
    denom   = -np.expm1(-jfh)
    sh      = (u_short * jfh) / (e_level * denom)
    s       = 1.0 - np.exp(-sh)
    s       = s.drop("2020-03-01", errors='ignore')

    # Quarterly averages
    s  = s.resample(freq).mean().dropna()
    jf = jf.resample(freq).mean().dropna()

    # ------------------------------------------------------------------ #
    # Labor share and real wage index                                      #
    # PRS85006173: labor share, nonfarm business.                         #
    # w = ls * lp is a real compensation index (index units, not level).  #
    # Appropriate for second-moment / growth-rate targets in estimation.  #
    # ------------------------------------------------------------------ #
    ls = fred.get_series('PRS85006173').resample(freq).mean()
    w  = ls * lp

    # ------------------------------------------------------------------ #
    # Business formation series                                            #
    # BAWBATOTALSAUS (Business Formation Statistics, high-propensity      #
    # business applications with planned wages) is preferred over         #
    # BFBF4QTOTALSAUS (4-quarter realized entry).  Reasons:              #
    #   1. Broader coverage and longer history (available from 2004Q3).  #
    #   2. Not subject to mechanical lags in conversion-to-entry during   #
    #      downturns, when entry is delayed but applications still occur. #
    #   3. Corresponds conceptually to the BF8 moment used in Bayesian   #
    #      estimation (high-propensity applicants ≈ latent entrants).    #
    # ------------------------------------------------------------------ #
    ba = fred.get_series("BAWBATOTALSAUS").resample(freq).mean().dropna() / pop

    # ------------------------------------------------------------------ #
    # Permanent establishment exit rate (BED Deaths, employment-weighted) #
    # BED dataclass=08 (Deaths): establishments absent for 4+ consecutive #
    # quarters — the same criterion used in construct_delta_instrument.py #
    # for the Bartik instrument.  Deaths ⊆ Closings; temporary shutdowns  #
    # are excluded.                                                        #
    #                                                                      #
    # delta = deaths_employment_t / payems_t                              #
    #   deaths_employment: employment at dying establishments (thousands), #
    #     BLS BED series BDS0000000000000000110008LQ5                     #
    #   payems: total nonfarm payroll employment (thousands), FRED PAYEMS #
    #                                                                      #
    # Model frequency note: the model is cast at monthly frequency.       #
    # BED Deaths are only available quarterly.  We represent delta at     #
    # monthly frequency by dividing the quarterly rate by 3 — a first-   #
    # order approximation equivalent to assuming deaths are distributed   #
    # uniformly within the quarter.  The same convention is applied when  #
    # computing calibration targets (delta_bar).                          #
    # ------------------------------------------------------------------ #
    deaths_emp = fetch_bed_deaths_national()
    payems     = fred.get_series('PAYEMS')          # total nonfarm, thousands
    delta_q    = (deaths_emp / payems).resample(freq).mean().dropna()
    # Convert quarterly rate to monthly equivalent (model period = 1 month)
    delta      = delta_q / 3

    # ------------------------------------------------------------------ #
    # Combine                                                              #
    # ------------------------------------------------------------------ #
    var_load_list = [c, u, v, theta, jf, lp, ls, s, delta, w, ba]
    dat = pd.concat(var_load_list, axis=1)
    dat.columns = ['c', 'u', 'v', 'theta', 'jf', 'lp', 'ls', 's', 'delta', 'w', 'ba']
    return dat

######################
# Main Script/Plots  #
######################

if __name__ == "__main__":

    # Date range and configuration
    init = '1951-01-01'
    final = '2025-10-01'
    freq = 'QS'
    load = False

    # Load or construct data
    if load:
        dat = pickle.load(open("var_load_list", "rb"))
    else:
        dat = construct_data(init, final, freq)
        save_object(dat, 'dat')

    # Extract relevant dates
    dat = dat.loc[init:final]

    # Optional: view separation series since BED Deaths begin (1992Q3)
    separation_panel = dat[["s", "delta"]].loc["1992":final]
    dat.to_pickle("raw_data.pkl")

    # Plot idiosyncratic match separation vs permanent exit rate
    # s measures total match separations; delta is establishment deaths.
    # Their difference approximates separations from continuing establishments.
    s_ind = dat.s - dat.delta
    fig, ax = plt.subplots()
    ax.plot(dat.delta.index, s_ind, label="Idiosyncratic separation (s − δ)")
    ax.plot(dat.delta.index, dat.delta, label="Permanent exit rate (δ, BED Deaths)")
    plt.legend()
    plt.savefig("separation_vs_exit.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[plot] saved: separation_vs_exit.png")
   
    
        
        
      
    
