
#######################
# Imports and Settings #
#######################

import io
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


# ------------------------------------------------------------------------ #
# 1994 CPS redesign correction (fixed Oct 28, 2026).                        #
# Since January 1994 the CPS no longer asks the unemployment duration of a  #
# worker unemployed in consecutive months; the BLS adds the elapsed weeks   #
# instead. Measured short-term unemployment (UEMPLT5) therefore FELL        #
# discontinuously in the six rotation groups that are not incoming.         #
# Shimer, "Reassessing the Ins and Outs of Unemployment" (working-paper     #
# version, Appendix A and its footnote 27): short-term unemployment was     #
# 41.6% of unemployment in the full CPS and 41.7% in the incoming rotation  #
# groups over 1976-93, but 37.9% and 44.2% over 1994-2005. He keeps the     #
# full sample through 1993 and corrects the data AFTER 1994; multiplying    #
# UEMPLT5 by a constant 1.1 after 1994 "delivers very similar results".     #
# We use that constant. (The code before Oct 28 multiplied PRE-1994 values  #
# by 1.16, which doubled the break instead of removing it.)                 #
# Published version: Review of Economic Dynamics 15(2), 2012, 127-148,      #
# doi:10.1016/j.red.2012.02.001; check page and footnote numbers there.     #
# ------------------------------------------------------------------------ #
CPS_REDESIGN_DATE   = '1994-01-01'
CPS_REDESIGN_FACTOR = 1.1


def shimer_flows(u_level, e_level, u_short, freq):
    """
    Monthly job-finding and separation probabilities by the Shimer (2012)
    continuous-time method, then quarterly averages. Inputs are the monthly
    FRED levels UNEMPLOY, CE16OV and UEMPLT5 (thousands). Returns (s, jf):
    s is the total separation probability tau used downstream.

    The redesign correction is applied to u_short before any formula.
    """
    u_short = u_short.copy()
    u_short.loc[u_short.index >= CPS_REDESIGN_DATE] *= CPS_REDESIGN_FACTOR

    # Flows are dated by the period [t, t+1) they cover, so both use the
    # stocks observed at t+1. Shimer eq. (4): F_t = 1 - (u_{t+1} - u^s_{t+1})/u_t.
    # Separation rate from eq. (2) with e held at e_t over the period:
    # u^s_{t+1} = e_t s_t (1 - exp(-f_t))/f_t.
    # (Fixed Oct 28, 2026: the old u_level[1:] kept its date labels, so pandas
    # aligned it with u_t and jf_t collapsed to u^s_t/u_t.)
    u_next  = u_level.shift(-1)
    us_next = u_short.shift(-1)
    jf      = 1 - (u_next - us_next) / u_level
    jfh     = -np.log(1 - jf)
    denom   = -np.expm1(-jfh)
    sh      = (us_next * jfh) / (e_level * denom)
    s       = 1.0 - np.exp(-sh)
    s       = s.drop(pd.Timestamp("2020-03-01"), errors='ignore')

    return s.resample(freq).mean().dropna(), jf.resample(freq).mean().dropna()


def fetch_fred_csv(series_id):
    """A FRED series from the public fredgraph CSV endpoint (no API key)."""
    url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={series_id}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text), index_col=0, parse_dates=True)
    return pd.to_numeric(df.iloc[:, 0], errors='coerce').dropna()


def refresh_flows(path, init, final, freq):
    """
    Rebuild option (b), Oct 28, 2026: recompute only s (tau) and jf in an
    existing raw_data.pkl from UNEMPLOY, CE16OV and UEMPLT5, leaving every
    other series at its earlier vintage. Prints the change.
    """
    dat = pd.read_pickle(path)
    s, jf = shimer_flows(*(fetch_fred_csv(k).loc[init:final]
                           for k in ('UNEMPLOY', 'CE16OV', 'UEMPLT5')), freq)
    for name, new in (('s', s), ('jf', jf)):
        old = dat[name].copy()
        dat[name] = new.reindex(dat.index)
        both = pd.concat([old, dat[name]], axis=1).dropna()
        for lo, hi in (('1951', '1993'), ('1994', '2019')):
            b = both.loc[lo:hi]
            print(f"  {name:2s} {lo}-{hi}: mean {b.iloc[:, 0].mean():.4f} -> "
                  f"{b.iloc[:, 1].mean():.4f}")
    return dat


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
    # PRS85006163: output per worker, nonfarm business sector, index      #
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
    # (see shimer_flows for the formulas and the 1994 CPS correction)     #
    # ------------------------------------------------------------------ #
    s, jf = shimer_flows(fred.get_series('UNEMPLOY').loc[init:final],
                         fred.get_series('CE16OV').loc[init:final],
                         fred.get_series('UEMPLT5').loc[init:final], freq)

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
    # Oct 28, 2026 (MS, rebuild option (b)): recompute only s and jf in the
    # existing raw_data.pkl after the CPS-redesign fix. Set to False for a full
    # rebuild, which re-downloads every series at the current vintage.
    flows_only = True

    # Load or construct data
    if flows_only:
        dat = refresh_flows("raw_data.pkl", init, final, freq)
    elif load:
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
   
    
        
        
      
    
