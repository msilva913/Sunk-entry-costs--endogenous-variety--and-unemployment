# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 09:33:25 2019

@author: Administrator
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.arima_model import ARMA
import statsmodels.api as sm

def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation. 
    Shifted data filled with NaNs 
    
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else: 
        return datax.corr(datay.shift(lag))

def plot_solution(dic_x, dic_y, sim_length=None, type='', suptitle ='', bins=20, freq='Quarterly'):
    " Automatically adjust subplots to number of variables "
    nvars = len(dic_y)
    if nvars < 4: 
        m, n = 1,3
    if nvars == 4:
        m, n = 2, 2
    if 5 <= nvars <= 6:
        m, n = 3, 2
    if 7 <= nvars <= 9:
        m, n = 3, 3
    if 10 <= nvars <= 12:
        m, n = 3, 4
    if 13 <= nvars <= 16:
        m, n = 4, 4
    fig = plt.figure(figsize=(12, 12))
    
    for index_x, key_x in enumerate(dic_x):
        for index, key in enumerate(dic_y):
            k, l = np.unravel_index(index, (m, n))
            ax = fig.add_subplot(m, n, index+1)
            if type == 'plot':
                ax.plot(dic_x[key_x], dic_y[key], label = key, lw=2)
                ax.set_xlabel(key_x, fontsize=12)
            elif type== 'scatter':
                ax.scatter(dic_x[key_x], dic_y[key], label = key, lw=2)
                ax.set_xlabel(key_x, fontsize=12)
            elif type== 'hist':
                ax.hist(dic_y[key], label = key, lw=2, bins=bins)
                ax.set_xlabel('Value', fontsize=12)
            else:
                if sim_length is None:
                    y = dic_y[key]
                else:
                    y = dic_y[key][0:sim_length]
                T = range(len(y))
                ax.plot(T, y, label = key, lw=2)
                ax.set_xlabel(freq, fontsize=12)
            ax.legend(loc='best')
            fig.suptitle(suptitle, fontsize=18)
    plt.show()


def time_series_plot(dataframe, labels, nrows=2, ncols=3, figsize=(14, 10), loc='upper center', 
                     fontsize=12):
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize= figsize)
    for index, x in enumerate(dataframe.columns):
        k, l = np.unravel_index(index, (nrows, ncols))
        ax[k,l].plot(dataframe[x], label= labels[index], lw=2)
        ax[k,l].legend(loc=loc, frameon=True, fontsize=fontsize)
    plt.show()


def fit_AR(x, p):
    mod = ARMA(x, order=(p, 0))
    res = mod.fit(trend='nc')
    rho = res.params[0]
    sigma = np.std(res.resid)
    return rho, sigma


def moments(dat, relative_std='Y', lab=['Y'], lags=[1,2]):
    
    def df_autocorr(df, lag, axis=0):
        """Compute full-sample column-wise autocorrelation for series"""
        return df.apply(lambda x: x.autocorr(lag), axis=axis)
    
    """
    Generates table of moments
    inputs: data, date range (init, final), index labels for variables, transformation_type, forecast gap h in Hamilton Method
    output: dataframe of moments
    Second, compute moments
    Output both dataframe and TeX code
    """
    
    corr_labels = []
    for lab_type in lab:
        corr_labels.append(f'Cor(x, {lab_type})')
        
    moments = np.zeros((len(dat.columns), 2+len(lab)+len(lags)))
    moments[:, 0] = dat.std()
    moments[:, 1] = dat.std()/dat.std()[relative_std]
    moments[:, 2:(2+len(lab))] = dat.corr()[lab]
    k = 2+len(lab)
    for lag in lags:
        moments[:, k] = np.asarray(df_autocorr(dat, lag))
        k += 1

    columns = ('SD(x)', 'RSD', *corr_labels)
    for lag in lags:
        columns = columns + (f'Cor(x,x_{-lag})',)
       
    moments = pd.DataFrame(moments, columns=columns, index=dat.columns)
    print(moments.to_string())
    return moments


def moments_dynamic(dat, dyn_lag=[4, -4]):
    
    " Moments augmented to include dynamic correlations "
    def df_autocorr(df, lag=1, axis=0):
        """Compute full-sample column-wise autocorrelation for series"""
        return df.apply(lambda x: x.autocorr(lag), axis=axis)
    
    dat_mod = dat.copy()
    " Lag 4 periods "
    for i, lag in enumerate(dyn_lag):
        y_shift = dat.iloc[:, 0].shift(lag)
        dat_mod.iloc[:, i] = y_shift

    """
    Generates table of moments
    inputs: data, date range (init, final), index labels for variables, transformation_type, forecast gap h in Hamilton Method
    output: dataframe of moments
    First, transform function and apply Hamilton filter
    Second, compute moments
    Output both dataframe and TeX code
    """
    # Log and Hamilton-filter u, d, V, p
    moments = np.vstack((dat.std(), dat.std()/dat.std()[0], dat.corr()[0:1],
                         dat_mod.corr()[0:len(dyn_lag)], df_autocorr(dat, lag=1))).T
    
    columns=('SD(x)', 'RSD', 'Cor(x,' + dat.columns[0]+')')
    
    for lag in dyn_lag:
        columns = columns + ('Cor(x,' + dat.columns[0]+f'_{-lag})',)
    
    columns = columns + ('Cor(x,x_{-1})',)
    
    moments = pd.DataFrame(moments, columns=columns, index=dat.columns)
    print(moments.to_latex())
    return moments

def hamilton_filter(x, h=8):
    """
    Applies the Hamilton regression filter to the input data.

    Parameters:
    x: DataFrame or Series that can be coerced into a DataFrame.
    h: Forecast horizon (default=8).

    Returns:
    residuals: Residuals of the linear projection.
    predictions: Predicted values from the regression.
    """
    # Convert input to DataFrame if it's not already
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame()

    # Lagged values
    x_h = x.shift(h)

    # Prepare the design matrix with a column of ones and lagged values
    X = pd.DataFrame(np.ones(len(x)), index=x.index, columns=['const'])  # Column of ones for intercept
    X = pd.concat([X, x_h, x_h.shift(1), x_h.shift(2), x_h.shift(3)], axis=1)
    X.columns = ['const', 'x_h', 'x_h-1', 'x_h-2', 'x_h-3']

    # Fit the model using OLS
    reg = sm.OLS(x, X, missing='drop')
    results = reg.fit()

    # Return residuals and predicted values
    #return results.resid, results.predict(X)
    return results.resid

def growth_filter(x):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    
    x = x.diff()
    #x = x-x.mean()
    return x


def linear_filter(x):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    T = np.arange(0, len(x.index), 1)
    X = pd.DataFrame(np.ones(len(x)), index=x.index)
    X['time'] = T
    reg = sm.OLS(x, exog=X, missing='drop')
    results = reg.fit()
    return results.resid


def linear_quadratic_filter(x):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    T = np.arange(0, len(x.index), 1)
    T2 = T**2
    X = pd.DataFrame(np.ones(len(x)), index=x.index)
    X['time'] = T
    X['time2'] = T2
    reg = sm.OLS(x, exog=X, missing='drop')
    results = reg.fit()
    return results.resid

def cubic_filter(x):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    T = np.arange(0, len(x.index), 1)
    T2 = T**2
    T3 = T**3
    X = pd.DataFrame(np.ones(len(x)), index=x.index)
    X['time'] = T
    X['time2'] = T2
    X['time3'] = T3
    reg = sm.OLS(x, exog=X, missing='drop')
    results = reg.fit()
    return results.resid


def hp_filter(x, lamb=1600):
    # Hodrick-Prescott filter
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    cycle, trend = sm.tsa.filters.hpfilter(x, lamb)
    return cycle

def bk_filter(x, P_L=6, P_H=32):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    cycle = sm.tsa.filters.bkfilter(x, P_L, P_H)
    return cycle

def cf_filter(x, P_L=6, P_H=32):
    if not isinstance(x, pd.DataFrame):
        x = pd.Series(x)
        x = x.to_frame() #convert to dataframe if not already
    cycle = sm.tsa.filters.cffilter(x, low=P_L, high=P_H)[0]
    return cycle

def filter_transform(x, init, final, transform_type = None, freq=None,
                     filter_type='hamilton', h=8, lamb=1600, P_L=6,
                     P_H=32, demean=False):
    """
    x: input series (pandas series)
    init: first data
    final: final data
    transform_type: transformation of data (log)
    h: forecast horizon (default=8)
    """
    #Transform to desired frequency
    if freq:
        x = x.resample(freq).mean()
    " Return non-null entries "
    x = x[x.notnull()]
    x = x[init:final]
    #Normalizing transformation of data (log, proportional deviations, or identity)
    if transform_type == 'log':
        transform = np.log
    elif transform_type == 'prop_dev':
        transform = lambda x: (x-np.mean(x))/np.mean(x)
    else:
        transform = lambda x: x
    
    z = transform(x)
    " Choice of filter "
    if filter_type == 'hamilton':
        cycle = hamilton_filter(z, h=h)
    elif filter_type == 'growth':
        cycle = growth_filter(z)
    elif filter_type == 'linear':
        cycle = linear_filter(z)
    elif  filter_type == 'quadratic':
        cycle = linear_quadratic_filter(z)
    elif  filter_type == 'cubic':
        cycle = cubic_filter(z)
    elif filter_type == 'hp_filter':
        cycle = hp_filter(z, lamb=lamb)
    elif filter_type == 'bkfilter':
        cycle = bk_filter(z, P_L=P_L, P_H=P_H)
    else:
        Exception('filter type is not found')
    " Subset data range "
    cycle = cycle.loc[init:final]
    " Demean (since using subset of data) "
    if demean:
        cycle = cycle-cycle.mean()
    
    return cycle


def monthly_to_quarterly(x):
    x = x.resample('Q').mean()
    x = 1 - (1-x)**3
    return x

def monthly_to_quarterly_datetime(sim, N=10000):
    
    capT= sim.shape[0]
    K = int(capT/N)
    date = pd.date_range(start='01-01-01', periods=K, freq='MS')
    df_array = pd.DataFrame()
    for i in range(N):
        df = sim.iloc[K*i:K*(i+1), :]
        df.index = date
        df = df.resample('QS').mean()
        df_array = pd.concat([df_array, df])
    
    return df_array

def quarterly_to_annual_datetime(sim, N=100):
    
    capT= sim.shape[0]
    K = int(capT/N)
    date = pd.date_range(start='01-01-01', periods=K, freq='QS')
    df_array = pd.DataFrame()
    for i in range(N):
        df = sim.iloc[K*i:K*(i+1), :]
        df.index = date
        df = df.resample('AS').mean()
        df_array = pd.concat([df_array, df])
    
    return df_array


def qnwmonomial1(vcv):
    n = vcv.shape[0]
    n_nodes = 2*n

    z1 = np.zeros((n_nodes, n))

    # In each node, random variable i takes value either 1 or -1, and
    # all other variables take value 0. For example, for N = 2,
    # z1 = [1 0; -1 0; 0 1; 0 -1]
    for i in range(n):
        z1[2*i:2*(i+1), i] = [1, -1]

    sqrt_vcv = np.linalg.cholesky(vcv)
    R = np.sqrt(n)*sqrt_vcv
    ϵj = z1 @ R
    ωj = np.ones(n_nodes) / n_nodes

    return ϵj, ωj


def qnwmonomial2(vcv):
    n = vcv.shape[0]
    assert n == vcv.shape[1], "Variance covariance matrix must be square"
    z0 = np.zeros((1, n))

    z1 = np.zeros((2*n, n))
    # In each node, random variable i takes value either 1 or -1, and
    # all other variables take value 0. For example, for N = 2,
    # z1 = [1 0; -1 0; 0 1; 0 -1]
    for i in range(n):
        z1[2*i:2*(i+1), i] = [1, -1]

    z2 = np.zeros((2*n*(n-1), n))
    i = 0

    # In each node, a pair of random variables (p,q) takes either values
    # (1,1) or (1,-1) or (-1,1) or (-1,-1), and all other variables take
    # value 0. For example, for N = 2, `z2 = [1 1; 1 -1; -1 1; -1 1]`
    for p in range(n-1):
        for q in range(p+1, n):
            z2[4*i:4*(i+1), p] = [1, -1, 1, -1]
            z2[4*i:4*(i+1), q] = [1, 1, -1, -1]
            i += 1

    sqrt_vcv = np.linalg.cholesky(vcv)
    R = np.sqrt(n+2)*sqrt_vcv
    S = np.sqrt((n+2)/2)*sqrt_vcv
    ϵj = np.row_stack([z0, z1 @ R, z2 @ S])

    ωj = np.concatenate([2/(n+2) * np.ones(z0.shape[0]),
                         (4-n)/(2*(n+2)**2) * np.ones(z1.shape[0]),
                         1/(n+2)**2 * np.ones(z2.shape[0])])
    return ϵj, ωj 


def dynamic_correlations(cycle, var, ylabel, nleads=12, nlags=12, title=None,
                         savefig=None):
    fig, ax = plt.subplots(figsize=(14, 5))
    # Loop over variables
    x = cycle[var[0]]
    y = cycle[var[1]]
    #rs = pd.Series([x.corr(y.shift(-lag))] for lag in range(-nlags, nleads))
    rs = pd.Series([crosscorr(cycle[var[0]], cycle[var[1]], -lag) for lag in range(-nlags, nleads)])
    ax.axhline(y=0.0, color="black", linestyle="--")
    ax.plot(rs, '-', alpha=0.7, linewidth=2.0)
    ax.axvline(np.argmax(abs(rs)),linestyle='--', color='red')
    ax.set_xlabel(r'$\Delta$ (Quarterly)', fontsize=14)
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

def stacked_moments(cycle, var_labels):
    cycle = cycle[var_labels]
    # Calculate standard deviations
    stds = cycle.std(axis=0)
    stds = pd.DataFrame(stds)
    stds.index = [f"std({label})" for label in var_labels]
    
    # Calculate correlation matrix and extract unique correlations
    corr_matrix = cycle.corr().values
    corr_array = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
    
    # Generate correlation labels
    cor_labels = []
    for i in range(len(var_labels)):
        for j in range(i+1, len(var_labels)):
            cor_labels.append(f"Cor({var_labels[i]}, {var_labels[j]})")
    
    cor_dat = pd.DataFrame(corr_array, index=cor_labels)
    
    # Calculate autocorrelations
    autocorr_dat = pd.DataFrame([cycle[x].autocorr() for x in var_labels])
    autocorr_dat.index = [f"Cor({label}, {label}_{{-1}})" for label in var_labels]
    
    # Combine all into a single DataFrame
    summ = pd.concat([stds, cor_dat, autocorr_dat])
    
    # Print formatted output
    print(summ.style.format(precision=3).to_latex())
    
    return summ

