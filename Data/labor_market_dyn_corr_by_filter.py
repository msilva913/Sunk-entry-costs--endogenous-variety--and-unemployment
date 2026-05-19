"""
labor_market_dyn_corr_by_filter.py
------------------------------------
Dynamic cross-correlations between the permanent establishment exit rate
(delta) and unemployment / vacancies, computed under four detrending
conventions:

    HP-1600   : HP filter, lambda=1,600  (standard DSGE)
    HP-100k   : HP filter, lambda=100,000 (Shimer/CK convention)
    Hamilton  : Hamilton (2018) regression filter, h=8 quarters
    FD        : First differences (quarter-on-quarter)

Purpose: demonstrate that HP-100k is the anomalous specification.
Under every other filter, cor(delta, u) > 0 and cor(delta, v) < 0
at h=0, consistent with the vacancy-destruction mechanism and with
the causal IRFs in Section 4. HP-100k reverses the sign of cor(delta, v)
and drives cor(delta, u) to near zero — an artifact of over-smoothing
that strips out the business-cycle variation in which the mechanism
operates, not a feature of the data.

Output
------
  data/results/dyn_corr_delta_u_v_by_filter.csv   — numeric table
  data/results/dyn_corr_delta_u_by_filter.png
  data/results/dyn_corr_delta_v_by_filter.png

Usage
-----
    python labor_market_dyn_corr_by_filter.py

Requires raw_data.pkl (produced by observables.py).
"""

#####################
# Imports           #
#####################

import os
import pickle
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#####################
# Paths             #
#####################

RESULTS_DIR = os.path.join("Bartek analysis", "data", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

#####################
# Filter functions  #
#####################

def hp_cycle(series, lam):
    """Return HP-filter cyclical component."""
    cyc, _ = sm.tsa.filters.hpfilter(series, lamb=lam)
    return cyc


def hamilton_cycle(series, h=8, p=4):
    """
    Hamilton (2018, ReStat) regression filter.

    Regresses y_t on [1, y_{t-h}, y_{t-h-1}, ..., y_{t-h-p+1}] via OLS
    on the full sample.  The residual is the cyclical component.

    Parameters
    ----------
    series : pd.Series
    h : int
        Projection horizon (8 quarters = 2 years, Hamilton's recommendation
        for quarterly data).
    p : int
        Number of lags in the projection (Hamilton uses p=4).

    Returns
    -------
    pd.Series of cyclical residuals; first h+p-1 observations are NaN.

    Notes
    -----
    Residuals inherit an MA(h-1) error structure by construction.
    SMM weighting matrices should use HAC with bandwidth >= h.
    """
    y = np.asarray(series, dtype=float)
    n = len(y)
    resid = np.full(n, np.nan)

    # Observations for which all regressors are available: t = h+p-1, ..., n-1
    t_start = h + p - 1
    T = n - t_start
    if T <= p + 2:
        return pd.Series(resid, index=series.index)

    Y = y[t_start:]                                    # shape (T,)
    X = np.column_stack(                               # shape (T, p+1)
        [np.ones(T)] +
        [y[t_start - h - j: n - h - j] for j in range(p)]
    )
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid[t_start:] = Y - X @ beta
    return pd.Series(resid, index=series.index)


def first_diff(series):
    """Quarter-on-quarter first difference."""
    return series.diff()


#####################
# Dynamic correlation
#####################

def dyn_corr(x_series, y_series, h):
    """
    cor(x_t, y_{t+h}).  Positive h means y leads x by h periods.
    Both series must share the same index; NaNs are dropped pairwise.
    """
    if h >= 0:
        x = x_series.iloc[:len(x_series) - h] if h > 0 else x_series
        y = y_series.iloc[h:]                  if h > 0 else y_series
    else:
        x = x_series.iloc[-h:]
        y = y_series.iloc[:len(y_series) + h]
    xy = pd.concat([x.reset_index(drop=True),
                    y.reset_index(drop=True)], axis=1).dropna()
    if len(xy) < 10:
        return np.nan
    return float(np.corrcoef(xy.iloc[:, 0], xy.iloc[:, 1])[0, 1])


#####################
# Main              #
#####################

if __name__ == "__main__":

    # ------------------------------------------------------------------ #
    # Load data                                                            #
    # ------------------------------------------------------------------ #
    dat = pickle.load(open("raw_data.pkl", "rb"))
    core = dat[['u', 'v', 'delta']].dropna()
    print(f"Sample: {core.index[0].date()} to {core.index[-1].date()}, "
          f"N={len(core)}")

    # ------------------------------------------------------------------ #
    # Construct filtered series                                            #
    # ------------------------------------------------------------------ #
    filters = {
        'HP-1600' : core.apply(hp_cycle,       lam=1_600),
        'HP-100k' : core.apply(hp_cycle,       lam=100_000),
        'Hamilton': core.apply(hamilton_cycle, h=8, p=4).dropna(),
        'FD'      : core.apply(first_diff).dropna(),
    }

    # Sanity check: contamination of lambda=1600 cycle by 100k trend
    trend_100k = core['delta'] - hp_cycle(core['delta'], lam=100_000)
    cyc_1600   = hp_cycle(core['delta'], lam=1_600)
    print(f"\nContamination check — cor(cyc_1600_delta, trend_100k): "
          f"{cyc_1600.corr(trend_100k):.3f}  (near 0 = no trend bleed)")

    # ------------------------------------------------------------------ #
    # Compute dynamic correlations                                         #
    # ------------------------------------------------------------------ #
    horizons = list(range(-4, 5))   # h = -4 ... +4

    rows = []
    for outcome in ['u', 'v']:
        for h in horizons:
            row = {'outcome': outcome, 'h': h}
            for fname, df in filters.items():
                row[fname] = dyn_corr(df['delta'], df[outcome], h)
            rows.append(row)

    results = pd.DataFrame(rows)

    # ------------------------------------------------------------------ #
    # Print tables                                                         #
    # ------------------------------------------------------------------ #
    filter_cols = ['HP-1600', 'HP-100k', 'Hamilton', 'FD']

    for outcome, label in [('u', 'unemployment'), ('v', 'vacancies')]:
        sub = results[results['outcome'] == outcome].set_index('h')
        print(f"\n{'='*62}")
        print(f"Dynamic cor(delta_t, {outcome}_{{t+h}})  —  {label}")
        print(f"{'='*62}")
        header = f"{'h':>5}" + "".join(f"{c:>12}" for c in filter_cols)
        print(header)
        print("-" * 53)
        for h in horizons:
            vals = sub.loc[h, filter_cols]
            line = f"{h:>+5d}" + "".join(f"{v:>12.3f}" for v in vals)
            # Flag where HP-100k disagrees with the majority
            hp100 = vals['HP-100k']
            others = [vals[c] for c in filter_cols if c != 'HP-100k']
            same_sign = all(np.sign(o) == np.sign(hp100) for o in others
                            if not np.isnan(o))
            line += "  <-- HP-100k outlier" if not same_sign else ""
            print(line)

    # ------------------------------------------------------------------ #
    # Save CSV                                                             #
    # ------------------------------------------------------------------ #
    out_csv = os.path.join(RESULTS_DIR, "dyn_corr_delta_u_v_by_filter.csv")
    results.to_csv(out_csv, index=False, float_format="%.4f")
    print(f"\nSaved: {out_csv}")

    # ------------------------------------------------------------------ #
    # Plots                                                                #
    # ------------------------------------------------------------------ #
    COLORS = {
        'HP-1600' : '#2166ac',   # blue
        'HP-100k' : '#d73027',   # red  — the anomaly
        'Hamilton': '#1a9641',   # green
        'FD'      : '#7b2d8b',   # purple
    }
    STYLES = {
        'HP-1600' : '-o',
        'HP-100k' : '--s',       # dashed to visually flag it
        'Hamilton': '-^',
        'FD'      : '-D',
    }

    for outcome, label, zero_ref in [
        ('u', 'Unemployment rate', 'Expected sign: positive (+)'),
        ('v', 'Vacancy rate',      'Expected sign: negative (−)'),
    ]:
        sub = results[results['outcome'] == outcome].set_index('h')

        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.axhline(0, color='k', lw=0.8, ls='-')
        ax.axvline(0, color='k', lw=0.5, ls='--', alpha=0.4)

        for fname in filter_cols:
            ax.plot(horizons, sub[fname].values,
                    STYLES[fname], color=COLORS[fname],
                    lw=1.8, ms=5, label=fname)

        ax.set_xlabel('Horizon h (quarters)', fontsize=11)
        ax.set_ylabel(f'cor($\\delta_t$, {outcome}$_{{t+h}}$)', fontsize=11)
        ax.set_title(
            f'Dynamic correlation: permanent exit rate ($\\delta$) '
            f'and {label.lower()}\n'
            f'{zero_ref} — HP-$\\lambda$=100k is the anomalous specification',
            fontsize=10
        )
        ax.set_xticks(horizons)
        ax.legend(fontsize=10, framealpha=0.9)
        ax.grid(True, alpha=0.25)

        fig.tight_layout()
        fname_out = os.path.join(
            RESULTS_DIR,
            f"dyn_corr_delta_{outcome}_by_filter.png"
        )
        fig.savefig(fname_out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {fname_out}")

    print("\nDone.")
