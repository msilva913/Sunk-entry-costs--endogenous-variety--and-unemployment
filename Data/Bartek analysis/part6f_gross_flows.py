"""
part6f_gross_flows.py -- separation and job-finding probabilities from CPS gross flows.

The baseline s and jf come from Shimer's duration method (observables.py, shimer_flows): both
are built from short-term unemployment u^s_{t+1}, so an error in u^s moves them together, and
inflows from non-participation are counted as separations (Fujita and Ramey 2009, p. 420).
This program builds the alternatives from the BLS labor force status flows (CPS, matched
households, monthly, seasonally adjusted, February 1990 on), as Fujita and Ramey (2009) do:
    EU_t = (E->U)_t / (E->E + E->U + E->N)_t      employed at t-1 who are unemployed at t
    UE_t = (U->E)_t / (U->E + U->U + U->N)_t      unemployed at t-1 who are employed at t
Discrete monthly probabilities, quarterly means (as observables.py averages its monthly rates).
Marginal inflows (entrants to the population) are left out of the denominators. No continuous-
time correction for time aggregation (the model is monthly and discrete).

Fetch: BLS series LNS17000000-LNS17900000 from DBnomics (api.db.nomics.world, dataset BLS/ln);
the BLS site and its anonymous API are not reachable from the cloud container. Cached once in
data/cache/bls_gross_flows_monthly.parquet; delete the file to refetch.

Output: data/cache/gross_flows_quarterly.parquet (tau_gf = EU, jf_gf = UE, quarterly).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import requests

CACHE_M = Path("data/cache/bls_gross_flows_monthly.parquet")
OUT_Q   = Path("data/cache/gross_flows_quarterly.parquet")
SERIES  = {"EE": "LNS17000000", "UE": "LNS17100000", "NE": "LNS17200000", "EU": "LNS17400000",
           "UU": "LNS17500000", "NU": "LNS17600000", "EN": "LNS17800000", "UN": "LNS17900000"}
URL     = "https://api.db.nomics.world/v22/series/BLS/ln/{}?observations=1&format=json"


def fetch_one(code):
    """One monthly series (thousands) from DBnomics; 'NA' entries become NaN."""
    resp = requests.get(URL.format(code), timeout=60)
    resp.raise_for_status()
    doc = resp.json()["series"]["docs"][0]
    values = pd.to_numeric(pd.Series(doc["value"]), errors="coerce").to_numpy()
    return pd.Series(values, index=pd.PeriodIndex(doc["period"], freq="M"), name=code)


def monthly_flows():
    if CACHE_M.exists():
        out = pd.read_parquet(CACHE_M)
        out.index = pd.PeriodIndex(out.index, freq="M")
        return out
    out = pd.DataFrame({name: fetch_one(code) for name, code in SERIES.items()})
    stored = out.copy()
    stored.index = stored.index.astype(str)
    stored.to_parquet(CACHE_M)
    return out


if __name__ == "__main__":
    m = monthly_flows()
    prob = pd.DataFrame({
        "tau_gf": m["EU"] / (m["EE"] + m["EU"] + m["EN"]),
        "jf_gf":  m["UE"] / (m["UE"] + m["UU"] + m["UN"]),
    })
    q = prob.groupby(prob.index.asfreq("Q")).mean()
    counts = prob.groupby(prob.index.asfreq("Q")).count()
    q = q[(counts == 3).all(axis=1)]                      # complete quarters only
    stored = q.copy()
    stored.index = stored.index.astype(str)
    stored.to_parquet(OUT_Q)
    print(f"Monthly flows {m.index[0]}-{m.index[-1]}; quarterly {q.index[0]}-{q.index[-1]} (n {len(q)})")
    print(q.describe().round(4).to_string())
    print(f"Saved: {OUT_Q}")
