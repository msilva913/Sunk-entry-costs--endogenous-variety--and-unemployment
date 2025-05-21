import pandas as pd
pd.set_option('display.max_columns', 10) 
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from fredapi import Fred


# Data
#download at https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html
#df = pd.read_csv('BDS_Extension.csv')
df = pd.read_csv('bds2022.csv')
fred = Fred(api_key="8d302fb5f121b2be4f7d6e795194bcde")
df['year'] = pd.to_datetime(df['year'], format='%Y')
df.set_index("year", inplace=True)

y = fred.get_series('GDPC1').resample('AS').mean()
pop = fred.get_series('CNP16OV').resample('AS').mean()
# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
# Population adjusted
y = y/pop
# Growth rate
y = np.log(y).diff()
y = y.to_frame()
y.rename(columns={0: "real_gdp_growth_rate"}, inplace=True)

df = pd.concat([df, y], axis=1)
df = df.loc["1978":"2022"]

#df.rename(columns={'firms': 'Firms', 'B': 'Column2'}, inplace=True)
# Episodes

df['firms_entry'] = df.firms.diff() + df.firmdeath_firms
df['firms_entry_rate'] = 100*df['firms_entry']/df['firms']
df['firms_exit_rate'] = 100*df['firmdeath_firms']/df['firms']

df.to_pickle("BDS_data_adj.pkl")
#df = pd.read_pickle("BDS_data_adj.pkl")





