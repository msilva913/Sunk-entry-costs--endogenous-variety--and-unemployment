import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from fredapi import Fred
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from BLSAPI import bls
from time_series_functions import (moments, stacked_moments, filter_transform)

fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')

# Define episodes
episodes = {
    'Great Recession': ('2007-12-01', '2009-07-01'),
    'Pandemic Covid': ('2020-03-01', '2022-02-01')
}

" Represent data in thousands "

# Total Business Formations within 4 Quarters
df1 = fred.get_series('BFBF4QTOTALSAUS').resample('QE').sum().dropna()/1000
# Business applications with planned wages
#df2 = fred.get_series("BFBF8QTOTALSAUS").resample("QE").sum().dropna()/1000
df2 = fred.get_series("BAWBATOTALSAUS").resample("QE").sum().dropna()/1000
# Establishment Births (Business Employment Dynamics)
df3 = bls.get_series("BDS0000000000000000120007LQ5").resample("QE").mean().dropna() #(Because raw levels in thousands)
pop = fred.get_series('CNP16OV').resample('QE').mean().dropna()

# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=1600)[1]
df1 = df1/pop
df2 = df2/pop
df3 = df3/pop

df = pd.concat([df1, df2, df3], axis=1)
df.columns = ["BF4", "BA", "EB"]

# Plot in levels (establishment per person)

fig, ax = plt.subplots(figsize=(12, 6))
ax.set_title("Business Formation")
ax.set_xlabel("Year")
ax.set_ylabel('Per Capita Business Formations/ Number of Establishments Birth')
ax.plot(df.index, df.BF4, label='BF4', color='blue')
ax.plot(df.index, df.BA, label='BA', color='red')
ax.plot(df.index, df.EB, label='EB', color='Purple')
ax.legend()
ax.grid(True)
plt.show()
# Plt.legend()


#plt.savefig('BFS_plot.pdf')

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

#init = "1993-06-30"
init = "2004-03-30"
final = "2019-12-31"
cycle_hp = df.apply(lambda x: filter_transform(x, init=init, final=final, transform_type='log',
                                      filter_type="hp_filter", lamb=100_000) )
cycle_ham = df.apply(lambda x: filter_transform(x, init=init, final=final, transform_type='log',
                                      filter_type="hamilton"))
mom_hp = moments(cycle_hp, relative_std="BF4", lab=["BF4", "EB"])
mom_ham = moments(cycle_ham, relative_std="BF4", lab=["BF4", "EB"])

df_log = np.log(df)
fig, ax = plt.subplots(figsize=(12, 6))
ax.set_title("Business application and formation")
ax.set_xlabel("Year")
ax.set_ylabel(' Cyclical component of log per Capita business formations/ number of establishments birth/application')
ax.plot(cycle_hp.index, cycle_hp.BF4, label='Business formation within 4 quarters', color='blue')
ax.plot(cycle_hp.index, cycle_hp.BA, label='Business applications with planned wages', color='red')
ax.plot(cycle_hp.index, cycle_hp.EB, label='Establishment births', color='Purple')
ax.legend()
ax.grid(True)
plt.show()

# For short
# df3_shorter = df3.loc[df1.index]
    
# plt.figure(figsize=(12, 6))
# plt.plot(df1.index, np.log(df1.values), label='SBF4', color='blue')
# plt.plot(df2.index, np.log(df2.values), label='SBF8', color='red')
# plt.plot(df3_shorter.index, np.log(df3_shorter.values), label='BED (Shorter)', color='purple')

# for start, end in episodes.values():
#     plt.axvspan(start, end, color='gray', alpha=0.3)

# plt.title('Business Formation')
# plt.xlabel('Year')
# plt.ylabel(' Log Per Capita Business Formations/ Number of Establishments birth')
# plt.legend(loc='upper left')
# plt.grid(True)
# plt.savefig('BFS_s_plot.pdf')
# plt.show()

# # Summary
# def create_summary_stats(df):
#     #full_stats = df.describe().rename(columns={'Value': 'Full'})
#     full_stats = df.describe()
#     episode_stats = []
#     for episode, (start, end) in episodes.items():
#         episode_data = df.loc[start:end]
#         #episode_stats.append(episode_data.describe().rename(columns={'episode_data': episode}))
#         episode_stats.append(episode_data.describe())
#     combined_stats = pd.concat([full_stats] + episode_stats, axis=1)
#     combined_stats.columns = ["Full", "Great Recession", "Pandemic"]
#     return combined_stats.round(2)

# sbf4_stats = create_summary_stats(df.BF4)
# sbf8_stats = create_summary_stats(df2)
# bed_stats = create_summary_stats(df3)

# print("\nSBF8 Summary Statistics:")
# display(sbf8_stats)
# print("\nSBF8 Summary Statistics:")
# display(sbf8_stats)
# print("\nBED Summary Statistics:")
# display(bed_stats)