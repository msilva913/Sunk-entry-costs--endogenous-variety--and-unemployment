import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from fredapi import Fred
import numpy as np
import statsmodels.api as sm
fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')
#df1 = pd.read_excel('SBF4.xlsx')
#df2 = pd.read_excel('SBF8.xlsx')

#define episodes
episodes = {
    'Great Recession': ('2007-12-01', '2009-07-01'),
    'Pandemic Covid': ('2020-03-01', '2022-02-01')
}

#df_melted1 = df1.melt(id_vars='Year', var_name='Month', value_name='Value')
#df_melted2 = df2.melt(id_vars='Year', var_name='Month', value_name='Value')

df1 = fred.get_series('BFBF4QTOTALSAUS').resample('MS').mean().dropna()
df2 = fred.get_series("BFBF8QTOTALSAUS").resample("MS").mean().dropna()
bawba = fred.get_series("BAWBATOTALSAUS").resample("MS").mean().dropna()
pop = fred.get_series('CNP16OV').resample('MS').mean().dropna()
# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
BF4 = df1/pop
BF8 = df2/pop
bawba = bawba/pop
df = pd.concat([BF4, BF8, bawba], axis=1)
df.columns = ["BF4", "BF8", "bawba"]
df.dropna(inplace=True)
# Divide business formation by population

#plot
df = np.log(df)
df = df - df.iloc[0, :]

plt.figure(figsize=(12, 6))
plt.plot(df.BF4, label='SBF4', color='blue')
plt.plot(df.BF8, label='SBF8', color='red')
plt.plot(df.bawba, label='Business applications: planned wages', color='magenta')
#plt.legend()

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel(' Log Per Capita Business application and formation')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('BFS_plot.pdf')
plt.show()

#summary
def create_summary_stats(df):
    #full_stats = df.describe().rename(columns={'Value': 'Full'})
    full_stats = df.describe()
    episode_stats = []
    for episode, (start, end) in episodes.items():
        episode_data = df.loc[start:end]
        #episode_stats.append(episode_data.describe().rename(columns={'episode_data': episode}))
        episode_stats.append(episode_data.describe())
    combined_stats = pd.concat([full_stats] + episode_stats, axis=1)
    combined_stats.columns = ["Full", "Great Recession", "Pandemic"]
    return combined_stats.round(2)

sbf4_stats = create_summary_stats(df1)
sbf8_stats = create_summary_stats(df2)

print("SBF4 Summary Statistics:")
display(sbf4_stats)
print("\nSBF8 Summary Statistics:")
display(sbf8_stats)