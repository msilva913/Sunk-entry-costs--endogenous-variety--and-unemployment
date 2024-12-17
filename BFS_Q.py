import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from fredapi import Fred
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from BLSAPI import bls

fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')

# Define episodes
episodes = {
    'Great Recession': ('2007-01-01', '2009-12-31'),
    'Pandemic Covid': ('2020-01-01', '2023-12-31')
}


df1 = fred.get_series('BFBF4QTOTALSAUS').resample('Q').mean().dropna()
df2 = fred.get_series("BFBF8QTOTALSAUS").resample("Q").mean().dropna()
df3 = bls.get_series("BDS0000000000000000120007LQ5").resample("Q").mean().dropna() * 1000 #(Because raw levels in thousands)
pop = fred.get_series('CNP16OV').resample('Q').mean().dropna()

# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=1600)[1]
df1 = df1/pop
df2 = df2/pop
df3 = df3/pop

df1.dropna(inplace=True)
df2.dropna(inplace=True)
df3.dropna(inplace=True)

plt.figure(figsize=(12, 6))
plt.plot(df1.index, np.log(df1.values), label='SBF4', color='blue')
plt.plot(df2.index, np.log(df2.values), label='SBF8', color='red')
plt.plot(df3.index, np.log(df3.values), label='BED', color='Purple')

# Plt.legend()

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel(' Log Per Capita Business Formations/ Number of Establishments Birth')
plt.legend(loc='upper left')
plt.grid(True)
plt.savefig('BFS_plot.pdf')
plt.show()

# For short
df3_shorter = df3.loc[df1.index]
    
plt.figure(figsize=(12, 6))
plt.plot(df1.index, np.log(df1.values), label='SBF4', color='blue')
plt.plot(df2.index, np.log(df2.values), label='SBF8', color='red')
plt.plot(df3_shorter.index, np.log(df3_shorter.values), label='BED (Shorter)', color='purple')

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel(' Log Per Capita Business Formations/ Number of Establishments birth')
plt.legend(loc='upper left')
plt.grid(True)
plt.savefig('BFS_s_plot.pdf')
plt.show()

# Summary
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
bed_stats = create_summary_stats(df3)

print("\nSBF8 Summary Statistics:")
display(sbf8_stats)
print("\nSBF8 Summary Statistics:")
display(sbf8_stats)
print("\nBED Summary Statistics:")
display(bed_stats)