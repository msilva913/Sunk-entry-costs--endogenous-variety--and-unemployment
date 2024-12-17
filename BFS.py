import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from fredapi import Fred
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.tsa.x13 import x13_arima_analysis
fred = Fred(api_key = 'd35aabd7dc07cd94481af3d1e2f0ecf3	')
#df1 = pd.read_excel('SBF4.xlsx')
#df2 = pd.read_excel('SBF8.xlsx')

#define episodes
episodes = {
    'Great Recession': ('2007-12-01', '2009-07-01'),
    'Pandemic Covid': ('2020-03-01', '2022-02-01')
}

df1 = fred.get_series('BFBF4QTOTALSAUS').resample('MS').mean().dropna()
df2 = fred.get_series("BFBF8QTOTALSAUS").resample("MS").mean().dropna()
df3 = fred.get_series("BFDUR4QTOTALNSAUS").resample("MS").mean().dropna()
bawba = fred.get_series("BAWBATOTALSAUS").resample("MS").mean().dropna()

# Seasonally adjusted series from X-13
x13_path = r"E:\文档\Github\Sunk-entry-costs--endogenous-variety--and-unemployment\x13as"
result = x13_arima_analysis(df3, freq='M', x12path=x13_path, outlier=True, print_stdout=True)
df3= result.seasadj  


pop = fred.get_series('CNP16OV').resample('MS').mean().dropna()
# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]

BF4 = df1/pop
BF8 = df2/pop
DU4 = df3
bawba = bawba/pop

df = pd.concat([BF4, BF8, bawba], axis=1)
df.columns = ["BF4", "BF8", "bawba"]
df.dropna(inplace=True)

df = np.log(df)
df = df - df.iloc[0, :]


#plot

plt.figure(figsize=(12, 6))
plt.plot(df.BF4, label='SBF4', color='blue')
plt.plot(df.BF8, label='SBF8', color='red')
plt.plot(df.bawba, label='Business applications: planned wages', color='magenta')

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel(' Log Per Capita Business application and formation')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('BFS_plot.pdf')
plt.show()

# DU4Q
plt.figure(figsize=(12, 6))
plt.plot(DU4, color='purple')
plt.xlabel('Year')
plt.ylabel('Average Duration in Quarters')
plt.grid(True)
plt.legend()
plt.title("DU4Q")
for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)
    
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
dur4_stats = create_summary_stats(df3)

print("SBF4 Summary Statistics:")
display(sbf4_stats)
print("\nSBF8 Summary Statistics:")
display(sbf8_stats)
print("\nDUR4 Summary Statistics:")
display(dur4_stats)