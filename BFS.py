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
    'Great Recession': ('2007-01-01', '2009-12-31'),
    'Pandemic Covid': ('2020-01-01', '2023-12-31')
}

#df_melted1 = df1.melt(id_vars='Year', var_name='Month', value_name='Value')
#df_melted2 = df2.melt(id_vars='Year', var_name='Month', value_name='Value')

df1 = fred.get_series('BFBF4QTOTALSAUS').resample('MS').mean().dropna()
df2 = fred.get_series("BFBF8QTOTALSAUS").resample("MS").mean().dropna()
df3 = fred.get_series("BFDUR4QTOTALNSAUS").resample("MS").mean().dropna()

x13_path = r"E:\文档\Github\Sunk-entry-costs--endogenous-variety--and-unemployment\x13as"
result = x13_arima_analysis(df3, freq='M', x12path=x13_path, outlier=True, print_stdout=True)
df3_adj= result.seasadj  # Seasonally adjusted series from X-13

plt.figure(figsize=(12, 6))
plt.plot(df3, label="Original Series", alpha=0.7)
plt.plot(df3_adj, label="X-13 Adjusted Series", color="red", alpha=0.7)
plt.legend()
plt.title("DU4Q")
plt.show()

pop = fred.get_series('CNP16OV').resample('MS').mean().dropna()
# Use HP-filtered trend for population to avoid discrete jumps around census dates
pop = sm.tsa.filters.hpfilter(pop, lamb=10_000)[1]
df1 = df1/pop
df2 = df2/pop
df3 = df3_adj

df1.dropna(inplace=True)
df2.dropna(inplace=True)
df3.dropna(inplace=True)
# Divide business formation by population

#mapping
#month_map = {
#    'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
#    'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
#}

#define df
# df_melted1['Month'] = df_melted1['Month'].map(month_map)
# df_melted2['Month'] = df_melted2['Month'].map(month_map)

# df_melted1['Year-Month'] = pd.to_datetime(df_melted1['Year'].astype(str) + '-' + df_melted1['Month'].astype(str))
# df_melted2['Year-Month'] = pd.to_datetime(df_melted2['Year'].astype(str) + '-' + df_melted2['Month'].astype(str))

# df_melted1.sort_values('Year-Month', inplace=True)
# df_melted2.sort_values('Year-Month', inplace=True)

# df_melted1.set_index('Year-Month', inplace=True)
# df_melted2.set_index('Year-Month', inplace=True)

# df_melted1.drop(['Year', 'Month'], axis=1, inplace=True)
# df_melted2.drop(['Year', 'Month'], axis=1, inplace=True)

#plot
plt.figure(figsize=(12, 6))
plt.plot(df1.index, np.log(df1.values), label='SBF4', color='blue')
plt.plot(df2.index, np.log(df2.values), label='SBF8', color='red')
#plt.legend()

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel(' Log Per Capita Business Formations')
plt.legend(loc='upper left')
plt.grid(True)
plt.savefig('BFS_plot.pdf')
plt.show()

plt.plot(df3.index, df3.values, label='DUR4Q', color='green')

for start, end in episodes.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)

plt.title('Business Formation')
plt.xlabel('Year')
plt.ylabel('Number of Quarters')
plt.ylim(0.4, 2.2)
plt.legend(loc='upper left')
plt.grid(True)
plt.savefig('DUR4_plot.pdf')
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