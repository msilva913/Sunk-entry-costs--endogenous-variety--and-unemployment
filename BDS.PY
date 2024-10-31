import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

df = pd.read_csv('BDS.csv')

# Episodes
periods = {
    'Volcker Rate Rise': (1979, 1987),
    'Great Recession': (2007, 2009),
    'Pandemic Covid': (2020, 2023)
}

def format_millions(ax):
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e6)}M'))

# Plot job creation and destruction
plt.figure(figsize=(12, 6))
ax = plt.gca()
plt.plot(df['year'], df['job_creation'], label='Job Creation', color='blue')
plt.plot(df['year'], df['job_destruction'], label='Job Destruction', color='red')
for start, end in periods.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)
plt.title('Job Creation and Destruction with Key Episodes')
plt.xlabel('Year')
plt.ylabel('Jobs (in millions)')
format_millions(ax)
plt.legend(loc='upper left')
plt.grid(True)
plt.show()

# Plot net job creation rate
plt.figure(figsize=(12, 6))
plt.plot(df['year'], df['net_job_creation_rate'], label='Net Job Creation Rate', color='green')
for start, end in periods.values():
    plt.axvspan(start, end, color='gray', alpha=0.3)
plt.title('Net Job Creation Rate with Key Episodes')
plt.xlabel('Year')
plt.ylabel('Net Job Creation Rate')
plt.legend(loc='upper left')
plt.grid(True)
plt.show()

# Create tables
df[['job_creation', 'job_destruction']] /= 1e6
def create_summary_table(column_name):
    stats = pd.DataFrame()
    stats['Full'] = df[column_name].describe()
    for period, (start, end) in periods.items():
        period_data = df[(df['year'] >= start) & (df['year'] <= end)]
        stats[period] = period_data[column_name].describe()
    return stats.round(2)

print("Summary Statistics for Job Creation:")
display(create_summary_table('job_creation'))

print("\nSummary Statistics for Job Destruction:")
display(create_summary_table('job_destruction'))

print("\nSummary Statistics for Net Job Creation Rate:")
display(create_summary_table('net_job_creation_rate'))