import pandas as pd
pd.set_option('display.max_columns', 10) 
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Data
#https://www.bls.gov/bdm/bdmover.htm#data BDM Overview
#download at https://data.bls.gov/PDQWeb/bd
"""
 Column labels
BDS0000000000000000110008RQ5 (1): Percentage employment lost from establishment deaths
BDS0000000000000000120008RQ5 (1): Percentage of establishments destroyed
BDS0000000000000000110007RQ5 (1): Percentage of employment gained from establishment births
BDS0000000000000000120007RQ5 (1): Percent of establishments with employment gained from births 
BDS0000000000000000110005RQ5 (1): Percent of employment lost from contractions
BDS0000000000000000120005RQ5 (1): Percent of establishments with employment lost from contractions for the total private sector in the U.S. (as a percent of total establishments in this sector)
"""

df = pd.read_excel('BED_data.xlsx', sheet_name='Data Import')

df['Series'] = pd.date_range(start="1992Q3", end="2024Q3", freq="QS")
df.set_index("Series", inplace=True)

df.rename(columns={"BDS0000000000000000120008RQ5": "estabs_exit_rate",
                   "BDS0000000000000000110008RQ5": "job_dest_rate"}, inplace=True)

recessions = [
    #(pd.Timestamp('1980-01-01'), pd.Timestamp('1980-07-31')),
    #(pd.Timestamp('1981-07-01'), pd.Timestamp('1982-11-30')),
    #(pd.Timestamp('1990-07-01'), pd.Timestamp('1991-03-31')),
    (pd.Timestamp('2001-03-01'), pd.Timestamp('2001-11-30')),
    (pd.Timestamp('2007-12-01'), pd.Timestamp('2009-06-30')),
    (pd.Timestamp('2020-02-01'), pd.Timestamp('2020-04-30')),
]


fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(df.estabs_exit_rate.loc["1992":"2019"], linewidth=2, alpha=0.7,
        label="Establishment exit:BED")
ax.plot(df.job_dest_rate.loc["1992":"2019"], linewidth=2, alpha=0.7,
        label="Job destruction:BED")
ax.set_xlabel("Time", fontsize=10)
ax.set_ylabel("Rate (%)", fontsize=10)
for start, end in recessions:
    ax.axvspan(start, end, color='gray', alpha=0.3)
ax.legend(loc="upper right", fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
#df.to_pickle("BDS_data_adj.pkl")
#df = pd.read_pickle("BDS_data_adj.pkl")

print("Annual BED dest rate =", df.estabs_exit_rate.loc["1992":"2019"].mean()*4)

"""
Purpose
Business Employment Dynamics data are quarterly series of gross job gains and 
gross job losses statistics for the entire economy. 
These data track changes in employment at the establishment level, 
and thus provide a picture of the dynamics underlying aggregate net employment growth statistics.

Top
Sources
The microdata used to construct the gross job gains and gross job losses
 statistics are from the
 Quarterly Census of Employment and Wages (QCEW) program. 
 These data include all establishments subject to 
 State unemployment insurance (UI) laws and 
 Federal agencies subject to the Unemployment Compensation for 
 Federal Employees program. 
 Each quarter, the State agencies edit and process the data and
 send the information to BLS in Washington, DC. 
 Major exclusions from UI coverage are the self-employed and
 certain nonprofit organizations. 
 Establishments report employment for the pay period including 
 the 12th of the month. The job flow estimates report employment changes between the third month of each quarter.
 """
 
"""
Concepts and Methodology
Business Employment Dynamics measure the net change in employment at the 
establishment level. These changes come about in one of four ways. 
A net increase in employment can come from either opening establishments or
 expanding establishments. A net decrease in employment can come from either closing establishments or contracting establishments. Gross job gains include the sum of all jobs added at either opening or expanding establishments. Gross job losses include the sum of all jobs lost in either closing or contracting establishments. The net change in employment is the difference between gross job gains and gross job losses.
The formal definitions of establishment-level employment changes are as follows:

Openings. These are either establishments with positive third month employment for the first time in the current quarter, with no links to the prior quarter, or with positive third month employment in the current quarter following zero employment in the previous quarter.

Expansions. These are establishments with positive employment in the third month in both the previous and current quarters, with a net increase in employment over this period.

Closings. These are either establishments with positive third month employment in the previous quarter, with no positive employment reported in the current quarter, or with positive third month employment in the previous quarter followed by zero employment in the current quarter.

Contractions. These are establishments with positive employment in the third month in both the previous and current quarters, with a net decrease in employment over this period.

All establishment-level employment changes are measured from the third month of each quarter. Not all establishments change their employment levels; these establishments are included in total employment, but do not affect counts of gross job gains and gross job losses.

Job flows are expressed as rates by dividing their levels by the average of employment in the current and previous quarters. This provides a symmetric growth rate. Job flows are calculated for the components of gross job gains and gross job losses and summed to form their respective totals. Job flow rates can be added and subtracted just as their levels can. For instance the difference between the gross job gains rate and the gross job loss rate is the net growth rate.
 """
 
 
"""
Comparison of BED to BDS
Scope and Coverage:

Business Employment Dynamics (BED): Focuses on employment changes in U.S. businesses. It provides data on job gains and losses from business openings, closures, expansions, and contractions. The BED is derived from the Quarterly Census of Employment and Wages (QCEW) and covers private-sector establishments.
Business Dynamics Statistics (BDS): Offers a broader view of business dynamics, including firm age and size, job creation and destruction, establishment births and deaths, and firm entry and exit. It uses data from the Census Bureau's Longitudinal Business Database, covering nearly all U.S. businesses.
Data Frequency:

BED: Provides quarterly data, allowing for more frequent analysis of employment trends over time.
BDS: Typically released annually, offering a comprehensive yearly snapshot of business dynamics.
Level of Detail:

BED: Focuses on employment dynamics at the establishment level, providing detailed insights into job flows within businesses.
BDS: Includes both establishment and firm-level data, offering insights into broader business dynamics, including firm age, size, and survival rates.
Methodology:

BED: Utilizes administrative data from unemployment insurance records, providing detailed employment change information.
BDS: Combines survey data with administrative records, allowing for more comprehensive analysis of business characteristics and dynamics.
"""