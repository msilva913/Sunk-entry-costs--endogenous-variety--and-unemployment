from fredapi import Fred
import pandas as pd

FRED_API_KEY='9c70445138df124be4928605b7e08bd4'
fred = Fred(api_key=FRED_API_KEY)

industries=['Nondurable Goods Manufacturing','Durable Goods Manufacturing','Construction','Trade, Transportation, and Utilities','Professional and Business Services','Leisure and Hospitality']

for ind in industries:
    print('\n=== Industry:', ind)
    try:
        res = fred.search('establishments '+ind, limit=50)
        # fredapi returns a pandas.DataFrame
        if isinstance(res, pd.DataFrame):
            for idx, row in res.iterrows():
                print(row['id'], '->', row['title'])
        else:
            print(res)
    except Exception as e:
        print('Search error:', e)
