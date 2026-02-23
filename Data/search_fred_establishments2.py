from fredapi import Fred
import pandas as pd

FRED_API_KEY='9c70445138df124be4928605b7e08bd4'
fred = Fred(api_key=FRED_API_KEY)

res = fred.search('establishments', limit=200)
print('Total results:', len(res))
for idx, row in res.iterrows():
    title = row['title']
    sid = row['id']
    if 'Establishments' in title or 'establishments' in title or 'Establishment' in title:
        print(sid, '->', title)
