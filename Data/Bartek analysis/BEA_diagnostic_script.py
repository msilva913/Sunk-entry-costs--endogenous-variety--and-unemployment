

import requests, json, os
key = "7878D119-FFC4-411E-969D-3B168F53A1E0"
#key = os.environ["BEA_API_KEY"]
r = requests.get("https://apps.bea.gov/api/data", params={
    "UserID": key, "method": "GetParameterValues",
    "DataSetName": "GDPbyIndustry", "ParameterName": "TableID",
    "ResultFormat": "JSON"
}, timeout=30)
data = r.json()
tables = data["BEAAPI"]["Results"]["ParamValue"]
for t in tables:
    print(t.get("Key"), t.get("Desc"))