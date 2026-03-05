
import requests

resp = requests.get(
    "https://download.bls.gov/pub/time.series/jt/jt.industry",
    timeout=30
)
for line in resp.text.splitlines():
    if "480" in line or "400" in line or "transport" in line.lower():
        print(line)


