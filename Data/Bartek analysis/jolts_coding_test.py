import requests, zipfile, io, pandas as pd

url = "https://data.bls.gov/cew/data/files/2006/csv/2006_annual_singlefile.zip"
resp = requests.get(url, timeout=300)
with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
    with zf.open([n for n in zf.namelist() if n.endswith(".csv")][0]) as f:
        raw = pd.read_csv(f, dtype=str, low_memory=False)

raw.columns = raw.columns.str.strip().str.lower().str.replace(" ", "_")
for col in ["area_fips","own_code","industry_code","agglvl_code"]:
    if col in raw.columns:
        raw[col] = raw[col].str.strip()

# Show every distinct (agglvl_code, industry_code) that appears for
# private-sector, national or state rows, where the industry code
# looks like it could be retail, transport, or utilities
mask = (
    (raw["own_code"] == "5") &
    (raw["area_fips"].str.endswith("000")) &          # state or national rows
    (raw["industry_code"].str.match(r'^[24][0-9]'))   # starts with 2x or 4x
)
print(raw.loc[mask, ["agglvl_code","industry_code"]]
      .drop_duplicates()
      .sort_values(["agglvl_code","industry_code"])
      .to_string(index=False))

