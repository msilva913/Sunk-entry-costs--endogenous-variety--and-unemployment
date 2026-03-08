import pandas as pd
from pathlib import Path

from construct_delta_instrument import _quarter_to_bds_year
cache = Path("data/cache")  # adjust if needed

bed = pd.read_parquet(cache / "bed_closings_national_12ind.parquet")
bds_path = list(cache.glob("bds_*.parquet"))
print("BDS cache files found:", bds_path)
print("\nBED closings sample (2010Q1):")
print(bed[bed["quarter_label"]=="2010Q1"][["industry_code","closings_nat"]])

import pandas as pd
bed = pd.read_parquet("data/cache/bed_closings_national_12ind.parquet")
bds = pd.read_parquet("data/cache/bds_exits_national_supersector.parquet")

# Show a representative year — 2010
bed_2010 = bed[bed["quarter_label"].str.startswith("2010")].groupby("industry_code")["closings_nat"].sum()
bds_2010 = bds[bds["bds_year"] == 2011]["bds_exits"]  # 2011 covers 2010Q2-2011Q1
print("BED annual sum 2010:\n", bed_2010.sort_index())
print("\nBDS exits 2011:\n", bds_2010.values)

bds = pd.read_parquet("data/cache/bds_job_destruction_deaths_supersector.parquet")
bed = pd.read_parquet("data/cache/bed_closings_national_12ind.parquet")
bed["bds_year"] = bed["quarter_label"].apply(_quarter_to_bds_year)
bed_sum = bed.groupby(["industry_code","bds_year"])["closings_nat"].sum()
print("Mining BDS exits (thousands):")
print(bds[bds["industry_code"]=="10"].set_index("bds_year")["bds_exits"].head(10))
print("\nMining BED annual sum:")
print(bed_sum["10"].head(10))


# Point at wherever your cache landed
cache = Path("data/cache")   # adjust if needed

bds = pd.read_parquet(cache / "bds_job_destruction_deaths_supersector.parquet")
print("BDS Mining (industry_code=10), first 5 years:")
print(bds[bds["industry_code"]=="10"].head(5)[["bds_year","bds_exits"]])

perm = pd.read_parquet("data/instruments/permanence_ratios_by_supersector.parquet")
print("\nPermanence Mining rows:")
print(perm[perm["industry_code"]=="10"].head(5)[["bds_year","pi","bds_exits","bed_closings_sum"]])
#########################################

perm = pd.read_parquet("data/instruments/permanence_ratios_by_supersector.parquet")

# BDS year 2021 covers 2020Q2–2021Q1 — this is the COVID window
print("π values for BDS year 2021:")
print(perm[perm["bds_year"]==2021][["industry_code","pi","bds_exits","bed_closings_sum"]]
      .sort_values("pi").to_string())


s = pd.read_parquet("data/instruments/shock_rates_1992Q3_2023Q1.parquet")
# 1. Ratio of 2020Q2 to 2019Q2 — how big is the spike after calibration?
m19 = s[s["quarter_label"]=="2019Q2"].groupby("industry_code")["g_delta_loo"].mean()
m20 = s[s["quarter_label"]=="2020Q2"].groupby("industry_code")["g_delta_loo"].mean()
print("2020Q2 / 2019Q2 ratio (calibrated):")
print((m20/m19).sort_values(ascending=False).round(2))

# 2. Is there a perm_adjusted column or is it baked into g_delta_loo?
print("\nColumns:", s.columns.tolist())
print("\n2020Q2 Leisure & hospitality (code 70):")
print(s[(s["quarter_label"]=="2020Q2") & (s["industry_code"]=="70")]["g_delta_loo"].describe())
    
