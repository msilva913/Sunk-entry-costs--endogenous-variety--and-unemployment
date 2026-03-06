import pandas as pd
from pathlib import Path
cache = Path("data/cache")  # adjust if needed

bed = pd.read_parquet(cache / "bed_closings_national_12ind.parquet")
bds_path = list(cache.glob("bds_*.parquet"))
print("BDS cache files found:", bds_path)
print("\nBED closings sample (2010Q1):")
print(bed[bed["quarter_label"]=="2010Q1"][["industry_code","closings_nat"]])
