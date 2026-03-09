import pandas as pd
import numpy as np
import os, Path   
 
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(Path.home() / "Documents" / "GitHub" / "Sunk_entry_costs_endogenous_variety_unemployment" / "Data" / "Bartek analysis")

from construct_delta_instrument import END_QUARTER as END_DELTA
from construct_s_instrument import END_QUARTER as END_S, Path
print(f"Delta instrument end: {END_DELTA}")
print(f"S instrument end:     {END_S}")

print(f"Delta instrument end: {END_DELTA}")
print(f"S instrument end:     {END_S}")

delta = pd.read_csv("data/instruments/delta_instrument_base2006.csv")
s     = pd.read_csv("data/instruments/s_instrument_base2006.csv")

both = delta.merge(s, on=["state", "state_fips", "quarter_label"])
pre  = both[both["quarter_label"] < "2020Q1"]

# Means
s_q     = pre["bartik_s"].mean()
delta_q = pre["bartik_delta"].mean()
tau_q   = s_q + delta_q * (1 - s_q)

s_m     = 1 - (1 - s_q)**(1/3)
delta_m = 1 - (1 - delta_q)**(1/3)
tau_m   = s_m + delta_m * (1 - s_m)

# Standard deviations
std_s_q     = pre["bartik_s"].std()
std_delta_q = pre["bartik_delta"].std()

# Delta method: d/dq [1-(1-q)^(1/3)] = (1/3)(1-q)^(-2/3)
std_s_m     = std_s_q     * (1/3) * (1 - s_q)**(-2/3)
std_delta_m = std_delta_q * (1/3) * (1 - delta_q)**(-2/3)

# Table 1: levels
table1 = pd.DataFrame(
    {
        "Quarterly (%)": [s_q * 100, delta_q * 100, tau_q * 100],
        "Monthly (%)":   [s_m * 100, delta_m * 100, tau_m * 100],
    },
    index=pd.Index(["s", "δ", "τ = s + δ(1−s)"], name="Parameter")
).round(3)

# Table 2: dispersion
table2 = pd.DataFrame(
    {
        "Quarterly (%)": [
            std_s_q * 100,
            std_delta_q * 100,
            (std_s_q / s_q) * 100,
            (std_delta_q / delta_q) * 100,
        ],
        "Monthly (%)": [
            std_s_m * 100,
            std_delta_m * 100,
            (std_s_m / s_m) * 100,
            (std_delta_m / delta_m) * 100,
        ],
    },
    index=pd.Index(
        ["Std(s)", "Std(δ)", "Std(s) / Mean(s)", "Std(δ) / Mean(δ)"],
        name="Parameter"
    )
).round(3)

print("--- Table 1: Steady-state rates ---")
print(table1.to_string())
print("\n--- Table 2: Cross-state and Time series dispersion (pre-COVID) ---")
print(table2.to_string())