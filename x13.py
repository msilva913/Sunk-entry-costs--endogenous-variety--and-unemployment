from statsmodels.tsa.x13 import x13_arima_analysis
from statsmodels.tsa.seasonal import seasonal_decompose
from fredapi import Fred
import matplotlib.pyplot as plt


fred = Fred(api_key='d35aabd7dc07cd94481af3d1e2f0ecf3')
df3 = fred.get_series("BFDUR4QTOTALNSAUS").resample("MS").mean().dropna()

# Path to X-13 (folder) executable
x13_path = r"E:\文档\Github\Sunk-entry-costs--endogenous-variety--and-unemployment\x13as"

result1 = x13_arima_analysis(df3, freq='M', x12path=x13_path, outlier=True, print_stdout=True)
result2 = seasonal_decompose(df3, model="additive", period=12)


x13_seasonally_adjusted = result1.seasadj  # Seasonally adjusted series from X-13
decompose_seasonally_adjusted = df3 - result2.seasonal  # Seasonally adjusted series from seasonal_decompose


plt.figure(figsize=(14, 8))
plt.plot(df3, label="Original Series", alpha=0.7)

plt.plot(x13_seasonally_adjusted, label="X-13 Seasonally Adjusted", color="red", alpha=0.7)
plt.plot(decompose_seasonally_adjusted, label="Seasonal Decompose Adjusted", color="green", alpha=0.7)

plt.legend(loc="best", fontsize=12)
plt.title("Comparison of Original, X-13, and Seasonal Decompose Adjusted Series", fontsize=16)
plt.xlabel("Date", fontsize=12)
plt.ylabel("DU4Q", fontsize=12)
plt.grid(True)


plt.tight_layout()
plt.show()