import requests
import pandas as pd
from datetime import datetime


class bls:
    api_key = "b1e1257b28b441b28a1ac602e7df7006"
    base_url = "https://api.bls.gov/publicAPI/v2/timeseries/data/"

    @classmethod
    def get_series(cls, series_id, start_year=None, end_year=None, freq="Q"):

        all_data_frames = []
        chunk_size = 20 

        if start_year is None:
            start_year = 1900
        if end_year is None:
            end_year = datetime.now().year

        for chunk_start in range(start_year, end_year + 1, chunk_size):
            chunk_end = min(chunk_start + chunk_size - 1, end_year)
            payload = {
                "seriesid": [series_id],
                "startyear": str(chunk_start),
                "endyear": str(chunk_end),
                "registrationkey": cls.api_key,
            }
            response = requests.post(cls.base_url, json=payload)
            if response.status_code != 200:
                raise Exception(f"API request failed with status code {response.status_code}: {response.text}")
            data = response.json()
            if data["status"] != "REQUEST_SUCCEEDED":
                raise Exception(f"BLS API request failed: {data.get('message', 'Unknown error')}")
            series_data = data.get("Results", {}).get("series", [])[0].get("data", [])
            if not series_data:
                continue
            df = pd.DataFrame(series_data)
            df["value"] = pd.to_numeric(df["value"], errors="coerce")
            df["year"] = pd.to_numeric(df["year"], errors="coerce")

            if df["period"].str.startswith("Q").any():  # Quarterly data
                df["quarter"] = df["period"].str[1:].astype(int)
                df["date"] = pd.PeriodIndex(df["year"].astype(str) + "Q" + df["quarter"].astype(str), freq="Q").to_timestamp()
            elif df["period"].str.startswith("M").any():  # Monthly data
                df["month"] = df["period"].str[1:].astype(int)
                df["date"] = pd.to_datetime(df[["year", "month"]].assign(day=1))
            elif df["period"].str.startswith("A").any():  # Annual data
                df["date"] = pd.to_datetime(df["year"].astype(str) + "-01-01")
            else:
                raise ValueError("Unsupported period format.")
            df = df.set_index("date").sort_index()
            all_data_frames.append(df[["value"]])

        if all_data_frames:
            combined_data = pd.concat(all_data_frames)

            if freq == "Q":  # Quarterly
                combined_data = combined_data.resample("Q").mean()
            elif freq == "M":  # Monthly
                combined_data = combined_data.resample("M").mean()
            elif freq == "A":  # Annual
                combined_data = combined_data.resample("A").mean()
            else:
                raise ValueError("Invalid frequency. Choose 'M', 'Q', or 'A'.")

            return combined_data["value"]
        else:
            raise ValueError("No data could be retrieved for the given series ID and date range.")