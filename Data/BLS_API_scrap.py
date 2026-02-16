# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 17:14:36 2026

@author: msilv
"""

# --- Configuration ---
# 1. Enter your BLS API key here.
API_KEY = '3acb5544952e47c09e0b4625d640ce3f' 

def fetch_bls_data(series_dict, start_year, end_year, api_key):
    """
    Fetches multiple data series from the BLS API and returns a clean, pivoted DataFrame.
    This version is robust against series that contain data points without a 'value' key.
    """
    if not api_key or api_key == 'YOUR_API_KEY_HERE':
        print("ERROR: A valid BLS API key is required. Please register for one at https://data.bls.gov/registrationEngine/")
        return None

    headers = {'Content-type': 'application/json'}
    data = json.dumps({
        "seriesid": list(series_dict.values()),
        "startyear": start_year,
        "endyear": end_year,
        "registrationkey": api_key,
        "catalog": False
    })
    
    try:
        response = requests.post('https://api.bls.gov/publicAPI/v2/timeseries/data/', data=data, headers=headers, timeout=30)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"ERROR: Network request failed. {e}")
        return None

    json_data = response.json()
    
    if json_data.get('status') != 'REQUEST_SUCCEEDED':
        messages = json_data.get('message', [])
        print(f"ERROR: BLS API returned a failure status.")
        for msg in messages:
            print(f"  - {msg}")
        # Let's print the full error JSON to see what's happening
        print("Full error response:")
        print(json_data)
        return None

    id_to_name_map = {v: k for k, v in series_dict.items()}
    all_series_dfs = []

    for series_data in json_data['Results']['series']:
        bls_series_id = series_data['seriesID']
        custom_series_name = id_to_name_map.get(bls_series_id)
        
        if not custom_series_name:
            print(f"Warning: Received data for an unexpected series ID: {bls_series_id}")
            continue

        if not series_data['data']:
            print(f"Info: No data points returned for series '{custom_series_name}' ({bls_series_id}). Skipping.")
            continue

        # >>>>>>>>>>>>>>>>>>>>>>>> THE FINAL, ROBUST FIX <<<<<<<<<<<<<<<<<<<<<<<<<<<
        # 1. Filter the list to include only dictionaries that contain a 'value' key.
        valid_data_points = [point for point in series_data['data'] if 'value' in point]

        # 2. Check if there are any valid data points left after filtering.
        if not valid_data_points:
            print(f"Info: Series '{custom_series_name}' had data points, but none contained a 'value' key. Skipping.")
            continue
        # >>>>>>>>>>>>>>>>>>>>>>>>>>> END OF FIX <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        # 3. Create the DataFrame from the *filtered* list. This is now guaranteed to work.
        df = pd.DataFrame(valid_data_points)
        
        # This line will now succeed without a KeyError
        df[custom_series_name] = pd.to_numeric(df['value'], errors='coerce')
        
        df['date'] = pd.to_datetime(df['year'] + '-' + (df['period'].str.replace('Q0', '').astype(int) * 3).astype(str))
        
        df = df[['date', custom_series_name]].set_index('date')
        all_series_dfs.append(df)

    if not all_series_dfs:
        print("Error: After processing, no valid data could be compiled from any series.")
        return None
    
    final_df = pd.concat(all_series_dfs, axis=1)
    final_df = final_df.sort_index().ffill()
    
    return final_df

# 3. Call the function
bed_data_bls = fetch_bls_data(bls_series_ids, '2000', '2023', API_KEY)