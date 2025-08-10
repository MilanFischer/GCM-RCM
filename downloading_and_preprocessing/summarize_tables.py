import os
import pandas as pd
import re
from collections import defaultdict

def parse_model_files(files, variable):
    """
    Parse files and group them by model and scenario (e.g., 'historical', 'ssp585', 'RCP85').
    """
    pattern = re.compile(
        rf"^{variable}_Amon_(.+?)_((?:historical)|(?:ssp585)|(?:RCP85)|(?:rcp85))_r\d+i\d+p\d+(?:f\d+)?_.*\.txt$"
    )

    file_info = defaultdict(dict)
    for file in files:
        match = pattern.match(file)
        if match:
            model, scenario = match.groups()
            file_info[model][scenario.lower()] = file  # normalize scenario names to lowercase
    return file_info

def extract_monthly_means(filepath):
    # Read header line
    with open(filepath) as f:
        headers = f.readline().strip().split()
    
    # Fix duplicate column names
    seen = {}
    clean_headers = []
    for h in headers:
        if h not in seen:
            seen[h] = 1
            clean_headers.append(h)
        else:
            seen[h] += 1
            clean_headers.append(f"{h}_{seen[h]}")

    # Read the actual data
    df = pd.read_csv(filepath, skiprows=1, sep=r"\s+", header=None, names=clean_headers)

    # Identify the correct columns
    date_col = next((col for col in df.columns if col.startswith("Date")), None)
    mean_col = next((col for col in df.columns if col.startswith("Mean")), None)

    if date_col is None or mean_col is None:
        raise ValueError(f"Missing 'Date' or 'Mean' column in {filepath}")

    df["Date"] = pd.to_datetime(df[date_col], errors="coerce")
    df["Month"] = df["Date"].dt.month
    df[mean_col] = pd.to_numeric(df[mean_col], errors="coerce")  # Force non-numeric to NaN

    monthly_means = df.groupby("Month")[mean_col].mean()
    return monthly_means

def process_climate_folder(path_in, variable, ensemble):
    files = [f for f in os.listdir(path_in) if f.endswith(".txt")]
    file_map = parse_model_files(files, variable)

    data_hist = pd.DataFrame(index=range(1, 13))
    data_fut = pd.DataFrame(index=range(1, 13))

    for model, scenarios in file_map.items():
        if 'historical' in scenarios and ('ssp585' in scenarios or 'rcp85' in scenarios):
            fut_scenario = 'ssp585' if 'ssp585' in scenarios else 'rcp85'
            print(f"Processing {model}: historical + {fut_scenario}")

            hist_means = extract_monthly_means(os.path.join(path_in, scenarios['historical']))
            fut_means = extract_monthly_means(os.path.join(path_in, scenarios[fut_scenario]))

            data_hist[model] = hist_means
            data_fut[model] = fut_means
        else:
            print(f"Skipping {model}: missing historical or future scenario.")

    # Save output
    out_base = os.path.basename(os.path.normpath(path_in))
    hist_file = os.path.join(path_in, f"{out_base}_1981-2005.csv")
    fut_file = os.path.join(path_in, f"{out_base}_2076-2100.csv")

    data_hist.to_csv(hist_file, index=False)
    data_fut.to_csv(fut_file, index=False)
    print(f"Saved: {hist_file}, {fut_file}")

# ▶ Example usage for CMIP6 or CMIP5
# EURO-CORDEX_EUR44X
process_climate_folder(
    path_in="/mnt/data/Other/GCM-RCM/climatology_tabs/EURO-CORDEX/EUR44",
    variable="vpd"
)

# CMIP5_4_EURO-CORDEX
process_climate_folder(
    path_in="/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP5_4_EURO-CORDEX",
    variable="vpd"
)

# CMIP5_4_EURO-CORDEX
process_climate_folder(
    path_in="/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP5_OTHERS",
    variable="vpd"
)

# CMIP6
process_climate_folder(
    path_in="/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP6",
    variable="vpd"
)
