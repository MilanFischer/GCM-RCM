import subprocess
import os
import re
import shutil
import xarray as xr
import numpy as np
from datetime import datetime

# Set working directory and prepare output structure
working_dir = "/mnt/data/Other/GCM-RCM"
os.makedirs(working_dir, exist_ok=True)
os.chdir(working_dir)

# Output base structure
output_base = os.path.join(working_dir, "GCM", "CMIP6")

# Scenarios and models
scenarios = ["historical", "ssp585"]
models = [
    "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1-HR",
    "EC-Earth3", "EC-Earth3-CC", "EC-EARTH3-Veg", "GFDL-CM4", "GFDL-ESM4", "HadGEM3-GC31-MM",
    "INM-CM5-0", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM", "TaiESM1"
]
# Regex to extract time structure
pattern = re.compile(r"_(r\d+i\d+p\d+f\d+)_([^_]+)_(\d{8}-\d{8})\.nc$")

# Function to auto-detect variable name
def get_variable(ds, possible_names):
    for name in possible_names:
        if name in ds:
            return ds[name]
    raise KeyError(f"None of the variable names {possible_names} found in dataset: {list(ds.data_vars)}")

# Function that checks whether a filename’s date range overlaps the target period.
def overlaps(target_start, target_end, file_range):
    file_start_str, file_end_str = file_range.split("-")
    file_start = datetime.strptime(file_start_str, "%Y%m%d")
    file_end = datetime.strptime(file_end_str, "%Y%m%d")
    return not (file_end < target_start or file_start > target_end)

# Loop over all models and scenarios
for model in models:
    for scenario in scenarios:
        print(f"\n=== Processing {model} - {scenario} ===")

        # Reconstruct full output paths that match input structure
        daily_output_dir = os.path.join(
            output_base, "daily_nc", "ScenarioMIP", scenario, model
        )
        monthly_output_dir = os.path.join(
            output_base, "monthly_nc", "ScenarioMIP", scenario, model
        )

        # Ensure the full nested paths exist
        os.makedirs(daily_output_dir, exist_ok=True)
        os.makedirs(monthly_output_dir, exist_ok=True)

        base_path = f"/mnt/k_drive/Data/GCM/CMIP6/daily_nc/ScenarioMIP/{scenario}/{model}/"
        if not os.path.isdir(base_path):
            print(f"  ➤ Directory not found: {base_path}, skipping.")
            continue

        all_files = os.listdir(base_path)

        # Index files by variant/grid/time
        file_index = {}
        for file in all_files:
            if any(var in file for var in ["tasmin", "tasmax", "hurs"]):
                match = pattern.search(file)
                if match:
                    variant = match.group(1)
                    grid = match.group(2)
                    timerange = match.group(3)

                    # Set scenario-specific time window
                    if scenario == "historical":
                        target_start = datetime(1981, 1, 1)
                        target_end = datetime(2005, 12, 31)
                    elif scenario.lower() in ["ssp585", "rcp85"]:
                        target_start = datetime(2076, 1, 1)
                        target_end = datetime(2100, 12, 31)
                    else:
                        continue  # Skip unrecognized scenarios
                    
                    # Check if file time range overlaps the target period
                    if not overlaps(target_start, target_end, timerange):
                        continue  # Skip this file
                    
                    key = (variant, grid, timerange)
                    var = file.split("_")[0]
                    file_index.setdefault(key, {})[var] = os.path.join(base_path, file)

        # Process complete triplets
        for (variant, grid, timerange), files in file_index.items():
            if all(var in files for var in ["tasmin", "tasmax", "hurs"]):
                print(f"  ➤ Processing time slice: {variant}, {grid}, {timerange}")

                # Copy input files to daily output directory
                local_files = {}
                for var in ["tasmin", "tasmax", "hurs"]:
                    src = files[var]
                    dst = os.path.join(daily_output_dir, os.path.basename(src))
                    shutil.copy2(src, dst)
                    local_files[var] = dst

                # Load data
                hurs_ds = xr.open_dataset(local_files["hurs"], chunks={"time": 365})
                tasmin_ds = xr.open_dataset(local_files["tasmin"], chunks={"time": 365})
                tasmax_ds = xr.open_dataset(local_files["tasmax"], chunks={"time": 365})

                hurs = get_variable(hurs_ds, ["hurs", "rh2m"])
                tasmin = get_variable(tasmin_ds, ["tasmin"]) - 273.15
                tasmax = get_variable(tasmax_ds, ["tasmax"]) - 273.15

                # Compute VPD
                es_tasmin = 0.6108 * np.exp((17.27 * tasmin) / (tasmin + 237.3))
                es_tasmax = 0.6108 * np.exp((17.27 * tasmax) / (tasmax + 237.3))
                es = (es_tasmin + es_tasmax) / 2
                ea = es * (hurs / 100)
                vpd = es - ea

                # Compute slope of the vapor pressure curve (kPa/K)
                tas = (tasmin + tasmax) / 2
                delta = (4098 * (0.6108 * np.exp((17.27 * tas) / (tas + 237.3)))) / ((tas + 237.3) ** 2)

                # Create Datasets
                vpd_ds = xr.Dataset(
                    {"vpd": vpd},
                    coords=hurs.coords,
                    attrs={
                        "description": f"VPD derived from {model} {scenario}",
                        "units": "kPa"
                    }
                )

                delta_ds = xr.Dataset(
                    {"delta": delta},
                    coords=hurs.coords,
                    attrs={
                        "description": f"Slope of vapor pressure curve derived from {model} {scenario}",
                        "units": "kPa/K"
                    }
                )

                # Output filenames
                common_suffix = f"{variant}_{grid}_{timerange}"
                vpd_daily_output_file =   os.path.join(daily_output_dir, f"vpd_day_{model}_{scenario}_{common_suffix}.nc")
                vpd_monthly_output_file = os.path.join(monthly_output_dir, f"vpd_Amon_{model}_{scenario}_{common_suffix}.nc")
                delta_daily_output_file = os.path.join(daily_output_dir, f"delta_day_{model}_{scenario}_{common_suffix}.nc")
                delta_monthly_output_file = os.path.join(monthly_output_dir, f"delta_Amon_{model}_{scenario}_{common_suffix}.nc")

                # Save daily file
                vpd_ds.to_netcdf(vpd_daily_output_file)
                print(f"    ✓ Daily VPD saved: {vpd_daily_output_file}")
                delta_ds.to_netcdf(delta_daily_output_file)
                print(f"    ✓ Daily slope saved: {delta_daily_output_file}")

                # Save monthly file using CDO
                cdo_command_vpd = [
                    "cdo", "-f", "nc4c", "-z", "zip", "monmean",
                    vpd_daily_output_file, vpd_monthly_output_file
                ]
                subprocess.run(cdo_command_vpd, check=True)
                print(f"    ✓ Monthly VPD saved: {vpd_monthly_output_file}")

                cdo_command_delta = [
                    "cdo", "-f", "nc4c", "-z", "zip", "monmean",
                    delta_daily_output_file, delta_monthly_output_file
                ]
                subprocess.run(cdo_command_delta, check=True)
                print(f"    ✓ Monthly slope saved: {delta_monthly_output_file}")

            else:
                print(f"    ✗ Incomplete time slice: {variant}, {grid}, {timerange}")
