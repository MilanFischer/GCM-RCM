import subprocess
import os
import re
import shutil
import xarray as xr
import cftime
import numpy as np
from datetime import datetime

# Set working directory and prepare output structure
working_dir = "/mnt/data/Other/GCM-RCM"
os.makedirs(working_dir, exist_ok=True)
os.chdir(working_dir)

# Output base structure
output_base = os.path.join(working_dir, "GCM", "CMIP5_OTHERS")

# Scenarios and models
scenarios = ["historical", "RCP85"]
models = [
    "CSIRO-Mk3-6-0",
    "GFDL-CM3",
    "GFDL-ESM2G",
    "GISS-E2-H",
    "GISS-E2-R",
    "IPSL-CM5B-LR",
    "MIROC5-ESM-CHEM",
    "MOHC-HadGEM2-AO",
    "MRI-ESM1"
]
# Regex to extract time structure
pattern = re.compile(r"_(r\d+i\d+p\d+)_(\d{8}-\d{8})\.nc$")

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

# Function to safely decode time using cftime
def safe_decode_time(ds):
    time_vals = ds.time.values.astype("float32")
    time_units = ds.time.attrs["units"]
    time_calendar = ds.time.attrs.get("calendar", "standard")

    decoded = cftime.num2date(
        time_vals,
        units=time_units,
        calendar=time_calendar,
        only_use_cftime_datetimes=True
    )
    ds["time"] = ("time", decoded)
    return ds

# Loop over all models and scenarios
for model in models:
    for scenario in scenarios:
        print(f"\n=== Processing {model} - {scenario} ===")

        # Reconstruct full output paths that match input structure
        daily_output_dir = os.path.join(
            output_base, "daily_nc", scenario, model
        )
        monthly_output_dir = os.path.join(
            output_base, "monthly_nc", scenario, model
        )

        # Ensure the full nested paths exist
        os.makedirs(daily_output_dir, exist_ok=True)
        os.makedirs(monthly_output_dir, exist_ok=True)

        base_path = f"/mnt/k_drive/Data/GCM/CMIP5_OTHERS/daily_nc/{scenario}/{model}/"
        if not os.path.isdir(base_path):
            print(f"  ➤ Directory not found: {base_path}, skipping.")
            continue

        all_files = os.listdir(base_path)

        # Index files by variant/grid/time
        file_index = {}
        for file in all_files:
            if any(var in file for var in ["tasmin", "tasmax", "huss", "psl"]):
                match = pattern.search(file)
                if match:
                    variant = match.group(1)
                    grid = ""  # CMIP5 doesn't use grid labels
                    timerange = match.group(2)

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
            if all(var in files for var in ["tasmin", "tasmax", "huss", "psl"]):
                print(f"  ➤ Processing time slice: {variant}, {grid}, {timerange}")

                # Copy input files to daily output directory
                local_files = {}
                for var in ["tasmin", "tasmax", "huss", "psl"]:
                    src = files[var]
                    dst = os.path.join(daily_output_dir, os.path.basename(src))
                    shutil.copy2(src, dst)
                    local_files[var] = dst

                # Load without decoding
                tasmin_ds = xr.open_dataset(local_files["tasmin"], decode_times=False)
                tasmax_ds = xr.open_dataset(local_files["tasmax"], decode_times=False)
                huss_ds = xr.open_dataset(local_files["huss"], decode_times=False)
                psl_ds = xr.open_dataset(local_files["psl"], decode_times=False)
                
                # Safely decode time
                tasmin_ds = safe_decode_time(tasmin_ds)
                tasmax_ds = safe_decode_time(tasmax_ds)
                huss_ds = safe_decode_time(huss_ds)
                psl_ds = safe_decode_time(psl_ds)

                huss = get_variable(huss_ds, ["huss"])
                psl = get_variable(psl_ds, ["psl"])
                tasmin = get_variable(tasmin_ds, ["tasmin"]) - 273.15
                tasmax = get_variable(tasmax_ds, ["tasmax"]) - 273.15

                # Mask invalid temperature values before computing anything
                tasmin = tasmin.where((tasmin > -100) & (tasmin < 100))
                tasmax = tasmax.where((tasmax > -100) & (tasmax < 100))

                # Compute VPD
                es_tasmin = 0.6108 * np.exp((17.27 * tasmin) / (tasmin + 237.3))
                es_tasmax = 0.6108 * np.exp((17.27 * tasmax) / (tasmax + 237.3))
                es = (es_tasmin + es_tasmax) / 2
                # Compute actual vapor pressure (kPa)
                ea = (psl * huss) / (0.378 * huss + 0.622) / 1000
                vpd = es - ea

                # Compute slope of the vapor pressure curve (kPa/K)
                tas = (tasmin + tasmax) / 2
                delta = (4098 * (0.6108 * np.exp((17.27 * tas) / (tas + 237.3)))) / ((tas + 237.3) ** 2)

                # Create Datasets
                vpd_ds = xr.Dataset(
                    {"vpd": vpd},
                    coords=huss.coords,
                    attrs={
                        "description": f"VPD derived from {model} {scenario}",
                        "units": "kPa"
                    }
                )

                delta_ds = xr.Dataset(
                    {"delta": delta},
                    coords=huss.coords,
                    attrs={
                        "description": f"Slope of vapor pressure curve derived from {model} {scenario}",
                        "units": "kPa/K"
                    }
                )

                # Output filenames
                # Clean suffix without extra underscores
                suffix_parts = [variant]
                if grid:  # only add if not empty
                    suffix_parts.append(grid)
                suffix_parts.append(timerange)
                common_suffix = "_".join(suffix_parts)
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
