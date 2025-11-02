import subprocess
import os
import re
import xarray as xr
import cftime
import numpy as np
from datetime import datetime

# Set working directory and prepare output structure
base_dir = "/mnt/data/Other/GCM-RCM/GCM/CMIP5"

# Scenarios and models
scenarios = ["historical", "RCP85"]
models = [
    "bcc-csm1-1",
    "CESM1-BGC",
    "CESM1-CAM5",
    "CMCC-CESM",
    "CMCC-CMS",
    "CNRM-CM5",
    "EC-EARTH",
    "GFDL-CM3",
    "GFDL-ESM2G",
    "HadGEM2-ES",
    "inmcm4",
    "IPSL-CM5A-MR",
    "MPI-ESM-LR",
    "MRI-ESM1",
    "NorESM1-M"
]
# Regex to extract time structure
pattern = re.compile(r"_(r\d+i\d+p\d+)_(\d{8}-\d{8})\.nc$")

# Constants for hypsometric equation
ELR = 0.0065                # K/m – environmenatal lapse rate (temperature drop per meter of altitude)
R = 8.31451                 # J/mol/K – universal gas constant (used in ideal gas law)
M_dry = 0.0289644           # kg/mol – molar mass of dry air
M_vap = 0.01801528          # kg/mol – molar mass of water vapor
g = 9.80665                 # m/s^2 – acceleration due to gravity (standard value at sea level)
R_dry = R / M_dry           # J/kg/K, specific gas constant for dry air
R_vap = R / M_vap           # J/kg/K, specific gas constant for dry air

def saturation_vapor_pressure(tas_C):
    """
    Compute saturation vapor pressure in kPa
    - over water for tas >= 0°C
    - over ice for tas < 0°C
    """
    # ECMWF (2018) IFS DOCUMENTATION – Cy45r1, Part IV: Physical Processes, Chapter 12
    es = np.where(
        tas_C >= 0,
        0.61121 * np.exp((17.502 * tas_C) / (tas_C + 240.96)),  # over water
        0.61121 * np.exp((22.587 * tas_C) / (tas_C + 273.85))   # over ice
    )
    return es  # in kPa

def slope_of_saturation_vapor_pressure_curve(tas_C):
    """
    Compute slope of saturation vapor pressure curve in kPa per degree Celsius
    - over water for tas >= 0°C
    - over ice for tas < 0°C
    """
    e_sat = saturation_vapor_pressure(tas_C)

    delta = np.where(
        tas_C >= 0,
        17.502 * 240.96 * e_sat/(240.96 + tas_C) ** 2,     # over water
        22.587 * 273.85 * e_sat/(273.85 + tas_C) ** 2      # over ice
    )
    return delta  # in kPa per degree Celsius

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
        daily_dir = os.path.join(
            base_dir, "daily_nc", scenario, model
        )
        monthly_dir = os.path.join(
            base_dir, "monthly_nc", scenario, model
        )

        # Ensure the full nested paths exist
        os.makedirs(monthly_dir, exist_ok=True)

        all_files = os.listdir(daily_dir)

        # Index files by variant/grid/time
        file_index = {}
        for file in all_files:
            if any(var in file for var in ["tasmin", "tasmax", "tas", "huss", "psl"]):
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
                    file_index.setdefault(key, {})[var] = os.path.join(daily_dir, file)

        # Process complete triplets
        for (variant, grid, timerange), files in file_index.items():
            if all(var in files for var in ["tasmin", "tasmax", "tas", "huss", "psl"]):
                print(f"  ➤ Processing time slice: {variant}, {grid}, {timerange}")

                # Load without decoding
                tasmin_ds = xr.open_dataset(files["tasmin"], decode_times=False)
                tasmax_ds = xr.open_dataset(files["tasmax"], decode_times=False)
                tas_ds = xr.open_dataset(files["tas"], decode_times=False)
                huss_ds = xr.open_dataset(files["huss"], decode_times=False)
                psl_ds = xr.open_dataset(files["psl"], decode_times=False)
                
                # Safely decode time
                tasmin_ds = safe_decode_time(tasmin_ds)
                tasmax_ds = safe_decode_time(tasmax_ds)
                tas_ds = safe_decode_time(tas_ds)
                huss_ds = safe_decode_time(huss_ds)
                psl_ds = safe_decode_time(psl_ds)
                
                huss = get_variable(huss_ds, ["huss"])
                psl = get_variable(psl_ds, ["psl"])
                tasmin = get_variable(tasmin_ds, ["tasmin"]) - 273.15
                tasmax = get_variable(tasmax_ds, ["tasmax"]) - 273.15
                tas = get_variable(tas_ds, ["tas"]) - 273.15

                # Mask invalid temperature values before computing anything
                tasmin = tasmin.where((tasmin > -100) & (tasmin < 100))
                tasmax = tasmax.where((tasmax > -100) & (tasmax < 100))
                tas = tas.where((tas > -100) & (tas < 100))

                # Load orography (assumes variable name is 'orog')
                orog_dir = os.path.join(
                    base_dir, "orography", model
                )
                ds_orog = xr.open_dataset(f'{orog_dir}/orog_fx_{model}_historical_r0i0p0.nc')
                orog = ds_orog['orog']  # in meters

                # Align dimensions
                psl, tas = xr.align(psl, tas)
                orog = orog.interp_like(tas, method='nearest')

                # Compute surface pressure (ps) from psl and orography
                ps = psl * (1 - (ELR * orog / (tas + 273.15))) ** (g / (R_dry * ELR))

                # Compute saturation vapor pressure (kPa) from tas
                es_tas = saturation_vapor_pressure(tas)

                # Compute actual vapor pressure (kPa)
                ea = (ps * huss) / ((1 -R_dry/R_vap) * huss + R_dry/R_vap) / 1000

                # Compute relative humidity (RH)
                RH = ea / es_tas * 100  # Relative Humidity
                RH =  RH.clip(min=0, max=100) # Clip RH to [0, 100]

                # Now simplify to use it consitently with CMIP6 where hurs is used
                # More precise would be to use the function taking into account ice and water
                # saturation_vapor_pressure(tasmin) and saturation_vapor_pressure(tasmax)
                # but this is skipped for simplicity now

                # Compute VPD
                es_tasmin = 0.6108 * np.exp((17.27 * tasmin) / (tasmin + 237.3))
                es_tasmax = 0.6108 * np.exp((17.27 * tasmax) / (tasmax + 237.3))
                es = (es_tasmin + es_tasmax) / 2

                ea = es * (RH / 100)
                vpd = es - ea               
                vpd = vpd.clip(min=0)  # Replace negative values with 0

                # Compute slope of the vapor pressure curve (kPa/K)
                # delta = slope_of_saturation_vapor_pressure_curve((tasmin + tasmax) / 2)
                tas = (tasmin + tasmax) / 2
                delta = (4098 * (0.6108 * np.exp((17.27 * tas) / (tas + 237.3)))) / ((tas + 237.3) ** 2)

                # Create Datasets
                RH_ds = xr.Dataset(
                    {"RH": RH},
                    coords=huss.coords,
                    attrs={
                        "description": f"RH derived from {model} {scenario}",
                        "units": "%"
                    }
                )

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
                RH_daily_output_file =   os.path.join(daily_dir, f"RH_day_{model}_{scenario}_{common_suffix}.nc")
                RH_monthly_output_file = os.path.join(monthly_dir, f"RH_Amon_{model}_{scenario}_{common_suffix}.nc")
                vpd_daily_output_file =   os.path.join(daily_dir, f"vpd_day_{model}_{scenario}_{common_suffix}.nc")
                vpd_monthly_output_file = os.path.join(monthly_dir, f"vpd_Amon_{model}_{scenario}_{common_suffix}.nc")
                delta_daily_output_file = os.path.join(daily_dir, f"delta_day_{model}_{scenario}_{common_suffix}.nc")
                delta_monthly_output_file = os.path.join(monthly_dir, f"delta_Amon_{model}_{scenario}_{common_suffix}.nc")

                # Save daily file
                RH_ds.to_netcdf(RH_daily_output_file)
                print(f"    ✓ Daily RH saved: {RH_daily_output_file}")
                vpd_ds.to_netcdf(vpd_daily_output_file)
                print(f"    ✓ Daily VPD saved: {vpd_daily_output_file}")
                delta_ds.to_netcdf(delta_daily_output_file)
                print(f"    ✓ Daily slope saved: {delta_daily_output_file}")

                # Save monthly file using CDO
                cdo_command_RH = [
                    "cdo", "-f", "nc4c", "-z", "zip", "monmean",
                    RH_daily_output_file, RH_monthly_output_file
                ]
                subprocess.run(cdo_command_RH, check=True)
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
