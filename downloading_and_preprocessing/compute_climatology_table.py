import os
import subprocess
import glob
import re

# ---------------------------------------------------------
# This script computes the climatology of a specified variable
# from multiple GCMs and RCMs, and saves the outputs to a .txt file.
# It uses CDO (Climate Data Operators) for processing NetCDF files.

# ----------------------------------------------------------
# Configuration

# "EURO-CORDEX_EUR44", "CMIP5_4_EURO-CORDEX", "CMIP5_OTHERS", "CMIP5", or "CMIP6"
ensemble = "EURO-CORDEX_EUR44"

# Variable to process - "vpd", "delta", "RH"
variable = "delta"

# Region bounding box
lonmin, lonmax = 9, 21
latmin, latmax = 47, 54

if ensemble == "EURO-CORDEX_EUR44":
    root = "/mnt/data/Other/GCM-RCM/RCM/EURO-CORDEX/EUR44/EUR44_monthly_nc"
    output_txt_dir = f"/mnt/data/Other/GCM-RCM/climatology_tabs/EURO-CORDEX/EUR44/{variable}"
    models = [
        "CCCma-CanESM2_r1i1p1_SMHI-RCA4_v1", "CNRM-CERFACS-CNRM-CM5_r1i1p1_CNRM-ALADIN53_v1", "CNRM-CERFACS-CNRM-CM5_r1i1p1_HMS-ALADIN52_v1",
        "CNRM-CERFACS-CNRM-CM5_r1i1p1_SMHI-RCA4_v1", "CSIRO-QCCCE-CSIRO-Mk3-6-0_r1i1p1_SMHI-RCA4_v1", "ICHEC-EC-EARTH_r12i1p1_SMHI-RCA4_v1",
        "ICHEC-EC-EARTH_r1i1p1_KNMI-RACMO22E_v1", "ICHEC-EC-EARTH_r3i1p1_DMI-HIRHAM5_v1", "IPSL-IPSL-CM5A-MR_r1i1p1_IPSL-INERIS-WRF331F_v1",
        "IPSL-IPSL-CM5A-MR_r1i1p1_SMHI-RCA4_v1", "MIROC-MIROC5_r1i1p1_SMHI-RCA4_v1", "MOHC-HadGEM2-ES_r1i1p1_KNMI-RACMO22E_v2",
        "MOHC-HadGEM2-ES_r1i1p1_SMHI-RCA4_v1", "MPI-M-MPI-ESM-LR_r1i1p1_CLMcom-CCLM4-8-17_v1", "MPI-M-MPI-ESM-LR_r1i1p1_MPI-CSC-REMO2009_v1",
        "MPI-M-MPI-ESM-LR_r1i1p1_SMHI-RCA4_v1", "MPI-M-MPI-ESM-LR_r2i1p1_MPI-CSC-REMO2009_v1", "NCC-NorESM1-M_r1i1p1_SMHI-RCA4_v1",
        "NOAA-GFDL-GFDL-ESM2M_r1i1p1_SMHI-RCA4_v1"
    ]

elif ensemble == "CMIP5_4_EURO-CORDEX":
    root = "/mnt/data/Other/GCM-RCM/GCM/CMIP5_4_EURO-CORDEX/monthly_nc"
    output_txt_dir = f"/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP5_4_EURO-CORDEX/{variable}"
    models = [
        "CanESM2", "CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "IPSL-CM5A-MR", "MPI-ESM-LR", "NorESM1-M"
    ]
elif ensemble == "CMIP5_OTHERS":
    root = "/mnt/data/Other/GCM-RCM/GCM/CMIP5_OTHERS/monthly_nc"
    output_txt_dir = f"/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP5_OTHERS/{variable}"
    models = [
        "CSIRO-Mk3-6-0", "GFDL-CM3", "GFDL-ESM2G", "GISS-E2-H", "GISS-E2-R",
        "IPSL-CM5B-LR", "MIROC5-ESM-CHEM", "MOHC-HadGEM2-AO", "MRI-ESM1"
    ]
elif ensemble == "CMIP5":
    root = "/mnt/data/Other/GCM-RCM/GCM/CMIP5/monthly_nc"
    output_txt_dir = f"/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP5/{variable}"
    models = [
        "bcc-csm1-1", "CESM1-BGC", "CESM1-CAM5", "CMCC-CESM", "CMCC-CMS",
        "CNRM-CM5", "EC-EARTH", "GFDL-CM3", "GFDL-ESM2G", "HadGEM2-ES",
        "inmcm4", "IPSL-CM5A-MR", "MPI-ESM-LR", "MRI-ESM1", "NorESM1-M"
    ]
elif ensemble == "CMIP6":
    root = "/mnt/data/Other/GCM-RCM/GCM/CMIP6/monthly_nc/ScenarioMIP"
    output_txt_dir = f"/mnt/data/Other/GCM-RCM/climatology_tabs/CMIP6/{variable}"
    models = [
        "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1-HR",
        "EC-Earth3", "EC-Earth3-CC", "EC-EARTH3-Veg", "GFDL-CM4", "GFDL-ESM4", "HadGEM3-GC31-MM",
        "INM-CM5-0", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM", "TaiESM1"
    ]
else:
    raise ValueError("Invalid ensemble type. Choose 'CMIP5_4_EURO-CORDEX', 'CMIP5_OTHERS', or 'CMIP6'.")

# Create output directory if it doesn't exist
os.makedirs(output_txt_dir, exist_ok=True)

scenarios = ["historical", "ssp585", "RCP85", "rcp85"]

# Time windows per scenario
time_windows = {
    "historical": ("1981-01-01", "2005-12-31"),
    "ssp585": ("2076-01-01", "2100-12-31"),
    "RCP85": ("2076-01-01", "2100-12-31"),
    "rcp85": ("2076-01-01", "2100-12-31")
}

# Processing loop
for scenario in scenarios:
    for model in models:
        model_dir = os.path.join(root, scenario, model)
        if not os.path.isdir(model_dir):
            print(f"⏩ Skipping: {model} {scenario} — directory not found")
            continue

        start_date, end_date = time_windows[scenario]

        # Match all variants
        if ensemble == "EURO-CORDEX_EUR44":
            pattern = f"{variable}_EUR-44_*_{scenario}_r*_Amon_*.nc"
        else:
            pattern = f"{variable}_Amon_{model}_{scenario}_r*.nc"

        input_files = sorted(glob.glob(os.path.join(model_dir, pattern)))
        if not input_files:
            print(f"⏩ No files for {model} {scenario}")
            continue

        print(f"▶ Processing: {model} - {scenario}")

        # Extract variant from the first file
        first_file = os.path.basename(input_files[0])
        if ensemble == "EURO-CORDEX_EUR44":
            variant_match = re.search(r"_(r\d+i\d+p\d+(?:f\d+)?)(?:_|\.)", first_file)
        else:
            variant_match = re.search(r"_(r\d+i\d+p\d+(f\d+)?)(?:_|\.nc)", first_file)

        variant = variant_match.group(1) if variant_match else "unknown"

        # Output file base name and .txt path
        if ensemble == "EURO-CORDEX_EUR44":
            gcm_rcm_parts = model.split("_")
            if len(gcm_rcm_parts) == 4:
                gcm, variant, rcm, rcm_ver = gcm_rcm_parts
                output_base = f"{variable}_Amon_{gcm}_{rcm}_{rcm_ver}_{scenario}_{variant}_{start_date.replace('-', '')}-{end_date.replace('-', '')}"
            else:
                print(f"⚠️ Unexpected model format: {model}")
                output_base = f"{variable}_Amon_{model}_{scenario}_{variant}_{start_date.replace('-', '')}-{end_date.replace('-', '')}"
        else:
            output_base = f"{variable}_Amon_{model}_{scenario}_{variant}_{start_date.replace('-', '')}-{end_date.replace('-', '')}"
        info_txt = os.path.join(output_txt_dir, output_base + ".txt")

        # Temp files
        merged = os.path.join(model_dir, "temp_merged.nc")
        selected_time = os.path.join(model_dir, "temp_selected.nc")
        selected_region = os.path.join(model_dir, "temp_region.nc")
        climatology = os.path.join(model_dir, "temp_climatology.nc")

        try:
            # Run CDO processing steps
            subprocess.run(["cdo", "mergetime"] + input_files + [merged], check=True)
            subprocess.run(["cdo", f"seldate,{start_date},{end_date}", merged, selected_time], check=True)
            subprocess.run(["cdo", f"sellonlatbox,{lonmin},{lonmax},{latmin},{latmax}", selected_time, selected_region], check=True)
            subprocess.run(["cdo", "ymonmean", selected_region, climatology], check=True)

            # Extract metadata
            result = subprocess.run(["cdo", "infon", climatology], capture_output=True, text=True, check=True)
            lines = [line for line in result.stdout.splitlines() if "average_DT" not in line]

            # Write to output .txt
            with open(info_txt, "w") as out_file:
                out_file.write("\n".join(lines))

            print(f"   ✅ Saved: {info_txt}")

        except subprocess.CalledProcessError as e:
            print(f"   ❌ CDO command failed for {model} {scenario}: {e}")

        finally:
            # Clean up temp files
            for f in [merged, selected_time, selected_region, climatology]:
                if os.path.exists(f):
                    os.remove(f)
