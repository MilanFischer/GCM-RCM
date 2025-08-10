import cdsapi
client = cdsapi.Client()
import requests
import os
os.chdir("/mnt/data/Other/GCM-RCM/CDS_downloads/RCP85")

# Downloading data using the CDS API

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Paste the following  code copied from the the online API request
dataset = "projections-cmip5-daily-single-levels"
request = {
    "experiment": "rcp_8_5",
    "variable": [
        "2m_temperature",
        "maximum_2m_temperature_in_the_last_24_hours",
        "mean_sea_level_pressure",
        "minimum_2m_temperature_in_the_last_24_hours",
        "near_surface_specific_humidity",
        "daily_near_surface_relative_humidity"
    ],
    "model": "noresm1_m",
    "ensemble_member": "r1i1p1",
    "period": ["20560101-21001231"]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#--------------------------------------------------------------------
import zipfile
import os
import shutil

# Step 1: Find and unzip the .zip file
zip_file = next((f for f in os.listdir() if f.endswith(".zip")), None)
if not zip_file:
    raise FileNotFoundError("No .zip file found.")

with zipfile.ZipFile(zip_file, 'r') as zip_ref:
    zip_ref.extractall()
    extracted_files = zip_ref.namelist()
print(f"✅ Unzipped: {zip_file}")

# Step 2: Delete the zip file
os.remove(zip_file)
print(f"🗑️ Deleted: {zip_file}")

# Step 3: Identify model name from extracted filename
# (Assumes consistent naming pattern)
if not extracted_files:
    raise RuntimeError("No files were extracted from the zip archive.")

model_name = extracted_files[0].split("_")[2]
target_dir = os.path.join(os.getcwd(), model_name)
os.makedirs(target_dir, exist_ok=True)

# Step 4: Move extracted files to the model-named directory
for fname in extracted_files:
    shutil.move(fname, os.path.join(target_dir, fname))
    print(f"📦 Moved {fname} → {target_dir}/")

print(f"✅ All files moved to folder: {model_name}")