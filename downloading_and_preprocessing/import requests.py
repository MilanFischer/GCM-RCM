import requests

# Step 1: Get the dataset's metadata
collection_url = "https://cds.climate.copernicus.eu/api/catalogue/v1/collections/projections-cmip5-daily-single-levels"
collection_metadata = requests.get(collection_url).json()

# Step 2: Extract the 'form' link (contains valid query options)
form_url = next(link["href"] for link in collection_metadata["links"] if link["rel"] == "form")

# Step 3: Fetch the form JSON with valid parameter values
form_data = requests.get(form_url).json()

# Step 4: Function to print clean options
def print_options(label, values):
    print(f"\n{label} ({len(values)}):")
    for v in values:
        print("  ", v.get("name", v))

# Step 5: Print out all relevant dimensions
print_options("Models", form_data.get("model", []))
print_options("Experiments", form_data.get("experiment", []))
print_options("Variables", form_data.get("variable", []))
print_options("Ensemble members", form_data.get("ensemble_member", []))
print_options("Periods", form_data.get("period", []))


# https://github.com/ecmwf/cdsapi/issues/80
import requests

pressure_level_address = "https://cds.climate.copernicus.eu/api/catalogue/v1/collections/reanalysis-era5-pressure-levels"

r = requests.get(pressure_level_address)
time_interval = r.json()["extent"]["temporal"]["interval"][0]
end_time = time_interval[1]









# Step 1: Get the dataset's metadata
collection_url = "https://cds.climate.copernicus.eu/api/catalogue/v1/collections/projections-cmip5-daily-single-levels"
collection_metadata = requests.get(collection_url).json()

# Step 2: Extract the 'form' link (contains valid query options)
form_url = next(link["href"] for link in collection_metadata["links"] if link["rel"] == "form")

# Step 3: Fetch the form JSON with valid parameter values
form_data = requests.get(form_url).json()