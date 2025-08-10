from pyesgf.search import SearchConnection
import re
import os
import urllib.request

os.environ["ESGF_PYCLIENT_NO_FACETS_STAR_WARNING"] = "1"

def download_esgf_variable(project, model, time_frequency, realm,
                           ensemble, experiment, variable, latest=True,
                           start_year=None, end_year=None,
                           search_url=None,
                           distrib=False,
                           download_dir="./esgf_data"):
    """
    Download NetCDF files for a specific variable from ESGF.

    Args:
        project (str): e.g., "CMIP5", "CMIP6", "CORDEX"
        model (str): e.g., "BCC-CSM1.1"
        time_frequency (str): e.g., "day"
        realm (str): e.g., "atmos"
        ensemble (str): e.g., "r1i1p1"
        experiment (str): e.g., "historical"
        variable (str): e.g., "tas"
        latest (bool): Only get latest dataset version
        search_url (str or None): ESGF search node; if None, tries fallback list
        distrib (bool): Search all ESGF federation (slower but comprehensive)
        download_dir (str): Where to save downloaded files
    """

    fallback_nodes = [
        "https://esgf.ceda.ac.uk/esg-search",
        "https://esgf-node.llnl.gov/esg-search",
        "https://esgf-data.dkrz.de/esg-search",
        "https://esgf-node.ipsl.upmc.fr/esg-search",
        "https://esgf-data.diasjp.net/esg-search",
        "https://esgf-node.cmcc.it/esg-search"
    ]

    # Try nodes until one works
    active_node = search_url
    if not search_url:
        for node in fallback_nodes:
            try:
                print(f"Trying ESGF node: {node} ...")
                conn = SearchConnection(node, distrib=distrib)
                _ = conn.new_context(project=project).search()  # quick test
                active_node = node
                print(f"✅ Connected to {node}")
                break
            except Exception as e:
                print(f"❌ Failed to connect: {node} ({e})")

        if not active_node:
            print("❌ No ESGF nodes available. Aborting.")
            return

    else:
        conn = SearchConnection(active_node, distrib=distrib)

    print(f"\nUsing ESGF node: {active_node}")
    print(f"Project: {project}, Model: {model}, Experiment: {experiment}, Variable: {variable}\n")

    # Prepare the search context
    ctx = conn.new_context(
        project=project,
        model=model,
        time_frequency=time_frequency,
        realm=realm,
        ensemble=ensemble,
        experiment=experiment,
        latest=latest
    )

    results = ctx.search()
    if not results:
        print("No datasets found.")
        return

    files = results[0].file_context().search()

    # Filter files for the specified variable
    filtered_files = []
    filenames = []
    var_pattern = rf"/{variable}(_|\.|/)"

    for f in files:
        url = f.download_url
        if not url or not re.search(var_pattern, url):
            continue
    
        fname = os.path.basename(url)
        match = re.search(r'(\d{4})\d{2}\d{2}-(\d{4})\d{2}\d{2}', fname)
        if match:
            start_file_year = int(match.group(1))
            end_file_year = int(match.group(2))

            # Include if overlaps with [start_year, end_year]
            if start_year and end_file_year < start_year:
                continue
            if end_year and start_file_year > end_year:
                continue
        
        filtered_files.append(f)
        filenames.append(fname)

    print(f"Found {len(filtered_files)} {variable} files.\n")

    if not filtered_files:
        print("No matching files found.")
        return

    print("Files to download:")
    for name in filenames:
        print(" -", name)

    # Create download directory
    os.makedirs(download_dir, exist_ok=True)

    # Download files with skipping and error logging
    for f in filtered_files:
        url = f.download_url
        filepath = os.path.join(download_dir, os.path.basename(url))

        if os.path.exists(filepath):
            print(f"Skipping (already exists): {filepath}")
            continue

        print(f"Downloading: {url}")
        try:
            urllib.request.urlretrieve(url, filepath)
            print(f"Saved to: {filepath}")
        except Exception as e:
            print(f"❌ Failed: {url}")
            with open("download_errors.log", "a") as log:
                log.write(f"{url} failed: {e}\n")

#--------------
# Example Usage

download_esgf_variable(
    project="CMIP5",
    model="EC-EARTH",
    time_frequency="day",
    realm="atmos",
    ensemble="r1i1p1",
    experiment="historical",
    variable="tas",
    start_year=1981,
    end_year=2005,
    search_url="https://esgf.ceda.ac.uk/esg-search",  # or any other node
    distrib=True,  # enables federation-wide search
    download_dir="./esgf_files"
)



"BCC-CSM1.1", "CESM1(BGC)", "CESM1(CAM5)", "CMCC-CESM", "CMCC-CMS", "CNRM-CM5", "EC-EARTH"