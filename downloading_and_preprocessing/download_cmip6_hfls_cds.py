"""
Download CMIP6 monthly 'hfls' (surface upward latent heat flux) from the CDS
into: <dst_root>/ScenarioMIP/<scenario>/<model>/

- dataset: 'projections-cmip6'
- variable: 'surface_upward_latent_heat_flux'
- experiment: 'historical' or 'ssp5_8_5'
- model: CDS-style id (e.g., 'awi_cm_1_1_mr')
"""

import argparse
import logging
import os
import shutil
import tempfile
import zipfile
from typing import List, Optional

import cdsapi

# ---------------------------
# Logging
# ---------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("download_cds_cmip6_hfls")

# ---------------------------
# Defaults
# ---------------------------
DEFAULT_SCENARIOS = ["historical", "ssp585"]
DEFAULT_MODELS = [
    "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2",
    "CNRM-CM6-1-HR", "EC-Earth3", "EC-Earth3-CC", "EC-EARTH3-Veg", "GFDL-CM4",
    "GFDL-ESM4", "HadGEM3-GC31-MM", "INM-CM5-0", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
    "MRI-ESM2-0", "NorESM2-MM", "TaiESM1"
]

# Your usual windows
DEFAULT_YEARS = {
    "historical": (1981, 2005),
    "ssp585": (2076, 2100),
}

# Scenario mapping to CDS "experiment"
EXPERIMENT_MAP = {
    "historical": "historical",
    "ssp585": "ssp5_8_5",  # CDS uses underscores
}

# Correct CDS variable id for CMIP6 HFLs
CDS_VARIABLE = "surface_upward_latent_heat_flux"

MONTHS_12 = [f"{m:02d}" for m in range(1, 13)]


def model_to_cds(model: str) -> str:
    """
    Convert your folder-style model name to CDS id heuristically:
    'EC-Earth3' -> 'ec_earth3', 'IPSL-CM6A-LR' -> 'ipsl_cm6a_lr'
    You can override per-model via --cds-models.
    """
    return model.replace("-", "_").replace(".", "_").lower()


def year_list(start: int, end: int) -> List[str]:
    return [str(y) for y in range(start, end + 1)]


def unzip_all(zip_path: str, out_dir: str):
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)


def move_netcdf(src_dir: str, dst_dir: str, overwrite: bool) -> int:
    os.makedirs(dst_dir, exist_ok=True)
    moved = 0
    for root, _, files in os.walk(src_dir):
        for fn in files:
            if not fn.endswith(".nc"):
                continue
            src = os.path.join(root, fn)
            dst = os.path.join(dst_dir, fn)
            if os.path.exists(dst) and not overwrite:
                L.info(f"exists, skip: {dst}")
                continue
            shutil.copy2(src, dst)
            L.info(f"saved: {dst}")
            moved += 1
    return moved


def retrieve_block(
    client: cdsapi.Client,
    experiment: str,
    cds_model: str,
    years: List[str],
    months: List[str],
    target_zip: str,
):
    """
    Make one CDS request for a list of years (monthly resolution).
    """
    req = {
        "temporal_resolution": "monthly",
        "experiment": experiment,
        "variable": CDS_VARIABLE,
        "model": cds_model,
        "year": years,
        "month": months,
    }
    L.info(f"CDS retrieve: exp={experiment} model={cds_model} years={years[0]}..{years[-1]}")
    client.retrieve("projections-cmip6", req, target_zip)


def download_model_scenario(
    dst_root: str,
    scenario: str,
    model_display: str,
    years: List[str],
    overwrite: bool,
    dry_run: bool,
    cds_model_name: Optional[str] = None,
    block_size: int = 10,  # split requests into <=10-year chunks for robustness
):
    experiment = EXPERIMENT_MAP.get(scenario.lower())
    if not experiment:
        L.warning(f"Unrecognized scenario '{scenario}', skipping.")
        return

    cds_model = cds_model_name or model_to_cds(model_display)
    out_dir = os.path.join(dst_root, "ScenarioMIP", scenario, model_display)
    os.makedirs(out_dir, exist_ok=True)

    if dry_run:
        L.info(f"[dry-run] would save to: {out_dir}")
        L.info(f"[dry-run] CDS: exp={experiment}, model={cds_model}, years={years[0]}..{years[-1]}")
        return

    c = cdsapi.Client()
    moved_total = 0

    with tempfile.TemporaryDirectory() as tmp:
        # chunk years (e.g., 10-year blocks)
        for i in range(0, len(years), block_size):
            yrs = years[i : i + block_size]
            zip_path = os.path.join(tmp, f"{cds_model}_{experiment}_{yrs[0]}_{yrs[-1]}.zip")
            try:
                retrieve_block(c, experiment, cds_model, yrs, MONTHS_12, zip_path)
            except Exception as e:
                L.exception(f"CDS request failed for {scenario}/{model_display} years {yrs[0]}–{yrs[-1]}: {e}")
                continue

            unpack_dir = os.path.join(tmp, f"unz_{yrs[0]}_{yrs[-1]}")
            os.makedirs(unpack_dir, exist_ok=True)
            try:
                unzip_all(zip_path, unpack_dir)
                moved = move_netcdf(unpack_dir, out_dir, overwrite)
                moved_total += moved
            except Exception as e:
                L.exception(f"Unzip/move failed for {zip_path}: {e}")

    if moved_total == 0:
        L.info(f"No files saved for {scenario}/{model_display}. Check CDS model id or availability.")
    else:
        L.info(f"Done {scenario}/{model_display}: {moved_total} files saved.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Download CMIP6 monthly HFLs from CDS into the CMIP6 tree."
    )
    p.add_argument("--dst-root",
                   default="/mnt/data/Other/GCM-RCM/GCM/CMIP6/monthly_nc",
                   help="Destination base; files go under <dst>/ScenarioMIP/<scenario>/<model>/")
    p.add_argument("--scenarios", nargs="+", default=DEFAULT_SCENARIOS,
                   help="historical, ssp585")
    p.add_argument("--models", nargs="+", default=DEFAULT_MODELS,
                   help="Target model folder names (used in destination).")
    p.add_argument("--start-year", type=int, help="Override start year (all scenarios).")
    p.add_argument("--end-year", type=int, help="Override end year (all scenarios).")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    p.add_argument("--dry-run", action="store_true", help="Show actions without downloading.")
    p.add_argument("--cds-models", nargs="+",
                   help="Optional explicit CDS model ids (1:1 with --models), e.g., ec_earth3 ipsl_cm6a_lr ...")
    p.add_argument("--block-size", type=int, default=10,
                   help="Years per CDS request (default 10).")
    return p.parse_args()


def main():
    args = parse_args()

    for scenario in args.scenarios:
        sy, ey = DEFAULT_YEARS.get(scenario.lower(), (None, None))
        start = args.start_year if args.start_year is not None else sy
        end = args.end_year if args.end_year is not None else ey
        if start is None or end is None:
            # Fallback to canonical CMIP6 windows if not provided
            start, end = (1850, 2014) if scenario.lower() == "historical" else (2015, 2100)

        yrs = year_list(start, end)
        L.info(f"=== {scenario}: {yrs[0]}–{yrs[-1]} ===")

        cds_models = args.cds_models or [None] * len(args.models)
        if len(cds_models) != len(args.models):
            raise SystemExit("--cds-models must have same length as --models")

        for model_display, cds_name in zip(args.models, cds_models):
            try:
                download_model_scenario(
                    dst_root=args.dst_root,
                    scenario=scenario,
                    model_display=model_display,
                    years=yrs,
                    overwrite=args.overwrite,
                    dry_run=args.dry_run,
                    cds_model_name=cds_name,
                    block_size=args.block_size,
                )
            except Exception as e:
                L.exception(f"Failed {scenario}/{model_display}: {e}")


if __name__ == "__main__":
    main()

# Example usage:
# python download_cds_cmip6_hfls.py --scenarios historical        # 1981–2005 by default
# python download_cds_cmip6_hfls.py --scenarios ssp585            # 2076–2100 by default
# python download_cds_cmip6_hfls.py --dry-run                     # Show planned requests only
# python download_cds_cmip6_hfls.py --overwrite                   # Replace existing files
# python download_cds_cmip6_hfls.py --scenarios historical ssp585 --models CESM2 IPSL-CM6A-LR MPI-ESM1-2-HR
# python download_cds_cmip6_hfls.py --scenarios ssp585 --models EC-Earth3 --cds-models ec_earth3
