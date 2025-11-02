import argparse
import logging
import os
import shutil
from typing import List

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("copy_hfls")

DEFAULT_SCENARIOS = ["historical", "ssp585"]

# The model lists you showed under your target dirs (you can override via --models)
DEFAULT_MODELS = [
    "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2",
    "CNRM-CM6-1-HR", "EC-Earth3", "EC-Earth3-CC", "EC-EARTH3-Veg", "GFDL-CM4",
    "GFDL-ESM4", "HadGEM3-GC31-MM", "INM-CM5-0", "IPSL-CM6A-LR", "MPI-ESM1-2-HR",
    "MRI-ESM2-0", "NorESM2-MM", "TaiESM1"
]

def copy_hfls_for_model(
    src_root: str,
    dst_root: str,
    scenario: str,
    model: str,
    overwrite: bool = False,
    dry_run: bool = False,
) -> int:
    """
    Copy hfls_*.nc files from:
      <src_root>/<scenario>/<model>/
    to:
      <dst_root>/<scenario>/<model>/
    Returns number of files copied/skipped (existing).
    """
    src_dir = os.path.join(src_root, scenario, model)
    if not os.path.isdir(src_dir):
        L.warning(f"Missing source dir, skipping: {src_dir}")
        return 0

    dst_dir = os.path.join(dst_root, scenario, model)
    os.makedirs(dst_dir, exist_ok=True)

    moved = 0
    for fname in sorted(os.listdir(src_dir)):
        if not (fname.startswith("hfls_") and fname.endswith(".nc")):
            continue
        src = os.path.join(src_dir, fname)
        dst = os.path.join(dst_dir, fname)

        if os.path.exists(dst) and not overwrite:
            L.info(f"exists, skip: {dst}")
            continue

        if dry_run:
            L.info(f"[dry-run] copy {src} -> {dst}")
        else:
            shutil.copy2(src, dst)
            L.info(f"copied: {dst}")
        moved += 1

    if moved == 0:
        L.info(f"No hfls files found for {scenario}/{model} (or all skipped).")
    return moved

def main():
    ap = argparse.ArgumentParser(
        description="Copy CMIP6 monthly hfls (latent heat flux) files into the GCM-RCM tree."
    )
    ap.add_argument(
        "--src-root",
        default="/mnt/k_drive/Data/GCM/CMIP6/monthly_nc/ScenarioMIP",
        help="Source root that contains <scenario>/<model>/hfls_*.nc",
    )
    ap.add_argument(
        "--dst-root",
        default="/mnt/data/Other/GCM-RCM/GCM/CMIP6/monthly_nc/ScenarioMIP",
        help="Destination root to place <scenario>/<model>/hfls_*.nc",
    )
    ap.add_argument(
        "--scenarios",
        nargs="+",
        default=DEFAULT_SCENARIOS,
        help="Scenarios to process (default: historical ssp585)",
    )
    ap.add_argument(
        "--models",
        nargs="+",
        default=DEFAULT_MODELS,
        help="Models to copy (default: your target set). If omitted, uses this list.",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing destination files.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be copied without writing.",
    )
    args = ap.parse_args()

    total = 0
    for scenario in args.scenarios:
        L.info(f"=== Scenario: {scenario} ===")
        for model in args.models:
            try:
                total += copy_hfls_for_model(
                    args.src_root, args.dst_root, scenario, model,
                    overwrite=args.overwrite, dry_run=args.dry_run
                )
            except Exception as e:
                L.exception(f"Error on {scenario}/{model}: {e}")

    L.info(f"Done. Files processed (copied or skipped due to overwrite policy): {total}")

if __name__ == "__main__":
    main()

# Example usage:
# python copy_hfls_CMIP6.py --scenarios historical # Copy historical scenario
# python copy_hfls_CMIP6.py --scenarios ssp585 # Copy ssp585 scenario
# python copy_hfls_CMIP6.py --dry-run # Show what would be copied without writing
# python copy_hfls_CMIP6.py --overwrite # Overwrite existing files
# python copy_hfls_CMIP6.py --scenarios historical ssp585 --models CESM2 IPSL-CM6A-LR MPI-ESM1-2-HR # Filter specific models