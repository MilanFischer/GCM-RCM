#!/usr/bin/env python3
import argparse
import logging
import os
import shutil
from typing import Dict, List, Optional

# ---------------------------
# Logging
# ---------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("copy_cmip5_hfls")

# ---------------------------
# Defaults / Paths
# ---------------------------
DEFAULT_SRC_ROOTS = [
    # Preference order: EURO-CORDEX first, then OTHERS
    "/mnt/k_drive/Data/GCM/CMIP5_4_EURO-CORDEX/monthly_nc",
    "/mnt/k_drive/Data/GCM/CMIP5_OTHERS/monthly_nc",
]
DEFAULT_DST_ROOT = "/mnt/data/Other/GCM-RCM/GCM/CMIP5/monthly_nc"

# Target scenarios and common aliases
SCENARIO_ALIASES = {
    "historical": "historical",
    "rcp85": "RCP85",
    "rcp8.5": "RCP85",
    "RCP85": "RCP85",
    "HISTORICAL": "historical",
}

# If --models is not provided, we’ll read existing model folder names under the destination
DEFAULT_FALLBACK_MODELS = [
    # A reasonable superset; actual list is loaded from target if present
    "bcc-csm1-1", "CESM1-BGC", "CESM1-CAM5", "CMCC-CESM", "CMCC-CMS", "CNRM-CM5",
    "EC-EARTH", "GFDL-CM3", "GFDL-ESM2G", "HadGEM2-ES", "inmcm4", "IPSL-CM5A-MR",
    "MPI-ESM-LR", "MRI-ESM1", "NorESM1-M",
]

def norm(s: str) -> str:
    """Normalize model names for case-insensitive matching (strip non-alnum)."""
    return "".join(ch.lower() for ch in s if ch.isalnum())

def discover_target_models(dst_root: str, scenario_dir: str) -> List[str]:
    """Read existing model subfolders from destination scenario directory, if it exists."""
    path = os.path.join(dst_root, scenario_dir)
    if not os.path.isdir(path):
        return []
    entries = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return sorted(entries)

def build_target_name_map(existing_models: List[str]) -> Dict[str, str]:
    """
    Map normalized name -> actual folder name in destination.
    This lets us place files in existing lowercase/cased dirs when names differ.
    """
    m: Dict[str, str] = {}
    for name in existing_models:
        m[norm(name)] = name
    return m

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def copy_hfls_folder(src_dir: str, dst_dir: str, overwrite: bool, dry_run: bool) -> int:
    """Copy all hfls_*.nc files from src_dir to dst_dir."""
    if not os.path.isdir(src_dir):
        return 0
    ensure_dir(dst_dir)
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
    return moved

def first_existing_dir(*candidates: str) -> Optional[str]:
    for d in candidates:
        if d and os.path.isdir(d):
            return d
    return None

def gather_source_models(src_root: str, scenario_dir: str) -> List[str]:
    path = os.path.join(src_root, scenario_dir)
    if not os.path.isdir(path):
        return []
    return sorted([d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))])

def main():
    ap = argparse.ArgumentParser(
        description="Copy CMIP5 monthly hfls (latent heat flux) files into the CMIP5 GCM-RCM tree."
    )
    ap.add_argument(
        "--src-roots",
        nargs="+",
        default=DEFAULT_SRC_ROOTS,
        help="Source roots searched in order of priority (default: EURO-CORDEX then OTHERS).",
    )
    ap.add_argument(
        "--dst-root",
        default=DEFAULT_DST_ROOT,
        help="Destination root containing <scenario>/<model>/ folders.",
    )
    ap.add_argument(
        "--scenarios",
        nargs="+",
        default=["historical", "RCP85"],
        help="Scenarios to process (aliases: rcp85/rcp8.5).",
    )
    ap.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="Models to copy. If omitted, we read existing model folders in the destination scenario dir; "
             "if none, we fall back to a generic list.",
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

    # Normalize scenarios via aliases
    scenarios = []
    for sc in args.scenarios:
        scenarios.append(SCENARIO_ALIASES.get(sc, sc))

    total = 0

    for scenario in scenarios:
        # Destination scenario dir, e.g., .../monthly_nc/historical or .../monthly_nc/RCP85
        dst_scenario_dir = os.path.join(args.dst_root, scenario)
        ensure_dir(dst_scenario_dir)

        # Determine models list
        if args.models:
            models = args.models
        else:
            existing = discover_target_models(args.dst_root, scenario)
            models = existing if existing else DEFAULT_FALLBACK_MODELS
            if existing:
                L.info(f"Using existing destination model folders for {scenario}: {', '.join(existing)}")
            else:
                L.info(f"No destination model folders found for {scenario}; using fallback model list.")

        # Build a case-insensitive mapping to existing model folder names (to preserve your casing)
        target_name_map = build_target_name_map(discover_target_models(args.dst_root, scenario))

        # Gather available source models per root to skip empties faster
        src_models_by_root = {root: set(gather_source_models(root, scenario)) for root in args.src_roots}

        L.info(f"=== Scenario: {scenario} ===")
        for model in models:
            # Decide the destination model folder name (preserve existing casing if present)
            dst_model_name = target_name_map.get(norm(model), model)
            dst_model_dir = os.path.join(dst_scenario_dir, dst_model_name)

            # Try each source root in priority order
            moved = 0
            for src_root in args.src_roots:
                src_models = src_models_by_root.get(src_root, set())
                # Prefer exact match; else try case-insensitive match against available source dirs
                if model in src_models:
                    src_model_name = model
                else:
                    # find any case-insensitive match
                    matches = [m for m in src_models if norm(m) == norm(model)]
                    src_model_name = matches[0] if matches else None

                if not src_model_name:
                    continue  # not in this root

                src_model_dir = os.path.join(src_root, scenario, src_model_name)
                L.info(f"{scenario}/{dst_model_name}: searching in {src_model_dir}")
                moved_here = copy_hfls_folder(src_model_dir, dst_model_dir, args.overwrite, args.dry_run)
                moved += moved_here

                # If we copied anything from the preferred root, we don't fall back
                if moved_here > 0:
                    break

            if moved == 0:
                L.info(f"No hfls files found for {scenario}/{model} in any source roots (or all skipped).")
            total += moved

    L.info(f"Done. Files processed (copied or skipped due to overwrite policy): {total}")

if __name__ == "__main__":
    main()

# Example usage:
# python copy_hfls_CMIP5.py --scenarios historical # Copy historical scenario
# python copy_hfls_CMIP5.py --scenarios rcp85 # Copy RCP85 scenario
# python copy_hfls_CMIP5.py --dry-run # Show what would be copied without writing
# python copy_hfls_CMIP5.py --overwrite # Overwrite existing files
# python copy_hfls_CMIP5.py --scenarios historical RCP85 --models EC-EARTH IPSL-CM5A-MR MPI-ESM-LR # Filter specific models
