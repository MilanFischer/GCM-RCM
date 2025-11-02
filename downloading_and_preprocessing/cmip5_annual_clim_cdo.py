#!/usr/bin/env python3
import argparse
import logging
import os
import glob
import subprocess
from typing import List

# -----------------------
# Logging
# -----------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("cmip5_annual_clim_cdo")

# -----------------------
# Paths & Defaults
# -----------------------
BASE_CMIP5 = "/mnt/data/Other/GCM-RCM/GCM/CMIP5"
SRC_ROOT = os.path.join(BASE_CMIP5, "monthly_nc")
DST_ROOT = os.path.join(BASE_CMIP5, "annual_nc")  # single-layer annual climatology files

SCENARIOS = {
    "historical": {"years": (1981, 2005)},
    "RCP85": {"years": (2076, 2100)},
}

VARS = ("hfls", "vpd")  # compute both if present

# -----------------------
# Helpers
# -----------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def run_cdo(args: List[str], **kwargs):
    cmd = ["cdo"] + args
    L.debug("RUN: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)

def find_var_files(model_dir: str, var: str, scenario: str) -> List[str]:
    """
    Find CMIP5 monthly files for a variable inside model_dir.
    Robust to: RCP85 vs rcp85, small naming differences, and missing scenario token.
    """
    scen_tokens = [scenario, scenario.lower(), scenario.upper()]
    patterns = []
    for scen in scen_tokens:
        patterns += [
            os.path.join(model_dir, f"{var}_Amon_*_{scen}_*.nc"),
            os.path.join(model_dir, f"{var}_*_{scen}_*.nc"),
            os.path.join(model_dir, f"*_{var}_*_{scen}_*.nc"),
        ]
    # Fallbacks (no scenario token)
    patterns += [
        os.path.join(model_dir, f"{var}_Amon_*.nc"),
        os.path.join(model_dir, f"{var}_*.nc"),
    ]
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))
    return sorted(set(files))

def output_name(var: str, model: str, scenario: str, yr0: int, yr1: int) -> str:
    # Single layer = mean of annual means over the period
    return f"{var}_AyrClim_{model}_{scenario}_{yr0}-{yr1}.nc"

def annual_clim_for_var(
    src_files: List[str],
    out_file: str,
    year0: int,
    year1: int,
    overwrite: bool,
    dry_run: bool,
):
    """
    Produce one-layer file:
      timmean(yearmean(sel(year0/year1, mergetime(files))))
    """
    if not src_files:
        L.info("  - no source files found (skip)")
        return
    if os.path.exists(out_file) and not overwrite:
        L.info(f"  - exists, skip: {out_file}")
        return

    workdir = os.path.dirname(src_files[0])
    merged = os.path.join(workdir, "tmp_merge.nc")
    sel = os.path.join(workdir, "tmp_sel.nc")
    yrmean = os.path.join(workdir, "tmp_yearmean.nc")

    try:
        if dry_run:
            L.info(f"  [dry-run] mergetime {len(src_files)} files")
            L.info(f"  [dry-run] selyear {year0}/{year1}")
            L.info(f"  [dry-run] yearmean then timmean -> {out_file}")
            return

        cdo_opts = ["-f", "nc4c", "-z", "zip"]

        # 1) merge time if needed
        if len(src_files) == 1:
            L.info("  - single file: skipping mergetime")
            run_cdo(cdo_opts + ["copy", src_files[0], merged])
        else:
            L.info(f"  - mergetime {len(src_files)} files -> {merged}")
            run_cdo(cdo_opts + ["mergetime"] + src_files + [merged])

        # 2) select years
        L.info(f"  - selyear {year0}/{year1} -> tmp_sel.nc")
        run_cdo(cdo_opts + [f"selyear,{year0}/{year1}", merged, sel])

        # 3) annual means (per year)
        L.info("  - yearmean -> tmp_yearmean.nc")
        run_cdo(cdo_opts + ["yearmean", sel, yrmean])

        # 4) mean over all annual means (single layer)
        L.info(f"  - timmean(yearmean) -> {out_file}")
        ensure_dir(os.path.dirname(out_file))
        run_cdo(cdo_opts + ["timmean", yrmean, out_file])

        # Optional sanity checks:
        # subprocess.run(["cdo", "-s", "-outputkey,ntime", "-info", out_file], check=True)  # should be 1

    finally:
        for f in (merged, sel, yrmean):
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception:
                pass

def main():
    ap = argparse.ArgumentParser(
        description="Compute CMIP5 single-layer annual climatology (mean of annual means) for hfls and vpd using CDO."
    )
    ap.add_argument("--src-root", default=SRC_ROOT, help="Source base: .../CMIP5/monthly_nc")
    ap.add_argument("--dst-root", default=DST_ROOT, help="Destination base: .../CMIP5/annual_nc (single-layer outputs)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    ap.add_argument("--dry-run", action="store_true", help="Show actions without running CDO")
    ap.add_argument("--scenarios", nargs="+", default=list(SCENARIOS.keys()),
                    help="Scenarios to process (default: historical RCP85)")
    ap.add_argument("--vars", nargs="+", default=list(VARS), help="Variables to process (default: hfls vpd)")
    args = ap.parse_args()

    for scenario in args.scenarios:
        if scenario not in SCENARIOS:
            L.warning(f"Unknown scenario '{scenario}', skipping.")
            continue

        yr0, yr1 = SCENARIOS[scenario]["years"]
        src_scen_dir = os.path.join(args.src_root, scenario)
        dst_scen_dir = os.path.join(args.dst_root, scenario)

        if not os.path.isdir(src_scen_dir):
            L.warning(f"Missing source scenario dir: {src_scen_dir}")
            continue

        models = sorted(
            d for d in os.listdir(src_scen_dir)
            if os.path.isdir(os.path.join(src_scen_dir, d))
        )
        L.info(f"=== {scenario}: {len(models)} models, years {yr0}-{yr1} ===")
        ensure_dir(dst_scen_dir)

        for model in models:
            L.info(f"- {scenario}/{model}")
            src_model_dir = os.path.join(src_scen_dir, model)
            dst_model_dir = os.path.join(dst_scen_dir, model)
            ensure_dir(dst_model_dir)

            for var in args.vars:
                files = find_var_files(src_model_dir, var, scenario)
                if not files:
                    L.info(f"  - {var}: no files found, skipping.")
                    continue

                out_file = os.path.join(dst_model_dir, output_name(var, model, scenario, yr0, yr1))
                try:
                    annual_clim_for_var(files, out_file, yr0, yr1, args.overwrite, args.dry_run)
                except subprocess.CalledProcessError as e:
                    L.error(f"  ❌ CDO failed for {scenario}/{model}/{var}: {e}")
                except Exception as e:
                    L.exception(f"  ❌ Unexpected error for {scenario}/{model}/{var}: {e}")

    L.info("Done.")

if __name__ == "__main__":
    main()

# Example usage:
# python cmip5_annual_clim_cdo.py                       # process both scenarios with defaults
# python cmip5_annual_clim_cdo.py --dry-run             # show planned operations only
# python cmip5_annual_clim_cdo.py --overwrite           # recompute and overwrite outputs
# python cmip5_annual_clim_cdo.py --scenarios historical --vars hfls  # only historical hfls
