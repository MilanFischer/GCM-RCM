#!/usr/bin/env python3
import argparse
import logging
import os
import glob
import re
import xarray as xr
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("cmip5_geff_from_annual_fixed")

GAMMA_KPA_PER_K = 0.066
RHO_AIR = 1.225
CP_AIR = 1004.67

BASE_CMIP5 = "/mnt/data/Other/GCM-RCM/GCM/CMIP5"
SRC_ROOT = os.path.join(BASE_CMIP5, "annual_nc")
DST_ROOT = os.path.join(BASE_CMIP5, "annual_nc")
SCENARIOS = ["historical", "RCP85"]

# ---------- helpers ----------
def ensure_dir(p): os.makedirs(p, exist_ok=True)

def find_file(model_dir: str, var: str, scenario: str):
    scen_tokens = [scenario, scenario.lower(), scenario.upper()]
    pats = []
    for scen in scen_tokens:
        pats += [
            os.path.join(model_dir, f"{var}_Amon_*_{scen}_yearmean_*.nc"),
            os.path.join(model_dir, f"{var}_*_{scen}_yearmean_*.nc"),
            os.path.join(model_dir, f"*_{var}_*_{scen}_*yearmean*.nc"),
        ]
    pats += [os.path.join(model_dir, f"{var}_*yearmean*.nc")]
    hits = []
    for p in pats: hits.extend(glob.glob(p))
    hits = sorted(set(hits))
    return hits[0] if hits else None

def pick_numeric_var(ds: xr.Dataset, prefer_tokens: list[str]) -> xr.DataArray:
    # 1) exact match
    for tok in prefer_tokens:
        if tok in ds.data_vars and np.issubdtype(ds[tok].dtype, np.number):
            return ds[tok]
    # 2) name contains token & is numeric
    for tok in prefer_tokens:
        for name in ds.data_vars:
            if tok in name.lower() and np.issubdtype(ds[name].dtype, np.number):
                if any(bad in name.lower() for bad in ["time", "bnds", "bounds", "average_t"]):
                    continue
                return ds[name]
    # 3) first numeric variable that isn’t a time/bounds helper
    for name, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.number) and not any(
            bad in name.lower() for bad in ["time", "bnds", "bounds", "average_t"]
        ):
            return da
    raise ValueError(f"No suitable numeric variable found. Vars: {list(ds.data_vars)}")

def maybe_convert_vpd_to_kpa(vpd: xr.DataArray) -> xr.DataArray:
    units = str(vpd.attrs.get("units", "")).lower()
    if units in ("pa", "pascal"):
        L.info("    vpd in Pa → converting to kPa")
        return (vpd.astype("float64") / 1000.0).assign_attrs(units="kPa")
    if "kpa" in units or units == "kilopascal":
        return vpd.astype("float64")
    return vpd.astype("float64")

def maybe_check_hfls_units(hfls: xr.DataArray) -> xr.DataArray:
    # Expect W m-2 (J s-1 m-2). Proceed but cast to float.
    return hfls.astype("float64")

def compute_geff(vpd: xr.DataArray, hfls: xr.DataArray) -> xr.DataArray:
    vpd, hfls = xr.align(vpd, hfls, join="inner")
    vpd = maybe_convert_vpd_to_kpa(vpd)
    hfls = maybe_check_hfls_units(hfls)

    rseff = (RHO_AIR * CP_AIR / GAMMA_KPA_PER_K) * (vpd / hfls)
    geff = (1.0 / rseff) * 1000.0
    geff.name = "geff"
    geff.attrs.update({
        "long_name": "effective conductance",
        "formula": "geff = 1000 / (rho_air * cp_air / gamma * vpd / hfls)",
        "units": "mm s-1",
        "rho_air": RHO_AIR, "cp_air": CP_AIR, "gamma": GAMMA_KPA_PER_K,
        "expects": "vpd in kPa, hfls in W m-2; annual means aligned on time",
    })
    return geff

def save_geff(geff: xr.DataArray, out_path: str):
    ensure_dir(os.path.dirname(out_path))
    comp = {geff.name: {"zlib": True, "complevel": 4, "_FillValue": None}}
    xr.Dataset({geff.name: geff}, coords=geff.coords).to_netcdf(out_path, encoding=comp)
    L.info(f"    ✓ saved: {out_path}")

def infer_year_window(path: str) -> str:
    """
    Extract 'YYYY-YYYY' from filename, else infer from time axis.
    Works for pandas.DatetimeIndex and cftime calendars.
    """
    base = os.path.basename(path)

    # strict: tail pattern ..._YYYY-YYYY.nc
    m = re.search(r'_(\d{4}-\d{4})(?:\.nc)?$', base)
    if m:
        return m.group(1)

    # loose: anywhere in the name
    m = re.search(r'(\d{4})-(\d{4})', base)
    if m:
        return f"{m.group(1)}-{m.group(2)}"

    # fallback: read from time axis
    try:
        with xr.open_dataset(path, decode_times=True) as ds:
            if "time" in ds.coords:
                idx = ds.indexes.get("time")
                if idx is not None and len(idx) > 0:
                    # idx may be pandas or cftime index; both expose .year
                    years = np.array([getattr(t, "year", None) for t in idx])
                    years = years[years != None]  # filter Nones
                    if years.size:
                        return f"{int(years.min())}-{int(years.max())}"
    except Exception:
        pass

    return "unknown"

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Compute geff (effective conductance) from CMIP5 annual vpd & hfls."
    )
    ap.add_argument("--src-root", default=SRC_ROOT, help=".../CMIP5/annual_nc")
    ap.add_argument("--dst-root", default=DST_ROOT, help=".../CMIP5/annual_nc_geff")
    ap.add_argument("--scenarios", nargs="+", default=SCENARIOS)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for scenario in args.scenarios:
        scen_src = os.path.join(args.src_root, scenario)
        if not os.path.isdir(scen_src):
            L.warning(f"Missing scenario dir: {scen_src}")
            continue

        models = sorted(d for d in os.listdir(scen_src) if os.path.isdir(os.path.join(scen_src, d)))
        L.info(f"=== {scenario}: {len(models)} models ===")

        for model in models:
            L.info(f"- {scenario}/{model}")
            mdir = os.path.join(scen_src, model)
            vpd_file = find_file(mdir, "vpd", scenario)
            hfls_file = find_file(mdir, "hfls", scenario)

            if not vpd_file or not hfls_file:
                if not vpd_file: L.info("    vpd: not found (skip)")
                if not hfls_file: L.info("    hfls: not found (skip)")
                continue

            # derive year window (prefer vpd name; fall back to hfls and then to 'unknown')
            yrwin = infer_year_window(vpd_file)
            if yrwin == "unknown":
                yrwin = infer_year_window(hfls_file)

            out_dir = os.path.join(args.dst_root, scenario, model)
            out_path = os.path.join(out_dir, f"geff_Ayr_{model}_{scenario}_{yrwin}.nc")
            if os.path.exists(out_path) and not args.overwrite:
                L.info(f"    exists, skip: {out_path}")
                continue
            if args.dry_run:
                L.info(f"    [dry-run] vpd:  {os.path.basename(vpd_file)}")
                L.info(f"    [dry-run] hfls: {os.path.basename(hfls_file)}")
                L.info(f"    [dry-run]  → {out_path}")
                continue

            try:
                with xr.open_dataset(vpd_file) as ds_v, xr.open_dataset(hfls_file) as ds_h:
                    vpd_da  = pick_numeric_var(ds_v, ["vpd"])
                    hfls_da = pick_numeric_var(ds_h, ["hfls","hfls_amon","surface_upward_latent_heat_flux"])
                    geff = compute_geff(vpd_da, hfls_da)
                    geff.attrs["vpd_source_file"] = os.path.basename(vpd_file)
                    geff.attrs["hfls_source_file"] = os.path.basename(hfls_file)
                    save_geff(geff, out_path)
            except Exception as e:
                L.exception(f"    ❌ failed {scenario}/{model}: {e}")

    L.info("Done.")

if __name__ == "__main__":
    main()

# Example usage:
# python cmip5_geff_from_annual_fixed.py
# python cmip5_geff_from_annual_fixed.py --dry-run
# python cmip5_geff_from_annual_fixed.py --overwrite
# python cmip5_geff_from_annual_fixed.py --scenarios RCP85
