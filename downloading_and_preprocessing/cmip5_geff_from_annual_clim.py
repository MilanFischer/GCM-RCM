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
L = logging.getLogger("cmip5_geff_from_annual_clim")

# ---- constants ----
GAMMA_KPA_PER_K = 0.066
RHO_AIR = 1.225
CP_AIR = 1004.67

# ---- paths ----
BASE_CMIP5 = "/mnt/data/Other/GCM-RCM/GCM/CMIP5"
SRC_ROOT = os.path.join(BASE_CMIP5, "annual_nc")   # source: single-layer climatologies
DST_ROOT = os.path.join(BASE_CMIP5, "annual_nc")   # write alongside
SCENARIOS = ["historical", "RCP85"]

# ---------- helpers ----------
def ensure_dir(p): os.makedirs(p, exist_ok=True)

def find_file(model_dir: str, var: str, scenario: str):
    """Prefer *_AyrClim_*; fallback to *_yearmean_*; robust to RCP85/rcp85."""
    scen_tokens = [scenario, scenario.lower(), scenario.upper()]

    prefer = []
    for s in scen_tokens:
        prefer += [
            os.path.join(model_dir, f"{var}_AyrClim_*_{s}_*.nc"),
            os.path.join(model_dir, f"{var}_AyrClim_{s}_*.nc"),
            os.path.join(model_dir, f"{var}_*_{s}_AyrClim_*.nc"),
        ]
    fb = []
    for s in scen_tokens:
        fb += [
            os.path.join(model_dir, f"{var}_Amon_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"{var}_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"*_{var}_*_{s}_*yearmean*.nc"),
        ]
    fb += [os.path.join(model_dir, f"{var}_*yearmean*.nc")]

    hits = []
    for p in prefer: hits.extend(glob.glob(p))
    if hits: return sorted(set(hits))[0]
    hits = []
    for p in fb: hits.extend(glob.glob(p))
    return sorted(set(hits))[0] if hits else None

def pick_numeric_var(ds: xr.Dataset, prefer_tokens: list[str]) -> xr.DataArray:
    for tok in prefer_tokens:
        if tok in ds.data_vars and np.issubdtype(ds[tok].dtype, np.number):
            return ds[tok]
    for tok in prefer_tokens:
        for name in ds.data_vars:
            nm = name.lower()
            if tok in nm and np.issubdtype(ds[name].dtype, np.number):
                if any(bad in nm for bad in ["time", "bnds", "bounds", "average_t"]):
                    continue
                return ds[name]
    for name, da in ds.data_vars.items():
        nm = name.lower()
        if np.issubdtype(da.dtype, np.number) and not any(
            bad in nm for bad in ["time", "bnds", "bounds", "average_t"]
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
    return hfls.astype("float64")

# --- grid harmonization helpers ---
def _standardize_coords(da: xr.DataArray) -> xr.DataArray:
    """Rename common coord names, drop bounds, sort, and round coords slightly."""
    # 1) rename latitude/longitude -> lat/lon (no errors kw on DataArray.rename)
    for c in list(da.coords):
        cl = c.lower()
        if cl == "latitude":
            da = da.rename({c: "lat"})
        elif cl == "longitude":
            da = da.rename({c: "lon"})

    # 2) drop typical bounds coords if present (quietly skip if not there)
    for c in list(da.coords):
        if any(k in c.lower() for k in ("bnds", "bounds")):
            try:
                da = da.drop_vars(c)
            except Exception:
                # Some xarray versions store bounds as non-dim coords; ignore if drop fails
                pass

    # 3) normalize longitude domain if needed (wrap -180..180 to 0..360 to match many CMIP5 files)
    if "lon" in da.coords:
        lon = da["lon"].values
        try:
            lon = lon.astype(float)
            # If lots of negative longitudes, shift to [0, 360)
            if np.nanmean(lon) < 0:
                lon = (lon + 360.0) % 360.0
            da = da.assign_coords(lon=lon)
        except Exception:
            pass

    # 4) sort and lightly round to avoid tiny FP grid mismatches
    if "lat" in da.coords:
        da = da.sortby("lat")
        try:
            da = da.assign_coords(lat=np.round(da["lat"].astype(float).values, 6))
        except Exception:
            pass
    if "lon" in da.coords:
        da = da.sortby("lon")
        try:
            da = da.assign_coords(lon=np.round(da["lon"].astype(float).values, 6))
        except Exception:
            pass

    return da

def _align_or_regrid(vpd: xr.DataArray, hfls: xr.DataArray, regrid_method: str = "nearest"):
    """Try inner align; if spatial overlap empty, regrid hfls to vpd grid and retry."""
    vpd_s = _standardize_coords(vpd)
    hfls_s = _standardize_coords(hfls)

    v_al, h_al = xr.align(vpd_s, hfls_s, join="inner")

    def _any_zero_dim(a: xr.DataArray):
        return any(a.sizes[d] == 0 for d in a.sizes)

    if _any_zero_dim(v_al) or _any_zero_dim(h_al):
        # Log context before regridding
        def _shape_info(a):
            lat = a.coords.get("lat")
            lon = a.coords.get("lon")
            return {
                "shape": tuple(a.shape),
                "dims": a.dims,
                "lat_range": (float(lat.min()), float(lat.max())) if isinstance(lat, xr.DataArray) else None,
                "lon_range": (float(lon.min()), float(lon.max())) if isinstance(lon, xr.DataArray) else None,
            }
        L.warning("    grids mismatch → regridding hfls to vpd grid (%s).", regrid_method)
        try:
            hfls_rg = hfls_s.interp_like(vpd_s, method=regrid_method)
        except Exception as e:
            L.warning("    interp_like failed with '%s', retrying with 'nearest'.", e)
            hfls_rg = hfls_s.interp_like(vpd_s, method="nearest")
        v_al, h_al = xr.align(vpd_s, hfls_rg, join="inner")

    return v_al, h_al

def compute_geff(vpd: xr.DataArray, hfls: xr.DataArray) -> xr.DataArray:
    # Align on spatial dims; fallback to regrid if needed
    vpd2, hfls2 = _align_or_regrid(vpd, hfls, regrid_method="nearest")
    vpd2 = maybe_convert_vpd_to_kpa(vpd2)
    hfls2 = maybe_check_hfls_units(hfls2)

    rseff = (RHO_AIR * CP_AIR / GAMMA_KPA_PER_K) * (vpd2 / hfls2)
    geff = (1.0 / rseff) * 1000.0
    geff.name = "geff"
    geff.attrs.update({
        "long_name": "effective conductance",
        "formula": "geff = 1000 / (rho_air * cp_air / gamma * vpd / hfls)",
        "units": "mm s-1",
        "rho_air": RHO_AIR, "cp_air": CP_AIR, "gamma": GAMMA_KPA_PER_K,
        "expects": "vpd in kPa, hfls in W m-2; single-layer climatologies or aligned annual means",
    })
    return geff

def save_geff(geff: xr.DataArray, out_path: str):
    ensure_dir(os.path.dirname(out_path))
    comp = {geff.name: {"zlib": True, "complevel": 4, "_FillValue": None}}
    xr.Dataset({geff.name: geff}, coords=geff.coords).to_netcdf(out_path, encoding=comp)
    L.info(f"    ✓ saved: {out_path}")

def infer_year_window(path: str) -> str:
    base = os.path.basename(path)
    m = re.search(r'_(\d{4}-\d{4})(?:\.nc)?$', base)
    if m: return m.group(1)
    m = re.search(r'(\d{4})-(\d{4})', base)
    if m: return f"{m.group(1)}-{m.group(2)}"
    # single-layer files usually have no time axis, but keep fallback
    try:
        with xr.open_dataset(path, decode_times=True) as ds:
            if "time" in ds.coords and len(ds["time"]) > 0:
                years = np.array([getattr(t, "year", None) for t in ds.indexes["time"]])
                years = years[years != None]
                if years.size:
                    return f"{int(years.min())}-{int(years.max())}"
    except Exception:
        pass
    return "unknown"

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Compute geff from CMIP5 single-layer annual climatologies of vpd & hfls."
    )
    ap.add_argument("--src-root", default=SRC_ROOT, help=".../CMIP5/annual_nc")
    ap.add_argument("--dst-root", default=DST_ROOT, help=".../CMIP5/annual_nc (outputs here)")
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

            yrwin = infer_year_window(vpd_file) or "unknown"
            if yrwin == "unknown":
                yrwin = infer_year_window(hfls_file)

            out_dir = os.path.join(args.dst_root, scenario, model)
            out_path = os.path.join(out_dir, f"geff_AyrClim_{model}_{scenario}_{yrwin}.nc")
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
# python cmip5_geff_from_annual_clim.py
# python cmip5_geff_from_annual_clim.py --dry-run
# python cmip5_geff_from_annual_clim.py --overwrite
# python cmip5_geff_from_annual_clim.py --scenarios RCP85
