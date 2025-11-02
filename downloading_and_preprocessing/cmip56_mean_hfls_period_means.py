#!/usr/bin/env python3
"""
Compute multi-model mean HFLS (surface upward latent heat flux; W m-2)
for two standard periods across CMIP5 and CMIP6:
  - historical  → 1981–2005
  - RCP85/ssp585 → 2076–2100

The script scans one or more CMIP roots that contain `annual_nc/.../<scenario>/<model>/`.
It finds either single-layer climatologies (preferred: *_AyrClim_*) or time series
of annual means (fallback: *_yearmean_*), harmonizes grids, and builds period-
specific multi-model means on a common reference grid.

Outputs a NetCDF containing variables:
  - mean_hfls_1981_2005 (if historical present)
  - mean_hfls_2076_2100 (if RCP85 and/or ssp585 present)

Example:
  python cmip56_mean_hfls_period_means.py \
    --roots /mnt/data/Other/GCM-RCM/GCM/CMIP5 /mnt/data/Other/GCM-RCM/GCM/CMIP6 \
    --scenarios historical RCP85 ssp585 \
    --regrid linear \
    --min-finite-frac 0.1 \
    --out /mnt/data/Other/GCM-RCM/outputs/mean_hfls_CMIP5_CMIP6.nc \
    --overwrite
"""
import argparse
import glob
import logging
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import xarray as xr

# ----------------------- Logging -----------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("mean_hfls_periods")

# ----------------------- Helpers -----------------------
def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def _standardize_coords(da: xr.DataArray) -> xr.DataArray:
    """Rename common coord names, drop bounds, coerce lon to 0..360 if needed,
    lightly round coords, sort, and ensure only lat/lon remain as dims."""
    # rename latitude/longitude → lat/lon
    for c in list(da.coords):
        cl = c.lower()
        if cl == "latitude":
            da = da.rename({c: "lat"})
        elif cl == "longitude":
            da = da.rename({c: "lon"})
    # drop bounds & stray scalar coords
    for c in list(da.coords):
        cl = c.lower()
        if any(k in cl for k in ("bnds", "bounds")):
            try:
                da = da.drop_vars(c)
            except Exception:
                pass
        if c not in ("lat", "lon") and c in da.coords and getattr(da[c], "ndim", 1) == 0:
            try:
                da = da.drop_vars(c)
            except Exception:
                pass
    # normalize lon to [0, 360)
    if "lon" in da.coords:
        try:
            lon = da["lon"].values.astype(float)
            if np.nanmean(lon) < 0:
                lon = (lon + 360.0) % 360.0
            da = da.assign_coords(lon=np.round(lon, 6))
        except Exception:
            pass
        da = da.sortby("lon")
    if "lat" in da.coords:
        try:
            da = da.assign_coords(lat=np.round(da["lat"].astype(float).values, 6))
        except Exception:
            pass
        da = da.sortby("lat")
    # drop non lat/lon dims if any slipped in
    for d in list(da.dims):
        if d not in ("lat", "lon"):
            try:
                da = da.isel({d: 0}, drop=True)
            except Exception:
                pass
    return da


def _regrid_like(src: xr.DataArray, tgt: xr.DataArray, method: str) -> xr.DataArray:
    src_s = _standardize_coords(src)
    tgt_s = _standardize_coords(tgt)
    try:
        return src_s.interp_like(tgt_s, method=method)
    except Exception as e:
        L.warning("interp_like '%s' failed (%s), falling back to nearest.", method, e)
        return src_s.interp_like(tgt_s, method="nearest")


def _find_single(pats: List[str]) -> Optional[str]:
    hits = []
    for p in pats:
        hits.extend(glob.glob(p))
    hits = sorted(set(hits))
    return hits[0] if hits else None


def _find_hfls_clim_or_series(model_dir: str, scenario: str) -> Tuple[Optional[str], Optional[str]]:
    scen = [scenario, scenario.lower(), scenario.upper()]
    pats_clim = []
    for s in scen:
        pats_clim += [
            os.path.join(model_dir, f"hfls_AyrClim_*_{s}_*.nc"),
            os.path.join(model_dir, f"hfls_*_{s}_AyrClim_*.nc"),
        ]
    hfls_clim = _find_single(pats_clim)
    if hfls_clim:
        return hfls_clim, None

    pats_series = []
    for s in scen:
        pats_series += [
            os.path.join(model_dir, f"hfls_Amon_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"hfls_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"*hfls*_{s}_*yearmean*.nc"),
        ]
    hfls_series = _find_single(pats_series)
    return None, hfls_series


def _to_single_layer_hfls(path: str) -> xr.DataArray:
    """If 'time' exists, mean over time; else return a suitable numeric var."""
    with xr.open_dataset(path) as ds:
        vname = None
        # try canonical
        if "hfls" in ds and np.issubdtype(ds["hfls"].dtype, np.number):
            vname = "hfls"
        # try any variable with hfls-like name
        if vname is None:
            for name, da in ds.data_vars.items():
                if np.issubdtype(da.dtype, np.number) and "hfls" in name.lower():
                    vname = name
                    break
        # last resort: first numeric non-time/bounds
        if vname is None:
            for name, da in ds.data_vars.items():
                nm = name.lower()
                if np.issubdtype(da.dtype, np.number) and not any(k in nm for k in ["time", "bnds", "bounds"]):
                    vname = name
                    break
        if vname is None:
            raise ValueError(f"No numeric HFLS variable in {path}")
        da = ds[vname]
        if "time" in da.dims and da.sizes.get("time", 0) > 0:
            da = da.mean("time", skipna=True)
        return da.squeeze(drop=True).astype("float64")


def _period_label(scen: str) -> Optional[str]:
    s = scen.lower()
    if s == "historical":
        return "1981_2005"
    if s in ("rcp85", "ssp585"):
        return "2076_2100"
    return None


def _scenario_dirs(root_annual_nc: str, scen: str) -> List[str]:
    """Return possible scenario directories under a given root annual_nc.
    - CMIP5: root/annual_nc/<scenario>/
    - CMIP6: root/annual_nc/ScenarioMIP/<scenario>/
    """
    cand = [
        os.path.join(root_annual_nc, scen),  # CMIP5
        os.path.join(root_annual_nc, "ScenarioMIP", scen),  # CMIP6
    ]
    return [d for d in cand if os.path.isdir(d)]


# ----------------------- Main -----------------------

def main():
    ap = argparse.ArgumentParser(
        description="Per-pixel multi-model mean HFLS (W m-2) for 1981–2005 and 2076–2100 across CMIP5 & CMIP6."
    )
    ap.add_argument(
        "--roots",
        nargs="+",
        required=True,
        help=(
            "Base roots; each must contain annual_nc/<scenario>/<model>/ (CMIP5) "
            "or annual_nc/ScenarioMIP/<scenario>/<model>/ (CMIP6)"
        ),
    )
    ap.add_argument("--scenarios", nargs="+", default=["historical", "RCP85", "ssp585"])  # what to scan
    ap.add_argument("--regrid", choices=["linear", "nearest"], default="linear")
    ap.add_argument(
        "--min-finite-frac",
        type=float,
        default=None,
        help="Skip a model if finite fraction below this (e.g., 0.1).",
    )
    ap.add_argument("--out", required=True, help="Output NetCDF path")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.overwrite:
        L.info("Output exists; use --overwrite: %s", args.out)
        return

    ref_grid = None
    canon_lat = None
    canon_lon = None

    # Accumulators for period-specific means (list of np arrays already on ref grid)
    mean_acc: Dict[str, List[np.ndarray]] = {
        "1981_2005": [],
        "2076_2100": [],
    }

    per_root_counts: Dict[str, int] = {}
    per_scen_counts: Dict[str, int] = {s: 0 for s in args.scenarios}
    per_scen_models: Dict[str, List[str]] = {s: [] for s in args.scenarios}
    skipped: Dict[str, Dict[str, List[str]]] = {s: {"hfls": [], "low_finite": []} for s in args.scenarios}

    for root in args.roots:
        annual_nc = os.path.join(root, "annual_nc")
        if not os.path.isdir(annual_nc):
            L.warning("Missing annual_nc under %s (skip).", root)
            continue
        before = sum(len(v) for v in mean_acc.values())

        for scen in args.scenarios:
            scen_dirs = _scenario_dirs(annual_nc, scen)
            if not scen_dirs:
                continue
            for scen_dir in scen_dirs:
                models = sorted(d for d in os.listdir(scen_dir) if os.path.isdir(os.path.join(scen_dir, d)))
                for model in models:
                    mdir = os.path.join(scen_dir, model)
                    hfls_clim, hfls_series = _find_hfls_clim_or_series(mdir, scen)
                    if not (hfls_clim or hfls_series):
                        skipped[scen]["hfls"].append(model)
                        continue

                    # read HFLS (clim or series→mean)
                    if hfls_clim:
                        with xr.open_dataset(hfls_clim) as ds:
                            # pick a numeric hfls-like variable
                            vname = "hfls" if "hfls" in ds.data_vars else None
                            if vname is None:
                                for name, da in ds.data_vars.items():
                                    if np.issubdtype(da.dtype, np.number) and "hfls" in name.lower():
                                        vname = name
                                        break
                            if vname is None:
                                # fallback: first numeric non-bounds
                                for name, da in ds.data_vars.items():
                                    nm = name.lower()
                                    if np.issubdtype(da.dtype, np.number) and not any(k in nm for k in ["time", "bnds", "bounds"]):
                                        vname = name
                                        break
                            if vname is None:
                                L.info("    no numeric hfls in %s", mdir)
                                skipped[scen]["hfls"].append(model)
                                continue
                            hfls = _standardize_coords(ds[vname].squeeze(drop=True).astype("float64"))
                    else:
                        hfls = _standardize_coords(_to_single_layer_hfls(hfls_series))

                    # data-quality filter
                    arr = np.asarray(hfls, dtype="float64")
                    ff = np.isfinite(arr).sum() / float(arr.size) if arr.size else 0.0
                    if args.min_finite_frac is not None and ff < args.min_finite_frac:
                        skipped[scen]["low_finite"].append(f"{model}(f={ff:.2f})")
                        continue

                    if ref_grid is None:
                        ref_grid = hfls
                        L.info("Reference grid set by %s/%s  (shape: %s)", model, scen, tuple(ref_grid.shape))
                        canon_lat = ref_grid["lat"]
                        canon_lon = ref_grid["lon"]

                    # Regrid if needed and enforce identical coords
                    hfls_rg = hfls if hfls.shape == ref_grid.shape else _regrid_like(hfls, ref_grid, method=args.regrid)
                    hfls_rg = hfls_rg.reindex_like(ref_grid).assign_coords(lat=canon_lat, lon=canon_lon).reset_coords(drop=True)

                    tag = _period_label(scen)
                    if tag is not None:
                        mean_acc[tag].append(hfls_rg.values)
                        per_scen_counts[scen] += 1
                        per_scen_models[scen].append(model)

        per_root_counts[root] = sum(len(v) for v in mean_acc.values()) - before

    if not (mean_acc["1981_2005"] or mean_acc["2076_2100"]):
        L.error("No HFLS inputs found — nothing to do.")
        return

    # Build output dataset
    data_vars = {}
    if mean_acc["1981_2005"]:
        mh = np.nanmean(np.stack(mean_acc["1981_2005"], axis=0), axis=0)
        data_vars["mean_hfls_1981_2005"] = (("lat", "lon"), mh, {"long_name": "multi-model mean HFLS (1981–2005)", "units": "W m-2"})
    if mean_acc["2076_2100"]:
        mh = np.nanmean(np.stack(mean_acc["2076_2100"], axis=0), axis=0)
        data_vars["mean_hfls_2076_2100"] = (("lat", "lon"), mh, {"long_name": "multi-model mean HFLS (2076–2100)", "units": "W m-2"})

    attrs = dict(
        title="Multi-model mean HFLS for 1981–2005 and 2076–2100 (CMIP5 & CMIP6)",
        note=(
            "Inputs are single-layer climatologies or annual means converted to a single layer; "
            "regridded to a canonical lat/lon via interp_like + reindex_like; output coords from reference grid."
        ),
        per_root_counts=str(per_root_counts),
        per_scenario_counts=str(per_scen_counts),
        per_scenario_models=str({k: sorted(set(v)) for k, v in per_scen_models.items()}),
        skipped=str(skipped),
        regrid_method=args.regrid,
        min_finite_frac=args.min_finite_frac if args.min_finite_frac is not None else "none",
    )

    ds_out = xr.Dataset(data_vars=data_vars, coords=dict(lat=canon_lat, lon=canon_lon), attrs=attrs)

    ensure_dir(os.path.dirname(args.out))
    comp = {v: {"zlib": True, "complevel": 4} for v in ds_out.data_vars}
    ds_out.to_netcdf(args.out, encoding=comp)
    L.info("\u2713 Saved: %s (vars: %s)", args.out, ", ".join(ds_out.data_vars))


if __name__ == "__main__":
    main()

#python cmip56_mean_hfls_period_means.py \
#  --roots /mnt/data/Other/GCM-RCM/GCM/CMIP5 /mnt/data/Other/GCM-RCM/GCM/CMIP6 \
#  --scenarios historical RCP85 ssp585 \
#  --regrid linear \
#  --min-finite-frac 0.1 \
#  --out /mnt/data/Other/GCM-RCM/outputs/mean_hfls_CMIP5_CMIP6.nc \
#  --overwrite