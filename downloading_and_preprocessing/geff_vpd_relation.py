#!/usr/bin/env python3
"""
Per-pixel regression of geff vs 1/sqrt(VPD) from annual climatologies,
plus period-specific multi-model means.

Outputs a NetCDF with layers:
  - a, b, r2  (from geff = a + b * 1/sqrt(VPD))
  - mean_geff_1981_2005, mean_vpd_1981_2005      (if historical present)
  - mean_geff_2076_2100, mean_vpd_2076_2100      (if RCP85/ssp585 present)

Example:
  python geff_vpd_relation.py \
    --roots /mnt/data/Other/GCM-RCM/GCM/CMIP5 /mnt/data/Other/GCM-RCM/GCM/CMIP6 \
    --scenarios historical RCP85 ssp585 \
    --regrid nearest \
    --min-finite-frac 0.1 \
    --vpd-floor 0.05 \
    --out /mnt/data/Other/GCM-RCM/outputs/geff_vpd_regression_CMIP5_CMIP6.nc \
    --overwrite
"""
import argparse
import glob
import logging
import os
from typing import List, Optional, Tuple, Dict

import numpy as np
import xarray as xr

# ----------------------- Logging -----------------------
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("stack_geff_vpd")

# ----------------------- Helpers -----------------------
def ensure_dir(p: str): os.makedirs(p, exist_ok=True)

def _standardize_coords(da: xr.DataArray) -> xr.DataArray:
    # rename common variants
    for c in list(da.coords):
        cl = c.lower()
        if cl == "latitude":
            da = da.rename({c: "lat"})
        elif cl == "longitude":
            da = da.rename({c: "lon"})
    # drop bounds & stray coords
    for c in list(da.coords):
        cl = c.lower()
        if any(k in cl for k in ("bnds", "bounds")):
            try: da = da.drop_vars(c)
            except Exception: pass
        if c not in ("lat", "lon") and c in da.coords and da[c].ndim == 0:
            try: da = da.drop_vars(c)
            except Exception: pass
    # normalize lon to 0..360 and sort
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
    for p in pats: hits.extend(glob.glob(p))
    hits = sorted(set(hits))
    return hits[0] if hits else None

def _find_geff_file(model_dir: str, scenario: str) -> Optional[str]:
    scen = [scenario, scenario.lower(), scenario.upper()]
    pats = []
    for s in scen:
        pats += [
            os.path.join(model_dir, f"geff_AyrClim_*_{s}_*.nc"),
            os.path.join(model_dir, f"geff_*_{s}_AyrClim_*.nc"),
        ]
    return _find_single(pats)

def _find_vpd_clim_or_series(model_dir: str, scenario: str) -> Tuple[Optional[str], Optional[str]]:
    scen = [scenario, scenario.lower(), scenario.upper()]
    pats_clim = []
    for s in scen:
        pats_clim += [
            os.path.join(model_dir, f"vpd_AyrClim_*_{s}_*.nc"),
            os.path.join(model_dir, f"vpd_*_{s}_AyrClim_*.nc"),
        ]
    vpd_clim = _find_single(pats_clim)
    if vpd_clim:
        return vpd_clim, None
    pats_series = []
    for s in scen:
        pats_series += [
            os.path.join(model_dir, f"vpd_Amon_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"vpd_*_{s}_yearmean_*.nc"),
            os.path.join(model_dir, f"*vpd*_{s}_*yearmean*.nc"),
        ]
    vpd_series = _find_single(pats_series)
    return None, vpd_series

def _to_single_layer_vpd(path: str) -> xr.DataArray:
    """If 'time' exists, mean over time; else return numeric var."""
    with xr.open_dataset(path) as ds:
        vname = None
        if "vpd" in ds and np.issubdtype(ds["vpd"].dtype, np.number):
            vname = "vpd"
        if vname is None:
            for name, da in ds.data_vars.items():
                if np.issubdtype(da.dtype, np.number) and "vpd" in name.lower():
                    vname = name; break
        if vname is None:
            for name, da in ds.data_vars.items():
                if np.issubdtype(da.dtype, np.number) and not any(k in name.lower() for k in ["time","bnds","bounds"]):
                    vname = name; break
        if vname is None:
            raise ValueError(f"No numeric VPD variable in {path}")
        da = ds[vname]
        if "time" in da.dims and da.sizes.get("time", 0) > 0:
            da = da.mean("time", skipna=True)
        return da.squeeze(drop=True)

def _vpd_to_kpa(vpd_da: xr.DataArray) -> xr.DataArray:
    units = str(vpd_da.attrs.get("units","")).lower()
    v = vpd_da.astype("float64")
    if units in ("pa","pascal"):
        v = v / 1000.0
    v.attrs["units"] = "kPa"
    return v

def _inv_sqrt_vpd(vpd_kpa: xr.DataArray, vpd_floor: Optional[float]) -> xr.DataArray:
    v = vpd_kpa.astype("float64")
    if vpd_floor is not None:
        v = xr.apply_ufunc(np.maximum, v, vpd_floor)
    inv = 1.0 / xr.apply_ufunc(np.sqrt, v)
    inv.name = "inv_sqrt_vpd"
    inv.attrs["units"] = "kPa^-0.5"
    return inv

def _finite_fraction(da: xr.DataArray) -> float:
    arr = np.asarray(da, dtype="float64")
    tot = arr.size
    if tot == 0: return 0.0
    return np.isfinite(arr).sum() / float(tot)

def _nan_safe_regression(x: np.ndarray, y: np.ndarray, axis=0):
    xbar = np.nanmean(x, axis=axis, keepdims=True)
    ybar = np.nanmean(y, axis=axis, keepdims=True)
    xm = x - xbar
    ym = y - ybar
    cov = np.nanmean(xm * ym, axis=axis)
    varx = np.nanmean(xm * xm, axis=axis)
    vary = np.nanmean(ym * ym, axis=axis)
    b = cov / varx
    a = (ybar.squeeze(axis=axis)) - b * (xbar.squeeze(axis=axis))
    r2 = (cov * cov) / (varx * vary)
    n_valid = np.sum(np.isfinite(x) & np.isfinite(y), axis=axis)
    bad = (n_valid < 2) | (varx <= 0) | (vary <= 0)
    b = np.where(bad, np.nan, b)
    a = np.where(bad, np.nan, a)
    r2 = np.where(bad, np.nan, r2)
    return a, b, r2

def _period_label(scen: str) -> Optional[str]:
    s = scen.lower()
    if s == "historical": return "1981_2005"
    if s in ("rcp85", "ssp585"): return "2076_2100"
    return None

def _scenario_dirs(root_annual_nc: str, scen: str) -> List[str]:
    """
    Return possible scenario directories under a given root annual_nc.
    - CMIP5: root/annual_nc/<scenario>/
    - CMIP6: root/annual_nc/ScenarioMIP/<scenario>/
    """
    cand = [
        os.path.join(root_annual_nc, scen),                     # CMIP5
        os.path.join(root_annual_nc, "ScenarioMIP", scen),      # CMIP6
    ]
    return [d for d in cand if os.path.isdir(d)]

# ----------------------- Main -----------------------
def main():
    ap = argparse.ArgumentParser(
        description="Per-pixel regression of geff vs 1/sqrt(VPD) + period means."
    )
    ap.add_argument("--roots", nargs="+", required=True,
                    help="Base roots; each must contain annual_nc/<scenario>/<model>/ (CMIP5) or annual_nc/ScenarioMIP/<scenario>/<model>/ (CMIP6)")
    ap.add_argument("--scenarios", nargs="+", default=["historical","RCP85","ssp585"])
    ap.add_argument("--regrid", choices=["linear","nearest"], default="linear")
    ap.add_argument("--min-finite-frac", type=float, default=None,
                    help="Skip model/period pairs if finite fraction of geff or VPD below this number (e.g., 0.1).")
    ap.add_argument("--vpd-floor", type=float, default=None,
                    help="Minimum VPD (kPa) before computing 1/sqrt(VPD); avoids huge leverage at tiny VPD.")
    ap.add_argument("--out", required=True, help="Output NetCDF path")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.overwrite:
        L.info("Output exists; use --overwrite: %s", args.out)
        return

    pairs = []   # (model, scenario, geff_rg, inv_rg, vpd_kpa_rg, period_tag)
    ref_grid = None
    canon_lat = None
    canon_lon = None

    # For period-specific means & accounting
    mean_acc: Dict[str, Dict[str, List[np.ndarray]]] = {
        "1981_2005": {"geff": [], "vpd": []},
        "2076_2100": {"geff": [], "vpd": []},
    }

    per_root_counts = {}
    per_scen_counts = {s: 0 for s in args.scenarios}
    per_scen_models  = {s: [] for s in args.scenarios}
    skipped_models = {s: {"geff": [], "vpd": [], "low_finite": []} for s in args.scenarios}

    for root in args.roots:
        annual_nc = os.path.join(root, "annual_nc")
        if not os.path.isdir(annual_nc):
            L.warning("Missing annual_nc under %s (skip).", root)
            continue
        before = len(pairs)
        for scen in args.scenarios:
            scen_dirs = _scenario_dirs(annual_nc, scen)
            if not scen_dirs:
                continue
            for scen_dir in scen_dirs:
                models = sorted(d for d in os.listdir(scen_dir) if os.path.isdir(os.path.join(scen_dir, d)))
                for model in models:
                    mdir = os.path.join(scen_dir, model)
                    geff_path = _find_geff_file(mdir, scen)
                    if not geff_path:
                        skipped_models[scen]["geff"].append(model)
                        continue
                    vpd_clim, vpd_series = _find_vpd_clim_or_series(mdir, scen)
                    if not (vpd_clim or vpd_series):
                        skipped_models[scen]["vpd"].append(model)
                        continue

                    # read geff
                    with xr.open_dataset(geff_path) as ds_g:
                        gvar = "geff" if "geff" in ds_g.data_vars else next(
                            (n for n, da in ds_g.data_vars.items()
                             if np.issubdtype(da.dtype, np.number) and "geff" in n.lower()), None)
                        if gvar is None:
                            gvar = next(n for n, da in ds_g.data_vars.items()
                                        if np.issubdtype(da.dtype, np.number))
                        geff = _standardize_coords(ds_g[gvar].squeeze(drop=True).astype("float64"))

                    # read vpd (clim or series->mean)
                    if vpd_clim:
                        with xr.open_dataset(vpd_clim) as ds_v:
                            vname = "vpd" if "vpd" in ds_v.data_vars else next(
                                (n for n, da in ds_v.data_vars.items()
                                 if np.issubdtype(da.dtype, np.number) and "vpd" in n.lower()), None)
                            if vname is None:
                                vname = next(n for n, da in ds_v.data_vars.items()
                                             if np.issubdtype(da.dtype, np.number))
                            vpd = _standardize_coords(ds_v[vname].squeeze(drop=True).astype("float64"))
                    else:
                        vpd = _standardize_coords(_to_single_layer_vpd(vpd_series))

                    # data-quality filters
                    ff_geff = _finite_fraction(geff)
                    ff_vpd  = _finite_fraction(vpd)
                    if args.min_finite_frac is not None and (ff_geff < args.min_finite_frac or ff_vpd < args.min_finite_frac):
                        skipped_models[scen]["low_finite"].append(f"{model}(g={ff_geff:.2f},v={ff_vpd:.2f})")
                        continue

                    vpd_kpa = _vpd_to_kpa(vpd)
                    inv = _inv_sqrt_vpd(vpd_kpa, args.vpd_floor)

                    if ref_grid is None:
                        ref_grid = geff
                        L.info("Reference grid set by %s/%s  (shape: %s)", model, scen, tuple(ref_grid.shape))
                        # Canonical coordinate objects (exact identity for labels)
                        canon_lat = ref_grid["lat"]
                        canon_lon = ref_grid["lon"]

                    # Regrid (if needed)
                    geff_rg    = geff    if geff.shape    == ref_grid.shape else _regrid_like(geff,    ref_grid, method=args.regrid)
                    vpd_kpa_rg = vpd_kpa if vpd_kpa.shape == ref_grid.shape else _regrid_like(vpd_kpa, ref_grid, method=args.regrid)
                    inv_rg     = inv     if inv.shape     == ref_grid.shape else _regrid_like(inv,     ref_grid, method=args.regrid)

                    # Enforce exact coordinate equality with canonical coords
                    geff_rg    = geff_rg.reindex_like(ref_grid)   .assign_coords(lat=canon_lat, lon=canon_lon)
                    vpd_kpa_rg = vpd_kpa_rg.reindex_like(ref_grid).assign_coords(lat=canon_lat, lon=canon_lon)
                    inv_rg     = inv_rg.reindex_like(ref_grid)    .assign_coords(lat=canon_lat, lon=canon_lon)

                    # drop stray coords to make concat robust
                    geff_rg    = geff_rg.reset_coords(drop=True)
                    vpd_kpa_rg = vpd_kpa_rg.reset_coords(drop=True)
                    inv_rg     = inv_rg.reset_coords(drop=True)

                    # add
                    tag = _period_label(scen)
                    pairs.append((model, scen, geff_rg, inv_rg, vpd_kpa_rg, tag))
                    per_scen_counts[scen] += 1
                    per_scen_models[scen].append(model)
                    if tag in mean_acc:
                        mean_acc[tag]["geff"].append(geff_rg.values)
                        mean_acc[tag]["vpd"].append(vpd_kpa_rg.values)

        per_root_counts[root] = len(pairs) - before

    total_pairs = len(pairs)
    if not pairs:
        L.error("No (geff, vpd) pairs found — nothing to do.")
        return

    # Summary to console
    L.info("Found %d model-scenario pairs.", total_pairs)
    for root, cnt in per_root_counts.items():
        L.info("  • From root %s: %d pairs", root, cnt)
    for scen in args.scenarios:
        if per_scen_counts[scen] > 0:
            uniq_models = ", ".join(sorted(set(per_scen_models[scen])))
            L.info("  • %s: %d pairs; models: %s", scen, per_scen_counts[scen], uniq_models)
    for scen in args.scenarios:
        miss_g = ", ".join(sorted(set(skipped_models[scen]["geff"]))) or "-"
        miss_v = ", ".join(sorted(set(skipped_models[scen]["vpd"]))) or "-"
        lowf   = ", ".join(skipped_models[scen]["low_finite"]) or "-"
        L.info("Missing inputs in %s -> geff: %s ; vpd: %s ; low_finite: %s", scen, miss_g, miss_v, lowf)

    # ---------- Build stacks & regress ----------
    G  = xr.concat([p[2] for p in pairs], dim="sample", coords="minimal").assign_coords(sample=np.arange(total_pairs))
    X  = xr.concat([p[3] for p in pairs], dim="sample", coords="minimal").assign_coords(sample=np.arange(total_pairs))

    # Use canonical coords (fixed  lat/lon labels and sizes)
    lat = canon_lat
    lon = canon_lon
    nlat = lat.size
    nlon = lon.size

    # Regression with numpy
    x = X.transpose("sample","lat","lon").values  # (S, nlat, nlon)
    y = G.transpose("sample","lat","lon").values  # (S, nlat, nlon)
    a, b, r2 = _nan_safe_regression(x, y, axis=0) # (nlat, nlon)

    # Safety: auto-fix if someone came out transposed
    if a.shape == (nlon, nlat):
        L.warning("Detected transposed regression outputs; auto-rotating to (lat, lon).")
        a = a.T; b = b.T; r2 = r2.T

    if a.shape != (nlat, nlon):
        raise ValueError(f"Regression outputs shape {a.shape} does not match stack grid {(nlat, nlon)}")

    # Period-specific means (already regridded → shapes match nlat×nlon)
    data_vars = dict(
        a=(("lat","lon"), a, {"long_name":"intercept", "units":"mm s-1"}),
        b=(("lat","lon"), b, {"long_name":"slope", "units":"mm s-1 per kPa^-0.5"}),
        r2=(("lat","lon"), r2, {"long_name":"coefficient of determination", "units":""}),
    )
    if mean_acc["1981_2005"]["geff"]:
        mg = np.nanmean(np.stack(mean_acc["1981_2005"]["geff"], axis=0), axis=0)
        mv = np.nanmean(np.stack(mean_acc["1981_2005"]["vpd"],  axis=0), axis=0)
        data_vars["mean_geff_1981_2005"] = (("lat","lon"), mg, {"long_name":"multi-model mean geff (1981–2005)", "units":"mm s-1"})
        data_vars["mean_vpd_1981_2005"]  = (("lat","lon"), mv, {"long_name":"multi-model mean VPD (1981–2005)",  "units":"kPa"})
    if mean_acc["2076_2100"]["geff"]:
        mg = np.nanmean(np.stack(mean_acc["2076_2100"]["geff"], axis=0), axis=0)
        mv = np.nanmean(np.stack(mean_acc["2076_2100"]["vpd"],  axis=0), axis=0)
        data_vars["mean_geff_2076_2100"] = (("lat","lon"), mg, {"long_name":"multi-model mean geff (2076–2100)", "units":"mm s-1"})
        data_vars["mean_vpd_2076_2100"]  = (("lat","lon"), mv, {"long_name":"multi-model mean VPD (2076–2100)",  "units":"kPa"})

    attrs = dict(
        title="Per-pixel OLS: geff = a + b * 1/sqrt(VPD); plus period-specific multi-model means",
        note="Regridded to canonical lat/lon via interp_like + reindex_like; output coords taken from canonical reference grid.",
        samples=total_pairs,
        per_root_counts=str(per_root_counts),
        per_scenario_counts=str(per_scen_counts),
        skipped=str(skipped_models),
        regrid_method=args.regrid,
        min_finite_frac=args.min_finite_frac if args.min_finite_frac is not None else "none",
        vpd_floor_kpa=args.vpd_floor if args.vpd_floor is not None else "none",
    )

    ds_out = xr.Dataset(data_vars=data_vars, coords=dict(lat=lat, lon=lon), attrs=attrs)

    ensure_dir(os.path.dirname(args.out))
    comp = {v: {"zlib": True, "complevel": 4} for v in ds_out.data_vars}
    ds_out.to_netcdf(args.out, encoding=comp)
    L.info("✓ Saved: %s (vars: %s)", args.out, ", ".join(ds_out.data_vars))


if __name__ == "__main__":
    main()
