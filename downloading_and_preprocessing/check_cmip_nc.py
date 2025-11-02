#!/usr/bin/env python3
import argparse, os, glob, logging
import numpy as np
import xarray as xr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
L = logging.getLogger("check_cmip_nc")

def pick_numeric_var(ds):
    # Prefer obvious names, then any numeric non-bounds
    prefs = ["geff","vpd","hfls","surface_upward_latent_heat_flux"]
    for p in prefs:
        if p in ds.data_vars and np.issubdtype(ds[p].dtype, np.number):
            return p
    for k, da in ds.data_vars.items():
        kl = k.lower()
        if any(b in kl for b in ["bnds","bounds","time_bnds","average_t","climatology_bounds"]):
            continue
        if np.issubdtype(da.dtype, np.number):
            return k
    return None

def stats_for(path):
    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            vname = pick_numeric_var(ds)
            if vname is None:
                return dict(status="no_numeric_var")
            da = ds[vname].squeeze(drop=True)
            # check dims non-empty
            for d, n in da.sizes.items():
                if n == 0:
                    return dict(status="zero_dim", var=vname, dims=dict(da.sizes))
            vals = da.astype("float64").values
            total = vals.size
            finite = np.isfinite(vals).sum()
            frac = (finite / total) if total else 0.0
            mn = float(np.nanmin(vals)) if finite else np.nan
            md = float(np.nanmedian(vals)) if finite else np.nan
            mx = float(np.nanmax(vals)) if finite else np.nan
            units = str(da.attrs.get("units",""))
            return dict(status="ok" if frac>0 else "all_nan",
                        var=vname, units=units, total=int(total),
                        finite=int(finite), frac=float(frac),
                        min=mn, med=md, max=mx,
                        dims={k:int(v) for k,v in da.sizes.items()})
    except Exception as e:
        return dict(status="error", error=str(e))

def main():
    ap = argparse.ArgumentParser(description="Scan CMIP files and flag empty/all-NaN datasets.")
    ap.add_argument("--root", required=True, help="e.g. /mnt/data/Other/GCM-RCM/GCM/CMIP6/annual_nc")
    ap.add_argument("--scenarios", nargs="+", default=["historical","ssp585"])
    ap.add_argument("--glob", default="**/*.nc", help="Glob pattern under scenario dirs")
    args = ap.parse_args()

    for scen in args.scenarios:
        base = os.path.join(args.root, scen)
        if not os.path.isdir(base):
            L.info("Skip missing scenario dir: %s", base)
            continue
        L.info("=== %s ===", scen)
        files = sorted(glob.glob(os.path.join(base, args.glob), recursive=True))
        bad = []
        for f in files:
            r = stats_for(f)
            if r["status"] != "ok":
                bad.append((f, r))
        if not bad:
            L.info("All OK in %s (%d files checked).", scen, len(files))
        else:
            L.info("Problems in %s (%d issues):", scen, len(bad))
            for f, r in bad:
                L.info(" - %s : %s", f, r)

if __name__ == "__main__":
    main()
