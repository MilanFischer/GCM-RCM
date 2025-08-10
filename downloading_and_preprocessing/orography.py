import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/mnt/data/Other/GCM-RCM/orography/bcc-csm1-1")

# Constants
R = 287.0     # J/kg/K, specific gas constant for dry air
g = 9.81      # m/s^2, gravity

# Load datasets
ds_ps = xr.open_dataset('ps_Amon_bcc-csm1-1_historical_r1i1p1_185001-201212.nc')
ds_psl = xr.open_dataset('psl_Amon_bcc-csm1-1_historical_r1i1p1_185001-201212.nc')
ds_tas = xr.open_dataset('tas_Amon_bcc-csm1-1_historical_r1i1p1_185001-201212.nc')

# Extract variables
ps = ds_ps['ps']     # Surface pressure [Pa]
psl = ds_psl['psl']  # Sea level pressure [Pa]
tas = ds_tas['tas']  # Near-surface temperature [K]

# Align time and spatial dimensions
ps, psl, tas = xr.align(ps, psl, tas)

# Compute elevation from hypsometric equation
z_est = (R * tas / g) * np.log(psl / ps)

# Compute time-mean and std deviation
z_mean = z_est.mean(dim='time')
z_std = z_est.std(dim='time')

# Plot estimated mean elevation
plt.figure(figsize=(10, 4))
z_mean.plot(cmap='terrain')
plt.title('Estimated Orography (Mean from ps, psl, tas)')
plt.show()

plt.savefig('estimated_orography_mean.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory

# Plot variability
plt.figure(figsize=(10, 4))
z_std.plot(cmap='magma')
plt.title('Estimated Orography Variability (Standard Deviation)')
plt.show()

plt.savefig('estimated_orography_sd.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory


lonmin, lonmax = 9, 21
latmin, latmax = 47, 54

# Subset the region
z_mean_region = z_mean.sel(lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))

# Compute mean elevation over that region
regional_mean_elevation = z_mean_region.mean().item()  # convert to scalar float

print(f"Mean estimated elevation in region ({latmin}-{latmax}N, {lonmin}-{lonmax}E): {regional_mean_elevation:.2f} meters")


# Assign variable name and attributes for clarity
z_mean.name = 'orog_estimated'
z_mean.attrs['units'] = 'meters'
z_mean.attrs['long_name'] = 'Estimated Orography from ps, psl, tas'

# Save to NetCDF
z_mean.to_netcdf('estimated_orography_mean.nc')
print("Saved z_mean to 'estimated_orography_mean.nc'")