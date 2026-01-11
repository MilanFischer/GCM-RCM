rm(list = ls())

library(terra)
library(rnaturalearth)
library(sf)

VPD_Tgt5C_hist <- rast("ERA5-land_VPD_mean_Tgt5C_1981-2005.nc")
Ta <- mean(rast("ERA5-land_T2m_daily_1981-2005_ymon.nc"))

hfls_ensemble <- rotate(rast("./mean_hfls_CMIP5_CMIP6.nc"))

plot(hfls_ensemble)
delta_hfls_over_hfls <- (hfls_ensemble$mean_hfls_2076_2100 - hfls_ensemble$mean_hfls_1981_2005) / hfls_ensemble$mean_hfls_1981_2005
delta_hfls_over_hfls <- ifel(delta_hfls_over_hfls>0.5, NA, delta_hfls_over_hfls)
delta_hfls_over_hfls <- ifel(delta_hfls_over_hfls< -0.5, NA, delta_hfls_over_hfls)
plot(delta_hfls_over_hfls)

hfls_hist <- hfls_ensemble$mean_hfls_1981_2005
hfls_fut <- hfls_ensemble$mean_hfls_2076_2100

lambda <- (2.501 * (10^6) - 2361 *  8.033441) / (10^6)
ET <- hfls_ensemble * 1 / lambda * 3600*24 / 10^6 *365.25

mean(values(ET))

plot(ifel(ET>300, ET, NA))

AI <- rast("../../era5_land_global/nc/ERA5-land_AI_1981-2005.nc")
low_veg_cov <- rast("../../era5_land_global/lai/low_vegetation_cover.tif")
low_veg_lai <- rast("../../era5_land_global/lai/low_vegetation_lai.tif")
high_veg_cov <- rast("../../era5_land_global/lai/high_vegetation_cover.tif")
high_veg_lai <- rast("../../era5_land_global/lai/high_vegetation_lai.tif")

LAI <- low_veg_cov*low_veg_lai + high_veg_cov*high_veg_lai

plot(ifel(LAI<1, NA, LAI))


# low_veg_lai_files <- list.files(
#   "../../era5_land_global/lai/LAI_low_veg",
#   pattern = "\\.(tif|nc)$",   # adjust as needed
#   full.names = TRUE
# )
# 
# low_veg_lai <- rast(low_veg_lai_files)
# low_veg_lai_avg <- mean(low_veg_lai)
# writeRaster(low_veg_lai_avg, "../../era5_land_global/lai/low_vegetation_lai.tif")
# 
# high_veg_lai_files <- list.files(
#   "../../era5_land_global/lai/LAI_high_veg",
#   pattern = "\\.(tif|nc)$",   # adjust as needed
#   full.names = TRUE
# )
# 
# high_veg_lai <- rast(high_veg_lai_files)
# high_veg_lai_avg <- mean(high_veg_lai)
# writeRaster(high_veg_lai_avg, "../../era5_land_global/lai/high_vegetation_lai.tif")

# --- build VPD_crit as you already do ---
Data <- rast("geff_vpd_regression_CMIP5_CMIP6.nc")

R2 <- Data$r2; a <- Data$a; b <- Data$b
VPD_crit <- (b^2) / (4 * a^2)
VPD_crit <- ifel(R2 < 0.75, NA, VPD_crit)
VPD_crit <- ifel(VPD_crit > 5, NA, VPD_crit)
VPD_hist <- rotate(Data$mean_vpd_1981_2005)
VPD_fut <- rotate(Data$mean_vpd_2076_2100)

geff_hist <- rotate(Data$mean_geff_1981_2005)
geff_fut <- rotate(Data$mean_geff_2076_2100)
geff_hist <- ifel(geff_hist<20, geff_hist, NA)


# --- put longitudes on -180..180 ---
# If your source is 0..360, rotate() is perfect:
VPD_crit <- rotate(VPD_crit)

# If you ever see -360..0 on the axis instead, do:
# VPD_crit <- wrap(VPD_crit, xmin = -180, xmax = 180)

# --- mask oceans (keep only land) ---
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land_v <- vect(land)
land_v <- project(land_v, crs(VPD_crit))        # match CRS
land_v <- crop(land_v, VPD_crit)                # optional speed-up

VPD_crit_land <- mask(VPD_crit, land_v)         # <- use the rotated raster

plot(ifel(VPD_crit_land>1,NA, VPD_crit_land))   # should show land only, -180..180

VPD_Tgt5C_hist_s <- resample(VPD_Tgt5C_hist, VPD_crit_land)
Ta_s <- resample(Ta, VPD_crit_land)
AI_s <- resample(AI, VPD_crit_land)
LAI_s <- resample(LAI, VPD_crit_land)

VPD_c <- crop(VPD_hist, ext(9,21,47,54))

VPD_vals <- values(VPD_c)

range(VPD_vals)

VPD_c <- crop(VPD_fut, ext(9,21,47,54))

VPD_vals <- values(VPD_c)


# The change
crop(VPD_hist, ext(9,21,47,54)) |> values() |> mean()
crop(VPD_crit_land, ext(9,21,47,54)) |> values() |> mean()
crop(VPD_fut, ext(9,21,47,54)) |> values() |> mean()

crop(VPD_hist, ext(9,21,47,54)) |> values() |> min()
crop(VPD_fut, ext(9,21,47,54)) |> values() |> max()


VPD_hist_range <- crop(VPD_hist, ext(9,21,47,54)) |> values() |> range(na.rm=TRUE)
AI_range <- crop(AI_s, ext(9,21,47,54)) |> values() |> range(na.rm=TRUE)
LAI_range <- crop(LAI_s, ext(9,21,47,54)) |> values() |> range(na.rm=TRUE)

0.327/0.41
0.60/0.41

range(VPD_vals)
quantile(VPD_vals, probs = c(0.25, 0.75), na.rm = TRUE)

range(values(crop(VPD_Tgt5C_hist_s, ext(9,21,47,54))))
quantile(values(crop(VPD_Tgt5C_hist_s, ext(9,21,47,54))), probs = c(0.05, 0.95), na.rm = TRUE)

# final_VPD <- ifel(VPD_hist<VPD_crit_land&VPD_fut>VPD_crit_land, VPD_crit_land, NA)
# final_VPD <- ifel(VPD_fut>1.05*VPD_crit_land, VPD_crit_land, NA)
# final_VPD <- ifel(VPD_hist<0.975*VPD_crit_land&VPD_fut>1.025*VPD_crit_land, VPD_crit_land, NA)
# final_VPD <- ifel(VPD_hist>0.22&VPD_hist<0.93, final_VPD, NA)
# final_VPD <- ifel(VPD_hist>0.2&VPD_hist<1, final_VPD, NA)
# final_VPD <- ifel(VPD_Tgt5C_hist_s>0.36&VPD_Tgt5C_hist_s<0.58, final_VPD, NA)
# final_VPD <- ifel(Ta_s>5, final_VPD, NA)
# final_VPD <- ifel(VPD_hist>VPD_hist_range[1]&VPD_hist<VPD_hist_range[2], final_VPD, NA)
# final_VPD <- ifel(AI_s>AI_range[1]&AI_s<AI_range[2], final_VPD, NA)
# final_VPD <- ifel(LAI_s>LAI_range[1]&LAI_s<LAI_range[2], final_VPD, NA)

# It goes over the peak
# final_VPD <- ifel(VPD_hist<0.95*VPD_crit_land&VPD_fut>1.05*VPD_crit_land, VPD_crit_land, NA)
final_VPD <- ifel(VPD_hist<VPD_crit_land&VPD_fut>VPD_crit_land, VPD_crit_land, NA)
# final_VPD <- ifel(VPD_hist<1.1*VPD_crit_land&VPD_fut>1.1*VPD_crit_land, VPD_crit_land, NA)

# It has reasonable vegetation
final_VPD <- ifel(LAI_s>1, final_VPD, NA)
plot(final_VPD, col = "red")
lines(land_v)


plot(ifel(hfls_fut<hfls_hist, 1, NA))
lines(land_v)

################################################################################
# Define required libraries
libs <- c(
  "elevatr", "terra", "tidyverse",
  "sf", "giscoR", "osmdata", "marmap", "tidyterra", "ggnewscale", "stringr"
)

# Install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == FALSE)) {
  install.packages(libs[!installed_libs])
}

# Load libraries
invisible(lapply(libs, library, character.only = TRUE))

VPD_ann_filt <- ifel(final_VPD > 0, 1, NA)

# Build a latitude raster on the same grid
lat_vals <- yFromCell(VPD_ann_filt, 1:ncell(VPD_ann_filt))   # latitude at each cell center
lat_r    <- setValues(rast(VPD_ann_filt), lat_vals)          # same geometry as VPD_ann_filt

# set cells sout/north of 66.5°N to NA
VPD_ann_filt[ abs(lat_r) > 66.5 ] <- NA

writeCDF(VPD_ann_filt, "AFFM.nc", overwrite=TRUE)

plot(VPD_ann_filt, col = "blue")
lines(land_v)

library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(terra)

# 0) Build mask polygons (value==1)
mask_poly <- terra::as.polygons(VPD_ann_filt, dissolve = TRUE, values = TRUE, na.out = TRUE)
mask_poly <- mask_poly[mask_poly[[1]] == 1, , drop = FALSE]
mask_sf   <- st_as_sf(mask_poly) |> st_make_valid()

# 1) LAND for clipping ONLY (not plotted)
land <- rnaturalearth::ne_download(scale = "medium", type = "land", category = "physical",
                                   returnclass = "sf") |>
  st_make_valid() |>
  st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"))
land_union <- st_union(land)

# Clip mask to land so oceans aren't hatched
mask_land <- suppressWarnings(st_intersection(mask_sf, land_union))
if (nrow(mask_land) == 0) stop("No land overlap in mask.")

# 2) COASTLINES for drawing borders (lines → no fill artifacts)
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") |>
  st_make_valid() |>
  st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"))

# 3) Work in a projected CRS (Equal Earth) for hatching + plot in that CRS
crs_proj <- 8857
mask_p   <- st_transform(mask_land, crs_proj)
coast_p  <- st_transform(coast,     crs_proj)

# ---- Neatline / projection frame (closed border) ----
mk_frame_ll <- function(xmin=-180, xmax=180, ymin=-85, ymax=85, step=0.25) {
  top    <- cbind(seq(xmin, xmax, by = step), ymax)
  right  <- cbind(xmax, seq(ymax, ymin, by = -step))
  bottom <- cbind(seq(xmax, xmin, by = -step), ymin)
  left   <- cbind(xmin, seq(ymin, ymax, by = step))
  ring   <- rbind(top, right, bottom, left, top[1, ])
  sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
}
frame_ll <- mk_frame_ll(step = 0.25)
frame_p  <- sf::st_transform(frame_ll, crs_proj)

# Padding so frame isn’t cropped
bb_fp   <- sf::st_bbox(frame_p)
pad_m   <- 200000
xlim_ex <- c(bb_fp["xmin"] - pad_m, bb_fp["xmax"] + pad_m)
ylim_ex <- c(bb_fp["ymin"] - pad_m, bb_fp["ymax"] + pad_m)

# ---- Graticule (lon/lat lines) → build in 4326, then project ----
grat_ll <- st_graticule(
  xlim = c(-180, 180),
  ylim = c(-85, 85),
  crs  = st_crs(4326),
  lon  = seq(-180, 180, by = 30),
  lat  = seq(-60, 60,  by = 15),
  ndiscr = 1e4
)
grat_p <- st_transform(grat_ll, crs_proj)

# Clip graticule to inset frame (e.g., 20 km inside)
frame_in <- sf::st_make_valid(sf::st_buffer(frame_p, dist = -20000))
grat_in  <- suppressWarnings(sf::st_intersection(grat_p, frame_in))

# -------- Graticule labels --------
lons <- seq(-150, 150, by = 30)
lats <- seq(-60,   60, by = 15)

y_bot <- -84.3
y_top <-  84.3
x_lft <- -178.5
x_rgt <-  178.5

fmt_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
fmt_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", ifelse(x > 0, "N", "")))

pts_lon_bot <- if (!is.null(y_bot)) sf::st_as_sf(data.frame(
  lon = lons, lat = y_bot, lab = fmt_lon(lons)
), coords = c("lon","lat"), crs = 4326) else NULL

pts_lon_top <- if (!is.null(y_top)) sf::st_as_sf(data.frame(
  lon = lons, lat = y_top, lab = fmt_lon(lons)
), coords = c("lon","lat"), crs = 4326) else NULL

pts_lat_lft <- if (!is.null(x_lft)) sf::st_as_sf(data.frame(
  lon = x_lft, lat = lats, lab = fmt_lat(lats)
), coords = c("lon","lat"), crs = 4326) else NULL

pts_lat_rgt <- if (!is.null(x_rgt)) sf::st_as_sf(data.frame(
  lon = x_rgt, lat = lats, lab = fmt_lat(lats)
), coords = c("lon","lat"), crs = 4326) else NULL

to_proj <- function(x) if (is.null(x)) NULL else sf::st_transform(x, crs_proj)
lab_lon_bot <- to_proj(pts_lon_bot)
lab_lon_top <- to_proj(pts_lon_top)
lab_lat_lft <- to_proj(pts_lat_lft)
lab_lat_rgt <- to_proj(pts_lat_rgt)

shift_xy <- function(sf_pts, dx = 0, dy = 0) {
  if (is.null(sf_pts)) return(NULL)
  crs_now <- sf::st_crs(sf_pts)
  m <- sf::st_coordinates(sf_pts)
  m[, 1] <- m[, 1] + dx
  m[, 2] <- m[, 2] + dy
  sf::st_as_sf(
    data.frame(sf_pts, X = m[, 1], Y = m[, 2]),
    coords = c("X", "Y"), crs = crs_now
  )
}

# inward label offsets (m)
off_tb <- -461000
off_lr <- -437000
lab_lon_top <- shift_xy(lab_lon_top, dy = -off_tb)
lab_lon_bot <- shift_xy(lab_lon_bot, dy =  off_tb)
lab_lat_lft <- shift_xy(lab_lat_lft, dx =  off_lr)
lab_lat_rgt <- shift_xy(lab_lat_rgt, dx = -off_lr)

# ---- Diagonal hatch pattern ----
bb <- st_bbox(mask_p)
spacing <- 180000

make_diag <- function(b, x1, x2) {
  st_linestring(matrix(c(x1, x1 + b, x2, x2 + b), ncol = 2, byrow = TRUE))
}
b_seq <- seq((bb["ymin"] - bb["xmax"]) - 5e6, (bb["ymax"] - bb["xmin"]) + 5e6, by = spacing)
lines_sfc <- st_sfc(lapply(b_seq, make_diag,
                           x1 = bb["xmin"] - 5e6, x2 = bb["xmax"] + 5e6),
                    crs = crs_proj)
bb_poly   <- st_as_sfc(bb) |> st_buffer(1e5)
lines_cut <- suppressWarnings(st_intersection(lines_sfc, bb_poly))
hatch_p   <- suppressWarnings(st_intersection(lines_cut, mask_p))

# ---- Split hatch lines by latitude bands ----
# Bands in lon/lat
trop_ll <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -23.5, ymax =  23.5), crs = 4326))
tempN_ll <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin =  23.5, ymax =  66.5), crs = 4326))
tempS_ll <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -66.5, ymax = -23.5), crs = 4326))
temp_ll  <- st_union(tempN_ll, tempS_ll)

# Reproject bands to map CRS
trop_p <- st_transform(trop_ll, crs_proj)
temp_p <- st_transform(temp_ll,  crs_proj)

# Intersect the already land-clipped hatch lines with bands
hatch_trop <- suppressWarnings(st_intersection(hatch_p, trop_p))
hatch_temp <- suppressWarnings(st_intersection(hatch_p, temp_p))

# Study domain – area of interest
AOI <- c(xmin = 9, xmax = 21, ymin = 47, ymax = 54)

# Create a high-resolution polygon for AOI
res <- 10
x_seq <- seq(AOI["xmin"], AOI["xmax"], length.out = res)
y_seq <- seq(AOI["ymin"], AOI["ymax"], length.out = res)
AOI_points <- rbind(
  cbind(x_seq, rep(AOI["ymin"], res)),  # Bottom side
  cbind(rep(AOI["xmax"], res), y_seq),  # Right side
  cbind(rev(x_seq), rep(AOI["ymax"], res)),  # Top side
  cbind(rep(AOI["xmin"], res), rev(y_seq))  # Left side
)

AOI_polygon <- st_as_sf(st_sfc(st_polygon(list(AOI_points))), crs = 4326)

# ---- Plot ----
p_h_AFFM_VPD <- ggplot() +
  geom_sf(data = grat_in, color = "grey55", linewidth = 0.3, linetype = "dotted") +
  # geom_sf_text(data = lab_lon_top, aes(label = lab), vjust = 0.9, size = 3, color = "grey20") +
  geom_sf_text(data = lab_lon_bot, aes(label = lab), vjust = 0.1, nudge_y = -120000, size = 3.2, color = "grey20") +
  geom_sf_text(data = lab_lat_lft, aes(label = lab), hjust = 0.95, nudge_x = -120000, size = 3.2, color = "grey20") +
  # geom_sf_text(data = lab_lat_rgt, aes(label = lab), hjust = 0.05, size = 3, color = "grey20") +
  # geom_sf(data = hatch_p, color = "#0a74da", linewidth = 0.6, alpha = 0.95, lineend = "round") +
  # Two-colored hatches by latitude band
  geom_sf(data = hatch_temp, color = "#0072B2", linewidth = 0.6, alpha = 0.95, lineend = "round") +
  geom_sf(data = hatch_trop, color = "#009E73", linewidth = 0.6, alpha = 0.95, lineend = "round") +
  geom_sf(data = coast_p, color = "grey40", linewidth = 0.32) +
  geom_sf(data = frame_p, fill = NA, color = "grey20", linewidth = 0.4) +
  # geom_sf(data = AOI_polygon, fill = "#FF4500", alpha = 0.1, color = NA) +    # faint translucent fill
  #geom_sf(data = AOI_polygon, fill = NA, color = "#FF4500", linewidth = 0.8) +   # inner stroke
  coord_sf(crs = sf::st_crs(crs_proj),
           xlim = xlim_ex, ylim = ylim_ex,
           expand = FALSE, clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("VPD_hatched_feedforward_mechanism.png", p_h_AFFM_VPD , width = 10, height = 4.9, dpi = 800, bg = "white")

# Save
saveRDS(p_h_AFFM_VPD, "p_h_AFFM_VPD.rds")

################################################################################

# ---- Plot (one-color filled mask) ----
p_h_AFFM_VPD <- ggplot() +
  geom_sf(data = grat_in, color = "grey65", linewidth = 0.25, linetype = "dotted") +
  
  # Mask as semi-transparent fill (no hatching)
  geom_sf(data = mask_p, fill = "#2C7FB8", color = NA, alpha = 0.35) +
  
  geom_sf(data = coast_p, color = "grey40", linewidth = 0.32) +
  geom_sf(data = frame_p, fill = NA, color = "grey20", linewidth = 0.4) +
  
  # Study domain box (optional but helpful)
  geom_sf(data = st_transform(AOI_polygon, crs_proj), fill = NA, color = "#8B1E1E", linewidth = 0.8) +
  
  coord_sf(crs = sf::st_crs(crs_proj),
           xlim = xlim_ex, ylim = ylim_ex,
           expand = FALSE, clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )


################################################################################
# ---- Prepare land for plotting (projected) ----
land_p <- st_transform(land_union, crs_proj)

# ---- Plot ----
p_h_AFFM_VPD <- ggplot() +
  
  # Ocean background = filled projection frame
  geom_sf(
    data = frame_p,
    fill = "#EEF3F8",
    color = NA
  ) +
  
  # Land fill (very light grey)
  geom_sf(
    data = land_p,
    fill = "#F2F2F2",
    color = NA
  ) +
  
  # Diagnostic mask (single-color, semi-transparent)
  geom_sf(
    data = mask_p,
    fill = "#2C7FB8",
    color = NA,
    alpha = 0.35
  ) +
  
  # Graticule (subtle)
  geom_sf(
    data = grat_in,
    color = "grey65",
    linewidth = 0.25,
    linetype = "dotted"
  ) +
  
  # Coastlines
  geom_sf(
    data = coast_p,
    color = "grey40",
    linewidth = 0.32
  ) +
  
  # Map frame
  geom_sf(
    data = frame_p,
    fill = NA,
    color = "grey20",
    linewidth = 0.4
  ) +
  
  # Central Europe study domain
  geom_sf(
    data = st_transform(AOI_polygon, crs_proj),
    fill = NA,
    color = "#8B1E1E",
    linewidth = 0.8
  ) +
  
  coord_sf(
    crs = sf::st_crs(crs_proj),
    xlim = xlim_ex,
    ylim = ylim_ex,
    expand = FALSE,
    clip = "off"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ---- Save ----
ggsave(
  "VPD_feedforward_diagnostic_map.png",
  p_h_AFFM_VPD,
  width = 10,
  height = 4.9,
  dpi = 800,
  bg = "white"
)

################################################################################
# ---- Plot: more visually appealing version ----

# IMPORTANT for polygon boolean ops with dateline-wrapped geometries
sf::sf_use_s2(FALSE)

# 1) Land in lon/lat, dateline-safe, unioned, valid
land_ll <- rnaturalearth::ne_download(
  scale = "medium", type = "land", category = "physical", returnclass = "sf"
) |>
  st_make_valid() |>
  st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")) |>
  st_make_valid()

land_union_ll <- st_union(land_ll) |> st_make_valid()

# 2) Your frame in lon/lat (you already have frame_ll from mk_frame_ll)
# Make sure it's valid too
frame_ll <- st_make_valid(frame_ll)

# 3) Compute ocean explicitly in lon/lat (this is the key change)
ocean_ll <- st_difference(frame_ll, land_union_ll) |> st_make_valid()

# 4) Now project both to Equal Earth (or whatever crs_proj is)
land_p  <- st_transform(land_union_ll, crs_proj)
ocean_p <- st_transform(ocean_ll,      crs_proj)


# 1
p_h_AFFM_VPD <- ggplot() +
  
  # Ocean background (clipped to globe)
  geom_sf(data = frame_p, fill = "#EEF3F8", color = NA) +
  
  # Land (very light)
  geom_sf(data = land_p, fill = "#F6F6F6", color = NA) +
  
  # Diagnostic mask (slightly stronger + optional subtle edge)
  geom_sf(data = mask_p, fill = "#3B82B8", color = NA, alpha = 0.60) +
  # Optional: a faint edge to help definition at small size
  # geom_sf(data = mask_p, fill = NA, color = "#3B82B8", linewidth = 0.12, alpha = 0.35) +
  
  # Graticule (lighter + thinner)
  geom_sf(data = grat_in, color = "grey75", linewidth = 0.18, linetype = "dotted") +
  
  # Coastlines (lighter and a bit thinner)
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  # Study domain: soft fill + crisp outline
  geom_sf(
    data = st_transform(AOI_polygon, crs_proj),
    fill = "#8B1E1E", alpha = 0.12,
    color = "#8B1E1E", linewidth = 0.85
  ) +
  
  # Frame outline (slightly lighter than before)
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  
  coord_sf(
    crs = sf::st_crs(crs_proj),
    xlim = xlim_ex, ylim = ylim_ex,
    expand = FALSE, clip = "off"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# 
# p_h_AFFM_VPD <- p_h_AFFM_VPD +
#   annotate("text", x = st_coordinates(st_centroid(st_transform(AOI_polygon, crs_proj)))[1],
#            y = st_coordinates(st_centroid(st_transform(AOI_polygon, crs_proj)))[2] + 500000,
#            label = "Study domain", size = 3.2, color = "#8B1E1E")

ggsave("VPD_feedforward_diagnostic_map_pretty_01.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 800, bg = "white")


# 2
p_h_AFFM_VPD <- ggplot() +
  geom_sf(data = frame_p, fill = "#F1F4F6", color = NA) +
  geom_sf(data = land_p,  fill = "#F5F3EF", color = NA) +
  
  geom_sf(data = mask_p,
          fill = "#4C8DAE",
          color = NA,
          alpha = 0.45) +
  
  geom_sf(data = grat_in, color = "grey78", linewidth = 0.18, linetype = "dotted") +
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  geom_sf(data = st_transform(AOI_polygon, crs_proj),
          fill = "#7A1F1F", alpha = 0.12,
          color = "#7A1F1F", linewidth = 0.85) +
  
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  coord_sf(crs = sf::st_crs(crs_proj),
           xlim = xlim_ex, ylim = ylim_ex,
           expand = FALSE, clip = "off") +
  theme_void()

ggsave("VPD_feedforward_diagnostic_map_pretty_02.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 800, bg = "white")

#3
p_h_AFFM_VPD <- ggplot() +
  geom_sf(data = frame_p, fill = "#EFEDE9", color = NA) +
  geom_sf(data = land_p,  fill = "#F3F6F8", color = NA) +

  geom_sf(data = mask_p,
          fill = "#5B7FA6",
          color = NA,
          alpha = 0.40) +
  
  geom_sf(data = grat_in, color = "grey78", linewidth = 0.18, linetype = "dotted") +
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  geom_sf(data = st_transform(AOI_polygon, crs_proj),
          fill = "#7A1F1F", alpha = 0.12,
          color = "#7A1F1F", linewidth = 0.85) +
  
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  coord_sf(crs = sf::st_crs(crs_proj),
           xlim = xlim_ex, ylim = ylim_ex,
           expand = FALSE, clip = "off") +
  theme_void()

ggsave("VPD_feedforward_diagnostic_map_pretty_03.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 800, bg = "white")

#4
p_h_AFFM_VPD <- ggplot() +
  geom_sf(data = frame_p, fill = "#F2F4F3", color = NA) +
  geom_sf(data = land_p,  fill = "#F4F2EE", color = NA) +

  # geom_sf(data = mask_p,
  #         fill = "#7A8F6B",
  #         color = NA,
  #         alpha = 0.60) +
  geom_sf(data = mask_p,
          fill = "#7A8F6B",
          color = "#7A8F6B",
          linewidth = 0.05,
          alpha = 0.45) +
  
  geom_sf(data = grat_in, color = "grey78", linewidth = 0.18, linetype = "dotted") +
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  geom_sf(data = st_transform(AOI_polygon, crs_proj),
          fill = "#6E1E1E", alpha = 0.12,
          color = "#6E1E1E", linewidth = 0.85) +
          # fill = "black", alpha = 0.12,
          # color = "black", linewidth = 0.85) +
  
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  coord_sf(crs = sf::st_crs(crs_proj),
           xlim = xlim_ex, ylim = ylim_ex,
           expand = FALSE, clip = "off") +
  theme_void()

ggsave("VPD_feedforward_diagnostic_map_pretty_04.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 800, bg = "white")

# 5
p_h_AFFM_VPD <- ggplot() +
  
  # Ocean background (clipped to globe)
  geom_sf(data = frame_p, fill = "#EEF3F8", color = NA) +
  
  # Land (very light)
  geom_sf(data = land_p, fill = "#F6F6F6", color = NA) +
  
  # Diagnostic mask (slightly stronger + optional subtle edge)
  geom_sf(data = mask_p, fill = "#7A8F6B", color = NA, alpha = 0.80) +
  # Optional: a faint edge to help definition at small size
  # geom_sf(data = mask_p, fill = NA, color = "#3B82B8", linewidth = 0.12, alpha = 0.35) +
  
  # Graticule (lighter + thinner)
  geom_sf(data = grat_in, color = "grey75", linewidth = 0.18, linetype = "dotted") +
  
  # Coastlines (lighter and a bit thinner)
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  # Study domain: soft fill + crisp outline
  geom_sf(
    data = st_transform(AOI_polygon, crs_proj),
    fill = "#8B1E1E", alpha = 0.12,
    color = "#8B1E1E", linewidth = 0.85
  ) +
  
  # Frame outline (slightly lighter than before)
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  
  coord_sf(
    crs = sf::st_crs(crs_proj),
    xlim = xlim_ex, ylim = ylim_ex,
    expand = FALSE, clip = "off"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# 
# p_h_AFFM_VPD <- p_h_AFFM_VPD +
#   annotate("text", x = st_coordinates(st_centroid(st_transform(AOI_polygon, crs_proj)))[1],
#            y = st_coordinates(st_centroid(st_transform(AOI_polygon, crs_proj)))[2] + 500000,
#            label = "Study domain", size = 3.2, color = "#8B1E1E")

ggsave("VPD_feedforward_diagnostic_map_pretty_05.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 1000, bg = "white")
################################################################################

# 5
p_h_AFFM_VPD <- ggplot() +
  
  # Ocean background (clipped to globe)
  geom_sf(data = frame_p, fill = "#EEF3F8", color = NA) +
  
  # Land (very light)
  geom_sf(data = land_p, fill = "#F6F6F6", color = NA) +
  
  # Diagnostic mask (slightly stronger + optional subtle edge)
  geom_sf(data = mask_p, fill = "#7A8F6B", color = NA, alpha = 0.80) +
  # Optional: a faint edge to help definition at small size
  # geom_sf(data = mask_p, fill = NA, color = "#3B82B8", linewidth = 0.12, alpha = 0.35) +
  
  # Graticule (lighter + thinner)
  geom_sf(data = grat_in, color = "grey75", linewidth = 0.18, linetype = "dotted") +
  
  # Coastlines (lighter and a bit thinner)
  geom_sf(data = coast_p, color = "grey45", linewidth = 0.26) +
  
  # Study domain: soft fill + crisp outline
  geom_sf(
    data = st_transform(AOI_polygon, crs_proj),
    fill = "#ff4500", alpha = 0.01,
    color = "#ff4500", linewidth = 0.4
  ) +
  
  # Frame outline (slightly lighter than before)
  geom_sf(data = frame_p, fill = NA, color = "grey35", linewidth = 0.35) +
  
  coord_sf(
    crs = sf::st_crs(crs_proj),
    xlim = xlim_ex, ylim = ylim_ex,
    expand = FALSE, clip = "off"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("VPD_feedforward_diagnostic_map_pretty_06.png", p_h_AFFM_VPD,
       width = 10, height = 4.9, dpi = 1000, bg = "white")
################################################################################












