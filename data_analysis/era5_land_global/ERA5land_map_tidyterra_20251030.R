# Clear workspace and garbage collect
rm(list = ls()); gc()

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

# Define constants for Clausius-Clapeyron equation parameters (Katul et al. 2012)
a <- 0.6111
b <- 17.5
c <- 240.97

Ta <- rast("./nc/ERA5-land_T2m_daily_1981-2005_ymon.nc")
RH <- rast("./nc/ERA5-land_RH_daily_1981-2005_ymon.nc")
e_sat <- a * exp(b * Ta / (Ta + c)) # kPa
e <- e_sat * RH / 100 # kPa
VPD <- e_sat - e # kPa

days_in_m <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# weights that sum to 1
w <- days_in_m / sum(days_in_m)

# annual mean weighted by days in month
VPD_ann <- sum(VPD * w, na.rm = TRUE)

VPD_ann_filt <- ifel(VPD_ann>0.35,NA, ifel(VPD_ann<0.25, NA, 1))

plot(VPD_ann_filt)
lines(continents)

library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)

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

# 2.5) Graticule (lon/lat lines) → build in 4326, then project to Equal Earth
grat_ll <- st_graticule(
  xlim = c(-180, 180),
  ylim = c(-85, 85),          # crop a bit to avoid polar clutter
  crs  = st_crs(4326),
  lon  = seq(-180, 180, by = 30),  # meridians every 30°
  lat  = seq(-60, 60,  by = 15),   # parallels every 15°
  ndiscr = 1e4
)
grat_p <- st_transform(grat_ll, crs_proj)

# Make an inside buffer of the frame (e.g., 80 km inset)
frame_in <- sf::st_buffer(frame_p, dist = -20000) |> sf::st_make_valid()

# Clip the graticule to that inset frame
grat_in <- suppressWarnings(sf::st_intersection(grat_p, frame_in))

# -------- Graticule labels (nice lon/lat labels) --------
# what to label
lons <- seq(-150, 150, by = 30)    # meridians
lats <- seq(-60,   60, by = 15)    # parallels

# a little inside the neatline so text doesn't touch the stroke
y_bot <- -84.3  # bottom labels' latitude
y_top <-  84.3  # top labels' latitude  (set NULL to hide top labels)
x_lft <- -178.5 # left labels' longitude
x_rgt <-  178.5 # right labels' longitude (set NULL to hide right labels)

# degree-label helpers
fmt_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
fmt_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", ifelse(x > 0, "N", "")))

# build label point sf
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

# project to Equal Earth
to_proj <- function(x) if (is.null(x)) NULL else sf::st_transform(x, crs_proj)
lab_lon_bot <- to_proj(pts_lon_bot)
lab_lon_top <- to_proj(pts_lon_top)
lab_lat_lft <- to_proj(pts_lat_lft)
lab_lat_rgt <- to_proj(pts_lat_rgt)

# helper: shift an sf POINT layer by (dx, dy) in projected units
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

# how far to push labels inward (meters, Equal Earth)
off_tb <- -461000  # top/bottom
off_lr <- -437000  # left/right

lab_lon_top <- shift_xy(lab_lon_top, dy = -off_tb)  # move DOWN
lab_lon_bot <- shift_xy(lab_lon_bot, dy =  off_tb)  # move UP
lab_lat_lft <- shift_xy(lab_lat_lft, dx =  off_lr)  # move RIGHT
lab_lat_rgt <- shift_xy(lab_lat_rgt, dx = -off_lr)  # move LEFT

# ---- Neatline / projection frame (closed border) ----
# Build a dense lon/lat rectangle (-180..180, -85..85) as a ring, then project
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
bb_fp   <- sf::st_bbox(frame_p)
pad_m   <- 200000  # 200 km pad so the stroke isn't cut
xlim_ex <- c(bb_fp["xmin"] - pad_m, bb_fp["xmax"] + pad_m)
ylim_ex <- c(bb_fp["ymin"] - pad_m, bb_fp["ymax"] + pad_m)

# Build diagonal hatch lines and clip to mask
bb <- st_bbox(mask_p)
spacing <- 180000  # meters between stripes (~180 km). Tweak as you like.

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

"#D1495B"
"#FCA311"
"#E26D5A"
"#E53935"

p <- ggplot() +
  geom_sf(data = grat_in, color = "grey55", linewidth = 0.3, linetype = "dotted") +
  
  # labels (directionally offset; no overlap with neatline)
  geom_sf_text(data = lab_lon_top, aes(label = lab), vjust = 0.9, size = 3, color = "grey20") +
  geom_sf_text(data = lab_lon_bot, aes(label = lab), vjust = 0.1, size = 3, color = "grey20") +
  geom_sf_text(data = lab_lat_lft, aes(label = lab), hjust = 0.95, size = 3, color = "grey20") +
  geom_sf_text(data = lab_lat_rgt, aes(label = lab), hjust = 0.05, size = 3, color = "grey20") +
  
  geom_sf(data = hatch_p, color = "#0a74da", linewidth = 0.6, alpha = 0.95, lineend = "round") +
  geom_sf(data = coast_p, color = "grey40", linewidth = 0.32) +
  geom_sf(data = frame_p, fill = NA, color = "grey20", linewidth = 0.4) +
  coord_sf(crs = sf::st_crs(crs_proj), expand = FALSE, clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(1, 10, 1, 10, "mm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("VPD_hatched_coastlines.png", p, width = 10, height = 4.9, dpi = 800, bg = "white")
