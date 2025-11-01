# Clear workspace and garbage collect
rm(list = ls()); gc()

# Define required libraries
libs <- c(
  "elevatr", "terra", "tidyverse",
  "sf", "giscoR", "osmdata", "marmap", "tidyterra", "ggnewscale", "stringr", "cowplot"
)

# Install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == FALSE)) {
  install.packages(libs[!installed_libs])
}

# Load libraries
invisible(lapply(libs, library, character.only = TRUE))

# Extra lib used later
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)

# Projection & plot extents (same as your elevation script)
CRS <- 3035  # set 4326 or 3035
if (CRS == 4326) {
  X_range <- c(-25, 60)
  Y_range <- c(30, 75)
} else if (CRS == 3035) {
  X_range <- c(2300000, 6900000)
  Y_range <- c(710000, 5450000)
}

# Filenames include these; set sensible defaults
elev_zoom <- 6
shp_res   <- "01"
DPI       <- 1200

# Bounding box for fetching/cropping borders (in target CRS)
bbox_sf <- st_bbox(
  c(xmin = X_range[1], xmax = X_range[2], ymin = Y_range[1], ymax = Y_range[2]),
  crs = st_crs(CRS)
)

# Borders (reuse cached shapefile if you have it; otherwise fetch from giscoR)
if (file.exists(paste0("./vect/borders_v-", shp_res, ".shp"))) {
  borders_sf <- st_read(paste0("./vect/borders_v-", shp_res, ".shp"), quiet = TRUE)
  borders_sf <- st_transform(borders_sf, CRS) |> st_crop(bbox_sf)
} else {
  borders_sf <- giscoR::gisco_get_countries(year = "2020", epsg = "4326", resolution = shp_res) |>
    st_transform(CRS) |>
    st_crop(bbox_sf)
}

# AOI in lon/lat then transform to target CRS (same AOI as before)
AOI <- c(xmin = 9, xmax = 21, ymin = 47, ymax = 54)
# Create a high-resolution polygon for AOI
res <- 1000
x_seq <- seq(AOI["xmin"], AOI["xmax"], length.out = res)
y_seq <- seq(AOI["ymin"], AOI["ymax"], length.out = res)
AOI_points <- rbind(
  cbind(x_seq, rep(AOI["ymin"], res)),  # Bottom side
  cbind(rep(AOI["xmax"], res), y_seq),  # Right side
  cbind(rev(x_seq), rep(AOI["ymax"], res)),  # Top side
  cbind(rep(AOI["xmin"], res), rev(y_seq))  # Left side
)

AOI_polygon <- st_as_sf(st_sfc(st_polygon(list(AOI_points))), crs = 4326)
AOI_transformed <- st_transform(AOI_polygon, crs = CRS)


# Graticule in target CRS over the plot window
graticule <- st_graticule(
  xlim = X_range, ylim = Y_range, crs = st_crs(CRS), ndiscr = 10000
)

# Custom legend builder (same as your elevation script)
build_variable_width_legend <- function(
    labels, fills, title = "Legend", widths = NULL, label_size = 5.2, title_size = 6.2,
    box_lwd = 0.3, bar_height = 0.03, bottom_gap = 0.08, label_offset = 0.08,
    title_gap = 0.05, use_mono = TRUE
) {
  stopifnot(length(labels) == length(fills))
  if (is.null(widths)) widths <- nchar(labels)
  widths <- pmax(as.numeric(widths), 1e-6)
  
  df <- tibble::tibble(label = labels, fill = fills, w = widths) |>
    dplyr::mutate(xmax = cumsum(w), xmin = dplyr::lag(xmax, default = 0), xmid = (xmin + xmax)/2)
  
  sumw <- sum(df$w); x_title <- -title_gap * sumw; x_min <- -(title_gap + 0.02) * sumw; x_max <- sumw
  
  ggplot(df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = bottom_gap, ymax = bottom_gap + bar_height, fill = label),
              color = "grey50", linewidth = box_lwd, show.legend = FALSE) +
    scale_fill_manual(values = stats::setNames(df$fill, df$label)) +
    geom_text(aes(x = xmid, y = bottom_gap + bar_height + label_offset, label = label),
              vjust = 0, hjust = 0.5, size = label_size, family = if (use_mono) "mono" else "") +
    annotate("text", x = x_title, y = bottom_gap + bar_height/2,
             label = title, hjust = 1, vjust = 0.5, fontface = "plain", size = title_size) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 0, b = 0, l = 0))
}

# Define constants for Clausius-Clapeyron equation parameters (Katul et al. 2012)
a <- 0.6111
b <- 17.5
c <- 240.97

# ---- EUROPE-ONLY VPD & AI MAPS (same design as elevation) --------------------

# Helper to normalize to [0,1]
get_vals01 <- function(x){
  rng <- range(x, finite = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# ---- Load/compute VPD annual (day-weighted) ----------------------------------
days_in_m <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

Ta <- rast("./nc/ERA5-land_T2m_daily_1981-2005_ymon.nc")        # °C
RH <- rast("./nc/ERA5-land_RH_daily_1981-2005_ymon.nc")          # %
e_sat <- a * exp(b * Ta / (Ta + c))                               # kPa
e     <- e_sat * RH / 100                                         # kPa
VPD   <- e_sat - e                                                # kPa

# weights sum to 1 (days per month)
w <- days_in_m / sum(days_in_m)
VPD_ann <- sum(VPD * w, na.rm = TRUE)

# ---- Load/compute AI = PET/P annual ------------------------------------------
# PET (Penman-Monteith, monthly → annual)
H   <- rast("./nc/ERA5-land_H_monthly_1981-2005_ymon.nc")  / (3600*24) * -1
LE  <- rast("./nc/ERA5-land_LE_monthly_1981-2005_ymon.nc") / (3600*24) * -1
ext(H)  <- ext(Ta); ext(LE) <- ext(Ta)
Rn_G <- H + LE

lambda   <- (2.501e6 - 2361 * Ta) / 1e6
W_to_mm  <- 1 / lambda * 3600*24 / 1e6
ET       <- W_to_mm * LE

u10   <- rast("./nc/ERA5-land_WS10m_daily_1981-2005_ymon.nc")
zo    <- (0.123 * 0.12)
u2    <- log(2/zo) / log(10/zo) * u10
gamma <- rast("./nc/ERA5-land_gamma_daily_1981-2005_ymon.nc")
SVP   <- b * c * e_sat/(c + Ta)^2
ra    <- 118 / u2
rs    <- 35
rhoAir <- 1.225
CpAir  <- 1004.67
VPD_d  <- rast("./nc/ERA5-land_VPD_daily_1981-2005_ymon.nc")

PET_d <- (SVP * Rn_G + rhoAir * CpAir * VPD_d / ra) /
  (SVP + gamma * (1 + rs / ra)) * W_to_mm
PET   <- PET_d * days_in_m
PET_ann <- sum(PET)

P      <- rast("./nc/ERA5-land_TP_daily_1981-2005_ymon.nc") * days_in_m
P_ann  <- sum(P)

AI <- PET_ann / P_ann

# ---- Project & crop to Europe (match your elevation map extents) -------------
target_crs <- paste0("EPSG:", CRS)
VPD_eu <- project(VPD_ann, target_crs, method = "bilinear") |> crop(ext(c(X_range, Y_range)))
AI_eu  <- project(AI,      target_crs, method = "bilinear") |> crop(ext(c(X_range, Y_range)))

# (Optional) align to elevation resolution for crisp overlay/consistent pixels
# VPD_eu <- resample(VPD_eu, elevation, method = "bilinear")
# AI_eu  <- resample(AI_eu,  elevation, method = "bilinear")

# ---- Color scales & legends ---------------------------------------------------
# VPD thresholds & palette (from your global map)
vpd_thresh  <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00, 3.50)
cols_vpd <- c(
  "#F7FBFF","#C6DBEF","#9ECAE1","#2CA25F","#9ec564",
  "#FEE08B","#FDAE61","#F46D43","#D53E4F","#B2182B","#8C1515","#67000D"
)
vals_vpd <- get_vals01(vpd_thresh); lims_vpd <- range(vpd_thresh)
vpd_labels <- sprintf("%.1f", vpd_thresh)
box_widths_equal <- rep(1, length(vpd_labels))

my_scale_fill_vpd <- function (...) {
  ggplot2::scale_fill_gradientn(
    colours = cols_vpd,
    breaks  = vpd_thresh,
    labels  = vpd_labels,
    values  = vals_vpd,
    limits  = lims_vpd,
    na.value = "grey85",   # <- ocean color here
    ...
  )
}

# AI thresholds & palette (from your AI global section; last label as >10)
AI_thresh  <- c(0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 4.818692e+03)
cols_AI <- c(
  "#F7FBFF","#C6DBEF","#9ECAE1","#6BAED6","#2CA25F",
  "#9ec564","#FEE08B","#FDAE61","#F46D43","#D53E4F","#8C1515","#67000D"
)
vals_AI <- get_vals01(AI_thresh); lims_AI <- range(AI_thresh)
AI_labels <- sprintf("%.1f", AI_thresh); AI_labels[length(AI_labels)] <- ">10"

my_scale_fill_AI <- function (...) {
  ggplot2::scale_fill_gradientn(
    colours = cols_AI,
    breaks  = AI_thresh,
    labels  = AI_labels,
    values  = vals_AI,
    limits  = lims_AI,
    na.value = "grey85",   # <- ocean color here
    ...
  )
}

# ---- Reusable map builder (same design as elevation) --------------------------
make_thematic_map <- function(r, scale_fun) {
  ggplot() +
    # # ocean/panel background
    # annotate("rect",
    #          xmin = X_range[1], xmax = X_range[2],
    #          ymin = Y_range[1], ymax = Y_range[2],
    #          fill = "grey85", colour = NA
    # ) +
    geom_spatraster(data = r, maxcell = 10^9, interpolate = FALSE, alpha = 1, show.legend = FALSE) +
    scale_fun(guide = "none") +
    xlim(X_range) + ylim(Y_range) +
    geom_sf(data = borders_sf, fill = NA, linewidth = 0.1) +
    geom_sf(data = graticule, color = "dodgerblue3", linewidth = 0.2, linetype = "solid") +
    geom_sf(data = AOI_transformed, fill = "#FF4500", alpha = 0.1, color = NA) +   # translucent fill
    geom_sf(data = AOI_transformed, fill = NA, color = "white",   linewidth = 1.4) + # casing
    geom_sf(data = AOI_transformed, fill = NA, color = "#FF4500", linewidth = 1.0) + # inner stroke
    coord_sf(crs = CRS, expand = FALSE) +
    labs(x = "", y = "", title = "", subtitle = "", caption = "") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
}

# ---- VPD map + legend ---------------------------------------------------------
vpd_map <- make_thematic_map(VPD_eu, my_scale_fill_vpd)

vpd_leg <- build_variable_width_legend(
  labels = vpd_labels,
  fills  = cols_vpd,
  title  = "VPD (kPa)",
  widths = box_widths_equal,
  bar_height  = 0.30,
  bottom_gap  = 0.30,
  label_size  = 5.3,
  title_size  = 5.7,
  label_offset = 0.16,
  title_gap    = 0.04,
  use_mono     = FALSE
)

vpd_legend_row <- cowplot::plot_grid(NULL, vpd_leg, NULL, ncol = 3, rel_widths = c(0.33, 1, 0.0))
vpd_final <- cowplot::plot_grid(vpd_map, vpd_legend_row, ncol = 1, rel_heights = c(1, 0.1))

ggsave(
  filename = paste0("./plots/VPD_Europe_r-", elev_zoom, "_v-", shp_res, "_DPI-", DPI, ".png"),
  width = 8, height = 8.9, dpi = DPI, device = "png",
  plot = vpd_final, bg = "white"
)

# ---- AI map + legend ----------------------------------------------------------
ai_map <- make_thematic_map(AI_eu, my_scale_fill_AI)

ai_leg <- build_variable_width_legend(
  labels = AI_labels,
  fills  = cols_AI,
  title  = " AI = PET / P",
  widths = rep(1, length(AI_labels)),
  bar_height  = 0.30,
  bottom_gap  = 0.30,
  label_size  = 5.3,
  title_size  = 5.7,
  label_offset = 0.16,
  title_gap    = 0.04,
  use_mono     = FALSE
)

ai_legend_row <- cowplot::plot_grid(NULL, ai_leg, NULL, ncol = 3, rel_widths = c(0.33, 1, 0.0))
ai_final <- cowplot::plot_grid(ai_map, ai_legend_row, ncol = 1, rel_heights = c(1, 0.1))

ggsave(
  filename = paste0("./plots/AI_Europe_r-", elev_zoom, "_v-", shp_res, "_DPI-", DPI, ".png"),
  width = 8, height = 8.9, dpi = DPI, device = "png",
  plot = ai_final, bg = "white"
)
