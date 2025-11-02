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

VPD_d <- rast("./nc/ERA5-land_VPD_daily_1981-2005_ymon.nc")

days_in_m <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

H <- rast("./nc/ERA5-land_H_monthly_1981-2005_ymon.nc") / (3600 * 24) * -1
LE <- rast("./nc/ERA5-land_LE_monthly_1981-2005_ymon.nc") / (3600 * 24) * -1

ext(H) <- ext(Ta)
ext(LE) <- ext(Ta)

Rn_G <- H + LE

lambda <- (2.501 * (10^6) - 2361 * Ta) / (10^6)
W_to_mm <- 1 / lambda * 3600*24 / 10^6

# Evapotranspiration
ET <- W_to_mm * LE

u10 <- rast("./nc/ERA5-land_WS10m_daily_1981-2005_ymon.nc")
zo <- (0.123 * 0.12)
u2 <- log(2 / zo) / log(10 / zo) * u10
gamma <- rast("./nc/ERA5-land_gamma_daily_1981-2005_ymon.nc")
SVP <- b * c * e_sat/(c + Ta)^2
ra <- 118 / u2
rs <- 35 # Original value in Allen et al. (2005) is 45 s/m, it was reduced in order to ensure that ET is always smaller than ETr
rhoAir <- 1.225
CpAir <- 1004.67
PET_d <- (SVP * Rn_G + rhoAir * CpAir * VPD_d / ra) /
  (SVP + gamma * (1 + rs / ra)) * W_to_mm
days_in_m <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

PET <- PET_d * days_in_m

PET_ann <- sum(PET)

P <- rast("./nc/ERA5-land_TP_daily_1981-2005_ymon.nc") * days_in_m
P_ann <- sum(P)

AI <- PET_ann / P_ann

writeCDF(AI, "./nc/ERA5-land_AI_1981-2005.nc")

# weights that sum to 1
w <- days_in_m / sum(days_in_m)

# annual mean weighted by days in month
VPD_ann <- sum(VPD * w, na.rm = TRUE)

library(terra)
library(tidyterra)
library(ggplot2)
library(viridisLite)
library(scico)
library(scales)
library(colorspace)

# AOI rectangle (lon 9–21, lat 47–54) in EPSG:4326
AOI <- st_as_sfc(st_bbox(c(xmin = 9, ymin = 47, xmax = 21, ymax = 54), crs = 4326)) |>
  st_as_sf()

# Helper to normalize to [0,1] like your get_alt01()
get_vals01 <- function(x){
  rng <- range(x, finite = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# ==== VPD scale (your custom breaks; feel free to tweak) ====
# thresholds
vpd_thresh <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00, 3.50)

# colors aligned to thresholds (blue ≤0.4, green at 0.4, red-ish >0.4)
cols_vpd <- c(
  "#F7FBFF", # 0.00
  "#C6DBEF", # 0.10
  "#9ECAE1", # 0.20
  "#2CA25F", # 0.30
  "#9ec564", # 0.40  <-- green pivot
  "#FEE08B", # 0.50
  "#FDAE61", # 0.60
  "#F46D43", # 0.80
  "#D53E4F", # 1.00
  "#B2182B", # 1.50
  "#8C1515", # 2.00
  "#67000D"  # 3.50
)


vpd_thresh <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00, 3.50)

# cols_vpd <- c(
#   "#F7FBFF", # 0.00  very light blue
#   "#DCEBFF", # 0.10
#   "#B8D8F6", # 0.20
#   "#66C2A4", # 0.30  light teal-green
#   "#2CA25F", # 0.40  stronger green pivot
#   "#FEE08B", # 0.50  yellow
#   "#E2BE6E", # 0.60  yellow-tan
#   "#C99856", # 0.80  tan
#   "#B27941", # 1.00  light brown
#   "#945A2A", # 1.50  medium brown
#   "#77431F", # 2.00  dark brown
#   "#532C16"  # 3.50  very dark brown
# )
# 
# cols_vpd <- c(
#   "#F7FBFF", # 0.00
#   "#DCEBFF", # 0.10
#   "#B8D8F6", # 0.20
#   "#5BC8AF", # 0.30 
#   "#1B9E4B", # 0.40
#   "#CFE66D", # 0.50
#   "#E6C15F", # 0.60
#   "#C99856", # 0.80
#   "#B27941", # 1.00
#   "#945A2A", # 1.50
#   "#77431F", # 2.00
#   "#532C16"  # 3.50
# )

vals_vpd <- get_vals01(vpd_thresh)
lims_vpd <- range(vpd_thresh)
vpd_labels <- sprintf("%.2f", vpd_thresh)

# Equal-width legend boxes
box_widths <- rep(1, length(vpd_labels))

# Scale function (like your my_scale_fill_etopo_pos)
my_scale_fill_vpd <- function (...) {
  ggplot2::scale_fill_gradientn(
    colours = cols_vpd,
    breaks  = vpd_thresh,
    labels  = vpd_labels,
    values  = vals_vpd,
    limits  = lims_vpd,
    ...
  )
}

# Variable-width legend builder (same as yours, just with a VPD title by default)
build_variable_width_legend <- function(
    labels, fills,
    title = "VPD (kPa)",
    widths = NULL,
    label_size = 5.2,
    title_size = 6.2,
    box_lwd = 0.3,
    bar_height = 0.03,
    bottom_gap = 0.08,
    label_offset = 0.08,
    title_gap = 0.05,
    use_mono = TRUE
) {
  stopifnot(length(labels) == length(fills))
  if (is.null(widths)) widths <- nchar(labels)
  widths <- pmax(as.numeric(widths), 1e-6)
  
  df <- tibble::tibble(label = labels, fill = fills, w = widths) |>
    dplyr::mutate(xmax = cumsum(w),
                  xmin = dplyr::lag(xmax, default = 0),
                  xmid = (xmin + xmax)/2)
  
  sumw    <- sum(df$w)
  x_title <- -title_gap * sumw
  x_min   <- -(title_gap + 0.02) * sumw
  x_max   <-  sumw
  
  ggplot(df) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax,
          ymin = bottom_gap, ymax = bottom_gap + bar_height, fill = label),
      color = "grey65", linewidth = box_lwd, show.legend = FALSE
    ) +
    scale_fill_manual(values = stats::setNames(df$fill, df$label)) +
    geom_text(
      aes(x = xmid, y = bottom_gap + bar_height + label_offset, label = label),
      vjust = 0, hjust = 0.5, size = label_size,
      family = if (use_mono) "mono" else ""
    ) +
    annotate(
      "text", x = x_title, y = bottom_gap + bar_height/2,
      label = title, hjust = 1, vjust = 0.5, fontface = "plain", size = title_size
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 0, b = 0, l = 0))
}

# Global graticule (no reprojection, no crop)
graticule <- st_graticule(
  xlim = c(-180, 180),
  ylim = c(-90, 90),
  crs  = st_crs(4326),
  ndiscr = 10000
)

# --- MAP: VPD over the whole globe in EPSG:4326 ---
get_vpd_map <- function() {
  ggplot() +
    geom_spatraster(data = VPD_ann, maxcell = 10^9, interpolate = FALSE, alpha = 1, show.legend = FALSE) +
    my_scale_fill_vpd(guide = "none") +   # hide built-in legend
    # AOI box (fill + white casing + orange stroke)
    # geom_sf(data = AOI, fill = "#FF4500", alpha = 0.10, color = NA) +
    # geom_sf(data = AOI, fill = NA, color = "white",   linewidth = 0.4) +
    # geom_sf(data = AOI, fill = NA, color = "#FF4500", linewidth = 0.3) +
    geom_sf(data = AOI, fill = "black", alpha = 0.10, color = NA) +
    geom_sf(data = AOI, fill = NA, color = "white",   linewidth = 0.6) +
    geom_sf(data = AOI, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = graticule, color = "grey70", linewidth = 0.2) +
    geom_sf(data = graticule, color = "grey70", linewidth = 0.2) +
    coord_sf(crs = "EPSG:4326", expand = FALSE,
             xlim = c(-180, 180), ylim = c(-90, 90)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "mm")
    )
}

region_map <- get_vpd_map()

# Custom legend with equal-sized boxes
custom_leg <- build_variable_width_legend(
  labels      = vpd_labels,
  fills       = cols_vpd,
  title       = "Mean annual VPD (kPa)",
  widths      = box_widths,   # <- equal box widths
  bar_height  = 0.30,
  bottom_gap  = 0.30,
  label_size  = 4,
  title_size  = 4.6,
  label_offset = 0.13,
  title_gap    = 0.04,
  use_mono     = FALSE
)

legend_row <- cowplot::plot_grid(NULL, custom_leg, NULL, ncol = 3,
                                 rel_widths = c(0.33, 1, 0.0))

final_plot <- cowplot::plot_grid(region_map, legend_row, ncol = 1,
                                 rel_heights = c(1, 0.12))

ggsave("VPD_global.png", plot = final_plot, width = 10, height = 6, dpi = 600, bg = "white")

################################################################################

# AOI rectangle (lon 9–21, lat 47–54) in EPSG:4326
AOI <- st_as_sfc(st_bbox(c(xmin = 9, ymin = 47, xmax = 21, ymax = 54), crs = 4326)) |>
  st_as_sf()

# Helper to normalize to [0,1] like your get_alt01()
get_vals01 <- function(x){
  rng <- range(x, finite = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# ==== VPD scale (your custom breaks; feel free to tweak) ====
# thresholds
AI_thresh <- c(0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 4.818692e+03)

# colors aligned to thresholds (blue ≤0.4, green at 0.4, red-ish >0.4)
# cols_AI <- c(
#   "#F7FBFF", # 0.00
#   "#C6DBEF", # 0.10
#   "#9ECAE1", # 0.20
#   "#2CA25F", # 0.30
#   "#9ec564", # 0.40  <-- green pivot
#   "#FEE08B", # 0.50
#   "#FDAE61", # 0.60
#   "#F46D43", # 0.80
#   "#D53E4F", # 1.00
#   "#B2182B", # 1.50
#   "#8C1515", # 2.00
#   "#67000D"  # 3.50
# )

cols_AI <- c(
  "#F7FBFF", # 0.00  very light blue
  "#C6DBEF", # 0.10  light blue
  "#9ECAE1", # 0.20  soft blue
  "#6BAED6", # 0.30  medium blue   (new)
  "#2CA25F", # 0.40  deeper blue   (new)
  "#9ec564", # 0.50  green
  "#FEE08B", # 0.60  green pivot
  "#FDAE61", # 0.60
  "#F46D43", # 0.80
  "#D53E4F", # 1.00
  "#8C1515", # 2.00
  "#67000D"  # 3.50
)

# cols_vpd <- c(
#   "#F7FBFF", # 0.00  very light blue
#   "#DCEBFF", # 0.10
#   "#B8D8F6", # 0.20
#   "#66C2A4", # 0.30  light teal-green
#   "#2CA25F", # 0.40  stronger green pivot
#   "#FEE08B", # 0.50  yellow
#   "#E2BE6E", # 0.60  yellow-tan
#   "#C99856", # 0.80  tan
#   "#B27941", # 1.00  light brown
#   "#945A2A", # 1.50  medium brown
#   "#77431F", # 2.00  dark brown
#   "#532C16"  # 3.50  very dark brown
# )
# 
# cols_vpd <- c(
#   "#F7FBFF", # 0.00
#   "#DCEBFF", # 0.10
#   "#B8D8F6", # 0.20
#   "#5BC8AF", # 0.30
#   "#1B9E4B", # 0.40
#   "#CFE66D", # 0.50
#   "#E6C15F", # 0.60
#   "#C99856", # 0.80
#   "#B27941", # 1.00
#   "#945A2A", # 1.50
#   "#77431F", # 2.00
#   "#532C16"  # 3.50
# )

vals_AI <- get_vals01(AI_thresh)
lims_AI <- range(AI_thresh)
AI_labels <- sprintf("%.2f", AI_thresh)
AI_labels[length(AI_labels)] <- ">10"

# Equal-width legend boxes
box_widths <- rep(1, length(AI_labels))

# Scale function (like your my_scale_fill_etopo_pos)
my_scale_fill_AI <- function (...) {
  ggplot2::scale_fill_gradientn(
    colours = cols_AI,
    breaks  = AI_thresh,
    labels  = AI_labels,
    values  = vals_AI,
    limits  = lims_AI,
    ...
  )
}

# Variable-width legend builder (same as yours, just with a VPD title by default)
build_variable_width_legend <- function(
    labels, fills,
    title = "AI = PET /P",
    widths = NULL,
    label_size = 5.2,
    title_size = 6.2,
    box_lwd = 0.3,
    bar_height = 0.03,
    bottom_gap = 0.08,
    label_offset = 0.08,
    title_gap = 0.05,
    use_mono = TRUE
) {
  stopifnot(length(labels) == length(fills))
  if (is.null(widths)) widths <- nchar(labels)
  widths <- pmax(as.numeric(widths), 1e-6)
  
  df <- tibble::tibble(label = labels, fill = fills, w = widths) |>
    dplyr::mutate(xmax = cumsum(w),
                  xmin = dplyr::lag(xmax, default = 0),
                  xmid = (xmin + xmax)/2)
  
  sumw    <- sum(df$w)
  x_title <- -title_gap * sumw
  x_min   <- -(title_gap + 0.02) * sumw
  x_max   <-  sumw
  
  ggplot(df) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax,
          ymin = bottom_gap, ymax = bottom_gap + bar_height, fill = label),
      color = "grey65", linewidth = box_lwd, show.legend = FALSE
    ) +
    scale_fill_manual(values = stats::setNames(df$fill, df$label)) +
    geom_text(
      aes(x = xmid, y = bottom_gap + bar_height + label_offset, label = label),
      vjust = 0, hjust = 0.5, size = label_size,
      family = if (use_mono) "mono" else ""
    ) +
    annotate(
      "text", x = x_title, y = bottom_gap + bar_height/2,
      label = title, hjust = 1, vjust = 0.5, fontface = "plain", size = title_size
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 2, r = 0, b = 0, l = 0))
}

# Global graticule (no reprojection, no crop)
graticule <- st_graticule(
  xlim = c(-180, 180),
  ylim = c(-90, 90),
  crs  = st_crs(4326),
  ndiscr = 10000
)

# --- MAP: AI over the whole globe in EPSG:4326 ---
get_AI_map <- function() {
  ggplot() +
    geom_spatraster(data = AI, maxcell = 10^9, interpolate = FALSE, alpha = 1, show.legend = FALSE) +
    my_scale_fill_AI(guide = "none") +   # hide built-in legend
    # AOI box (fill + white casing + orange stroke)
    # geom_sf(data = AOI, fill = "#FF4500", alpha = 0.10, color = NA) +
    # geom_sf(data = AOI, fill = NA, color = "white",   linewidth = 0.4) +
    # geom_sf(data = AOI, fill = NA, color = "#FF4500", linewidth = 0.3) +
    geom_sf(data = AOI, fill = "black", alpha = 0.10, color = NA) +
    geom_sf(data = AOI, fill = NA, color = "white",   linewidth = 0.6) +
    geom_sf(data = AOI, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = graticule, color = "grey70", linewidth = 0.2) +
    geom_sf(data = graticule, color = "grey70", linewidth = 0.2) +
    coord_sf(crs = "EPSG:4326", expand = FALSE,
             xlim = c(-180, 180), ylim = c(-90, 90)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "mm")
    )
}

region_map <- get_AI_map()

# Custom legend with equal-sized boxes
custom_leg <- build_variable_width_legend(
  labels      = AI_labels,
  fills       = cols_AI,
  title       = "Mean annual AI = PET / P",
  widths      = box_widths,   # <- equal box widths
  bar_height  = 0.30,
  bottom_gap  = 0.30,
  label_size  = 4,
  title_size  = 4.6,
  label_offset = 0.13,
  title_gap    = 0.04,
  use_mono     = FALSE
)

legend_row <- cowplot::plot_grid(NULL, custom_leg, NULL, ncol = 3,
                                 rel_widths = c(0.33, 1, 0.0))

final_plot <- cowplot::plot_grid(region_map, legend_row, ncol = 1,
                                 rel_heights = c(1, 0.12))

ggsave("AI_global.png", plot = final_plot, width = 10, height = 6, dpi = 600, bg = "white")


# ------------------------------------------------------------------------------
AI_val <- values(AI)
VPD_val <- values(VPD_ann)

AI_val_s <- AI_val[VPD_val<1.5]
VPD_val_s <- VPD_val[VPD_val<1.5]
plot(VPD_val_s, AI_val_s)

Fit <- lm(AI_val_s ~ VPD_val_s)
abline(Fit, col = "red")

vpd_thresh * coef(Fit)[2] + coef(Fit)[1]
