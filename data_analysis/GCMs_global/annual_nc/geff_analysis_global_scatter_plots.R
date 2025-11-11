rm(list = ls())

library(tidyverse)
library(cowplot)
library(terra)

# # Path to your CSV
# csv_path <- "./geff_vpd_land.csv"
# 
# # Read (specify types so huge files load faster/reliably)
# dat <- readr::read_csv(
#   file = csv_path,
#   col_types = cols(
#     VPD = col_double(),
#     geff = col_double(),
#     model = col_character(),
#     period = col_character(),
#     continent = col_character(),
#     longitude = col_double(),
#     latitude = col_double()
#   )
# )
# 
# glimpse(dat)
# 
# eu <- dat %>% 
#   filter(continent == "Europe")
# 
# ggplot(eu, aes(x = VPD, y = geff)) +
#   geom_point(alpha = 0.1, size = 0.6) +
#   labs(
#     x = "VPD (kPa)",
#     y = "geff (mm s⁻¹)",
#     title = "geff vs VPD over Europe (all models & periods)"
#   ) +
#   xlim(0, 5) +
#   ylim(0, 20) +
#   theme_minimal()
# 
# 
# glob_no_antarctica <- dat |>  
#   filter(continent != "Antarctica")
# 
# ggplot(glob_no_antarctica, aes(x = VPD, y = geff)) +
#   geom_point(alpha = 0.1, size = 0.6) +
#   labs(
#     x = "VPD (kPa)",
#     y = "geff (mm s⁻¹)",
#     title = "geff vs VPD over globe (all models & periods)"
#   ) +
#   xlim(0, 5) +
#   ylim(0, 20) +
#   theme_minimal()
# 
# library(hexbin)
# ggplot(dat, aes(VPD, geff)) +
#   geom_bin2d(bins = 200) +
#   scale_fill_viridis_c() +
#   labs(
#     x = expression("VPD)"~(kPa)),
#     y = expression("geff (mm s"^{-1}*")"),
#     fill = "Count"
#   ) +
#   xlim(0, 5) +
#   ylim(0, 50) +
#   theme_minimal()
# 
# ggplot(dat, aes(VPD, geff)) +
#   geom_bin2d(
#     bins = 1000,
#     aes(fill = after_stat(ifelse(count < 10, NA, count)))  # hide low-count bins
#   ) +
#   scale_fill_viridis_c(na.value = NA) +
#   labs(
#     x = expression("VPD (kPa)"),
#     y = expression("geff (mm s"^{-1}*")"),
#     fill = "Count"
#   ) +
#   xlim(0, 5) +
#   ylim(0, 20) +
#   theme_minimal()


################################################################################
# 1) Read CSV
dat <- read_csv("./geff_vpd_land.csv",
                col_types = cols(
                  VPD = col_double(),
                  geff = col_double(),
                  model = col_character(),
                  period = col_character(),
                  continent = col_character(),
                  longitude = col_double(),
                  latitude = col_double()
                )
)

# 2) Match raster lon convention (-180..180 for your SpatRaster)
dat <- dat |> 
  mutate(lon_180 = if_else(longitude > 180, longitude - 360, longitude),
         lat = latitude)

# 3) Extract raster values at points (works across terra versions)
coords  <- as.matrix(dat[, c("lon_180", "lat")])

AFFM <- rast("./AFFM.nc")


#---------------
# Central Europe

AFFM_CE <- crop(AFFM, ext(9,21,47,54))

vals_df <- terra::extract(AFFM_CE , coords) |> as.data.frame()

# drop ID if present
if ("ID" %in% names(vals_df)) vals_df$ID <- NULL

# find the extracted column name (some versions use the layer name, others use a generic name)
layer_name <- names(AFFM )[1]
if (!layer_name %in% names(vals_df)) {
  # fallback to the first column in vals_df
  layer_name <- names(vals_df)[1]
}

# 4) Bind and filter where raster == 1
dat_sel <- bind_cols(dat, vals_df) |> 
  filter(.data[[layer_name]] == 1)

## --- Model B: geff ~ 1/sqrt(VPD) ---
## --- Ensure no bad values (esp. VPD <= 0 for 1/sqrt(VPD)) ---
dat_ok <- dat_sel |> 
  filter(is.finite(VPD), is.finite(geff)) |> 
  filter(VPD > 0) |>                # avoid Inf
  mutate(inv_sqrt_VPD = 1/sqrt(VPD))

fit_inv <- lm(geff ~ inv_sqrt_VPD, data = dat_ok)
sum_inv <- summary(fit_inv)
a_inv <- unname(coef(fit_inv)[1])           # intercept
b_inv <- unname(coef(fit_inv)[2])           # slope for 1/sqrt(VPD)
r2_inv <- unname(sum_inv$r.squared)

lab_inv <- sprintf("geff = %.2f + %.2f·(1/sqrt(VPD))\nR² = %.3f", a_inv, b_inv, r2_inv)

# Scatter + regression line + annotation
p_scatter <- ggplot(dat_ok, aes(x = inv_sqrt_VPD, y = geff)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(
    x = expression("1 / sqrt(VPD)"~(kPa^-0.5)),
    y = expression("geff (mm s"^{-1}*")"),
    title = "geff vs 1/sqrt(VPD) over Europe"
  ) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 30)) +
  theme_minimal() +
  annotate(
    "text",
    x = Inf, y = Inf, label = lab_inv,
    hjust = 1.02, vjust = 1.1, size = 3.6
  )

p_scatter

# Air density (kg m-3)
rhoAir <- 1.225

# Specific heat of air (J K-1kg-1),
CpAir <- 1004.67

# Scatter + regression line + annotation
p_scatter <- ggplot(dat_ok, aes(x = VPD, y = rhoAir * CpAir / 0.066 * VPD / 1000 * geff * 0.408 * 3600*24 / 10^6 *365.25)) +
  geom_point(alpha = 0.1, size = 0.6) +
  labs(
    x = expression("VPD"~(kPa)),
    y = expression("ET")
  ) +
  # coord_cartesian(xlim = c(0, 3), ylim = c(0, 10)) +
  theme_minimal() +
  annotate(
    "text",
    x = Inf, y = Inf, label = lab_inv,
    hjust = 1.02, vjust = 1.1, size = 3.6
  )

p_scatter

#-----------------------------------
# Tropics and temperate zone

AFFM_trop <- crop(AFFM, ext(xmin(AFFM), xmax(AFFM), -23.5, 23.5))
AFFM_tempN <- crop(AFFM, ext(xmin(AFFM), xmax(AFFM), 23.5, 66.5))
AFFM_tempS <- crop(AFFM, ext(xmin(AFFM), xmax(AFFM), -66.5, -23.5))
# put them on the same geometry and merge (NA outside each crop)
AFFM_temp  <- mosaic(AFFM_tempN, AFFM_tempS, fun = "max")  # or "mean", it won’t matter for a mask/value raster

################################################################################
# Temperate zone 

vals_df_temp <- terra::extract(AFFM_temp, coords) |> as.data.frame()
if ("ID" %in% names(vals_df_temp)) vals_df_temp$ID <- NULL
layer_name <- names(AFFM_temp)[1]
if (!(layer_name %in% names(vals_df_temp))) layer_name <- names(vals_df_temp)[1]

dat_sel_temp <- bind_cols(dat, vals_df_temp) |>
  filter(.data[[layer_name]] == 1)

dat_ok_temp <- dat_sel_temp |>
  filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
  mutate(inv_sqrt_VPD = 1/sqrt(VPD),
         zone = "Temperate")

# Tropical zone 
vals_df_trop <- terra::extract(AFFM_trop, coords) |> as.data.frame()
if ("ID" %in% names(vals_df_trop)) vals_df_trop$ID <- NULL
layer_name <- names(AFFM_trop)[1]
if (!(layer_name %in% names(vals_df_trop))) layer_name <- names(vals_df_trop)[1]

dat_sel_trop <- bind_cols(dat, vals_df_trop) |>
  filter(.data[[layer_name]] == 1)

dat_ok_trop <- dat_sel_trop |>
  filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
  mutate(inv_sqrt_VPD = 1/sqrt(VPD),
         zone = "Tropical")

# --- Combine
dat_both <- bind_rows(dat_ok_temp, dat_ok_trop)

fits_df <- dat_both |> 
  group_by(zone) |> 
  do({
    fit <- lm(geff ~ inv_sqrt_VPD, data = .)
    sm  <- summary(fit)
    tibble(
      a = unname(coef(fit)[1]),
      b = unname(coef(fit)[2]),
      r2 = sm$r.squared
    )
  }) |> 
  mutate(label = sprintf("y = %.1f + %.1f·x\nR² = %.2f", a, b, r2),
         x = Inf, y = Inf,
         hjust = 1.02,
         vjust = c(1.1, 2.4)[match(zone, unique(zone))]) # stagger labels

fits_df <- dat_both |> 
  group_by(zone) |> 
  do({
    fit <- lm(geff ~ inv_sqrt_VPD, data = .)
    sm  <- summary(fit)
    tibble(
      a = unname(coef(fit)[1]),
      b = unname(coef(fit)[2]),
      r2 = sm$r.squared
    )
  }) |> 
  mutate(label = sprintf("y = %.1f + %.1f·x\nR² = %.2f", a, b, r2),
         x = Inf, y = Inf,
         hjust = 1.02,
         vjust = c(1.1, 2.4)[match(zone, unique(zone))]) # stagger labels

dat_plot <- dat_both |> 
  dplyr::group_by(zone) |> 
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), floor(0.01 * nrow(.x))))) |> 
  dplyr::ungroup()


n_slices <- 50000

dat_slices_width <- dat_both |> 
  filter(is.finite(VPD), is.finite(geff), VPD > 0) |> 
  group_by(zone) |> 
  mutate(bin = cut_interval(VPD, n = n_slices)) |>    # equal-width bins within each zone
  group_by(zone, bin) |> 
  summarise(
    VPD   = mean(VPD,  na.rm = TRUE),
    VPD_sd   = mean(VPD,  na.rm = TRUE),
    inv_sqrt_VPD  = mean(inv_sqrt_VPD,  na.rm = TRUE),
    inv_sqrt_VPD_sd  = sd(inv_sqrt_VPD,  na.rm = TRUE),
    geff  = mean(geff, na.rm = TRUE),
    geff_sd    = sd(geff,   na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  )


# 1) Lock order of groups (edit strings to match yours exactly)
zone_order <- c("Temperate", "Tropical")
dat_slices_width <- dat_slices_width |> 
  mutate(zone = factor(zone, levels = zone_order))

# 2) Recompute fits (if needed)
fits_df <- dat_both |> 
  mutate(zone = factor(zone, levels = zone_order)) |> 
  group_by(zone) |> 
  do({
    fit <- lm(geff ~ I(1/sqrt(VPD)), data = .)
    sm  <- summary(fit)
    tibble(a = coef(fit)[1], b = coef(fit)[2], r2 = sm$r.squared)
  }) |>  ungroup()

# 3) FLIP positions: North -> bottom-right, South+Tropics -> upper-left
pos_df <- tibble::tibble(
  zone  = zone_order,
  x     = c( Inf, -Inf),   # flipped
  y     = c(-Inf,  Inf),   # flipped
  hjust = c(1.1, -0.1),
  vjust = c(-0.22,  1.2)
)

# 4) Two-line plain-text labels (two decimals)
fits_lbl <- fits_df |> 
  left_join(pos_df, by = "zone") |> 
  mutate(label = sprintf("y = %.1f + %.1fx\nR\u00B2 = %.2f", a, b, r2))

# define your two colors
cols <- c("#0072B2", "#009E73")

# name them to match your zone levels (so colors stick to labels)
zlv  <- unique(dat_slices_width$zone)          # keep current order
pal  <- setNames(cols[seq_along(zlv)], zlv)    # works if you have 1–2 zones
# If you ever have >2 zones, add more colors or subset zlv to the two you want

# 1) geff ~ 1/sqrt(VPD)
p_inv <- ggplot(dat_slices_width, aes(x = inv_sqrt_VPD, y = geff, color = zone)) +
  geom_point(alpha = 0.2, size = 0.5, shape = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = pal, name = NULL) +
  labs(
    x = bquote("1 / "~sqrt(VPD) ~ "(" * kPa^{-0.5} * ")"),
    y = bquote(g[eff] ~ "(" * mm ~ s^{-1} * ")")
  ) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 40)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()   # ← no grid lines
  ) +
  geom_text(
    data = fits_lbl,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust, color = zone),
    inherit.aes = FALSE, size = 3, lineheight = 1.05
  )

# 2) geff ~ VPD
p_vpd <- ggplot(dat_slices_width, aes(x = VPD, y = geff, color = zone)) +
  geom_point(alpha = 0.2, size = 0.5, shape = 3) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = bquote("VPD (kPa)"), y =  bquote(g[eff] ~ "(" * mm ~ s^{-1} * ")")) +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 40)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()   # ← no grid lines
  )

# 3) ET ~ VPD
rhoAir <- 1.225
CpAir  <- 1004.67

dat_slices_width <- dat_slices_width |> 
  mutate(ET = rhoAir * CpAir / 0.066 * VPD / 1000 * geff * 0.408 * 3600*24 / 10^6 * 365.25)


dat_plot <- dat_plot |> 
  mutate(ET = rhoAir * CpAir / 0.066 * VPD / 1000 * geff * 0.408 * 3600*24 / 10^6 * 365.25)

p_et <- ggplot(dat_slices_width, aes(x = VPD, y = ET, color = zone)) +
  geom_point(alpha = 0.2, size = 0.1, shape = 3) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = bquote("VPD (kPa)"), y = bquote(ET ~ "(" * mm ~ yr^{-1} * ")")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank()   # ← no grid lines
  )

ggsave(filename = "./geff_VPD.png", plot = p_vpd, width = 80, height = 80, dpi = 600, units = "mm")
ggsave(filename = "./geff_inv_sqrt_VPD.png", plot = p_inv, width = 80, height = 80, dpi = 600, units = "mm")
ggsave(filename = "./ET_VPD.png", plot = p_et, width = 80, height = 80, dpi = 600, units = "mm")

# Load the map
p_h_AFFM_VPD <- readRDS("p_h_AFFM_VPD.rds")

# Tiny, consistent margins
tight_bot <- theme(plot.margin = margin(3, 3, 2, 3))   # for the map (small bottom margin)
tight_top <- theme(plot.margin = margin(2, 3, 3, 3))   # for the bottom panels (small top margin)

p_map <- p_h_AFFM_VPD + tight_bot
p1    <- p_vpd        + tight_top
p2    <- p_inv        + tight_top
p3    <- p_et         + tight_top

# Align (keeps panel sizes/axes consistent)
aligned <- cowplot::align_plots(p_map, p1, p2, p3, align = "hv", axis = "tblr")

# 3) Top row: just the map (label a)
top_row <- cowplot::plot_grid(
  aligned[[1]],
  ncol = 1,
  labels = "a)",
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# 4) Bottom row: the three scatter/hex plots with small gaps
bottom_row <- cowplot::plot_grid(
  aligned[[2]], NULL, aligned[[3]], NULL, aligned[[4]],
  ncol = 5,
  rel_widths = c(1, 0.05, 1, 0.05, 1),
  labels = c("b)", "", "c)", "", "d)"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# 5) Stack rows with vertical spacing (tune rel_heights as you like)
combined <- cowplot::plot_grid(
  NULL,
  top_row,
  NULL,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.05, 1.7, 0, 1)  # map taller than the bottom row
)

combined_annot <- ggdraw(combined) +
  draw_label("Temperate", x = 0.05, y = 0.95,
             hjust = 0, vjust = 1, size = 12,
             fontface = "bold.italic", color = "#0072B2") +
  draw_label("Tropical",  x = 0.05, y = 0.91,
             hjust = 0, vjust = 1, size = 12,
             fontface = "bold.italic", color = "#009E73")

# 6) Save
ggsave("geff_vpd_panel_map_top_3bottom.png",
       combined_annot, width = 240, height = 180, units = "mm", dpi = 600, bg = "white")


################################################################################
# Bin averages
# -------- parameters --------

# ===============================
# Parameters (pick ONE approach)
# ===============================
bin_frac <- 0.05   # e.g., ~1% of rows per bin WITHIN EACH zone (set n_bins <- NULL)
n_bins   <- NULL   # OR set e.g. 25 and set bin_frac <- NULL

pal <- c("Temperate"="#0072B2", "Tropical"="#009E73")

# dat_both must have: VPD, geff, zone ("Temperate"/"Tropical")
stopifnot(all(c("VPD","geff","zone") %in% names(dat_both)))

# Derived variables
rhoAir <- 1.225; CpAir <- 1004.67
dat_both <- dat_both %>%
  filter(is.finite(VPD), is.finite(geff)) %>%
  mutate(
    inv_sqrt_VPD = 1/sqrt(VPD),
    ET = rhoAir * CpAir / 0.066 * VPD / 1000 * geff * 0.408 * 3600*24 / 1e6 * 365.25
  )

# ============================================================
# Helper: equal-count bins per group (zone), RELATIVE bin size
# - If bin_frac is given: bin_size_g = max(1, floor(n_g * bin_frac))
# - If n_bins   is given: bin_size_g = ceiling(n_g / n_bins)
# - Returns median + 25th/75th (x & y) and n per bin
# ============================================================
summarise_equalcount_rel <- function(df, x_var, y_var,
                                     bin_frac = NULL, n_bins = NULL,
                                     min_n = 20) {
  stopifnot(!is.null(bin_frac) || !is.null(n_bins))
  
  df0 <- df %>%
    mutate(.x = {{ x_var }}, .y = {{ y_var }}) %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    select(zone, .x, .y)
  
  # split by group, process each independently, then bind
  df_list <- df0 %>%
    group_split(zone, keep = TRUE)
  
  out <- map_df(df_list, function(gdf) {
    z <- gdf$zone[1]
    gdf <- arrange(gdf, .x)
    n_g <- nrow(gdf)
    
    bin_size_g <- if (!is.null(bin_frac)) {
      max(1L, floor(n_g * bin_frac))
    } else {
      ceiling(n_g / n_bins)
    }
    
    gdf <- gdf %>%
      mutate(idx = row_number(),
             bin_id = ceiling(idx / bin_size_g))
    
    gdf %>%
      group_by(bin_id) %>%
      summarise(
        zone  = z,
        x_med = median(.x, na.rm = TRUE),
        x_q25 = quantile(.x, 0.25, na.rm = TRUE, type = 7),
        x_q75 = quantile(.x, 0.75, na.rm = TRUE, type = 7),
        y_med = median(.y, na.rm = TRUE),
        y_q25 = quantile(.y, 0.25, na.rm = TRUE, type = 7),
        y_q75 = quantile(.y, 0.75, na.rm = TRUE, type = 7),
        n     = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(x_plot = x_med) %>%
      arrange(bin_id)
  })
  
  out
}

# Build summaries for each relationship
bin_inv <- summarise_equalcount_rel(dat_both, x_var = inv_sqrt_VPD, y_var = geff,
                                    bin_frac = bin_frac, n_bins = n_bins)
bin_vpd <- summarise_equalcount_rel(dat_both, x_var = VPD,           y_var = geff,
                                    bin_frac = bin_frac, n_bins = n_bins)
bin_et  <- summarise_equalcount_rel(dat_both, x_var = VPD,           y_var = ET,
                                    bin_frac = bin_frac, n_bins = n_bins)

# Plot function: median point + IQR bars on both axes
plot_binned_q <- function(d, x_lab, y_lab, xlim = NULL, ylim = NULL) {
  ggplot(d, aes(x = x_plot, y = y_med, color = zone)) +
    # horizontal IQR in x
    ggstance::geom_errorbarh(aes(xmin = x_q25, xmax = x_q75, y = y_med),
                             height = 0, size = 0.2, alpha = 0.9) +
    # vertical IQR in y
    geom_errorbar(aes(ymin = y_q25, ymax = y_q75),
                  width = 0, size = 0.2, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    { if (!is.null(xlim) || !is.null(ylim)) coord_cartesian(xlim = xlim, ylim = ylim) else coord_cartesian() } +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

# Build plots
p_inv_med <- plot_binned_q(
  bin_inv,
  x_lab = expression("1 / " * sqrt(VPD) ~ "(" * kPa^-0.5 * ")"),
  y_lab = expression(g[eff] ~ "(" * mm ~ s^-1 * ")"),
  xlim  = c(0, 4), ylim = c(0, 30)
)
p_vpd_med <- plot_binned_q(
  bin_vpd,
  x_lab = expression("VPD (kPa)"),
  y_lab = expression(g[eff] ~ "(" * mm ~ s^-1 * ")"),
  xlim  = c(0, 2.5), ylim = c(0, 30)
)
p_et_med <- plot_binned_q(
  bin_et,
  x_lab = expression("VPD (kPa)"),
  y_lab = expression(ET ~ "(" * mm ~ yr^-1 * ")"),
  xlim  = c(0, 2.5), ylim = c(300, 1700)
)

# Save if you want
ggsave("./geff_inv_sqrt_VPD_equalcountREL.png", p_inv_med, width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./geff_VPD_equalcountREL.png",         p_vpd_med, width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./ET_VPD_equalcountREL.png",           p_et_med,  width = 80, height = 80, dpi = 600, units = "mm", bg = "white")

# Load the map
p_h_AFFM_VPD <- readRDS("p_h_AFFM_VPD.rds")

# Tiny, consistent margins
tight_bot <- theme(plot.margin = margin(3, 3, 2, 3))   # for the map (small bottom margin)
tight_top <- theme(plot.margin = margin(2, 3, 3, 3))   # for the bottom panels (small top margin)

p_map <- p_h_AFFM_VPD + tight_bot
p1    <- p_vpd_med        + tight_top
p2    <- p_inv_med        + tight_top
p3    <- p_et_med         + tight_top

# Align (keeps panel sizes/axes consistent)
aligned <- cowplot::align_plots(p_map, p1, p2, p3, align = "hv", axis = "tblr")

# 3) Top row: just the map (label a)
top_row <- cowplot::plot_grid(
  aligned[[1]],
  ncol = 1,
  labels = "a)",
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# 4) Bottom row: the three scatter/hex plots with small gaps
bottom_row <- cowplot::plot_grid(
  aligned[[2]], NULL, aligned[[3]], NULL, aligned[[4]],
  ncol = 5,
  rel_widths = c(1, 0.05, 1, 0.05, 1),
  labels = c("b)", "", "c)", "", "d)"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# 5) Stack rows with vertical spacing (tune rel_heights as you like)
combined <- cowplot::plot_grid(
  NULL,
  top_row,
  NULL,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.05, 1.7, 0, 1)  # map taller than the bottom row
)

combined_annot <- ggdraw(combined) +
  draw_label("Temperate", x = 0.05, y = 0.95,
             hjust = 0, vjust = 1, size = 12,
             fontface = "bold.italic", color = "#0072B2") +
  draw_label("Tropical",  x = 0.05, y = 0.91,
             hjust = 0, vjust = 1, size = 12,
             fontface = "bold.italic", color = "#009E73")

# 6) Save
ggsave("geff_vpd_panel_map_top_3bottom_bin_med.png",
       combined_annot, width = 240, height = 180, units = "mm", dpi = 600, bg = "white")


################################################################################
# Normalize VPD and geff per each grid

# ===============================
# Colors and parameters
# ===============================
pal <- c("Temperate"="#0072B2", "Tropical"="#009E73")

# Binning controls (pick ONE style)
bin_frac <- 0.05   # relative bin size within each zone (e.g., 5% per bin)
n_bins   <- NULL   # OR set an integer (e.g., 25) and set bin_frac <- NULL

# ===============================
# Helper: equal-count bins per group (zone), RELATIVE bin size
# - If bin_frac is given: bin_size_g = max(1, floor(n_g * bin_frac))
# - If n_bins   is given: bin_size_g = ceiling(n_g / n_bins)
# - Returns median + 25th/75th (x & y) and n per bin
# ===============================
summarise_equalcount_rel <- function(df, x_var, y_var,
                                     bin_frac = NULL, n_bins = NULL,
                                     min_n = 20) {
  stopifnot(!is.null(bin_frac) || !is.null(n_bins))
  
  df0 <- df %>%
    mutate(.x = {{ x_var }}, .y = {{ y_var }}) %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    select(zone, .x, .y)
  
  df_list <- df0 %>%
    group_split(zone, keep = TRUE)
  
  out <- map_df(df_list, function(gdf) {
    z <- gdf$zone[1]
    gdf <- arrange(gdf, .x)
    n_g <- nrow(gdf)
    
    bin_size_g <- if (!is.null(bin_frac)) {
      max(1L, floor(n_g * bin_frac))
    } else {
      ceiling(n_g / n_bins)
    }
    
    gdf <- gdf %>%
      mutate(idx = row_number(),
             bin_id = ceiling(idx / bin_size_g))
    
    gdf %>%
      group_by(bin_id) %>%
      summarise(
        zone  = z,
        x_med = median(.x, na.rm = TRUE),
        x_q25 = quantile(.x, 0.25, na.rm = TRUE, type = 7),
        x_q75 = quantile(.x, 0.75, na.rm = TRUE, type = 7),
        y_med = median(.y, na.rm = TRUE),
        y_q25 = quantile(.y, 0.25, na.rm = TRUE, type = 7),
        y_q75 = quantile(.y, 0.75, na.rm = TRUE, type = 7),
        n     = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(x_plot = x_med) %>%
      arrange(bin_id)
  })
  
  out
}

# ===============================
# Extract + per-grid normalization
# ===============================

# ── Temperate zone: keep grid id from raster cell
vals_df_temp <- terra::extract(AFFM_temp, coords, cells = TRUE) |> 
  as.data.frame() |>
  dplyr::rename(grid_id = cell)       # keep the cell index as grid id

# (keep the AFFM layer name logic)
if ("ID" %in% names(vals_df_temp)) vals_df_temp$ID <- NULL
layer_name <- names(AFFM_temp)[1]
if (!(layer_name %in% names(vals_df_temp))) layer_name <- names(vals_df_temp)[1]

dat_sel_temp <- dplyr::bind_cols(dat, vals_df_temp) |>
  dplyr::filter(.data[[layer_name]] == 1)

# Physical constants for ET
rhoAir <- 1.225
CpAir  <- 1004.67
K <- rhoAir * CpAir / 0.066 * 0.408 * 3600*24 / 1e6 * 365.25   # overall factor

dat_ok_temp <- dat_sel_temp |>
  dplyr::filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
  dplyr::group_by(grid_id) |>
  dplyr::mutate(
    # per-grid maxima for normalization
    VPD_max   = max(VPD,  na.rm = TRUE),
    geff_max  = max(geff, na.rm = TRUE),
    # normalized inputs
    VPD_norm  = dplyr::if_else(VPD_max  > 0, VPD  / VPD_max,  NA_real_),
    geff_norm = dplyr::if_else(geff_max > 0, geff / geff_max, NA_real_),
    # physics: ET from ORIGINALS (units), then per-grid normalized ET
    ET        = K * (VPD/1000) * geff,
    ET_max    = max(ET, na.rm = TRUE),
    ET_norm   = dplyr::if_else(ET_max > 0, ET / ET_max, NA_real_),
    # convenience: inverse sqrt on normalized VPD
    inv_sqrt_VPD_norm = 1 / sqrt(VPD_norm),
    zone = "Temperate"
  ) |>
  dplyr::ungroup()

# ── Tropical zone: same treatment
vals_df_trop <- terra::extract(AFFM_trop, coords, cells = TRUE) |> 
  as.data.frame() |>
  dplyr::rename(grid_id = cell)

if ("ID" %in% names(vals_df_trop)) vals_df_trop$ID <- NULL
layer_name <- names(AFFM_trop)[1]
if (!(layer_name %in% names(vals_df_trop))) layer_name <- names(vals_df_trop)[1]

dat_sel_trop <- dplyr::bind_cols(dat, vals_df_trop) |>
  dplyr::filter(.data[[layer_name]] == 1)

dat_ok_trop <- dat_sel_trop |>
  dplyr::filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
  dplyr::group_by(grid_id) |>
  dplyr::mutate(
    VPD_max   = max(VPD,  na.rm = TRUE),
    geff_max  = max(geff, na.rm = TRUE),
    VPD_norm  = dplyr::if_else(VPD_max  > 0, VPD  / VPD_max,  NA_real_),
    geff_norm = dplyr::if_else(geff_max > 0, geff / geff_max, NA_real_),
    ET        = K * (VPD/1000) * geff,
    ET_max    = max(ET, na.rm = TRUE),
    ET_norm   = dplyr::if_else(ET_max > 0, ET / ET_max, NA_real_),
    inv_sqrt_VPD_norm = 1 / sqrt(VPD_norm),
    zone = "Tropical"
  ) |>
  dplyr::ungroup()

# ── Combine
dat_both <- dplyr::bind_rows(dat_ok_temp, dat_ok_trop)

# Quick sanity
stopifnot(all(c("VPD","geff","ET","VPD_norm","geff_norm","ET_norm",
                "inv_sqrt_VPD_norm","zone") %in% names(dat_both)))

# ===============================
# Build binned summaries
# ===============================
bin_inv <- summarise_equalcount_rel(
  dat_both, x_var = inv_sqrt_VPD_norm, y_var = geff_norm,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_vpd <- summarise_equalcount_rel(
  dat_both, x_var = VPD_norm,          y_var = geff_norm,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_et  <- summarise_equalcount_rel(
  dat_both, x_var = VPD,               y_var = ET,        # physical ET vs physical VPD
  bin_frac = bin_frac, n_bins = n_bins
)
# If you want normalized ET vs normalized VPD instead, swap to:
# bin_et  <- summarise_equalcount_rel(dat_both, x_var = VPD_norm, y_var = ET_norm,
#                                     bin_frac = bin_frac, n_bins = n_bins)

# ===============================
# Plot helper: median point + IQR bars
# ===============================
plot_binned_q <- function(d, x_lab, y_lab, xlim = NULL, ylim = NULL) {
  ggplot(d, aes(x = x_plot, y = y_med, color = zone)) +
    ggstance::geom_errorbarh(aes(xmin = x_q25, xmax = x_q75, y = y_med),
                             height = 0, size = 0.2, alpha = 0.9) +
    geom_errorbar(aes(ymin = y_q25, ymax = y_q75),
                  width = 0, size = 0.2, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

# ===============================
# Build plots
# ===============================
p_inv_med <- plot_binned_q(
  bin_inv,
  x_lab = expression("1 / " * sqrt(VPD[norm]) ~ "(" * kPa^-0.5 * ")"),
  y_lab = expression(g["eff norm"]),
  xlim  = c(0, 4), ylim = c(0, 1)
)
p_vpd_med <- plot_binned_q(
  bin_vpd,
  x_lab = expression("VPD"[norm]),
  y_lab = expression(g["eff norm"]),
  xlim  = c(0, 1), ylim = c(0, 1)
)
p_et_med <- plot_binned_q(
  bin_et,
  x_lab = expression("VPD (kPa)"),
  y_lab = expression(ET ~ "(" * mm ~ yr^-1 * ")"),
  xlim  = c(0, quantile(dat_both$VPD, 0.99, na.rm = TRUE)),
  ylim  = c(quantile(dat_both$ET, 0.01, na.rm = TRUE),
            quantile(dat_both$ET, 0.99, na.rm = TRUE))
)

# ===============================
# Save
# ===============================
ggsave("./geff_inv_sqrt_VPD_equalcountREL_norm.png", p_inv_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./geff_VPD_equalcountREL_norm.png",         p_vpd_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./ET_VPD_equalcountREL_physical.png",       p_et_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")

################################################################################

################################################################################
# Normalize VPD and geff per each grid — SINGLE ZONE
################################################################################

# ===============================
# Zone selection (edit these two)
# ===============================
zone_name <- "Global"   # or "Temperate" / "Tropical" / any label you prefer
AFFM_one  <- AFFM  # choose the mask for your single zone (e.g., AFFM_temp or AFFM_trop)

# ===============================
# Colors and parameters
# ===============================
pal <- c(setNames("#0072B2", zone_name))  # single color palette keyed by your zone

# Binning controls (pick ONE style)
bin_frac <- 0.05   # relative bin size (e.g., 5% of rows per bin)
n_bins   <- NULL   # OR set an integer (e.g., 25) and set bin_frac <- NULL

# ===============================
# Helper: equal-count bins per group (zone), RELATIVE bin size
# - If bin_frac is given: bin_size_g = max(1, floor(n_g * bin_frac))
# - If n_bins   is given: bin_size_g = ceiling(n_g / n_bins)
# - Returns median + 25th/75th (x & y) and n per bin
# ===============================
summarise_equalcount_rel <- function(df, x_var, y_var,
                                     bin_frac = NULL, n_bins = NULL,
                                     min_n = 20) {
  stopifnot(!is.null(bin_frac) || !is.null(n_bins))
  
  df0 <- df %>%
    mutate(.x = {{ x_var }}, .y = {{ y_var }}) %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    select(zone, .x, .y)
  
  df_list <- df0 %>%
    group_split(zone, keep = TRUE)
  
  out <- map_df(df_list, function(gdf) {
    z <- gdf$zone[1]
    gdf <- arrange(gdf, .x)
    n_g <- nrow(gdf)
    
    bin_size_g <- if (!is.null(bin_frac)) {
      max(1L, floor(n_g * bin_frac))
    } else {
      ceiling(n_g / n_bins)
    }
    
    gdf <- gdf %>%
      mutate(idx = row_number(),
             bin_id = ceiling(idx / bin_size_g))
    
    gdf %>%
      group_by(bin_id) %>%
      summarise(
        zone  = z,
        x_med = median(.x, na.rm = TRUE),
        x_q25 = quantile(.x, 0.25, na.rm = TRUE, type = 7),
        x_q75 = quantile(.x, 0.75, na.rm = TRUE, type = 7),
        y_med = median(.y, na.rm = TRUE),
        y_q25 = quantile(.y, 0.25, na.rm = TRUE, type = 7),
        y_q75 = quantile(.y, 0.75, na.rm = TRUE, type = 7),
        n     = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(x_plot = x_med) %>%
      arrange(bin_id)
  })
  
  out
}

# ===============================
# Extract + per-grid normalization (single zone)
# ===============================

# Keep grid id (raster cell index) to normalize within each grid
vals_df <- terra::extract(AFFM_one, coords, cells = TRUE) |> 
  as.data.frame() |>
  dplyr::rename(grid_id = cell)

# Keep the AFFM layer name logic
if ("ID" %in% names(vals_df)) vals_df$ID <- NULL
layer_name <- names(AFFM_one)[1]
if (!(layer_name %in% names(vals_df))) layer_name <- names(vals_df)[1]

# Apply the mask (keep points where the mask == 1)
dat_sel <- dplyr::bind_cols(dat, vals_df) |>
  dplyr::filter(.data[[layer_name]] == 1)

# Physical constants for ET
rhoAir <- 1.225
CpAir  <- 1004.67
K <- rhoAir * CpAir / 0.066 * 0.408 * 3600*24 / 1e6 * 365.25   # overall factor

# Per-grid normalization + ET (physical) + ET_norm
dat_both <- dat_sel |>
  dplyr::filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
  dplyr::group_by(grid_id) |>
  dplyr::mutate(
    # per-grid maxima for normalization
    VPD_max   = max(VPD,  na.rm = TRUE),
    geff_max  = max(geff, na.rm = TRUE),
    # normalized inputs
    VPD_norm  = dplyr::if_else(VPD_max  > 0, VPD  / VPD_max,  NA_real_),
    geff_norm = dplyr::if_else(geff_max > 0, geff / geff_max, NA_real_),
    # physics: ET from ORIGINALS (units), then per-grid normalized ET
    ET_norm       = K * (VPD_norm/1000) * geff_norm,
    # convenience: inverse sqrt on normalized VPD
    inv_sqrt_VPD_norm = 1 / sqrt(VPD_norm),
    zone = zone_name
  ) |>
  dplyr::ungroup()

# Quick sanity
stopifnot(all(c("VPD","geff","VPD_norm","geff_norm","ET_norm",
                "inv_sqrt_VPD_norm","zone") %in% names(dat_both)))

# ===============================
# Build binned summaries
# ===============================
bin_inv <- summarise_equalcount_rel(
  dat_both, x_var = inv_sqrt_VPD_norm, y_var = geff_norm,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_vpd <- summarise_equalcount_rel(
  dat_both, x_var = VPD_norm,          y_var = geff_norm,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_et  <- summarise_equalcount_rel(
  dat_both, x_var = VPD_norm,               y_var = ET_norm,        # physical ET vs physical VPD
  bin_frac = bin_frac, n_bins = n_bins
)
# If you want normalized ET vs normalized VPD instead, swap to:
# bin_et  <- summarise_equalcount_rel(dat_both, x_var = VPD_norm, y_var = ET_norm,
#                                     bin_frac = bin_frac, n_bins = n_bins)

# ===============================
# Plot helper: median point + IQR bars
# ===============================
plot_binned_q <- function(d, x_lab, y_lab, xlim = NULL, ylim = NULL) {
  ggplot(d, aes(x = x_plot, y = y_med, color = zone)) +
    ggstance::geom_errorbarh(aes(xmin = x_q25, xmax = x_q75, y = y_med),
                             height = 0, size = 0.2, alpha = 0.9) +
    geom_errorbar(aes(ymin = y_q25, ymax = y_q75),
                  width = 0, size = 0.2, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

# ===============================
# Build plots
# ===============================
p_inv_med <- plot_binned_q(
  bin_inv,
  x_lab = expression("1 / " * sqrt(VPD[norm]) ~ "(" * kPa^-0.5 * ")"),
  y_lab = expression(g["eff norm"]),
  xlim  = c(0, 4), ylim = c(0, 1)
)
p_vpd_med <- plot_binned_q(
  bin_vpd,
  x_lab = expression("VPD"[norm]),
  y_lab = expression(g["eff norm"]),
  xlim  = c(0, 1), ylim = c(0, 1)
)
p_et_med <- plot_binned_q(
  bin_et,
  x_lab = expression("VPD (kPa)"),
  y_lab = expression(ET_norm ~ "(" * mm ~ yr^-1 * ")"),
  xlim  = c(0, 1), ylim = c(0, 60)
)

# ===============================
# Save
# ===============================
ggsave("./geff_inv_sqrt_VPD_equalcountREL_norm_single.png", p_inv_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./geff_VPD_equalcountREL_norm_single.png",         p_vpd_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./ET_VPD_equalcountREL_physical_single.png",       p_et_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")


################################################################################
# Normalize g_eff by predicted g_eff at VPD_ref = 0.1 kPa (per grid)
################################################################################

# -------------------------------
# User choices
# -------------------------------
zone_name <- "Global"   # label for this zone
AFFM_one  <- AFFM_temp  # pick the mask raster for this zone
VPD_ref   <- 0.1        # kPa, reference VPD for normalization
pal       <- c(setNames("#0072B2", zone_name))

# Binning controls (pick ONE)
bin_frac <- 0.05
n_bins   <- NULL

# -------------------------------
# Helper: equal-count bins per zone
# -------------------------------
summarise_equalcount_rel <- function(df, x_var, y_var,
                                     bin_frac = NULL, n_bins = NULL,
                                     min_n = 20) {
  stopifnot(!is.null(bin_frac) || !is.null(n_bins))
  df0 <- df %>%
    mutate(.x = {{ x_var }}, .y = {{ y_var }}) %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    select(zone, .x, .y)
  df_list <- df0 %>% group_split(zone, keep = TRUE)
  
  map_df(df_list, function(gdf) {
    z <- gdf$zone[1]
    gdf <- arrange(gdf, .x)
    n_g <- nrow(gdf)
    bin_size_g <- if (!is.null(bin_frac)) max(1L, floor(n_g * bin_frac)) else ceiling(n_g / n_bins)
    
    gdf %>%
      mutate(idx = row_number(), bin_id = ceiling(idx / bin_size_g)) %>%
      group_by(bin_id) %>%
      summarise(
        zone  = z,
        x_med = median(.x, na.rm = TRUE),
        x_q25 = quantile(.x, 0.25, na.rm = TRUE, type = 7),
        x_q75 = quantile(.x, 0.75, na.rm = TRUE, type = 7),
        y_med = median(.y, na.rm = TRUE),
        y_q25 = quantile(.y, 0.25, na.rm = TRUE, type = 7),
        y_q75 = quantile(.y, 0.75, na.rm = TRUE, type = 7),
        n     = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(x_plot = x_med) %>%
      arrange(bin_id)
  })
}

# -------------------------------
# Extract grid_id + apply mask (single zone)
# -------------------------------
vals_df <- terra::extract(AFFM_one, coords, cells = TRUE) |>
  as.data.frame() |>
  dplyr::rename(grid_id = cell)

if ("ID" %in% names(vals_df)) vals_df$ID <- NULL
layer_name <- names(AFFM_one)[1]
if (!(layer_name %in% names(vals_df))) layer_name <- names(vals_df)[1]

dat_sel <- dplyr::bind_cols(dat, vals_df) |>
  dplyr::filter(.data[[layer_name]] == 1)

# -------------------------------
# Prepare data, fit per-grid linear model, get g_ref at VPD_ref
# -------------------------------
dat_base <- dat_sel %>%
  filter(is.finite(VPD), is.finite(geff), VPD > 0) %>%
  mutate(inv_sqrt_VPD = 1/sqrt(VPD),
         zone = zone_name)

# Fit geff ~ inv_sqrt_VPD per grid to get a and b
grid_fits <- dat_base %>%
  group_by(grid_id) %>%
  group_modify(~{
    df <- .x
    # need at least 3 points and some variation
    ok <- nrow(df) >= 3 && sd(df$inv_sqrt_VPD, na.rm = TRUE) > 0 && sd(df$geff, na.rm = TRUE) > 0
    if (!ok) return(tibble(a = NA_real_, b = NA_real_))
    fit <- try(lm(geff ~ inv_sqrt_VPD, data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble(a = NA_real_, b = NA_real_))
    co <- coef(fit)
    tibble(a = unname(co[1]), b = unname(co[2]))
  }) %>%
  ungroup()

# Evaluate g_ref = a + b * (1/sqrt(VPD_ref))
inv_ref <- 1/sqrt(VPD_ref)
grid_fits <- grid_fits %>%
  mutate(g_ref = a + b * inv_ref)

# Join g_ref back and build the ratio
dat_both <- dat_base %>%
  left_join(grid_fits %>% select(grid_id, g_ref), by = "grid_id") %>%
  mutate(
    geff_ratio = dplyr::if_else(is.finite(g_ref) & g_ref > 0, geff / g_ref, NA_real_)
  )

# Sanity check
stopifnot(all(c("VPD","geff","geff_ratio","zone") %in% names(dat_both)))

# -------------------------------
# Binned summaries: geff/geff(0.1) vs VPD
# -------------------------------
bin_ratio <- summarise_equalcount_rel(
  dat_both, x_var = VPD, y_var = geff_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)

# -------------------------------
# Plot helper
# -------------------------------
plot_binned_q <- function(d, x_lab, y_lab, xlim = NULL, ylim = NULL) {
  ggplot(d, aes(x = x_plot, y = y_med, color = zone)) +
    ggstance::geom_errorbarh(aes(xmin = x_q25, xmax = x_q75, y = y_med),
                             height = 0, size = 0.2, alpha = 0.9) +
    geom_errorbar(aes(ymin = y_q25, ymax = y_q75),
                  width = 0, size = 0.2, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

# ===============================
# (Optional) ET ratios using the same reference
# ===============================
dat_both <- dat_both %>%
  mutate(
    # Physical ET (if not computed yet, uncomment the next 2 lines)
    # rhoAir = 1.225; CpAir = 1004.67;
    # K = rhoAir * CpAir / 0.066 * 0.408 * 3600*24 / 1e6 * 365.25,
    ET_ref    = K * (VPD_ref/1000) * g_ref,                         # per-grid ET at VPD_ref
    ET_ratio  = dplyr::if_else(is.finite(ET_ref) & ET_ref > 0,
                               (VPD / VPD_ref) * geff_ratio, NA_real_)
  )

# ===============================
# Build binned summaries for the three plots
# ===============================
bin_inv <- summarise_equalcount_rel(
  dat_both, x_var = 1/sqrt(VPD), y_var = geff_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_vpd <- summarise_equalcount_rel(
  dat_both, x_var = VPD,          y_var = geff_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_et  <- summarise_equalcount_rel(
  dat_both, x_var = VPD,          y_var = ET_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)

# ===============================
# Plot helper (reuses your existing helper)
# ===============================
# plot_binned_q(...) already defined above

# Axis ranges from data (robust to extremes)
x_vpd_max  <- quantile(dat_both$VPD, 0.99, na.rm = TRUE)
y_ratio_hi <- quantile(dat_both$geff_ratio, 0.99, na.rm = TRUE)
y_etr_hi   <- quantile(dat_both$ET_ratio,   0.99, na.rm = TRUE)

# ===============================
# Build plots
# ===============================
p_inv_med <- plot_binned_q(
  bin_inv,
  x_lab = expression("1 / " * sqrt(VPD) ~ "(" * kPa^-0.5 * ")"),
  y_lab = bquote(g[eff] / g[eff] ~ (VPD==.(VPD_ref))),
  xlim  = c(0, quantile(1/sqrt(dat_both$VPD), 0.99, na.rm = TRUE)),
  ylim  = c(0, max(1.2, y_ratio_hi, na.rm = TRUE))
)

p_vpd_med <- plot_binned_q(
  bin_vpd,
  x_lab = expression("VPD (kPa)"),
  y_lab = bquote(g[eff] / g[eff] ~ (VPD==.(VPD_ref))),
  xlim  = c(0, x_vpd_max),
  ylim  = c(0, max(1.2, y_ratio_hi, na.rm = TRUE))
)

p_et_med <- plot_binned_q(
  bin_et,
  x_lab = expression("VPD (kPa)"),
  y_lab = bquote(ET / ET ~ (VPD==.(VPD_ref))),
  xlim  = c(0, x_vpd_max),
  ylim  = c(0, max(1.5, y_etr_hi, na.rm = TRUE))
)

# ===============================
# Save
# ===============================
ggsave("./geff_ratio_vs_inv_sqrtVPD_equalcount_single.png", p_inv_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./geff_ratio_vs_VPD_equalcount_single.png",        p_vpd_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./ET_ratio_vs_VPD_equalcount_single.png",          p_et_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")


################################################################################
# Normalize g_eff by predicted g_eff at VPD_ref = 0.1 kPa, per grid, for 2 zones
################################################################################

# ===============================
# Colors and parameters
# ===============================
pal <- c("Temperate"="#0072B2", "Tropical"="#009E73")

# Binning controls (pick ONE style)
bin_frac <- 0.05   # relative bin size within each zone (e.g., 5% per bin)
n_bins   <- NULL   # OR set an integer (e.g., 25) and set bin_frac <- NULL

# Reference VPD for normalization
VPD_ref <- 0.1   # kPa

# ===============================
# Helper: equal-count bins per group (zone), RELATIVE bin size
# ===============================
summarise_equalcount_rel <- function(df, x_var, y_var,
                                     bin_frac = NULL, n_bins = NULL,
                                     min_n = 20) {
  stopifnot(!is.null(bin_frac) || !is.null(n_bins))
  
  df0 <- df %>%
    mutate(.x = {{ x_var }}, .y = {{ y_var }}) %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    select(zone, .x, .y)
  
  df_list <- df0 %>%
    group_split(zone, keep = TRUE)
  
  map_df(df_list, function(gdf) {
    z <- gdf$zone[1]
    gdf <- arrange(gdf, .x)
    n_g <- nrow(gdf)
    
    bin_size_g <- if (!is.null(bin_frac)) {
      max(1L, floor(n_g * bin_frac))
    } else {
      ceiling(n_g / n_bins)
    }
    
    gdf %>%
      mutate(idx = row_number(),
             bin_id = ceiling(idx / bin_size_g)) %>%
      group_by(bin_id) %>%
      summarise(
        zone  = z,
        x_med = median(.x, na.rm = TRUE),
        x_q25 = quantile(.x, 0.25, na.rm = TRUE, type = 7),
        x_q75 = quantile(.x, 0.75, na.rm = TRUE, type = 7),
        y_med = median(.y, na.rm = TRUE),
        y_q25 = quantile(.y, 0.25, na.rm = TRUE, type = 7),
        y_q75 = quantile(.y, 0.75, na.rm = TRUE, type = 7),
        n     = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(x_plot = x_med) %>%
      arrange(bin_id)
  })
}

# ===============================
# Extract + model-based normalization helpers
# ===============================

# Physical constants for ET
rhoAir <- 1.225
CpAir  <- 1004.67
K <- rhoAir * CpAir / 0.066 * 0.408 * 3600*24 / 1e6 * 365.25   # overall factor

# Small utility to process a single zone block
process_zone <- function(mask_rast, zone_label) {
  # Keep grid id (raster cell index)
  vals_df <- terra::extract(mask_rast, coords, cells = TRUE) |> 
    as.data.frame() |>
    dplyr::rename(grid_id = cell)
  
  if ("ID" %in% names(vals_df)) vals_df$ID <- NULL
  layer_name <- names(mask_rast)[1]
  if (!(layer_name %in% names(vals_df))) layer_name <- names(vals_df)[1]
  
  dat_sel <- dplyr::bind_cols(dat, vals_df) |>
    dplyr::filter(.data[[layer_name]] == 1)
  
  dat_base <- dat_sel |>
    dplyr::filter(is.finite(VPD), is.finite(geff), VPD > 0) |>
    dplyr::mutate(inv_sqrt_VPD = 1/sqrt(VPD),
                  zone = zone_label)
  
  # Fit geff ~ inv_sqrt_VPD PER GRID to get a,b
  grid_fits <- dat_base |>
    dplyr::group_by(grid_id) |>
    dplyr::group_modify(~{
      df <- .x
      ok <- nrow(df) >= 3 &&
        stats::sd(df$inv_sqrt_VPD, na.rm = TRUE) > 0 &&
        stats::sd(df$geff, na.rm = TRUE) > 0
      if (!ok) return(tibble(a = NA_real_, b = NA_real_))
      fit <- try(stats::lm(geff ~ inv_sqrt_VPD, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) return(tibble(a = NA_real_, b = NA_real_))
      co <- stats::coef(fit)
      tibble(a = unname(co[1]), b = unname(co[2]))
    }) |>
    dplyr::ungroup()
  
  # Predicted reference at VPD_ref
  inv_ref <- 1/sqrt(VPD_ref)
  grid_fits <- grid_fits |> dplyr::mutate(g_ref = a + b * inv_ref)
  
  # Join back, compute ratios and ET variants
  dat_out <- dat_base |>
    dplyr::left_join(grid_fits |> dplyr::select(grid_id, g_ref), by = "grid_id") |>
    dplyr::mutate(
      geff_ratio = dplyr::if_else(is.finite(g_ref) & g_ref > 0, geff / g_ref, NA_real_),
      ET         = K * (VPD/1000) * geff,
      ET_ref     = K * (VPD_ref/1000) * g_ref,
      ET_ratio   = dplyr::if_else(is.finite(ET_ref) & ET_ref > 0,
                                  (VPD / VPD_ref) * geff_ratio, NA_real_)  # algebraic simplification
    )
  
  dat_out
}

# ── Temperate
dat_ok_temp <- process_zone(AFFM_temp, "Temperate")

# ── Tropical
dat_ok_trop <- process_zone(AFFM_trop, "Tropical")

# ── Combine
dat_both <- dplyr::bind_rows(dat_ok_temp, dat_ok_trop)

# Quick sanity
stopifnot(all(c("VPD","geff","geff_ratio","ET","ET_ratio","zone") %in% names(dat_both)))

# ===============================
# Build binned summaries (ratios)
# ===============================
bin_inv <- summarise_equalcount_rel(
  dat_both, x_var = 1/sqrt(VPD), y_var = geff_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_vpd <- summarise_equalcount_rel(
  dat_both, x_var = VPD,          y_var = geff_ratio,
  bin_frac = bin_frac, n_bins = n_bins
)
bin_et  <- summarise_equalcount_rel(
  # dat_both, x_var = VPD,          y_var = ET_ratio,
  dat_both, x_var = VPD,          y_var = ET,
  bin_frac = bin_frac, n_bins = n_bins
)

# ===============================
# Plot helper: median point + IQR bars
# ===============================
plot_binned_q <- function(d, x_lab, y_lab, xlim = NULL, ylim = NULL) {
  ggplot(d, aes(x = x_plot, y = y_med, color = zone)) +
    ggstance::geom_errorbarh(aes(xmin = x_q25, xmax = x_q75, y = y_med),
                             height = 0, size = 0.2, alpha = 0.9) +
    geom_errorbar(aes(ymin = y_q25, ymax = y_q75),
                  width = 0, size = 0.2, alpha = 0.9) +
    geom_point(size = 1.8) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

# Robust axis ranges
x_vpd_max   <- quantile(dat_both$VPD, 0.99, na.rm = TRUE)
x_inv_max   <- quantile(1/sqrt(dat_both$VPD), 0.99, na.rm = TRUE)
y_ratio_hi  <- quantile(dat_both$geff_ratio, 0.99, na.rm = TRUE)
y_etr_hi    <- quantile(dat_both$ET_ratio,   0.99, na.rm = TRUE)

# ===============================
# Build plots (ratios)
# ===============================
p_inv_med <- plot_binned_q(
  bin_inv,
  x_lab = expression("1 / " * sqrt(VPD) ~ "(" * kPa^-0.5 * ")"),
  y_lab = bquote(g[eff] / g[eff] ~ (VPD==.(VPD_ref))),
  xlim  = c(0, x_inv_max), ylim = c(0, max(1.2, y_ratio_hi, na.rm = TRUE))
)

p_vpd_med <- plot_binned_q(
  bin_vpd,
  x_lab = expression("VPD (kPa)"),
  y_lab = bquote(g[eff] / g[eff] ~ (VPD==.(VPD_ref))),
  xlim  = c(0, x_vpd_max), ylim = c(0, max(1.2, y_ratio_hi, na.rm = TRUE))
)

# p_et_med <- plot_binned_q(
#   bin_et,
#   x_lab = expression("VPD (kPa)"),
#   y_lab = bquote(ET / ET ~ (VPD==.(VPD_ref))),
#   xlim  = c(0, x_vpd_max), ylim = c(0, max(1.5, y_etr_hi, na.rm = TRUE))
# )

p_et_med <- plot_binned_q(
  bin_et,
  x_lab = expression("VPD (kPa)"),
  y_lab = bquote(ET~(mm~yr^"-1")),
  xlim  = c(0, x_vpd_max), ylim = c(0, 2000)
)

# ===============================
# Save
# ===============================
ggsave("./geff_ratio_vs_inv_sqrtVPD_equalcount_twozones.png", p_inv_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./geff_ratio_vs_VPD_equalcount_twozones.png",        p_vpd_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
ggsave("./ET_ratio_vs_VPD_equalcount_twozones.png",          p_et_med,
       width = 80, height = 80, dpi = 600, units = "mm", bg = "white")
