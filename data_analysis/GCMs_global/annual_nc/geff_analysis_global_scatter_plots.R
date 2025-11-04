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
