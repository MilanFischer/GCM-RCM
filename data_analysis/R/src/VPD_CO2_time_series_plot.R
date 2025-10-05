library(tidyverse)

Data <- read_csv("../inputs/CE_VPD_annual_domain_mean.csv")

# Keep only 2015–2100 (works even if you don't have all years)
plot_df <- Data |> 
  mutate(date = as.numeric(date)) |> 
  filter(date >= 2015, date <= 2100) |> 
  arrange(date)

# ---- CO2 time series (ppm) ----
p_co2 <- ggplot(plot_df, aes(x = date, y = CO2_SSP585)) +
  geom_line(linewidth = 1, color = "grey30") +
  labs(
    title = "Atmosferická koncentrace CO₂ (SSP5-8.5)",
    subtitle = "2015–2100",
    x = "Rok",
    y = "CO₂ (ppm)"
  ) +
  scale_x_continuous(breaks = seq(2015, 2100, by = 10), limits = c(2015, 2100)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/co2_timeseries.png", p_co2, width = 6, height = 4.5, units = "in", dpi = 600)

# ---- VPD time series (kPa) ----
p_vpd <- ggplot(plot_df, aes(x = date, y = VPD)) +
  geom_line(linewidth = 1, color = "steelblue") +
  labs(
    title = "Sytostní doplněk (D)",
    subtitle = "2015–2100",
    x = "Rok",
    y = "D (kPa)"
  ) +
  scale_x_continuous(breaks = seq(2015, 2100, by = 10), limits = c(2015, 2100)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("../plots/vpd_timeseries.png", p_vpd, width = 6, height = 4.5, units = "in", dpi = 600)

# ---- Stomatal conductance ----

# --- Inputs ---
VPD_seq <- seq(0.01, 5, 0.2)
K_ET <- (rhoAir * CpAir / gamma) * VPD_seq * (1/1000) *
  1 / ((2.501e6 - 2361 * 15) / 1e6) * (3600 * 24) / 1e6

gs_ref1 <- 6.5
gs_ref2 <- 6.2
gs_ref3 <- 4.5
m1 <- 0.25
m2 <- 0.42
m3 <- 0.60

# --- Stomatal conductance models ---
gs_1 <- pmin(gs_ref1 * (1 - m1 * log(VPD_seq)), gs_ref1)
gs_2 <- pmin(gs_ref2 * (1 - m2 * log(VPD_seq)), gs_ref2)
gs_3 <- pmin(gs_ref3 * (1 - m3 * log(VPD_seq)), gs_ref3)

# --- ET series ---
ET_1 <- gs_ref1 * K_ET
ET_2 <- pmin(gs_2, gs_ref2) * K_ET
ET_3 <- pmin(gs_3, gs_ref3) * K_ET

# --- Data pro loess a vyhlazení na jemném gridu ---
df_raw <- tibble(D = VPD_seq,
                 ref = ET_1,
                 mid = ET_2,
                 strong = ET_3)

# Jemnější grid jen pro zobrazovaný rozsah D (0–4 kPa)
grid <- tibble(D = seq(0, 5, by = 0.01))

# Loess fit pro každou sérii (stejné span jako dřív)
fit_ref    <- loess(ref    ~ D, data = df_raw, span = 2)
fit_mid    <- loess(mid    ~ D, data = df_raw, span = 2)
fit_strong <- loess(strong ~ D, data = df_raw, span = 2)

# Predikce na gridu
smooth_ref    <- predict(fit_ref,    newdata = grid)
smooth_mid    <- predict(fit_mid,    newdata = grid)
smooth_strong <- predict(fit_strong, newdata = grid)

# --- CAPPING: oříznout regulované křivky podle vyhlazené reference ---
smooth_mid_capped    <- pmin(smooth_mid,    smooth_ref)
smooth_strong_capped <- pmin(smooth_strong, smooth_ref)

# --- Dlouhý formát pro ggplot ---
plot_df <- grid |>
  mutate(`Bez regulace` = smooth_ref,
         `Střední regulace`        = smooth_mid_capped,
         `Silná regulace`          = smooth_strong_capped) |>
  pivot_longer(-D, names_to = "Scénář", values_to = "ET") |>
  mutate(
    Scénář = factor(Scénář,
                    levels = c("Bez regulace", "Střední regulace", "Silná regulace"))
  )

# --- Vykreslení: nejprve regulace, pak reference NAVRCHU ---
p_et <- ggplot() +
  # ostatní křivky nejdřív
  geom_line(
    data = filter(plot_df, Scénář != "Bez regulace"),
    aes(x = D, y = ET, color = Scénář),
    linewidth = 1.2
  ) +
  # referenční křivka nakonec (navrchu, trochu silnější)
  geom_line(
    data = filter(plot_df, Scénář == "Bez regulace"),
    aes(x = D, y = ET, color = Scénář),
    linewidth = 1.6
  ) +
  labs(
    title = "Závislost evapotranspirace (ET)\nna sytostním doplňku (D)\na regulaci průduchů",
    x = "D (kPa)",
    y = "ET (mm / den)"
  ) +
  # Okabe–Ito paleta
  scale_color_manual(values = c(
    "Bez regulace" = "#0072B2",
    "Střední regulace"        = "#009E73",
    "Silná regulace"          = "#D55E00"
  ), breaks = c("Bez regulace", "Střední regulace", "Silná regulace")) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 10)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", lineheight = 1.1),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), color = "grey30"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("../plots/ET_vs_VPD.png", p_et, width = 4.6, height = 5.4, units = "in", dpi = 600)