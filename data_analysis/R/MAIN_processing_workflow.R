# Clear the workspace
rm(list=ls()); gc()

library(DEoptim)
# library(ggplot2)
library(grid)
library(patchwork)
library(gridExtra)
library(ggrepel)
library(cowplot)
# library(scales)
# library(png)
library(magick)
# library(gtable)
library(MASS)
library(ggrepel)
library(numDeriv)
library(tidyverse)

# To set the RStudio theme like VS Code (my personal preference) use:
# https://github.com/anthonynorth/rscodeio
# https://github.com/anthonynorth/rscodeio/blob/master/inst/resources/rscodeio.rstheme
# If the installation prompts do not work, copy the file rscodeio.rstheme into %AppData%\RStudio\themes\ ion Windows or ~/.config/rstudio/themes/ Linux
# Then use Tools → Global Options → Appearance → Add… and select the file

# Preprocess
source("./src/general_functions.R")
source("./src/make_plot_tidy.R")
source("./src/colors.R")
source("./src/resize_plot_elements.R")
source("./src/Budyko_curve_tidy.R")
source("./src/pair_periods.R")

# VPD from daily values is in average about 15% higher, so constant 1.15 can be applied

# Plotting setup
# COL_CMIP5 <- "#117733"
# COL_CMIP6 <- "#AA4499"
# COL_RCMs <- "#6699CC"

COL_CMIP5 <- "#1c91df"
COL_CMIP6 <- "#4e79a7"
COL_RCMs <- "#ff4500"
BORDER_COLOR <- "#2b2b2b"

c("#1f77b4", "#2a9df4", "#4e79a7", "#5f9ea0", "#1c91df", "#4682b4",
  "#0073cf", "#87cefa", "#2e8b57", "#0a74da", "#d62728", "#ff6347",
  "#e74c3c", "#c94c4c", "#ff4500", "#3cb371", "#2e8b57", "#66cdaa")

Pl_width <- 120
Pl_height <- 120
oma_mm <- c(20, 20, 0, 0)
mar_mm <- c(0, 0, 0, 0) + 2
n_row <- 1
n_col <- 1

RES = 600

# Specify the list of driving GCM models
driving_GCMs <- c("MPI", "EARTH", "HADGEM", "CNRM", "IPSL")

# Load datasets from CSV files and add an 'ENSEMBLE' column to each
CMIP5_data <- read_csv("../inputs/CMIP5.csv") |> 
  mutate(ENSEMBLE = "CMIP5")

CMIP6_data <- read_csv("../inputs/CMIP6.csv") |> 
  mutate(ENSEMBLE = "CMIP6")

RCMs_data <- read_csv("../inputs/EUR-44.csv") |> 
  mutate(ENSEMBLE = "EUR-44")

# Define constants for Clausius-Clapeyron equation parameters (Katul et al. 2012)
a <- 0.6111
b <- 17.5
c <- 240.97

# Allen et al. (1998)
# a <- 0.6108
# b <- 17.27
# c <- 237.3

# Psychrometric constant (considering the average elevation of the domain)
gamma <- 0.063

# Air density (kg m-3)
rhoAir <- 1.225

# Specific heat of air (J K-1kg-1),
CpAir <- 1004.67

# Atmospheric pressure (kPa)
AP <- 97.73785 # From ERA5-land for the domain in 1981–2005
# lambda <- (2.501 * (10^6) - 2361 *  8.033441) / (10^6)
# gamma <- CpAir*AP/(0.622*lambda*10^6)

# Based on the "SUPPLEMENT_DataTables_Meinshausen_6May2020.xlsx" for the CMIP6
CO2_2076_2100_RCP85_CMIP5 <- 828.38 # CMIP5, https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=download
CO2_2076_2100_RCP85_CMIP6 <- 979.42 # ppm, 972.28 ppm for whole world, 979.42 for northern hemisphere
CO2_1981_2005 <- 359.47

# dgs/gs predicted based on optimality theory
# 0.7–0.9
s <- 0.7

Kc25 <- 300
Ko25 <- 300
Coa <- 210
Ta_ST <- 25
Kc <- Kc25 * exp(0.074 * Ta_ST -25)
Ko <- Kc25 * exp(0.015 * Ta_ST -25)
a2 <- Kc * (1 + Coa / Ko)
a2
# a2 <- 510 
-1/((a2 / (s * CO2_1981_2005)) + 1) * (CO2_2076_2100_RCP85_CMIP6 - CO2_1981_2005) / CO2_1981_2005

# The correction by Yang et al. (2018) https://doi.org/10.1038/s41558-018-0361-0
S_rs_CO2 <- 0.09 # ppm^-1

rs_CO2 <- function(rs_1981_2005 = 70, S_rs_CO2 = 0.09, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP6){
  # S_rs_CO2 <- 0.09 % ppm^-1
  rs <- rs_1981_2005 * (1 + S_rs_CO2 / 100 * (CO2_2076_2100_RCP85 - CO2_1981_2005))
  return(rs)
}

# dgs/gs predicted by Yang et al. (2018)
(1 / rs_CO2(70) * 1000 - (1 / 70 * 1000)) / (1 / 70 * 1000)

# Inverse of the correction
rs_CO2_inv <- function(rs_2076_2100 = 109.0568, S_rs_CO2 = 0.09, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP6){
  # S_rs_CO2 <- 0.09 % ppm^-1
  rs <-  rs_2076_2100 / (1 + S_rs_CO2 / 100 * (CO2_2076_2100_RCP85 - CO2_1981_2005))
  return(rs)
}

# Surface resistance of the grass or alfalfa reference surface (s/m)
rs_ref_grass <- 70
rs_ref_alfalfa <- 35 # Original value in Allen et al. (2005) is 45 s/m, it was reduced in order to ensure that ET is always smaller than ETr
rs_ref_grass_2076_2100_RCP85_CMIP5 <- rs_CO2(rs_1981_2005 = rs_ref_grass, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP5)
rs_ref_grass_2076_2100_RCP85_CMIP6 <- rs_CO2(rs_1981_2005 = rs_ref_grass, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP6)
rs_ref_alfalfa_2076_2100_RCP85_CMIP5 <- rs_CO2(rs_ref_alfalfa, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP5)
rs_ref_alfalfa_2076_2100_RCP85_CMIP6 <- rs_CO2(rs_ref_alfalfa, CO2_2076_2100_RCP85 = CO2_2076_2100_RCP85_CMIP6)

# Gap-filling of variable using mean values for the purpose of the Penman-Monteith PET
GF_u <- TRUE
GF_RH <- TRUE

# Combine all datasets into a single dataframe
Data <- bind_rows(CMIP5_data, CMIP6_data, RCMs_data)

# Reorder columns to place key identifiers first
Data <- Data |> 
  select(MONTH, VAR, PERIOD, ENSEMBLE, everything()) |> 
  # Transform from wide to long format for easier manipulation
  pivot_longer(
    cols = -c(MONTH, VAR, PERIOD, ENSEMBLE), # Exclude these columns from pivoting
    names_to = "MODEL",                      # New column for model names
    values_to = "VALUE"                      # New column for corresponding values
  )

# Make ERA5 as independent ensemble
Data <- Data |>
  mutate(ENSEMBLE = if_else(MODEL == "ERA5", "ERA5", ENSEMBLE))

# Function for gapfilling of u and RH
gap_fill_tidy <- function(data, var, period) {
  data |>
    group_by(ENSEMBLE, MONTH) |>
    mutate(
      VALUE = if_else(
        VAR == var & PERIOD == period & is.na(VALUE),
        # mean(VALUE[VAR == var & PERIOD == period & MODEL != "ERA5"], na.rm = TRUE),
        mean(VALUE[VAR == var & PERIOD == period], na.rm = TRUE),
        VALUE
      )
    ) |>
    ungroup()
}

if (GF_u == TRUE) {
  Data <- gap_fill_tidy(Data, "u10", "1981_2005")
  Data <- gap_fill_tidy(Data, "u10", "2076_2100")
}

if (GF_RH == TRUE) {
  Data <- gap_fill_tidy(Data, "RH", "1981_2005")
  Data <- gap_fill_tidy(Data, "RH", "2076_2100")
}

# Further data transformations
Data <- Data |> 
  mutate(
    # Convert Roman numeral months to numeric
    MONTH = as.numeric(as.roman(MONTH)),
    # Add full month names based on numeric month
    MONTH_NAME = month.name[MONTH],
    # Assign the number of days in each month
    DAYS_IN_MONTH = case_when(
      MONTH == 1  ~ 31,
      MONTH == 2  ~ 28.25,  # Average considering leap years
      MONTH == 3  ~ 31,
      MONTH == 4  ~ 30,
      MONTH == 5  ~ 31,
      MONTH == 6  ~ 30,
      MONTH == 7  ~ 31,
      MONTH == 8  ~ 31,
      MONTH == 9  ~ 30,
      MONTH == 10 ~ 31,
      MONTH == 11 ~ 30,
      MONTH == 12 ~ 31,
      TRUE ~ NA_real_
    )
  ) |> 
  # Reorder columns to place new date-related columns appropriately
  select(MONTH, MONTH_NAME, DAYS_IN_MONTH, everything())

# Identify whether each model is a driving GCM
Data <- Data |> 
  mutate(
    Driving_GCM = if_else(MODEL %in% driving_GCMs, "yes", "no")
  ) |> 
  # Ensure 'VALUE' column is positioned at the end for readability
  relocate(VALUE, .after = last_col())

# Drop lines when VAR is NA
Data <- Data |> 
  filter(!is.na(VAR))

# Remove the duplicates
Data <- Data |> 
  distinct(MONTH, VAR, PERIOD, ENSEMBLE, MODEL, .keep_all = TRUE)

# Pivot data to a wider format to facilitate calculations across variables
Data_wide <- Data |> 
  pivot_wider(names_from = VAR, values_from = VALUE) |> 
  mutate(
    # Correct the sign for 'LW_net' (there is an unusual convention in the input data)
    LW_net = -LW_net,
    # Calculate Net Radiation
    Rn = SW_net + LW_net,
    # Compute Saturated Vapor Pressure using Clausius-Clapeyron equation
    e_sat = a * exp(b * Ta / (Ta + c)),
    # Actual vapor pressure
    e = (RH / 100 * e_sat),
    # Determine Vapor Pressure Deficit
    VPD = e_sat - e,
    # Slope of the saturated vapor pressure curve
    SVP = b * c * e_sat/(c + Ta)^2,
    # Latent heat of vaporization (MJ/kg)
    lambda = (2.501 * (10^6) - 2361 * Ta) / (10^6),
    # Unit conversions
    W_to_mm = 1 / lambda * 3600*24 / 10^6,
    # W_to_mm = 0.408 * 3600*24 / 10^6,
    # Evapotranspiration
    ET = W_to_mm * DAYS_IN_MONTH * LE,
    # Runoff
    RO = P - ET
  )

vpd_mean_by_ensemble <- Data_wide |>
  filter(PERIOD == "1981_2005") |>
  group_by(ENSEMBLE) |>
  summarise(mean_VPD = mean(VPD, na.rm = TRUE), .groups = "drop")

vpd_mean_by_ensemble
# 
# # Now drop ERA5
# Data <- Data |> 
#   filter(ENSEMBLE != "ERA5")
# 
# Data_wide <- Data_wide |> 
#   filter(ENSEMBLE != "ERA5")

# --- helper to fit slope-only lm(y ~ x - 1) per ENSEMBLE using lm ---
fit_slope_per_ensemble <- function(df, x, y, slope_name) {
  df |> 
    filter(!is.na(ENSEMBLE), !is.na(.data[[x]]), !is.na(.data[[y]])) |> 
    nest_by(ENSEMBLE) |> 
    mutate(!!slope_name := {
      m <- lm(reformulate(x, y, intercept = FALSE), data = data)
      unname(coef(m)[[x]])
    }) |> 
    ungroup() |> 
    select(ENSEMBLE, !!slope_name)
}

# Fit the slopes using linear regression
slopes_vpd   <- fit_slope_per_ensemble(Data_wide, "VPD",  "VPD_d",   "slope_vpd")
slopes_delta <- fit_slope_per_ensemble(Data_wide, "SVP", "delta_d", "slope_delta")

slopes_vpd <- slopes_vpd |>
  add_row(ENSEMBLE = "ERA5", slope_vpd = 1.16)
slopes_delta <- slopes_delta |>
  add_row(ENSEMBLE = "ERA5", slope_delta = 1.01)

# --- Join and fill both variables ---
Data_wide <- Data_wide |>
  left_join(slopes_vpd,  by = "ENSEMBLE") |>
  left_join(slopes_delta, by = "ENSEMBLE") |>
  mutate(
    VPD_d  = if_else(is.na(VPD_d)  & !is.na(VPD) & !is.na(slope_vpd),   slope_vpd   * VPD,  VPD_d),
    delta_d = if_else(is.na(delta_d) & !is.na(SVP) & !is.na(slope_delta), slope_delta * SVP, delta_d)
  ) |>
  select(-slope_vpd, -slope_delta)

# Continue with your downstream calculations
Data_wide <- Data_wide |> 
  mutate(
    # Equilibrium evaporation
    ET_eq = delta_d*(H + LE) / (delta_d + gamma) * W_to_mm * DAYS_IN_MONTH,
    # Priestley–Taylor potential evaporation
    PET_PT_alpha_1.26 = 1.26 * ET_eq,
    # Soil heat flux
    G = Rn - H - LE,
    # Available energy
    A = H + LE,
    # Specific humidity (g/kg)
    q = 0.622 * e / (AP - 0.378 * e) * 1000,
    # ETo FAO56 grass (with inline temporary vars)
    ETo_FAO56_grass = {
      zo = 0.123 * 0.12
      u2 = log(2 / zo) / log(10 / zo) * u10
      ra = 208 / u2
      rs = rs_ref_grass
      (delta_d * A + rhoAir * CpAir * VPD_d / ra) /
        (delta_d + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 grass with correction for CO2 effect on stomatal conductance for GCMs
    # Only when air temperature is above 5 deg. C
    ETo_FAO56_grass_GCM_CO2_corr = {
      zo = 0.123 * 0.12
      u2 = log(2 / zo) / log(10 / zo) * u10
      ra = 208 / u2
      rs = case_when(
        PERIOD == "2076_2100" & ENSEMBLE == "CMIP5" & Ta > 5 ~ rs_ref_grass_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ENSEMBLE == "CMIP6" & Ta > 5 ~ rs_ref_grass_2076_2100_RCP85_CMIP6,
        TRUE ~ rs_ref_grass       # <- default fallback value
      )
      (delta_d * A + rhoAir * CpAir * VPD_d / ra) /
        (delta_d + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 alfalfa (with inline temporary vars)
    ETo_FAO56_alfalfa = {
      zo = (0.123 * 0.12)
      u2 = log(2 / zo) / log(10 / zo) * u10
      # u2 = u10 * 4.87 / log(67.8 * 10 - 5.42)
      ra = 118 / u2
      rs = rs_ref_alfalfa
      (delta_d * A + rhoAir * CpAir * VPD_d / ra) /
        (delta_d + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 alfalfa with correction for CO2 effect on stomatal conductance for GCMs
    # Only when air temperature is above 5 deg. C
    ETo_FAO56_alfalfa_GCM_CO2_corr = {
      zo = (0.123*0.12)
      u2 = log(2 / zo) / log(10 / zo) * u10
      # u2 = u10 * 4.87 / log(67.8 * 10 - 5.42)
      ra = 118 / u2
      rs = case_when(
        PERIOD == "2076_2100" & ENSEMBLE == "CMIP5" & Ta > 5 ~ rs_ref_alfalfa_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ENSEMBLE == "CMIP6" & Ta > 5 ~ rs_ref_alfalfa_2076_2100_RCP85_CMIP6,
        TRUE ~ rs_ref_alfalfa       # <- default fallback value
      )
      (delta_d * A + rhoAir * CpAir * VPD_d / ra) /
        (delta_d + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    }
  )

# Drop variables which are not needed anymore
Data_wide <- Data_wide |> 
  select(-RH_d, -VPD_d, -delta_d)

# Calculate and assign the maximum ET/ET_eq ratio within each MODEL and ENSEMBLE
Data_wide <- Data_wide |> 
  group_by(MODEL, ENSEMBLE) |>
  mutate(
    max_ET_ratio = if (all(is.na(ET[MONTH %in% 4:8 & PERIOD == "1981_2005"])) ||
                       all(is.na(ET_eq[MONTH %in% 4:8 & PERIOD == "1981_2005"]))) {
      NA_real_
    } else {
      max(
        if_else(MONTH %in% 4:8 & PERIOD == "1981_2005" & !is.na(ET) & !is.na(ET_eq),
                ET / ET_eq, NA_real_),
        na.rm = TRUE
      )
    }
  ) |> 
  ungroup()

# # Remove ERA5 land 
# Data_wide <- Data_wide |> 
#   filter(MODEL != "ERA5")

# Create the boxplot
p <- ggplot(Data_wide |> filter(ENSEMBLE != "ERA5"), aes(x = ENSEMBLE, y = max_ET_ratio, fill = ENSEMBLE)) +
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2, coef = 3) +    # Bars width
  geom_boxplot(coef = 3) +
  scale_fill_manual(values = c(COL_CMIP5, COL_CMIP6, COL_RCMs)) +
  labs(y = bquote('Priestley-Taylor'~alpha), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

# Save the plot to a file
ggsave(filename = "../plots/ggplot2/Priestley-Taylor_alpha_ggplot2_TIDY.png", plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")

# Pivot data back to long format for subsequent analyses
Data_long <- Data_wide |> 
  pivot_longer(
    cols = -c(MONTH, MONTH_NAME, DAYS_IN_MONTH, PERIOD, ENSEMBLE, MODEL, Driving_GCM),
    names_to = "VAR",
    values_to = "VALUE",
    cols_vary = "slowest"  # Attempt to preserve original row order
  )

# Display unique variable names to verify transformations
Data_long |> 
  pull(VAR) |> 
  unique()

# Examine distinct models identified as driving GCMs
glimpse(
  Data_long |> 
    filter(Driving_GCM == "yes") |> 
    distinct(MODEL)
)

# Calculate annual statistics for each variable, period, ensemble, and model
annual_stats <- Data_long |> 
  group_by(VAR, PERIOD, ENSEMBLE, MODEL) |> 
  reframe(
    # "^" anchors the match to the beginning of the string
    annual_stat = if (grepl("^(ET|P|RO|PET|ETo)", unique(VAR))){
      sum(VALUE, na.rm = FALSE)  # Simple sum
    } else {
      sum(VALUE * DAYS_IN_MONTH, na.rm = FALSE) / sum(DAYS_IN_MONTH, na.rm = FALSE) # Weighted mean for other variables
    }
  )

# Verify unique variables in annual statistics
annual_stats |> 
  pull(VAR) |> 
  unique()

# Compute mean values across models for each variable, period, and ensemble
mean_values <- annual_stats |> 
  group_by(VAR, PERIOD, ENSEMBLE) |> 
  summarize(
    mean_value = mean(annual_stat, na.rm = TRUE),
    .groups = "drop"
  )

# Visualize mean values using a bar plot
p <- mean_values |>
  filter(ENSEMBLE != "ERA5") |> 
  ggplot(aes(x = PERIOD, y = mean_value, fill = ENSEMBLE)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ VAR, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "CMIP5" = COL_CMIP5,  # Color for CMIP5
      "CMIP6" = COL_CMIP6,  # Color for CMIP6
      "EUR-44" = COL_RCMs   # Color for EUR-44 (RCMs)
    )
  ) +
  labs(
    title = "Mean Values by Ensemble and Variable",
    y = "Mean Value",
    x = "Period"
  ) +
  theme_bw() + theme(panel.grid = element_blank())

# Save the plot to a file
ggsave(filename = "../plots/ggplot2/summary_ggplot2_TIDY.png", plot = p, width = Pl_width * 4, height = Pl_height * 2, dpi = RES, units = "mm")

################################################################################
# Prepare wide-format data to compute derived variables (rs_eff and g_eff)
# Remove irrelevant columns, ensure Ta is present, and calculate variables
annual_stats_wide <- annual_stats |> 
  pivot_wider(names_from = VAR, values_from = annual_stat) |> 
  select(-"W_to_mm", -"lambda") |> 
  drop_na(Ta) |>  # Ensure Ta is available for g_eff and rs_eff calculation
  mutate(
    rs_eff = rhoAir * CpAir / gamma * VPD / LE,
    g_eff = 1 / rs_eff * 1000,
    Bo = H / LE,
    # Integrated water vapor (mm)
    # Ruckstuhl et al. (2007) Observed relationship between surface specific humidity,
    # integrated water vapor, and longwave downward radiation at different altitude
    IWV = 2.67 * q + 0.68,
    # Precipitation efficiency (%)
    PE = P / (IWV * 365.25) * 100,
    # PE = P/e
  )

# Add aridity index and evaporative index and ET to PET ratio
annual_stats_wide <- annual_stats_wide |> 
  mutate(
    AI_FAO56_alfalfa = ETo_FAO56_alfalfa / P,
    AI_ETo_FAO56_alfalfa_GCM_CO2_corr = ETo_FAO56_alfalfa_GCM_CO2_corr / P,
    EI = ET / P,
    ET_over_ETo_FAO56_alfalfa = ETo_FAO56_alfalfa / ET,
    ET_over_ETo_FAO56_alfalfa_GCM_CO2_corr = ETo_FAO56_alfalfa_GCM_CO2_corr / ET
  )

# Convert derived variables back to long format and append to the main dataset
annual_stats <- bind_rows(
  annual_stats,
  annual_stats_wide |>
    select(PERIOD, ENSEMBLE, MODEL, rs_eff, g_eff, Bo, IWV, PE,
           AI_FAO56_alfalfa, AI_ETo_FAO56_alfalfa_GCM_CO2_corr,
           EI, ET_over_ETo_FAO56_alfalfa, ET_over_ETo_FAO56_alfalfa_GCM_CO2_corr) |>
    pivot_longer(cols = c(rs_eff, g_eff, Bo, IWV, PE,
                          AI_FAO56_alfalfa, AI_ETo_FAO56_alfalfa_GCM_CO2_corr,
                          EI, ET_over_ETo_FAO56_alfalfa, ET_over_ETo_FAO56_alfalfa_GCM_CO2_corr),
                 names_to = "VAR", values_to = "annual_stat")
)

# Fit a linear model
fit <- lm(ETo_FAO56_alfalfa ~ ETo_FAO56_grass, data = annual_stats_wide)
slope <- sprintf("%.2f", coef(fit)[2])
intercept <- sub("([+-])", "\\1 ", sprintf("%+.0f", coef(fit)[1]))
rsq <- sprintf("%.2f", summary(fit)$r.squared)

eq_label <- bquote(
  italic(y) == .(slope) * italic(x) ~ .(intercept) * ";" ~~
    R^2 == .(rsq)
)
p <- ggplot(annual_stats_wide, aes(x = ETo_FAO56_grass, y = ETo_FAO56_alfalfa)) +
  geom_point(color = "#1c91df", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = as.expression(eq_label),
           hjust = 1.1, vjust = -0.5, parse = TRUE, size = 4, color = "black") +
  labs(
    x = bquote("ETo"["FAO56 grass"]~"(mm yr"^"-1"*")"),
    y = bquote("ETo"["FAO56 alfalfa"]~"(mm yr"^"-1"*")"),
    title = "FAO56 vs VPD-adjusted ETo"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Save the plot to a file
ggsave(filename = "../plots/ggplot2/FAO_ET_ggplot2_TIDY.png", plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")

annual_stats_wide |> mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> pull(ET_diff) |> plot()
annual_stats_wide |> mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> pull(ET_diff) |> min()
annual_stats_wide |> mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> filter(ET_diff < 0)

annual_stats_wide |> filter(ENSEMBLE == "ERA5") |> View()
annual_stats_wide |>
  filter(ENSEMBLE == "EUR-44", PERIOD == "1981_2005") |>
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE))
  ) |>  View()

annual_stats_wide |> filter(ENSEMBLE == "ERA5") |> select(VPD, EI, AI_FAO56_alfalfa) |> View()

# Calculate the normalized difference between periods
norm_diff <- annual_stats |>
  group_by(VAR, ENSEMBLE, MODEL) |>
  reframe(
    d_annual_stat = if (VAR[1] == "Bo") {
      (annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"]) /
        (1 + annual_stat[PERIOD == "1981_2005"])
    } else {
      (annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"]) /
        annual_stat[PERIOD == "1981_2005"]
    },
    VAR = if (VAR[1] == "Bo") "d_Bo_adj" else paste0("d_", VAR[1], "_over_", VAR[1])
  )

# Calculate the difference between periods
diff <- annual_stats |> 
  group_by(VAR, ENSEMBLE, MODEL) |> 
  reframe(
    d_annual_stat = annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"],
    VAR = paste0("d_", VAR[1])
  )

# Combine the two data frames
all_diff <- bind_rows(norm_diff, diff)

# Reshape data: pivot wider to have P and ET in separate columns
Data_to_plot <- list()
Data_to_plot[["diffs"]] <- all_diff |> 
  pivot_wider(names_from = VAR, values_from = d_annual_stat) |> 
  mutate(
    # Assign colors based on ENSEMBLE
    color = case_when(
      ENSEMBLE == "CMIP5" ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" ~ COL_CMIP6,
      ENSEMBLE == "EUR-44" ~ COL_RCMs,
      TRUE ~ NA_character_
    ),
    # Assign fill colors based on conditions
    fill = case_when(
      ENSEMBLE == "CMIP5" & MODEL %in% driving_GCMs ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" & MODEL %in% driving_GCMs ~ COL_CMIP6,
      ENSEMBLE == "EUR-44" ~ COL_RCMs,
      TRUE ~ NA_character_
    ),
    # Assign border colors based on conditions
    border = case_when(
      (ENSEMBLE %in% c("CMIP5", "CMIP6") & MODEL %in% driving_GCMs) | ENSEMBLE == "EUR-44" ~ BORDER_COLOR,
      ENSEMBLE == "CMIP5" ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" ~ COL_CMIP6,
      TRUE ~ NA_character_
    ),
    # Set shape and linetype
    shape = 21,
    linetype = "solid"
  ) |>
  filter(ENSEMBLE != "ERA5") |> 
  rename(ensemble = ENSEMBLE, model = MODEL) |> 
  mutate(model = gsub("_", "-", model))

# --------------
# Start plotting

#----------------
# dP/P vs. dET/ET
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.88, 0.83, 0.78),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

LM_eq_labels <- tibble(
  x = c(0.7, 0.7, 0.7),
  y = c(0.14, 0.09, 0.04)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_ET_over_ET, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_ET_over_ET", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'ET/'*'ET'),  y_lab = bquote(Delta*'P/'*'P'),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.87, 0.95), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels, plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_P_over_P_vs_delta_ET_over_ET_ggplot2_TIDY.png", save_ggplot2_obj_as="p1")

#----------------
# dRO/RO vs. dET/ET
LM_eq_labels <- data.frame(
  x = c(0.27, 0.27, 0.27),
  y = c(0.22, 0.17, 0.12)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_ET_over_ET, d_RO_over_RO, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_ET_over_ET", y = "d_RO_over_RO"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'ET/'*'ET'),  y_lab = bquote(Delta*'R'[o]*'/'*'R'[o]),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.80, 0.96), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_RO_over_RO_vs_delta_ET_over_ET_ggplot2_TIDY.png", save_ggplot2_obj_as="p2")

#----------------
# dET/ET vs. dP/P
LM_eq_labels <- tibble(
  x = c(0.7, 0.7, 0.7),
  y = c(0.14, 0.09, 0.04)
)

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.95, 0.90, 0.85),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_P_over_P, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_P_over_P", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.79), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels, plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_ET_over_ET_vs_delta_P_over_P_ggplot2_TIDY.png", save_ggplot2_obj_as="p3")

#----------------
# dET/ET vs. dRO/RO
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_P_over_P, d_RO_over_RO, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_P_over_P", y = "d_RO_over_RO"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'R'[o]*'/'*'R'[o]),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.75, 0.97), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_RO_over_RO_vs_delta_P_over_P_ggplot2_TIDY.png", save_ggplot2_obj_as="p4")


#---------------
#  Combined plot

# Define colors
cols <- c(CMIP5 = COL_CMIP5, CMIP6 = COL_CMIP6, `EUR-44` = COL_RCMs)

#---- ET and RO ratios for p1 ----
et_ratios <- annual_stats |>
  filter(PERIOD == "1981_2005", VAR %in% c("ET", "RO")) |>
  pivot_wider(names_from = VAR, values_from = annual_stat) |>
  mutate(
    et_ratio = ET / (ET + RO),
    ro_ratio = RO / (ET + RO)
  ) |>
  group_by(ENSEMBLE) |>
  summarise(
    et = sprintf("%.2f", mean(et_ratio, na.rm = TRUE)),
    ro = sprintf("%.2f", mean(ro_ratio, na.rm = TRUE)),
    .groups = "drop"
  )

# Add annotations to p1
for (i in seq_len(nrow(et_ratios))) {
  ens <- et_ratios$ENSEMBLE[i]
  et_val <- et_ratios$et[i]
  ro_val <- et_ratios$ro[i]
  
  p1 <- p1 +
    annotation_custom(
      textGrob(
        bquote('ET/(ET + R'[o]*') = '~.(et_val)),
        x = unit(0.27, "npc"), y = unit(0.93 - 0.05 * (i - 1), "npc"),
        hjust = 0, gp = gpar(fontsize = 8, col = cols[ens])
      )
    ) +
    annotation_custom(
      textGrob(
        bquote('R'[o]*'/(ET + R'[o]*') = '~.(ro_val)),
        x = unit(0.27, "npc"), y = unit(0.76 - 0.05 * (i - 1), "npc"),
        hjust = 0, gp = gpar(fontsize = 8, col = cols[ens])
      )
    )
}

#---- P and RO ratios for p3 ----
p_ratios <- annual_stats |>
  filter(PERIOD == "1981_2005", VAR %in% c("P", "RO")) |>
  pivot_wider(names_from = VAR, values_from = annual_stat) |>
  mutate(
    p_ratio = P / (P - RO),
    ro_ratio = RO / (P - RO)
  ) |>
  group_by(ENSEMBLE) |>
  summarise(
    p = sprintf("%.2f", mean(p_ratio, na.rm = TRUE)),
    ro = sprintf("%.2f", mean(ro_ratio, na.rm = TRUE)),
    .groups = "drop"
  )

# Add annotations to p3
for (i in seq_len(nrow(p_ratios))) {
  ens <- p_ratios$ENSEMBLE[i]
  p_val <- p_ratios$p[i]
  ro_val <- p_ratios$ro[i]
  
  p3 <- p3 +
    annotation_custom(
      textGrob(
        bquote('P/(P - R'[o]*') = '~.(p_val)),
        x = unit(0.7, "npc"), y = unit(0.55 - 0.05 * (i - 1), "npc"),
        hjust = 0, gp = gpar(fontsize = 8, col = cols[ens])
      )
    ) +
    annotation_custom(
      textGrob(
        bquote('R'[o]*'/(P - R'[o]*') = '~.(ro_val)),
        x = unit(0.7, "npc"), y = unit(0.38 - 0.05 * (i - 1), "npc"),
        hjust = 0, gp = gpar(fontsize = 8, col = cols[ens])
      )
    )
}

#---- Subplot labels
plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
labels <- c("a", "b", "c", "d")

plots <- Map(function(plot, label) {
  plot + annotation_custom(
    textGrob(label, x = unit(0.05, "npc"), y = unit(0.96, "npc"),
             gp = gpar(fontsize = 12))
  )
}, plots, labels)

p1 <- plots$p1; p2 <- plots$p2; p3 <- plots$p3; p4 <- plots$p4

#---- Combine and save
panel_figure <- (p1 | p2) / (p3 | p4)

ggsave('../plots/ggplot2/panel_fig_hydrological_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#-------------------------
# Delta_T versus delta_P/P
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

# Manual adjustment
Y_range_man <- c(-0.12, 0.32)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p1")

plot_data <- ggplot_build(p1)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

# Test of Clausius-Clapeyron
# Central Europe with Ta 8 °C and globe with Ta 15 °C
(-b*8/(c+8)^2+b/(c+8))
(-b*15/(c+15)^2+b/(c+15))

b*c/(c+8)^2

b/(c+8)
b/(c+15)

-b*8/(c+8)^2
-b*15/(c+15)^2

b/(c+8)-b*8/(c+8)^2
b/(c+15)-b*15/(c+15)^2

# Clausius-Clapeyron
X <- seq(0, 10, 0.1)

# Mean air temperature form CMIP5, CMIP6 and RCMs
T_avg <- annual_stats |> 
  filter(VAR == "Ta", PERIOD == "1981_2005", ENSEMBLE %in% c("CMIP5", "CMIP6", "EUR-44")) |> 
  group_by(ENSEMBLE) |> 
  summarise(T_avg = mean(annual_stat, na.rm = TRUE)) |>  
  pull(T_avg) |> 
  mean()  # Overall mean across the 3 ensemble means

rand_int_I <- 0
d_P_over_P_CC_linear <- b/(c+T_avg)*X

d_P_over_P_CC_nonlinear <- (-b*T_avg/(c+T_avg)^2+b/(c+T_avg))*X

# Create data frames for the lines
DF1_CC_linear <- data.frame(X = X, Y = rand_int_I + d_P_over_P_CC_linear)
DF1_CC_nonlinear <- data.frame(X = X, Y = rand_int_I + d_P_over_P_CC_nonlinear)

# Normalize the X and Y values
DF1_CC_linear$X_norm <- (DF1_CC_linear$X - X_range[1]) / (X_range[2] - X_range[1])
DF1_CC_linear$Y_norm <- (DF1_CC_linear$Y - Y_range[1]) / (Y_range[2] - Y_range[1])

# Calculate the differences (deltas) for the normalized values
dX_norm <- diff(DF1_CC_linear$X_norm)
dY_norm <- diff(DF1_CC_linear$Y_norm)

# Calculate the angles in radians
angles_radians <- atan2(dY_norm, dX_norm)

# Convert the angles to degrees
angle_degrees <- mean(angles_radians * (180 / pi))

p1 <- p1 + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.02, X_range), y = rescale(0.77, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

# Save the plot
ggsave('../plots/ggplot2/delta_P_over_P_vs_Ta_ggplot2_TIDY.png', plot = p1, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Delta_T versus delta_ET/ET

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p2")

p2 <- p2 + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.02, X_range), y = rescale(0.77, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

# Save the plot
ggsave('../plots/ggplot2/delta_ET_over_ET_vs_Ta_ggplot2_TIDY.png', plot = p2, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Delta_T versus d_VPD/VPD

Plot_labels <- tibble(
  x = c(0.04, 0.04, 0.04),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_VPD_over_VPD, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_VPD_over_VPD"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'VPD/'*'VPD'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p3")

plot_data <- ggplot_build(p3)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

rand_int_I <- 0
CC_linear <- b/(c+T_avg)*X
CC_nonlinear <- (-b*T_avg/(c+T_avg)^2+b/(c+T_avg))*X

# Create data frames for the lines
DF1_CC_linear <- data.frame(X = X, Y = rand_int_I + CC_linear)
DF1_CC_nonlinear <- data.frame(X = X, Y = rand_int_I + CC_nonlinear)

# Normalize the X and Y values
DF1_CC_linear$X_norm <- (DF1_CC_linear$X - X_range[1]) / (X_range[2] - X_range[1])
DF1_CC_linear$Y_norm <- (DF1_CC_linear$Y - Y_range[1]) / (Y_range[2] - Y_range[1])

# Calculate the differences (deltas) for the normalized values
dX_norm <- diff(DF1_CC_linear$X_norm)
dY_norm <- diff(DF1_CC_linear$Y_norm)

# Calculate the angles in radians
angles_radians <- atan2(dY_norm, dX_norm)

# Convert the angles to degrees
angle_degrees <- mean(angles_radians * (180 / pi))


p3 <- p3 + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.6, X_range), y = rescale(0.28, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

# Save the plot
ggsave('../plots/ggplot2/delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p3, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Delta_T versus d_geff/geff

# Account also for the interactive term
# (d_ET_over_ET_CMIP5_a - d_VPD_over_VPD_CMIP5_a)/(1 + d_VPD_over_VPD_CMIP5_a)
# (d_ET_over_ET_CMIP6_a - d_VPD_over_VPD_CMIP6_a)/(1 + d_VPD_over_VPD_CMIP6_a)
# (d_ET_over_ET_RCMs_a - d_VPD_over_VPD_RCMs_a)/(1 + d_VPD_over_VPD_RCMs_a)

# --- Compute CO2-corrected d_g_eff for EUR-44 ---

g_eff_corr_RCMs <- annual_stats_wide |>
  filter(PERIOD %in% c("1981_2005", "2076_2100"), ENSEMBLE == "EUR-44") |>
  select(PERIOD, MODEL, rs_eff, g_eff) |>
  pivot_wider(names_from = PERIOD, values_from = c(rs_eff, g_eff)) |>
  mutate(
    g_eff_corr = (1 / (rs_CO2(rs_eff_2076_2100)) ) * 1000,
    d_g_eff_CO2_corr = (g_eff_corr - g_eff_1981_2005) / g_eff_1981_2005
  ) |>
  select(MODEL, d_g_eff_CO2_corr) |> 
  mutate(MODEL = gsub("_", "-", MODEL))

# --- Prepare plotting data ---

# Original CMIP5, CMIP6, and EUR-44 points
points_base <- Data_to_plot$diffs |>
  select(d_Ta, d_g_eff_over_g_eff, ensemble, color, fill, border, shape, model, linetype) |>
  rename(x = d_Ta, y = d_g_eff_over_g_eff)

# CO2-corrected EUR-44 points
points_corr <- Data_to_plot$diffs |>
  filter(ensemble == "EUR-44") |>
  left_join(g_eff_corr_RCMs, by = c("model" = "MODEL")) |>
  mutate(
    ensemble = "EUR-44_CO2_corr",
    x = d_Ta,
    y = d_g_eff_CO2_corr,
    fill = "#2e8b57",
    color = "#2e8b57",
    linetype = "dashed"
  ) |>
  select(x, y, ensemble, color, fill, border, shape, model, linetype)

# Combine all points
points_all <- bind_rows(points_base, points_corr)

make_scatter_plot(data = points_all,
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'g'['eff']~'/'~'g'['eff']),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p4")

plot_data <- ggplot_build(p4)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

p4 <- p4 + annotate("text", x = rescale(0.66, X_range), y = rescale(0.96, Y_range),
                    label = as.expression(bquote(bold("EUR-44 CO"["2"]~"corr"))), hjust = 0, color = "#2e8b57", size=4)

# Save the plot
ggsave('../plots/ggplot2/delta_g_eff_over_g_eff_vs_Ta_ggplot2_TIDY.png', plot = p4, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# #---- Subplot labels
# plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
# labels <- c("a", "b", "c", "d")
# 
# plots <- Map(function(plot, label) {
#   plot + annotation_custom(
#     textGrob(label, x = unit(0.05, "npc"), y = unit(0.96, "npc"),
#              gp = gpar(fontsize = 12))
#   )
# }, plots, labels)
# 
# p1 <- plots$p1; p2 <- plots$p2; p3 <- plots$p3; p4 <- plots$p4
# 
# #---- Combine and save
# panel_figure <- (p1 | p2) / (p3 | p4)
# 
# ggsave('../plots/ggplot2/panel_fig_Ta_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

# notes: make the labels to be read from csv file, avoid overlapping with lines, allow to not show labels for selected data

#-------------------------------------------------------------------------------

#---------------------------------------
# Effective conductance for double check
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

Data_to_plot$diffs <- Data_to_plot$diffs |> 
  mutate(d_g_eff_over_g_eff_CHECK = (d_Rn_over_Rn - d_Bo_adj - d_VPD_over_VPD) / (1 + d_VPD_over_VPD)
  )

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_g_eff_over_g_eff_CHECK, d_g_eff_over_g_eff, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_g_eff_over_g_eff_CHECK", y = "d_g_eff_over_g_eff"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote('('*Delta*'R'['n']*'/'*'R'['n']~'-'~Delta*'B'['o']*'/(1 + '*'B'['o']*')'~'-'~Delta*'VPD'*'/VPD'*') / ('*1 + ~Delta*'VPD'*'/VPD)'),
                  y_lab = bquote(Delta*'g'['eff']~'/'~'g'['eff']),
                  hline = FALSE, vline = FALSE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.93, 0.84), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_g_eff_star_over_g_eff_star_check_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

# Reshape data: pivot wider to have P and ET in separate columns
Data_to_plot[["abs"]] <- annual_stats |> 
  pivot_wider(names_from = VAR, values_from = annual_stat) |> 
  mutate(
    # Assign colors based on ENSEMBLE
    color = case_when(
      ENSEMBLE == "CMIP5" ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" ~ COL_CMIP6,
      ENSEMBLE == "EUR-44" ~ COL_RCMs,
      TRUE ~ NA_character_
    ),
    # Assign fill colors based on conditions
    fill = case_when(
      ENSEMBLE == "CMIP5" & MODEL %in% driving_GCMs ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" & MODEL %in% driving_GCMs ~ COL_CMIP6,
      ENSEMBLE == "EUR-44" ~ COL_RCMs,
      TRUE ~ NA_character_
    ),
    # Assign border colors based on conditions
    border = case_when(
      (ENSEMBLE %in% c("CMIP5", "CMIP6") & MODEL %in% driving_GCMs) | ENSEMBLE == "EUR-44" ~ BORDER_COLOR,
      ENSEMBLE == "CMIP5" ~ COL_CMIP5,
      ENSEMBLE == "CMIP6" ~ COL_CMIP6,
      TRUE ~ NA_character_
    ),
    # Set shape and linetype
    shape = case_when(
      PERIOD == "1981_2005" ~ 21,
      PERIOD == "2076_2100" ~ 24,
      TRUE ~ NA_real_
    ),
    linetype = "solid"
  ) |>
  filter(ENSEMBLE != "ERA5") |> 
  rename(ensemble = ENSEMBLE, model = MODEL) |> 
  mutate(model = gsub("_", "-", model))

attributes(Data_to_plot$abs)

attr(Data_to_plot$abs, "description") <- "This data frame was created to plot not only the differences but also the actual values within the two examined periods."
#---------------------------
# VPD versus g_eff

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD, y = g_eff, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  Y_range_man = c(2, 17.2),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote("g"["eff"]~"(mm s"^"-1"*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

add_manual_legend <- function(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70){
  plot_data <- ggplot_build(p)
  X_range <- plot_data$layout$panel_params[[1]]$x.range
  Y_range <- plot_data$layout$panel_params[[1]]$y.range
  
  # Manually add two points and their labels to the top right
  p <- p +
    annotate("point", x = rescale(x1, X_range), y = rescale(y1, Y_range), shape = 21, size = 2, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y1, Y_range), label = "1981–2005", size = 3, hjust = 0)
  
  p <- p +
    annotate("point", x = rescale(x1, X_range), y = rescale(y2, Y_range), shape = 24, size = 2, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y2, Y_range), label = "2076–2100 RCP 8.5", size = 3, hjust = 0)
  
  return(p)
}

p <- add_manual_legend(p, x1 = 0.06, x2 = 0.1, y1 = 0.1, y2 = 0.05)

ggsave('../plots/ggplot2/g_eff_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# Remove layers that are geom_smooth (fit lines and ribbons)
p$layers <- lapply(p$layers, function(layer) {
  if (inherits(layer$geom, c("GeomSmooth", "GeomRibbon"))) {
    NULL
  } else {
    layer
  }
})

# Drop NULLs
p$layers <- Filter(Negate(is.null), p$layers)

# Now plot without the smooth fit lines or confidence belts
ggsave('../plots/ggplot2/g_eff_versus_VPD_ggplot2_no_fit_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# Normalized geff vs. VPD
LM_eq_labels <- tibble(
  x = c(0.53, 0.53, 0.53),
  y = c(0.94, 0.86, 0.78)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_g_eff_over_g_eff, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_g_eff_over_g_eff"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*g[eff]~"/"~g[eff]),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  save_ggplot2_obj_as="p_geff_VPD_norm")

p_geff_VPD_norm <- p_geff_VPD_norm +
  theme(
    panel.background = element_rect(fill = NA, colour = NA),  # panel area transparent
    plot.background  = element_rect(fill = NA, colour = NA)   # outer background transparent
  )


# Shrink text size in geom_text_repel and annotate("text")
p_geff_VPD_norm <- resize_plot_elements(
  p_geff_VPD_norm,
  repel_size = 1,
  text_size  = 3
)

# Combining plots
combined_plot <- p + 
  inset_element(
    p_geff_VPD_norm,
    left   = 0.3,  # relative x-position (0–1)
    bottom = 0.3,  # relative y-position (0–1)
    right  = 1.0,  # relative width (0–1)
    top    = 1.0   # relative height (0–1)
  )

# Save the plot
ggsave('../plots/ggplot2/g_eff_versus_VPD_ggplot2_no_fit_and_delta_geff_over_geff_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = combined_plot, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for GCMs

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & ensemble %in% c("CMIP5", "CMIP6"),
        1 / (rs_CO2_inv(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff,
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for RCMs

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & ensemble %in% c("EUR-44"),
        1 / (rs_CO2(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff,
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------------------------------------
# VPD versus g_eff normalized by available energy

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD, y = g_eff / (H + LE), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_AE_norm_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#----------------------------------------------------------------------
# VPD versus g_eff GCM CO2 corrected and normalized by available energy

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        1 / (rs_CO2_inv(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (H + LE),
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_AE_norm_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff RCM CO2 corrected and normalized by available energy

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        1 / (rs_CO2(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (H + LE),
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_AE_norm_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------------------------------------
# VPD versus g_eff normalized by global radiation

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD, y = g_eff / (Rg), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_Rg_norm_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff GCM CO2 corrected and normalized by global radiation

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        1 / (rs_CO2_inv(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (Rg),
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_Rg_norm_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff RCM CO2 corrected and normalized by global radiation

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        1 / (rs_CO2(rs_eff)) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (Rg),
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/g_eff_Rg_norm_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------
# VPD versus rs_eff

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD, y = rs_eff, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('r'['s eff']~(s/m)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/rs_eff_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for GCMs

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      rs_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        rs_CO2_inv(rs_eff),
        rs_eff
      ),
      x = VPD,
      y = rs_eff,
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/rs_eff_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for RCMs
make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(
      rs_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        rs_CO2(rs_eff),
        rs_eff
      ),
      x = VPD,
      y = rs_eff,
      ensemble = interaction(ensemble, PERIOD, drop = TRUE)
    ) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/rs_eff_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

# Assessing whether aerodynamics or thermodynamics dominate the relative change in ET

make_scatter_plot(data = Data_to_plot$diffs |>
                    mutate(OMEGA = d_g_eff_over_g_eff / d_VPD_over_VPD) |>
                    select(d_Ta, OMEGA, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = d_Ta, y = OMEGA),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'g'['eff']~'/'~'g'['eff']~'/'~'('*Delta*'VPD'~'/'~'VPD)'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_path = "../plots/ggplot2", plot_name = "Omega_versus_delta_Ta_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

#-------------------------------------------------
# Delta_Rn and Bowen ratio term versus delta_ET/ET
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    mutate(d_A_minus_Bo_adj = d_A_over_A - d_Bo_adj) |>
                    select(d_A_minus_Bo_adj, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = d_A_minus_Bo_adj, y = d_ET_over_ET),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G) - Delta * B[o] / (1 + B[o])),
                  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.89), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_Rn&delta_Bo_versus_delta_ET_over_ET_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

#---------------------------
# Delta_T versus d_e/e
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_e_over_e, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_e_over_e"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'e'~'/'~'e'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p")

plot_data <- ggplot_build(p)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

rand_int_I <- 0
CC_linear <- b/(c+T_avg)*X
CC_nonlinear <- (-b*T_avg/(c+T_avg)^2+b/(c+T_avg))*X

# Create data frames for the lines
DF1_CC_linear <- data.frame(X = X, Y = rand_int_I + CC_linear)
DF1_CC_nonlinear <- data.frame(X = X, Y = rand_int_I + CC_nonlinear)

# Normalize the X and Y values
DF1_CC_linear$X_norm <- (DF1_CC_linear$X - X_range[1]) / (X_range[2] - X_range[1])
DF1_CC_linear$Y_norm <- (DF1_CC_linear$Y - Y_range[1]) / (Y_range[2] - Y_range[1])

# Calculate the differences (deltas) for the normalized values
dX_norm <- diff(DF1_CC_linear$X_norm)
dY_norm <- diff(DF1_CC_linear$Y_norm)

# Calculate the angles in radians
angles_radians <- atan2(dY_norm, dX_norm)

# Convert the angles to degrees
angle_degrees <- mean(angles_radians * (180 / pi))

p <- p + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.7, X_range), y = rescale(0.7, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

ggsave('../plots/ggplot2/delta_e_over_e_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#-----------------------
# Delta_T versus d_PE/PE
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_PE_over_PE, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_PE_over_PE"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'PE'~'/'~'PE'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_PE_over_PE_vs_Ta_ggplot2_TIDY.png", save_ggplot2_obj_as="p")

annual_stats_wide |>
  reframe(PE_mean = mean(PE, na.rm = TRUE), .by = c(ENSEMBLE, PERIOD))

#-------------------------------------------------------------------------------

#-----------------------------
# Delta Rg versus delta SW_net
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.88, 0.83, 0.78),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Rg_over_Rg, d_SW_net_over_SW_net, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_SW_net_over_SW_net"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'SW'['net']*'/'*'SW'['net']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.94, 0.92), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p5")

# Save the plot
ggsave('../plots/ggplot2/delta_Rg_versus_delta_SW_net_ggplot2_TIDY.png', plot = p5, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#-----------------------------
# Delta Rg versus delta LW_net

# Note that the mean LW_net is a negative number
# CMIP5 has slightly more negative values in the future
# RCMs has less negative values in the future
# Because the denominator is negative, the d_LW_net_over_LW_net_RCMs_a with delta being positive becomes a negative number

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Rg_over_Rg, d_LW_net_over_LW_net, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_LW_net_over_LW_net"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'LW'['net']*'/'*'LW'['net']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.02, 0.25), robust_regression = TRUE,
                  save_ggplot2_obj_as="p6")

# Save the plot
ggsave('../plots/ggplot2/delta_Rg_versus_delta_LW_net_ggplot2_TIDY.png', plot = p6, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------
# Delta Rg versus delta Rn
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Rg_over_Rg, d_Rn_over_Rn, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_Rn_over_Rn"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'R'['n']*'/'*'R'['n']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.94, 0.75), robust_regression = TRUE,
                  save_ggplot2_obj_as="p7")

# Save the plot
ggsave('../plots/ggplot2/delta_Rg_versus_delta_Rn_ggplot2_TIDY.png', plot = p7, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------
# Delta_Rn versus delta_ET/ET

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_A_over_A, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_A_over_A", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G)),  y_lab = bquote(Delta * ET / ET %~~% Delta * (R[n] - G) / (R[n] - G) - Delta * B[o] / (1 + B[o])),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.91), robust_regression = TRUE,
                  save_ggplot2_obj_as="p8")

# Save the plot
ggsave('../plots/ggplot2/delta_Rn_versus_delta_ET_over_ET_ggplot2_TIDY.png', plot = p8, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---- Subplot labels
plots <- list(p5 = p5, p6 = p6, p7 = p7, p8 = p8)
labels <- c("a", "b", "c", "d")

plots <- Map(function(plot, label) {
  x_pos <- if (label == "d") 0.08 else 0.05  # Move "d" to the right
  plot + annotation_custom(
    textGrob(label, x = unit(x_pos, "npc"), y = unit(0.96, "npc"),
             gp = gpar(fontsize = 12))
  )
}, plots, labels)

p5 <- plots$p5; p6 <- plots$p6; p7 <- plots$p7; p8 <- plots$p8

#---- Combine and save
panel_figure <- (p5 | p6) / (p7 | p8)

#---- Combine and save
panel_figure <- (p5 | p6) / (p7 | p8)

ggsave('../plots/ggplot2/panel_fig_Rg&Rn_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#------------------------------------
# Delta_Rn and delta Bowen ratio term

Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_A_over_A, d_Bo_adj, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_A_over_A", y = "d_Bo_adj"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G)),
                  y_lab = bquote(Delta * B[o] / (1 + B[o])),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.59, 0.96), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../plots/ggplot2", plot_name = "delta_Rn_versus_delta_Bo_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

# Focus on complementary relation in Budyko space

# Budyko curves

#-------------------------------------------------------------------------------
# Fitting Budyko curve
Data_to_plot$abs <- Data_to_plot$abs |>
  mutate(PET = ETo_FAO56_alfalfa) |>
  rowwise() |>
  mutate(
    n = if (all(!is.na(c(PET, ET, P)))) Budyko_curve_optim(PET, ET, P) else NA_real_,
    ω = if (!is.na(n)) n + 0.72 else NA_real_
  ) |>
  ungroup() |>
  relocate(PET, n, ω, .before = color)

# Create the boxplot
# p <- ggplot(Data_to_plot$abs, aes(x = PERIOD, y = ω, fill = ensemble)) +
#   stat_boxplot(geom = "errorbar", width = 0.2, coef = 3,
#                position = position_dodge(width = 0.8)) +
#   geom_boxplot(coef = 3, position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = c(COL_CMIP5, COL_CMIP6, COL_RCMs)) +
#   scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
#   labs(y = "ω", x = "Period") +
#   theme_bw() +
#   theme(
#     panel.grid = element_blank(),
#     legend.position = c(0.95, 0.95),   # <- top-right corner
#     legend.justification = c("right", "top"), # anchor legend box
#     legend.background = element_rect(fill = "white", colour = "black") # optional box
#   )


p_omega_box <- ggplot(Data_to_plot$abs, aes(x = PERIOD, y = ω, fill = ensemble)) +
  stat_boxplot(geom = "errorbar", width = 0.2, coef = 3,
               position = position_dodge(width = 0.8)) +
  geom_boxplot(coef = 3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(COL_CMIP5, COL_CMIP6, COL_RCMs)) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  # relabel the x axis
  scale_x_discrete(labels = c(`1981_2005` = "1981–2005",
                              `2076_2100` = "2076–2100")) +
  labs(y = "ω", x = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"   # remove legend
  )

# Save the plot to a file
ggsave(filename = "../plots/ggplot2/ω_alpha_ggplot2_TIDY.png", plot = p_omega_box, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")

# Make a simple plot with linear regression and equation
simple_scatter_plot <- function(data, x, y, xpos = Inf, ypos = Inf,
                                hjust = 1.1, vjust = 2) {
  x_expr <- rlang::enquo(x)
  y_expr <- rlang::enquo(y)
  
  # fit the model
  fit <- lm(rlang::eval_tidy(y_expr, data) ~ rlang::eval_tidy(x_expr, data))
  
  # build label
  label <- paste0(
    "y = ", round(coef(fit)[2], 3), "x + ", round(coef(fit)[1], 3),
    "   R² = ", round(summary(fit)$r.squared, 3)
  )
  
  # plot
  ggplot(data, aes(x = !!x_expr, y = !!y_expr)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = xpos, y = ypos,
             hjust = hjust, vjust = vjust, label = label) +
    theme_bw()
}

# Omega and VPD
simple_scatter_plot(Data_to_plot$abs, VPD, ω)

# Omega and P - ET
simple_scatter_plot(Data_to_plot$abs, P - ET, ω)

# Omega and PET / P
simple_scatter_plot(Data_to_plot$abs, PET / P, ω)

# Omega and LAI
simple_scatter_plot(Data_to_plot$abs, LAI, ω)

# Omega and g_eff
simple_scatter_plot(Data_to_plot$abs, g_eff, ω,
                    xpos = -Inf, ypos = Inf, hjust = -0.1, vjust = 2)

# Omega^-1 and g_eff
simple_scatter_plot(Data_to_plot$abs, g_eff, 1 / ω)

# Omega and ET / PET
simple_scatter_plot(Data_to_plot$abs, ET / PET, ω,
                    xpos = -Inf, ypos = Inf, hjust = -0.1, vjust = 2)

# ET / PET and VPD
simple_scatter_plot(Data_to_plot$abs, ET / PET, VPD)


#######################
# Omega and ET / PET

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = ET / PET, y = n, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("ET / PET"),  
  y_lab = bquote("ω"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

add_manual_legend <- function(p, x1 = 0.05, x2 = 0.09, y1 = 0.75, y2 = 0.70){
  plot_data <- ggplot_build(p)
  X_range <- plot_data$layout$panel_params[[1]]$x.range
  Y_range <- plot_data$layout$panel_params[[1]]$y.range
  
  # Manually add two points and their labels to the top right
  p <- p +
    annotate("point", x = rescale(x1, X_range), y = rescale(y1, Y_range), shape = 21, size = 2, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y1, Y_range), label = "History", size = 3, hjust = 0)
  
  p <- p +
    annotate("point", x = rescale(x1, X_range), y = rescale(y2, Y_range), shape = 24, size = 2, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y2, Y_range), label = "Scenario", size = 3, hjust = 0)
  
  return(p)
}

p <- add_manual_legend(p, x1 = 0.05, x2 = 0.09, y1 = 0.75, y2 = 0.70)

ggsave('../plots/ggplot2/ω_versus_ET_over_PET_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

###############
# Omega and VPD

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD, y = n + 0.72, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("VPD (kPa)"),  
  y_lab = bquote("ω"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p"
)

ggsave('../plots/ggplot2/ω_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

p_omega_box <- p_omega_box +
  theme(
    panel.background = element_rect(fill = NA, colour = NA),  # panel area transparent
    plot.background  = element_rect(fill = NA, colour = NA)   # outer background transparent
  )


# Shrink text size in geom_text_repel and annotate("text")
p_omega_box <- resize_plot_elements(
  p_omega_box,
  repel_size = 1,
  text_size  = 3
)

# Combining plots
combined_plot <- p + 
  inset_element(
    p_omega_box,
    left   = 0.3,   # relative x-position (0–1)
    bottom = 0.35,  # relative y-position (0–1)
    right  = 0.99,  # 1.0,  # relative width (0–1)
    top    = 0.99   # 1.0,  # relative height (0–1)
  )

ggsave('../plots/ggplot2/ω_versus_VPD_and_omega_boxplot_ggplot2_TIDY.png', plot = combined_plot, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# For further use
omega_versus_VPD_and_omega_boxplot <- combined_plot

################
# Omega and geff

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = g_eff, y = (n + 0.72), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("g"["eff"]~"(mm s"^"-1"*")"),  
  y_lab = bquote("ω"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p"
)

ggsave('../plots/ggplot2/ω_versus_g_eff_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

##################
# Omega and ET/PET

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = ET / ETo_FAO56_alfalfa, y = (n + 0.72), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("ET / PET"),  
  y_lab = bquote("ω"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p"
)

ggsave('../plots/ggplot2/ω_versus_ET_PET_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#################
# Omega and PET/P

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = ETo_FAO56_alfalfa / P, y = (n + 0.72), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("AI = PET / P"),  
  y_lab = bquote("ω"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p"
)

ggsave('../plots/ggplot2/ω_versus_AI_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------
# Remove ERA5 land
annual_stats_wide <- annual_stats_wide |>
  filter(ENSEMBLE != "ERA5")

# ------------------
# Plot Budyko curves
out_BC_FAO56_alfalfa <- Budyko_curve(annual_stats_wide,
                                     pet_col = "ETo_FAO56_alfalfa", et_col = "ET", p_col = "P",
                                     X_range = c(0.0, 2.2), Y_range = c(0.0, 1.06),
                                     Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
                                     plot = TRUE, plot_name = "../plots/ggplot2/Budyko_curve_ggplot2.png")

out_BC_FAO56_alfalfa_CO2_corr <- Budyko_curve(annual_stats_wide,
                                              pet_col = "ETo_FAO56_alfalfa_GCM_CO2_corr", et_col = "ET", p_col = "P",
                                              X_range = c(0.0, 2.2), Y_range = c(0.0, 1.06),
                                              Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
                                              plot = TRUE, plot_name = "../plots/ggplot2/Budyko_curve_CO2_ggplot2.png")

# Add label "a" to the first plot
p9 <- out_BC_FAO56_alfalfa$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "a", hjust = 1.5, vjust = 1.5, size = 6)

# Add label "b" to the second plot
p10 <- out_BC_FAO56_alfalfa_CO2_corr$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "b", hjust = 1.5, vjust = 1.5, size = 6)

# Combine the two plots side by side using cowplot
panel_figure <- (p9 | p10) 

# Save the combined plot
ggsave('../plots/ggplot2/BC_PET_PM_and_PM_corr_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height, dpi = RES, units = 'mm')

################################################################################

# Add label "a" to the first plot
p9 <- out_BC_FAO56_alfalfa$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "a", hjust = 1.5, vjust = 1.5, size = 6)

# Add label "b" to the second plot
p10 <- out_BC_FAO56_alfalfa_CO2_corr$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "b", hjust = 1.5, vjust = 1.5, size = 6)

# Combine the two plots side by side using cowplot
panel_figure <- (p9 | p10) 

# Save the combined plot
ggsave('../plots/ggplot2/BC_PET_PM_and_PM_corr.png', plot = panel_figure, width = Pl_width*2, height = Pl_height, dpi = RES, units = 'mm')

################################################################################

# Runoff from mHM
Q_EOBS_mHM <- 224.6725
P_EOBS <- 689.7999
ET_EOBS_mHM <- 466.8086

#_______________________________________________________________________________
################################################################################
# Plot of the domain and the main hydrological and atmospheric state variables #
################################################################################
#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

## --- helpers ---
summarize_var <- function(df, var) {
  df |> 
    mutate(
      Period = case_when(
        str_detect(PERIOD, "1981_?2005") ~ "Baseline",
        str_detect(PERIOD, "2076_?2100") ~ "Future",
        TRUE ~ PERIOD
      ),
      Period = factor(Period, levels = c("Baseline", "Future"))
    ) |> 
    group_by(ModelGroup = ENSEMBLE, Period)  |> 
    summarise(
      Mean = mean(.data[[var]], na.rm = TRUE),
      SD   = sd(.data[[var]],   na.rm = TRUE),
      .groups = "drop"
    )
}

create_barplot <- function(data, ylabel, y_min = NULL, y_max = NULL) {
  pd <- position_dodge(width = 0.7)
  p <- ggplot(data, aes(x = Period, y = Mean, fill = ModelGroup)) +
    geom_col(position = pd, width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD, color = ModelGroup),
      position = pd, width = 0.2, size = 0.15
    ) +
    scale_fill_manual(values = c(CMIP5 = COL_CMIP5, CMIP6 = COL_CMIP6, "EUR-44" = COL_RCMs)) +
    scale_color_manual(values = c(CMIP5 = "#333333", CMIP6 = "#333333", "EUR-44" = "#333333")) +
    labs(x = NULL, y = ylabel) +
    scale_x_discrete(
      labels = c(Baseline = "1981–2005", Future = "2076–2100"),
      expand = expansion(mult = 0, add = 0.04)
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 4)),
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 10, margin = margin(r = 2)),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  if (!is.null(y_min) && !is.null(y_max)) p <- p + coord_cartesian(ylim = c(y_min, y_max))
  p
}

add_tb_padding <- function(g, top = 0.05, bottom = 0.005) {
  grid.arrange(nullGrob(), g, nullGrob(), ncol = 1, heights = c(top, 1, bottom))
}

# label_panel <- function(g, lab = "a", col = "#333333", pad_pt = 1) {
#   if (inherits(g, "gg")) g <- ggplotGrob(g)
#   stopifnot(grid::is.grob(g))
#   title_g <- grid::textGrob(
#     lab, x = grid::unit(0.01, "npc"), y = grid::unit(1, "npc"),
#     just = c("left", "top"),
#     gp = grid::gpar(col = col, fontsize = 16, fontface = "bold")
#   )
#   h <- grid::grobHeight(title_g) + grid::unit(pad_pt, "pt")  # enough room for the text
#   gridExtra::arrangeGrob(title_g, g, ncol = 1,
#                          heights = grid::unit.c(h, grid::unit(1, "null")))
# }

label_panel <- function(g, lab = "a", col = "#333333", pad_pt = 1,
                        fontsize = 16, fontface = "bold") {
  if (inherits(g, "gg")) g <- ggplotGrob(g)
  stopifnot(grid::is.grob(g))
  
  title_g <- grid::textGrob(
    lab,
    x = grid::unit(0.01, "npc"), y = grid::unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(col = col, fontsize = fontsize, fontface = fontface)
  )
  
  h <- grid::grobHeight(title_g) + grid::unit(pad_pt, "pt")
  gridExtra::arrangeGrob(title_g, g, ncol = 1,
                         heights = grid::unit.c(h, grid::unit(1, "null")))
}

## --- summaries ---
Data_P   <- summarize_var(annual_stats_wide, "P")
Data_ET  <- summarize_var(annual_stats_wide, "ET")
Data_RO  <- summarize_var(annual_stats_wide, "RO")
Data_Ta  <- summarize_var(annual_stats_wide, "Ta")
Data_RH  <- summarize_var(annual_stats_wide, "RH")
Data_VPD <- summarize_var(annual_stats_wide, "VPD")

## --- panel a (map) ---
# map_image <- image_read("../elevation_map/plots/Elevation_and_domain_map_r-6_v-01_DPI-300.png")
map_image <- image_read("../era5_land_global/plots/VPD_Europe_r-6_v-01_DPI-1200.png")
crop_image_by_percentage <- function(image, left_percent, right_percent, top_percent, bottom_percent) {
  info <- image_info(image); w <- info$width; h <- info$height
  left <- round((left_percent/100) * w);  right <- round((right_percent/100) * w)
  top  <- round((top_percent/100) * h);   bottom <- round((bottom_percent/100) * h)
  cw <- w - left - right; ch <- h - top - bottom
  image_crop(image, geometry = paste0(cw, "x", ch, "+", left, "+", top))
}
# panel_a <- crop_image_by_percentage(map_image, 0, 3.1, 5, 2.6) |> # left_percent, right_percent, top_percent, bottom_percent
panel_a <- crop_image_by_percentage(map_image, 0, 3.1, 5, 2.55) |> # left_percent, right_percent, top_percent, bottom_percent
  grDevices::as.raster() |> rasterGrob(interpolate = TRUE) |>
  label_panel("a", col = "#333333", fontsize = 12) |>
  add_tb_padding(top = 0.02, bottom = 0.005)

## --- plots (barplots) ---
plot_P   <- create_barplot(Data_P,  expression(P ~ (mm ~ yr^{-1})),    y_min = 700, y_max = 1100)
plot_ET  <- create_barplot(Data_ET, expression(ET ~ (mm ~ yr^{-1})),   y_min = 400, y_max = 800)
plot_RO  <- create_barplot(Data_RO, expression(R[o] ~ (mm ~ yr^{-1})), y_min = 100, y_max = 500)

plot_Ta  <- create_barplot(Data_Ta, expression("T"["a"] ~ "(°C)"),     y_min = 5,   y_max = 15)
plot_RH  <- create_barplot(Data_RH, expression("RH" ~ "(%)"),          y_min = 65,  y_max = 85)
plot_VPD <- create_barplot(Data_VPD, expression("VPD" ~ "(kPa)"),      y_min = 0.2, y_max = 0.7)

# Adding the headers to barplots
overlay_header <- function(g, left_lab = "", title = "",
                           top_pad = 0.06,
                           lab_col = "#333333", title_col = "#333333",
                           lab_size = 16, title_size = 12,
                           title_x = 0.5, title_hjust = 0.5,  # <-- center controls
                           lab_x = 0.01, lab_hjust = 0) {
  if (inherits(g, "gg")) g <- ggplotGrob(g)
  y_lab <- 1 - top_pad/2
  cowplot::ggdraw() +
    cowplot::draw_grob(g, x = 0, y = 0, width = 1, height = 1 - top_pad) +
    cowplot::draw_label(left_lab, x = lab_x,   y = y_lab, hjust = lab_hjust, vjust = 0.5,
                        fontface = "bold", size = lab_size, color = lab_col) +
    cowplot::draw_label(title,    x = title_x, y = y_lab, hjust = title_hjust, vjust = 0.5,
                        fontface = "bold", size = title_size, color = title_col)
}

# build the stacked columns as before
panel_hydro_core  <- grid.arrange(plot_P, plot_ET, plot_RO,  ncol = 1, heights = c(1,1,1))
panel_thermo_core <- grid.arrange(plot_Ta, plot_RH, plot_VPD, ncol = 1, heights = c(1,1,1))

# Leave ~7% of the height for the header; tweak 0.04–0.08 to taste
panel_hydro  <- overlay_header(panel_hydro_core,  left_lab = "b", title = "Water balance", 
                               top_pad = 0.07, lab_size = 12, title_x = 0.6, title_hjust = 0.5)
panel_thermo <- overlay_header(panel_thermo_core, left_lab = "",   title = "Atmospheric states",
                               top_pad = 0.07, title_x = 0.6, title_hjust = 0.5)

# convert to grobs and add tiny outer padding if you want
panel_hydro  <- ggplotGrob(panel_hydro)  |> add_tb_padding(top = 0.002, bottom = 0.005)
panel_thermo <- ggplotGrob(panel_thermo) |> add_tb_padding(top = 0.002, bottom = 0.005)

lg_rect <- function(col) rectGrob(width = unit(6, "mm"), height = unit(3, "mm"),
                                  gp = gpar(fill = col, col = NA))
lg_text <- function(txt) textGrob(txt, gp = gpar(fontsize = 11))

legend_row <- gridExtra::arrangeGrob(
  lg_rect(COL_CMIP5), lg_text("CMIP5"),
  lg_rect(COL_CMIP6), lg_text("CMIP6"),
  lg_rect(COL_RCMs), lg_text("EUR-44"),
  ncol = 6, widths = c(1, 2, 1, 2, 1, 2)
)
legend_centered <- gridExtra::arrangeGrob(nullGrob(), legend_row, nullGrob(),
                                          ncol = 3, widths = c(0.3, 1, 0.3))

b_with_legend <- grid.arrange(
  panel_hydro, panel_thermo, legend_centered,
  layout_matrix = rbind(c(1, 2),
                        c(3, 3)),
  widths  = c(1, 1),
  heights = c(1, 0.10)
)

## --- top row: a | spacer | (b+legend) ---
top_block <- grid.arrange(panel_a, nullGrob(), b_with_legend, ncol = 3, widths = c(1.92, 0.13, 1.95))

# <<----------------- Start of the Budyko curve ------------------------------>>

out_BC_FAO56_alfalfa <- Budyko_curve(annual_stats_wide,
                                     pet_col = "ETo_FAO56_alfalfa", et_col = "ET", p_col = "P",
                                     X_range = c(0.0, 2.2), Y_range = c(0.0, 1.2),
                                     Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
                                     boundary_line_col = "gray45",
                                     boundary_line_type = "solid",
                                     boundary_line_size = 0.2,
                                     plot = TRUE, plot_name = "../plots/ggplot2/Budyko_curve_ggplot2.png")
Budyko_plot <- out_BC_FAO56_alfalfa$out_plot
plot_data <- ggplot_build(Budyko_plot)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

Budyko_plot <- Budyko_plot +
  # Text label
  annotate("text",
           x = rescale(0.3, X_range),
           y = rescale(0.9, Y_range),
           label = as.expression(bquote(bold("Energy") ~ italic("vs.") ~ bold("Water")~"limitation")),
           hjust = 0, color = "#333333", size = 4) +
  
  # Left-pointing arrow above "Energy"
  annotate("segment",
           x = rescale(0.42, X_range),  # adjust to match text
           xend = rescale(0.30, X_range),
           y = rescale(0.95, Y_range),
           yend = rescale(0.95, Y_range),
           arrow = arrow(type = "closed", length = unit(0.2, "cm")),
           color = "#333333") +
  
  # Right-pointing arrow above "water"
  annotate("segment",
           x = rescale(0.51, X_range),
           xend = rescale(0.60, X_range),
           y = rescale(0.95, Y_range),
           yend = rescale(0.95, Y_range),
           arrow = arrow(type = "closed", length = unit(0.2, "cm")),
           color = "#2b2b2b")

# Save the plot
ggsave('../plots/ggplot2/Budyko_curve_ggplot2.png', plot = Budyko_plot, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# <<----------------- End of the Budyko curve   ------------------------------>>


## --- bottom row: c | spacer | d ---
panel_c <- Budyko_plot |>
  label_panel("c", col = "#333333", fontsize = 12) |>
  add_tb_padding(top = 0.01, bottom = 0.005)

# <<----------------- Start of dTa vs. dVPD/VPD ------------------------------>>

Plot_labels <- tibble(
  x = c(0.04, 0.04, 0.04),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_VPD_over_VPD, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_VPD_over_VPD"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'VPD/'*'VPD'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_Ta_VPDn")

plot_data <- ggplot_build(p_Ta_VPDn)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

rand_int_I <- 0
CC_linear <- b/(c+T_avg)*X
CC_nonlinear <- (-b*T_avg/(c+T_avg)^2+b/(c+T_avg))*X

# Create data frames for the lines
DF1_CC_linear <- data.frame(X = X, Y = rand_int_I + CC_linear)
DF1_CC_nonlinear <- data.frame(X = X, Y = rand_int_I + CC_nonlinear)

# Normalize the X and Y values
DF1_CC_linear$X_norm <- (DF1_CC_linear$X - X_range[1]) / (X_range[2] - X_range[1])
DF1_CC_linear$Y_norm <- (DF1_CC_linear$Y - Y_range[1]) / (Y_range[2] - Y_range[1])

# Calculate the differences (deltas) for the normalized values
dX_norm <- diff(DF1_CC_linear$X_norm)
dY_norm <- diff(DF1_CC_linear$Y_norm)

# Calculate the angles in radians
angles_radians <- atan2(dY_norm, dX_norm)

# Convert the angles to degrees
angle_degrees <- mean(angles_radians * (180 / pi))


p_Ta_VPDn <- p_Ta_VPDn + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.6, X_range), y = rescale(0.28, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

# <<----------------- End of dTa vs. dVPD/VPD -------------------------------->>

panel_d <- p_Ta_VPDn |>
  label_panel("d", col = "#333333", fontsize = 12) |>
  add_tb_padding(top = 0.01, bottom = 0.005)

bottom_block <- grid.arrange(panel_c, nullGrob(), panel_d, ncol = 3, widths = c(2, 0.1, 2))

## --- final stack ---
final_figure <- grid.arrange(top_block, bottom_block, ncol = 1, heights = c(1, 0.8))

# grid.newpage(); grid.draw(final_figure)
# ggsave("../plots/ggplot2/Elevation_water_balance&atmosphere.png",
ggsave("../plots/ggplot2/Domain_VPD_water_balance&atmosphere.png",
       #plot = final_figure, width = Pl_width*3*0.65, height = Pl_height*3.1*0.65,
       plot = final_figure, width = 240, height = 248,
       # dpi = RES, units = "mm")
       dpi = 1200, units = "mm", bg = "white")

graphics.off()
#__________________________________________________________
###########################################################
# Budyko curve components against VPD perturbation analysis
###########################################################
#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

Fit_ET_norm  <- lm(Data_to_plot$diffs$d_ET_over_ET~Data_to_plot$diffs$d_VPD_over_VPD)
Fit_PET_norm <- lm(Data_to_plot$diffs$d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa~Data_to_plot$diffs$d_VPD_over_VPD)
Fit_P_norm   <- lm(Data_to_plot$diffs$d_P_over_P~Data_to_plot$diffs$d_VPD_over_VPD)
Fit_AI_norm  <- lm(Data_to_plot$diffs$d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa~Data_to_plot$diffs$d_VPD_over_VPD)
Fit_EI_norm  <- lm(Data_to_plot$diffs$d_EI_over_EI~Data_to_plot$diffs$d_VPD_over_VPD)

plot(Data_to_plot$abs$VPD, Data_to_plot$abs$ET)
plot(Data_to_plot$abs$VPD, Data_to_plot$abs$ETo_FAO56_alfalfa)
plot(Data_to_plot$abs$VPD, Data_to_plot$abs$P)
plot(Data_to_plot$abs$VPD, Data_to_plot$abs$AI_FAO56_alfalfa)
plot(Data_to_plot$abs$VPD, Data_to_plot$abs$EI)

# See the distribution of VPD
ggplot(Data_to_plot$abs |> 
         filter(!is.na(VPD)),  # remove missing VPD values
       aes(x = VPD, fill = ensemble, color = ensemble)) +
  geom_density(alpha = 0.3, adjust = 2) +
  facet_wrap(~ PERIOD, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density plot of VPD across ensembles and periods",
    x = "Vapor Pressure Deficit (VPD)",
    y = "Density",
    fill = "Model",
    color = "Model"
  )

#######################
# Normalized AI vs. VPD

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'AI / AI = '*Delta*'(PET/P) / (PET/P)'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_AI_VPD_norm")

# Save the plot
ggsave('../plots/ggplot2/delta_AI_over_AI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_AI_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#######################
# Normalized EI vs. VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_EI_over_EI, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_EI_over_EI"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'EI / EI = '*Delta*'(ET/P) / (ET/P)'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_EI_over_EI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_EI_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Normalized PET vs. VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'PET/'*'PET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_PET_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_PET_over_PET_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_PET_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Normalized P vs. VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_P_over_P_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_P_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Normalized ET vs. VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_ET_over_ET_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_ET_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Normalized geff vs. VPD
Plot_labels <- tibble(
  x = c(0.82, 0.82, 0.82),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

LM_eq_labels <- tibble(
  x = c(0.04, 0.04, 0.04),
  y = c(0.14, 0.09, 0.04)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_g_eff_over_g_eff, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_g_eff_over_g_eff"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*g[eff]~"/"~g[eff]),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels, plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_geff_VPD_norm")

# Save the plot
ggsave('../plots/ggplot2/delta_geff_over_geff_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_geff_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Optional but helps: tiny, identical plot margins
tight <- theme(plot.margin = margin(3,3,3,3))
pA <- p_AI_VPD_norm + tight
pB <- p_EI_VPD_norm + tight
pC <- p_PET_VPD_norm + tight
pD <- p_P_VPD_norm   + tight
pE <- p_ET_VPD_norm  + tight

# (2) Align *all five* first; keep all four sides ('tblr')
aligned <- cowplot::align_plots(pA, pB, pC, pD, pE, align = "hv", axis = "tblr")

# (3) Build rows from the aligned grobs
# --- top row with spacing between A and B ---
top <- cowplot::plot_grid(
  aligned[[1]], NULL, aligned[[2]],
  ncol = 3,
  rel_widths = c(1, 0.05, 1),
  labels = c("a", "", "b"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# --- bottom row with spacing between C, D, E ---
bottom <- cowplot::plot_grid(
  aligned[[3]], NULL, aligned[[4]], NULL, aligned[[5]],
  ncol = 5,
  rel_widths = c(1, 0.05, 1, 0.05, 1),
  labels = c("c", "", "d", "", "e"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# --- add vertical spacing between top and bottom ---
combined <- cowplot::plot_grid(
  NULL, top, NULL, bottom,
  ncol = 1,
  rel_heights = c(0.1, 3, 0.1, 1.9)  # the middle '0.1' is vertical gap
)

# (5) Save (use any size you like; e.g. 240×200 mm fits a 3:2 row split nicely)
ggsave("../plots/ggplot2/delta_AI,EI,PET,P,ET_norm_vs_delta_VPD_norm_ggplot2_TIDY.png", combined,
       width = 240, height = 192, units = "mm", dpi = RES, bg = "white")

#-------------------------------------------------------------------------------
# Normalized EI vs. VPD check

tmp <- Data_to_plot$diffs
tmp <- tmp |> 
  mutate(EI_1st_order_predict = d_ET_over_ET - d_P_over_P)

make_scatter_plot(data = tmp |>
                    select(d_EI_over_EI, EI_1st_order_predict, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_EI_over_EI", y = "EI_1st_order_predict"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'EI / EI = '*Delta*'(ET/P) / (ET/P)'),  y_lab = bquote(Delta*'EI / EI = '*Delta*'(ET/P) / (ET/P)'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_EI_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_EI_over_EI_vs_delta_EI_over_EI_ggplot2_TIDY.png', plot = p_EI_EI_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

rm(tmp)

#####################################
# Normalized AI_CO2 corrected vs. VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_AI_ETo_FAO56_alfalfa_GCM_CO2_corr_over_AI_ETo_FAO56_alfalfa_GCM_CO2_corr, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_AI_ETo_FAO56_alfalfa_GCM_CO2_corr_over_AI_ETo_FAO56_alfalfa_GCM_CO2_corr"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*AI[CO[2]]/AI),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_CO2corr_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_AI_CO2_corr_over_AI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_AI_CO2corr_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#______________________________________
#######################################
# Plot of the P, ET and VPD relations #
#######################################
#‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

###########
# P vs. VPD
make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = VPD, y = P, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("P (mm/yr)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD")
# Save the plot
ggsave('../plots/ggplot2/P_vs_VPD_ggplot2_TIDY.png', plot = p_P_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


###########
# ET vs. VPD
make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = VPD, y = ET, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE,
                  xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  # plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_ET_VPD")
# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_ggplot2_TIDY.png', plot = p_ET_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

##############
# geff vs. VPD
VPD_g_eff_data <- Data_to_plot$abs |> select(VPD, g_eff, PERIOD)
VPD_g_eff_data <- VPD_g_eff_data |>
  mutate(log_VPD = log(VPD),
         sqrt_VPD = 1/sqrt(VPD))

plot(VPD_g_eff_data |> select(VPD, g_eff))
plot(VPD_g_eff_data |> select(log_VPD, g_eff))
plot(VPD_g_eff_data |> select(sqrt_VPD, g_eff))

VPD_g_eff_data |> select(VPD, g_eff) |> plot()

fit <- lm(g_eff ~ log_VPD, data = VPD_g_eff_data)
summary(fit)

# Fit models
fit_log  <- lm(g_eff ~ log_VPD,  data = VPD_g_eff_data)
fit_sqrt <- lm(g_eff ~ sqrt_VPD, data = VPD_g_eff_data)

# Extract coefficients and R^2
coef_log  <- coef(fit_log)
r2_log    <- summary(fit_log)$r.squared

coef_sqrt <- coef(fit_sqrt)
r2_sqrt   <- summary(fit_sqrt)$r.squared

# Build equation labels
eq_log  <- sprintf("y = %.3f + %.3f*x,  R² = %.2f", coef_log[1], coef_log[2], r2_log)
eq_sqrt <- sprintf("y = %.3f + %.3f*x,  R² = %.2f", coef_sqrt[1], coef_sqrt[2], r2_sqrt)

# Plot 1: g_eff ~ log(VPD)
p1 <- ggplot(VPD_g_eff_data, aes(x = log_VPD, y = g_eff)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text", x = max(VPD_g_eff_data$log_VPD, na.rm = TRUE),
           y = max(VPD_g_eff_data$g_eff, na.rm = TRUE),
           label = eq_log, hjust = 1, vjust = 1) +
  labs(x = "log(VPD)", y = "g_eff")

# Plot 2: g_eff ~ sqrt(VPD)
p2 <- ggplot(VPD_g_eff_data, aes(x = sqrt_VPD, y = g_eff)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text", x = max(VPD_g_eff_data$sqrt_VPD, na.rm = TRUE),
           y = max(VPD_g_eff_data$g_eff, na.rm = TRUE),
           label = eq_sqrt, hjust = 1, vjust = 1) +
  labs(x = bquote(1/sqrt(VPD)),   # √VPD
       y = "g_eff")

# Plot 3: g_eff ~ sqrt(VPD)
fit_sqrt_p3 <- lm(g_eff ~ sqrt_VPD, data = VPD_g_eff_data |> filter(PERIOD=="1981_2005"))
coef_sqrt_p3 <- coef(fit_sqrt_p3)
r2_sqrt_p3   <- summary(fit_sqrt_p3)$r.squared
eq_sqrt_p3 <- sprintf("y = %.3f + %.3f*x,  R² = %.2f", coef_sqrt_p3[1], coef_sqrt_p3[2], r2_sqrt_p3)

p3 <- ggplot(VPD_g_eff_data |> filter(PERIOD=="1981_2005"), aes(x = sqrt_VPD, y = g_eff)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate("text", x = max(VPD_g_eff_data$sqrt_VPD, na.rm = TRUE),
           y = max(VPD_g_eff_data$g_eff, na.rm = TRUE),
           label = eq_sqrt_p3, hjust = 1, vjust = 1) +
  labs(x = bquote(1/sqrt(VPD)),   # √VPD
       y = "g_eff")

# Save
ggsave("../plots/ggplot2/g_eff_vs_logVPD.png", p1, width = 4, height = 4, dpi = 300)
ggsave("../plots/ggplot2/g_eff_vs_sqrtVPD.png", p2, width = 4, height = 4, dpi = 300)

# 1) Fit models (if not already fit)
fit_log  <- lm(g_eff ~ log_VPD,  data = VPD_g_eff_data)
fit_sqrt <- lm(g_eff ~ sqrt_VPD, data = VPD_g_eff_data)

peak_log <- exp((-coef(fit_log)[1] - coef(fit_log)[2]) / coef(fit_log)[2])
peak_sqrt <- (coef(fit_sqrt)[2] / (-2 * (coef(fit_sqrt)[1])))^2

# 2) Prediction grid
new_VPD <- seq(0, 1, by = 0.001)

# --- Log fit line (drop VPD == 0 since log(0) = -Inf) ---
new_log <- tibble(VPD = new_VPD) |> 
  mutate(log_VPD = log(VPD)) |> 
  filter(is.finite(log_VPD)) |> 
  mutate(
    predicted_g_eff = predict(fit_log, newdata = pick(everything())),
    predicted_rs_eff = 1 / (predicted_g_eff / 1000),
    predicted_LE = (rhoAir * CpAir / gamma) * VPD / predicted_rs_eff
  )

lambda = (2.501 * (10^6) - 2361 * 8) / (10^6)
# Unit conversions
W_to_mm = 1 / lambda * 3600*24 / 10^6

new_log <- new_log |> 
  mutate(predicted_ET = predicted_LE * W_to_mm * 365.25) |> 
  transmute(VPD, predicted_ET, fit = "log(VPD)")

# --- Sqrt fit line (VPD=0 is fine here) ---
new_sqrt <- tibble(VPD = new_VPD) |> 
  mutate(sqrt_VPD = 1/sqrt(VPD)) |> 
  mutate(
    predicted_g_eff = predict(fit_sqrt, newdata = cur_data()),
    predicted_rs_eff = 1 / (predicted_g_eff / 1000),
    predicted_LE = (rhoAir * CpAir / gamma) * VPD / predicted_rs_eff
  ) |> 
  mutate(predicted_ET = predicted_LE * W_to_mm * 365.25) |> 
  transmute(VPD, predicted_ET, fit = "sqrt(VPD)")

# 3) Combine both lines
lines_df <- bind_rows(new_log, new_sqrt)

# 4) Add both lines to the existing ggplot
p_ET_VPD <- p_ET_VPD +
  
  geom_line(
    data = lines_df,
    aes(x = VPD, y = predicted_ET, linetype = fit),
    color = "black"
  ) +
  geom_vline(xintercept = peak_log, 
             colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = peak_sqrt, 
             colour = "grey35", linetype = "solid",  linewidth = 0.25) +
  scale_linetype_manual(
    values = c("log(VPD)" = "dashed", "sqrt(VPD)" = "solid"),
    labels = c("log(VPD)" = bquote(log(VPD)),
               "sqrt(VPD)" = bquote(1/sqrt(VPD)))
  ) +
  labs(linetype = NULL) +   # remove legend title
  theme(
    legend.position = c(0.86, 0.99),        # inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, color = NA)  # no square
  )

p_ET_VPD <- p_ET_VPD +
  annotate("point", x = 0.65, y = 715, shape = 21, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.69, y = 715, label = "1981–2005", size = 3, hjust = 0)

p_ET_VPD <- p_ET_VPD +
  annotate("point", x = 0.65, y = 695, shape = 24, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.69, y = 695, label = "2076–2100 RCP 8.5", size = 3, hjust = 0)

# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_ET_VPD$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_ET_VPD$layers <- c(p_ET_VPD$layers[-idx_pts_txt], p_ET_VPD$layers[idx_pts_txt])
}

# Save again
ggsave(
  '../plots/ggplot2/ET_vs_VPD_ggplot2_TIDY.png',
  plot = p_ET_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm'
)

#######################
# ET vs. VPD arrow plot
plot_data <- ggplot_build(p_ET_VPD)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

# Build arrow start/end positions per (model, label) between the two periods
arrow_data <- Data_to_plot$abs  |> 
  filter(PERIOD %in% c("1981_2005", "2076_2100"),
         !is.na(VPD), !is.na(ET)) |>
  mutate(PERIOD = factor(PERIOD, levels = c("1981_2005", "2076_2100"))) |>
  group_by(ensemble, model) |>
  arrange(PERIOD, .by_group = TRUE) |>
  summarise(
    x     = first(VPD),         # start x (1981–2005)
    y     = first(ET),          # start y (1981–2005)
    xend  = last(VPD),          # end x (2076–2100)
    yend  = last(ET),           # end y (2076–2100)
    color = first(color),       # line colour comes from tibble
    .groups = "drop"
  ) |>
  filter(!is.na(x), !is.na(y), !is.na(xend), !is.na(yend)) |>
  mutate(
    Scenario = ensemble,  # used only for thickness rule below
    # “thick” if EUR-44 or model mentions as driving GCM
    driving_gcm = str_detect(toupper(model), paste0(driving_GCMs, collapse = "|")),
    arrow_type = if_else(Scenario == "RCMs" | driving_gcm, "thick", "thin"),
    line_size  = if_else(arrow_type == "thick", 0.8, 0.3)  # mapped via scale_size_manual below
  )

# Base plot:
# - Black fit lines are added FIRST so they appear behind the coloured arrows
# - Vertical reference lines next
# - Arrows last (on top)
p_ET_VPD_arrow <- ggplot(arrow_data) +
  geom_line(
    data = lines_df,
    aes(x = VPD, y = predicted_ET, linetype = fit),
    color = "black"
  ) +
  geom_vline(xintercept = peak_log, 
             colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = peak_sqrt, 
             colour = "grey35", linetype = "solid",  linewidth = 0.25) +
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend,
        colour = color, size = arrow_type),
    lineend = "round",
    arrow = arrow(length = unit(2.5, "mm"), type = "closed")
  ) +
  # Identity scales: use colours/shapes/linetypes as provided in the data
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  # Map “thin/thick” to numeric linewidths (no legend for size)
  scale_size_manual(values = c(thin = 0.3, thick = 0.8), guide = "none") +
  labs(
    x = "VPD (kPa)",
    y = bquote("ET (mm yr"^"-1"*")")
  ) +
  # Keep your theme stack; theme_bw() overrides theme_minimal()
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Styling for the two fitted lines (legend text and line types)
p_ET_VPD_arrow <- p_ET_VPD_arrow +
  scale_linetype_manual(
    values = c("log(VPD)" = "dashed", "sqrt(VPD)" = "solid"),
    labels = c("log(VPD)" = bquote(log(VPD)),
               "sqrt(VPD)" = bquote(1/sqrt(VPD)))
  ) +
  labs(linetype = NULL) +   # remove legend title
  theme(
    legend.position = c(0.96, 0.99),               # inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, color = NA)  # no legend box
  ) +
  # Fix the visible window to match an existing plot’s ranges
  scale_x_continuous(limits = X_range, expand = c(0, 0)) +
  scale_y_continuous(limits = Y_range, expand = c(0, 0))

# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_arrow_ggplot2_TIDY.png',
       plot = p_ET_VPD_arrow, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

##########################
# g_eff vs. VPD single fit

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = 1/sqrt(VPD), y = g_eff, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, 
  xy_round = 0.05, xy_offset = 0.04, X_range_man = c(0.904, 2.496),
  x_lab = bquote(1/sqrt(VPD)~"(kPa"^"-0.5"*")"),  
  y_lab = bquote('g'['eff']~"(mm s"^"-1"*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p_geff_sqrt_VPD_single_fit"
)

# 1) Extract the first point layer from `p`
idx_pt <- which(vapply(p_geff_sqrt_VPD_single_fit$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
stopifnot(!is.na(idx_pt))  # ensure we found a point layer

pts <- layer_data(p_geff_sqrt_VPD_single_fit, idx_pt)
pts <- pts[is.finite(pts$x) & is.finite(pts$y), ]  # keep valid points

# 2) Fit OLS on the plotted coordinates (same space as the plot)
fit <- lm(y ~ x, data = pts)
b0  <- unname(coef(fit)[1])
b1  <- unname(coef(fit)[2])
r2  <- summary(fit)$r.squared

# Build a pretty, signed label with 2 decimals
sign <- if (b1 >= 0) "+" else "−"
# eq_label <- sprintf("y = %.2f %s %.2fx; R\u00B2 = %.2f",
#                     b0, sign, abs(b1), r2)
eq_label <- sprintf(
  "y = %.2f %s %.2f x\nR\u00B2 = %.2f",
  b0, sign, abs(b1), r2
)

# 3) Add line + 95% CI and the equation label (top-left)
plot_data <- ggplot_build(p_geff_sqrt_VPD_single_fit)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

p_geff_sqrt_VPD_single_fit <- p_geff_sqrt_VPD_single_fit +
  geom_smooth(
    data = pts,
    aes(x = x, y = y),
    method = "lm", formula = y ~ x,
    se = TRUE,
    inherit.aes = FALSE,
    color = "black", fill = "grey40", alpha = 0.25, linewidth = 0.6
  ) +
  annotate("text",
           # x = rescale(0.54, X_range), y = rescale(0.05, Y_range), label = eq_label,
           x = rescale(0.15, X_range), y = rescale(0.60, Y_range), label = eq_label,
           hjust = 0, vjust = 0.5, size = 3.6)

add_manual_legend <- function(p_geff_sqrt_VPD_single_fit, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70){
  plot_data <- ggplot_build(p_geff_sqrt_VPD_single_fit)
  X_range <- plot_data$layout$panel_params[[1]]$x.range
  Y_range <- plot_data$layout$panel_params[[1]]$y.range
  
  # Manually add two points and their labels to the top right
  p_geff_sqrt_VPD_single_fit <- p_geff_sqrt_VPD_single_fit +
    annotate("point", x = rescale(x1, X_range), y = rescale(y1, Y_range), shape = 21, size = 2.5, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y1, Y_range), label = "1981–2005", size = 3.5, hjust = 0)
  
  p_geff_sqrt_VPD_single_fit <- p_geff_sqrt_VPD_single_fit +
    annotate("point", x = rescale(x1, X_range), y = rescale(y2, Y_range), shape = 24, size = 2.5, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y2, Y_range), label = "2076–2100 RCP 8.5", size = 3.5, hjust = 0)
  
  return(p_geff_sqrt_VPD_single_fit)
}

p_geff_sqrt_VPD_single_fit <- add_manual_legend(p_geff_sqrt_VPD_single_fit, x1 = 0.06, x2 = 0.11, y1 = 0.78, y2 = 0.73)

# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_geff_sqrt_VPD_single_fit$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_geff_sqrt_VPD_single_fit$layers <- c(p_geff_sqrt_VPD_single_fit$layers[-idx_pts_txt], p_geff_sqrt_VPD_single_fit$layers[idx_pts_txt])
}

ggsave('../plots/ggplot2/g_eff_versus_sqrt_VPD_single_fit_ggplot2_TIDY.png', plot = p_geff_sqrt_VPD_single_fit, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


###################################
# g_eff vs. VPD ensemble-specif fit
LM_eq_labels <- tibble(
  x = rep(0.03, 6),
  y = rep(1, 6) - seq(from = 0.12, to = 0.5, length.out = 6)
)

make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = 1/sqrt(VPD), y = g_eff, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, LM_eq_labels = LM_eq_labels,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote(1/sqrt(VPD)~"(kPa"^"-0.5"*")"),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p_geff_sqrt_VPD_ensemble_fit"
)

ggsave('../plots/ggplot2/g_eff_versus_sqrt_VPD_ensemble_fit_ggplot2_TIDY.png', plot = p_geff_sqrt_VPD_ensemble_fit, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

FITS <- Data_to_plot$abs |>
  mutate(
    x = 1 / sqrt(VPD),
    y = g_eff,
    ensemble = interaction(ensemble, PERIOD, drop = TRUE)
  ) |>
  group_by(model) |>
  summarise(
    intercept = coef(lm(y ~ x))[1],
    slope     = coef(lm(y ~ x))[2],
    r2        = summary(lm(y ~ x))$r.squared,
    .groups = "drop"
  )


# 2) Rebuild the data with x,y,model and join the coefs
Data_with_pred <- Data_to_plot$abs |> 
  mutate(
    x = 1/sqrt(VPD),
    y = g_eff,
    ensemble = interaction(ensemble, PERIOD, drop = TRUE)
  ) |> 
  left_join(FITS, by = "model") |> 
  mutate(g_eff_predicted = intercept + slope * x)

# Result:
Data_with_pred |> select(g_eff, g_eff_predicted) |> plot()

Data_with_pred <- Data_with_pred |> 
  mutate(ET_predicted = (rhoAir * CpAir / gamma) * VPD /
           (1 / (g_eff_predicted / 1000)) *                                       # Surface resistance
           1 / ((2.501 * (10^6) - 2361 * Ta) / (10^6)) * 3600*24 / 10^6 * 365.25  # From LE to annual ET
  )

# Predicted vs. original ET
make_scatter_plot(data = Data_with_pred  |>
                    select(ET, ET_predicted, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "ET", y = "ET_predicted"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(ET~"(mm/yr)"),  y_lab = bquote(ET[predicted]~"(mm/yr)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_ET")

# Save the plot
ggsave('../plots/ggplot2/ET_predicted_vs_ET.png', plot = ET_predicted_vs_ET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

################################################################################
# --- VPD grid (avoid 0 because 1/sqrt(0) = Inf) ---
# If you want to match your data range: use range(Data_to_plot$abs$VPD, na.rm = TRUE)
VPD_grid <- tibble(VPD = seq(1e-4, 1, by = 0.001))

# --- Physical constants (use your existing values/objects) ---
# rhoAir, CpAir, gamma already defined in your environment
# keep your same ET conversion (daily-to-annual in your earlier code)
ET_conv <- W_to_mm   * 365.25

# --- Build lines for all 6 fits (one per model) ---
lines_df <- tidyr::crossing(FITS, VPD_grid) |>
  mutate(
    x = 1/sqrt(VPD),
    predicted_g_eff = intercept + slope * x,
    predicted_rs_eff = 1 / (predicted_g_eff / 1000),
    predicted_LE = (rhoAir * CpAir / gamma) * VPD / predicted_rs_eff,
    predicted_ET = predicted_LE * ET_conv
  ) |>
  # keep tidy columns; rename 'model' to 'fit' if you want to match your prior naming
  transmute(VPD, predicted_ET, fit = model) |>
  # optional: drop non-finite or negative-g_eff induced artifacts
  filter(is.finite(predicted_ET))

# Now lines_df contains one predicted ET curve per 'fit' (i.e., per model)
# Example to plot:
ggplot(lines_df, aes(VPD, predicted_ET, color = fit)) +
  geom_line() +
  coord_cartesian(ylim = c(400, 700)) +
  theme_bw() +
  theme(panel.grid = element_blank())

#---------------------------------------------------
# Check the residuals of the ensemble specific model
Data_with_pred <- Data_with_pred |> 
  mutate(g_eff_resids = g_eff - g_eff_predicted)

Data_with_pred |> pull(g_eff_resids) |> plot()

make_scatter_plot(data = Data_with_pred  |>
                    select(VPD, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_VPD.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Global radiation
make_scatter_plot(data = Data_with_pred |>
                    select(Rg, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Rg", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(R[g]),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_Rg.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Precipitation
make_scatter_plot(data = Data_with_pred |>
                    select(P, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "P", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(P),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_P.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Temperature
make_scatter_plot(data = Data_with_pred |>
                    select(Ta, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Ta", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(P),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_Ta.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

################################################################################
# Jarvis type model

# Predicted geff vs. observed
make_scatter_plot(data = Data_to_plot$abs |>
                    select(ETo_FAO56_alfalfa, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "ETo_FAO56_alfalfa", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(352, 1648), Y_range_man = c(384, 816),
                  x_lab = bquote("PET (mm yr"^"-1"*")"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_vs_PET")

# Save the plot
ggsave('../plots/ggplot2/ET_vs_PET.png', plot = ET_vs_PET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

names(Data_to_plot$abs)

Data_to_plot$abs |>
  select(PET, ET) |>
  ggplot(aes(x = PET, y = ET)) +
  geom_point()


plot(1/Data_to_plot$abs$rs_eff * 1000, Data_to_plot$abs$g_eff)

g_eff_corr_test <- 1 / (rs_CO2(Data_to_plot$abs$rs_eff)) * 1000
g_eff_test <- Data_to_plot$abs$g_eff

Fit_test <- lm(g_eff_corr_test~g_eff_test -1)
plot(g_eff_test, g_eff_corr_test)
abline(Fit_test)

df_test <- data.frame(
  g_eff_test = g_eff_test,
  g_eff_corr_test = g_eff_corr_test
)

test_plot <- ggplot(df_test, aes(x = g_eff_test, y = g_eff_corr_test)) +
  geom_point(color = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey30", size = 1) +
  annotate("text", x = 14, y = 13.2, label = "1:1", color = "grey30") +
  labs(
    x = expression(g[eff] ~ "(mm s"^{-1}*")"),
    y = expression(g[eff[corr]] ~ "(mm s"^{-1}*")")
  ) +
  coord_cartesian(xlim = c(-0.6, 15.6), ylim = c(-0.6, 15.6)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8) # bounding box
  )

ggsave('../plots/ggplot2/geff_corr_vs_geff.png', plot = test_plot, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm', bg = "white")

#------------------------------
# Fitting Jarvis reduction model
Fit_Jarvis <- FALSE
if(Fit_Jarvis){
  source("./src/Jarvis_model_min_max.R")
  list2env(jarvis_bundle, envir = .GlobalEnv)
}else{
  load("./RData/20251226_jarvis_objects.RData")
  list2env(jarvis_bundle, envir = .GlobalEnv)
}

Fit_RF <- FALSE
if(Fit_RF){
  source("./src/RF_hybrid_vpd_constrained_geff_et.R")
  list2env(rf_hybrid_bundle, envir = .GlobalEnv)
}else{
  load("./RData/20251226_RF_objects.RData")
  list2env(rf_hybrid_bundle, envir = .GlobalEnv)
}

source("./src/Jarvis&RF_permutation_importance.R")

run_permutation_importance(
  bundle_path = "./RData/20251226_jarvis_objects.RData",
  perm_B = 1000,
  q_lo = 0.25,
  q_hi = 0.75
)

run_permutation_importance(
  bundle_path = "./RData/20251226_RF_objects.RData",
  perm_B = 1000,
  q_lo = 0.25,
  q_hi = 0.75
)

source("./src/Jarvis&RF_ALE_plots.R")

res_jarvis <- run_ALE_plots(
  bundle_path = "./RData/20251226_jarvis_objects.RData",
  n_bins = 50,
  smooth_span = 0.5
)

res_rf <- run_ALE_plots(
  bundle_path = "./RData/20251226_RF_objects.RData",
  n_bins = 50,
  smooth_span = 0.5
)

source("./src/Jarvis&RF_PDP_plots.R")

res_j <- run_PDP_plots("./RData/20251226_jarvis_objects.RData")
res_r <- run_PDP_plots("./RData/20251226_RF_objects.RData")

pdp_jarvis <- res_j$pdp
pdp_rf     <- res_r$pdp


# Extract data frames from the Jarvis and RF bundle
jarvis_out <- jarvis_bundle$output_df
rf_out <- rf_hybrid_bundle$output_df

# # Ensure, there are no underscores in the model names
# jarvis_out <- jarvis_out |> mutate(model = gsub("_", "-", model))
# rf_out <- rf_out |> mutate(model = gsub("_", "-", model))

# Predicted ET vs. observed RF
make_scatter_plot(data = rf_out  |>
                    mutate(x = ET, y = ET_pred, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(400-16, 800+16), Y_range_man = c(400-16, 800+16),
                  x_lab = bquote(ET),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_ET")

# Save the plot
ggsave('../plots/ggplot2/ET_predicted_vs_ET_RF.png', plot = ET_predicted_vs_ET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
cor.test(rf_out$ET, rf_out$ET_pred)

# Predicted ET vs. observed
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = ET, y = ET_pred, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(400-16, 800+16), Y_range_man = c(400-16, 800+16),
                  x_lab = bquote(ET),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_ET")

# Save the plot
ggsave('../plots/ggplot2/ET_predicted_vs_ET.png', plot = ET_predicted_vs_ET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
cor.test(jarvis_out$ET, jarvis_out$ET_pred)

# Predicted geff vs. observed
make_scatter_plot(data = rf_out |>
                    select(g_eff, g_eff_pred, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "g_eff", y = "g_eff_pred"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(1-0.52, 14+0.52), Y_range_man = c(1-0.52, 14+0.52),
                  x_lab = bquote(g_eff),  y_lab = bquote("g"["eff predicted"]),
                  hline = FALSE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="geff_predicted_vs_geff_RF")

# Save the plot
ggsave('../plots/ggplot2/geff_predicted_vs_geff_RF.png', plot = geff_predicted_vs_geff_RF, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
cor.test(jarvis_out$g_eff, jarvis_out$g_eff_pred)

# Predicted geff vs. observed
make_scatter_plot(data = jarvis_out |>
                    select(g_eff, g_eff_pred, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "g_eff", y = "g_eff_pred"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(1-0.52, 14+0.52), Y_range_man = c(1-0.52, 14+0.52),
                  x_lab = bquote(g_eff),  y_lab = bquote("g"["eff predicted"]),
                  hline = FALSE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="geff_predicted_vs_geff")

# Save the plot
ggsave('../plots/ggplot2/geff_predicted_vs_geff.png', plot = geff_predicted_vs_geff, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
cor.test(jarvis_out$g_eff, jarvis_out$g_eff_pred)


# Predicted ET vs. VPD
make_scatter_plot(data = rf_out |>
                    select(VPD, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = c(384, 816),
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_VPD_AFFM")

VPD_vect <- seq(min(rf_out$VPD, na.rm = TRUE), max(rf_out$VPD, na.rm = TRUE), 0.001)


K_ET_T_avg <- (rhoAir * CpAir / gamma) * VPD_vect * (1/1000) *
  1 / ((2.501e6 - 2361 * T_avg) / 1e6) * (3600 * 24) / 1e6 * 365.25

AFFM <- K_ET_T_avg * (b0 + b1 / sqrt(VPD_vect))

# VPD critical:
# -- numerically
VPD_vect[which.max(AFFM)]
# -- analytically
VPD_at_AFFM_peak <- b1^2 / (4 * b0^2)

ET_predicted_vs_VPD_AFFM <- ET_predicted_vs_VPD_AFFM +
  geom_vline(xintercept = VPD_at_AFFM_peak, 
             colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_line(
    data = tibble(VPD = VPD_vect, ET = AFFM),
    aes(x = VPD, y = ET),
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

ET_predicted_vs_VPD_AFFM <- put_line_behind(ET_predicted_vs_VPD_AFFM)

# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_AFFM.png', plot = ET_predicted_vs_VPD_AFFM, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# ------------------------------------------------------------------------------
# Extract ALE data
# ale_all_sm <- res_jarvis$ale_smoothed
ale_all_sm <- res_rf$ale_smoothed

# Predicted ET vs. VPD
make_scatter_plot(data = jarvis_out |>
                    select(VPD, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = c(384, 816),
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_VPD")

VPD_at_ET_peak <- ale_all_sm |> 
  filter(var == "VPD") |> 
  slice_max(ALE_ET_abs_s, n = 1) |> pull(x)

ET_predicted_vs_VPD <- ET_predicted_vs_VPD +
  geom_vline(xintercept = VPD_at_ET_peak, 
             colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_smooth(
    data = ale_all_sm |>  filter(var == "VPD") |>  select(x, ALE_ET_abs_s),
    aes(x = x, y = ALE_ET_abs_s),
    color = "black",
    linetype = "solid",
    linewidth = 0.5,
    se = FALSE,        # no confidence interval shading
    method = "loess",  # or "gam" if you prefer
    span = 0.6         # adjust smoothness (default ~0.75)
  )

ET_predicted_vs_VPD <- put_line_behind(ET_predicted_vs_VPD)

# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_ALE.png', plot = ET_predicted_vs_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#--------------------
# Predicted ET vs. Rg
make_scatter_plot(data = jarvis_out |>
                    select(Rg, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Rg", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(70, 170), Y_range_man = c(384, 816),
                  x_lab = bquote("R"[g]~"(W m"^"-2"*")"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_Rg")


ET_predicted_vs_Rg <- ET_predicted_vs_Rg +
  geom_smooth(
    data = ale_all_sm |>  filter(var == "Rg") |>  select(x, ALE_ET_abs_s),
    aes(x = x, y = ALE_ET_abs_s),
    color = "black",
    linetype = "solid",
    linewidth = 0.5,
    se = FALSE,        # no confidence interval shading
    method = "loess",  # or "gam" if you prefer
    span = 0.6         # adjust smoothness (default ~0.75)
  )

ET_predicted_vs_Rg <- put_line_behind(ET_predicted_vs_Rg)

# Save the plot
ggsave('../plots/ggplot2/ET_vs_Rg_ALE.png', plot = ET_predicted_vs_Rg, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#--------------------
# Predicted ET vs. Ta
make_scatter_plot(data = jarvis_out |>
                    select(Ta, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Ta", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(5-0.48, 17.48), Y_range_man = c(384, 816),
                  x_lab = bquote("T"[a]~"(°C)"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_Ta")


ET_predicted_vs_Ta <- ET_predicted_vs_Ta +
  geom_smooth(
    data = ale_all_sm |>  filter(var == "Ta") |>  select(x, ALE_ET_abs_s),
    aes(x = x, y = ALE_ET_abs_s),
    color = "black",
    linetype = "solid",
    linewidth = 0.5,
    se = FALSE,        # no confidence interval shading
    method = "loess",  # or "gam" if you prefer
    span = 0.6         # adjust smoothness (default ~0.75)
  )

ET_predicted_vs_Ta <- put_line_behind(ET_predicted_vs_Ta)

# Save the plot
ggsave('../plots/ggplot2/ET_vs_Ta_ALE.png', plot = ET_predicted_vs_Ta, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------
# Predicted ET vs. P
make_scatter_plot(data = jarvis_out |>
                    select(P, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "P", y = "ET"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(576, 1224), Y_range_man = c(384, 816),
                  x_lab = bquote("P (mm yr"^"-1"*")"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_P")


ET_predicted_vs_P <- ET_predicted_vs_P +
  geom_smooth(
    data = ale_all_sm |>  filter(var == "P") |>  select(x, ALE_ET_abs_s),
    aes(x = x, y = ALE_ET_abs_s),
    color = "black",
    linetype = "solid",
    linewidth = 0.5,
    se = FALSE,        # no confidence interval shading
    method = "loess",  # or "gam" if you prefer
    span = 0.6         # adjust smoothness (default ~0.75)
  )

ET_predicted_vs_P <- put_line_behind(ET_predicted_vs_P)

# Save the plot
ggsave('../plots/ggplot2/ET_vs_P_ALE.png', plot = ET_predicted_vs_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# ------------------------------------------------------------------------------

source("./src/make_perm_plot.R")

inset_geff <- make_perm_plot(
  csv_path = "../outputs/perm_rf_hybrid_vpd_constrained_fit_summary.csv",
  response = "geff"
)

inset_ET <- make_perm_plot(
  csv_path = "../outputs/perm_rf_hybrid_vpd_constrained_fit_summary.csv",
  response = "ET"
)

# inset_ET <- make_perm_plot_ET_vpd_channels(
#   summary_csv      = "../outputs/perm_rf_hybrid_vpd_constrained_fit_summary.csv",
#   vpd_channels_csv = "../outputs/perm_rf_hybrid_vpd_constrained_VPD_channels_ET_summary.csv"
# )

#--------------
# Combined plot

pA <- p_geff_sqrt_VPD_single_fit + tight
pB <- ET_predicted_vs_VPD_AFFM   + tight
pC <- ET_predicted_vs_Rg         + tight
pD <- ET_predicted_vs_Ta         + tight
pE <- ET_predicted_vs_P          + tight

# 1) Align the base plots (no insets yet)
aligned <- cowplot::align_plots(pA, pB, pC, pD, pE, align = "hv", axis = "tblr")

# 2) Add insets AFTER alignment (wrap the already-aligned grobs)
A_inset <- ggdraw(aligned[[1]]) +
  draw_plot(inset_geff, x = 0.685, y = 0.12, width = 0.30, height = 0.30)

B_inset <- ggdraw(aligned[[2]]) +
  draw_plot(inset_ET,   x = 0.685, y = 0.68, width = 0.30, height = 0.30)

# 3) Build rows (use A_inset and B_inset, and the already-aligned others)
top <- plot_grid(
  A_inset, NULL, B_inset,
  ncol = 3,
  rel_widths = c(1, 0.05, 1),
  labels = c("a", "", "b"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

bottom <- plot_grid(
  aligned[[3]], NULL, aligned[[4]], NULL, aligned[[5]],
  ncol = 5,
  rel_widths = c(1, 0.05, 1, 0.05, 1),
  labels = c("c", "", "d", "", "e"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = -0.1
)

combined <- plot_grid(
  NULL, top, NULL, bottom,
  ncol = 1,
  rel_heights = c(0.1, 3, 0.1, 1.9)
)

ggsave(
  "../plots/ggplot2/g_eff_vs_sqrt_VPD,ET_vs_VPD,Rg,Ta,P_ggplot2_TIDY.png",
  combined, width = 240, height = 192, units = "mm", dpi = RES, bg = "white"
)

# Sort by VPD
jarvis_out |> 
  select(VPD, ET, P, Ta, model, ensemble, PERIOD) |>
  na.omit() |>
  arrange(VPD) |>
  write.csv("../outputs/models_rank.csv", row.names = FALSE)

################################################################################

Rg_min    <- jarvis_bundle$par_hat["Rg_min"]
Rg_max    <- jarvis_bundle$par_hat["Rg_max"]
Ta_min    <- jarvis_bundle$par_hat["Ta_min"]
Ta_max    <- jarvis_bundle$par_hat["Ta_max"]
P_min     <- jarvis_bundle$par_hat["P_min"]
P_max     <- jarvis_bundle$par_hat["P_max"]
k_CO2     <- jarvis_bundle$par_hat["k_CO2"]
b0_VPD    <- jarvis_bundle$par_hat["b0_VPD"]
b1_VPD    <- jarvis_bundle$par_hat["b1_VPD"]

jarvis_out <- jarvis_out |> 
  mutate(
    fRg = pmin(1, pmax(0, (Rg - Rg_min) / (Rg_max - Rg_min))),
    fTa = pmin(1, pmax(0, (Ta - Ta_min) / (Ta_max - Ta_min))),
    fP  = pmin(1, pmax(0, (P  - P_min ) / (P_max - P_min))),
    fVPD = b0_VPD + b1_VPD * 1 / (sqrt(VPD)),
    g_eff_recovered = fP * fRg * fTa * fVPD,
    ET_no_recovered = K_ET * g_eff_recovered
  )

plot(jarvis_out$ET_pred, jarvis_out$ET_no_recovered)

# Now hold all but one functions at their mean and see the response of ET
fP_mean <- mean(jarvis_out$fP, na.rm = TRUE)
fRg_mean <- mean(jarvis_out$fRg, na.rm = TRUE)
fTa_mean <- mean(jarvis_out$fTa, na.rm = TRUE)
fVPD_mean <- mean(jarvis_out$fVPD, na.rm = TRUE)

plot(jarvis_out$g_eff,
     jarvis_out$g_eff_predicted)

plot(jarvis_out$ET,
     jarvis_out$ET_predicted)

plot(jarvis_out$VPD,
     with(jarvis_out, ET))

plot(jarvis_out$VPD,
     with(jarvis_out, K_ET * fP * fRg * fTa * fVPD))

plot(jarvis_out$VPD,
     with(jarvis_out, K_ET * (fP*0 + mean(fP, na.rm = TRUE)) * fRg * fTa * fVPD),
     xlab = "VPD", ylab = "ET")

# Predicted ET vs. observed
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP * fRg * fTa * fVPD, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(0, 1), Y_range_man = c(400-16, 800+16),
                  x_lab = bquote(ET),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_VPD")

make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * (fP*0 + mean(fP, na.rm = TRUE)) * fRg * fTa * fVPD, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04, X_range_man = c(0, 1), Y_range_man = c(400-16, 800+16),
                  x_lab = bquote(ET),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_mean_fP_vs_VPD")

make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP * fRg * fTa * (fVPD*0 + mean(fVPD, na.rm = TRUE)), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_mean_fVPD_vs_VPD")

# VPD effect
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP_mean * fRg_mean * fTa_mean * fVPD, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_VPD_effect_vs_VPD")

# Precipitation effect
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP * fRg_mean * fTa_mean * fVPD_mean, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_P_effect_vs_VPD")

# Temperature effect
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP_mean * fRg_mean * fTa * fVPD_mean, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_Ta_effect_vs_VPD")


# Radiation effect
make_scatter_plot(data = jarvis_out  |>
                    mutate(x = VPD, y = K_ET * fP_mean * fRg * fTa_mean * fVPD_mean, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["eff predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_Rg_effect_vs_VPD")


print(ET_VPD_effect_vs_VPD)
print(ET_P_effect_vs_VPD)
print(ET_Ta_effect_vs_VPD)
print(ET_Rg_effect_vs_VPD)

# ET - VPD effect
plot(jarvis_out$VPD,
     with(jarvis_out, ET - K_ET * fP_mean * fRg_mean * fTa_mean * fVPD),
     title = "VPD effect")

################################################################################

# helper: normalized difference
.normdiff <- function(base, fut, var_name = "") {
  if (is.na(base) || is.na(fut)) return(NA_real_)
  if (var_name == "Bo") return((fut - base) / (1 + base))
  if (base == 0) return(NA_real_)
  (fut - base) / base
}

# 1) Collapse to one row per (label, model) with *_hist and *_fut
num_vars <- names(jarvis_out)[sapply(jarvis_out, is.numeric)]

jarvis_diffs <- jarvis_out |> 
  group_by(model, ensemble) |> 
  summarise(
    across(all_of(num_vars), ~ .[PERIOD == "1981_2005"][1], .names = "{.col}_hist"),
    across(all_of(num_vars), ~ .[PERIOD == "2076_2100"][1], .names = "{.col}_fut"),
    .groups = "drop"
  )

# 2) Add d<var> (absolute) and d<var>_norm (normalized)
for (v in num_vars) {
  hist_col <- paste0(v, "_hist")
  fut_col  <- paste0(v, "_fut")
  d_abs    <- paste0("d", v)
  d_norm   <- paste0("d", v, "_norm")
  
  jarvis_diffs[[d_abs]]  <- jarvis_diffs[[fut_col]] - jarvis_diffs[[hist_col]]
  jarvis_diffs[[d_norm]] <- mapply(.normdiff,
                                   jarvis_diffs[[hist_col]],
                                   jarvis_diffs[[fut_col]],
                                   MoreArgs = list(var_name = v))
}

jarvis_diffs

names(jarvis_diffs)

# Full predicted normalized ET against VPD

# This is a necessary step because of g_eff needs temperature for conversion to ET
a_Ta <- 2.501
b_Ta <- 0.002361

jarvis_diffs <- jarvis_diffs |> 
  mutate(
    eT      = (b_Ta * Ta_hist) / (a_Ta - b_Ta * Ta_hist),
    pT      = 2 * ((b_Ta * Ta_hist) / (a_Ta - b_Ta * Ta_hist))^2
  )

jarvis_diffs <- jarvis_diffs |>
  mutate(
    dET_norm_from_VPD_g_eff_Ta =
      dVPD_norm +                   # First order term
      dg_eff_norm +                 # First order term
      eT * dTa_norm +               # First order term
      dg_eff_norm * dVPD_norm +     # Second order term
      eT * dVPD_norm * dTa_norm +   # Second order term
      eT * dg_eff_norm * dTa_norm + # Second order term
      0.5 * pT * dTa_norm^2         # Second order term
  )

# Now plot normalized ET vs. VPD
plot(jarvis_diffs$dVPD_norm,
     with(jarvis_diffs, dET_norm) )

# It is equivalent because temperature dependence is considered for unit conversion
points(jarvis_diffs$dVPD_norm,
       with(jarvis_diffs, dET_norm_from_VPD_g_eff_Ta),
       col = "red", pch = 3, cex = 0.8)

# If you ignore the temperature dependence, there is more but maybe insignificant scatter
points(jarvis_diffs$dVPD_norm,
       with(jarvis_diffs, dVPD_norm + dg_eff_norm + dg_eff * dVPD / (g_eff_hist * VPD_hist)),
       col = "grey30", pch = 3, cex = 0.8)

# This is of course equivalent
points(jarvis_diffs$dVPD_norm,
       with(jarvis_diffs, dVPD_norm + dg_eff_norm + dg_eff_norm * dVPD_norm),
       col = "blue", pch = 3, cex = 0.3)

plot(jarvis_diffs$dET_pred_norm,
     with(jarvis_diffs,
          dVPD_norm + dg_eff_pred_norm + dg_eff_pred_norm * dVPD_norm))

plot(jarvis_diffs$dET_pred_norm,
     with(jarvis_diffs,
          dVPD_norm + dg_eff_pred_norm + dg_eff_pred_norm * dVPD_norm +
            eT * dTa_norm +               # First order term
            eT * dVPD_norm * dTa_norm +   # Second order term
            eT * dg_eff_norm * dTa_norm + # Second order term
            0.5 * pT * dTa_norm^2) )

################################################################################
# Predicting the conductance and ET

# First order
dg_eff_norm_first_order <- function(dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm) {
  # first order
  FO <- dfRg_norm + dfTa_norm + dfP_norm + dfVPD_norm
  
  FO
}

# Second order cross-terms
dg_eff_norm_second_order <- function(dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm) {
  # first order
  FO <- dfRg_norm + dfTa_norm + dfP_norm + dfVPD_norm
  
  # second-order cross terms (all pairwise products, no squares)
  SO <- 
    dfRg_norm*dfTa_norm +
    dfRg_norm*dfP_norm +
    dfRg_norm*dfVPD_norm +
    
    dfTa_norm*dfP_norm +
    dfTa_norm*dfVPD_norm +
    
    dfP_norm*dfVPD_norm
  
  FO + SO
}

dg_eff_norm_pred_from_FO <- with(jarvis_diffs,
                                      dg_eff_norm_first_order(dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm))

dg_eff_norm_pred_from_SO <- with(jarvis_diffs,
                                      dg_eff_norm_second_order(dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm))

plot(jarvis_diffs$dg_eff_pred_norm, dg_eff_norm_pred_from_FO)
plot(jarvis_diffs$dg_eff_pred_norm, dg_eff_norm_pred_from_SO)

# Exact match
mean(jarvis_diffs$dg_eff_pred_norm - dg_eff_norm_pred_from_SO, na.rm = TRUE)

plot(jarvis_diffs$dET_norm, jarvis_diffs$dVPD_norm + jarvis_diffs$dg_eff_norm + jarvis_diffs$dVPD_norm * jarvis_diffs$dg_eff_norm)
plot(jarvis_diffs$dET_pred_norm, jarvis_diffs$dVPD_norm + jarvis_diffs$dg_eff_pred_norm + jarvis_diffs$dVPD_norm * jarvis_diffs$dg_eff_pred_norm)

plot(jarvis_diffs$dET_pred_norm, jarvis_diffs$dVPD_norm + dg_eff_norm_pred_from_SO + jarvis_diffs$dVPD_norm * dg_eff_norm_pred_from_SO)

dET_norm_third_order <- function(dVPD_norm, dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm) {
  # 1st + 2nd order
  FO <- dVPD_norm + dfRg_norm + dfTa_norm + dfP_norm + dfVPD_norm
  SO <- dVPD_norm*dfRg_norm + dVPD_norm*dfTa_norm + dVPD_norm*dfP_norm + dVPD_norm*dfVPD_norm +
    dfRg_norm*dfTa_norm + dfRg_norm*dfP_norm + dfRg_norm*dfVPD_norm +
    dfTa_norm*dfP_norm + dfTa_norm*dfVPD_norm +
    dfP_norm*dfVPD_norm
  
  # 3rd order: all distinct triples (C(5,3) = 10)
  TO <- dVPD_norm*dfRg_norm*dfTa_norm +
    dVPD_norm*dfRg_norm*dfP_norm +
    dVPD_norm*dfRg_norm*dfVPD_norm +
    dVPD_norm*dfTa_norm*dfP_norm +
    dVPD_norm*dfTa_norm*dfVPD_norm +
    dVPD_norm*dfP_norm*dfVPD_norm +
    dfRg_norm*dfTa_norm*dfP_norm +
    dfRg_norm*dfTa_norm*dfVPD_norm +
    dfRg_norm*dfP_norm*dfVPD_norm +
    dfTa_norm*dfP_norm*dfVPD_norm
  
  FO + SO + TO
}

dET_norm_pred_from_TO <- with(jarvis_diffs,
                                   dET_norm_third_order(dVPD_norm, dfRg_norm, dfTa_norm, dfP_norm, dfVPD_norm))

plot(jarvis_diffs$dET_pred_norm, dET_norm_pred_from_TO)


################################################################################
# Now remove the modifiers
jarvis_diffs <- jarvis_diffs |>
  arrange(ensemble, model)

jarvis_diffs <- jarvis_diffs |>
  mutate(dETnorm_pred_from_TO = dET_norm_third_order(1*dVPD_norm, 1*dfRg_norm, 1*dfTa_norm, 1*dfP_norm, 1*dfVPD_norm)) |> 
  mutate(ET_fut_pred_from_TO = dETnorm_pred_from_TO * ET_pred_hist + ET_pred_hist) |> 
  mutate(ET_fut_pred_from_TO_and_hist = dETnorm_pred_from_TO * ET_hist + ET_hist)

plot(jarvis_diffs$ET_pred_fut, jarvis_diffs$ET_fut_pred_from_TO)
plot(jarvis_diffs$ET_pred_fut, jarvis_diffs$ET_fut_pred_from_TO_and_hist)

part_1 <- jarvis_diffs |> 
  select(
    ensemble, model,
    ET_pred_from_TO_and_hist = ET_hist,
    P = P_hist
  ) |> 
  mutate(PERIOD = "1981_2005")

part_2 <- jarvis_diffs |> 
  select(
    ensemble, model,
    ET_pred_from_TO_and_hist = ET_fut_pred_from_TO_and_hist,
    P = P_fut
  ) |> 
  mutate(PERIOD = "2076_2100")

merged <- bind_rows(part_1, part_2)

# normalize keys (trim spaces) on both sides, keep only needed column from Data_to_plot$abs
dtp_join <- Data_to_plot$abs |> 
  mutate(across(c(PERIOD, ensemble, model), ~ str_squish(as.character(.)))) |> 
  select(PERIOD, ensemble, model, ETo_FAO56_alfalfa) |> 
  distinct()

merged <- merged |> 
  mutate(across(c(PERIOD, ensemble, model), ~ str_squish(as.character(.)))) |> 
  left_join(dtp_join, by = c("PERIOD", "ensemble", "model")) |> 
  mutate(ENSEMBLE = ensemble, MODEL = model)


out_BC_FAO56_alfalfa_ET_semipredicted <- Budyko_curve(merged,
                                                      pet_col = "ETo_FAO56_alfalfa", et_col = "ET_pred_from_TO_and_hist", p_col = "P",
                                                      X_range = c(0.0, 2.2), Y_range = c(0.0, 1.2),
                                                      Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
                                                      boundary_line_col = "gray45",
                                                      boundary_line_type = "solid",
                                                      boundary_line_size = 0.2,
                                                      plot = TRUE, plot_name = "../plots/ggplot2/Budyko_curve_ET_semipredicted_ggplot2.png")

################################################################################
# Now remove fVPD
dETnorm_pred_from_TO_no_VPD <- with(jarvis_diffs,
                                         dET_norm_third_order(dVPD_norm, dfRg_norm, dfTa_norm, dfP_norm,
                                                              0 * dfVPD_norm))

plot(jarvis_diffs$dET_pred_norm, dETnorm_pred_from_TO_no_VPD)


ET_future_no_fVPD <- dETnorm_pred_from_TO_no_VPD *
  jarvis_diffs$ET_pred_hist + jarvis_diffs$ET_pred_hist

plot(jarvis_diffs$ET_pred_fut, ET_future_no_fVPD)

plot(jarvis_diffs$VPD_fut, ET_future_no_fVPD)

jarvis_diffs$ET_future_no_fVPD <- ET_future_no_fVPD

test_data <- Data_to_plot$abs |> filter(PERIOD == "2076_2100")

test_data_sorted <- test_data |> 
  arrange(ensemble, model)

jarvis_diffs_sorted <- jarvis_diffs |> 
  arrange(ensemble, model)

plot(test_data_sorted$ET, jarvis_diffs_sorted$ET_fut)

plot(test_data_sorted$VPD, test_data_sorted$ETo_FAO56_alfalfa)
points(test_data_sorted$VPD, jarvis_diffs_sorted$ET_future_no_fVPD, col = "red")

plot(test_data_sorted$ETo_FAO56_alfalfa,
     jarvis_diffs_sorted$ET_future_no_fVPD)

plot(test_data_sorted$ETo_FAO56_grass,
     jarvis_diffs_sorted$ET_future_no_fVPD)

lm(jarvis_diffs_sorted$ET_future_no_fVPD~test_data_sorted$ETo_FAO56_alfalfa - 1)
lm(jarvis_diffs_sorted$ET_fut~test_data_sorted$ETo_FAO56_alfalfa - 1)

plot(test_data_sorted$ETo_FAO56_alfalfa / test_data_sorted$P,
     jarvis_diffs_sorted$ET_future_no_fVPD / test_data_sorted$P)

# Predicted ET vs. VPD
make_scatter_plot(data = jarvis_out |>
                    select(VPD, ET_pred, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "ET_pred"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD),  y_lab = bquote("ET"["predicted"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_VPD")

# Save the plot
ggsave('../plots/ggplot2/ET_predicted_vs_VPD.png', plot = ET_predicted_vs_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# 1) Keep just the keys + required variables from each table
part_1 <- annual_stats_wide |> 
  select(PERIOD, ENSEMBLE, MODEL, ETo_FAO56_alfalfa = ETo_FAO56_alfalfa, ET = ET, P = P)

part_2 <- jarvis_out |> 
  select(PERIOD, ensemble, model, ET_pred = ET_pred, ET = ET, P = P)

# 2) Inner-join on the matching keys:
part_join <- part_1 |> 
  inner_join(part_2,
             by = c("PERIOD" = "PERIOD",
                    "ENSEMBLE" = "ensemble",
                    "MODEL"    = "model"))

View(part_join)

out_BC_FAO56_alfalfa_ET_predicted <- Budyko_curve(part_join,
                                                  pet_col = "ETo_FAO56_alfalfa", et_col = "ET_pred", p_col = "P.x",
                                                  X_range = c(0.0, 2.2), Y_range = c(0.0, 1.2),
                                                  Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
                                                  boundary_line_col = "gray45",
                                                  boundary_line_type = "solid",
                                                  boundary_line_size = 0.2,
                                                  plot = TRUE, plot_name = "../plots/ggplot2/Budyko_curve_ET_predicted_ggplot2.png")

###############
# Combined plot
###############

#---- Subplot labels
plots <- list(p1 = p_AI_VPD_norm, p2 = p_EI_VPD_norm, p3 = p_geff_sqrt_VPD_single_fit, p4 = p_ET_VPD)
labels <- c("a", "b", "c", "d")

plots <- Map(function(plot, label) {
  plot + annotation_custom(
    textGrob(label, x = unit(0.05, "npc"), y = unit(0.96, "npc"),
             gp = gpar(fontsize = 12))
  )
}, plots, labels)

p1 <- plots$p1; p2 <- plots$p2; p3 <- plots$p3; p4 <- plots$p4

#---- Combine and save
panel_figure <- (p1 | p2) / (p3 | p4)

ggsave('../plots/ggplot2/panel_fig_EI,AI,geff&VPD_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

############
# P vs. VPD
make_scatter_plot(data = Data_to_plot$abs |>
                    select(VPD, P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "P"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("P (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  # plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_P_VPD")
# Save the plot
ggsave('../plots/ggplot2/P_vs_VPD_ggplot2_TIDY.png', plot = p_P_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


############
# EI vs. VPD
make_scatter_plot(data = Data_to_plot$abs |>
                    select(VPD, EI, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "EI"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("EI = ET/P"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_VPD")

# Normalize by mean precipitation form CMIP5, CMIP6 and RCMs
P_avg <- annual_stats |> 
  filter(VAR == "P", PERIOD == "1981_2005", ENSEMBLE %in% c("CMIP5", "CMIP6", "EUR-44")) |> 
  group_by(ENSEMBLE) |> 
  summarise(P_avg = mean(annual_stat, na.rm = TRUE)) |>  
  pull(P_avg) |> 
  mean()  # Overall mean across the 3 ensemble means

lines_df <- lines_df |> 
  mutate(EI_mean_P = predicted_ET / P_avg)

p_EI_VPD <- p_EI_VPD +
  geom_line(
    data = lines_df,
    mapping = aes(x = VPD, y = EI_mean_P, group = fit, linetype = fit),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = peak_log,  colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = peak_sqrt, colour = "grey35", linetype = "solid",  linewidth = 0.25) +
  # Make sure the panel shows everything
  expand_limits(
    x = range(lines_df$VPD, na.rm = TRUE),
    y = range(lines_df$EI_mean_P, na.rm = TRUE)
  ) +
  # If the base plot also uses linetype, make sure this scale includes *both* sets of values:
  scale_linetype_manual(
    values = c("log(VPD)" = "dashed", "sqrt(VPD)" = "solid"),
    labels = c("log(VPD)" = bquote(log(VPD)),
               "sqrt(VPD)" = bquote(1/sqrt(VPD)))
  ) +
  labs(linetype = NULL) +   # remove legend title
  theme(
    legend.position = c(0.86, 1.01),        # inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, color = NA)  # no square
  )

p_EI_VPD <- p_EI_VPD +
  annotate("point", x = 0.66, y = 0.85, shape = 21, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.70, y = 0.85, label = "1981–2005", size = 3, hjust = 0)

p_EI_VPD <- p_EI_VPD +
  annotate("point", x = 0.66, y = 0.82, shape = 24, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.70, y = 0.82, label = "2076–2100 RCP8.5", size = 3, hjust = 0)


# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_EI_VPD$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_EI_VPD$layers <- c(p_EI_VPD$layers[-idx_pts_txt], p_EI_VPD$layers[idx_pts_txt])
}

# Save the plot
ggsave('../plots/ggplot2/EI_vs_VPD_ggplot2_TIDY.png', plot = p_EI_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#######################
# EI vs. VPD arrow plot
arrow_data <- Data_to_plot$abs  |> 
  filter(PERIOD %in% c("1981_2005", "2076_2100"),
         !is.na(VPD), !is.na(EI)) |>
  mutate(PERIOD = factor(PERIOD, levels = c("1981_2005", "2076_2100"))) |>
  group_by(ensemble, model) |>
  arrange(PERIOD, .by_group = TRUE) |>
  summarise(
    x     = first(VPD),
    y     = first(EI),
    xend  = last(VPD),
    yend  = last(EI),
    color = first(color),  # color comes from your tibble
    .groups = "drop"
  ) |>
  filter(!is.na(x), !is.na(y), !is.na(xend), !is.na(yend)) |>
  mutate(
    Scenario = ensemble,
    is_driven_by_key_gcm = str_detect(toupper(model), paste0(driving_GCMs, collapse = "|")),
    arrow_type = if_else(Scenario == "RCMs" | is_driven_by_key_gcm, "thick", "thin"),
    line_size  = if_else(arrow_type == "thick", 0.8, 0.3)
  )

p_EI_VPD_arrow <- ggplot(arrow_data) +
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend,
        colour = color, size = arrow_type),
    lineend = "round",
    arrow = arrow(length = unit(3.5, "mm"), type = "closed")
  ) +
  # identity scales applied once, as in your style
  scale_color_identity() +
  scale_fill_identity() +
  scale_linetype_identity() +
  scale_shape_identity() +
  scale_size_manual(values = c(thin = 0.3, thick = 0.8), guide = "none") +
  labs(
    x = "VPD (kPa)",
    y = "EI = ET/P",
  ) +
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Save the plot
ggsave('../plots/ggplot2/EI_vs_VPD_arrow_ggplot2_TIDY.png', plot = p_EI_VPD_arrow, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-----------------------------------------------
# Examining the scatter in g_eff to VPD relation

Data_to_plot$abs |> names()

Data_to_plot$abs <- Data_to_plot$abs |> 
  mutate(g_eff_pred = coef(fit_sqrt)[1] + coef(fit_sqrt)[2] * 1 / sqrt(VPD))

# Data_to_plot$abs <- Data_to_plot$abs |> 
#   mutate(g_eff_pred = coef(fit_log)[1] + coef(fit_log)[2] * log(VPD))

# Data_to_plot$abs <- Data_to_plot$abs |> 
#   mutate(g_eff_pred = g_eff)

# Predicted vs. original ET
make_scatter_plot(data = Data_to_plot$abs |>
                    select(g_eff, g_eff_pred, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "g_eff", y = "g_eff_pred"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(g["eff"]~"(mm/s)"),  y_lab = bquote(g["eff predicted"]~"(mm/s)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="geff_predicted_vs_geff")

# Save the plot
ggsave('../plots/ggplot2/geff_predicted_vs_geff.png', plot = geff_predicted_vs_geff, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

Data_to_plot$abs <- Data_to_plot$abs |> 
  mutate(ET_pred = (rhoAir * CpAir / gamma) * VPD /
           (1 / (g_eff_pred / 1000)) *                                       # Surface resistance
           1 / ((2.501 * (10^6) - 2361 * Ta) / (10^6)) * 3600*24 / 10^6 * 365.25  # From LE to annual ET
  )

# Predicted vs. original ET
make_scatter_plot(data = Data_to_plot$abs |>
                    select(ET, ET_pred, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "ET", y = "ET_pred"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(ET~"(mm/yr)"),  y_lab = bquote(ET[predicted]~"(mm/yr)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="ET_predicted_vs_ET")

# Save the plot
ggsave('../plots/ggplot2/ET_predicted_vs_ET.png', plot = ET_predicted_vs_ET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#--------------------
# Check the residuals
Data_to_plot$abs <- Data_to_plot$abs |> 
  mutate(g_eff_resids = g_eff - g_eff_pred)

Data_to_plot$abs |> pull(g_eff_resids) |> plot()

make_scatter_plot(data = Data_to_plot$abs |>
                    select(VPD, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_VPD.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Global radiation
make_scatter_plot(data = Data_to_plot$abs |>
                    select(Rg, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Rg", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(R[g]),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_Rg.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Air temperature
make_scatter_plot(data = Data_to_plot$abs |>
                    select(Ta, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "Ta", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(T[a]),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_Ta.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Precipitation
make_scatter_plot(data = Data_to_plot$abs |>
                    select(P, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "P", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(P),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_P.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------
# Merging tibles

keys <- c("PERIOD", "ensemble", "model")
aesthetic_cols <- c("color", "fill", "border", "shape", "linetype")

# 1) Normalize key types
Data_to_plot$abs <- Data_to_plot$abs |> 
  mutate(across(all_of(keys), as.character))

Data_to_plot_right <- Data_to_plot$diffs |> 
  mutate(PERIOD = "2076_2100") |>               # <- add the constant period
  mutate(across(all_of(keys), as.character)) |> 
  select(-any_of(aesthetic_cols)) |>            # avoid collisions with left-side aesthetics
  distinct(across(all_of(keys)), .keep_all = TRUE)

# 2) Guard against duplicate keys on the right
dups <- Data_to_plot_right |> 
  count(across(all_of(keys))) |> 
  filter(n > 1)
if (nrow(dups)) {
  stop("`Data_to_plot` has duplicate (PERIOD, model, label) keys; resolve before joining.")
}

# 3) Left-join: only matches 2076_2100 rows; others remain as-is
merged_df <- Data_to_plot$abs |> 
  left_join(Data_to_plot_right, by = keys)

# 4) Sanity checks
stopifnot(nrow(merged_df) == nrow(Data_to_plot$abs))

# Which (PERIOD, model, label) in the left didn’t match?
unmatched <- Data_to_plot$abs |> 
  distinct(across(all_of(keys))) |> 
  anti_join(Data_to_plot_right |>  distinct(across(all_of(keys))), by = keys)
if (nrow(unmatched)) {
  message("Unmatched keys on right (won't be filled):\n",
          paste(capture.output(print(unmatched)), collapse = "\n"))
}

#--------------
# Precipitation
make_scatter_plot(data = merged_df |>
                    select(P, g_eff_resids, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "P", y = "g_eff_resids"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(P),  y_lab = bquote("g"["eff predicted"] - "g"[eff]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="g_eff_resids_relations")

# Save the plot
ggsave('../plots/ggplot2/g_eff_resids_and_P.png', plot = g_eff_resids_relations, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Normalized reference g_eff
a0 = coef(fit_sqrt)[1]
b0 = coef(fit_sqrt)[2]
merged_df <- merged_df |> 
  mutate(d_g_eff_ref_over_g_eff_ref = d_ET_over_ET -
           1 / (2 * VPD) * (b0 + 2 * a0 * sqrt(VPD)) / (b0 + a0 *sqrt(VPD)))

#---------------------------------------
# Normalized aridity index and g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa,
                           d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'AI/'*'AI'),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_d_AI_over_AI")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_d_AI_over_AI.png', plot = d_g_eff_ref_g_eff_ref_vs_d_AI_over_AI, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------------------
# Normalized precipitation and g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_P_over_P, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_P_over_P", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_P_over_P")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_P_over_P.png', plot = d_g_eff_ref_g_eff_ref_vs_P_over_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------
# Precipitation and g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(P, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "P", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_P")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_P.png', plot = d_g_eff_ref_g_eff_ref_vs_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------------------
# Normalized net radiation and g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_Rn_over_Rn, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rn_over_Rn", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"R"[n]*" / R"[n]),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_Rn_over_Rn")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_Rn_over_Rn.png', plot = d_g_eff_ref_g_eff_ref_vs_Rn_over_Rn, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------------------------------
# Normalized global radiation and g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_Rg_over_Rg, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"R"[g]*" / R"[g]),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_Rg_over_Rg")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_Rg_over_Rg.png', plot = d_g_eff_ref_g_eff_ref_vs_Rg_over_Rg, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------------------------------
# Change in air temperature and normalized g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_Ta, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"T"[a]*"(°C)"),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_d_Ta")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_d_Ta.png', plot = d_g_eff_ref_g_eff_ref_vs_d_Ta, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------
# Change in vapor pressure deficit and normalized g_eff_ref
make_scatter_plot(data = merged_df |>
                    select(d_VPD_over_VPD, d_g_eff_ref_over_g_eff_ref, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_g_eff_ref_over_g_eff_ref") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"T"[a]*"(°C)"),  y_lab = bquote(Delta*'g'["eff ref"]*" / g"["eff ref"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_g_eff_ref_g_eff_ref_vs_d_VPD_over_VPD")

# Save the plot
ggsave('../plots/ggplot2/d_g_eff_ref_g_eff_ref_vs_d_VPD_over_VPD.png', plot = d_g_eff_ref_g_eff_ref_vs_d_VPD_over_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------------------
# Normalized net radiation and precipitation
make_scatter_plot(data = merged_df |>
                    select(d_Rn_over_Rn, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rn_over_Rn", y = "d_P_over_P") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"R"[n]*" / R"[n]),  y_lab = bquote(Delta*"P / P"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_P_over_P_vs_Rn_over_Rn")

# Save the plot
ggsave('../plots/ggplot2/d_P_over_P_vs_Rn_over_Rn.png', plot = d_P_over_P_vs_Rn_over_Rn, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------
# Normalized global radiation and precipitation
make_scatter_plot(data = merged_df |>
                    select(d_Rg_over_Rg, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_P_over_P") |> 
                    mutate(shape = 21),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"R"[g]*" / R"[g]),  y_lab = bquote(Delta*"P / P"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_P_over_P_vs_Rg_over_Rg")

# Save the plot
ggsave('../plots/ggplot2/d_P_over_P_vs_Rg_over_Rg.png', plot = d_P_over_P_vs_Rg_over_Rg, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

############
# EI vs. VPD
make_scatter_plot(data = Data_to_plot$abs |>
                    select(VPD, EI, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "EI"),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("EI = ET/P"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_VPD")

# Normalize by mean precipitation form CMIP5, CMIP6 and RCMs
P_avg <- annual_stats |> 
  filter(VAR == "P", PERIOD == "1981_2005", ENSEMBLE %in% c("CMIP5", "CMIP6", "EUR-44")) |> 
  group_by(ENSEMBLE) |> 
  summarise(P_avg = mean(annual_stat, na.rm = TRUE)) |>  
  pull(P_avg) |> 
  mean()  # Overall mean across the 3 ensemble means

lines_df <- lines_df |> 
  mutate(EI_mean_P = predicted_ET / P_avg)

p_EI_VPD <- p_EI_VPD +
  geom_line(
    data = lines_df,
    mapping = aes(x = VPD, y = EI_mean_P, group = fit, linetype = fit),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = peak_log,  colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = peak_sqrt, colour = "grey35", linetype = "solid",  linewidth = 0.25) +
  # Make sure the panel shows everything
  expand_limits(
    x = range(lines_df$VPD, na.rm = TRUE),
    y = range(lines_df$EI_mean_P, na.rm = TRUE)
  ) +
  # If the base plot also uses linetype, make sure this scale includes *both* sets of values:
  scale_linetype_manual(
    values = c("log(VPD)" = "dashed", "sqrt(VPD)" = "solid"),
    labels = c("log(VPD)" = bquote(log(VPD)),
               "sqrt(VPD)" = bquote(1/sqrt(VPD)))
  ) +
  labs(linetype = NULL) +   # remove legend title
  theme(
    legend.position = c(0.86, 1.01),        # inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, color = NA)  # no square
  )

p_EI_VPD <- p_EI_VPD +
  annotate("point", x = 0.66, y = 0.85, shape = 21, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.70, y = 0.85, label = "1981–2005", size = 3, hjust = 0)

p_EI_VPD <- p_EI_VPD +
  annotate("point", x = 0.66, y = 0.82, shape = 24, size = 2, stroke = 0.8,
           color = "#2b2b2b", fill = NA) +
  annotate("text", x = 0.70, y = 0.82, label = "2076–2100 RCP8.5", size = 3, hjust = 0)


# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_EI_VPD$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_EI_VPD$layers <- c(p_EI_VPD$layers[-idx_pts_txt], p_EI_VPD$layers[idx_pts_txt])
}

# Save the plot
ggsave('../plots/ggplot2/EI_vs_VPD_ggplot2_TIDY.png', plot = p_EI_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

################################################################################

# With normalized VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_ET_over_ET_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_ET_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# With normalized VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'PET/'*'PET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_PET_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_PET_over_PET_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_PET_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# ---------------------------------------------------------
# Budyko curve components against VPD perturbation analysis

# With normalized VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_EI_over_EI, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_EI_over_EI"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'EI/'*'EI'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_EI_over_EI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_EI_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# With normalized VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_AI_FAO56_alfalfa_over_AI_FAO56_alfalfa"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'AI/'*'AI'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_AI_over_AI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_AI_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_AI_ETo_FAO56_alfalfa_GCM_CO2_corr_over_AI_ETo_FAO56_alfalfa_GCM_CO2_corr, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_AI_ETo_FAO56_alfalfa_GCM_CO2_corr_over_AI_ETo_FAO56_alfalfa_GCM_CO2_corr"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*AI[CO[2]]/AI),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_CO2corr_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_AI_CO2_corr_over_AI_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_AI_CO2corr_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(dET_to_PET_ratio_over_ET_to_PET_ratio = (d_ET_over_ET - d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa) / 
           (1 + d_ETo_FAO56_alfalfa_over_ETo_FAO56_alfalfa))

# With normalized VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, dET_to_PET_ratio_over_ET_to_PET_ratio, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "dET_to_PET_ratio_over_ET_to_PET_ratio"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'VPD/VPD'),  y_lab = bquote(Delta*'(ET/PET) / (ET/PET)'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_over_PET_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_ET_to_PET_ratio_over_ET_to_PET_ratio_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_ET_over_PET_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD_over_VPD, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD_over_VPD", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD/'*'VPD'),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD_norm")
# Save the plot
ggsave('../plots/ggplot2/delta_P_over_P_vs_delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p_P_VPD_norm, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# With absolute values of VPD
make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD (kPa)'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_norm_VPD")
# Save the plot
ggsave('../plots/ggplot2/delta_ET_over_ET_vs_delta_VPD_ggplot2_TIDY.png', plot = p_ET_norm_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_VPD, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_VPD", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'VPD (kPa)'),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_norm_VPD")
# Save the plot
ggsave('../plots/ggplot2/delta_P_over_P_vs_delta_VPD_ggplot2_TIDY.png', plot = p_P_norm_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


# With absolute values of VPD in the historical run
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(VPD = d_VPD / d_VPD_over_VPD)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD, d_ET_over_ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_norm_VPD_historical")
# Save the plot
ggsave('../plots/ggplot2/delta_ET_over_ET_vs_VPD_historical_ggplot2_TIDY.png', plot = p_ET_norm_VPD_historical, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD, d_P_over_P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_norm_VPD_historical")
# Save the plot
ggsave('../plots/ggplot2/delta_P_over_P_vs_VPD_historical_ggplot2_TIDY.png', plot = p_P_norm_VPD_historical, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# With absolute values of VPD, ET and P in the historical run
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(VPD = d_VPD / d_VPD_over_VPD)
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(P = d_P / d_P_over_P)
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(ET = d_ET / d_ET_over_ET)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD, ET, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(ET["1981–2005"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_VPD_historical")
# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_historical_ggplot2_TIDY.png', plot = p_ET_VPD_historical, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD, P, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD", y = "P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(P["1981–2005"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD_historical")
# Save the plot
ggsave('../plots/ggplot2/P_vs_VPD_historical_ggplot2_TIDY.png', plot = p_P_VPD_historical, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# With absolute values of VPD, ET and P in the RCP run
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(VPD_RCP = d_VPD + VPD)
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(P_RCP = d_P + P)
Data_to_plot$diffs <- Data_to_plot$diffs |>
  mutate(ET_RCP = d_ET + ET)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD_RCP, ET_RCP, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD_RCP", y = "ET_RCP"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["2076–2100"]~"(kPa)"),  y_lab = bquote(ET["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/ET_vs_VPD_RCP_ggplot2_TIDY.png', plot = p_ET_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = Data_to_plot$diffs |>
                    select(VPD_RCP, P_RCP, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "VPD_RCP", y = "P_RCP"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["2076–2100"]~"(kPa)"),  y_lab = bquote(P["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/P_vs_VPD_RCP_ggplot2_TIDY.png', plot = p_P_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

################################################################################
# Emergent constraint
plot(jarvis_diffs$VPD_hist, jarvis_diffs$VPD_fut)
plot(jarvis_diffs$VPD_hist, (jarvis_diffs$VPD_fut - jarvis_diffs$VPD_hist))
plot(jarvis_diffs$VPD_hist, (jarvis_diffs$VPD_fut - jarvis_diffs$VPD_hist) / jarvis_diffs$VPD_hist)

plot(jarvis_diffs$ET_hist / jarvis_diffs$P_hist, jarvis_diffs$ET_fut / jarvis_diffs$P_fut)
plot(jarvis_diffs$ET_hist, jarvis_diffs$ET_fut)
plot(jarvis_diffs$P_hist, jarvis_diffs$P_fut)
plot(jarvis_diffs$VPD_hist, jarvis_diffs$P_fut)
plot(jarvis_diffs$VPD_hist, jarvis_diffs$g_eff_fut)
plot(jarvis_diffs$VPD_hist, (jarvis_diffs$g_eff_fut - jarvis_diffs$g_eff_hist) / jarvis_diffs$g_eff_hist)
plot(jarvis_diffs$VPD_hist, jarvis_diffs$ET_fut)


plot(jarvis_diffs$Ta_hist, jarvis_diffs$Ta_fut)
plot(jarvis_diffs$P_hist, jarvis_diffs$P_fut)
plot(jarvis_diffs$ET_hist, jarvis_diffs$ET_fut)

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = EI_hist, y = EI_fut, ensemble = interaction(ensemble, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(EI["1981–2005"]),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_EI_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/EI_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_EI_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "VPD"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_VPD_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "VPD", "P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = 1000 * VPD_hist / P_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("1000 * VPD / P"["1981–2005"]~"(kPa / mm yr)"),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_P_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_P_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_VPD_P_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "e"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = e_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(e["1981–2005"]~"(kPa)"),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_e_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/e_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_e_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "e_sat"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = e_sat_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(e["sat 1981–2005"]~"(kPa)"),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_e_sat_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/e_sat_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_e_sat_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "ET"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(ET["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/ET_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_ET_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI", "ET", "ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist/ETo_FAO56_alfalfa_hist , y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET / PET"["1981–2005"]),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_EI_RCP")
# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_EI_RCP_ggplot2_TIDY.png', plot = p_ET2PET_EI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("AI_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = AI_FAO56_alfalfa_hist, y = AI_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(AI["1981–2005"]),  y_lab = bquote(AI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_AI_RCP")
# Save the plot
ggsave('../plots/ggplot2/AI_hist_vs_AI_RCP_ggplot2_TIDY.png', plot = p_AI_AI_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "AI_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = AI_FAO56_alfalfa_fut - AI_FAO56_alfalfa_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]),  y_lab = bquote(Delta*AI),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_dAI_VPD_hist")
# Save the plot
ggsave('../plots/ggplot2/delta_AI_vs_VPD_hist_ggplot2_TIDY.png', plot = p_dAI_VPD_hist, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "EI"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = EI_fut - EI_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]),  y_lab = bquote(Delta*EI),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_dEI_VPD_hist")
# Save the plot
ggsave('../plots/ggplot2/delta_EI_vs_VPD_hist_ggplot2_TIDY.png', plot = p_dEI_VPD_hist, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = ET_fut/ ETo_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("ET/PET"["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_ET2PET_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_ET2PET_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_ET2PET_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = ET_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("ET"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_ET_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_ET_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_ET_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa", "P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = P_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("P"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_P_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_P_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_P_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa", "EI"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = EI_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("EI"["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_EI_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_EI_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_EI_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa", "VPD"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = VPD_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("VPD"["2076–2100"]~"(kPa)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_VPD_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_VPD_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_VPD_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa", "Ta"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist / ETo_FAO56_alfalfa_hist, y = Ta_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET/PET"["1981–2005"]),  y_lab = bquote("T"["a 2076–2100"]~"(°C)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_hist_Ta_fut")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_hist_vs_Ta_fut_ggplot2_TIDY.png', plot = p_ET2PET_hist_Ta_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ETo_FAO56_alfalfa_hist, y = ETo_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("PET"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote("PET"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_PET_hist_PET_fut")

# Save the plot
ggsave('../plots/ggplot2/PET_hist_vs_PET_RCP_ggplot2_TIDY.png', plot = p_PET_hist_PET_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ETo_FAO56_alfalfa_hist, y = (ETo_FAO56_alfalfa_fut - ETo_FAO56_alfalfa_hist) / ETo_FAO56_alfalfa_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("PET"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote(Delta*'PET/'*'PET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_PET_hist_delta_PET")

# Save the plot
ggsave('../plots/ggplot2/PET_hist_vs_delta_PET_ggplot2_TIDY.png', plot = p_PET_hist_delta_PET, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = P_hist, y = (P_fut - P_hist) / P_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("P"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_hist_delta_P")

# Save the plot
ggsave('../plots/ggplot2/P_hist_vs_delta_P_ggplot2_TIDY.png', plot = p_P_hist_delta_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = (P_fut - P_hist) / P_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD"["1981–2005"]~"(kPa)"),  y_lab = bquote(Delta*'P/'*'P'),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_hist_delta_P")

# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_delta_P_ggplot2_TIDY.png', plot = p_VPD_hist_delta_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist, y = ET_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(ET["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote(ET["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_ET_RCP")
# Save the plot
ggsave('../plots/ggplot2/ET_hist_vs_ET_RCP_ggplot2_TIDY.png', plot = p_ET_ET_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = P_hist, y = P_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(P["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote(P["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_P_RCP")
# Save the plot
ggsave('../plots/ggplot2/P_hist_vs_P_RCP_ggplot2_TIDY.png', plot = p_P_P_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = VPD_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(VPD["2076–2100"]~(kPa)),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_VPD_RCP_ggplot2_TIDY.png', plot = p_VPD_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs  |> filter(ensemble  %in% c("CMIP5", "EUR-44")),
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = VPD_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(VPD["2076–2100"]~(kPa)),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_VPD_RCP_CMIP&EUR-44-ggplot2_TIDY.png', plot = p_VPD_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "ET"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = ET_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(ET["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_ET_RCP_ggplot2_TIDY.png', plot = p_ET_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "AI_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = AI_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(AI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_AI_RCP_ggplot2_TIDY.png', plot = p_AI_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = P_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(P["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_P_VPD_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_P_RCP_ggplot2_TIDY.png', plot = p_P_VPD_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("g_eff"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = g_eff_hist, y = g_eff_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(g["eff 1981–2005"]~"(mm s"^"-1"*")"),  y_lab = bquote(g["eff 2076–2100"]~"(mm s"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_g_eff_g_eff_RCP")
# Save the plot
ggsave('../plots/ggplot2/g_eff_hist_vs_g_eff_RCP_ggplot2_TIDY.png', plot = p_g_eff_g_eff_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "g_eff"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = g_eff_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~"(kPa)"),  y_lab = bquote(g["eff 2076–2100"]~"(mm s"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_VPD_g_eff_RCP")
# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_g_eff_RCP_ggplot2_TIDY.png', plot = p_VPD_g_eff_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("Ta"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = Ta_hist, y = Ta_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(T["a 1981–2005"]~"(°C)"),  y_lab = bquote(T["a 2076–2100"]~"(°C)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_Ta_Ta_RCP")
# Save the plot
ggsave('../plots/ggplot2/Ta_hist_vs_Ta_RCP_ggplot2_TIDY.png', plot = p_Ta_Ta_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("RH"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = RH_hist, y = RH_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(RH["1981–2005"]~"(%)"),  y_lab = bquote(RH["2076–2100"]~"(%)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = FALSE,
                  save_ggplot2_obj_as="p_RH_RH_RCP")
# Save the plot
ggsave('../plots/ggplot2/RH_hist_vs_RH_RCP_ggplot2_TIDY.png', plot = p_RH_RH_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# VPD vs. Delta VPD over VPD
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "n"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist,
                           y = (VPD_fut - VPD_hist) / VPD_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = "VPD (kPa)",  y_lab = bquote(Delta*"VPD / VPD"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="VPD_hist_vs_d_VPD_over_VPD")

# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_d_VPD_over_VPD_ggplot2_TIDY.png', plot = VPD_hist_vs_d_VPD_over_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# VPD vs. Delta VPD
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "n"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist,
                           y = (VPD_fut - VPD_hist), model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = "VPD (kPa)",  y_lab = bquote(Delta*"VPD (kPa)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="VPD_hist_vs_d_VPD")

# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_d_VPD_ggplot2_TIDY.png', plot = VPD_hist_vs_d_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Test log(VPD_fut/VPD_hist) vs Delta Ta

LM_eq_labels <- tibble(
  x = rep(0.03, 6),
  y = rep(1, 6) - seq(from = 0.12, to = 0.5, length.out = 6)
)
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "Ta"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = Ta_fut - Ta_hist, y = log(VPD_fut / VPD_hist), model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  LM_eq_labels = LM_eq_labels, force_origin = FALSE,
                  x_lab = bquote(Delta*T["a"]~"(°C)"),  y_lab = bquote("ln(VPD/VPD)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_log_VPD_fut_VPD_hist_delta_Ta")
# Save the plot
ggsave('../plots/ggplot2/log_VPD_fut_VPD_hist_delta_Ta_ggplot2_TIDY.png', plot = p_log_VPD_fut_VPD_hist_delta_Ta, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

LM_eq_labels <- tibble(
  x = rep(0.03, 6),
  y = rep(1, 6) - seq(from = 0.12, to = 0.5, length.out = 6)
)
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "Ta"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = Ta_fut - Ta_hist, y = (VPD_fut - VPD_hist) / VPD_hist, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  LM_eq_labels = LM_eq_labels, force_origin = FALSE,
                  x_lab = bquote(Delta*T["a"]~"(°C)"),  y_lab = bquote(Delta*"VPD / VPD"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_delta_VPD_norm_delta_Ta")
# Save the plot
ggsave('../plots/ggplot2/delta_VPD_norm_delta_Ta_ggplot2_TIDY.png', plot = p_delta_VPD_norm_delta_Ta, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("Bo"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = Bo_hist, y = Bo_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Bo["1981–2005"]),  y_lab = bquote(Bo["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_Bo_Bo_RCP")
# Save the plot
ggsave('../plots/ggplot2/Bo_hist_vs_Bo_RCP_ggplot2_TIDY.png', plot = p_Bo_Bo_RCP, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
################################################################################
# Emergent constraint – single fit

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  ensemble= c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("EI"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = EI_hist, y = EI_fut, ensemble = interaction(ensemble, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(EI["1981–2005"]),  y_lab = bquote(EI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p_EI_EI_RCP_single_fit")

# 1) Extract the first point layer from `p`
idx_pt <- which(vapply(p_EI_EI_RCP_single_fit$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
stopifnot(!is.na(idx_pt))  # ensure we found a point layer

pts <- layer_data(p_EI_EI_RCP_single_fit, idx_pt)
pts <- pts[is.finite(pts$x) & is.finite(pts$y), ]  # keep valid points

# 2) Fit OLS on the plotted coordinates (same space as the plot)
fit <- lm(y ~ x, data = pts)
b0  <- unname(coef(fit)[1])
b1  <- unname(coef(fit)[2])
r2  <- summary(fit)$r.squared

# Build a pretty, signed label with 2 decimals
sign <- if (b1 >= 0) "+" else "−"
eq_label <- sprintf("y = %.2f %s %.2fx; R\u00B2 = %.2f",
                    b0, sign, abs(b1), r2)

# 3) Add line + 95% CI and the equation label (top-left)
plot_data <- ggplot_build(p_EI_EI_RCP_single_fit)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

p_EI_EI_RCP_single_fit <- p_EI_EI_RCP_single_fit +
  geom_smooth(
    data = pts,
    aes(x = x, y = y),
    method = "lm", formula = y ~ x,
    se = TRUE,
    inherit.aes = FALSE,
    color = "black", fill = "grey40", alpha = 0.25, linewidth = 0.6
  ) +
  annotate("text",
           x = rescale(0.54, X_range), y = rescale(0.05, Y_range), label = eq_label,
           hjust = 0, vjust = 0.5, size = 3.6)

# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_EI_EI_RCP_single_fit$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_EI_EI_RCP_single_fit$layers <- c(p_EI_EI_RCP_single_fit$layers[-idx_pts_txt], p_EI_EI_RCP_single_fit$layers[idx_pts_txt])
}

# Save the plot
ggsave('../plots/ggplot2/EI_hist_vs_EI_RCP_single_fit_ggplot2_TIDY.png', plot = p_EI_EI_RCP_single_fit, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "AI_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = VPD_hist, y = AI_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(VPD["1981–2005"]~(kPa)),  y_lab = bquote(AI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_VPD_RCP_single_fit")

# 1) Extract the first point layer from `p`
idx_pt <- which(vapply(p_AI_VPD_RCP_single_fit$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
stopifnot(!is.na(idx_pt))  # ensure we found a point layer

pts <- layer_data(p_AI_VPD_RCP_single_fit, idx_pt)
pts <- pts[is.finite(pts$x) & is.finite(pts$y), ]  # keep valid points

# 2) Fit OLS on the plotted coordinates (same space as the plot)
fit <- lm(y ~ x, data = pts)
b0  <- unname(coef(fit)[1])
b1  <- unname(coef(fit)[2])
r2  <- summary(fit)$r.squared

# Build a pretty, signed label with 2 decimals
sign <- if (b1 >= 0) "+" else "−"
eq_label <- sprintf("y = %.2f %s %.2fx; R\u00B2 = %.2f",
                    b0, sign, abs(b1), r2)

# 3) Add line + 95% CI and the equation label (top-left)
plot_data <- ggplot_build(p_AI_VPD_RCP_single_fit)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

p_AI_VPD_RCP_single_fit <- p_AI_VPD_RCP_single_fit +
  geom_smooth(
    data = pts,
    aes(x = x, y = y),
    method = "lm", formula = y ~ x,
    se = TRUE,
    inherit.aes = FALSE,
    color = "black", fill = "grey40", alpha = 0.25, linewidth = 0.6
  ) +
  annotate("text",
           x = rescale(0.54, X_range), y = rescale(0.05, Y_range), label = eq_label,
           hjust = 0, vjust = 0.5, size = 3.6)

# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_AI_VPD_RCP_single_fit$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_AI_VPD_RCP_single_fit$layers <- c(p_AI_VPD_RCP_single_fit$layers[-idx_pts_txt], p_AI_VPD_RCP_single_fit$layers[idx_pts_txt])
}

# Save the plot
ggsave('../plots/ggplot2/VPD_hist_vs_AI_RCP_single_fit_ggplot2_TIDY.png', plot = p_AI_VPD_RCP_single_fit, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("AI_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = AI_FAO56_alfalfa_hist, y = AI_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(AI["1981–2005"]),  y_lab = bquote(AI["2076–2100"]),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_AI_AI_RCP_single_fit")

# 1) Extract the first point layer from `p`
idx_pt <- which(vapply(p_AI_AI_RCP_single_fit$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
stopifnot(!is.na(idx_pt))  # ensure we found a point layer

pts <- layer_data(p_AI_AI_RCP_single_fit, idx_pt)
pts <- pts[is.finite(pts$x) & is.finite(pts$y), ]  # keep valid points

# 2) Fit OLS on the plotted coordinates (same space as the plot)
fit <- lm(y ~ x, data = pts)
b0  <- unname(coef(fit)[1])
b1  <- unname(coef(fit)[2])
r2  <- summary(fit)$r.squared

# Build a pretty, signed label with 2 decimals
sign <- if (b1 >= 0) "+" else "−"
eq_label <- sprintf("y = %.2f %s %.2fx; R\u00B2 = %.2f",
                    b0, sign, abs(b1), r2)

# 3) Add line + 95% CI and the equation label (top-left)
plot_data <- ggplot_build(p_AI_AI_RCP_single_fit)
X_range <- plot_data$layout$panel_params[[1]]$x.range
Y_range <- plot_data$layout$panel_params[[1]]$y.range

p_AI_AI_RCP_single_fit <- p_AI_AI_RCP_single_fit +
  geom_smooth(
    data = pts,
    aes(x = x, y = y),
    method = "lm", formula = y ~ x,
    se = TRUE,
    inherit.aes = FALSE,
    color = "black", fill = "grey40", alpha = 0.25, linewidth = 0.6
  ) +
  annotate("text",
           x = rescale(0.54, X_range), y = rescale(0.05, Y_range), label = eq_label,
           hjust = 0, vjust = 0.5, size = 3.6)

# Identify layers that are GeomPoint or GeomTextRepel
idx_pts_txt <- which(vapply(p_AI_AI_RCP_single_fit$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p_AI_AI_RCP_single_fit$layers <- c(p_AI_AI_RCP_single_fit$layers[-idx_pts_txt], p_AI_AI_RCP_single_fit$layers[idx_pts_txt])
}

# Save the plot
ggsave('../plots/ggplot2/AI_hist_vs_AI_RCP_single_fit_ggplot2_TIDY.png', plot = p_AI_AI_RCP_single_fit, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

with(pts, cor(x, y, use = "complete.obs", method = "pearson"))

################################################################################
# Complementary relations
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist + ETo_FAO56_alfalfa_hist, y = ET_fut + ETo_FAO56_alfalfa_fut , model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET + PET"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote("ET + PET"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_PET_hist_ET_PET_fut")

# Save the plot
ggsave('../plots/ggplot2/p_ET_PET_hist_vs_ET_PET_fut_ggplot2_TIDY.png', plot = p_ET_PET_hist_ET_PET_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "P"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist + P_hist, y = ET_fut + P_fut , model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET + P"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote("ET + P"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_P_hist_ET_P_fut")

# Save the plot
ggsave('../plots/ggplot2/p_ET_P_hist_vs_ET_P_fut_ggplot2_TIDY.png', plot = p_ET_P_hist_ET_P_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("ET", "P", "ETo_FAO56_alfalfa"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = ET_hist + P_hist + ETo_FAO56_alfalfa_hist,
                           y = ET_fut + P_fut + ETo_FAO56_alfalfa_fut, model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("ET + P + PET"["1981–2005"]~"(mm yr"^"-1"*")"),  y_lab = bquote("ET + P + PET"["2076–2100"]~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = TRUE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_P_PET_hist_ET_P_PET_fut")

# Save the plot
ggsave('../plots/ggplot2/p_ET_P_PET_hist_vs_ET_P_PET_fut_ggplot2_TIDY.png', plot = p_ET_P_PET_hist_ET_P_PET_fut, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# Delta n over n vs. Delta VPD over VPD
make_scatter_plot(data = pair_periods(df   = Data_to_plot$abs,
                                      hist = "1981_2005",
                                      fut  = "2076_2100",
                                      vars = c("VPD", "n"),
                                      by   = c("ensemble","model")) |>
                    mutate(x = (VPD_fut - VPD_hist) / VPD_hist,
                           y = (n_fut + 0.72 - (n_hist + 0.72)) / (n_hist + 0.72), model = interaction(model, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*"VPD / VPD"),  y_lab = bquote(Delta*omega~"/"~omega),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="d_omega_over_omega_vs_d_VPD_over_VPD")

# Save the plot
ggsave('../plots/ggplot2/d_omega_over_omega_vs_d_VPD_over_VPD_ggplot2_TIDY.png', plot = d_omega_over_omega_vs_d_VPD_over_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

LM_eq_labels <- tibble(
  x = rep(0.03, 6),
  y = rep(1, 6) - seq(from = 0.12, to = 0.5, length.out = 6)
)

make_scatter_plot(data = Data_to_plot$diffs |>
                    select(d_Ta, d_VPD_over_VPD, ensemble, color, fill, border, shape, model, linetype) |> 
                    rename(x = "d_Ta", y = "d_VPD_over_VPD"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  LM_eq_labels = LM_eq_labels,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'VPD/'*'VPD'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p3")

################################################################################

# Pick variables to pair (keep only those present & numeric)
vars <- c(
  "A","ET","ET_eq","ETo_FAO56_alfalfa","ETo_FAO56_alfalfa_GCM_CO2_corr",
  "ETo_FAO56_grass","ETo_FAO56_grass_GCM_CO2_corr","G","H","LAI","LE",
  "LW_net","P","PET_PT_alpha_1.26","RH","RO","Rg","Rn","SVP","SW_net",
  "Ta","VPD","W_to_mm","e","e_sat","lambda","max_ET_ratio","q","u10",
  "rs_eff","g_eff","Bo","IWV","PE","AI_FAO56_alfalfa",
  "AI_ETo_FAO56_alfalfa_GCM_CO2_corr","EI","ET_over_ETo_FAO56_alfalfa",
  "ET_over_ETo_FAO56_alfalfa_GCM_CO2_corr","PET","n","ω"
)

paired <- pair_periods(
  df   = Data_to_plot$abs,
  hist = "1981_2005",
  fut  = "2076_2100",
  vars = vars,
  by   = c("ensemble","model")
)

hist_cols <- paste0(vars, "_hist")
fut_cols  <- paste0(vars, "_fut")

# Keep only columns that actually exist & are numeric
hist_cols <- intersect(hist_cols, names(paired))
fut_cols  <- intersect(fut_cols,  names(paired))

X <- as.matrix(paired[, hist_cols, drop = FALSE])
Y <- as.matrix(paired[, fut_cols,  drop = FALSE])

# Optional: drop columns with <2 finite values or zero variance
ok_hist <- colSums(is.finite(X)) > 1 & apply(X, 2, sd, na.rm = TRUE) > 0
ok_fut  <- colSums(is.finite(Y)) > 1 & apply(Y, 2, sd, na.rm = TRUE) > 0
X <- X[, ok_hist, drop = FALSE]
Y <- Y[, ok_fut,  drop = FALSE]

# Cross correlation matrix (pairwise NA handling -> matches vector-by-vector)
R <- cor(X, Y, use = "pairwise.complete.obs", method = "pearson")

# Tidy long table
cross_cor2 <- as_tibble(R, rownames = "var_hist") |>
  pivot_longer(-var_hist, names_to = "var_fut", values_to = "r")


top_all <- cross_cor2 |> 
  filter(!is.na(r)) |> 
  arrange(desc(abs(r)))   # highest absolute correlations first

print(top_all, n = 20)    # top 20

filter(top_all, var_fut  == "AI_FAO56_alfalfa_fut") |> View()

################################################################################
raw <- pair_periods(
  df   = Data_to_plot$abs,
  hist = "1981_2005",
  fut  = "2076_2100",
  vars = "AI_FAO56_alfalfa",
  by   = c("ensemble","model")
) |>
  transmute(x = AI_FAO56_alfalfa_hist,
            y = AI_FAO56_alfalfa_fut) |>
  filter(is.finite(x), is.finite(y))

cor(raw$x, raw$y, use = "complete.obs", method = "pearson")

################################################################################
make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = AI_FAO56_alfalfa, y = ET / ETo_FAO56_alfalfa, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("AI = PET / P"),  y_lab = bquote("ET / PET"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_AI")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_vs_AI_ggplot2_TIDY.png', plot = p_ET2PET_AI, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = AI_FAO56_alfalfa, y = ET + ETo_FAO56_alfalfa, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("AI = PET / P"),  y_lab = bquote("ET + PET"~"(mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_and_PET_AI")

# Save the plot
ggsave('../plots/ggplot2/ET&PET_vs_AI_ggplot2_TIDY.png', plot = p_ET_and_PET_AI, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = AI_FAO56_alfalfa, y = ET / (ET + ETo_FAO56_alfalfa), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("AI = PET / P"),  y_lab = bquote("ET / (ET + PET)"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2ET_and_PET_AI")

# Save the plot
ggsave('../plots/ggplot2/ET2ET&PET_vs_AI_ggplot2_TIDY.png', plot = p_ET2ET_and_PET_AI, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = VPD, y = ET / ETo_FAO56_alfalfa, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("VPD (kPa)"),  y_lab = bquote("ET / PET"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_VPD")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_vs_VPD_ggplot2_TIDY.png', plot = p_ET2PET_VPD, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = P, y = ET / ETo_FAO56_alfalfa, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("P (mm yr"^"-1"*")"),  y_lab = bquote("ET / PET"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET2PET_P")

# Save the plot
ggsave('../plots/ggplot2/ET2PET_vs_P_ggplot2_TIDY.png', plot = p_ET2PET_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

make_scatter_plot(data = Data_to_plot$abs |>
                    mutate(x = P, y = ET, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
                    select(ensemble, model, color, fill, border, shape, linetype, x, y),
                  FIT = FALSE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote("P (mm yr"^"-1"*")"),  y_lab = bquote("ET (mm yr"^"-1"*")"),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p_ET_vs_P")

# Save the plot
ggsave('../plots/ggplot2/ET2_vs_P_ggplot2_TIDY.png', plot = p_ET_vs_P, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------
# Fit one Budyko curve
make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = ETo_FAO56_alfalfa / P, y = ET / P, ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("AI = PET/P"),
  y_lab = bquote("EI = ET/P"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

BC_general <- Data_to_plot$abs |> 
  filter(!is.na(ETo_FAO56_alfalfa), !is.na(ET), !is.na(P)) |> 
  summarise(
    n = Budyko_curve_optim(ETo_FAO56_alfalfa, ET, P),
    ω = n + 0.72
  )

n_scalar <- pull(BC_general, n)

pb <- ggplot_build(p)
pp <- pb$layout$panel_params[[1]]
xr <- if (!is.null(pp$x.range)) pp$x.range else pp$x$range
yr <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range

# draw the curve over the full x-range
p <- p +
  stat_function(
    fun  = function(x) 1 / (1 + (1 / x)^n_scalar)^(1 / n_scalar),
    xlim = xr,
    inherit.aes = FALSE,
    linewidth = 1,
    color = "grey40"
  ) +
  coord_cartesian(xlim = xr, ylim = yr) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0))


source("./src/delta_derivations.R")

p <- p +
  geom_line(
    data = pred_df |> arrange(AI_hat),
    mapping = aes(x = AI_hat, y = EI_hat_J),  # mapping, not data
    linetype = "dashed",
    linewidth = 1,
    inherit.aes = FALSE
  )

idx_pts_txt <- which(vapply(p$layers,
                            function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                            logical(1)))

# Move them to the end so they draw last (on top of lines etc.)
if (length(idx_pts_txt)) {
  p$layers <- c(p$layers[-idx_pts_txt], p$layers[idx_pts_txt])
}

p <- add_manual_legend(p, x1 = 0.59, x2 = 0.64, y1 = 0.10, y2 = 0.05)

# ggsave('../plots/ggplot2/BC_general_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

# p_BC <- p

# Budyko curve
budyko_obs <- tibble::tibble(
  x = seq(1.8, 2.15, length.out = 80)
) |>
  dplyr::mutate(
    y = 1 / (1 + (1 / x)^n_scalar)^(1 / n_scalar)
  )

# Apparent feedforward response
ff_obs <- pred_df |>
  dplyr::arrange(AI_hat) |>
  dplyr::filter(AI_hat > 1.45, AI_hat < 2.0) |>
  dplyr::transmute(x = AI_hat, y = EI_hat_J)

new_pts <- dplyr::bind_rows(
  budyko_obs,
  budyko_obs |> dplyr::mutate(y = y + 0.01),
  ff_obs,
  ff_obs |> dplyr::mutate(y = y + 0.01)
)

# add obstacle points (visible or invisible)
p <- p +
  geom_point(
    data = new_pts,
    aes(x, y),
    inherit.aes = FALSE,
    alpha = 0  # set to 0 if you want them invisible
  )

repel_idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))
stopifnot(length(repel_idx) == 1)

L <- p$layers[[repel_idx]]

# label column name used by repel layer (robust)
lab_var <- rlang::as_name(L$mapping$label)

# 1) remove any previously-added obstacle rows (prevents growth/duplicates)
lab_data <- L$data %>%
  filter(!is.na(.data[[lab_var]]), .data[[lab_var]] != "")

# 2) grab the per-row colour vector for JUST the label rows
col_vec <- L$aes_params$colour
stopifnot(length(col_vec) == nrow(lab_data))

# 3) rebuild repel data cleanly: labels + NEW obstacles
repel_data <- bind_rows(
  lab_data,
  new_pts %>% mutate(!!lab_var := "")
)

# 4) store colours as a column so length matches new data
repel_data$.repel_col <- c(col_vec, rep(NA, nrow(new_pts)))

# 5) update the existing repel layer in-place (no re-adding!)
L$data <- repel_data
L$aes_params$colour <- NULL
L$mapping$colour <- rlang::sym(".repel_col")

# --- tweak repel behavior ---
# L$geom_params$point.padding       <- 0.04   # 1e-6; distance between labels and points
# L$geom_params$box.padding         <- 0.25    # 0.25; spacing around label boxes
# L$geom_params$force               <- 1.5    # 1; repulsion strength
# L$geom_params$force_pull          <- 1.2    # 1; pull toward anchor point
# L$geom_params$max.overlaps        <- 10000     # 10; overlaps allowed before dropping labels
# L$geom_params$max.iter            <- 20000  # 10000; max optimization iterations
# L$geom_params$max.time            <- 1      # 0.5; max optimization time (seconds)
# 
# # leader line (segment) control
# L$geom_params$min.segment.length  <- 5    # 0.5; suppress short leader lines
# L$geom_params$segment.size        <- 0.5    # 0.5; thinner connector lines
# L$geom_params$segment.alpha       <- 1    # 1; make connector lines subtle
# ---------------------------

p$layers[[repel_idx]] <- L

ggsave('../plots/ggplot2/BC_general_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

library(geomtextpath)

budyko_df <- tibble::tibble(
  AI = seq(from = xr[1], to = xr[2], length.out = 400)
) |>
  dplyr::mutate(
    EI = 1 / (1 + (1 / AI)^n_scalar)^(1 / n_scalar)
  )

p <- p +
  geom_textpath(
    data = budyko_df |> dplyr::filter(AI > 1.8, AI < 2.3),
    aes(x = AI, y = EI, label = "Budyko curve"),
    linewidth = 0,
    inherit.aes = FALSE,
    color = "grey40",
    size = 3.0,
    vjust = 0,              # push label slightly above line
    hjust = 0.5,               # centered along path
    offset = unit(1, "mm")   # extra separation from the line
  )

p <- p +
  geom_textpath(
    data = pred_df |>
      dplyr::arrange(AI_hat) |>
      dplyr::filter(AI_hat > 1.45, AI_hat < 2.1),
    aes(x = AI_hat, y = EI_hat_J, label = "Apparent feedforward response"),
    linewidth = 0,
    inherit.aes = FALSE,
    size = 3.0,
    hjust = 0.5,
    vjust = 0,
    offset = unit(0.7, "mm")   # ABOVE the dashed curve

  )

ggsave('../plots/ggplot2/BC_general_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

p_BC <- p

#-------------------------------------------------------------------------------
make_scatter_plot(
  data = Data_to_plot$abs |>
    mutate(x = VPD,
           y = (ET / P) - (1 / (1 + (1 / (ETo_FAO56_alfalfa / P))^n_scalar)^(1 / n_scalar)), ensemble = interaction(ensemble, PERIOD, drop = TRUE)) |>
    select(ensemble, model, color, fill, border, shape, linetype, x, y),
  FIT = FALSE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote("VPD"),  
  y_lab = bquote("EI – EI"[omega]),
  hline = TRUE, vline = FALSE, one_to_one_line = FALSE,
  save_ggplot2_obj_as="p"
)

ggsave('../plots/ggplot2/EI_resids_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

p_EI_resids <- p
#-------------------------------------------------------------------------------

#--------------
# Combined plot

# Optional but helps: tiny, identical plot margins
tight <- theme(plot.margin = margin(3,3,3,3))
pA <- p_BC + tight
pB <- p_EI_resids + tight
pC <- omega_versus_VPD_and_omega_boxplot + tight
pD <- d_omega_over_omega_vs_d_VPD_over_VPD   + tight


# (2) Align *all five* first; keep all four sides ('tblr')
aligned <- cowplot::align_plots(pA, pB, pC, pD, align = "hv", axis = "tblr")

# (3) Build rows from the aligned grobs
# --- top row with spacing between A and B ---
top <- cowplot::plot_grid(
  aligned[[1]], NULL, aligned[[2]],
  ncol = 3,
  rel_widths = c(1, 0.05, 1),
  labels = c("a", "", "b"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = 0.1
)

# --- bottom row with spacing between C, D ---
bottom <- cowplot::plot_grid(
  aligned[[3]], NULL, aligned[[4]],
  ncol = 3,
  rel_widths = c(1, 0.05, 1),
  labels = c("c", "", "d"),
  label_colour = "#333333",
  label_size = 12,
  hjust = -1, vjust = -0.1
)

# --- add vertical spacing between top and bottom ---
combined <- cowplot::plot_grid(
  NULL, top, NULL, bottom,
  ncol = 1,
  rel_heights = c(0.1, 3, 0.1, 3)  # the middle '0.1' is vertical gap
)

# (5) Save (use any size you like; e.g. 240×200 mm fits a 3:2 row split nicely)
ggsave("../plots/ggplot2/combined_BC,EI,omega_ggplot2_TIDY.png", combined,
       width = 240, height = 239, units = "mm", dpi = RES, bg = "white")

#-------------------------------------------------------------------------------

# What next
# - can be BC with unchanged omega considered water-limited only and what is above is the effect of VPD
# - what if one single omega is fitted and the difference is in ET or EI are analysed in repsonse to VPD?
# - determine what should be a reduction of rs in PM that would ensure that all models are following the same BC
# - check in BC why the in the energy-limited region the gray line does not origin in zero-zero

