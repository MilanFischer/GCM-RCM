# Clear the workspace
rm(list=ls()); gc()

library(DEoptim)
# library(ggplot2)
library(grid)
library(patchwork)
# library(gridExtra)
library(ggrepel)
# library(scales)
# library(png)
# library(magick)
# library(gtable)
library(MASS)
library(ggrepel)
library(tidyverse)

# Preprocess
source("./src/general_functions.R")
source('./src/make_plot.R')
source("./src/colors.R")
source("./src/Budyko_curve_tidy.R")

# VPD from daily values and then averaged is about 15% higher, so constant 1.15 can be apllied

# Plotting setup
# COL_CMIP5="#117733"
# COL_CMIP6="#AA4499"
# COL_RCMs="#6699CC"

COL_CMIP5 <- "#1c91df"
COL_CMIP6 <- "#4e79a7"
COL_RCMs <- "#ff4500"
BORDER_COLOR <- "#2b2b2b"

c("#1f77b4", "#2a9df4", "#4e79a7", "#5f9ea0", "#1c91df", "#4682b4",
  "#0073cf", "#87cefa", "#2e8b57", "#0a74da", "#d62728", "#ff6347",
  "#e74c3c", "#c94c4c", "#ff4500", "#3cb371", "#2e8b57", "#66cdaa")

Pl_width = 120
Pl_height = 120
oma_mm = c(20, 20, 0, 0)
mar_mm = c(0, 0, 0, 0) + 2
n_row = 1
n_col = 1

RES = 600

# Specify the list of driving GCM models
driving_GCMs <- c("MPI", "EARTH", "HADGEM", "CNRM", "IPSL")

# Load datasets from CSV files and add an 'ENSEMBLE' column to each
CMIP5_data <- read_csv("../Inputs/CMIP5.csv") |> 
  mutate(ENSEMBLE = "CMIP5")

CMIP6_data <- read_csv("../Inputs/CMIP6.csv") |> 
  mutate(ENSEMBLE = "CMIP6")

RCMs_data <- read_csv("../Inputs/EUR-44.csv") |> 
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

# Based on the "SUPPLEMENT_DataTables_Meinshausen_6May2020.xlsx"
CO2_2076_2100_RCP85 <- 979.42 # ppm
CO2_1981_2005 <- 359.47
# The correction by Yang et al. (2018) https://doi.org/10.1038/s41558-018-0361-0
d_rs <- 0.05*(CO2_2076_2100_RCP85 - CO2_1981_2005)
d_g_eff <- 1/d_rs * 1000

# Surface resistance of the grass or alfalfa reference surface (s/m)
rs_ref_grass <- 70
rs_ref_alfalfa <- 40 # Original value in Allen et al. (2005) is 45 s/m
rs_ref_grass_2076_2100_RCP85 <- rs_ref_grass + 0.05*(CO2_2076_2100_RCP85 - CO2_1981_2005)
rs_ref_alfalfa_2076_2100_RCP85 <- rs_ref_alfalfa + 0.05*(CO2_2076_2100_RCP85 - CO2_1981_2005)

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

# Function for gapfilling of u and RH
gap_fill_tidy <- function(data, var, period) {
  data |>
    group_by(ENSEMBLE, MONTH) |>
    mutate(
      VALUE = if_else(
        VAR == var & PERIOD == period & is.na(VALUE),
        mean(VALUE[VAR == var & PERIOD == period & MODEL != "ERA5"], na.rm = TRUE),
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
    VPD_adj = VPD * 1.15, # Linearly adjusted to account for the impact of temporal averaging
    # Slope of the saturated vapor pressure curve
    SVP = b * c * e_sat/(c + Ta)^2,
    # Latent heat of vaporization (MJ/kg)
    lambda = (2.501 * (10^6) - 2361 * Ta) / (10^6),
    # Unit conversions
    # W_to_mm = 1 / lambda * 3600*24 / 10^6,
    W_to_mm = 0.408 * 3600*24 / 10^6,
    # Evapotranspiration
    ET = W_to_mm * DAYS_IN_MONTH * LE,
    # Runoff
    RO = P - ET,
    # Equilibrium evaporation
    ET_eq = SVP*(H + LE) / (SVP + gamma) * W_to_mm * DAYS_IN_MONTH,
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
      (SVP * A + rhoAir * CpAir * VPD_adj / ra) /
        (SVP + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 grass with correction for CO2 effect on stomatal conductance for GCMs
    # Only when air temperature is above 5 deg. C
    ETo_FAO56_grass_GCM_CO2_corr = {
      zo = 0.123 * 0.12
      u2 = log(2 / zo) / log(10 / zo) * u10
      ra = 208 / u2
      rs = if_else(
        PERIOD == "2076_2100" & ENSEMBLE %in% c("CMIP5", "CMIP6") & Ta > 5,
        rs_ref_grass_2076_2100_RCP85,
        rs_ref_grass
      )
      (SVP * A + rhoAir * CpAir * VPD_adj / ra) /
        (SVP + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 alfalfa (with inline temporary vars)
    ETo_FAO56_alfalfa = {
      zo = (0.123*0.12)
      u2 = log(2 / zo) / log(10 / zo) * u10
      # u2 = u10 * 4.87 / log(67.8 * 10 - 5.42)
      ra = 118 / u2
      rs = rs_ref_alfalfa
      (SVP * A + rhoAir * CpAir * VPD_adj / ra) /
        (SVP + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    },
    # ETo FAO56 alfalfa with correction for CO2 effect on stomatal conductance for GCMs
    # Only when air temperature is above 5 deg. C
    ETo_FAO56_alfalfa_GCM_CO2_corr = {
      zo = (0.123*0.12)
      u2 = log(2 / zo) / log(10 / zo) * u10
      # u2 = u10 * 4.87 / log(67.8 * 10 - 5.42)
      ra = 118 / u2
      rs = if_else(
        PERIOD == "2076_2100" & ENSEMBLE %in% c("CMIP5", "CMIP6") & Ta > 5,
        rs_ref_alfalfa_2076_2100_RCP85,
        rs_ref_alfalfa
      )
      (SVP * A + rhoAir * CpAir * VPD_adj / ra) /
        (SVP + gamma * (1 + rs / ra)) * W_to_mm * DAYS_IN_MONTH
    }
  )

# Test <- Data_wide |> 
#   group_by(MODEL, ENSEMBLE) |>
#   select(MONTH, MONTH_NAME, PERIOD, ENSEMBLE, MODEL, ET, ET_eq) |> 
#   filter(MODEL %in% c("CNRM", "EARTH"), ENSEMBLE == "CMIP5", MONTH %in% 4:8) |> 
#   mutate(max_ET_ratio = max(ET / ET_eq)) |> 
#   ungroup()

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

write_csv(Data_wide, "test.csv")
# View(Data_wide |> 
#        filter(PERIOD == "1981_2005", ENSEMBLE == "CMIP5", MODEL == "EARTH") |> 
#        select(MONTH, MONTH_NAME, PERIOD, ENSEMBLE, MODEL, ET))

# 
# Data_wide <- Data_wide |> 
#   select(MONTH, MONTH_NAME, PERIOD, ENSEMBLE, MODEL, max_ET_ratio)


# View(filter(Data_wide, MODEL %in% c("CNRM"), ENSEMBLE == "CMIP5"))


Data_wide <- Data_wide |> 
  filter(MODEL != "ERA5")

# Data_wide <- Data_wide |> 
#   filter(PERIOD != "2076_2100")

# Create the boxplot
p <- ggplot(Data_wide, aes(x = ENSEMBLE, y = max_ET_ratio, fill = ENSEMBLE)) +
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2, coef = 3) +    # Bars width
  geom_boxplot(coef = 3) +
  scale_fill_manual(values = c(COL_CMIP5, COL_CMIP6, COL_RCMs)) +
  labs(y = bquote('Priestley-Taylor'~alpha), x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

# Save the plot to a file
ggsave(filename = "../Plots/ggplot2/Priestley-Taylor_alpha_ggplot2_TIDY.png", plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")

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

# View(annual_stats |> 
#        filter(PERIOD == "1981_2005", ENSEMBLE == "CMIP5",
#               MODEL == "EARTH", VAR == "ET"))



# annual_stats <- Data_long |> 
#   group_by(VAR, PERIOD, ENSEMBLE, MODEL) |> 
#   summarise(
#     annual_stat = if (unique(VAR) == "P") {
#       sum(VALUE, na.rm = FALSE)  # Sum for precipitation
#     } else {
#       sum(VALUE * DAYS_IN_MONTH, na.rm = FALSE) / sum(DAYS_IN_MONTH, na.rm = FALSE) # Weighted mean for other variables
#     },
#     .groups = "drop"
#   )

# Verify unique variables in annual statistics
annual_stats |> 
  pull(VAR) |> 
  unique()

# Inspect annual statistics specifically for precipitation
# annual_stats |> 
#   filter(VAR == "P")

# Compute mean values across models for each variable, period, and ensemble
mean_values <- annual_stats |> 
  group_by(VAR, PERIOD, ENSEMBLE) |> 
  summarize(
    mean_value = mean(annual_stat, na.rm = TRUE),
    .groups = "drop"
  )

# Visualize mean values using a bar plot
mean_values |> 
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
  theme_minimal()


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
    Bo = H/LE,
    # Integrated water vapor (mm)
    # Ruckstuhl et al. (2007) Observed relationship between surface specific humidity,
    # integrated water vapor, and longwave downward radiation at different altitude
    IWV = 2.67 * q + 0.68,
    # Precipitation efficiency (%)
    PE = P / (IWV * 365.25) * 100,
    # PE = P/e
  )

# Convert derived variables back to long format and append to the main dataset
annual_stats <- bind_rows(
  annual_stats,
  annual_stats_wide |>
    select(PERIOD, ENSEMBLE, MODEL, rs_eff, g_eff, Bo, IWV, PE) |>
    pivot_longer(cols = c(rs_eff, g_eff, Bo, IWV, PE), names_to = "VAR", values_to = "annual_stat")
)

# Fit a linear model
fit <- lm(ETo_FAO56_alfalfa ~ ETo_FAO56_grass, data = annual_stats_wide)

# Extract coefficients and R²
coefs <- coef(fit)
r_squared <- summary(fit)$r.squared

# Create equation label
eq_label <- bquote(
  italic(y) == .(format(coefs[2], digits = 3)) %.% italic(x) + .(format(coefs[1], digits = 2)) ~~
    R^2 == .(format(r_squared, digits = 3))
)

ggplot(annual_stats_wide, aes(x = ETo_FAO56_grass, y = ETo_FAO56_alfalfa)) +
  geom_point(color = "#1c91df", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = as.expression(eq_label),
           hjust = 1.1, vjust = -0.5, parse = TRUE, size = 4, color = "black") +
  labs(
    x = expression("ETo"["FAO56 grass"]~"(mm)"),
    y = expression("ETo"["FAO56 alfalfa"]~"(mm)"),
    title = "FAO56 vs VPD-adjusted ETo"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

annual_stats_wide |> mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> pull(ET_diff) |> plot()

annual_stats_wide |> mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> pull(ET_diff) |> min()

annual_stats_wide |>  mutate(ET_diff = ETo_FAO56_alfalfa - ET) |> filter(ET_diff < 0)


# Calculate the normalized difference between periods
# norm_diff <- annual_stats |> 
#   group_by(VAR, ENSEMBLE, MODEL) |> 
#   summarize(
#     d_annual_stat = (annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"]) / annual_stat[PERIOD == "1981_2005"],
#     .groups = 'drop'
#   ) |> 
#   mutate(VAR = paste0("d_", VAR, "_over_", VAR))

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
# diff <- annual_stats |> 
#   group_by(VAR, ENSEMBLE, MODEL) |> 
#   summarize(
#     d_annual_stat = annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"],
#     .groups = 'drop'
#   ) |> 
#   mutate(VAR = paste0("d_", VAR))
diff <- annual_stats |> 
  group_by(VAR, ENSEMBLE, MODEL) |> 
  reframe(
    d_annual_stat = annual_stat[PERIOD == "2076_2100"] - annual_stat[PERIOD == "1981_2005"],
    VAR = paste0("d_", VAR[1])
  )

# Combine the two data frames
all_diff <- bind_rows(norm_diff, diff)


# Reshape data: pivot wider to have P and ET in separate columns
Data_to_plot <- all_diff |> 
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
  rename(model = ENSEMBLE, label = MODEL)

#----------------
# dP/P vs. dET/ET
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.88, 0.83, 0.78),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

LM_eq_labels <- tibble(
  x = c(0.7, 0.7, 0.7),
  y = c(0.14, 0.09, 0.04)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_ET_over_ET, d_P_over_P, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_ET_over_ET", y = "d_P_over_P"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'ET/'*'ET'),  y_lab = bquote(Delta*'P/'*'P'),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.87, 0.95), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels, plot_labels = Plot_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_P_over_P_vs_delta_ET_over_ET_ggplot2_TIDY.png", save_ggplot2_obj_as="p1")

#----------------
# dRO/RO vs. dET/ET
LM_eq_labels <- data.frame(
  x = c(0.27, 0.27, 0.27),
  y = c(0.22, 0.17, 0.12)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_ET_over_ET, d_RO_over_RO, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_ET_over_ET", y = "d_RO_over_RO"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'ET/'*'ET'),  y_lab = bquote(Delta*'R'[o]*'/'*'R'[o]),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.80, 0.96), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_RO_over_RO_vs_delta_ET_over_ET_ggplot2_TIDY.png", save_ggplot2_obj_as="p2")

#----------------
# dET/ET vs. dP/P
LM_eq_labels <- tibble(
  x = c(0.7, 0.7, 0.7),
  y = c(0.14, 0.09, 0.04)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_P_over_P, d_ET_over_ET, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_P_over_P", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.79), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_ET_over_ET_vs_delta_P_over_P_ggplot2_TIDY.png", save_ggplot2_obj_as="p3")

#----------------
# dET/ET vs. dRO/RO
make_scatter_plot(data = Data_to_plot |>
                    select(d_P_over_P, d_RO_over_RO, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_P_over_P", y = "d_RO_over_RO"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'P/'*'P'),  y_lab = bquote(Delta*'R'[o]*'/'*'R'[o]),
                  one_to_one_line = TRUE, one_to_one_line_lab = c(0.75, 0.97), robust_regression = TRUE,
                  LM_eq_labels = LM_eq_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_RO_over_RO_vs_delta_P_over_P_ggplot2_TIDY.png", save_ggplot2_obj_as="p4")


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
labels <- c("a)", "b)", "c)", "d)")

plots <- Map(function(plot, label) {
  plot + annotation_custom(
    textGrob(label, x = unit(0.05, "npc"), y = unit(0.96, "npc"),
             gp = gpar(fontsize = 12))
  )
}, plots, labels)

p1 <- plots$p1; p2 <- plots$p2; p3 <- plots$p3; p4 <- plots$p4

#---- Combine and save
panel_figure <- (p1 | p2) / (p3 | p4)

ggsave('../Plots/ggplot2/panel_fig_hydrological_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#-------------------------
# Delta_T versus delta_P/P
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

# Manual adjustment
Y_range_man <- c(-0.12, 0.32)

make_scatter_plot(data = Data_to_plot |>
                    select(d_Ta, d_P_over_P, model, color, fill, border, shape, label, linetype) |> 
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
ggsave('../Plots/ggplot2/delta_P_over_P_vs_Ta_ggplot2_TIDY.png', plot = p1, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Delta_T versus delta_ET/ET

make_scatter_plot(data = Data_to_plot |>
                    select(d_Ta, d_ET_over_ET, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Ta", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, Y_range_man = Y_range_man,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  save_ggplot2_obj_as="p2")

p2 <- p2 + geom_line(data = DF1_CC_nonlinear, aes(x = X, y = Y), linetype = 2) +
  # geom_line(data = DF1_CC_linear, aes(x = X, y = Y), linetype = 2, color = "darkred") +
  annotate("text", x = rescale(0.02, X_range), y = rescale(0.77, Y_range), label = "Clausius–Clapeyron", hjust = 0, color = "#2b2b2b", angle=angle_degrees, size=4)

# Save the plot
ggsave('../Plots/ggplot2/delta_ET_over_ET_vs_Ta_ggplot2_TIDY.png', plot = p2, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---------------------------
# Delta_T versus d_VPD/VPD

make_scatter_plot(data = Data_to_plot |>
                    select(d_Ta, d_VPD_over_VPD, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Ta", y = "d_VPD_over_VPD"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'VPD/'*'VPD'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
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
ggsave('../Plots/ggplot2/delta_VPD_over_VPD_ggplot2_TIDY.png', plot = p3, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

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
    g_eff_corr = (1 / (rs_eff_2076_2100 + d_rs)) * 1000,
    d_g_eff_CO2_corr = (g_eff_corr - g_eff_1981_2005) / g_eff_1981_2005
  ) |>
  select(MODEL, d_g_eff_CO2_corr)

# --- Prepare plotting data ---

# Original CMIP5, CMIP6, and EUR-44 points
points_base <- Data_to_plot |>
  select(d_Ta, d_g_eff_over_g_eff, model, color, fill, border, shape, label, linetype) |>
  rename(x = d_Ta, y = d_g_eff_over_g_eff)

# CO2-corrected EUR-44 points
points_corr <- Data_to_plot |>
  filter(model == "EUR-44") |>
  left_join(g_eff_corr_RCMs, by = c("label" = "MODEL")) |>
  mutate(
    model = "EUR-44_CO2_corr",
    x = d_Ta,
    y = d_g_eff_CO2_corr,
    fill = "#2e8b57",
    color = "#2e8b57",
    linetype = "dashed"
  ) |>
  select(x, y, model, color, fill, border, shape, label, linetype)

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

p4 <- p4 + annotate("text", x = rescale(0.66, X_range), y = rescale(0.96, Y_range), label = bquote(bold("EUR-44 CO"["2"]~"corr")), hjust = 0, color = "#2e8b57", size=4)

# Save the plot
ggsave('../Plots/ggplot2/delta_g_eff_over_g_eff_vs_Ta_ggplot2_TIDY.png', plot = p4, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#---- Subplot labels
plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
labels <- c("a)", "b)", "c)", "d)")

plots <- Map(function(plot, label) {
  plot + annotation_custom(
    textGrob(label, x = unit(0.05, "npc"), y = unit(0.96, "npc"),
             gp = gpar(fontsize = 12))
  )
}, plots, labels)

p1 <- plots$p1; p2 <- plots$p2; p3 <- plots$p3; p4 <- plots$p4

#---- Combine and save
panel_figure <- (p1 | p2) / (p3 | p4)

ggsave('../Plots/ggplot2/panel_fig_Ta_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

# notes: make the labels to be read from csv file, avoid overlapping with lines, allow to not show labels for selected data

#-------------------------------------------------------------------------------

#---------------------------------------
# Effective conductance for double check
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

Data_to_plot <- Data_to_plot |> 
  mutate(d_g_eff_over_g_eff_CHECK = (d_Rn_over_Rn - d_Bo_adj - d_VPD_over_VPD) / (1 + d_VPD_over_VPD)
  )

make_scatter_plot(data = Data_to_plot |>
                    select(d_g_eff_over_g_eff_CHECK, d_g_eff_over_g_eff, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_g_eff_over_g_eff_CHECK", y = "d_g_eff_over_g_eff"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote('('*Delta*'R'['n']*'/'*'R'['n']~'-'~Delta*'B'['o']*'/(1 + '*'B'['o']*')'~'-'~Delta*'VPD'*'/VPD'*') / ('*1 + ~Delta*'VPD'*'/VPD)'),
                  y_lab = bquote(Delta*'g'['eff']~'/'~'g'['eff']),
                  hline = FALSE, vline = FALSE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.93, 0.84), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_g_eff_star_over_g_eff_star_check_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

# Reshape data: pivot wider to have P and ET in separate columns
Data_to_plot_II <- annual_stats |> 
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
  rename(model = ENSEMBLE, label = MODEL)

#---------------------------
# VPD versus g_eff

Plot_labels <- tibble(
  x = c(0.785, 0.785, 0.785),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(x = VPD, y = g_eff, model = interaction(model, PERIOD, drop = TRUE)) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
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
    annotate("text", x = rescale(x2 , X_range), y = rescale(y1, Y_range), label = "History", size = 3, hjust = 0)
  
  p <- p +
    annotate("point", x = rescale(x1, X_range), y = rescale(y2, Y_range), shape = 24, size = 2, stroke = 0.8,
             color = "#2b2b2b", fill = NA) +
    annotate("text", x = rescale(x2 , X_range), y = rescale(y2, Y_range), label = "Scenario", size = 3, hjust = 0)
  
  return(p)
}

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for GCMs

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        1 / (rs_eff - d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff,
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#-------------------------------
# VPD versus g_eff_corr for RCMs

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        1 / (rs_eff + d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff,
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------------------------------------
# VPD versus g_eff normalized by available energy

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(x = VPD, y = g_eff / (H + LE), model = interaction(model, PERIOD, drop = TRUE)) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_AE_norm_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#----------------------------------------------------------------------
# VPD versus g_eff GCM CO2 corrected and normalized by available energy

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        1 / (rs_eff - d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (H + LE),
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_AE_norm_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff RCM CO2 corrected and normalized by available energy

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        1 / (rs_eff + d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (H + LE),
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_AE_norm_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------------------------------------
# VPD versus g_eff normalized by global radiation

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(x = VPD, y = g_eff / (Rg), model = interaction(model, PERIOD, drop = TRUE)) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_Rg_norm_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff GCM CO2 corrected and normalized by global radiation

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        1 / (rs_eff - d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (Rg),
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_Rg_norm_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------------------------------------------------
# VPD versus g_eff RCM CO2 corrected and normalized by global radiation

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      g_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        1 / (rs_eff + d_rs) * 1000,
        g_eff
      ),
      x = VPD,
      y = g_eff / (Rg),
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~"(kPa"^-1*")"),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.80, x2 = 0.84, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/g_eff_Rg_norm_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#------------------
# VPD versus rs_eff

Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)


make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(x = VPD, y = rs_eff, model = interaction(model, PERIOD, drop = TRUE)) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('r'['s eff']~(s/m)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/rs_eff_versus_VPD_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------
# VPD versus g_eff_corr for GCMs

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      rs_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("CMIP5", "CMIP6"),
        rs_eff - d_rs,
        rs_eff
      ),
      x = VPD,
      y = rs_eff,
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/rs_eff_versus_VPD_GCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#-------------------------------
# VPD versus g_eff_corr for RCMs

make_scatter_plot(
  data = Data_to_plot_II |>
    mutate(
      rs_eff = if_else(
        PERIOD == "2076_2100" & model %in% c("EUR-44"),
        rs_eff + d_rs,
        rs_eff
      ),
      x = VPD,
      y = rs_eff,
      model = interaction(model, PERIOD, drop = TRUE)
    ) |>
    select(model, label, color, fill, border, shape, linetype, x, y),
  FIT = TRUE, robust_regression = TRUE,
  xy_round = 0.05, xy_offset = 0.04,
  x_lab = bquote('VPD (kPa)'),  
  y_lab = bquote('g'['eff']~(mm/s)),
  hline = FALSE, vline = FALSE, one_to_one_line = FALSE,
  plot_labels = Plot_labels,
  save_ggplot2_obj_as="p"
)

p <- add_manual_legend(p, x1 = 0.04, x2 = 0.08, y1 = 0.75, y2 = 0.70)

ggsave('../Plots/ggplot2/rs_eff_versus_VPD_RCM_corr_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

# Assessing whether aerodynamics or thermodynamics dominate the relative change in ET

make_scatter_plot(data = Data_to_plot |>
                    mutate(OMEGA = d_g_eff_over_g_eff / d_VPD_over_VPD) |>
                    select(d_Ta, OMEGA, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = d_Ta, y = OMEGA),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'g'['eff']~'/'~'g'['eff']~'/'~'('*Delta*'VPD'~'/'~'VPD)'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_path = "../Plots/ggplot2", plot_name = "Omega_versus_delta_Ta_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

#-------------------------------------------------
# Delta_Rn and Bowen ratio term versus delta_ET/ET
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot |>
                    mutate(d_A_minus_Bo_adj = d_A_over_A - d_Bo_adj) |>
                    select(d_A_minus_Bo_adj, d_ET_over_ET, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = d_A_minus_Bo_adj, y = d_ET_over_ET),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G) - Delta * B[o] / (1 + B[o])),
                  y_lab = bquote(Delta*'ET/'*'ET'),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.89), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_Rn&delta_Bo_versus_delta_ET_over_ET_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

#---------------------------
# Delta_T versus d_e/e
make_scatter_plot(data = Data_to_plot |>
                    select(d_Ta, d_e_over_e, model, color, fill, border, shape, label, linetype) |> 
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

ggsave('../Plots/ggplot2/delta_e_over_e_ggplot2_TIDY.png', plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#-----------------------
# Delta_T versus d_PE/PE
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_Ta, d_PE_over_PE, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Ta", y = "d_PE_over_PE"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'T'['a']~'(°C)'),  y_lab = bquote(Delta*'PE'~'/'~'PE'),
                  hline = FALSE, vline = FALSE, one_to_one_line = FALSE, robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_PE_over_PE_vs_Ta_ggplot2_TIDY.png", save_ggplot2_obj_as="p")

annual_stats_wide |>
  reframe(PE_mean = mean(PE, na.rm = TRUE), .by = c(ENSEMBLE, PERIOD))

#-------------------------------------------------------------------------------

#-----------------------------
# Delta Rg versus delta SW_net
Plot_labels <- tibble(
  x = c(0.03, 0.03, 0.03),
  y = c(0.88, 0.83, 0.78),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_Rg_over_Rg, d_SW_net_over_SW_net, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_SW_net_over_SW_net"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'SW'['net']*'/'*'SW'['net']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.94, 0.92), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  save_ggplot2_obj_as="p5")

# Save the plot
ggsave('../Plots/ggplot2/delta_Rg_versus_delta_SW_net_ggplot2_TIDY.png', plot = p5, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')


#-----------------------------
# Delta Rg versus delta LW_net

# Note that the mean LW_net is a negative number
# CMIP5 has slightly more negative values in the future
# RCMs has less negative values in the future
# Because the denominator is negative, the d_LW_net_over_LW_net_RCMs_a with delta being positive becomes a negative number

make_scatter_plot(data = Data_to_plot |>
                    select(d_Rg_over_Rg, d_LW_net_over_LW_net, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_LW_net_over_LW_net"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'LW'['net']*'/'*'LW'['net']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.02, 0.25), robust_regression = TRUE,
                  save_ggplot2_obj_as="p6")

# Save the plot
ggsave('../Plots/ggplot2/delta_Rg_versus_delta_LW_net_ggplot2_TIDY.png', plot = p6, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#-------------------------
# Delta Rg versus delta Rn
Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_Rg_over_Rg, d_Rn_over_Rn, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_Rg_over_Rg", y = "d_Rn_over_Rn"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta*'R'['g']*'/'*'R'['g']),  y_lab = bquote(Delta*'R'['n']*'/'*'R'['n']),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.94, 0.75), robust_regression = TRUE,
                  save_ggplot2_obj_as="p7")

# Save the plot
ggsave('../Plots/ggplot2/delta_Rg_versus_delta_Rn_ggplot2_TIDY.png', plot = p7, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#----------------------------
# Delta_Rn versus delta_ET/ET

make_scatter_plot(data = Data_to_plot |>
                    select(d_A_over_A, d_ET_over_ET, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_A_over_A", y = "d_ET_over_ET"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G)),  y_lab = bquote(Delta * ET / ET %~~% Delta * (R[n] - G) / (R[n] - G) - Delta * B[o] / (1 + B[o])),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.92, 0.91), robust_regression = TRUE,
                  save_ggplot2_obj_as="p8")

# Save the plot
ggsave('../Plots/ggplot2/delta_Rn_versus_delta_ET_over_ET_ggplot2_TIDY.png', plot = p8, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')

#---- Subplot labels
plots <- list(p5 = p5, p6 = p6, p7 = p7, p8 = p8)
labels <- c("a)", "b)", "c)", "d)")

plots <- Map(function(plot, label) {
  x_pos <- if (label == "d)") 0.08 else 0.05  # Move "d)" to the right
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

ggsave('../Plots/ggplot2/panel_fig_Rg&Rn_relations_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height*2, dpi = RES, units = 'mm')

#-------------------------------------------------------------------------------

#------------------------------------
# Delta_Rn and delta Bowen ratio term

Plot_labels <- tibble(
  x = c(0.84, 0.84, 0.84),
  y = c(0.96, 0.91, 0.86),
  model = c("CMIP5", "CMIP6", "EUR-44"),
  color = c(COL_CMIP5, COL_CMIP6, COL_RCMs)
)

make_scatter_plot(data = Data_to_plot |>
                    select(d_A_over_A, d_Bo_adj, model, color, fill, border, shape, label, linetype) |> 
                    rename(x = "d_A_over_A", y = "d_Bo_adj"),
                  FIT = TRUE, xy_round = 0.05, xy_offset = 0.04,
                  x_lab = bquote(Delta * (R[n] - G) / (R[n] - G)),
                  y_lab = bquote(Delta * B[o] / (1 + B[o])),
                  hline = TRUE, vline = TRUE, one_to_one_line = TRUE, one_to_one_line_lab = c(0.59, 0.96), robust_regression = TRUE,
                  plot_labels = Plot_labels,
                  plot_path = "../Plots/ggplot2", plot_name = "delta_Rn_versus_delta_Bo_ggplot2_TIDY.png")

#-------------------------------------------------------------------------------

# Focus on complementary relation in Budyko space

# Budyko curves

#-------------------------------------------------------------------------------
# Fitting Budyko curve

# Compute Budyko parameter n per MODEL and PERIOD
n_estimates <- annual_stats_wide |> 
  group_by(ENSEMBLE, MODEL, PERIOD) |> 
  # Safely wrap the function to prevent it from failing on bad data
  reframe(n = purrr::possibly(Budyko_curve_optim, otherwise = NA_real_)(PET = ETo_FAO56_alfalfa, ET = ET, P = P))


# Ensure PERIOD is ordered if needed
n_estimates$PERIOD <- factor(n_estimates$PERIOD, levels = c("1981_2005", "2076_2100"), labels = c("Historical", "Future"))

# Boxplot
ggplot(n_estimates, aes(x = ENSEMBLE, y = n, fill = PERIOD)) +
  geom_boxplot(position = position_dodge(0.7), width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = PERIOD), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.6, size = 1.5) +
  labs(
    title = "Budyko parameter n by Ensemble and Period",
    x = "Ensemble",
    y = "Budyko parameter (n)",
    fill = "Period",
    color = "Period"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


# Run paired t-test and extract relevant summary
test_n_diff <- function(ens) {
  df <- n_estimates |> 
    filter(ENSEMBLE == ens) |> 
    select(MODEL, PERIOD, n) |> 
    pivot_wider(names_from = PERIOD, values_from = n) |>
    drop_na()
  
  ttest <- t.test(df[[2]], df[[3]], paired = TRUE)
  
  tibble(
    ENSEMBLE = ens,
    n_models = nrow(df),
    t_stat = ttest$statistic,
    p_value = ttest$p.value,
    mean_diff = mean(df[[2]] - df[[3]], na.rm = TRUE),
    CI_lower = ttest$conf.int[1],
    CI_upper = ttest$conf.int[2]
  )
}

# Run the test for all ensembles
ensembles <- unique(n_estimates$ENSEMBLE)

test_results <- map_dfr(ensembles, test_n_diff)

print(test_results)

# View or join this back to your main dataset if needed
View(n_estimates)

n_estimates |> pull(n) |> max()
n_estimates |> pull(n) |> hist()


out_BC_FAO56_alfalfa <- Budyko_curve(annual_stats_wide,
             pet_col = "ETo_FAO56_alfalfa", et_col = "ET", p_col = "P",
             X_range = c(0.0, 2.2), Y_range = c(0.0, 1.06),
             Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
             plot = TRUE, plot_name = "Budyko_curve_ggplot2.png")

out_BC_FAO56_alfalfa_CO2_corr <- Budyko_curve(annual_stats_wide,
             pet_col = "ETo_FAO56_alfalfa_GCM_CO2_corr", et_col = "ET", p_col = "P",
             X_range = c(0.0, 2.2), Y_range = c(0.0, 1.06),
             Xin_range = c(0.8, 1.6), Yin_range = c(0.60, 0.8),
             plot = TRUE, plot_name = "Budyko_curve_CO2_ggplot2.png")

# Add label "a)" to the first plot
p9 <- out_BC_FAO56_alfalfa$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "a)", hjust = 1.5, vjust = 1.5, size = 6)

# Add label "b)" to the second plot
p10 <- out_BC_FAO56_alfalfa_CO2_corr$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "b)", hjust = 1.5, vjust = 1.5, size = 6)

# Combine the two plots side by side using cowplot
panel_figure <- (p9 | p10) 

# Save the combined plot
ggsave('../Plots/ggplot2/BC_PET_PM_and_PM_corr_TIDY.png', plot = panel_figure, width = Pl_width*2, height = Pl_height, dpi = RES, units = 'mm')






# Recode period to match annual_stats_wide
n_estimates_fixed <- n_estimates |>
  mutate(PERIOD = recode(PERIOD,
                         "Historical" = "1981_2005",
                         "Future" = "2076_2100"))

# Now join with annual_stats_wide to get LAI
n_lai_data <- n_estimates_fixed |>
  left_join(
    annual_stats_wide |> select(PERIOD, ENSEMBLE, MODEL, LAI),
    by = c("PERIOD", "ENSEMBLE", "MODEL")
  ) |>
  drop_na(n, LAI)

# Scatter plot of n vs. LAI
ggplot(n_lai_data, aes(x = LAI, y = n, color = ENSEMBLE, shape = PERIOD)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(
    title = "Budyko Parameter (n) vs. Leaf Area Index (LAI)",
    x = "Leaf Area Index (LAI)",
    y = "Budyko Parameter (n)",
    color = "Ensemble",
    shape = "Period"
  ) +
  theme_minimal()





# Join LAI to Budyko `n`
n_lai_data <- n_estimates |>
  mutate(PERIOD = recode(PERIOD,
                         "Historical" = "1981_2005",
                         "Future" = "2076_2100")) |>
  left_join(
    annual_stats_wide |> select(PERIOD, ENSEMBLE, MODEL, LAI),
    by = c("PERIOD", "ENSEMBLE", "MODEL")
  )

# Pivot wider to get both periods side-by-side
n_lai_wide <- n_lai_data |>
  pivot_wider(names_from = PERIOD, values_from = c(n, LAI), names_sep = "_")

# Compute normalized differences
norm_diff_n_lai <- n_lai_wide |>
  mutate(
    d_n   = (n_2076_2100 - n_1981_2005) / n_1981_2005,
    d_LAI = (LAI_2076_2100 - LAI_1981_2005) / LAI_1981_2005
  )


library(broom)  # <- required for tidy()

# Safe linear model wrapper
safe_lm <- possibly(function(df) lm(d_n ~ d_LAI, data = df), otherwise = NULL)

# Fit model for each ensemble
ensemble_fits <- norm_diff_n_lai  |> 
  group_by(ENSEMBLE) |>
  nest() |>
  mutate(
    model = map(data, safe_lm),
    summary = map(model, function(m) {
      if (!is.null(m)) broom::tidy(m) else NULL
    })
  ) |>
  unnest(summary) |>
  filter(term == "d_LAI") |>
  select(ENSEMBLE, estimate, std.error, statistic, p.value)

fit_lines <- norm_diff_n_lai %>%
  group_by(ENSEMBLE) %>%
  nest() %>%
  mutate(
    model = map(data, safe_lm),
    fitted = map2(model, data, ~ if (!is.null(.x)) {
      predict(.x, newdata = .y)
    } else {
      rep(NA, nrow(.y))
    })
  ) %>%
  unnest(c(data, fitted)) %>%
  rename(fitted_d_n = fitted)

# Plot with actual points and fitted lines
ggplot(fit_lines, aes(x = d_LAI, y = d_n, color = ENSEMBLE)) +
  geom_point() +
  geom_line(aes(y = fitted_d_n)) +
  theme_minimal() +
  labs(
    title = "Normalized ΔLAI vs Δn with Fitted Lines",
    x = "Normalized Change in LAI (ΔLAI)",
    y = "Normalized Change in n (Δn)"
  )


























Data_wide



# Here I ended 2025-06-10

Data_to_plot |>
  select(d_Ta, d_g_eff_over_g_eff, model, color, fill, border, shape, label, linetype) |> 
  filter(model == "CMIP5")



Data_wide |> View()





################################################################################

# Focus on complementary relation in Budyko space

# Budyko curves

out_BC_alpha <- Budyko_curve(PET = "PET_PT", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_alpha1.26 <- Budyko_curve(PET = "PET_PT_1.26", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                 X_range = c(0.3, 1.5), Y_range = c(0.4, 1), Xin_range = c(0.55, 0.951), Yin_range = c(0.6, 0.80))

out_BC_Zhou <- Budyko_curve(PET = "PET_PT_Zhou", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_Morton <- Budyko_curve(PET = "PET_PT_Morton", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                              X_range = c(0.4, 1.8), Y_range = c(0.4, 1), Xin_range = c(0.85, 1.15), Yin_range = c(0.7, 0.80))

out_BC_Makkink_Hansen <- Budyko_curve(PET = "PET_Makkink_Hansen", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_deBruin_2016 <- Budyko_curve(PET = "PET_deBruin_2016", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                    Xin_range = c(0.7, 1.1), Yin_range = c(0.6, 0.80))

out_BC_deBruin_1979 <- Budyko_curve(PET = "PET_deBruin_1979", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_FAO56_grass <- Budyko_curve(PET = "ETo_FAO56_grass", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_FAO56_alfalfa <- Budyko_curve(PET = "ETo_FAO56_alfalfa", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                     X_range = c(0.0, 2.0), Y_range = c(0.0, 1.06), Xin_range = c(0.9, 1.50), Yin_range = c(0.6, 0.8))

out_BC_FAO56_alfalfa_CO2_corr <- Budyko_curve(PET = "ETo_FAO56_alfalfa_CO2_corr", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                              X_range = c(0.0, 2.0), Y_range = c(0.0, 1.06), Xin_range = c(0.9, 1.50), Yin_range = c(0.6, 0.8))

out_BC_ET_PM_rs_shift_CO2_corr <- Budyko_curve(PET = "ETo_FAO56_alfalfa", ET = "ET_PM_rs_shift_CO2_corr", trace = FALSE, print_n = TRUE, plot=TRUE,
                                              X_range = c(0.0, 2.0), Y_range = c(0.0, 1.06), Xin_range = c(0.9, 1.50), Yin_range = c(0.55, 0.8))

out_BC_FAO56_alfalfa

# Add label "a)" to the first plot
p9 <- out_BC_FAO56_alfalfa$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "a)", hjust = 1.5, vjust = 1.5, size = 6)

# Add label "b)" to the second plot
p10 <- out_BC_FAO56_alfalfa_CO2_corr$out_plot + 
  annotate("text", x = Inf, y = Inf, label = "b)", hjust = 1.5, vjust = 1.5, size = 6)

# Combine the two plots side by side using cowplot
panel_figure <- (p9 | p10) 

# Save the combined plot
ggsave('../Plots/ggplot2/BC_PET_PM_and_PM_corr.png', plot = panel_figure, width = Pl_width*2, height = Pl_height, dpi = RES, units = 'mm')


# Complementary hypothesis
out_BC_alpha_CH <- Budyko_curve(PET = "PET_PT_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                X_range = c(0.3, 1.3), Y_range = c(0.4, 1), Xin_range = c(0.6, 1.3), Yin_range = c(0.6, 0.80))

out_BC_alpha1.26_CH <- Budyko_curve(PET = "PET_PT_1.26_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE)

out_BC_Morton_CH <- Budyko_curve(PET = "PET_PT_Morton_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                 X_range = c(0.4, 1.8), Y_range = c(0.4, 1), Xin_range = c(0.8, 1.6), Yin_range = c(0.6, 0.80))

out_BC_Zhou_CH <- Budyko_curve(PET = "PET_PT_Zhou_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                X_range = c(0.3, 1.3), Y_range = c(0.4, 1), Xin_range = c(0.6, 1.3), Yin_range = c(0.6, 0.80))

out_BC_Makkink_Hansen_CH <- Budyko_curve(PET = "PET_Makkink_Hansen_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                               X_range = c(0.3, 1.3), Y_range = c(0.4, 1), Xin_range = c(0.6, 1.4), Yin_range = c(0.6, 0.80))

out_BC_deBruin_2016_CH <- Budyko_curve(PET = "PET_deBruin_2016_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                         X_range = c(0.3, 1.5), Y_range = c(0.3, 1), Xin_range = c(0.8, 1.5), Yin_range = c(0.6, 0.80))

out_BC_deBruin_1979_CH <- Budyko_curve(PET = "PET_deBruin_1979_CH", ET = "ET", trace = FALSE, print_n = TRUE, plot=TRUE,
                                       X_range = c(0.3, 1.3), Y_range = c(0.4, 1), Xin_range = c(0.6, 1.3), Yin_range = c(0.6, 0.80))

################################################################################

# Runoff

# Slope and intercept of the delta P over P to delta ET over ET relation



Q_EOBS_mHM <- 224.6725
P_EOBS <- 689.7999
ET_EOBS_mHM <- 466.8086


################################################################################
mean(P_CMIP5_1981_2005_a)
mean(ET_CMIP5_1981_2005_a)
mean(RO_CMIP5_1981_2005_a)

mean(P_CMIP6_1981_2005_a)
mean(ET_CMIP6_1981_2005_a)
mean(RO_CMIP6_1981_2005_a)

mean(P_RCMs_1981_2005_a)
mean(ET_RCMs_1981_2005_a)
mean(RO_RCMs_1981_2005_a)

mean(P_CMIP5_2076_2100_a)
mean(ET_CMIP5_2076_2100_a)
mean(RO_CMIP5_2076_2100_a)

mean(P_CMIP6_2076_2100_a)
mean(ET_CMIP6_2076_2100_a)
mean(RO_CMIP6_2076_2100_a)

mean(P_RCMs_2076_2100_a)
mean(ET_RCMs_2076_2100_a)
mean(RO_RCMs_2076_2100_a)
     
# Create the data frames with the mean and standard deviation values for each dataset
Data_P <- data.frame(
  ModelGroup = rep(c("CMIP5", "CMIP6", "RCMs"), 2),
  Period = rep(c("Baseline", "Future"), each = 3),
  Mean = c(mean(P_CMIP5_1981_2005_a),
           mean(P_CMIP6_1981_2005_a),
           mean(P_RCMs_1981_2005_a),
           mean(P_CMIP5_2076_2100_a),
           mean(P_CMIP6_2076_2100_a),
           mean(P_RCMs_2076_2100_a)),
  SD = c(sd(P_CMIP5_1981_2005_a),
         sd(P_CMIP6_1981_2005_a),
         sd(P_RCMs_1981_2005_a),
         sd(P_CMIP5_2076_2100_a),
         sd(P_CMIP6_2076_2100_a),
         sd(P_RCMs_2076_2100_a))
)

Data_ET <- data.frame(
  ModelGroup = rep(c("CMIP5", "CMIP6", "RCMs"), 2),
  Period = rep(c("Baseline", "Future"), each = 3),
  Mean = c(mean(ET_CMIP5_1981_2005_a),
           mean(ET_CMIP6_1981_2005_a),
           mean(ET_RCMs_1981_2005_a),
           mean(ET_CMIP5_2076_2100_a),
           mean(ET_CMIP6_2076_2100_a),
           mean(ET_RCMs_2076_2100_a)),
  SD = c(sd(ET_CMIP5_1981_2005_a),
         sd(ET_CMIP6_1981_2005_a),
         sd(ET_RCMs_1981_2005_a),
         sd(ET_CMIP5_2076_2100_a),
         sd(ET_CMIP6_2076_2100_a),
         sd(ET_RCMs_2076_2100_a))
)

Data_RO <- data.frame(
  ModelGroup = rep(c("CMIP5", "CMIP6", "RCMs"), 2),
  Period = rep(c("Baseline", "Future"), each = 3),
  Mean = c(mean(RO_CMIP5_1981_2005_a),
           mean(RO_CMIP6_1981_2005_a),
           mean(RO_RCMs_1981_2005_a),
           mean(RO_CMIP5_2076_2100_a),
           mean(RO_CMIP6_2076_2100_a),
           mean(RO_RCMs_2076_2100_a)),
  SD = c(sd(RO_CMIP5_1981_2005_a),
         sd(RO_CMIP6_1981_2005_a),
         sd(RO_RCMs_1981_2005_a),
         sd(RO_CMIP5_2076_2100_a),
         sd(RO_CMIP6_2076_2100_a),
         sd(RO_RCMs_2076_2100_a))
)

# Function to create bar plot
create_barplot <- function(data, ylabel, y_min = NULL, y_max = NULL) {
  p <- ggplot(data, aes(x = Period, y = Mean, fill = ModelGroup)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD, color = ModelGroup),  # T-shaped error bars (above bars only)
      position = position_dodge(0.7),  # Ensure error bars are aligned with the bars
      width = 0.2,  # Width of the horizontal part of the T
      size = 0.15  # Thickness of the error bars
    ) +
    scale_fill_manual(values = c(CMIP5 = "#1c91df", CMIP6 = "#4e79a7", RCMs = "#ff4500")) +
    scale_color_manual(values = c(CMIP5 = "#333333", CMIP6 = "#333333", RCMs = "#333333")) +
    xlab(NULL) +
    ylab(ylabel) +
    scale_x_discrete(labels = c("1981–2005", "2076–2100")) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 15)),  # Move y-axis label to the left
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",  # No legend
      panel.grid.major.x = element_blank()
    )
  
  # Apply custom y-axis limits if provided
  if (!is.null(y_min) & !is.null(y_max)) {
    p <- p + ylim(y_min, y_max)
  }
  
  return(p)
}

# Create each of the three plots with custom y-axis limits
plot_P <- create_barplot(Data_P, expression(P ~ (mm ~ yr^{-1})), y_min = 0, y_max = 1200)
plot_ET <- create_barplot(Data_ET, expression(ET ~ (mm ~ yr^{-1})), y_min = 0, y_max = 800)
plot_RO <- create_barplot(Data_RO, expression(R[o] ~ (mm ~ yr^{-1})), y_min = 0, y_max = 500)

# Load the map PNG
map_image <- image_read("../Elevation_map/plots/Elevation_and_domain_map_r-6_v-01_DPI-300.png")

# Function to crop image based on percentage cropping from each side
crop_image_by_percentage <- function(image, left_percent, right_percent, top_percent, bottom_percent) {
  # Get the dimensions of the image
  image_info <- image_info(image)
  image_width <- image_info$width
  image_height <- image_info$height
  
  # Calculate the pixel values for each side based on the percentages
  left_crop <- (left_percent / 100) * image_width
  right_crop <- (right_percent / 100) * image_width
  top_crop <- (top_percent / 100) * image_height
  bottom_crop <- (bottom_percent / 100) * image_height
  
  # Calculate the new width and height after cropping
  cropped_width <- image_width - left_crop - right_crop
  cropped_height <- image_height - top_crop - bottom_crop
  
  # Round pixel values to integers for the cropping geometry
  left_crop <- round(left_crop)
  top_crop <- round(top_crop)
  cropped_width <- round(cropped_width)
  cropped_height <- round(cropped_height)
  
  # Crop the image using the calculated pixel values
  cropped_image <- image_crop(image, geometry = paste0(cropped_width, "x", cropped_height, "+", left_crop, "+", top_crop))
  
  # Return the cropped image
  return(cropped_image)
}

cropped_map <- crop_image_by_percentage(map_image, 6.5, 9.5, 5.6, 3)
# image_write(cropped_map, path = "../Plots/ggplot2/cropped_test.png")

# Convert the cropped image to raster format
cropped_map_raster <- as.raster(cropped_map)

# Convert the raster to a grob
cropped_map_grob <- rasterGrob(cropped_map_raster, interpolate = TRUE)

# Create a custom legend
legend_grob <- grid::grobTree(
  grid::rectGrob(gp = gpar(fill = "#1c91df"), width = 0.05, height = 0.3, x = 0.20, y = 0.2),
  grid::textGrob("CMIP5", x = 0.25, y = 0.2, hjust = 0, gp = gpar(fontsize = 12)),
  
  grid::rectGrob(gp = gpar(fill = "#4e79a7"), width = 0.05, height = 0.3, x = 0.50, y = 0.2),
  grid::textGrob("CMIP6", x = 0.55, y = 0.2, hjust = 0, gp = gpar(fontsize = 12)),
  
  grid::rectGrob(gp = gpar(fill = "#ff4500"), width = 0.05, height = 0.3, x = 0.75, y = 0.2),
  grid::textGrob("EUR-44", x = 0.80, y = 0.2, hjust = 0, gp = gpar(fontsize = 12))
)

# Combine the three plots vertically and add the custom legend below them
combined_plots_with_legend <- grid.arrange(
  plot_P, plot_ET, plot_RO,
  legend_grob,
  ncol = 1,
  heights = c(1, 1, 1, 0.2)
)

# Add empty space (nullGrob) above, below, and to the right of the combined plots
combined_with_padding <- grid.arrange(
  nullGrob(),  # Empty space above
  combined_plots_with_legend,
  nullGrob(),  # Empty space below
  ncol = 1,  # Stack everything vertically (1 column)
  heights = c(0.05, 1, 0.025)  # Maintain padding for the rows
)

# Add extra empty space to the right
combined_with_padding_and_right_margin <- grid.arrange(
  combined_with_padding,  # Main plot with top and bottom padding
  nullGrob(),  # Empty space to the right
  ncol = 2,  # 2 columns (plot and right margin)
  widths = c(0.95, 0.05)  # Adjust width proportion, keeping the plot at 95% width
)

# Add empty space (nullGrob) above, below, and to the right of the map
map_with_padding <- grid.arrange(
  nullGrob(),  # Empty space above
  cropped_map_grob,
  nullGrob(),  # Empty space below
  ncol = 1,  # Stack everything vertically (1 column)
  heights = c(0.05, 1, 0.025)  # Maintain padding for the rows
)

# Add labels to the map and plots
map_with_label <- grobTree(
  map_with_padding,
  grid.text("a)", x = 0.05, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
)

combined_with_padding_and_label <- grobTree(
  combined_with_padding_and_right_margin,
  grid.text("b)", x = -0.02, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
)

# Final arrangement: map on the left, plots on the right
final_figure <- grid.arrange(
  map_with_label, combined_with_padding_and_label, 
  ncol = 2,
  widths = c(2, 1)
)

# Save or display the final figure
ggsave("../Plots/ggplot2/Elevation_and_hydroclimate.png", plot = final_figure, width = Pl_width*2.55, height = Pl_height*1.8, , dpi = RES, units = "mm")

################################################################################




















CMIP5_out <- data.frame(
  d_P_over_P=d_P_over_P_CMIP5_a,
  P_mm=annSums(P_CMIP5_1981_2005[,CMIP5_range]),
  d_ET_over_ET=d_ET_over_ET_CMIP5_a,
  ET_mm=annSums(ET_CMIP5_1981_2005[,CMIP5_range]),
  d_Ro_over_Ro=((annSums(P_CMIP5_2076_2100[,CMIP5_range])-annSums(ET_CMIP5_2076_2100[,CMIP5_range])) -
  (annSums(P_CMIP5_1981_2005[,CMIP5_range])-annSums(ET_CMIP5_1981_2005[,CMIP5_range]))) / (annSums(P_CMIP5_1981_2005[,CMIP5_range])-annSums(ET_CMIP5_1981_2005[,CMIP5_range])),
  Ro_mm=annSums(P_CMIP5_1981_2005[,CMIP5_range])-annSums(ET_CMIP5_1981_2005[,CMIP5_range]),
  PET=PET_PT_CMIP5_a_1981_2005,
  d_PET=PET_PT_CMIP5_a_2076_2100-PET_PT_CMIP5_a_1981_2005,
  alpha=ifelse(as.vector(t(alpha_PT_CMIP5_1981_2005[1,CMIP5_range]))<alpha_PT,alpha_PT,as.vector(t(alpha_PT_CMIP5_1981_2005[1,CMIP5_range]))),
  SVP=annMean(SVP_CMIP5_1981_2005[,CMIP5_range]),
  d_SVP=annMean(SVP_CMIP5_2076_2100[,CMIP5_range])-annMean(SVP_CMIP5_1981_2005[,CMIP5_range]),
  Bo=Bo_CMIP5_1981_2005_a,
  d_Bo=Bo_CMIP5_2076_2100_a-Bo_CMIP5_1981_2005_a,
  d_Rn_mm=(annMean(Rn_CMIP5_2076_2100[,CMIP5_range])-annMean(Rn_CMIP5_1981_2005[,CMIP5_range]))*3600*24*365.25/10^6*0.408,
  omega=rep(out_BC_alpha$omega_CMIP5_1981_2005,length(CMIP5_range)),
  d_omega=rep(out_BC_alpha$omega_CMIP5_2076_2100-out_BC_alpha$omega_CMIP5_1981_2005,length(CMIP5_range))
)

CMIP6_out <- data.frame(
  d_P_over_P=d_P_over_P_CMIP6_a,
  P_mm=annSums(P_CMIP6_1981_2005[,CMIP6_range]),
  d_ET_over_ET=d_ET_over_ET_CMIP6_a,
  ET_mm=annSums(ET_CMIP6_1981_2005[,CMIP6_range]),
  d_Ro_over_Ro=((annSums(P_CMIP6_2076_2100[,CMIP6_range])-annSums(ET_CMIP6_2076_2100[,CMIP6_range])) -
                  (annSums(P_CMIP6_1981_2005[,CMIP6_range])-annSums(ET_CMIP6_1981_2005[,CMIP6_range]))) / (annSums(P_CMIP6_1981_2005[,CMIP6_range])-annSums(ET_CMIP6_1981_2005[,CMIP6_range])),
  Ro_mm=annSums(P_CMIP6_1981_2005[,CMIP6_range])-annSums(ET_CMIP6_1981_2005[,CMIP6_range]),
  PET=PET_PT_CMIP6_a_1981_2005,
  d_PET=PET_PT_CMIP6_a_2076_2100-PET_PT_CMIP6_a_1981_2005,
  alpha=ifelse(as.vector(t(alpha_PT_CMIP6_1981_2005[1,CMIP6_range]))<alpha_PT,alpha_PT,as.vector(t(alpha_PT_CMIP6_1981_2005[1,CMIP6_range]))),
  SVP=annMean(SVP_CMIP6_1981_2005[,CMIP6_range]),
  d_SVP=annMean(SVP_CMIP6_2076_2100[,CMIP6_range])-annMean(SVP_CMIP6_1981_2005[,CMIP6_range]),
  Bo=Bo_CMIP6_1981_2005_a,
  d_Bo=Bo_CMIP6_2076_2100_a-Bo_CMIP6_1981_2005_a,
  d_Rn_mm=(annMean(Rn_CMIP6_2076_2100[,CMIP6_range])-annMean(Rn_CMIP6_1981_2005[,CMIP6_range]))*3600*24*365.25/10^6*0.408,
  omega=rep(out_BC_alpha$omega_CMIP6_1981_2005,length(CMIP6_range)),
  d_omega=rep(out_BC_alpha$omega_CMIP6_2076_2100-out_BC_alpha$omega_CMIP6_1981_2005,length(CMIP6_range))
)

RCMs_out <- data.frame(
  d_P_over_P=d_P_over_P_RCMs_a,
  P_mm=annSums(P_RCMs_1981_2005[,RCMs_range]),
  d_ET_over_ET=d_ET_over_ET_RCMs_a,
  ET_mm=annSums(ET_RCMs_1981_2005[,RCMs_range]),
  d_Ro_over_Ro=((annSums(P_RCMs_2076_2100[,RCMs_range])-annSums(ET_RCMs_2076_2100[,RCMs_range])) -
                  (annSums(P_RCMs_1981_2005[,RCMs_range])-annSums(ET_RCMs_1981_2005[,RCMs_range]))) / (annSums(P_RCMs_1981_2005[,RCMs_range])-annSums(ET_RCMs_1981_2005[,RCMs_range])),
  Ro_mm=annSums(P_RCMs_1981_2005[,RCMs_range])-annSums(ET_RCMs_1981_2005[,RCMs_range]),
  PET=PET_PT_RCMs_a_1981_2005,
  d_PET=PET_PT_RCMs_a_2076_2100-PET_PT_RCMs_a_1981_2005,
  alpha=ifelse(as.vector(t(alpha_PT_RCMs_1981_2005[1,RCMs_range]))<alpha_PT,alpha_PT,as.vector(t(alpha_PT_RCMs_1981_2005[1,RCMs_range]))),
  SVP=annMean(SVP_RCMs_1981_2005[,RCMs_range]),
  d_SVP=annMean(SVP_RCMs_2076_2100[,RCMs_range])-annMean(SVP_RCMs_1981_2005[,RCMs_range]),
  Bo=Bo_RCMs_1981_2005_a,
  d_Bo=Bo_RCMs_2076_2100_a-Bo_RCMs_1981_2005_a,
  d_Rn_mm=(annMean(Rn_RCMs_2076_2100[,RCMs_range])-annMean(Rn_RCMs_1981_2005[,RCMs_range]))*3600*24*365.25/10^6*0.408,
  omega=rep(out_BC_alpha$omega_RCMs_1981_2005,length(RCMs_range)),
  d_omega=rep(out_BC_alpha$omega_RCMs_2076_2100-out_BC_alpha$omega_RCMs_1981_2005,length(RCMs_range))
)

write.table(data.frame('Model'=rownames(CMIP5_out),CMIP5_out), 'CMIP5_out.csv', sep=',',row.names=FALSE)
write.table(data.frame('Model'=rownames(CMIP6_out),CMIP6_out), 'CMIP6_out.csv', sep=',',row.names=FALSE)
write.table(data.frame('Model'=rownames(RCMs_out),RCMs_out), 'RCMs_out.csv', sep=',',row.names=FALSE)

################################################################################

# Runoff





################################################################################

mean((annSums(d_P_CMIP5[,CMIP5_range]) - annSums(d_ET_CMIP5[,CMIP5_range]))/annSums(P_CMIP5_1981_2005[,CMIP5_range]))

mean((annSums(d_P_RCMs[,RCMs_range]) - annSums(d_ET_RCMs[,RCMs_range]))/annSums(P_RCMs_1981_2005[,RCMs_range]))

mean(annSums(P_CMIP5_1981_2005[,CMIP5_range]))
mean(annSums(P_RCMs_1981_2005[,RCMs_range]))

mean(annSums(P_CMIP5_2076_2100[,CMIP5_range]))
mean(annSums(P_RCMs_2076_2100[,RCMs_range]))

mean(annSums(ET_CMIP5_1981_2005[,CMIP5_range]))
mean(annSums(ET_RCMs_1981_2005[,RCMs_range]))

mean(annSums(ET_CMIP5_2076_2100[,CMIP5_range]))
mean(annSums(ET_RCMs_2076_2100[,RCMs_range]))

mean(annSums(PET_PT_CMIP5_1981_2005[,CMIP5_range]),na.rm=TRUE)
mean(annSums(PET_PT_RCMs_1981_2005[,RCMs_range]),na.rm=TRUE)

mean(annSums(PET_PT_CMIP5_2076_2100[,CMIP5_range]),na.rm=TRUE)
mean(annSums(PET_PT_RCMs_2076_2100[,RCMs_range]),na.rm=TRUE)

mean(annMean(Rn_CMIP5_1981_2005[,CMIP5_range]),na.rm=TRUE)
mean(annMean(Rn_RCMs_1981_2005[,RCMs_range]),na.rm=TRUE)

mean(annMean(Rn_CMIP5_2076_2100[,CMIP5_range]),na.rm=TRUE)
mean(annMean(Rn_RCMs_2076_2100[,RCMs_range]),na.rm=TRUE)

mean(annMean(Ta_CMIP5_1981_2005[,CMIP5_range]),na.rm=TRUE)
mean(annMean(Ta_RCMs_1981_2005[,RCMs_range]),na.rm=TRUE)

mean(annMean(Ta_CMIP5_2076_2100[,CMIP5_range]),na.rm=TRUE)
mean(annMean(Ta_RCMs_2076_2100[,RCMs_range]),na.rm=TRUE)

mean(annMean(SVP_CMIP5_1981_2005[,CMIP5_range]),na.rm=TRUE)
mean(annMean(SVP_RCMs_1981_2005[,RCMs_range]),na.rm=TRUE)

mean(annMean(SVP_CMIP5_2076_2100[,CMIP5_range]),na.rm=TRUE)
mean(annMean(SVP_RCMs_2076_2100[,RCMs_range]),na.rm=TRUE)

mean(annMean(LW_net_CMIP5_1981_2005[,CMIP5_range]),na.rm=TRUE)
mean(annMean(LW_net_RCMs_1981_2005[,RCMs_range]),na.rm=TRUE)

mean(annMean(LW_net_CMIP5_2076_2100[,CMIP5_range]),na.rm=TRUE)
mean(annMean(LW_net_RCMs_2076_2100[,RCMs_range]),na.rm=TRUE)
