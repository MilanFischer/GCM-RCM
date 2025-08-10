#-------------------------------------------------------------------------------
# Function to estimate the Budyko curve parameter 'n' using optimization
# 
# Inputs:
# PET     - Potential evapotranspiration (numeric vector)
# ET      - Actual evapotranspiration (numeric vector)
# P       - Precipitation (numeric vector)
# maxiter - Maximum number of iterations for the DEoptim algorithm (default = 500)
# lower   - Lower bound for the 'n' parameter (default = 0)
# upper   - Upper bound for the 'n' parameter (default = 1000)
#
# This function computes the Budyko 'n' parameter by minimizing the sum of squared 
# errors between observed evaporative index (ET/P) and the theoretical Budyko function 
# proposed by Choudhury (1999). Optimization is performed using differential evolution (DEoptim).
#
# Returns:
# Estimated value of the Budyko 'n' parameter that best fits the observed data.
#-------------------------------------------------------------------------------

Budyko_curve_optim <- function(PET, ET, P, maxiter = 500, lower = 0, upper = 1000){
  
  # Aridity index
  X <- PET / P
  
  # Evaporative index
  Y <- ET / P
  
  # Cost function
  Prediction <- function(pars, X, Y){
    n <- pars[1]
    Choudhury_EI <- 1 / (1 + (1 / X)^n)^(1 / n)
    sum((Y - Choudhury_EI)^2, na.rm = TRUE)
  }
  
  # Suppress warnings from DEoptim internally
  Out <- suppressWarnings(
    DEoptim(
      fn = Prediction,
      lower = lower,
      upper = upper,
      DEoptim.control(itermax = maxiter, trace = FALSE),
      X, Y
    )
  )
  
  return(Out$optim$bestmem)
}


#-------------------------------------------------------------------------------
# Budyko_curve: Estimates Budyko 'n' parameters across model ensembles and visualizes hydroclimatic shifts
#
# Description:
# This function computes the Budyko curve parameter 'n' for different climate model ensembles 
# and time periods, using annual PET, ET, and precipitation data. It then creates a comprehensive 
# diagnostic plot comparing historical and future shifts in evaporative index (ET/P) versus 
# aridity index (PET/P), including both model-level and ensemble-level trajectories.
#
# Arguments:
# - annual_stats_wide : Data frame with annual statistics, including PET, ET, and P columns
# - pet_col, et_col, p_col : Strings specifying column names for PET, ET, and P (default = "ETo_FAO56_alfalfa", "ET", "P")
# - X_range, Y_range : Axis ranges for the main plot (defaults cover [0,2] for PET/P and [0,1.06] for ET/P)
# - Xin_range, Yin_range : Axis ranges for the inset plot showing ensemble means
# - plot : Logical; whether to generate and save the plot (default = TRUE)
# - plot_name : Filename to save the main plot (default = "Budyko_curve_ggplot2.png")
#
# Returns:
# A list containing:
# - Budyko_pars : Data frame with estimated 'n' and derived 'omega' values per ENSEMBLE–PERIOD combination
# - out_plot    : If plot = TRUE, the ggplot object is returned invisibly; otherwise returns NA
#
# Notes:
# - The Budyko 'n' parameter is estimated via the DEoptim algorithm.
# - Colors and line styles are customized by ENSEMBLE and PERIOD grouping.
# - The function suppresses direct plot printing and writes to file silently.
# - Assumes presence of globals: COL_CMIP5, COL_CMIP6, COL_RCMs, Pl_width, Pl_height, RES.
#-------------------------------------------------------------------------------

Budyko_curve <- function(annual_stats_wide,
                         pet_col = "ETo_FAO56_alfalfa", et_col = "ET", p_col = "P",
                         X_range = c(0, 2), Y_range = c(0, 1.06),
                         Xin_range = c(0.9, 1.5), Yin_range = c(0.6, 0.8),
                         plot = TRUE, plot_name = "Budyko_curve_ggplot2.png") {
  
  # -- Axis buffer extensions
  X_range <- X_range + c(-diff(X_range) * 0.04, diff(X_range) * 0.04)
  Y_range <- Y_range + c(-diff(Y_range) * 0.04, diff(Y_range) * 0.04)
  
  # -- Symbolic column access
  pet_sym <- sym(pet_col)
  et_sym  <- sym(et_col)
  p_sym   <- sym(p_col)
  
  # -- Estimate Budyko 'n' parameter per ENSEMBLE and PERIOD
  n_estimates_ens <- annual_stats_wide |>
    group_by(ENSEMBLE, PERIOD) |>
    reframe(n = purrr::possibly(Budyko_curve_optim, otherwise = NA_real_)(
      PET = !!pet_sym, ET = !!et_sym, P = !!p_sym
    )) |>
    mutate(
      ENSEMBLE = recode(ENSEMBLE, "EUR-44" = "RCMs"),
      Group = paste0(ENSEMBLE, "_", PERIOD)
    )
  
  # -- Compute PET/P and ET/P ratios
  pet_et_ratios <- annual_stats_wide |>
    mutate(
      ENSEMBLE = recode(ENSEMBLE, "EUR-44" = "RCMs"),
      PET_over_P = !!pet_sym / !!p_sym,
      ET_over_P  = !!et_sym  / !!p_sym
    ) |>
    filter(PERIOD %in% c("1981_2005", "2076_2100")) |>
    select(ENSEMBLE, PERIOD, MODEL, PET_over_P, ET_over_P)
  
  X_vals <- seq(from = X_range[1], to = X_range[2], by = 0.001)
  
  lines_data <- n_estimates_ens |>
    mutate(
      LineType = if_else(grepl("1981_2005", Group), "solid", "dashed"),
      Color = case_when(
        grepl("^CMIP5", Group) ~ COL_CMIP5,
        grepl("^CMIP6", Group) ~ COL_CMIP6,
        grepl("^RCMs",  Group) ~ COL_RCMs
      )
    ) |>
    rowwise() |>
    mutate(
      X = list(X_vals),
      Y = list(1 / (1 + (1 / X_vals)^n)^(1 / n))
    ) |>
    unnest(c(X, Y)) |>
    rename(Scenario = Group)
  
  # -- Filter valid model entries for both periods and clean ensemble names
  valid_models <- pet_et_ratios |>
    count(ENSEMBLE, MODEL) |>
    filter(n == 2) |>
    select(ENSEMBLE, MODEL)
  
  Data <- pet_et_ratios |> semi_join(valid_models, by = c("ENSEMBLE", "MODEL"))
  
  # -- Construct model-wise shift data
  arrow_data <- Data |>
    pivot_wider(names_from = PERIOD, values_from = c(PET_over_P, ET_over_P)) |>
    mutate(
      x = PET_over_P_1981_2005,
      y = ET_over_P_1981_2005,
      xend = PET_over_P_2076_2100,
      yend = ET_over_P_2076_2100,
      Scenario = ENSEMBLE,
      Color = case_when(
        Scenario == "CMIP5" ~ COL_CMIP5,
        Scenario == "CMIP6" ~ COL_CMIP6,
        Scenario == "RCMs"  ~ COL_RCMs
      ),
      arrow_type = if_else(
        Scenario == "RCMs" | MODEL %in% driving_GCMs,
        "thick", "thin"
      ),
      line_size = if_else(arrow_type == "thick", 0.6, 0.2)
    ) |>
    select(x, y, xend, yend, Scenario, Color, Model = MODEL, arrow_type, line_size)
  
  # -- Ensemble means (all models)
  summarize_ensemble <- function(df, sel) {
    df |>
      group_by(ENSEMBLE, PERIOD) |>
      summarise(
        PET_over_P = mean(PET_over_P, na.rm = TRUE),
        ET_over_P  = mean(ET_over_P,  na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(Group = paste0(ENSEMBLE, "_", PERIOD), SEL = sel)
  }
  
  Data_ens_all   <- summarize_ensemble(Data, sel = "all")
  Data_ens_drive <- summarize_ensemble(
    Data |> filter(ENSEMBLE %in% c("CMIP5", "CMIP6"), MODEL %in% driving_GCMs),
    sel = "drive"
  )
  
  Data_ens <- bind_rows(Data_ens_all, Data_ens_drive)
  
  # -- Ensemble-level arrows
  arrow_data_ens <- tibble(
    x = Data_ens$PET_over_P[c(1,3,5,7,9)],
    y = Data_ens$ET_over_P[c(1,3,5,7,9)],
    xend = Data_ens$PET_over_P[c(2,4,6,8,10)],
    yend = Data_ens$ET_over_P[c(2,4,6,8,10)],
    Scenario = c("CMIP5", "CMIP6", "RCMs", "CMIP5", "CMIP6"),
    Color = c(COL_CMIP5, COL_CMIP6, COL_RCMs, COL_CMIP5, COL_CMIP6)
  )
  
  arrow_data_thick <- arrow_data |> filter(arrow_type == "thick")
  arrow_data_thin  <- arrow_data |> filter(arrow_type == "thin")
  
  # -- Text annotation helpers
  get_n <- function(ens, per) {
    n_estimates_ens |> filter(ENSEMBLE == ens, PERIOD == per) |> pull(n)
  }
  
  # -- Generate annotation labels for omega and delta omega values
  # These are mathematical expressions to be displayed in the plot.
  # `bquote()` creates the expression, and `deparse()` converts it to a character string
  # so it can be used with ggplot2::annotate() and `parse = TRUE` without triggering warnings.
  txt_CMIP5_1981_2005 <- deparse(bquote(omega == .(sprintf("%.2f", get_n("CMIP5", "1981_2005") + 0.72))))
  txt_CMIP6_1981_2005 <- deparse(bquote(omega == .(sprintf("%.2f", get_n("CMIP6", "1981_2005") + 0.72))))
  txt_RCMs_1981_2005  <- deparse(bquote(omega == .(sprintf("%.2f", get_n("RCMs", "1981_2005") + 0.72))))
  
  txt_CMIP5_delta_omega <- deparse(bquote(Delta * omega == .(sprintf("%.2f", get_n("CMIP5", "2076_2100") - get_n("CMIP5", "1981_2005")))))
  txt_CMIP6_delta_omega <- deparse(bquote(Delta * omega == .(sprintf("%.2f", get_n("CMIP6", "2076_2100") - get_n("CMIP6", "1981_2005")))))
  txt_RCMs_delta_omega  <- deparse(bquote(Delta * omega == .(sprintf("%.2f", get_n("RCMs", "2076_2100") - get_n("RCMs", "1981_2005")))))
  
  # -- Main plot
  if (plot) {
    boundary_line <- tibble(x = c(X_range[1], 1, X_range[2]), y = c(Y_range[1], 1, 1))
    p <- ggplot() +
      
      # Add continuous boundary line from (0,0) to (1,1) and then horizontally at y = 1
      geom_line(data = boundary_line, aes(x = x, y = y), color = "gray", linetype = "solid", size = 0.1, na.rm = TRUE) +
      
      # Add continuous Budyko curves
      geom_line(data = lines_data, aes(x = X, y = Y, group = Scenario, color = Scenario, linetype = LineType), na.rm = TRUE) +
      
      # Add arrows
      geom_segment(data = arrow_data_thin, aes(x = x, y = y, xend = xend, yend = yend, color = Scenario, size = line_size),
                   arrow = arrow(type = "closed", length = unit(0.03, "inches")), lineend = "round") +
      geom_segment(data = arrow_data_thick, aes(x = x, y = y, xend = xend, yend = yend, color = Scenario, size = line_size),
                   arrow = arrow(type = "closed", length = unit(0.09, "inches")), lineend = "round") +
      scale_color_manual(values = c(
        "CMIP5_1981_2005" = COL_CMIP5, "CMIP6_1981_2005" = COL_CMIP6, "RCMs_1981_2005" = COL_RCMs,
        "CMIP5_2076_2100" = COL_CMIP5, "CMIP6_2076_2100" = COL_CMIP6, "RCMs_2076_2100" = COL_RCMs,
        "RCMs" = COL_RCMs, "CMIP5" = COL_CMIP5, "CMIP6" = COL_CMIP6
      )) +
      scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
      scale_size_continuous(range = c(0.2, 0.6)) +
      labs(x = expression(E[p] * "/" * P), y = expression("ET/P")) +
      scale_x_continuous(limits = X_range, expand = c(0, 0)) +
      scale_y_continuous(limits = Y_range, expand = c(0, 0)) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none")
    
    
    # Adding annotations
    p <- p +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.96, Y_range), label = "CMIP5", hjust = 0, color = COL_CMIP5, fontface = "bold") +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.91, Y_range), label = "CMIP6", hjust = 0, color = COL_CMIP6, fontface = "bold") +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.86, Y_range), label = "EUR-44", hjust = 0, color = COL_RCMs, fontface = "bold") +
      annotate("text", x = rescale(0.054, X_range), y = rescale(0.75, Y_range), label = txt_CMIP5_1981_2005, parse = TRUE, hjust = 0, color = COL_CMIP5, fontface = "bold") +
      annotate("text", x = rescale(0.054, X_range), y = rescale(0.69, Y_range), label = txt_CMIP6_1981_2005, parse = TRUE, hjust = 0, color = COL_CMIP6, fontface = "bold") +
      annotate("text", x = rescale(0.054, X_range), y = rescale(0.63, Y_range), label = txt_RCMs_1981_2005, parse = TRUE, hjust = 0, color = COL_RCMs, fontface = "bold") +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.57, Y_range), label = txt_CMIP5_delta_omega, parse = TRUE, hjust = 0, color = COL_CMIP5, fontface = "bold") +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.51, Y_range), label = txt_CMIP6_delta_omega, parse = TRUE, hjust = 0, color = COL_CMIP6, fontface = "bold") +
      annotate("text", x = rescale(0.03, X_range), y = rescale(0.45, Y_range), label = txt_RCMs_delta_omega, parse = TRUE, hjust = 0, color = COL_RCMs, fontface = "bold")
    
    # -- Inset plot
    inset_plot <- ggplot() +
      
      # Add arrows
      geom_segment(data = arrow_data_ens[1:2, ], aes(x = x, y = y, xend = xend, yend = yend, color = Scenario, size = 0.01),
                   arrow = arrow(type = "closed", length = unit(0.03, "inches")), lineend = "round") +
      geom_segment(data = arrow_data_ens[3:5, ], aes(x = x, y = y, xend = xend, yend = yend, color = Scenario, size = 0.03),
                   arrow = arrow(type = "closed", length = unit(0.09, "inches")), lineend = "round") +
      scale_color_manual(values = c(
        "CMIP5" = COL_CMIP5, "CMIP6" = COL_CMIP6, "RCMs" = COL_RCMs,
        "CMIP5" = COL_CMIP5, "CMIP6" = COL_CMIP6)) +
      scale_size_continuous(range = c(0.2, 0.6)) +
      labs(x = expression(E[p] * "/" * P), y = expression("ET/P")) +
      scale_x_continuous(limits = Xin_range, expand = c(0, 0)) +
      scale_y_continuous(limits = Yin_range, expand = c(0, 0)) +
      ggtitle("Ensemble means") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),  # Adjust the plot title size and center it
        axis.title.x = element_text(size = 8),  # Adjust x-axis title size
        axis.title.y = element_text(size = 8),  # Adjust y-axis title size
        axis.text.x = element_text(size = 8),   # Adjust x-axis text size
        axis.text.y = element_text(size = 8),   # Adjust y-axis text size
        legend.position = "none"
      ) +
      coord_fixed(ratio = (X_range[2]-X_range[1])/(Y_range[2]-Y_range[1])) # This is important to preserve the ratio of x and y axes
    
    # Combining plots
    combined_plot <- p +
      annotation_custom(
        grob = ggplotGrob(inset_plot),
        xmin = rescale(0.40, X_range), xmax = rescale(0.98, X_range), 
        ymin = rescale(0, Y_range), ymax = rescale(0.4, Y_range)
      )
    
    # Save silently
    ggsave(plot_name, plot = combined_plot, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")
  }
  
  # -- Return results
  list(
    Budyko_pars = n_estimates_ens |>
      mutate(ENSEMBLE = recode(ENSEMBLE, "RCMs" = "EUR-44"), omega = n + 0.72) |>
      select(-Group),
    out_plot = if (plot) invisible(combined_plot) else NA
  )
}
