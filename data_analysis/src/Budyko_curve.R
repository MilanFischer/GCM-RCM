################################################################################
# Fitting Budyko curve

Prediction <- function(pars, X, Y){
  n <- pars[1]
  Choudhury_EI <- 1/(1+(1/X)^n)^(1/n)
  sum((Y-Choudhury_EI)^2, na.rm=TRUE)
}

Budyko_curve_optim <- function(model, period, PET, ET, maxiter = 200, lower = 0, upper = 50, trace = FALSE, print_n = FALSE){
  
  # Aridity index
  X <- get(paste0(PET, "_over_P_", model, "_", period, "_a"))
  
  # Evaporative index
  Y <- get(paste0(ET, "_over_P_", model, "_", period, "_a"))
  
  Out <- DEoptim(fn=Prediction,
                 lower = lower,
                 upper = upper,
                 DEoptim.control(itermax = maxiter, trace = trace), X, Y)
  
  assign(x=paste0("n_", PET, "_over_", ET, "_", model, "_", period, "_a"), value=as.numeric(Out$optim$bestmem), envir = parent.frame())
  
  if(print_n == TRUE){
    cat(paste0("\nn_", PET, "_over_", ET, "_", model, "_", period, "_a = ", round(as.numeric(Out$optim$bestmem), 3), "\n"))
  }
}

################################################################################

Budyko_curve <- function(PET, ET, maxiter=200, lower = 0, upper = 50, trace = FALSE, print_n = FALSE, plot=FALSE,
                         X_range = c(0.4, 1.4), Y_range = c(0.4, 1), Xin_range = c(0.65, 1.05), Yin_range = c(0.64, 0.78)){
  for(M in models){
    for(P in periods){
      Budyko_curve_optim(model = M, PET=PET, ET, period = P, maxiter=maxiter, lower = lower, upper = upper, trace = trace, print_n = print_n)
    }
  }
  
  PET_over_P_CMIP5_1981_2005_a <- get(paste0(PET, "_over_P_CMIP5_1981_2005_a"))
  PET_over_P_CMIP6_1981_2005_a <- get(paste0(PET, "_over_P_CMIP6_1981_2005_a"))
  PET_over_P_RCMs_1981_2005_a <- get(paste0(PET, "_over_P_RCMs_1981_2005_a"))
  PET_over_P_CMIP5_2076_2100_a <- get(paste0(PET, "_over_P_CMIP5_2076_2100_a"))
  PET_over_P_CMIP6_2076_2100_a <- get(paste0(PET, "_over_P_CMIP6_2076_2100_a"))
  PET_over_P_RCMs_2076_2100_a <- get(paste0(PET, "_over_P_RCMs_2076_2100_a"))
  
  ET_over_P_CMIP5_1981_2005_a <- get(paste0(ET, "_over_P_CMIP5_1981_2005_a"))
  ET_over_P_CMIP6_1981_2005_a <- get(paste0(ET, "_over_P_CMIP6_1981_2005_a"))
  ET_over_P_RCMs_1981_2005_a <- get(paste0(ET, "_over_P_RCMs_1981_2005_a"))
  ET_over_P_CMIP5_2076_2100_a <- get(paste0(ET, "_over_P_CMIP5_2076_2100_a"))
  ET_over_P_CMIP6_2076_2100_a <- get(paste0(ET, "_over_P_CMIP6_2076_2100_a"))
  ET_over_P_RCMs_2076_2100_a <- get(paste0(ET, "_over_P_RCMs_2076_2100_a"))
  
  n_CMIP5_1981_2005_a <- get(paste0("n_", PET, "_over_", ET, "_CMIP5_1981_2005_a"))
  n_CMIP6_1981_2005_a <- get(paste0("n_", PET, "_over_", ET, "_CMIP6_1981_2005_a"))
  n_RCMs_1981_2005_a <- get(paste0("n_", PET, "_over_", ET, "_RCMs_1981_2005_a"))
  n_CMIP5_2076_2100_a <- get(paste0("n_", PET, "_over_", ET, "_CMIP5_2076_2100_a"))
  n_CMIP6_2076_2100_a <- get(paste0("n_", PET, "_over_", ET, "_CMIP6_2076_2100_a"))
  n_RCMs_2076_2100_a <- get(paste0("n_", PET, "_over_", ET, "_RCMs_2076_2100_a"))
  
  ##############################################################################
  # ggplot figure
  
  if(plot == TRUE){
    plot_name <- paste0("../Plots/ggplot2/Budyko_", PET, "&", ET, "_ggplot2.png")
    
    # Generate the vector of aridity index
    X <- seq(from = 0.001, to = 5, by = 0.001)
    
    # Generate the vectors of evaporative index
    EI_CMIP5_1981_2005_a <- 1/(1+(1/X)^n_CMIP5_1981_2005_a)^(1/n_CMIP5_1981_2005_a)
    EI_CMIP6_1981_2005_a <- 1/(1+(1/X)^n_CMIP6_1981_2005_a)^(1/n_CMIP6_1981_2005_a)
    EI_RCMs_1981_2005_a <- 1/(1+(1/X)^n_RCMs_1981_2005_a)^(1/n_RCMs_1981_2005_a)
    EI_CMIP5_2076_2100_a <- 1/(1+(1/X)^n_CMIP5_2076_2100_a)^(1/n_CMIP5_2076_2100_a)
    EI_CMIP6_2076_2100_a <- 1/(1+(1/X)^n_CMIP6_2076_2100_a)^(1/n_CMIP6_2076_2100_a)
    EI_RCMs_2076_2100_a <- 1/(1+(1/X)^n_RCMs_2076_2100_a)^(1/n_RCMs_2076_2100_a)
    
    # Create a data frame for lines
    lines_data <- data.frame(
      X = rep(X, times = 6),
      Y = c(EI_CMIP5_1981_2005_a, EI_CMIP6_1981_2005_a, EI_RCMs_1981_2005_a,
            EI_CMIP5_2076_2100_a, EI_CMIP6_2076_2100_a, EI_RCMs_2076_2100_a),
      Scenario = rep(c("CMIP5_1981_2005", "CMIP6_1981_2005", "RCMs_1981_2005",
                       "CMIP5_2076_2100", "CMIP6_2076_2100", "RCMs_2076_2100"), each = length(X)),
      LineType = rep(c("solid", "solid", "solid", "dashed", "dashed", "dashed"), each = length(X)),
      Color = rep(c(COL_CMIP5, COL_CMIP6, COL_RCMs, COL_CMIP5, COL_CMIP6, COL_RCMs), each = length(X))
    )
    
    Data <- data.frame(
      PET_over_P = c(PET_over_P_RCMs_1981_2005_a, PET_over_P_CMIP5_1981_2005_a, PET_over_P_CMIP6_1981_2005_a,
                     PET_over_P_RCMs_2076_2100_a, PET_over_P_CMIP5_2076_2100_a, PET_over_P_CMIP6_2076_2100_a),
      ET_over_P = c(ET_over_P_RCMs_1981_2005_a, ET_over_P_CMIP5_1981_2005_a, ET_over_P_CMIP6_1981_2005_a,
                    ET_over_P_RCMs_2076_2100_a, ET_over_P_CMIP5_2076_2100_a, ET_over_P_CMIP6_2076_2100_a),
      Group = factor(rep(c("RCMs_1981_2005", "CMIP5_1981_2005", "CMIP6_1981_2005",
                           "RCMs_2076_2100", "CMIP5_2076_2100", "CMIP6_2076_2100"), 
                         times = c(length(PET_over_P_RCMs_1981_2005_a), length(PET_over_P_CMIP5_1981_2005_a), 
                                   length(PET_over_P_CMIP6_1981_2005_a), length(PET_over_P_RCMs_2076_2100_a), 
                                   length(PET_over_P_CMIP5_2076_2100_a), length(PET_over_P_CMIP6_2076_2100_a))))
    )
    
    arrow_data <- data.frame(
      x = c(PET_over_P_RCMs_1981_2005_a, PET_over_P_CMIP5_1981_2005_a, PET_over_P_CMIP6_1981_2005_a),
      y = c(ET_over_P_RCMs_1981_2005_a, ET_over_P_CMIP5_1981_2005_a, ET_over_P_CMIP6_1981_2005_a),
      xend = c(PET_over_P_RCMs_2076_2100_a, PET_over_P_CMIP5_2076_2100_a, PET_over_P_CMIP6_2076_2100_a),
      yend = c(ET_over_P_RCMs_2076_2100_a, ET_over_P_CMIP5_2076_2100_a, ET_over_P_CMIP6_2076_2100_a),
      Scenario = rep(c("RCMs", "CMIP5", "CMIP6"), times = c(length(PET_over_P_RCMs_1981_2005_a), length(PET_over_P_CMIP5_1981_2005_a), length(PET_over_P_CMIP6_1981_2005_a))),
      Color = rep(c(COL_RCMs, COL_CMIP5, COL_CMIP6), times = c(length(PET_over_P_RCMs_1981_2005_a), length(PET_over_P_CMIP5_1981_2005_a), length(PET_over_P_CMIP6_1981_2005_a))),
      Model = c(names(PET_over_P_RCMs_1981_2005_a),names(PET_over_P_CMIP5_1981_2005_a), names(PET_over_P_CMIP6_1981_2005_a))
    )
    
    arrow_data$arrow_type <- ifelse(arrow_data$Model %in% c("MPI", "EARTH", "HADGEM", "CNRM", "IPSL"), "thick", "thin")
    arrow_data$arrow_type <- ifelse(arrow_data$Scenario == "RCMs", "thick", arrow_data$arrow_type)
    arrow_data$line_size <- ifelse(arrow_data$arrow_type == "thick", 0.6, 0.2)
    
    # Make a subset
    arrow_data_thick <- arrow_data[arrow_data$arrow_type == "thick",]
    arrow_data_thin <- arrow_data[arrow_data$arrow_type != "thick",]
    
    # Compute the shifts in ensembles
    Data_ens <- data.frame(aggregate(. ~ Group, FUN = mean, data = Data), SEL = "all")
    
    # Only for driving model in CMIP5 and CMIP6
    Data_driv <- Data
    Data_driv$Model <-  c(names(PET_over_P_RCMs_1981_2005_a), names(PET_over_P_CMIP5_1981_2005_a), names(PET_over_P_CMIP6_1981_2005_a))
    Data_driv <- Data_driv[which(Data_driv$Model %in% c("MPI", "EARTH", "HADGEM", "CNRM", "IPSL")), ]
    Data_ens <- rbind(Data_ens, data.frame(aggregate(. ~ Group, FUN = mean, data = Data_driv[,1:3]), SEL = "drive"))
    
    arrow_data_ens <- data.frame(
      x = Data_ens$PET_over_P[c(1,3,5,7,9)],
      y = Data_ens$ET_over_P[c(1,3,5,7,9)],
      xend = Data_ens$PET_over_P[c(2,4,6,8,10)],
      yend =Data_ens$ET_over_P[c(2,4,6,8,10)],
      Scenario = c("CMIP5", "CMIP6", "RCMs", "CMIP5", "CMIP6"),
      Color = c(COL_CMIP5, COL_CMIP6, COL_RCMs, COL_CMIP5, COL_CMIP6)
    )
    
    x_offs <- (X_range[2]-X_range[1])*0.04
    y_offs <- (Y_range[2]-Y_range[1])*0.04
    
    X_range <- X_range+c(-x_offs, x_offs)
    Y_range <- Y_range+c(-y_offs, y_offs)
    
    txt_CMIP5_1981_2005 <- bquote(omega == .(sprintf("%.2f",n_CMIP5_1981_2005_a+0.72)))
    txt_CMIP6_1981_2005 <- bquote(omega == .(sprintf("%.2f",n_CMIP6_1981_2005_a+0.72)))
    txt_RCMs_1981_2005 <- bquote(omega == .(sprintf("%.2f",n_RCMs_1981_2005_a+0.72)))
    txt_CMIP5_delta_omega <- bquote(delta*omega == .(sprintf("%.2f",n_CMIP5_2076_2100_a-n_CMIP5_1981_2005_a)))
    txt_CMIP6_delta_omega <- bquote(delta*omega == .(sprintf("%.2f",n_CMIP6_2076_2100_a-n_CMIP6_1981_2005_a)))
    txt_RCMs_delta_omega <- bquote(delta*omega == .(sprintf("%.2f",n_RCMs_2076_2100_a-n_RCMs_1981_2005_a)))
    
    
    # Create data for the continuous line
    continuous_line_data <- data.frame(
      x = c(min(X_range), 1, max(X_range)),  # x coordinates: 0 to 1, and then to the max of X_range
      y = c(min(Y_range), 1, 1)              # y coordinates: 0 to 1 (1:1 line), and then stays at 1 (horizontal line)
    )
    
    # Now create the ggplot
    p <- ggplot() +
      
      # Add continuous line from (0,0) to (1,1) and then horizontally at y = 1
      geom_line(data = continuous_line_data, aes(x = x, y = y), color = "gray", linetype = "solid", size = 0.1) +
      
      # Add continuous lines from original data
      geom_line(data = lines_data, aes(x = X, y = Y, group = Scenario, color = Scenario, linetype = LineType)) +
      
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
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.75, Y_range), label = txt_CMIP5_1981_2005, hjust = 0, color = COL_CMIP5, fontface = "bold") +
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.69, Y_range), label = txt_CMIP6_1981_2005, hjust = 0, color = COL_CMIP6, fontface = "bold") +
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.63, Y_range), label = txt_RCMs_1981_2005, hjust = 0, color = COL_RCMs, fontface = "bold") +
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.57, Y_range), label = txt_CMIP5_delta_omega, hjust = 0, color = COL_CMIP5, fontface = "bold") +
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.51, Y_range), label = txt_CMIP6_delta_omega, hjust = 0, color = COL_CMIP6, fontface = "bold") +
      annotate("text", x = rescale(0.04, X_range), y = rescale(0.45, Y_range), label = txt_RCMs_delta_omega, hjust = 0, color = COL_RCMs, fontface = "bold")
    
    # Inset plot
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
    
    # Saving the plot
    ggsave(plot_name, plot = combined_plot, width = Pl_width, height = Pl_height, dpi = RES, units = "mm")
  }
  
  ##############################################################################
  out <- list(omega_CMIP5_1981_2005=n_CMIP5_1981_2005_a+0.72,
              omega_CMIP6_1981_2005=n_CMIP6_1981_2005_a+0.72,
              omega_RCMs_1981_2005=n_RCMs_1981_2005_a+0.72,
              omega_CMIP5_2076_2100=n_CMIP5_2076_2100_a+0.72,
              omega_CMIP6_2076_2100=n_CMIP6_2076_2100_a+0.72,
              omega_RCMs_2076_2100=n_RCMs_2076_2100_a+0.72,
              out_plot = if(exists("combined_plot")){
                combined_plot
              }else{
                NA
              } 
  )
  
  return(out)
}
