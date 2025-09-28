# Define a custom method function for robust linear regression
robust_rlm <- function(formula, data, ...) {
  rlm(formula, data = data, ...)
}

LM_equation <- function(slope, intercept = NULL, digits = 2, zero_tol = 1e-8) {
  s0 <- is.na(slope) || is.null(slope) || abs(slope) < zero_tol
  i0 <- is.null(intercept) || is.na(intercept) || abs(intercept) < zero_tol
  
  s <- formatC(slope, format = "f", digits = digits)
  # If slope ~ 0, just show the intercept (or 0)
  if (s0) {
    if (i0) return("y = 0")
    i <- formatC(intercept, format = "f", digits = digits)
    return(paste0("y = ", i))
  }
  
  if (i0) {
    # no intercept term shown
    return(paste0("y = ", s, "x"))
  } else {
    i_abs <- formatC(abs(intercept), format = "f", digits = digits)
    sign  <- if (intercept >= 0) " + " else " - "
    return(paste0("y = ", s, "x", sign, i_abs))
  }
}

# Define custom function for scatter plots
make_scatter_plot <- function(data,
                              FIT = TRUE, xy_round = 0.05, xy_offset = 0.04, X_range_man = FALSE, Y_range_man = FALSE,
                              x_lab = "x_lab", y_lab = "y_lab", one_to_one_line = FALSE, one_to_one_line_lab = FALSE,
                              vline = TRUE, hline = TRUE, robust_regression = FALSE, force_origin = FALSE,
                              LM_eq_labels = FALSE, plot_labels = FALSE,
                              plot_path = "./", plot_name = FALSE, save_ggplot2_obj_as = FALSE) {
  
  # Get the model variants
  models <- unique(data$model)
  
  # Perform regression if FIT is TRUE
  if(FIT){
    
    if(robust_regression == TRUE){
      lin_reg <- rlm
    }else{
      lin_reg <- lm
    }
    
    # List to store fits
    fits <- list()
    
    # Iterate over each dataset in data.frame1
    for (i in seq_along(models)) {
      XY <- subset(data, model == models[i], select = c(x, y))
      
      # Calculate fits
      if(force_origin == TRUE){
        fits[[i]] <- lin_reg(XY$y ~ XY$x -1)
      }else{
        fits[[i]] <- lin_reg(XY$y ~ XY$x)
      }
      
    }
  }
  
  # Set x and y range based on all provided data
  XY_range(x = data$x, y = data$y, round = xy_round, offset = xy_offset)
  
  # Overwrite ranges if manual ranges are provided
  if (!identical(X_range_man, FALSE)) X_range <- X_range_man
  if (!identical(Y_range_man, FALSE)) Y_range <- Y_range_man
  
  # Base ggplot setup
  p <- ggplot() +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # Add base lines if enabled
  if (hline) p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey70")
  if (vline) p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")
  if (one_to_one_line) p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#2b2b2b")
  
  # Add fit lines for each dataset in `models` (so they appear behind points)
  if (FIT) {
    for (i in seq_along(models)) {
      
      # Subset the data for the current model
      df <- subset(data, model == models[i], select = c(x, y, color, linetype))
      
      # Add a regression line for the current subset
      p <- p + geom_smooth(
        method = lin_reg,
        formula = (if (isTRUE(force_origin)) (y ~ x - 1) else (y ~ x)),
        data = df,
        aes(x = x, y = y, color = color, fill = color, linetype = linetype),
        se = TRUE
      )
      
      # Add regression equation label if LM_eq_labels is provided and positions are specified
      if (!isFALSE(LM_eq_labels) && !is.na(LM_eq_labels$x[i]) && !is.na(LM_eq_labels$y[i])) {
        # Extract x and y positions for the label
        x_label_pos <- LM_eq_labels$x[i]
        y_label_pos <- LM_eq_labels$y[i]
        
        # Compute the label text first
        label_text <- if (isTRUE(force_origin)) {
          slope <- unname(coef(fits[[i]])[1])
          LM_equation(slope = slope, intercept = 0)  # prints "y = ax" (no + 0)
        } else {
          b <- coef(fits[[i]])
          LM_equation(slope = unname(b[2]), intercept = unname(b[1]))
        }
        
        # Ensure x and y are numeric before rescaling
        if (is.numeric(x_label_pos) && is.numeric(y_label_pos)) {
          p <- p + annotate("text",
                            x = rescale(x_label_pos, X_range),
                            y = rescale(y_label_pos, Y_range),
                            label = label_text,
                            hjust = 0,
                            color = df$color[1])
        } else {
          warning(paste("Non-numeric values found in LM_eq_labels for dataset", i, "- skipping annotation"))
        }
      }
    }
  }
  
  # Add scatter points after fit lines, so points appear on top
  p <- p + geom_point(data = data, aes(x = x, y = y, fill = fill, color = border, shape = shape), size = 3)
  
  # Apply identity scales only once after all layers have been added
  p <- p +
    scale_color_identity() +
    scale_fill_identity() +
    scale_linetype_identity() +
    scale_shape_identity() 

  # Add plot labels if specified in plot_labels
  if (!isFALSE(plot_labels)) {
    
    for (i in seq_along(plot_labels$model)) {
      if(!is.na(plot_labels$x[i]) && !is.na(plot_labels$y[i]) && !is.na(plot_labels$model[i])){

                # Set defaults for label attributes
        label_hjust <- ifelse(!is.null(plot_labels$hjust[i]), plot_labels$hjust[i], 0)
        label_color <- ifelse(!is.null(plot_labels$color[i]), plot_labels$color[i], "#2b2b2b")
        label_fontface <- ifelse(!is.null(plot_labels$fontface[i]), plot_labels$fontface[i], "bold")
        label_size <- ifelse(!is.null(plot_labels$size[i]), plot_labels$size[i], 4)
        
        p <- p + annotate("text",
                          x = rescale(plot_labels$x[i], X_range),
                          y = rescale(plot_labels$y[i], Y_range),
                          label = plot_labels$model[i],
                          hjust = label_hjust,
                          color = label_color,
                          fontface = label_fontface,
                          size = label_size)
      }
    }
  }
  
  # Adding 1:1 line annotation
  if(one_to_one_line == TRUE && length(one_to_one_line_lab) == 2 &&
     is.numeric(one_to_one_line_lab[1]) && is.numeric(one_to_one_line_lab[2])) {
    p <- p + annotate("text", x = rescale(one_to_one_line_lab[1], X_range), y = rescale(one_to_one_line_lab[2], Y_range),
                      label = "1:1", hjust = 0, color = "grey25")
  }
  
  # Add labels using ggrepel
  p <- p + geom_text_repel(
    data = data,
    aes(x = x, y = y, label = label),
    size = 1.5,
    color = ifelse(data$color == COL_RCMs, darken_col(COL_RCMs, 0.6), data$color)
  )
  
  # Adjust x and y limits
  p <- p + scale_x_continuous(limits = X_range, expand = c(0, 0)) +
    scale_y_continuous(limits = Y_range, expand = c(0, 0))
  
  # Save the plot
  if (plot_name != FALSE) ggsave(paste0(plot_path, "/", plot_name), plot = p, width = Pl_width, height = Pl_height, dpi = RES, units = 'mm')
  
  # Return or save ggplot object if needed
  if (save_ggplot2_obj_as != FALSE) assign(x = save_ggplot2_obj_as, value = p, envir = parent.frame())
}

# Function allowing to put symbols and labels above the lines
# Useful for example when line is fitted after the plot is already crated
put_line_behind <- function(ggplot2_plot){
  # Identify layers that are GeomPoint or GeomTextRepel
  idx_pts_txt <- which(vapply(ggplot2_plot$layers,
                              function(l) inherits(l$geom, c("GeomPoint", "GeomTextRepel", "GeomLabelRepel")),
                              logical(1)))
  
  # Move them to the end so they draw last (on top of lines etc.)
  if (length(idx_pts_txt)) {
    ggplot2_plot$layers <- c(ggplot2_plot$layers[-idx_pts_txt], ggplot2_plot$layers[idx_pts_txt])
  }
  
  ggplot2_plot
}

