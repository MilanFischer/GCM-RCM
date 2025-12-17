# Permutation insets
make_perm_plot <- function(csv_path,
                            response = c("geff","ET"),
                            order = c("VPD","P","Ta","Rg","CO2_term"),
                            q_lo = 0.25, q_hi = 0.75) {
  
  response <- match.arg(response)
  
  df <- readr::read_csv(csv_path, show_col_types = FALSE)
  
  # -----------------------
  # Pick delta columns
  # -----------------------
  if (response == "geff") {
    med <- df$med_delta_geff
    lo  <- df$lo_geff
    hi  <- df$hi_geff
  } else {
    med <- df$med_delta_ET
    lo  <- df$lo_ET
    hi  <- df$hi_ET
  }
  
  # -----------------------
  # Option A normalization (share of total ΔRMSE)
  # -----------------------
  den <- sum(pmax(med, 0), na.rm = TRUE)
  
  df2 <- df |>
    dplyr::mutate(
      feature = factor(feature, levels = order),
      share    = 100 * pmax(med, 0) / den,
      share_lo = 100 * pmax(lo,  0) / den,
      share_hi = 100 * pmax(hi,  0) / den
    ) |>
    dplyr::filter(!is.na(feature)) |>
    dplyr::arrange(dplyr::desc(share))
  
  # -----------------------
  # Expression-style labels (IMPORTANT PART)
  # -----------------------
  label_map <- c(
    VPD      = "VPD",
    P        = "P",
    Ta       = "T[a]",
    Rg       = "R[g]",
    CO2_term = "CO[2]"
  )
  
  df2$lab <- unname(label_map[as.character(df2$feature)])
  
  # -----------------------
  # Plot
  # -----------------------
  ggplot2::ggplot(df2, ggplot2::aes(x = reorder(lab, share), y = share)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = share_lo, ymax = share_hi),
      width = 0.15, linewidth = 0.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) +  # <-- enables subscripts
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7, colour = "#222222"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(2,2,2,2)
    ) +
    ggplot2::labs(y = NULL, x = NULL)
}

make_perm_plot_ET_vpd_channels <- function(summary_csv,
                                            vpd_channels_csv,
                                            order = c("VPD[demand]","VPD[stomatal]","P","Ta","Rg","CO2_term"),
                                            q_lo = 0.25, q_hi = 0.75) {
  
  # -----------------------
  # Read CSVs
  # -----------------------
  s <- readr::read_csv(summary_csv, show_col_types = FALSE)
  v <- readr::read_csv(vpd_channels_csv, show_col_types = FALSE)
  
  # -----------------------
  # Non-VPD features
  # -----------------------
  nonvpd <- s |>
    dplyr::filter(feature %in% c("P","Ta","Rg","CO2_term")) |>
    dplyr::transmute(
      key = feature,
      med = med_delta_ET,
      lo  = lo_ET,
      hi  = hi_ET
    )
  
  # -----------------------
  # VPD channels
  # -----------------------
  vpdch <- v |>
    dplyr::transmute(
      key = channel,          # "VPD[demand]" / "VPD[stomatal]"
      med = med_delta_ET,
      lo  = lo_ET,
      hi  = hi_ET
    )
  
  # -----------------------
  # Combine
  # -----------------------
  df <- dplyr::bind_rows(vpdch, nonvpd) |>
    dplyr::mutate(
      key = factor(as.character(key), levels = order)
    ) |>
    dplyr::filter(!is.na(key))
  
  # -----------------------
  # Option A normalization (share of total ΔRMSE)
  # -----------------------
  den <- sum(pmax(df$med, 0), na.rm = TRUE)
  
  df <- df |>
    dplyr::mutate(
      share    = 100 * pmax(med, 0) / den,
      share_lo = 100 * pmax(lo,  0) / den,
      share_hi = 100 * pmax(hi,  0) / den
    ) |>
    dplyr::arrange(dplyr::desc(share))
  
  # -----------------------
  # Expression-style labels (THIS IS THE IMPORTANT PART)
  # -----------------------
  label_map <- c(
    "VPD[demand]"   = "VPD[demand]",
    "VPD[stomatal]" = "VPD[stomatal]",
    P               = "P",
    Ta              = "T[a]",
    Rg              = "R[g]",
    CO2_term        = "CO[2]"
  )
  
  df$lab <- unname(label_map[as.character(df$key)])
  
  # -----------------------
  # Plot
  # -----------------------
  ggplot2::ggplot(df, ggplot2::aes(x = reorder(lab, share), y = share)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = share_lo, ymax = share_hi),
      width = 0.15, linewidth = 0.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) +  # <-- parses subscripts
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7, colour = "#222222"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(2,2,2,2)
    )
}
