# ============================================================
# Jarvis_permutation_importance.R
# (universal: works with jarvis_bundle OR rf_hybrid_bundle)
#
# CHANGE vs original:
#   - Plotting is now "Option A" relative importance:
#       share (%) of total ΔRMSE across permuted features
#   - All CSV outputs + return objects remain unchanged (still ΔRMSE)
# ============================================================

run_permutation_importance <- function(
    bundle_path,
    out_dir  = "../outputs",
    plot_dir = "../Plots/ggplot2",
    perm_B   = 1000L,
    seed_perm = 1990L,
    q_lo = 0.25,
    q_hi = 0.75,
    # controls
    vars_perm_order = c("VPD", "P", "Ta", "Rg", "CO2_term"),
    do_vpd_channels = TRUE
) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ranger)
  })
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -----------------------
  # Load bundle from .RData into isolated env, detect which object is present
  # -----------------------
  e0 <- new.env(parent = emptyenv())
  load(bundle_path, envir = e0)
  
  bundle <-
    if (exists("jarvis_bundle", envir = e0, inherits = FALSE)) get("jarvis_bundle", envir = e0) else
      if (exists("rf_hybrid_bundle", envir = e0, inherits = FALSE)) get("rf_hybrid_bundle", envir = e0) else
        stop("No recognized bundle object found in: ", bundle_path,
             "\nExpected: jarvis_bundle or rf_hybrid_bundle")
  
  stopifnot(is.list(bundle))
  stopifnot(is.function(bundle$pred_fun))
  stopifnot(!is.null(bundle$df_opt))
  
  model_id <- bundle$model_id %||% tools::file_path_sans_ext(basename(bundle_path))
  df_opt   <- bundle$df_opt
  eps      <- bundle$eps %||% 1e-6
  
  message("Loaded bundle: ", model_id)
  message("Rows in df_opt: ", nrow(df_opt))
  
  # -----------------------
  # Basic column checks (fail early if data not as expected)
  # -----------------------
  required_cols <- c("ET", "g_eff")
  miss_req <- setdiff(required_cols, names(df_opt))
  if (length(miss_req) > 0) {
    stop("df_opt is missing required columns: ", paste(miss_req, collapse = ", "))
  }
  
  # -----------------------
  # Unified predictor wrapper: pred_fun(newdata) -> list(ET=..., geff=...)
  # -----------------------
  pred_fun <- function(newdata) {
    out <- bundle$pred_fun(newdata)
    
    ET   <- out$ET   %||% out$ET_pred   %||% out$ET_hat
    geff <- out$geff %||% out$geff_pred %||% out$g_eff %||% out$g_eff_pred
    
    if (is.null(ET) || is.null(geff)) {
      stop("bundle$pred_fun(newdata) must return ET and geff (or compatible names).")
    }
    if (length(ET) != nrow(newdata) || length(geff) != nrow(newdata)) {
      stop("Prediction lengths do not match nrow(newdata).")
    }
    
    list(ET = as.numeric(ET), geff = as.numeric(geff))
  }
  
  # -----------------------
  # Select variables to permute
  # - If Jarvis switches exist: respect them
  # - Else if rf_predictors exists: use those + VPD
  # - Else: fallback to fixed order (only those that exist)
  # -----------------------
  vars_perm <- character(0)
  
  has_switches <- any(c("use_VPD","use_P","use_Ta","use_Rg","use_CO2") %in% names(bundle))
  if (isTRUE(has_switches)) {
    if (isTRUE(bundle$use_VPD)) vars_perm <- c(vars_perm, "VPD")
    if (isTRUE(bundle$use_P))   vars_perm <- c(vars_perm, "P")
    if (isTRUE(bundle$use_Ta))  vars_perm <- c(vars_perm, "Ta")
    if (isTRUE(bundle$use_Rg))  vars_perm <- c(vars_perm, "Rg")
    if (isTRUE(bundle$use_CO2)) vars_perm <- c(vars_perm, "CO2_term")
  } else if (!is.null(bundle$rf_predictors)) {
    vars_perm <- unique(c("VPD", bundle$rf_predictors))
  } else {
    vars_perm <- vars_perm_order
  }
  
  vars_perm <- intersect(vars_perm_order, intersect(vars_perm, names(df_opt)))
  if (length(vars_perm) < 1) stop("No permutation variables found in df_opt.")
  
  # -----------------------
  # Label map for subscripts (parse-able expressions)
  # -----------------------
  label_map <- c(
    VPD      = "VPD",
    P        = "P",
    Ta       = "T[a]",
    Rg       = "R[g]",
    CO2_term = "CO[2]"
  )
  label_expr <- function(v) {
    x <- unname(label_map[v])          # NOTE: single-bracket lookup is safe
    if (is.na(x) || length(x) == 0) v else x
  }
  
  # ============================================================
  # Permutation importance with draws (UNCHANGED)
  # ============================================================
  perm_importance_model_draws <- function(data, vars, pred_fun, B = 200L, seed = 2024,
                                          q_lo = 0.25, q_hi = 0.75) {
    set.seed(seed)
    
    pred0 <- pred_fun(data)
    base_rmse_ET   <- sqrt(mean((data$ET    - pred0$ET  )^2, na.rm = TRUE))
    base_rmse_geff <- sqrt(mean((data$g_eff - pred0$geff)^2, na.rm = TRUE))
    
    draws <- bind_rows(lapply(vars, function(v) {
      dET   <- numeric(B)
      dGeff <- numeric(B)
      
      for (b in seq_len(B)) {
        dfp <- data
        dfp[[v]] <- sample(dfp[[v]])
        
        # prevent stale derived cols if present
        if ("inv_sqrt_VPD" %in% names(dfp)) dfp$inv_sqrt_VPD <- NULL
        if ("K_ET"         %in% names(dfp)) dfp$K_ET         <- NULL
        
        predp <- pred_fun(dfp)
        
        rmse_ET_p   <- sqrt(mean((data$ET    - predp$ET  )^2, na.rm = TRUE))
        rmse_geff_p <- sqrt(mean((data$g_eff - predp$geff)^2, na.rm = TRUE))
        
        dET[b]   <- rmse_ET_p   - base_rmse_ET
        dGeff[b] <- rmse_geff_p - base_rmse_geff
      }
      
      tibble(feature = v, b = seq_len(B), dET = dET, dGeff = dGeff)
    }))
    
    summary <- draws %>%
      group_by(feature) %>%
      summarise(
        med_delta_ET   = median(dET,   na.rm = TRUE),
        lo_ET          = quantile(dET,   q_lo, na.rm = TRUE),
        hi_ET          = quantile(dET,   q_hi, na.rm = TRUE),
        
        med_delta_geff = median(dGeff, na.rm = TRUE),
        lo_geff        = quantile(dGeff, q_lo, na.rm = TRUE),
        hi_geff        = quantile(dGeff, q_hi, na.rm = TRUE),
        
        B = dplyr::n(),
        .groups = "drop"
      )
    
    list(draws = draws, summary = summary)
  }
  
  perm_res <- perm_importance_model_draws(
    data     = df_opt,
    vars     = vars_perm,
    pred_fun = pred_fun,
    B        = perm_B,
    seed     = seed_perm,
    q_lo     = q_lo,
    q_hi     = q_hi
  )
  
  perm_draws   <- perm_res$draws
  perm_summary <- perm_res$summary
  
  print(perm_summary)
  
  # CSV outputs UNCHANGED
  write_csv(perm_summary, file.path(out_dir, paste0("perm_", model_id, "_fit_summary.csv")))
  write_csv(perm_draws,   file.path(out_dir, paste0("perm_", model_id, "_fit_draws.csv")))
  
  # ============================================================
  # VPD pathway decomposition (ET) WITH DRAWS (optional) (UNCHANGED)
  # ============================================================
  recompute_K_ET <- bundle$recompute_K_ET %||% NULL
  can_vpd_channels <- isTRUE(do_vpd_channels) &&
    ("VPD" %in% vars_perm) &&
    is.function(recompute_K_ET)
  
  perm_vpd_channels_draws <- NULL
  perm_vpd_channels_summary <- NULL
  plot_ET_combined <- NULL
  
  if (isTRUE(can_vpd_channels)) {
    
    perm_VPD_channels_draws <- function(data, pred_fun, recompute_K_ET, B = 200L, seed = 2024,
                                        q_lo = 0.25, q_hi = 0.75) {
      set.seed(seed)
      
      pred0 <- pred_fun(data)
      base_rmse_ET <- sqrt(mean((data$ET - pred0$ET)^2, na.rm = TRUE))
      
      K0    <- recompute_K_ET(data)
      geff0 <- pred0$geff
      
      d_stom   <- numeric(B)
      d_demand <- numeric(B)
      
      for (b in seq_len(B)) {
        dfp <- data
        dfp$VPD <- sample(dfp$VPD)
        
        if ("inv_sqrt_VPD" %in% names(dfp)) dfp$inv_sqrt_VPD <- NULL
        if ("K_ET"         %in% names(dfp)) dfp$K_ET         <- NULL
        
        # STOMATAL only: permute VPD -> recompute g, keep K fixed
        geff_p  <- pred_fun(dfp)$geff
        ET_stom <- K0 * geff_p
        rmse_stom <- sqrt(mean((data$ET - ET_stom)^2, na.rm = TRUE))
        d_stom[b] <- rmse_stom - base_rmse_ET
        
        # DEMAND only: permute VPD -> recompute K, keep g fixed
        K_p    <- recompute_K_ET(dfp)
        ET_dem <- K_p * geff0
        rmse_dem <- sqrt(mean((data$ET - ET_dem)^2, na.rm = TRUE))
        d_demand[b] <- rmse_dem - base_rmse_ET
      }
      
      draws <- tibble(
        b = seq_len(B),
        dET_stom   = d_stom,
        dET_demand = d_demand
      ) %>%
        pivot_longer(-b, names_to = "channel", values_to = "dET") %>%
        mutate(channel = recode(
          channel,
          dET_stom   = "VPD[stomatal]",
          dET_demand = "VPD[demand]"
        ))
      
      summary <- draws %>%
        group_by(channel) %>%
        summarise(
          med_delta_ET = median(dET, na.rm = TRUE),
          lo_ET        = quantile(dET, q_lo, na.rm = TRUE),
          hi_ET        = quantile(dET, q_hi, na.rm = TRUE),
          B = dplyr::n(),
          .groups = "drop"
        )
      
      list(draws = draws, summary = summary)
    }
    
    vpd_res <- perm_VPD_channels_draws(
      data           = df_opt,
      pred_fun       = pred_fun,
      recompute_K_ET = recompute_K_ET,
      B             = perm_B,
      seed          = seed_perm,
      q_lo          = q_lo,
      q_hi          = q_hi
    )
    
    perm_vpd_channels_draws   <- vpd_res$draws
    perm_vpd_channels_summary <- vpd_res$summary
    
    print(perm_vpd_channels_summary)
    
    # CSV outputs UNCHANGED
    write_csv(perm_vpd_channels_summary,
              file.path(out_dir, paste0("perm_", model_id, "_VPD_channels_ET_summary.csv")))
    write_csv(perm_vpd_channels_draws,
              file.path(out_dir, paste0("perm_", model_id, "_VPD_channels_ET_draws.csv")))
  } else {
    message("Skipping VPD channel analysis (need VPD + bundle$recompute_K_ET).")
  }
  
  # ============================================================
  # NEW: helper to compute "share of total ΔRMSE" from draws
  #   - uses pmax(delta, 0) to avoid negative shares
  #   - returns median + quantile whiskers in % share
  # ============================================================
  .share_summary_from_draws <- function(draws_long, value_col, group_col,
                                        q_lo = 0.25, q_hi = 0.75) {
    value_col <- rlang::ensym(value_col)
    group_col <- rlang::ensym(group_col)
    
    d <- draws_long %>%
      mutate(.delta_pos = pmax(!!value_col, 0)) %>%
      group_by(b) %>%
      mutate(.den = sum(.delta_pos, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(.share = if_else(is.finite(.den) & .den > 0, 100 * .delta_pos / .den, NA_real_))
    
    out <- d %>%
      group_by(!!group_col) %>%
      summarise(
        med_share = median(.share, na.rm = TRUE),
        lo_share  = quantile(.share, q_lo, na.rm = TRUE),
        hi_share  = quantile(.share, q_hi, na.rm = TRUE),
        .groups = "drop"
      )
    
    out
  }
  
  # ============================================================
  # Plots (3) — NOW RELATIVE (share % of total ΔRMSE)
  #   File names remain UNCHANGED
  # ============================================================
  subtitle_txt <- sprintf("Bars: median; whiskers: %.0f–%.0f%% quantiles (B=%d)",
                          100*q_lo, 100*q_hi, perm_B)
  
  # ---- (1) g_eff plot (share %) ----
  geff_share <- perm_draws %>%
    filter(feature %in% vars_perm_order) %>%
    .share_summary_from_draws(value_col = dGeff, group_col = feature, q_lo = q_lo, q_hi = q_hi) %>%
    mutate(label = map_chr(as.character(feature), label_expr))
  
  p_perm_geff <- ggplot(geff_share,
                        aes(x = reorder(label, med_share), y = med_share)) +
    geom_col() +
    geom_errorbar(aes(ymin = lo_share, ymax = hi_share), width = 0.2) +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = NULL,
      y = "Relative importance (% of total ΔRMSE from permutation)",
      title = bquote("Permutation importance — " * g[eff] * " (" * .(model_id) * ")"),
      subtitle = subtitle_txt
    )
  
  # ---- (2) ET plot (share %) ----
  et_share <- perm_draws %>%
    filter(feature %in% vars_perm_order) %>%
    .share_summary_from_draws(value_col = dET, group_col = feature, q_lo = q_lo, q_hi = q_hi) %>%
    mutate(label = map_chr(as.character(feature), label_expr))
  
  p_perm_ET <- ggplot(et_share,
                      aes(x = reorder(label, med_share), y = med_share)) +
    geom_col() +
    geom_errorbar(aes(ymin = lo_share, ymax = hi_share), width = 0.2) +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = NULL,
      y = "Relative importance (% of total ΔRMSE from permutation)",
      title = paste0("Permutation importance — ET (", model_id, ")"),
      subtitle = subtitle_txt
    )
  
  print(p_perm_geff)
  print(p_perm_ET)
  
  ggsave(file.path(plot_dir, paste0("perm_", model_id, "_geff_median_quantiles.png")),
         p_perm_geff, width = 7, height = 4.8, dpi = 300)
  
  ggsave(file.path(plot_dir, paste0("perm_", model_id, "_ET_median_quantiles.png")),
         p_perm_ET, width = 7, height = 4.8, dpi = 300)
  
  # ---- (3) ET plot with VPD channels (share %) ----
  p_perm_ET_channels <- NULL
  if (isTRUE(can_vpd_channels)) {
    
    # Build combined DRAWS table for shares:
    #  - P, Ta, Rg, CO2_term from perm_draws (dET)
    #  - VPD stom/demand from perm_vpd_channels_draws (dET)
    nonvpd_draws <- perm_draws %>%
      filter(feature %in% c("P", "Ta", "Rg", "CO2_term")) %>%
      transmute(b = b, label = feature, dET = dET)
    
    vpd_draws <- perm_vpd_channels_draws %>%
      transmute(b = b, label = channel, dET = dET)
    
    combined_draws <- bind_rows(nonvpd_draws, vpd_draws) %>%
      mutate(
        label = as.character(label),
        label = if_else(label %in% names(label_map), unname(label_map[label]), label)
      )
    
    combined_share <- .share_summary_from_draws(
      draws_long = combined_draws,
      value_col  = dET,
      group_col  = label,
      q_lo = q_lo, q_hi = q_hi
    )
    
    plot_ET_combined <- combined_share %>%
      transmute(
        label,
        med_delta_ET = med_share,
        lo_ET = lo_share,
        hi_ET = hi_share
      ) %>%
      arrange(desc(med_delta_ET))
    
    p_perm_ET_channels <- ggplot(plot_ET_combined,
                                 aes(x = reorder(label, med_delta_ET), y = med_delta_ET)) +
      geom_col() +
      geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
      coord_flip() +
      theme_bw() +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      labs(
        x = NULL,
        y = "Relative importance (% of total ΔRMSE from permutation)",
        title = paste0("Permutation importance — ET (", model_id, ") with VPD pathways"),
        subtitle = subtitle_txt
      )
    
    print(p_perm_ET_channels)
    
    ggsave(file.path(plot_dir, paste0("perm_", model_id, "_ET_VPDchannels_median_quantiles.png")),
           p_perm_ET_channels, width = 7, height = 5.2, dpi = 300)
    
    # CSV output name UNCHANGED; content now matches the plotted quantity (shares).
    write_csv(plot_ET_combined,
              file.path(out_dir, paste0("perm_", model_id, "_ET_combined_median_quantiles.csv")))
  }
  
  invisible(list(
    model_id = model_id,
    bundle_path = bundle_path,
    vars_perm = vars_perm,
    perm_summary = perm_summary,
    perm_draws = perm_draws,
    vpd_channels_summary = perm_vpd_channels_summary,
    vpd_channels_draws = perm_vpd_channels_draws,
    plots = list(geff = p_perm_geff, ET = p_perm_ET, ET_channels = p_perm_ET_channels),
    plot_ET_combined = plot_ET_combined
  ))
}
