# ============================================================
# Jarvis&RF_ALE_plots.R (universal)
#   - Option A: compute_ALE = TRUE  -> ALE (smoothed over full observed range)
#   - Option B: compute_ALE = FALSE -> 1-D prediction slices (others held fixed)
# Saves ONLY TWO PLOTS (anchored + observed overlay), WITH model_id in filename:
#   - ALE_ET_smoothed_abs_<model_id>.png   (or SLICES_ET_abs_<model_id>.png)
#   - ALE_geff_smoothed_abs_<model_id>.png (or SLICES_geff_abs_<model_id>.png)
# ============================================================

run_ALE_plots <- function(
    bundle_path,
    compute_ALE = TRUE,
    vars = c("P","Rg","Ta","VPD"),
    # Option B settings (slices)
    hold_others = c("median","mean"),
    slice_grid_points = 300,
    # Option A settings (ALE)
    n_bins = 50,
    smooth_span = 0.5,
    smooth_grid_points = 300,
    # overlays
    show_rug = TRUE,
    overlay_on_relative = FALSE,   # optional (not saved by default)
    obs_max_points = 5000,
    # output
    plot_dir = "../Plots/ggplot2",
    width = 9,
    height = 6.5,
    dpi = 300
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(purrr)
  })
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  hold_others <- match.arg(hold_others)
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -----------------------
  # Load bundle (isolated env) and detect object
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
  
  # ensure response columns exist for overlays
  req_cols <- c("ET","g_eff")
  miss <- setdiff(req_cols, names(df_opt))
  if (length(miss)) stop("df_opt missing required columns: ", paste(miss, collapse = ", "))
  
  vars <- intersect(vars, names(df_opt))
  if (length(vars) < 1) stop("None of requested vars exist in df_opt.")
  
  # -----------------------
  # Unified predictor wrapper: pred_fun(newdata)-> list(ET, geff)
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
  # Shared utilities
  # -----------------------
  var_labels <- c(
    P = "P (precipitation)",
    Rg = "Rg (radiation)",
    Ta = "Ta (°C)",
    VPD = "VPD (kPa)",
    CO2_term = "CO₂ term"
  )
  
  make_obs_pts <- function(y_name) {
    out <- bind_rows(lapply(vars, function(v) {
      tibble(var = v, x = df_opt[[v]], y = df_opt[[y_name]])
    })) %>% filter(is.finite(x), is.finite(y))
    if (nrow(out) > obs_max_points) out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
    out
  }
  
  obs_pts_et   <- make_obs_pts("ET")
  obs_pts_geff <- make_obs_pts("g_eff")
  
  # Helper: sanitize model_id for filenames
  safe_id <- gsub("[^A-Za-z0-9_\\-]+", "_", model_id)
  
  # -----------------------
  # OPTION A: ALE
  # -----------------------
  if (isTRUE(compute_ALE)) {
    
    ale_one <- function(var, data, n_bins) {
      z <- data[[var]]
      zf <- z[is.finite(z)]
      edges <- as.numeric(quantile(zf, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
      edges[1] <- min(zf, na.rm = TRUE)
      edges[length(edges)] <- max(zf, na.rm = TRUE)
      
      dET <- numeric(n_bins)
      dG  <- numeric(n_bins)
      cnt <- integer(n_bins)
      
      for (j in seq_len(n_bins)) {
        lo <- edges[j]; hi <- edges[j + 1]
        idx <- which(z >= lo & (z < hi | (j == n_bins & z <= hi)))
        cnt[j] <- length(idx)
        if (cnt[j] == 0) { dET[j] <- 0; dG[j] <- 0; next }
        
        df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
        df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
        df_both <- rbind(df_lo, df_hi)
        
        pred_both <- pred_fun(df_both)
        m <- cnt[j]
        
        et_lo <- pred_both$ET[seq_len(m)]
        et_hi <- pred_both$ET[(m + 1L):(2L * m)]
        g_lo  <- pred_both$geff[seq_len(m)]
        g_hi  <- pred_both$geff[(m + 1L):(2L * m)]
        
        dET[j] <- mean(et_hi - et_lo, na.rm = TRUE)
        dG[j]  <- mean(g_hi  - g_lo,  na.rm = TRUE)
      }
      
      ale_ET <- cumsum(dET)
      ale_G  <- cumsum(dG)
      mids   <- (edges[-1] + edges[-length(edges)]) / 2
      
      N <- sum(cnt)
      if (N > 0) {
        ale_ET <- ale_ET - sum(ale_ET * cnt) / N
        ale_G  <- ale_G  - sum(ale_G  * cnt) / N
      }
      
      tibble(var = var, x = mids, ALE_ET = ale_ET, ALE_geff = ale_G, n_in_bin = cnt)
    }
    
    ale_all <- bind_rows(lapply(vars, ale_one, data = df_opt, n_bins = n_bins))
    
    # anchor to mean predictions
    pred0 <- pred_fun(df_opt)
    mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
    mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)
    
    .fit_loess_safe <- function(y, x, w, span) {
      try(loess(y ~ x, weights = w, span = span,
                control = loess.control(surface = "direct")), silent = TRUE)
    }
    .predict_at <- function(fit, x_train, y_train, x_new) {
      if (!inherits(fit, "try-error")) {
        as.numeric(predict(fit, data.frame(x = x_new)))
      } else {
        as.numeric(approx(x = x_train, y = y_train, xout = x_new, rule = 2, ties = mean)$y)
      }
    }
    
    smooth_one_var <- function(vname, ale_mid_df, span, grid_n) {
      dfv <- ale_mid_df %>% filter(var == vname)
      x_min <- min(df_opt[[vname]], na.rm = TRUE)
      x_max <- max(df_opt[[vname]], na.rm = TRUE)
      x_grid <- seq(x_min, x_max, length.out = grid_n)
      
      fit_ET <- .fit_loess_safe(dfv$ALE_ET,   dfv$x, dfv$n_in_bin, span)
      fit_G  <- .fit_loess_safe(dfv$ALE_geff, dfv$x, dfv$n_in_bin, span)
      
      y_ET_grid <- .predict_at(fit_ET, dfv$x, dfv$ALE_ET,   x_grid)
      y_G_grid  <- .predict_at(fit_G,  dfv$x, dfv$ALE_geff, x_grid)
      
      y_ET_mid <- .predict_at(fit_ET, dfv$x, dfv$ALE_ET,   dfv$x)
      y_G_mid  <- .predict_at(fit_G,  dfv$x, dfv$ALE_geff, dfv$x)
      
      off_ET <- weighted.mean(y_ET_mid, w = dfv$n_in_bin, na.rm = TRUE)
      off_G  <- weighted.mean(y_G_mid,  w = dfv$n_in_bin, na.rm = TRUE)
      
      tibble(
        var = vname,
        x = x_grid,
        ALE_ET_abs_s   = (y_ET_grid - off_ET) + mu_ET_pred,
        ALE_geff_abs_s = (y_G_grid  - off_G ) + mu_geff_pred,
        ALE_ET_s       = (y_ET_grid - off_ET),
        ALE_geff_s     = (y_G_grid  - off_G )
      )
    }
    
    ale_all_sm <- map_dfr(vars, smooth_one_var,
                          ale_mid_df = ale_all,
                          span = smooth_span,
                          grid_n = smooth_grid_points)
    
    # ---- TWO anchored plots, with model_id in title ----
    p_et_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_ET_abs_s)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
      labs(
        x = NULL, y = "ET (mm/yr) — anchored",
        title = paste0("ALE (ET) — smoothed, absolute (anchored) | ", model_id)
      ) +
      theme_bw() +
      geom_point(data = obs_pts_et, aes(x = x, y = y),
                 inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    
    p_geff_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_geff_abs_s)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
      labs(
        x = NULL, y = "g_eff (mm s^-1) — anchored",
        title = paste0("ALE (g_eff) — smoothed, absolute (anchored) | ", model_id)
      ) +
      theme_bw() +
      geom_point(data = obs_pts_geff, aes(x = x, y = y),
                 inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    
    print(p_et_abs); print(p_geff_abs)
    
    # ---- filenames WITH model_id ----
    f1 <- file.path(plot_dir, paste0("ALE_ET_smoothed_abs_", safe_id, ".png"))
    f2 <- file.path(plot_dir, paste0("ALE_geff_smoothed_abs_", safe_id, ".png"))
    
    ggsave(f1, p_et_abs,   width = width, height = height, dpi = dpi)
    ggsave(f2, p_geff_abs, width = width, height = height, dpi = dpi)
    
    invisible(list(
      model_id = model_id,
      option = "ALE",
      plots = list(ET_abs = p_et_abs, geff_abs = p_geff_abs),
      files = c(ET_abs = f1, geff_abs = f2),
      ale_raw = ale_all,
      ale_smoothed = ale_all_sm
    ))
    
  } else {
    
    # -----------------------
    # OPTION B: 1-D SLICES
    # -----------------------
    .mode <- function(x) {
      x <- x[!is.na(x)]
      if (!length(x)) return(NA)
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    ref_row <- df_opt[1, , drop = FALSE]
    num_cols <- names(df_opt)[vapply(df_opt, is.numeric, logical(1))]
    
    if (hold_others == "mean") {
      ref_row[1, num_cols] <- lapply(df_opt[num_cols], function(col) mean(col, na.rm = TRUE))
    } else {
      ref_row[1, num_cols] <- lapply(df_opt[num_cols], function(col) median(col, na.rm = TRUE))
    }
    
    non_num_cols <- setdiff(names(df_opt), num_cols)
    if (length(non_num_cols)) {
      ref_modes <- lapply(df_opt[non_num_cols], .mode)
      for (nm in non_num_cols) ref_row[[nm]] <- ref_modes[[nm]]
    }
    
    make_slice <- function(v) {
      x_min <- min(df_opt[[v]], na.rm = TRUE)
      x_max <- max(df_opt[[v]], na.rm = TRUE)
      x_grid <- seq(x_min, x_max, length.out = slice_grid_points)
      
      df <- ref_row[rep(1, length(x_grid)), , drop = FALSE]
      df[[v]] <- x_grid
      
      preds <- pred_fun(df)
      
      tibble(
        var = v,
        x = x_grid,
        ET_pred   = preds$ET,
        geff_pred = preds$geff
      )
    }
    
    slices <- map_dfr(vars, make_slice)
    
    p_et_abs <- ggplot(slices, aes(x = x, y = ET_pred)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
      labs(
        x = NULL, y = "ET (mm/yr)",
        title = paste0("1-D slices for ET — absolute predictions | ", model_id)
      ) +
      theme_bw() +
      geom_point(data = obs_pts_et, aes(x = x, y = y),
                 inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    
    p_geff_abs <- ggplot(slices, aes(x = x, y = geff_pred)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
      labs(
        x = NULL, y = "g_eff (mm s^-1)",
        title = paste0("1-D slices for g_eff — absolute predictions | ", model_id)
      ) +
      theme_bw() +
      geom_point(data = obs_pts_geff, aes(x = x, y = y),
                 inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    
    print(p_et_abs); print(p_geff_abs)
    
    # ---- filenames WITH model_id ----
    f1 <- file.path(plot_dir, paste0("SLICES_ET_abs_", safe_id, ".png"))
    f2 <- file.path(plot_dir, paste0("SLICES_geff_abs_", safe_id, ".png"))
    
    ggsave(f1, p_et_abs,   width = width, height = height, dpi = dpi)
    ggsave(f2, p_geff_abs, width = width, height = height, dpi = dpi)
    
    invisible(list(
      model_id = model_id,
      option = "SLICES",
      plots = list(ET_abs = p_et_abs, geff_abs = p_geff_abs),
      files = c(ET_abs = f1, geff_abs = f2),
      slices = slices,
      ale_raw = ale_all,
      ale_smoothed = ale_all_sm
    ))
  }
}
