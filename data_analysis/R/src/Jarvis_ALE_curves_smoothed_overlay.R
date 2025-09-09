# ============================================================
# Jarvis model interpretation:
# - Option A (compute_ALE = TRUE): ALE curves (ET and g_eff), smoothed over full observed range
# - Option B (compute_ALE = FALSE): 1-D prediction slices with other variables held fixed
# In both options, observed points are overlaid on ANCHORED plots.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
})

# -----------------------
# Requirements (must exist)
# -----------------------
req <- c("jarvis_pars","jarvis_df_opt","jarvis_fVPD_mode","rhoAir","CpAir","gamma","jarvis_model_ET_components")
miss <- setdiff(req, ls(envir = .GlobalEnv))
if (length(miss)) stop("Missing in environment: ", paste(miss, collapse = ", "))

# -----------------------
# User settings (tune here)
# -----------------------
compute_ALE       <- TRUE      # TRUE = ALE; FALSE = 1-D slices (others held fixed)
hold_others       <- "median"  # for slices: "median" or "mean" for numeric; non-numeric -> mode
slice_grid_points <- 300       # for slices: points per variable for 1-D lines

if (!exists("eps", inherits = TRUE)) eps <- 1e-6  # numeric guard

# --- ALE binning
n_bins        <- 50    # equal-frequency ALE bins per variable
min_per_bin   <- 200   # (optional diagnostic usage; not enforced below)

# --- LOESS smoothing for ALE
smooth_span        <- 0.5   # LOESS span (↑ smoother, ↓ more detail)
smooth_grid_points <- 300   # x-points per variable for the smoothed ALE line

# --- Plotting & overlays
show_rug            <- TRUE    # show rug under relative plots
overlay_on_relative <- FALSE   # also overlay obs on relative plots (centered)
obs_max_points      <- 5000    # cap for observed scatter points per plot

# -----------------------
# Helpers: derive cols Jarvis needs + predict ET & g_eff
# -----------------------
rederive_cols <- function(df) {
  df |>
    mutate(
      inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
      K_ET = (rhoAir * CpAir / gamma) * VPD * (1/1000) *
        1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
    )
}

predict_ET_geff <- function(df) {
  df2 <- rederive_cols(df)
  comp <- jarvis_model_ET_components(par = jarvis_pars, input_df = df2, fVPD_mode = jarvis_fVPD_mode)
  list(ET = as.numeric(comp$ET_pred), geff = as.numeric(comp$g_eff_pred))
}

# -----------------------
# ALE for one variable (Molnar)
# -----------------------
ale_one <- function(var, data = jarvis_df_opt, n_bins = 20) {
  z <- data[[var]]
  edges <- as.numeric(quantile(z[is.finite(z)], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  edges[1] <- min(z, na.rm = TRUE); edges[length(edges)] <- max(z, na.rm = TRUE)
  
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
    
    pred_both <- predict_ET_geff(df_both)
    
    m <- cnt[j]
    et_lo <- pred_both$ET[seq_len(m)]
    et_hi <- pred_both$ET[(m + 1L):(2L * m)]
    g_lo  <- pred_both$geff[seq_len(m)]
    g_hi  <- pred_both$geff[(m + 1L):(2L * m)]
    
    dET[j] <- mean(et_hi - et_lo, na.rm = TRUE)
    dG[j]  <- mean(g_hi  - g_lo,  na.rm = TRUE)
  }
  
  # accumulate increments and center (standard ALE)
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

# -----------------------
# Utilities shared by both options
# -----------------------
vars <- c("P","Rg","Ta","VPD")
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")

make_obs_pts <- function(y_name) {
  out <- dplyr::bind_rows(lapply(vars, function(v) {
    tibble(var = v, x = jarvis_df_opt[[v]], y = jarvis_df_opt[[y_name]])
  })) %>% dplyr::filter(is.finite(x), is.finite(y))
  if (nrow(out) > obs_max_points) {
    out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
  }
  out
}

# -----------------------
# Compute plots (ALE vs Slices)
# -----------------------
if (isTRUE(compute_ALE)) {
  # ======== OPTION A: ALE (smoothed over full observed range) ========
  ale_list <- lapply(vars, ale_one, data = jarvis_df_opt, n_bins = n_bins)
  ale_all  <- dplyr::bind_rows(ale_list)  # centered, at midpoints
  
  # Mean predictions (anchoring absolute-scale curves)
  pred0 <- predict_ET_geff(jarvis_df_opt)
  mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
  mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)
  
  # Dense-grid LOESS smoothing per variable over full observed range
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
    x_min <- min(jarvis_df_opt[[vname]], na.rm = TRUE)
    x_max <- max(jarvis_df_opt[[vname]], na.rm = TRUE)
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
      ALE_ET_s   = y_ET_grid - off_ET,
      ALE_geff_s = y_G_grid  - off_G,
      ALE_ET_abs_s   = (y_ET_grid - off_ET) + mu_ET_pred,
      ALE_geff_abs_s = (y_G_grid  - off_G)  + mu_geff_pred
    )
  }
  ale_all_sm <- map_dfr(vars, smooth_one_var,
                        ale_mid_df = ale_all,
                        span = smooth_span,
                        grid_n = smooth_grid_points)
  
  # Relative (centered)
  plot_ale_et_rel <- ggplot(ale_all_sm, aes(x = x, y = ALE_ET_s)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "ALE of ET (mm/yr) — centered",
         title = "Accumulated Local Effects (ET) — smoothed over full observed range") +
    theme_bw()
  
  plot_ale_geff_rel <- ggplot(ale_all_sm, aes(x = x, y = ALE_geff_s)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "ALE of g_eff (mm s^-1) — centered",
         title = "Accumulated Local Effects (g_eff) — smoothed over full observed range") +
    theme_bw()
  
  if (isTRUE(show_rug)) {
    plot_ale_et_rel  <- plot_ale_et_rel  + geom_rug(data = ale_all, aes(x = x), sides = "b", inherit.aes = FALSE, alpha = 0.2)
    plot_ale_geff_rel<- plot_ale_geff_rel+ geom_rug(data = ale_all, aes(x = x), sides = "b", inherit.aes = FALSE, alpha = 0.2)
  }
  
  # Anchored (absolute)
  plot_ale_et_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_ET_abs_s)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "ET (mm/yr) — anchored",
         title = "Accumulated Local Effects (ET) — smoothed, absolute (anchored to mean)") +
    theme_bw()
  
  plot_ale_geff_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_geff_abs_s)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "g_eff (mm s^-1) — anchored",
         title = "Accumulated Local Effects (g_eff) — smoothed, absolute (anchored to mean)") +
    theme_bw()
  
  # ---- Overlay observed points (anchored plots) ----
  obs_pts_et   <- make_obs_pts("ET")
  obs_pts_geff <- make_obs_pts("g_eff")
  plot_ale_et_abs   <- plot_ale_et_abs   +
    geom_point(data = obs_pts_et,   aes(x = x, y = y), inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  plot_ale_geff_abs <- plot_ale_geff_abs +
    geom_point(data = obs_pts_geff, aes(x = x, y = y), inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  
  # Optional: overlay observations on RELATIVE panels (center to model mean)
  if (isTRUE(overlay_on_relative)) {
    obs_pts_et_rel   <- dplyr::mutate(obs_pts_et,   y = y - mu_ET_pred)
    obs_pts_geff_rel <- dplyr::mutate(obs_pts_geff, y = y - mu_geff_pred)
    plot_ale_et_rel  <- plot_ale_et_rel  + geom_point(data = obs_pts_et_rel,   aes(x = x, y = y),
                                                      inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    plot_ale_geff_rel<- plot_ale_geff_rel+ geom_point(data = obs_pts_geff_rel, aes(x = x, y = y),
                                                      inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  }
  
  # Save (ALE)
  ggsave("../Plots/ggplot2/ALE_ET_smoothed_abs.png",   plot_ale_et_abs,   width = 9, height = 6.5, dpi = 300)
  ggsave("../Plots/ggplot2/ALE_geff_smoothed_abs.png", plot_ale_geff_abs, width = 9, height = 6.5, dpi = 300)
  # Optionally:
  # ggsave("../Plots/ggplot2/ALE_ET_smoothed_rel.png",   plot_ale_et_rel,   width = 9, height = 6.5, dpi = 300)
  # ggsave("../Plots/ggplot2/ALE_geff_smoothed_rel.png", plot_ale_geff_rel, width = 9, height = 6.5, dpi = 300)
  
  # Expose
  jarvis_ale_raw      <- ale_all
  jarvis_ale_smoothed <- ale_all_sm
  
} else {
  # ======== OPTION B: 1-D SLICES (others held fixed) ========
  message("compute_ALE == FALSE: Plotting 1-D prediction slices with other variables held fixed.")
  
  # Helper: most frequent level for non-numerics
  .mode <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Build a single reference row:
  ref_row <- jarvis_df_opt[1, , drop = FALSE]
  num_cols <- names(jarvis_df_opt)[vapply(jarvis_df_opt, is.numeric, logical(1))]
  if (hold_others == "mean") {
    ref_row[1, num_cols] <- lapply(jarvis_df_opt[num_cols], function(col) mean(col, na.rm = TRUE))
  } else {
    ref_row[1, num_cols] <- lapply(jarvis_df_opt[num_cols], function(col) median(col, na.rm = TRUE))
  }
  non_num_cols <- setdiff(names(jarvis_df_opt), num_cols)
  if (length(non_num_cols)) {
    ref_modes <- lapply(jarvis_df_opt[non_num_cols], .mode)
    for (nm in non_num_cols) ref_row[[nm]] <- ref_modes[[nm]]
  }
  
  # Reference predictions at the fixed point (for centering)
  ref_preds <- predict_ET_geff(ref_row)
  ref_ET    <- as.numeric(ref_preds$ET)
  ref_geff  <- as.numeric(ref_preds$geff)
  
  # Sweep each target variable over its observed range, others fixed to ref_row
  make_slice <- function(v) {
    x_min <- min(jarvis_df_opt[[v]], na.rm = TRUE)
    x_max <- max(jarvis_df_opt[[v]], na.rm = TRUE)
    x_grid <- seq(x_min, x_max, length.out = slice_grid_points)
    df <- ref_row[rep(1, length(x_grid)), , drop = FALSE]
    df[[v]] <- x_grid
    preds <- predict_ET_geff(df)
    tibble(
      var = v,
      x = x_grid,
      ET_pred   = preds$ET,
      geff_pred = preds$geff
    )
  }
  slices <- map_dfr(vars, make_slice)
  
  # Relative (centered at the fixed-point prediction)
  slices <- slices %>%
    mutate(ET_rel = ET_pred - ref_ET,
           geff_rel = geff_pred - ref_geff)
  
  # ---- Build plots ----
  plot_ale_et_rel <- ggplot(slices, aes(x = x, y = ET_rel)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "ΔET vs fixed point (mm/yr)",
         title = "1-D model slices for ET — others held fixed") +
    theme_bw()
  
  plot_ale_geff_rel <- ggplot(slices, aes(x = x, y = geff_rel)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "Δg_eff vs fixed point (mm s^-1)",
         title = "1-D model slices for g_eff — others held fixed") +
    theme_bw()
  
  if (isTRUE(show_rug)) {
    plot_ale_et_rel  <- plot_ale_et_rel  + geom_rug(data = slices %>% distinct(var, x), aes(x = x), sides = "b", inherit.aes = FALSE, alpha = 0.2)
    plot_ale_geff_rel<- plot_ale_geff_rel+ geom_rug(data = slices %>% distinct(var, x), aes(x = x), sides = "b", inherit.aes = FALSE, alpha = 0.2)
  }
  
  # Anchored (absolute predictions)
  plot_ale_et_abs <- ggplot(slices, aes(x = x, y = ET_pred)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "ET (mm/yr)",
         title = "1-D model slices for ET — absolute predictions") +
    theme_bw()
  
  plot_ale_geff_abs <- ggplot(slices, aes(x = x, y = geff_pred)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(x = NULL, y = "g_eff (mm s^-1)",
         title = "1-D model slices for g_eff — absolute predictions") +
    theme_bw()
  
  # ---- Overlay observed points (anchored plots) ----
  obs_pts_et   <- make_obs_pts("ET")
  obs_pts_geff <- make_obs_pts("g_eff")
  plot_ale_et_abs   <- plot_ale_et_abs   +
    geom_point(data = obs_pts_et,   aes(x = x, y = y), inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  plot_ale_geff_abs <- plot_ale_geff_abs +
    geom_point(data = obs_pts_geff, aes(x = x, y = y), inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  
  # Optional: overlay observations on RELATIVE panels (center to fixed-point prediction)
  if (isTRUE(overlay_on_relative)) {
    obs_pts_et_rel   <- dplyr::mutate(obs_pts_et,   y = y - ref_ET)
    obs_pts_geff_rel <- dplyr::mutate(obs_pts_geff, y = y - ref_geff)
    plot_ale_et_rel  <- plot_ale_et_rel  + geom_point(data = obs_pts_et_rel,   aes(x = x, y = y),
                                                      inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
    plot_ale_geff_rel<- plot_ale_geff_rel+ geom_point(data = obs_pts_geff_rel, aes(x = x, y = y),
                                                      inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  }
  
  # Save (Slices)
  ggsave("../Plots/ggplot2/SLICES_ET_abs.png",   plot_ale_et_abs,   width = 9, height = 6.5, dpi = 300)
  ggsave("../Plots/ggplot2/SLICES_geff_abs.png", plot_ale_geff_abs, width = 9, height = 6.5, dpi = 300)
  # Optionally:
  # ggsave("../Plots/ggplot2/SLICES_ET_rel.png",   plot_ale_et_rel,   width = 9, height = 6.5, dpi = 300)
  # ggsave("../Plots/ggplot2/SLICES_geff_rel.png", plot_ale_geff_rel, width = 9, height = 6.5, dpi = 300)
  
  # Expose (named analogously)
  jarvis_ale_raw      <- NULL
  jarvis_ale_smoothed <- slices
}
