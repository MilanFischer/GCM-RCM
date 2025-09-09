# ============================================================
# ALE curves for Jarvis model (ET and g_eff) — split per variable
# Each variable gets its own standalone plot => y-axis ticks visible on all
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------
# Requirements (must exist)
# -----------------------
req <- c("jarvis_pars","jarvis_df_opt","jarvis_fVPD_mode","rhoAir","CpAir","gamma","jarvis_model_ET_components")
miss <- setdiff(req, ls(envir = .GlobalEnv))
if (length(miss)) stop("Missing in environment: ", paste(miss, collapse = ", "))

# -----------------------
# Settings
# -----------------------
if (!exists("eps", inherits = TRUE)) eps <- 1e-6  # numeric guard
n_bins              <- 20     # ALE equal-frequency bins per variable
smooth_span         <- 0.5    # LOESS span (↑ smoother, ↓ more detail)
show_rug            <- TRUE   # show data-support rug on relative plots
overlay_on_relative <- FALSE # also overlay obs on relative (centered) plots
obs_max_points      <- 5000   # downsample cap for observed scatter per plot

out_dir <- "../Plots/ggplot2"  # change if you like
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Helpers: derive cols Jarvis needs + predict ET & g_eff
# -----------------------
rederive_cols <- function(df) {
  df %>%
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
# ALE for one variable (Molnar). IMPORTANT: evaluate lo & hi TOGETHER
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
# Compute ALE (raw)
# -----------------------
vars <- c("P","Rg","Ta","VPD")
ale_list <- lapply(vars, ale_one, data = jarvis_df_opt, n_bins = n_bins)
ale_all  <- dplyr::bind_rows(ale_list)

# -----------------------
# Smooth ALE & anchor to mean prediction
# -----------------------
pred0 <- predict_ET_geff(jarvis_df_opt)
mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)

.loess_pred <- function(y, x, w, span) {
  fit <- try(loess(y ~ x, weights = w, span = span,
                   control = loess.control(surface = "direct")), silent = TRUE)
  if (inherits(fit, "try-error")) return(y)
  as.numeric(predict(fit, data.frame(x = x)))
}

ale_all_sm <- ale_all %>%
  group_by(var) %>%
  mutate(
    ALE_ET_s   = .loess_pred(ALE_ET,   x, n_in_bin, smooth_span),
    ALE_geff_s = .loess_pred(ALE_geff, x, n_in_bin, smooth_span)
  ) %>%
  mutate(
    ALE_ET_s   = ALE_ET_s   - weighted.mean(ALE_ET_s,   w = n_in_bin, na.rm = TRUE),
    ALE_geff_s = ALE_geff_s - weighted.mean(ALE_geff_s, w = n_in_bin, na.rm = TRUE),
    ALE_ET_abs_s   = ALE_ET_s   + mu_ET_pred,
    ALE_geff_abs_s = ALE_geff_s + mu_geff_pred
  ) %>%
  ungroup()

# Pretty names
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")

# -----------------------
# Observed points (for anchored overlays)
# -----------------------
make_obs_pts <- function(y_name) {
  out <- dplyr::bind_rows(lapply(vars, function(v) {
    tibble(var = v, x = jarvis_df_opt[[v]], y = jarvis_df_opt[[y_name]])
  })) %>% dplyr::filter(is.finite(x), is.finite(y))
  if (nrow(out) > obs_max_points) {
    out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
  }
  out
}
obs_pts_et   <- make_obs_pts("ET")
obs_pts_geff <- make_obs_pts("g_eff")

# -----------------------
# Plot builders (SPLIT per variable)
# -----------------------
theme_base <- theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.ticks.y = element_line()
  )

build_rel_plot <- function(v, what = c("ET","geff")) {
  what <- match.arg(what)
  dat  <- dplyr::filter(ale_all_sm, var == v)
  y    <- if (what == "ET") "ALE_ET_s" else "ALE_geff_s"
  ylab <- if (what == "ET") "ALE of ET (mm/yr) — centered" else "ALE of g_eff (mm s^-1) — centered"
  ttl  <- paste0("Accumulated Local Effects (", if (what=="ET") "ET" else "g_eff",
                 ") — smoothed, relative to mean  |  ", var_labels[[v]])
  
  p <- ggplot(dat, aes(x = x, y = .data[[y]])) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.4) +
    geom_line(size = 0.9) +
    labs(x = NULL, y = ylab, title = ttl) +
    scale_y_continuous(breaks = waiver()) +
    theme_base
  
  if (isTRUE(show_rug)) {
    # rug needs bin midpoints from the raw ALE for that var
    p <- p + geom_rug(data = dplyr::filter(ale_all, var == v),
                      aes(x = x), sides = "b", inherit.aes = FALSE, alpha = 0.2)
  }
  
  if (isTRUE(overlay_on_relative)) {
    mu <- if (what == "ET") mu_ET_pred else mu_geff_pred
    obs <- if (what == "ET") obs_pts_et else obs_pts_geff
    obs <- dplyr::filter(obs, var == v) %>% dplyr::mutate(y = y - mu)
    p <- p + geom_point(data = obs, aes(x = x, y = y),
                        inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25, color = "black")
  }
  
  p
}

build_abs_plot <- function(v, what = c("ET","geff")) {
  what <- match.arg(what)
  dat  <- dplyr::filter(ale_all_sm, var == v)
  y    <- if (what == "ET") "ALE_ET_abs_s" else "ALE_geff_abs_s"
  ylab <- if (what == "ET") "ET (mm/yr) — anchored" else "g_eff (mm s^-1) — anchored"
  ttl  <- paste0("Accumulated Local Effects (", if (what=="ET") "ET" else "g_eff",
                 ") — smoothed, absolute (anchored to mean)  |  ", var_labels[[v]])
  
  p <- ggplot(dat, aes(x = x, y = .data[[y]])) +
    geom_line(size = 0.9) +
    labs(x = NULL, y = ylab, title = ttl) +
    scale_y_continuous(breaks = waiver()) +
    theme_base
  
  # overlay observed points
  obs <- if (what == "ET") obs_pts_et else obs_pts_geff
  p <- p + geom_point(
    data = dplyr::filter(obs, var == v),
    mapping = aes(x = x, y = y),
    inherit.aes = FALSE,
    shape = 3, size = 0.5, alpha = 0.25, color = "black"
  )
  
  p
}

# -----------------------
# Create & save FOUR plots per set (P, Rg, Ta, VPD)
# -----------------------
for (v in vars) {
  # Relative ET
  g <- build_rel_plot(v, "ET")
  ggsave(file.path(out_dir, sprintf("ALE_ET_smoothed_rel_%s.png", v)),
         g, width = 6, height = 4.5, dpi = 300)
  
  # Relative g_eff
  g <- build_rel_plot(v, "geff")
  ggsave(file.path(out_dir, sprintf("ALE_geff_smoothed_rel_%s.png", v)),
         g, width = 6, height = 4.5, dpi = 300)
  
  # Anchored ET
  g <- build_abs_plot(v, "ET")
  ggsave(file.path(out_dir, sprintf("ALE_ET_smoothed_abs_%s.png", v)),
         g, width = 6, height = 4.5, dpi = 300)
  
  # Anchored g_eff
  g <- build_abs_plot(v, "geff")
  ggsave(file.path(out_dir, sprintf("ALE_geff_smoothed_abs_%s.png", v)),
         g, width = 6, height = 4.5, dpi = 300)
}

# -----------------------
# Expose results for interactive use
# -----------------------
jarvis_ale_raw      <- ale_all        # un-smoothed (centered)
jarvis_ale_smoothed <- ale_all_sm     # smoothed + anchored variants
