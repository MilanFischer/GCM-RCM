# ============================================================
# ALE curves for Jarvis model (ET and g_eff) — smoothed & anchored
# Save as: ./src/Jarvis_ALE_curves_smoothed.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------
# Requirements from the fit step (must exist in the session)
# -----------------------
req <- c("par_hat","df_opt","fVPD_mode","rhoAir","CpAir","gamma","model_ET_components")
miss <- setdiff(req, ls(envir = .GlobalEnv))
if (length(miss)) stop("Missing in environment: ", paste(miss, collapse = ", "))

# -----------------------
# Settings
# -----------------------
if (!exists("eps", inherits = TRUE)) eps <- 1e-6  # numeric guard
n_bins      <- 20    # ALE equal-frequency bins
smooth_span <- 0.5   # LOESS span (↑ smoother, ↓ more detail)

# -----------------------
# Helpers: derive columns Jarvis needs + predict ET & g_eff
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
  comp <- model_ET_components(par = par_hat, input_df = df2, fVPD_mode = fVPD_mode)
  list(ET = as.numeric(comp$ET_pred), geff = as.numeric(comp$g_eff_pred))
}

# -----------------------
# ALE for a single variable (Molnar definition)
# KEY: evaluate lo & hi TOGETHER so fVPD normalization is shared
# -----------------------
ale_one <- function(var, data = df_opt, n_bins = 20) {
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
  
  # accumulate and center (standard ALE)
  ale_ET <- cumsum(dET)
  ale_G  <- cumsum(dG)
  mids   <- (edges[-1] + edges[-length(edges)]) / 2
  
  N <- sum(cnt)
  if (N > 0) {
    mu_ET <- sum(ale_ET * cnt) / N
    mu_G  <- sum(ale_G  * cnt) / N
    ale_ET <- ale_ET - mu_ET
    ale_G  <- ale_G  - mu_G
  }
  
  tibble(var = var, x = mids, ALE_ET = ale_ET, ALE_geff = ale_G, n_in_bin = cnt)
}

# -----------------------
# Compute ALE (raw)
# -----------------------
vars <- c("P","Rg","Ta","VPD")
ale_list <- lapply(vars, ale_one, data = df_opt, n_bins = n_bins)
ale_all  <- dplyr::bind_rows(ale_list)

# -----------------------
# Smooth ALE & build anchored (absolute-scale) versions
# -----------------------
# Baseline mean predictions for anchoring
pred0 <- predict_ET_geff(df_opt)
mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)

# Safe LOESS predict helper (returns original if loess fails)
.loess_pred <- function(y, x, span) {
  fit <- try(loess(y ~ x, span = span), silent = TRUE)
  if (inherits(fit, "try-error")) return(y)
  as.numeric(predict(fit, data.frame(x = x)))
}

ale_all_sm <- ale_all %>%
  group_by(var) %>%
  mutate(
    ALE_ET_s   = .loess_pred(ALE_ET,   x, span = smooth_span),
    ALE_geff_s = .loess_pred(ALE_geff, x, span = smooth_span)
  ) %>%
  # re-center smoothed curves (keep "relative to mean" meaning)
  mutate(
    ALE_ET_s   = ALE_ET_s   - weighted.mean(ALE_ET_s,   w = n_in_bin, na.rm = TRUE),
    ALE_geff_s = ALE_geff_s - weighted.mean(ALE_geff_s, w = n_in_bin, na.rm = TRUE),
    # anchored absolute versions (add back mean prediction)
    ALE_ET_abs_s   = ALE_ET_s   + mu_ET_pred,
    ALE_geff_abs_s = ALE_geff_s + mu_geff_pred
  ) %>%
  ungroup()

# -----------------------
# Plotting
# -----------------------
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")

# Relative (centered) ALE
plot_ale_et_rel <- ggplot(ale_all_sm, aes(x = x, y = ALE_ET_s)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = bquote("ET (mm  yr"^"-1"*")"),
       title = "Accumulated Local Effects (ET) — smoothed, relative to mean") +
  theme_bw()

plot_ale_geff_rel <- ggplot(ale_all_sm, aes(x = x, y = ALE_geff_s)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = bquote("ALE of " * g[eff] * " (mm s"^"-1"*")"~"— centered"),
       title = bquote("Accumulated Local Effects (" * g[eff] * ") — smoothed, relative to mean")) +
  theme_bw()

# Anchored (absolute-scale) ALE
plot_ale_et_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_ET_abs_s)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = bquote("ET (mm  yr"^"-1"*")"),
       title = "Accumulated Local Effects (ET) — smoothed, absolute (anchored to mean)") +
  theme_bw()

plot_ale_geff_abs <- ggplot(ale_all_sm, aes(x = x, y = ALE_geff_abs_s)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = bquote(g[eff]~"(mm s"^"-1"*")"),
       title = bquote("Accumulated Local Effects (" * g[eff] * ") — smoothed, absolute (anchored to mean)")) +
  theme_bw()

# -----------------------
# Save (adjust paths as desired)
# -----------------------
ggsave("../Plots/ggplot2/ALE_ET_smoothed_rel.png",   plot_ale_et_rel,   width = 9, height = 6.5, dpi = 300)
ggsave("../Plots/ggplot2/ALE_geff_smoothed_rel.png", plot_ale_geff_rel, width = 9, height = 6.5, dpi = 300)
ggsave("../Plots/ggplot2/ALE_ET_smoothed_abs.png",   plot_ale_et_abs,   width = 9, height = 6.5, dpi = 300)
ggsave("../Plots/ggplot2/ALE_geff_smoothed_abs.png", plot_ale_geff_abs, width = 9, height = 6.5, dpi = 300)

# -----------------------
# Expose results for interactive use
# -----------------------
jarvis_ale_raw      <- ale_all        # un-smoothed (centered)
jarvis_ale_smoothed <- ale_all_sm     # smoothed + anchored variants
