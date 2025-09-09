suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(ranger)
  library(DEoptim)
  library(patchwork)   # only used to arrange a couple plots; safe to remove
})

# ===============================
# 0) Guards & small conveniences
# ===============================
# numeric guard for divisions etc.
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# ensure df_opt has inv_sqrt_VPD for optimization
if (!("inv_sqrt_VPD" %in% names(df_opt))) {
  stopifnot("VPD" %in% names(df_opt))
  df_opt <- df_opt |> mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))
}

# Small helper: recompute K_ET row-wise (no reliance on df$K_ET column)
recompute_K_ET <- function(d) {
  vVPD <- pmax(d$VPD, eps)
  vTa  <- d$Ta
  (rhoAir * CpAir / gamma) * vVPD * (1/1000) *
    1 / ((2.501e6 - 2361 * vTa) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

# ============================================
# 1) Fit parametric VPD-only g_eff model (DE)
# ============================================
obj_fun_geff_vpd <- function(par, input_df) {
  b0 <- par[[1]]; b1 <- par[[2]]
  g_hat <- b0 + b1 * input_df$inv_sqrt_VPD
  # keep predictions non-negative by penalization
  pen <- 1e6 * mean(pmax(-g_hat, 0)^2, na.rm = TRUE)
  rmse <- sqrt(mean((input_df$g_eff - g_hat)^2, na.rm = TRUE))
  rmse + pen
}

set.seed(if (exists("rng_seed")) rng_seed else 123)
de_geff <- DEoptim::DEoptim(
  fn       = obj_fun_geff_vpd,
  lower    = c(-10, -1e3),
  upper    = c( 10,  1e3),
  input_df = df_opt,
  control  = DEoptim::DEoptim.control(
    NP = 60L, itermax = 2000L, reltol = 1e-8,
    trace = if (exists("trace_flag")) trace_flag else FALSE
  )
)

b_hat <- de_geff$optim$bestmem
names(b_hat) <- c("b0_VPD","b1_VPD")
message(sprintf("g_eff VPD-only fit: b0=%.6f, b1=%.6f", b_hat[1], b_hat[2]))

# parametric predictions on full df
df <- df |> mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))
g_eff_vpd_hat <- with(df, b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD)

# ================================================
# 2) Random Forest residual model (flexible inputs)
# ================================================
# We train RF on residuals vs available predictors among {Rg, Ta, P, CO2_term}.
rf_candidate_vars <- c("Rg", "Ta", "P", "CO2_term")
rf_vars_present   <- intersect(rf_candidate_vars, names(df_opt))
stopifnot(length(rf_vars_present) >= 1)

df_rf <- df_opt |>
  mutate(resid_vpd = g_eff - (b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD)) |>
  select(resid_vpd, all_of(rf_vars_present)) |>
  filter(if_all(where(is.numeric), ~ is.finite(.x)))
stopifnot(nrow(df_rf) > 0)

set.seed(if (exists("rng_seed")) rng_seed else 123)
rf_fit <- ranger(
  resid_vpd ~ .,
  data = df_rf,
  num.trees = 800,
  mtry = min(3, length(rf_vars_present)),
  min.node.size = 5,
  importance = "impurity",
  seed = if (exists("rng_seed")) rng_seed else 123,
  oob.error = TRUE
)

oob_rmse <- sqrt(rf_fit$prediction.error)
message(sprintf("RF residual model OOB RMSE = %.3f", oob_rmse))

# Helper to predict RF residuals safely on any data.frame with columns available
.predict_rf_resid <- function(newdata) {
  vars <- intersect(rf_vars_present, colnames(newdata))
  # ranger needs at least one column
  if (length(vars) == 0) return(rep(0, nrow(newdata)))
  as.numeric(predict(rf_fit, data = newdata[, vars, drop = FALSE])$predictions)
}

# RF residual predictions for all rows
rf_pred <- .predict_rf_resid(df)

# ================================
# 3) Final predictions & diagnostics
# ================================
g_eff_final <- g_eff_vpd_hat
ok <- is.finite(rf_pred)
g_eff_final[ok] <- g_eff_vpd_hat[ok] + rf_pred[ok]

# ET prediction using recomputed K_ET row-wise
ET_pred_final <- recompute_K_ET(df) * g_eff_final

output_df_geff <- df |>
  mutate(
    g_eff_vpd_hat = g_eff_vpd_hat,
    resid_vpd     = g_eff - g_eff_vpd_hat,
    resid_rf_hat  = rf_pred,
    g_eff_final   = g_eff_final,
    g_eff_resids  = g_eff - g_eff_final,
    ET_pred_final = ET_pred_final
  )

# Fit metrics
RMSE_vpd   <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_eff_vpd_hat)^2, na.rm = TRUE))
R2_vpd   <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_eff_vpd_hat, use="complete.obs")^2)
RMSE_final <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_eff_final)^2, na.rm = TRUE))
R2_final   <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_eff_final, use="complete.obs")^2)

message(sprintf("Stage 1 (VPD-only): RMSE = %.3f, R^2 = %.3f", RMSE_vpd, R2_vpd))
message(sprintf("Stage 2 (VPD + RF): RMSE = %.3f, R^2 = %.3f", RMSE_final, R2_final))

# Quick observed vs predicted scatter (g_eff, ET)
p_geff <- output_df_geff |>
  transmute(g_eff_obs = g_eff, g_eff_pred = g_eff_final) |>
  filter(is.finite(g_eff_obs), is.finite(g_eff_pred)) |>
  ggplot(aes(g_eff_obs, g_eff_pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff])),
    title = expression(paste("Observed vs Predicted ", g[eff])),
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", RMSE_final, R2_final)
  )

p_et <- output_df_geff |>
  transmute(ET_obs = ET, ET_pred = ET_pred_final) |>
  filter(is.finite(ET_obs), is.finite(ET_pred)) |>
  ggplot(aes(ET_obs, ET_pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(x = "Observed ET", y = "Predicted ET", title = "Observed vs Predicted ET")

print(p_geff + p_et)

# ====================================================
# 4) ALE (Accumulated Local Effects) for g_eff and ET
#     — avoids unrealistic combos typical for PDP
# ====================================================

# Model prediction (g_eff & ET) for an arbitrary newdata
predict_geff_ET <- function(newdata) {
  inv_sqrt_VPD <- 1 / sqrt(pmax(newdata$VPD, eps))
  g_param      <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD
  rf_resid     <- .predict_rf_resid(newdata)
  g_eff_hat    <- g_param + rf_resid
  K_ET_new     <- recompute_K_ET(newdata)
  ET_hat       <- K_ET_new * g_eff_hat
  list(geff = as.numeric(g_eff_hat), ET = as.numeric(ET_hat))
}

# 1D ALE for a single variable (equal-frequency bins)
ale_one <- function(var, data = df, n_bins = 50) {
  stopifnot(var %in% names(data))
  z_all <- data[[var]]
  z <- z_all[is.finite(z_all)]
  if (!length(z)) stop("Variable ", var, " has no finite values.")
  
  # bin edges (equal frequency), over full observed range
  edges <- as.numeric(quantile(z, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  edges[1] <- min(z, na.rm = TRUE)
  edges[length(edges)] <- max(z, na.rm = TRUE)
  
  dET <- numeric(n_bins)
  dG  <- numeric(n_bins)
  cnt <- integer(n_bins)
  
  for (j in seq_len(n_bins)) {
    lo <- edges[j]; hi <- edges[j + 1]
    idx <- which(z_all >= lo & (z_all < hi | (j == n_bins & z_all <= hi)))
    cnt[j] <- length(idx)
    if (cnt[j] == 0) { dET[j] <- 0; dG[j] <- 0; next }
    
    df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
    df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
    
    pred_lo <- predict_geff_ET(df_lo)
    pred_hi <- predict_geff_ET(df_hi)
    
    dET[j] <- mean(pred_hi$ET   - pred_lo$ET,   na.rm = TRUE)
    dG[j]  <- mean(pred_hi$geff - pred_lo$geff, na.rm = TRUE)
  }
  
  # accumulate local effects and center (ALE definition)
  ale_ET <- cumsum(dET)
  ale_G  <- cumsum(dG)
  mids   <- (edges[-1] + edges[-length(edges)]) / 2
  N      <- sum(cnt)
  
  if (N > 0) {
    ale_ET <- ale_ET - sum(ale_ET * cnt, na.rm = TRUE) / N
    ale_G  <- ale_G  - sum(ale_G  * cnt, na.rm = TRUE) / N
  }
  
  tibble(var = var, x = mids, ALE_ET = ale_ET, ALE_geff = ale_G, n_in_bin = cnt)
}

# Weighted loess with safe fallback to linear interpolation
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

# Choose variables to interpret (four predictors requested)
vars_for_ALE <- intersect(c("VPD","Rg","Ta","P"), names(df))
stopifnot(length(vars_for_ALE) >= 1)

# Compute raw ALE (centered, at bin midpoints)
n_bins <- 50L
ale_raw_list <- lapply(vars_for_ALE, ale_one, data = df, n_bins = n_bins)
ale_raw <- bind_rows(ale_raw_list)

# Anchor ALE to absolute predictions (adds model mean back)
pred0 <- predict_geff_ET(df)
mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)

# Smooth ALE to a dense grid (per var), then re-center within var (preserve ALE meaning)
smooth_span <- 0.5
smooth_grid_points <- 300L
var_labels <- c(VPD="VPD (kPa)", Rg="Rg (radiation)", Ta="Ta (°C)", P="P (precip.)")

smooth_one_var <- function(vname) {
  dv <- ale_raw |> filter(var == vname)
  x_min <- min(df[[vname]], na.rm = TRUE)
  x_max <- max(df[[vname]], na.rm = TRUE)
  x_grid <- seq(x_min, x_max, length.out = smooth_grid_points)
  
  fit_ET <- .fit_loess_safe(dv$ALE_ET,   dv$x, dv$n_in_bin, smooth_span)
  fit_G  <- .fit_loess_safe(dv$ALE_geff, dv$x, dv$n_in_bin, smooth_span)
  
  y_ET_grid <- .predict_at(fit_ET, dv$x, dv$ALE_ET,   x_grid)
  y_G_grid  <- .predict_at(fit_G,  dv$x, dv$ALE_geff, x_grid)
  
  # re-center smoothed curves to zero-mean (weighted) and anchor to absolute predictions
  y_ET_mid <- .predict_at(fit_ET, dv$x, dv$ALE_ET,   dv$x)
  y_G_mid  <- .predict_at(fit_G,  dv$x, dv$ALE_geff, dv$x)
  off_ET <- weighted.mean(y_ET_mid, w = dv$n_in_bin, na.rm = TRUE)
  off_G  <- weighted.mean(y_G_mid,  w = dv$n_in_bin, na.rm = TRUE)
  
  tibble(
    var = vname,
    x = x_grid,
    ALE_ET_rel   = y_ET_grid - off_ET,
    ALE_geff_rel = y_G_grid  - off_G,
    ALE_ET_abs   = (y_ET_grid - off_ET) + mu_ET_pred,
    ALE_geff_abs = (y_G_grid  - off_G)  + mu_geff_pred
  )
}
ale_smoothed <- map_dfr(vars_for_ALE, smooth_one_var)

# Observed points for overlay on anchored panels (helps sanity-check alignment)
obs_max_points <- 5000L
make_obs_pts <- function(y_name) {
  out <- bind_rows(lapply(vars_for_ALE, function(v) {
    tibble(var = v, x = df[[v]], y = df[[y_name]])
  })) |> filter(is.finite(x), is.finite(y))
  if (nrow(out) > obs_max_points) {
    out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
  }
  out
}
obs_pts_ET   <- make_obs_pts("ET")
obs_pts_geff <- make_obs_pts("g_eff")

# -----------------
# ALE PLOTS
# -----------------

# Centered (relative) ALE — shows *effect shape* per variable
plot_ale_et_rel <- ggplot(ale_smoothed, aes(x = x, y = ALE_ET_rel)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = "ALE(ET) — centered", title = "Accumulated Local Effects for ET") +
  theme_bw()

plot_ale_geff_rel <- ggplot(ale_smoothed, aes(x = x, y = ALE_geff_rel)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = expression(ALE~of~g[eff]~"— centered"),
       title = expression("Accumulated Local Effects for " * g[eff])) +
  theme_bw()

# Anchored (absolute) ALE — adds model mean; overlay observed points
plot_ale_et_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_ET_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = "ET (mm/yr) — anchored to model mean",
       title = "ALE (ET): absolute scale") +
  theme_bw() +
  geom_point(data = obs_pts_ET, aes(x = x, y = y),
             inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25)

plot_ale_geff_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_geff_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*") — anchored to model mean"),
       title = expression("ALE (" * g[eff] * "): absolute scale")) +
  theme_bw() +
  geom_point(data = obs_pts_geff, aes(x = x, y = y),
             inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25)

# Show ALE plots
print(plot_ale_et_rel);  print(plot_ale_geff_rel)
print(plot_ale_et_abs);  print(plot_ale_geff_abs)

# ======================
# (Optional) Save plots
# ======================
# ggsave("../Plots/ggplot2/obs_vs_pred_geff_RF.png", p_geff, width = 6, height = 5, dpi = 300)
# ggsave("../Plots/ggplot2/obs_vs_pred_ET_RF.png",   p_et,   width = 6, height = 5, dpi = 300)
# ggsave("../Plots/ggplot2/ALE_ET_rel.png",   plot_ale_et_rel,   width = 9, height = 6.5, dpi = 300)
# ggsave("../Plots/ggplot2/ALE_geff_rel.png", plot_ale_geff_rel, width = 9, height = 6.5, dpi = 300)
# ggsave("../Plots/ggplot2/ALE_ET_abs.png",   plot_ale_et_abs,   width = 9, height = 6.5, dpi = 300)
# ggsave("../Plots/ggplot2/ALE_geff_abs.png", plot_ale_geff_abs, width = 9, height = 6.5, dpi = 300)

# ==============================
# (Optional) Save model bundle
# ==============================
geff_vpd_rf_bundle <- list(
  metadata   = if (exists("metadata")) paste0(metadata, " | g_eff VPD-only + RF residuals + ALE.") else "g_eff VPD-only + RF residuals + ALE",
  b0_VPD     = b_hat[["b0_VPD"]],
  b1_VPD     = b_hat[["b1_VPD"]],
  rf_fit     = rf_fit,
  rf_vars    = rf_vars_present,
  df_opt     = df_opt,
  out_df     = output_df_geff,
  RMSE_vpd   = RMSE_vpd,
  RMSE_final = RMSE_final,
  R2_final   = R2_final,
  rng_seed   = if (exists("rng_seed")) rng_seed else NA_integer_,
  ale_raw    = ale_raw,
  ale_smooth = ale_smoothed
)
if (exists("out_file")) {
  save(geff_vpd_rf_bundle, file = sub("\\.RData$", "_geff_vpd_rf_ALE.RData", out_file))
}
