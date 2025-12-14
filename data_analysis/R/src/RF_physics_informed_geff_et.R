################################################################################
## RF-only g_eff + ET model (one-stage), with physics-style VPD transform,
## CO2 regime monotonicity gate, permutation importance, and PDPs.
##
## Model:
##   g_eff ~ RF( inv_sqrt_VPD, Ta, P, Rg, CO2_term )
##   ET_hat = K_ET(VPD, Ta, ...) * g_eff_hat
##
## CO2 handling:
##   - Model uses CO2_term (dimensionless, e.g. log(CO2_ppm / CO2_1981_2005))
##   - We also create CO2_ppm (for plotting/PDP axis), including EUR-44 (RCMs)
##
## CO2 gate:
##   - Keep CO2 only if predicted g_eff is monotone non-increasing across regimes
##     (EUR44 -> CMIP5 -> CMIP6) for at least "majority_frac" of rows AND mean Δ < tol.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(ranger)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 0) Knobs & basic checks
# ------------------------------------------------------------------------------

eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6
perm_B <- if (exists("perm_B", inherits = TRUE)) perm_B else 200L

rf_num_trees <- if (exists("rf_num_trees", inherits = TRUE)) rf_num_trees else 1200L
rf_tune_mode <- if (exists("rf_tune_mode", inherits = TRUE)) rf_tune_mode else "full"

# CO2 gate knobs
co2_neg_majority_frac <- if (exists("co2_neg_majority_frac", inherits = TRUE)) co2_neg_majority_frac else 0.60
co2_neg_total_tol     <- if (exists("co2_neg_total_tol",     inherits = TRUE)) co2_neg_total_tol     else -1e-6

stopifnot(exists("df_opt"), exists("df"))
stopifnot(all(c("VPD","Ta","P","Rg","ET","g_eff") %in% names(df_opt)))
stopifnot(all(c("VPD","Ta","P","Rg","ET","g_eff") %in% names(df)))

# Optional: harmonize ensemble column name if you sometimes have ENSEMBLE
if (!"ensemble" %in% names(df_opt) && "ENSEMBLE" %in% names(df_opt)) df_opt <- df_opt %>% mutate(ensemble = ENSEMBLE)
if (!"ensemble" %in% names(df)     && "ENSEMBLE" %in% names(df))     df     <- df     %>% mutate(ensemble = ENSEMBLE)

# Optional: harmonize PERIOD column name if you sometimes have Period
if (!"PERIOD" %in% names(df_opt) && "Period" %in% names(df_opt)) df_opt <- df_opt %>% mutate(PERIOD = Period)
if (!"PERIOD" %in% names(df)     && "Period" %in% names(df))     df     <- df     %>% mutate(PERIOD = Period)

# Ensure VPD positive (needed for 1/sqrt(VPD))
df_opt <- df_opt %>% mutate(VPD = pmax(VPD, eps))
df     <- df     %>% mutate(VPD = pmax(VPD, eps))

# ------------------------------------------------------------------------------
# 0b) Build CO2_ppm (including EUR44) and CO2_term if needed
# ------------------------------------------------------------------------------

# We interpret EUR-44 (RCMs) as "no CO2 physiology", but for plotting on ppm axis,
# we place them at reference CO2 (CO2_1981_2005). This is only for axis convenience.
if (exists("CO2_1981_2005", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE) &&
    "PERIOD" %in% names(df_opt) && "ensemble" %in% names(df_opt) &&
    "PERIOD" %in% names(df)     && "ensemble" %in% names(df)) {
  
  df_opt <- df_opt %>%
    mutate(
      CO2_ppm = case_when(
        grepl("EUR", ensemble, ignore.case = TRUE) ~ CO2_1981_2005,
        PERIOD == "1981_2005" ~ CO2_1981_2005,
        PERIOD == "2076_2100" & ensemble == "CMIP5" ~ CO2_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ensemble == "CMIP6" ~ CO2_2076_2100_RCP85_CMIP6,
        TRUE ~ NA_real_
      )
    )
  
  df <- df %>%
    mutate(
      CO2_ppm = case_when(
        grepl("EUR", ensemble, ignore.case = TRUE) ~ CO2_1981_2005,
        PERIOD == "1981_2005" ~ CO2_1981_2005,
        PERIOD == "2076_2100" & ensemble == "CMIP5" ~ CO2_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ensemble == "CMIP6" ~ CO2_2076_2100_RCP85_CMIP6,
        TRUE ~ NA_real_
      )
    )
} else {
  # If you don't have those scalars/columns, we still run the RF without CO2_ppm.
  if (!"CO2_ppm" %in% names(df_opt)) df_opt$CO2_ppm <- NA_real_
  if (!"CO2_ppm" %in% names(df))     df$CO2_ppm     <- NA_real_
}

# If CO2_term doesn't exist, create it from CO2_ppm (log ratio to CO2_1981_2005),
# but only where CO2_ppm is available.
if (!"CO2_term" %in% names(df_opt)) df_opt$CO2_term <- NA_real_
if (!"CO2_term" %in% names(df))     df$CO2_term     <- NA_real_

if (exists("CO2_1981_2005", inherits = TRUE)) {
  df_opt <- df_opt %>%
    mutate(CO2_term = if_else(is.finite(CO2_term),
                              CO2_term,
                              if_else(is.finite(CO2_ppm), log(CO2_ppm / CO2_1981_2005), NA_real_)))
  df <- df %>%
    mutate(CO2_term = if_else(is.finite(CO2_term),
                              CO2_term,
                              if_else(is.finite(CO2_ppm), log(CO2_ppm / CO2_1981_2005), NA_real_)))
}

# ------------------------------------------------------------------------------
# 1) Create physics-style VPD feature: inv_sqrt_VPD = 1/sqrt(VPD)
# ------------------------------------------------------------------------------

df_opt <- df_opt %>% mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))
df     <- df     %>% mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))

# ------------------------------------------------------------------------------
# 2) Fit RF for g_eff (one-stage)
# ------------------------------------------------------------------------------

# Candidate predictors (VPD enters only via inv_sqrt_VPD)
rf_vars_full <- c("inv_sqrt_VPD", "Ta", "P", "Rg", "CO2_term")
rf_vars_full <- intersect(rf_vars_full, names(df_opt))

# Build RF training table (drop non-finite numeric rows)
df_rf_geff <- df_opt %>%
  select(g_eff, all_of(rf_vars_full)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

if (!nrow(df_rf_geff)) stop("No usable rows in df_rf_geff after filtering finiteness.")
if (length(setdiff(names(df_rf_geff), "g_eff")) < 1) stop("No predictors available for RF g_eff.")

# --- Train RF (tuned or fixed) ---
set.seed(if (exists("rng_seed")) rng_seed else 123)

if (rf_tune_mode == "none") {
  rf_geff <- fit_rf_fixed(
    df_rf_geff,
    target_col = "g_eff",
    num.trees  = rf_num_trees,
    seed       = if (exists("rng_seed")) rng_seed else 123
  )
  best_params_geff <- list(
    mtry            = rf_geff$mtry,
    min.node.size   = rf_geff$min.node.size,
    sample.fraction = 1.0,
    splitrule       = "variance",
    max.depth       = 0L
  )
} else {
  if (rf_tune_mode == "quick") {
    tuned_geff <- tune_rf_oob(
      df_rf_geff,
      target_col = "g_eff",
      num.trees       = rf_num_trees,
      min.node.size   = c(3L,5L,10L),
      sample.fraction = c(0.8,1.0),
      splitrule       = c("variance","extratrees"),
      max.depth       = c(0L,12L)
    )
  } else {
    tuned_geff <- tune_rf_oob(
      df_rf_geff,
      target_col = "g_eff",
      num.trees  = rf_num_trees
    )
  }
  rf_geff          <- tuned_geff$best_fit
  best_params_geff <- tuned_geff$best_params
}

oob_rmse_geff <- sqrt(rf_geff$prediction.error)
rf_vars_present <- setdiff(names(df_rf_geff), "g_eff")

message(sprintf("RF g_eff OOB RMSE = %.3f", oob_rmse_geff))
message(sprintf("RF params: mtry=%s, min.node.size=%s, sample.fraction=%.2f, splitrule=%s, max.depth=%s",
                best_params_geff$mtry, best_params_geff$min.node.size,
                best_params_geff$sample.fraction,
                as.character(best_params_geff$splitrule),
                best_params_geff$max.depth))

# Prediction helper for g_eff
.predict_rf_geff <- function(newdata) {
  vars <- intersect(rf_vars_present, names(newdata))
  if (!length(vars)) stop("No RF predictors found in newdata.")
  as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
}

# ------------------------------------------------------------------------------
# 3) Full prediction function: g_eff_hat and ET_hat
# ------------------------------------------------------------------------------

predict_geff_ET_RF <- function(newdata) {
  geff_hat <- .predict_rf_geff(newdata)
  K_ET     <- recompute_K_ET(newdata)     # must exist (your physics conversion)
  ET_hat   <- K_ET * geff_hat
  list(geff = as.numeric(geff_hat), ET = as.numeric(ET_hat))
}

# ------------------------------------------------------------------------------
# 4) CO2 regime monotonicity gate (EUR -> CMIP5 -> CMIP6)
#    Keep CO2 only if g_eff is monotone non-increasing in most rows AND mean Δ < tol
# ------------------------------------------------------------------------------

co2_regime_gate <- function(data, co2_levels,
                            majority_frac = 0.60,
                            tol = -1e-6) {
  stopifnot("CO2_term" %in% names(data))
  stopifnot(length(co2_levels) >= 2)
  
  preds <- lapply(co2_levels, function(val) {
    d <- data
    d$CO2_term <- val
    predict_geff_ET_RF(d)$geff
  })
  
  deltas <- Map(function(a, b) b - a, preds[-length(preds)], preds[-1])
  
  # fraction of row-wise comparisons (across all consecutive regime steps) that are <= tol
  frac_neg <- mean(unlist(lapply(deltas, function(d) d <= tol)), na.rm = TRUE)
  
  # mean overall change from lowest to highest regime
  mean_total <- mean(preds[[length(preds)]] - preds[[1]], na.rm = TRUE)
  
  list(frac_neg = frac_neg,
       mean_total = mean_total,
       pass = (frac_neg >= majority_frac && mean_total <= tol))
}

# Only attempt the gate if CO2_term is in the RF predictors AND regime constants exist
if ("CO2_term" %in% rf_vars_present &&
    exists("CO2_1981_2005", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
  
  co2_levels <- c(
    0,  # EUR-44 (no physiology)
    log(CO2_2076_2100_RCP85_CMIP5 / CO2_1981_2005),
    log(CO2_2076_2100_RCP85_CMIP6 / CO2_1981_2005)
  )
  
  gchk <- co2_regime_gate(df_opt, co2_levels,
                          majority_frac = co2_neg_majority_frac,
                          tol = co2_neg_total_tol)
  
  message(sprintf("CO2 gate: frac negative = %.1f%%, mean Δg_eff(hi-low) = %.4g",
                  100 * gchk$frac_neg, gchk$mean_total))
  
  if (!isTRUE(gchk$pass)) {
    message("Dropping CO2_term due to failed monotonicity gate, then refitting RF.")
    
    # Refit RF without CO2_term
    rf_vars_present2 <- setdiff(rf_vars_present, "CO2_term")
    
    df_rf_geff2 <- df_opt %>%
      select(g_eff, all_of(rf_vars_present2)) %>%
      filter(if_all(where(is.numeric), ~ is.finite(.x)))
    
    set.seed(if (exists("rng_seed")) rng_seed else 123)
    
    if (rf_tune_mode == "none") {
      rf_geff <- fit_rf_fixed(df_rf_geff2, target_col = "g_eff",
                              num.trees = rf_num_trees,
                              seed = if (exists("rng_seed")) rng_seed else 123)
    } else {
      # keep same tuning strategy
      if (rf_tune_mode == "quick") {
        tuned2 <- tune_rf_oob(df_rf_geff2, "g_eff", num.trees = rf_num_trees,
                              min.node.size = c(3L,5L,10L),
                              sample.fraction = c(0.8,1.0),
                              splitrule = c("variance","extratrees"),
                              max.depth = c(0L,12L))
      } else {
        tuned2 <- tune_rf_oob(df_rf_geff2, "g_eff", num.trees = rf_num_trees)
      }
      rf_geff <- tuned2$best_fit
    }
    
    rf_vars_present <- setdiff(names(df_rf_geff2), "g_eff")
    oob_rmse_geff <- sqrt(rf_geff$prediction.error)
    message(sprintf("Refit RF (no CO2_term): OOB RMSE = %.3f", oob_rmse_geff))
    
    # Update predictor helper to use new predictor set
    .predict_rf_geff <- function(newdata) {
      vars <- intersect(rf_vars_present, names(newdata))
      if (!length(vars)) stop("No RF predictors found in newdata.")
      as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
    }
    
    # Update full predictor
    predict_geff_ET_RF <- function(newdata) {
      geff_hat <- .predict_rf_geff(newdata)
      K_ET     <- recompute_K_ET(newdata)
      ET_hat   <- K_ET * geff_hat
      list(geff = as.numeric(geff_hat), ET = as.numeric(ET_hat))
    }
  } else {
    message("CO2_term passes monotonicity gate → keeping CO2_term.")
  }
}

# ------------------------------------------------------------------------------
# 5) Final predictions on df + metrics
# ------------------------------------------------------------------------------

pred_full <- predict_geff_ET_RF(df)

output_df <- df %>%
  mutate(
    g_eff_RF = pred_full$geff,
    ET_RF    = pred_full$ET
  )

rmse_geff <- sqrt(mean((output_df$g_eff - output_df$g_eff_RF)^2, na.rm = TRUE))
r2_geff   <- suppressWarnings(cor(output_df$g_eff, output_df$g_eff_RF, use = "complete.obs")^2)

rmse_ET <- sqrt(mean((output_df$ET - output_df$ET_RF)^2, na.rm = TRUE))
r2_ET   <- suppressWarnings(cor(output_df$ET, output_df$ET_RF, use = "complete.obs")^2)

message(sprintf("Final model: g_eff RMSE=%.3f, R^2=%.3f | ET RMSE=%.3f, R^2=%.3f",
                rmse_geff, r2_geff, rmse_ET, r2_ET))

# ------------------------------------------------------------------------------
# 6) Observed vs Predicted plots (g_eff and ET)
# ------------------------------------------------------------------------------

p_obs_pred_geff <- output_df %>%
  filter(is.finite(g_eff), is.finite(g_eff_RF)) %>%
  ggplot(aes(x = g_eff, y = g_eff_RF)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff])),
    title = expression(paste("Observed vs Predicted ", g[eff], " (RF-only, 1/sqrt(VPD))")),
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", rmse_geff, r2_geff)
  )

p_obs_pred_ET <- output_df %>%
  filter(is.finite(ET), is.finite(ET_RF)) %>%
  ggplot(aes(x = ET, y = ET_RF)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(
    x = "Observed ET",
    y = "Predicted ET",
    title = "Observed vs Predicted ET (ET = K_ET * g_eff_RF)",
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", rmse_ET, r2_ET)
  )

print(p_obs_pred_geff + p_obs_pred_ET)

# ------------------------------------------------------------------------------
# 7) Permutation importance for BOTH targets (g_eff and ET)
#    Importance is ΔRMSE (in target units) caused by permuting ONE predictor.
# ------------------------------------------------------------------------------

perm_importance_geff_ET <- function(data, vars, B = 200L, seed = 2024) {
  vars <- intersect(vars, names(data))
  if (!length(vars)) stop("No requested vars exist in data.")
  
  pred0 <- predict_geff_ET_RF(data)
  base_rmse_geff <- sqrt(mean((data$g_eff - pred0$geff)^2, na.rm = TRUE))
  base_rmse_ET   <- sqrt(mean((data$ET    - pred0$ET  )^2, na.rm = TRUE))
  
  purrr::map_dfr(vars, function(v) {
    set.seed(seed)
    dG <- numeric(B)
    dE <- numeric(B)
    
    for (b in seq_len(B)) {
      dp <- data
      dp[[v]] <- sample(dp[[v]])
      
      pr <- predict_geff_ET_RF(dp)
      rmseG <- sqrt(mean((dp$g_eff - pr$geff)^2, na.rm = TRUE))
      rmseE <- sqrt(mean((dp$ET    - pr$ET  )^2, na.rm = TRUE))
      
      dG[b] <- rmseG - base_rmse_geff
      dE[b] <- rmseE - base_rmse_ET
    }
    
    tibble(
      feature = v,
      mean_delta_geff = mean(dG, na.rm = TRUE),
      lo_geff = quantile(dG, 0.05, na.rm = TRUE),
      hi_geff = quantile(dG, 0.95, na.rm = TRUE),
      mean_delta_ET = mean(dE, na.rm = TRUE),
      lo_ET = quantile(dE, 0.05, na.rm = TRUE),
      hi_ET = quantile(dE, 0.95, na.rm = TRUE),
      B = B
    )
  })
}

# Use predictors that are actually present in the final RF
vars_for_perm <- intersect(c("inv_sqrt_VPD","Ta","P","Rg","CO2_term"), rf_vars_present)
perm_full <- perm_importance_geff_ET(df_opt, vars_for_perm, B = perm_B)

print(perm_full %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

p_imp_geff <- ggplot(perm_full, aes(x = reorder(feature, mean_delta_geff), y = mean_delta_geff)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_geff, ymax = hi_geff), width = 0.2) +
  coord_flip() + theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta*" RMSE_"*g[eff]*" (permute)"),
    title = expression("Permutation importance for "*g[eff]),
    subtitle = "RF-only model (VPD enters as 1/sqrt(VPD))"
  )

p_imp_ET <- ggplot(perm_full, aes(x = reorder(feature, mean_delta_ET), y = mean_delta_ET)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
  coord_flip() + theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta*" RMSE_ET (permute)"),
    title = "Permutation importance for ET",
    subtitle = "ET derived from K_ET * g_eff_RF"
  )

print(p_imp_geff + p_imp_ET)

# ------------------------------------------------------------------------------
# 8) Partial dependence plots (PDP) for g_eff and ET
#    - PDP over CO2 uses ppm on x-axis, internally converted to CO2_term
# ------------------------------------------------------------------------------

partial_dependence <- function(data, var, pred_fun, grid_size = 60) {
  stopifnot(var %in% names(data))
  
  if (var == "CO2_ppm") {
    if (!exists("CO2_1981_2005", inherits = TRUE)) stop("CO2_1981_2005 is required for CO2_ppm PDP.")
    x <- data$CO2_ppm
    xg <- seq(quantile(x, 0.02, na.rm = TRUE), quantile(x, 0.98, na.rm = TRUE), length.out = grid_size)
    
    pd <- sapply(xg, function(ppm) {
      d <- data
      d$CO2_term <- log(ppm / CO2_1981_2005)
      mean(pred_fun(d), na.rm = TRUE)
    })
  } else {
    x <- data[[var]]
    xg <- seq(quantile(x, 0.02, na.rm = TRUE), quantile(x, 0.98, na.rm = TRUE), length.out = grid_size)
    
    pd <- sapply(xg, function(v) {
      d <- data
      d[[var]] <- v
      mean(pred_fun(d), na.rm = TRUE)
    })
  }
  
  tibble(var = var, x = xg, pd = pd)
}

pred_geff_only <- function(d) predict_geff_ET_RF(d)$geff
pred_ET_only   <- function(d) predict_geff_ET_RF(d)$ET

# PDP variables: show VPD on original scale (more intuitive), but note:
# the model actually uses inv_sqrt_VPD internally.
vars_pdp <- intersect(c("VPD","Ta","P","Rg","CO2_ppm"), names(df))

# For VPD PDP we overwrite VPD and also update inv_sqrt_VPD consistently
partial_dependence_VPD <- function(data, pred_fun, grid_size = 60) {
  x <- data$VPD
  xg <- seq(quantile(x, 0.02, na.rm = TRUE), quantile(x, 0.98, na.rm = TRUE), length.out = grid_size)
  
  pd <- sapply(xg, function(vpd) {
    d <- data
    d$VPD <- pmax(vpd, eps)
    d$inv_sqrt_VPD <- 1 / sqrt(pmax(d$VPD, eps))
    mean(pred_fun(d), na.rm = TRUE)
  })
  
  tibble(var = "VPD", x = xg, pd = pd)
}

# Build PDP tables
pdp_geff_list <- list()
pdp_ET_list   <- list()

if ("VPD" %in% vars_pdp) {
  pdp_geff_list <- c(pdp_geff_list, list(partial_dependence_VPD(df, pred_geff_only, 60)))
  pdp_ET_list   <- c(pdp_ET_list,   list(partial_dependence_VPD(df, pred_ET_only,   60)))
}
for (v in setdiff(vars_pdp, "VPD")) {
  pdp_geff_list <- c(pdp_geff_list, list(partial_dependence(df, v, pred_geff_only, 60)))
  pdp_ET_list   <- c(pdp_ET_list,   list(partial_dependence(df, v, pred_ET_only,   60)))
}

pdp_geff <- bind_rows(pdp_geff_list)
pdp_ET   <- bind_rows(pdp_ET_list)

p_pdp_geff <- ggplot(pdp_geff, aes(x = x, y = pd)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ var, scales = "free_x") +
  theme_bw() +
  labs(
    x = NULL,
    y = expression(paste("Partial dependence of ", g[eff])),
    title = expression("Partial dependence plots for " * g[eff]),
    subtitle = "RF-only model; VPD effect enters through 1/sqrt(VPD)"
  )

p_pdp_ET <- ggplot(pdp_ET, aes(x = x, y = pd)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ var, scales = "free_x") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Partial dependence of ET",
    title = "Partial dependence plots for ET",
    subtitle = "ET derived from K_ET(VPD,Ta,...) × g_eff_RF"
  )

print(p_pdp_geff)
print(p_pdp_ET)

# ------------------------------------------------------------------------------
# 9) (Optional) Quick RF impurity importance for g_eff model (unitless)
# ------------------------------------------------------------------------------

if (!is.null(rf_geff$variable.importance)) {
  vi <- tibble(feature = names(rf_geff$variable.importance),
               importance = as.numeric(rf_geff$variable.importance)) %>%
    arrange(desc(importance))
  
  p_vip <- ggplot(vi, aes(x = reorder(feature, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = NULL, y = "Impurity importance (unitless)",
         title = "RF (g_eff) impurity importance (quick ranking)")
  print(p_vip)
}


