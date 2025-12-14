################################################################################
## Hybrid g_eff + ET model with physics constraint, CO2 regime gate, importance,
## and partial dependence plots.
##
## OVERVIEW (what this script does)
##  1) Fit a parametric stomatal/“effective conductance” baseline:
##         g_eff ≈ b0 + b1 / sqrt(VPD)
##     This enforces the desired 1/sqrt(VPD) scaling (Jarvis-like behavior).
##
##  2) Fit a Random Forest to the *residuals* of g_eff after removing the VPD law:
##         resid_geff = g_eff - g_param(VPD)
##     using Ta, P, Rg, and (optionally) CO2_term as predictors.
##
##  3) Predict full g_eff and ET using:
##         g_eff_hat = g_param(VPD) + RF_residual(Ta,P,Rg,CO2_term)
##         ET_hat    = K_ET(VPD,Ta) * g_eff_hat
##
##  4) Apply a CO2 “regime monotonicity” gate:
##     we only keep CO2_term if increasing CO2 (0 → CMIP5 → CMIP6) produces
##     a monotonic *decrease* in predicted g_eff for most samples.
##
##  5) Produce diagnostics (RMSE, R²), observed-vs-predicted plots,
##     permutation importance for the FULL hybrid model (incl. VPD),
##     and partial dependence plots (PDP) for g_eff and ET.
##
## REQUIREMENTS (assumed to exist)
##   - df_opt : training frame
##   - df     : prediction/evaluation frame
##   - recompute_K_ET(df) : function returning K_ET per row (units consistent with ET)
##   - tune_rf_oob(), fit_rf_fixed() : your RF tuning helpers (as in your earlier code)
##   - CO2_1981_2005, CO2_2076_2100_RCP85_CMIP5, CO2_2076_2100_RCP85_CMIP6 : scalars
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(DEoptim)
  library(ranger)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 0) Global knobs & basic sanity checks
# ------------------------------------------------------------------------------

# "eps" prevents division-by-zero and sqrt-of-zero issues (important for VPD)
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# Number of repeats for permutation importance (higher = smoother CIs, slower)
perm_B <- if (exists("perm_B", inherits = TRUE)) perm_B else 200L

# Random Forest training knobs (do not overwrite if already set in your session)
rf_num_trees <- if (exists("rf_num_trees", inherits = TRUE)) rf_num_trees else 1200L
rf_tune_mode <- if (exists("rf_tune_mode", inherits = TRUE)) rf_tune_mode else "full"

# We require these columns for both training and prediction.
stopifnot(all(c("VPD","Ta","P","Rg","CO2_term","ET","g_eff") %in% names(df_opt)))
stopifnot(all(c("VPD","Ta","P","Rg","CO2_term","ET","g_eff") %in% names(df)))

# Guard: VPD must be positive for 1/sqrt(VPD); clip tiny/negative values
df_opt <- df_opt %>% mutate(VPD = pmax(VPD, eps))
df     <- df     %>% mutate(VPD = pmax(VPD, eps))

# ------------------------------------------------------------------------------
# Add real CO2 concentration (ppm) for interpretation & PDPs ONLY
# Model fitting continues to use CO2_term
# ------------------------------------------------------------------------------

df_opt <- df_opt %>%
  mutate(
    CO2_ppm = case_when(
      # EUR-44: no physiology → represent on ppm axis as reference CO2
      grepl("EUR", ensemble, ignore.case = TRUE) ~ CO2_1981_2005,
      
      # Historical reference period
      PERIOD == "1981_2005" ~ CO2_1981_2005,
      
      # Future scenarios
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

# CO2 gate knobs:
# - majority_frac: fraction of samples that must show monotone negative response
# - tol: small negative tolerance (<= tol counts as "negative enough")
co2_neg_majority_frac <- 0.60
co2_neg_total_tol     <- -1e-6

# ------------------------------------------------------------------------------
# 1) Parametric VPD-only baseline for g_eff
#    g_eff ≈ b0 + b1 / sqrt(VPD)
# ------------------------------------------------------------------------------

# Objective = RMSE(g_eff_obs - g_eff_hat) for given (b0,b1).
# This enforces the "g_eff scales with 1/sqrt(VPD)" structure.
obj_fun_geff_vpd <- function(par, input_df) {
  b0 <- par[[1]]
  b1 <- par[[2]]
  inv_sqrt_vpd <- 1 / sqrt(pmax(input_df$VPD, eps))
  g_hat <- b0 + b1 * inv_sqrt_vpd
  
  # Optional: add penalty if you want to constrain g_hat >= 0.
  # pen <- 1e6 * mean(pmax(-g_hat, 0)^2, na.rm = TRUE)
  pen <- 0
  
  sqrt(mean((input_df$g_eff - g_hat)^2, na.rm = TRUE)) + pen
}

# Optimize (b0,b1) using DEoptim (global optimizer; robust for nonconvex surfaces)
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

# Store fitted parameters (for readability later)
b_hat <- de_geff$optim$bestmem
names(b_hat) <- c("b0_VPD", "b1_VPD")
message(sprintf("Parametric g_eff fit: b0=%.6f, b1=%.6f", b_hat[1], b_hat[2]))

# Compute parametric baseline g_eff for training (df_opt) and evaluation (df)
inv_sqrt_VPD_opt <- 1 / sqrt(df_opt$VPD)
inv_sqrt_VPD     <- 1 / sqrt(df$VPD)

g_param_opt  <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD_opt
g_param_full <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD

# Convert baseline g_eff into baseline ET using your physics factor K_ET(VPD, Ta).
# This is a diagnostic baseline (Stage 1), not the final ET prediction.
K_ET_opt  <- recompute_K_ET(df_opt)
K_ET_full <- recompute_K_ET(df)

ET_vpd_opt <- K_ET_opt  * g_param_opt
ET_vpd_hat <- K_ET_full * g_param_full

# ------------------------------------------------------------------------------
# 2) Random Forest on residual g_eff (no VPD in RF predictors)
#    resid_geff = g_eff - g_param(VPD)
# ------------------------------------------------------------------------------

# RF predictors: exclude VPD intentionally (VPD effect already captured by parametric term)
rf_vars_full <- c("Ta", "P", "Rg", "CO2_term")

# Training table for RF: target + candidate predictors; filter invalid numeric rows
df_rf_geff <- df_opt %>%
  mutate(resid_geff = g_eff - g_param_opt) %>%
  select(resid_geff, all_of(rf_vars_full)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

target_col_geff <- "resid_geff"

# Train RF (tuned or fixed, depending on rf_tune_mode)
set.seed(if (exists("rng_seed")) rng_seed else 123)

if (rf_tune_mode == "none") {
  # Fixed/default RF settings via your helper
  rf_geff <- fit_rf_fixed(
    df_rf_geff,
    target_col = target_col_geff,
    num.trees  = rf_num_trees,
    seed       = if (exists("rng_seed")) rng_seed else 123
  )
  
  # Keep a record of params (useful for reproducibility/logging)
  best_params_geff <- list(
    mtry            = rf_geff$mtry,
    min.node.size   = rf_geff$min.node.size,
    sample.fraction = 1.0,
    splitrule       = "variance",
    max.depth       = 0L
  )
  
} else {
  # Hyperparameter tuning via your OOB tuner helper
  if (rf_tune_mode == "quick") {
    tuned_geff <- tune_rf_oob(
      df_rf_geff,
      target_col = target_col_geff,
      num.trees       = rf_num_trees,
      min.node.size   = c(3L, 5L, 10L),
      sample.fraction = c(0.8, 1.0),
      splitrule       = c("variance", "extratrees"),
      max.depth       = c(0L, 12L)
    )
  } else {
    tuned_geff <- tune_rf_oob(
      df_rf_geff,
      target_col = target_col_geff,
      num.trees  = rf_num_trees
    )
  }
  
  rf_geff          <- tuned_geff$best_fit
  best_params_geff <- tuned_geff$best_params
}

# RF out-of-bag error (ranger stores MSE as prediction.error for regression)
base_oob_rmse_geff <- sqrt(rf_geff$prediction.error)
rf_vars_present_geff <- setdiff(names(df_rf_geff), target_col_geff)

message(sprintf("RF_g_eff (residual) OOB RMSE = %.3f", base_oob_rmse_geff))
message(sprintf(
  "RF_g_eff params: mtry=%s, min.node.size=%s, sample.fraction=%.2f, splitrule=%s, max.depth=%s",
  best_params_geff$mtry, best_params_geff$min.node.size,
  best_params_geff$sample.fraction,
  as.character(best_params_geff$splitrule),
  best_params_geff$max.depth
))

# Predict residual g_eff on new data using the currently active RF predictors.
# If predictors are missing, return 0 (meaning "no residual correction").
.predict_rf_resid_geff <- function(newdata) {
  vars <- intersect(rf_vars_present_geff, colnames(newdata))
  if (!length(vars)) return(rep(0, nrow(newdata)))
  as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
}

# ------------------------------------------------------------------------------
# 3) Full hybrid prediction function for g_eff and ET
# ------------------------------------------------------------------------------

predict_geff_ET_full <- function(newdata) {
  # 1) physics/constrained baseline: g_param(VPD)
  vVPD <- pmax(newdata$VPD, eps)
  inv_sqrt_vpd <- 1 / sqrt(vVPD)
  g_param <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_vpd
  
  # 2) meteorological conversion factor (ET demand) from your helper
  K_ET <- recompute_K_ET(newdata)
  
  # 3) data-driven residual correction: RF(Ta, P, Rg, CO2_term)
  g_resid <- .predict_rf_resid_geff(newdata)
  
  # 4) recombine
  g_eff_hat <- g_param + g_resid
  ET_hat    <- K_ET * g_eff_hat
  
  list(geff = as.numeric(g_eff_hat), ET = as.numeric(ET_hat))
}

# ------------------------------------------------------------------------------
# 4) CO2 regime monotonicity gate (physics sanity check)
#    Keep CO2_term only if increasing CO2 causes monotone *decrease* in g_eff.
#
#    Why we do this:
#      - With climate model ensembles and time slices, CO2 can be confounded with
#        other shifts. RF may exploit spurious correlation.
#      - We enforce a sign-consistent physiological prior: higher CO2 ⇒ lower g_eff.
#
#    How this gate works:
#      - Hold climate fixed (same rows), change only CO2_term between regimes:
#           0 (EUR-44) → CMIP5 → CMIP6
#      - Check that g_eff decreases in most samples and on average.
# ------------------------------------------------------------------------------

co2_regime_gate <- function(data, var = "CO2_term",
                            co2_levels,
                            majority_frac = 0.60,
                            tol = -1e-6) {
  stopifnot(var %in% names(data))
  stopifnot(length(co2_levels) >= 2)
  
  # Predicted g_eff under each CO2 regime, with all other columns unchanged
  preds <- lapply(co2_levels, function(val) {
    d <- data
    d[[var]] <- val
    predict_geff_ET_full(d)$geff
  })
  
  # Consecutive differences: should be <= tol for monotone non-increasing response
  deltas <- Map(function(a, b) b - a, preds[-length(preds)], preds[-1])
  
  # Share of comparisons that are "negative enough"
  frac_neg <- mean(unlist(lapply(deltas, function(d) d <= tol)), na.rm = TRUE)
  
  # Mean overall change from lowest to highest CO2 regime
  mean_total <- mean(preds[[length(preds)]] - preds[[1]], na.rm = TRUE)
  
  list(
    frac_neg   = frac_neg,
    mean_total = mean_total,
    pass       = (frac_neg >= majority_frac && mean_total <= tol)
  )
}

# Only attempt the gate if CO2_term currently exists in the RF predictor set
if ("CO2_term" %in% rf_vars_present_geff) {
  
  # CO2 regime constants must exist in the environment for this test
  if (!exists("CO2_1981_2005", inherits = TRUE) ||
      !exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) ||
      !exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
    
    warning("CO2 regime values not found in environment; skipping CO2 sign gate.")
    
  } else {
    
    # Use *actual* CO2 regimes:
    #   - 0: EUR-44 (RCMs ignore CO2 physiology)
    #   - log(CO2_future/CO2_hist): CMIP5 / CMIP6 end-century forcing
    co2_levels <- c(
      0,
      log(CO2_2076_2100_RCP85_CMIP5 / CO2_1981_2005),
      log(CO2_2076_2100_RCP85_CMIP6 / CO2_1981_2005)
    )
    
    # Gate evaluated on training data (df_opt) to avoid "peeking" at evaluation set
    gchk <- co2_regime_gate(
      df_opt,
      var           = "CO2_term",
      co2_levels    = co2_levels,
      majority_frac = co2_neg_majority_frac,
      tol           = co2_neg_total_tol
    )
    
    message(sprintf(
      "CO2 regime gate: frac(monotone negative)=%.1f%%, mean Δg_eff(hi-low)=%.4g",
      100 * gchk$frac_neg, gchk$mean_total
    ))
    
    # If the model implies non-physiological CO2 behavior, drop CO2_term and refit RF
    if (!isTRUE(gchk$pass)) {
      
      message("Decision: CO2_term fails regime monotonicity gate → dropping CO2_term and refitting RF.")
      
      rf_vars_full2 <- setdiff(rf_vars_full, "CO2_term")
      
      df_rf_geff2 <- df_opt %>%
        mutate(resid_geff = g_eff - g_param_opt) %>%
        select(resid_geff, all_of(rf_vars_full2)) %>%
        filter(if_all(where(is.numeric), ~ is.finite(.x)))
      
      set.seed(if (exists("rng_seed")) rng_seed else 123)
      
      # Refit RF using the same tuning strategy as before
      if (rf_tune_mode == "none") {
        rf_geff <- fit_rf_fixed(
          df_rf_geff2,
          target_col = "resid_geff",
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
          tuned_geff2 <- tune_rf_oob(
            df_rf_geff2,
            "resid_geff",
            num.trees       = rf_num_trees,
            min.node.size   = c(3L,5L,10L),
            sample.fraction = c(0.8,1.0),
            splitrule       = c("variance","extratrees"),
            max.depth       = c(0L,12L)
          )
        } else {
          tuned_geff2 <- tune_rf_oob(
            df_rf_geff2,
            "resid_geff",
            num.trees = rf_num_trees
          )
        }
        rf_geff          <- tuned_geff2$best_fit
        best_params_geff <- tuned_geff2$best_params
      }
      
      # Update tracking vars after refit
      base_oob_rmse_geff <- sqrt(rf_geff$prediction.error)
      rf_vars_present_geff <- setdiff(names(df_rf_geff2), "resid_geff")
      
      message(sprintf("Refit RF_g_eff (without CO2_term) OOB RMSE = %.3f", base_oob_rmse_geff))
      
      # IMPORTANT: update predictor helper so the rest of the script uses the new RF vars
      .predict_rf_resid_geff <- function(newdata) {
        vars <- intersect(rf_vars_present_geff, colnames(newdata))
        if (!length(vars)) return(rep(0, nrow(newdata)))
        as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
      }
      
    } else {
      message("Decision: CO2_term passes regime monotonicity gate → keeping CO2_term.")
    }
  }
}

# ------------------------------------------------------------------------------
# 5) Final predictions on df (evaluation / full dataset)
# ------------------------------------------------------------------------------

pred_full <- predict_geff_ET_full(df)

output_df_geff <- df %>%
  mutate(
    # Stage-1 (parametric only)
    g_param     = g_param_full,
    ET_vpd_hat  = ET_vpd_hat,
    
    # Full hybrid prediction (parametric + RF residuals)
    g_eff_final = pred_full$geff,
    ET_final    = pred_full$ET
  )

# ------------------------------------------------------------------------------
# 6) Diagnostics (Stage1 vs Full): RMSE & R²
# ------------------------------------------------------------------------------

RMSE_param_geff <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_param)^2, na.rm = TRUE))
R2_param_geff   <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_param, use = "complete.obs")^2)

RMSE_final_geff <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_eff_final)^2, na.rm = TRUE))
R2_final_geff   <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_eff_final, use = "complete.obs")^2)

RMSE_param_ET   <- sqrt(mean((output_df_geff$ET - output_df_geff$ET_vpd_hat)^2, na.rm = TRUE))
R2_param_ET     <- suppressWarnings(cor(output_df_geff$ET, output_df_geff$ET_vpd_hat, use = "complete.obs")^2)

RMSE_final_ET   <- sqrt(mean((output_df_geff$ET - output_df_geff$ET_final)^2, na.rm = TRUE))
R2_final_ET     <- suppressWarnings(cor(output_df_geff$ET, output_df_geff$ET_final, use = "complete.obs")^2)

message(sprintf(
  "g_eff — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Full: RMSE=%.3f, R^2=%.3f",
  RMSE_param_geff, R2_param_geff, RMSE_final_geff, R2_final_geff
))
message(sprintf(
  "ET    — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Full: RMSE=%.3f, R^2=%.3f",
  RMSE_param_ET, R2_param_ET, RMSE_final_ET, R2_final_ET
))

# ------------------------------------------------------------------------------
# 7) Variable importance of the *residual RF* (excludes VPD by construction)
#    NOTE: ranger impurity importance is unitless and can be biased; we mainly
#    use it as a quick “within-RF” ranking of residual drivers.
# ------------------------------------------------------------------------------

p_vip_geff <- NULL
if (!is.null(rf_geff$variable.importance)) {
  
  vi_geff <- tibble(
    feature    = names(rf_geff$variable.importance),
    importance = as.numeric(rf_geff$variable.importance)
  ) %>% arrange(desc(importance))
  
  p_vip_geff <- ggplot(vi_geff, aes(x = reorder(feature, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      x = NULL,
      y = "RF importance (impurity; unitless)",
      title = expression("g"[eff] * " residual RF: Variable Importance"),
      subtitle = "Explains residual g_eff beyond the 1/sqrt(VPD) parametric term"
    )
  
  print(p_vip_geff)
}

# ------------------------------------------------------------------------------
# 8) Observed vs predicted (full hybrid model)
# ------------------------------------------------------------------------------

p_obs_pred_geff_full <- output_df_geff %>%
  transmute(g_eff_obs = g_eff, geff_pred = g_eff_final) %>%
  filter(is.finite(g_eff_obs), is.finite(geff_pred)) %>%
  ggplot(aes(g_eff_obs, geff_pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff], " (hybrid)")),
    title = expression(paste("Observed vs Predicted ", g[eff], " (hybrid model)")),
    subtitle = sprintf("Stage1: RMSE=%.3f R^2=%.3f   |   Full: RMSE=%.3f R^2=%.3f",
                       RMSE_param_geff, R2_param_geff, RMSE_final_geff, R2_final_geff)
  )

p_obs_pred_ET_full <- output_df_geff %>%
  transmute(ET_obs = ET, ET_pred = ET_final) %>%
  filter(is.finite(ET_obs), is.finite(ET_pred)) %>%
  ggplot(aes(ET_obs, ET_pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Observed ET",
    y = "Predicted ET (hybrid)",
    title = "Observed vs Predicted ET (hybrid model)",
    subtitle = sprintf("Stage1: RMSE=%.3f R^2=%.3f   |   Full: RMSE=%.3f R^2=%.3f",
                       RMSE_param_ET, R2_param_ET, RMSE_final_ET, R2_final_ET)
  )

print(p_obs_pred_geff_full + p_obs_pred_ET_full)

# ------------------------------------------------------------------------------
# 9) Full-model permutation importance (includes VPD)
#    Why permutation importance?
#      - model-agnostic
#      - answers: “how much does prediction skill degrade if X’s information is destroyed?”
#      - uses ΔRMSE, which is interpretable in target units.
# ------------------------------------------------------------------------------

perm_importance_full_model <- function(
    data,
    vars = c("VPD", "Ta", "P", "Rg", "CO2_term"),
    B = 100L,
    seed = 2024
) {
  vars <- intersect(vars, names(data))
  if (!length(vars)) stop("None of the requested vars found in data.")
  
  # Baseline prediction errors using the current hybrid model
  pred0 <- predict_geff_ET_full(data)
  base_rmse_ET   <- sqrt(mean((data$ET    - pred0$ET  )^2, na.rm = TRUE))
  base_rmse_geff <- sqrt(mean((data$g_eff - pred0$geff)^2, na.rm = TRUE))
  
  purrr::map_dfr(vars, function(v) {
    set.seed(seed)
    dET   <- numeric(B)
    dGeff <- numeric(B)
    
    for (b in seq_len(B)) {
      # Permute ONE predictor across rows: breaks its relationship to the target
      dfp <- data
      dfp[[v]] <- sample(dfp[[v]])
      
      pred_p <- predict_geff_ET_full(dfp)
      
      rmse_ET_p   <- sqrt(mean((dfp$ET    - pred_p$ET  )^2, na.rm = TRUE))
      rmse_geff_p <- sqrt(mean((dfp$g_eff - pred_p$geff)^2, na.rm = TRUE))
      
      # Importance = increase in RMSE caused by destroying predictor v
      dET[b]   <- rmse_ET_p   - base_rmse_ET
      dGeff[b] <- rmse_geff_p - base_rmse_geff
    }
    
    tibble(
      feature        = v,
      
      mean_delta_ET  = mean(dET,   na.rm = TRUE),
      lo_ET          = quantile(dET,   0.05, na.rm = TRUE),
      hi_ET          = quantile(dET,   0.95, na.rm = TRUE),
      rel_ET_pct     = 100 * mean(dET,   na.rm = TRUE) / base_rmse_ET,
      
      mean_delta_geff = mean(dGeff, na.rm = TRUE),
      lo_geff         = quantile(dGeff, 0.05, na.rm = TRUE),
      hi_geff         = quantile(dGeff, 0.95, na.rm = TRUE),
      rel_geff_pct    = 100 * mean(dGeff, na.rm = TRUE) / base_rmse_geff,
      
      B = B
    )
  }) %>%
    arrange(desc(mean_delta_geff))
}

perm_full <- perm_importance_full_model(
  df_opt,
  vars = c("VPD", "Ta", "P", "Rg", "CO2_term"),
  B = perm_B
)

print(perm_full %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# ------------------------------------------------------------------------------
# 10) Importance barplots
# ------------------------------------------------------------------------------

p_imp_full_geff <- ggplot(perm_full, aes(x = reorder(feature, mean_delta_geff), y = mean_delta_geff)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_geff, ymax = hi_geff), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta * " RMSE_" * g[eff] * " (permute)"),
    title = expression("Full-model variable importance for " * g[eff]),
    subtitle = "Hybrid model (1/sqrt(VPD) parametric + RF residuals; includes VPD)"
  )

p_imp_full_ET <- ggplot(perm_full, aes(x = reorder(feature, mean_delta_ET), y = mean_delta_ET)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta * " RMSE_ET (permute)"),
    title = "Full-model variable importance for ET",
    subtitle = "Hybrid model (1/sqrt(VPD) parametric + RF residuals; includes VPD)"
  )

# Safe printing: p_vip_geff can be NULL if ranger didn't return importance
if (!is.null(p_vip_geff)) {
  print(p_vip_geff + p_imp_full_geff)
} else {
  print(p_imp_full_geff)
}
print(p_imp_full_geff + p_imp_full_ET)

# ------------------------------------------------------------------------------
# 11) Partial dependence plots (PDP)
#    PDP shows the *average* response of the model when a variable is swept
#    over a grid, while all other variables are kept as observed.
#
#    IMPORTANT CAVEAT:
#      PDP assumes predictors can vary independently; with correlated climate
#      variables and regime structure, interpret PDP shapes qualitatively.
# ------------------------------------------------------------------------------

partial_dependence <- function(data, var, pred_fun, grid_size = 50) {
  stopifnot(var %in% names(data))
  
  # Special handling for real CO2 (ppm):
  # We sweep CO2_ppm on the x-axis, but internally convert to CO2_term
  if (var == "CO2_ppm") {
    
    xg <- seq(
      quantile(data$CO2_ppm, 0.02, na.rm = TRUE),
      quantile(data$CO2_ppm, 0.98, na.rm = TRUE),
      length.out = grid_size
    )
    
    pd <- sapply(xg, function(ppm) {
      d <- data
      d$CO2_term <- log(ppm / CO2_1981_2005)  # convert ppm → log forcing
      mean(pred_fun(d), na.rm = TRUE)
    })
    
  } else {
    
    x <- data[[var]]
    xg <- seq(
      quantile(x, 0.02, na.rm = TRUE),
      quantile(x, 0.98, na.rm = TRUE),
      length.out = grid_size
    )
    
    pd <- sapply(xg, function(v) {
      d <- data
      d[[var]] <- v
      mean(pred_fun(d), na.rm = TRUE)
    })
  }
  
  tibble(var = var, x = xg, pd = pd)
}

# Convenience wrappers to reuse your hybrid predictor for a single target
pred_geff_only <- function(d) predict_geff_ET_full(d)$geff
pred_ET_only   <- function(d) predict_geff_ET_full(d)$ET

# Only compute PDPs for variables that exist (in case CO2_term got dropped upstream)
vars_pdp <- intersect(c("VPD", "Ta", "P", "Rg", "CO2_ppm"), names(df))

# PDP for g_eff
pdp_geff <- bind_rows(lapply(vars_pdp, partial_dependence,
                             data = df, pred_fun = pred_geff_only, grid_size = 60))

p_pdp_geff <- ggplot(pdp_geff, aes(x = x, y = pd)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ var, scales = "free_x") +
  theme_bw() +
  labs(
    x = NULL,
    y = expression(paste("Partial dependence of ", g[eff])),
    title = expression("Partial dependence plots for " * g[eff]),
    subtitle = "Hybrid model: parametric VPD + RF residuals"
  )
print(p_pdp_geff)

# PDP for ET
pdp_ET <- bind_rows(lapply(vars_pdp, partial_dependence,
                           data = df, pred_fun = pred_ET_only, grid_size = 60))

p_pdp_ET <- ggplot(pdp_ET, aes(x = x, y = pd)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ var, scales = "free_x") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Partial dependence of ET",
    title = "Partial dependence plots for ET",
    subtitle = "Hybrid model"
  )
print(p_pdp_ET)
