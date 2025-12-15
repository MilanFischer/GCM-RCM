################################################################################
## Hybrid g_eff + ET model with physics constraint + CO2 regime gate (FIT ONLY)
## - Keeps informative messages (CO2 gate debug, RF params, diagnostics)
## - NO plots / NO permutation / NO PDP / NO ALE (moved to separate scripts)
## - Saves a UNIVERSAL bundle for downstream scripts
##
## REQUIREMENTS (assumed to exist)
##   - df_opt : training frame
##   - df     : prediction/evaluation frame
##   - recompute_K_ET(df) : function returning K_ET per row (units consistent with ET)
##   - tune_rf_oob(), fit_rf_fixed() : your RF tuning helpers
##   - CO2_1981_2005, CO2_2076_2100_RCP85_CMIP5, CO2_2076_2100_RCP85_CMIP6 : scalars
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DEoptim)
  library(ranger)
})

source("./src/Jarvis&RF_preprocessing.R")

# ---- output ----
out_file <- "./RData/20251214_RF_objects.RData"

# --------------------------------------------------------------------------
# 0) Global knobs & basic sanity checks (keep behavior)
# --------------------------------------------------------------------------

eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

rf_num_trees <- if (exists("rf_num_trees", inherits = TRUE)) rf_num_trees else 1200L
rf_tune_mode <- if (exists("rf_tune_mode", inherits = TRUE)) rf_tune_mode else "full"
rng_seed     <- if (exists("rng_seed", inherits = TRUE)) rng_seed else 1990L

stopifnot(exists("df_opt", inherits = TRUE))
stopifnot(exists("df", inherits = TRUE))
stopifnot(exists("recompute_K_ET", inherits = TRUE))
stopifnot(is.function(recompute_K_ET))

stopifnot(exists("tune_rf_oob", inherits = TRUE))
stopifnot(exists("fit_rf_fixed", inherits = TRUE))

stopifnot(all(c("VPD","Ta","P","Rg","CO2_term","ET","g_eff") %in% names(df_opt)))
stopifnot(all(c("VPD","Ta","P","Rg","CO2_term","ET","g_eff") %in% names(df)))

# Guard: VPD must be positive
df_opt <- df_opt |> mutate(VPD = pmax(VPD, eps))
df     <- df     |> mutate(VPD = pmax(VPD, eps))

# Guard required metadata columns for CO2_ppm mapping (optional, for later scripts)
.has_co2_meta_opt <- all(c("ensemble", "PERIOD") %in% names(df_opt))
.has_co2_meta_df  <- all(c("ensemble", "PERIOD") %in% names(df))

if (!.has_co2_meta_opt) message("Note: df_opt lacks ensemble/PERIOD → CO2_ppm will not be created for df_opt.")
if (!.has_co2_meta_df)  message("Note: df lacks ensemble/PERIOD → CO2_ppm will not be created for df.")

# Optional: store CO2_ppm for interpretation (won't affect fitting)
if (.has_co2_meta_opt &&
    exists("CO2_1981_2005", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
  
  df_opt <- df_opt |>
    mutate(
      CO2_ppm = case_when(
        grepl("EUR", ensemble, ignore.case = TRUE) ~ CO2_1981_2005,
        PERIOD == "1981_2005" ~ CO2_1981_2005,
        PERIOD == "2076_2100" & ensemble == "CMIP5" ~ CO2_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ensemble == "CMIP6" ~ CO2_2076_2100_RCP85_CMIP6,
        TRUE ~ NA_real_
      )
    )
}

if (.has_co2_meta_df &&
    exists("CO2_1981_2005", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) &&
    exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
  
  df <- df |>
    mutate(
      CO2_ppm = case_when(
        grepl("EUR", ensemble, ignore.case = TRUE) ~ CO2_1981_2005,
        PERIOD == "1981_2005" ~ CO2_1981_2005,
        PERIOD == "2076_2100" & ensemble == "CMIP5" ~ CO2_2076_2100_RCP85_CMIP5,
        PERIOD == "2076_2100" & ensemble == "CMIP6" ~ CO2_2076_2100_RCP85_CMIP6,
        TRUE ~ NA_real_
      )
    )
}

# CO2 gate knobs
co2_neg_majority_frac <- 0.60
co2_neg_total_tol     <- -1e-6

# --------------------------------------------------------------------------
# 1) Parametric VPD-only baseline for g_eff: b0 + b1 / sqrt(VPD)
# --------------------------------------------------------------------------

obj_fun_geff_vpd <- function(par, input_df) {
  b0 <- par[[1]]
  b1 <- par[[2]]
  inv_sqrt_vpd <- 1 / sqrt(pmax(input_df$VPD, eps))
  g_hat <- b0 + b1 * inv_sqrt_vpd
  pen <- 0
  sqrt(mean((input_df$g_eff - g_hat)^2, na.rm = TRUE)) + pen
}

set.seed(rng_seed)
de_geff <- DEoptim::DEoptim(
  fn       = obj_fun_geff_vpd,
  lower    = c(-10, -1e3),
  upper    = c( 10,  1e3),
  input_df = df_opt,
  control  = DEoptim::DEoptim.control(
    NP = 60L, itermax = 2000L, reltol = 1e-8,
    trace = if (exists("trace_flag", inherits = TRUE)) trace_flag else FALSE
  )
)

b_hat <- de_geff$optim$bestmem
names(b_hat) <- c("b0_VPD", "b1_VPD")
message(sprintf("Parametric g_eff fit: b0=%.6f, b1=%.6f", b_hat[1], b_hat[2]))

g_param_opt  <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * (1 / sqrt(df_opt$VPD))
g_param_full <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * (1 / sqrt(df$VPD))

K_ET_opt  <- recompute_K_ET(df_opt)
K_ET_full <- recompute_K_ET(df)

ET_vpd_opt <- K_ET_opt  * g_param_opt
ET_vpd_hat <- K_ET_full * g_param_full

# --------------------------------------------------------------------------
# 2) Random Forest on residual g_eff (no VPD in RF predictors)
# --------------------------------------------------------------------------

rf_vars_full <- c("Ta", "P", "Rg", "CO2_term")

df_rf_geff <- df_opt |>
  mutate(resid_geff = g_eff - g_param_opt) |>
  select(resid_geff, all_of(rf_vars_full)) |>
  filter(if_all(where(is.numeric), is.finite))

target_col_geff <- "resid_geff"

set.seed(rng_seed)

if (rf_tune_mode == "none") {
  
  rf_geff <- fit_rf_fixed(
    df_rf_geff,
    target_col = target_col_geff,
    num.trees  = rf_num_trees,
    seed       = rng_seed
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

.predict_rf_resid_geff <- function(newdata) {
  vars <- intersect(rf_vars_present_geff, names(newdata))
  if (!length(vars)) return(rep(0, nrow(newdata)))
  as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
}

# --------------------------------------------------------------------------
# 3) Full hybrid predictor
# --------------------------------------------------------------------------

predict_geff_ET_full <- function(newdata) {
  vVPD <- pmax(newdata$VPD, eps)
  inv_sqrt_vpd <- 1 / sqrt(vVPD)
  
  g_param <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_vpd
  K_ET    <- recompute_K_ET(newdata)
  g_resid <- .predict_rf_resid_geff(newdata)
  
  g_eff_hat <- g_param + g_resid
  ET_hat    <- K_ET * g_eff_hat
  
  list(geff = as.numeric(g_eff_hat), ET = as.numeric(ET_hat))
}

# --------------------------------------------------------------------------
# 4) CO2 regime monotonicity gate (messages preserved)
# --------------------------------------------------------------------------

co2_gate_info <- list(
  majority_frac = co2_neg_majority_frac,
  tol           = co2_neg_total_tol,
  passed        = NA,
  frac_negative = NA,
  mean_delta    = NA
)

co2_regime_gate <- function(data, var = "CO2_term",
                            co2_levels,
                            majority_frac = 0.60,
                            tol = -1e-6) {
  stopifnot(var %in% names(data))
  stopifnot(length(co2_levels) >= 2)
  
  preds <- lapply(co2_levels, function(val) {
    d <- data
    d[[var]] <- val
    predict_geff_ET_full(d)$geff
  })
  
  g0 <- preds[[1]]
  deltas_vs0 <- lapply(preds[-1], function(gk) gk - g0)
  
  frac_neg <- mean(unlist(lapply(deltas_vs0, function(d) d <= tol)), na.rm = TRUE)
  mean_total <- mean(unlist(deltas_vs0), na.rm = TRUE)
  
  list(
    frac_neg   = frac_neg,
    mean_total = mean_total,
    pass       = (frac_neg >= majority_frac && mean_total <= tol)
  )
}

if ("CO2_term" %in% rf_vars_present_geff) {
  
  if (!exists("CO2_1981_2005", inherits = TRUE) ||
      !exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) ||
      !exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
    
    warning("CO2 regime values not found in environment; skipping CO2 sign gate.")
    
  } else {
    
    co2_levels <- c(
      0,
      log(CO2_2076_2100_RCP85_CMIP5 / CO2_1981_2005),
      log(CO2_2076_2100_RCP85_CMIP6 / CO2_1981_2005)
    )
    
    gchk <- co2_regime_gate(
      df_opt,
      var           = "CO2_term",
      co2_levels    = co2_levels,
      majority_frac = co2_neg_majority_frac,
      tol           = co2_neg_total_tol
    )
    
    message(sprintf(
      "CO2 gate: frac(Δg_eff<=tol)=%.1f%%, mean Δg_eff(vs 0)=%.4g",
      100 * gchk$frac_neg, gchk$mean_total
    ))
    
    preds_dbg <- lapply(co2_levels, function(val) {
      d <- df_opt
      d$CO2_term <- val
      predict_geff_ET_full(d)$geff
    })
    
    d5 <- preds_dbg[[2]] - preds_dbg[[1]]
    d6 <- preds_dbg[[3]] - preds_dbg[[1]]
    
    message(sprintf("CO2 debug: CMIP5 vs 0  frac(d<=tol)=%.1f%%  mean(d)=%.4g",
                    100 * mean(d5 <= co2_neg_total_tol, na.rm=TRUE),
                    mean(d5, na.rm=TRUE)))
    message(sprintf("CO2 debug: CMIP6 vs 0  frac(d<=tol)=%.1f%%  mean(d)=%.4g",
                    100 * mean(d6 <= co2_neg_total_tol, na.rm=TRUE),
                    mean(d6, na.rm=TRUE)))
    
    print(stats::quantile(d5, c(.05,.25,.5,.75,.95), na.rm=TRUE))
    print(stats::quantile(d6, c(.05,.25,.5,.75,.95), na.rm=TRUE))
    
    co2_gate_info$passed        <- gchk$pass
    co2_gate_info$frac_negative <- gchk$frac_neg
    co2_gate_info$mean_delta    <- gchk$mean_total
    
    if (!isTRUE(gchk$pass)) {
      
      message("Decision: CO2_term fails regime monotonicity gate → dropping CO2_term and refitting RF.")
      
      rf_vars_full2 <- setdiff(rf_vars_full, "CO2_term")
      
      df_rf_geff2 <- df_opt |>
        mutate(resid_geff = g_eff - g_param_opt) |>
        select(resid_geff, all_of(rf_vars_full2)) |>
        filter(if_all(where(is.numeric), is.finite))
      
      set.seed(rng_seed)
      
      if (rf_tune_mode == "none") {
        rf_geff <- fit_rf_fixed(
          df_rf_geff2,
          target_col = "resid_geff",
          num.trees  = rf_num_trees,
          seed       = rng_seed
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
            target_col = "resid_geff",
            num.trees       = rf_num_trees,
            min.node.size   = c(3L,5L,10L),
            sample.fraction = c(0.8,1.0),
            splitrule       = c("variance","extratrees"),
            max.depth       = c(0L,12L)
          )
        } else {
          tuned_geff2 <- tune_rf_oob(
            df_rf_geff2,
            target_col = "resid_geff",
            num.trees  = rf_num_trees
          )
        }
        rf_geff          <- tuned_geff2$best_fit
        best_params_geff <- tuned_geff2$best_params
      }
      
      base_oob_rmse_geff <- sqrt(rf_geff$prediction.error)
      rf_vars_present_geff <- setdiff(names(df_rf_geff2), "resid_geff")
      
      message(sprintf("Refit RF_g_eff (without CO2_term) OOB RMSE = %.3f", base_oob_rmse_geff))
      
      .predict_rf_resid_geff <- function(newdata) {
        vars <- intersect(rf_vars_present_geff, names(newdata))
        if (!length(vars)) return(rep(0, nrow(newdata)))
        as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
      }
      
    } else {
      message("Decision: CO2_term passes regime monotonicity gate → keeping CO2_term.")
    }
  }
}

# --------------------------------------------------------------------------
# 5) Final predictions + diagnostics (messages preserved)
# --------------------------------------------------------------------------

pred_full <- predict_geff_ET_full(df)

output_df <- df |>
  mutate(
    # keep your components/diagnostics
    g_param    = g_param_full,
    ET_vpd_hat = ET_vpd_hat,
    
    # Jarvis-style harmonized names (downstream scripts expect these)
    g_eff_pred = pred_full$geff,
    ET_pred    = pred_full$ET,
    ET_resid   = ET - ET_pred
  )

RMSE_param_geff <- sqrt(mean((output_df$g_eff - output_df$g_param)^2, na.rm = TRUE))
R2_param_geff   <- suppressWarnings(cor(output_df$g_eff, output_df$g_param, use = "complete.obs")^2)

RMSE_pred_geff <- sqrt(mean((output_df$g_eff - output_df$g_eff_pred)^2, na.rm = TRUE))
R2_pred_geff   <- suppressWarnings(cor(output_df$g_eff, output_df$g_eff_pred, use = "complete.obs")^2)

RMSE_param_ET <- sqrt(mean((output_df$ET - output_df$ET_vpd_hat)^2, na.rm = TRUE))
R2_param_ET   <- suppressWarnings(cor(output_df$ET, output_df$ET_vpd_hat, use = "complete.obs")^2)

RMSE_pred_ET <- sqrt(mean((output_df$ET - output_df$ET_pred)^2, na.rm = TRUE))
R2_pred_ET   <- suppressWarnings(cor(output_df$ET, output_df$ET_pred, use = "complete.obs")^2)

message(sprintf(
  "g_eff — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Full: RMSE=%.3f, R^2=%.3f",
  RMSE_param_geff, R2_param_geff, RMSE_pred_geff, R2_pred_geff
))
message(sprintf(
  "ET    — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Full: RMSE=%.3f, R^2=%.3f",
  RMSE_param_ET, R2_param_ET, RMSE_pred_ET, R2_pred_ET
))

# --------------------------------------------------------------------------
# 6) UNIVERSAL bundle
# --------------------------------------------------------------------------

pred_fun_universal <- local({
  
  eps_local <- eps
  b0_local  <- b_hat[["b0_VPD"]]
  b1_local  <- b_hat[["b1_VPD"]]
  
  rf_geff_local <- rf_geff
  rf_vars_local <- rf_vars_present_geff
  recompute_K_ET_local <- recompute_K_ET
  
  .predict_rf_resid_geff_local <- function(newdata) {
    vars <- intersect(rf_vars_local, names(newdata))
    if (!length(vars)) return(rep(0, nrow(newdata)))
    as.numeric(predict(rf_geff_local, data = newdata[, vars, drop = FALSE])$predictions)
  }
  
  function(newdata) {
    vVPD <- pmax(newdata$VPD, eps_local)
    inv_sqrt_vpd <- 1 / sqrt(vVPD)
    
    g_param <- b0_local + b1_local * inv_sqrt_vpd
    K_ET    <- recompute_K_ET_local(newdata)
    g_resid <- .predict_rf_resid_geff_local(newdata)
    
    g_eff_hat <- g_param + g_resid
    ET_hat    <- K_ET * g_eff_hat
    
    list(ET = as.numeric(ET_hat), geff = as.numeric(g_eff_hat))
  }
})

rf_hybrid_bundle <- list(
  
  metadata = "Hybrid g_eff model: parametric 1/sqrt(VPD) + RF residuals; ET = K_ET × g_eff; CO2 regime monotonicity gate. FIT-ONLY universal bundle.",
  model_id = "rf_hybrid_vpd_constrained",
  
  # universal interface (for permutation/plots scripts)
  df_opt         = df_opt,
  df             = df,
  pred_fun       = pred_fun_universal,
  recompute_K_ET = recompute_K_ET,
  eps            = eps,
  
  # fitted VPD law
  b0_VPD = b_hat[["b0_VPD"]],
  b1_VPD = b_hat[["b1_VPD"]],
  
  # RF residual model
  rf_geff_fit    = rf_geff,
  rf_predictors  = rf_vars_present_geff,
  rf_params      = best_params_geff,
  rf_oob_rmse    = base_oob_rmse_geff,
  
  # CO2 gate
  co2_gate = co2_gate_info,
  
  # outputs + diagnostics
  output_df       = output_df,
  RMSE_param_geff = RMSE_param_geff,
  R2_param_geff   = R2_param_geff,
  RMSE_pred_geff = RMSE_pred_geff,
  R2_pred_geff   = R2_pred_geff,
  RMSE_param_ET   = RMSE_param_ET,
  R2_param_ET     = R2_param_ET,
  RMSE_pred_ET   = RMSE_pred_ET,
  R2_pred_ET     = R2_pred_ET,
  
  # reproducibility
  rng_seed     = rng_seed,
  rf_num_trees = rf_num_trees,
  rf_tune_mode = rf_tune_mode
)

save(rf_hybrid_bundle, file = out_file)

message("RF hybrid model fitted and universal bundle saved to: ", out_file)
message("RF predictors used: ", paste(rf_vars_present_geff, collapse = ", "))
