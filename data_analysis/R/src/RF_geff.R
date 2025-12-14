# ------------------------------------
# 1) Define predictors for g_eff RF
#    (same as ET RF)
# ------------------------------------
rf_vars_full_geff <- c("VPD", "Ta", "P", "Rg", "CO2_term")

df_rf_geff <- df_opt %>%
  select(g_eff, all_of(rf_vars_full_geff)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

target_col_geff <- "g_eff"

# ------------------------------------
# 2) Fit tuned RF on g_eff directly
#    (reuse your tuning helpers)
# ------------------------------------
set.seed(if (exists("rng_seed")) rng_seed else 123)

if (rf_tune_mode == "none") {
  rf_geff <- fit_rf_fixed(
    df_rf_geff,
    target_col = target_col_geff,
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
      target_col_geff,
      num.trees       = rf_num_trees,
      min.node.size   = c(3L,5L,10L),
      sample.fraction = c(0.8,1.0),
      splitrule       = c("variance","extratrees"),
      max.depth       = c(0L,12L)
    )
  } else {
    tuned_geff <- tune_rf_oob(
      df_rf_geff,
      target_col_geff,
      num.trees = rf_num_trees
    )
  }
  rf_geff          <- tuned_geff$best_fit
  best_params_geff <- tuned_geff$best_params
}

base_oob_rmse_geff <- sqrt(rf_geff$prediction.error)
rf_vars_present_geff <- setdiff(names(df_rf_geff), target_col_geff)

message(sprintf("RF_g_eff OOB RMSE = %.3f", base_oob_rmse_geff))
message(sprintf("RF_g_eff params: mtry=%s, min.node.size=%s, sample.fraction=%.2f, splitrule=%s, max.depth=%s",
                best_params_geff$mtry, best_params_geff$min.node.size,
                best_params_geff$sample.fraction,
                as.character(best_params_geff$splitrule), best_params_geff$max.depth))

# ------------------------------------
# 3) Variable importance for g_eff RF
# ------------------------------------
if (!is.null(rf_geff$variable.importance)) {
  vi_geff <- tibble(
    feature    = names(rf_geff$variable.importance),
    importance = as.numeric(rf_geff$variable.importance)
  ) %>% arrange(desc(importance))
  
  p_vip_geff <- ggplot(vi_geff,
                       aes(x = reorder(feature, importance), y = importance)) +
    geom_col() + coord_flip() + theme_bw() +
    labs(
      x = NULL,
      y = "RF importance (impurity)",
      title = expression("g"[eff]*" RF: Variable Importance (VPD, Ta, P, Rg, CO"[2]*")")
    )
  print(p_vip_geff)
}

# ------------------------------------
# 4) Safe prediction helper for g_eff
# ------------------------------------
.predict_rf_geff <- function(newdata) {
  vars <- intersect(rf_vars_present_geff, colnames(newdata))
  if (length(vars) == 0) return(rep(NA_real_, nrow(newdata)))
  as.numeric(predict(rf_geff, data = newdata[, vars, drop = FALSE])$predictions)
}

g_eff_rf_hat <- .predict_rf_geff(df)

df_obs_pred_geff <- df %>%
  mutate(g_eff_pred_RF = g_eff_rf_hat)

rmse_geff_RF <- sqrt(mean((df_obs_pred_geff$g_eff - df_obs_pred_geff$g_eff_pred_RF)^2, na.rm = TRUE))
r2_geff_RF   <- suppressWarnings(
  cor(df_obs_pred_geff$g_eff, df_obs_pred_geff$g_eff_pred_RF, use = "complete.obs")^2
)

p_obs_pred_geff_RF <- ggplot(df_obs_pred_geff,
                             aes(x = g_eff, y = g_eff_pred_RF)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff], " (RF-only)")),
    title = expression(paste("Observed vs Predicted ", g[eff], " (RF-only model)")),
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", rmse_geff_RF, r2_geff_RF)
  )

print(p_obs_pred_geff_RF)

library(pdp)
library(patchwork)

pdp_plots <- lapply(rf_vars_present_geff, function(v) {
  pd <- partial(
    object   = rf_geff,
    pred.var = v,
    train    = df_rf_geff,
    grid.resolution = 50,
    progress = "none"
  )
  
  ggplot(pd, aes_string(x = v, y = "yhat")) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(
      x = v,
      y = expression(hat(g)[eff]),
      title = paste("Partial dependence of", v)
    )
})

wrap_plots(pdp_plots, ncol = 2)

direction_summary <- lapply(rf_vars_present_geff, function(v) {
  pd <- partial(
    rf_geff,
    pred.var = v,
    train = df_rf_geff,
    grid.resolution = 50,
    progress = "none"
  )
  
  slope <- diff(pd$yhat) / diff(pd[[v]])
  
  tibble(
    predictor = v,
    mean_slope = mean(slope, na.rm = TRUE),
    direction  = case_when(
      mean_slope > 0  ~ "positive",
      mean_slope < 0  ~ "negative",
      TRUE            ~ "neutral"
    )
  )
})

bind_rows(direction_summary)


