perm_importance_full_model <- function(
    data,
    vars = c("VPD", "Ta", "P", "Rg", "CO2_term"),
    B = 100,
    seed = 2024
) {
  # Only use vars that actually exist
  vars <- intersect(vars, names(data))
  if (!length(vars)) {
    stop("None of the requested vars found in data.")
  }
  
  # Baseline RMSE for ET and g_eff using full model
  pred0 <- predict_geff_ET(data)
  base_rmse_ET   <- sqrt(mean((data$ET   - pred0$ET  )^2, na.rm = TRUE))
  base_rmse_geff <- sqrt(mean((data$g_eff - pred0$geff)^2, na.rm = TRUE))
  
  purrr::map_dfr(vars, function(v) {
    set.seed(seed)
    dET   <- numeric(B)
    dGeff <- numeric(B)
    
    for (b in seq_len(B)) {
      dfp <- data
      dfp[[v]] <- sample(dfp[[v]])  # permute v
      
      pred_p <- predict_geff_ET(dfp)
      
      rmse_ET_p   <- sqrt(mean((dfp$ET   - pred_p$ET  )^2, na.rm = TRUE))
      rmse_geff_p <- sqrt(mean((dfp$g_eff - pred_p$geff)^2, na.rm = TRUE))
      
      dET[b]   <- rmse_ET_p   - base_rmse_ET
      dGeff[b] <- rmse_geff_p - base_rmse_geff
    }
    
    tibble::tibble(
      feature        = v,
      mean_delta_ET   = mean(dET,   na.rm = TRUE),
      lo_ET           = stats::quantile(dET,   0.05, na.rm = TRUE),
      hi_ET           = stats::quantile(dET,   0.95, na.rm = TRUE),
      rel_ET_pct      = 100 * mean(dET,   na.rm = TRUE) / base_rmse_ET,
      mean_delta_geff = mean(dGeff, na.rm = TRUE),
      lo_geff         = stats::quantile(dGeff, 0.05, na.rm = TRUE),
      hi_geff         = stats::quantile(dGeff, 0.95, na.rm = TRUE),
      rel_geff_pct    = 100 * mean(dGeff, na.rm = TRUE) / base_rmse_geff,
      B               = B
    )
  }) %>%
    arrange(desc(mean_delta_ET))  # or desc(mean_delta_geff) if you prefer
}








perm_full <- perm_importance_full_model(
  df_opt,
  vars = c("VPD", "Ta", "P", "Rg", "CO2_term"),
  B = perm_B  # you already defined perm_B
)

print(perm_full %>% mutate(across(where(is.numeric), ~ round(.x, 3))))







# ET importance
p_imp_ET <- ggplot(perm_full,
                   aes(x = reorder(feature, mean_delta_ET), y = mean_delta_ET)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
  coord_flip() + theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta*" RMSE_ET (permute)"),
    title = "Full-model variable importance for ET",
    subtitle = "Permutation importance on two-stage ET model (incl. VPD)"
  )

# g_eff importance
p_imp_geff <- ggplot(perm_full,
                     aes(x = reorder(feature, mean_delta_geff), y = mean_delta_geff)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_geff, ymax = hi_geff), width = 0.2) +
  coord_flip() + theme_bw() +
  labs(
    x = NULL,
    y = expression(Delta*" RMSE_"*g[eff]*" (permute)"),
    title = expression("Full-model variable importance for " * g[eff]),
    subtitle = "Permutation importance on two-stage g_eff model (incl. VPD)"
  )

print(p_imp_ET + p_imp_geff)


# ================================
# Observed vs Predicted (full model)
# ================================
pred_full <- predict_geff_ET(df)

df_obs_pred <- df %>%
  mutate(
    ET_pred_full   = pred_full$ET,
    geff_pred_full = pred_full$geff
  )

rmse_ET_full   <- sqrt(mean((df_obs_pred$ET   - df_obs_pred$ET_pred_full  )^2, na.rm = TRUE))
r2_ET_full     <- suppressWarnings(cor(df_obs_pred$ET,   df_obs_pred$ET_pred_full,   use = "complete.obs")^2)
rmse_geff_full <- sqrt(mean((df_obs_pred$g_eff - df_obs_pred$geff_pred_full)^2, na.rm = TRUE))
r2_geff_full   <- suppressWarnings(cor(df_obs_pred$g_eff, df_obs_pred$geff_pred_full, use = "complete.obs")^2)

p_obs_pred_geff_full <- ggplot(
  df_obs_pred,
  aes(x = g_eff, y = geff_pred_full)
) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff], " (full model)")),
    title = expression(paste("Observed vs Predicted ", g[eff], " (full model)")),
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", rmse_geff_full, r2_geff_full)
  )

p_obs_pred_ET_full <- ggplot(
  df_obs_pred,
  aes(x = ET, y = ET_pred_full)
) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Observed ET",
    y = "Predicted ET (full model)",
    title = "Observed vs Predicted ET (full model)",
    subtitle = sprintf("RMSE = %.3f, R^2 = %.3f", rmse_ET_full, r2_ET_full)
  )

# If you want them side by side (patchwork loaded already):
print(p_obs_pred_geff_full + p_obs_pred_ET_full)

