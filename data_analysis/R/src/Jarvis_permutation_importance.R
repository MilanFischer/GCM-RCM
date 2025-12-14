# ================================
# Permutation importance — Jarvis
# (run AFTER fitting; loads bundle)
# ================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- user settings ----
out_dir <- "../outputs"
plot_dir <- "../Plots/ggplot2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

perm_B_jarvis <- 1000L
seed_perm     <- 1990  # or jarvis_bundle$rng_seed

# ---- load model bundle ----
stopifnot(exists("jarvis_bundle"))

# pull things we need (avoid relying on global env)
par_hat   <- jarvis_bundle$jarvis_pars
fVPD_mode <- jarvis_bundle$jarvis_fVPD_mode
df_opt    <- jarvis_bundle$jarvis_df_opt

rhoAir <- jarvis_bundle$rhoAir
CpAir  <- jarvis_bundle$CpAir
gamma  <- jarvis_bundle$gamma
eps    <- jarvis_bundle$eps %||% 1e-6

use_P   <- jarvis_bundle$use_P
use_Rg  <- jarvis_bundle$use_Rg
use_Ta  <- jarvis_bundle$use_Ta
use_CO2 <- jarvis_bundle$use_CO2
use_VPD <- jarvis_bundle$use_VPD

# function handles
model_ET_components <- jarvis_bundle$jarvis_model_ET_components

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------
# Unified predictor (critical!)
# -----------------------
predict_jarvis <- function(newdata) {
  nd <- newdata %>%
    mutate(
      VPD = pmax(VPD, eps),
      inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
      K_ET = (rhoAir * CpAir / gamma) * VPD * (1/1000) *
        1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
    )
  
  comps <- model_ET_components(par_hat, nd, fVPD_mode = fVPD_mode)
  
  list(
    ET   = as.numeric(comps$ET_pred),
    geff = as.numeric(comps$g_eff_pred)
  )
}

# -----------------------
# Vars to permute (respect switches)
# -----------------------
vars_perm <- c()
if (isTRUE(use_VPD)) vars_perm <- c(vars_perm, "VPD")
if (isTRUE(use_Ta))  vars_perm <- c(vars_perm, "Ta")
if (isTRUE(use_P))   vars_perm <- c(vars_perm, "P")
if (isTRUE(use_Rg))  vars_perm <- c(vars_perm, "Rg")
if (isTRUE(use_CO2)) vars_perm <- c(vars_perm, "CO2_term")
vars_perm <- intersect(vars_perm, names(df_opt))
stopifnot(length(vars_perm) >= 1)

# -----------------------
# Permutation importance (model-agnostic)
# -----------------------
perm_importance_model <- function(data, vars, pred_fun, B = 200L, seed = 2024) {
  set.seed(seed)
  
  pred0 <- pred_fun(data)
  base_rmse_ET   <- sqrt(mean((data$ET    - pred0$ET  )^2, na.rm = TRUE))
  base_rmse_geff <- sqrt(mean((data$g_eff - pred0$geff)^2, na.rm = TRUE))
  
  bind_rows(lapply(vars, function(v) {
    dET   <- numeric(B)
    dGeff <- numeric(B)
    
    for (b in seq_len(B)) {
      dfp <- data
      dfp[[v]] <- sample(dfp[[v]])
      
      # prevent stale derived cols
      dfp$inv_sqrt_VPD <- NULL
      dfp$K_ET <- NULL
      
      predp <- pred_fun(dfp)
      
      rmse_ET_p   <- sqrt(mean((data$ET    - predp$ET  )^2, na.rm = TRUE))
      rmse_geff_p <- sqrt(mean((data$g_eff - predp$geff)^2, na.rm = TRUE))
      
      dET[b]   <- rmse_ET_p   - base_rmse_ET
      dGeff[b] <- rmse_geff_p - base_rmse_geff
    }
    
    tibble(
      feature = v,
      mean_delta_ET   = mean(dET, na.rm = TRUE),
      lo_ET           = quantile(dET, 0.05, na.rm = TRUE),
      hi_ET           = quantile(dET, 0.95, na.rm = TRUE),
      mean_delta_geff = mean(dGeff, na.rm = TRUE),
      lo_geff         = quantile(dGeff, 0.05, na.rm = TRUE),
      hi_geff         = quantile(dGeff, 0.95, na.rm = TRUE),
      B = B
    )
  }))
}

perm_jarvis <- perm_importance_model(
  data     = df_opt,
  vars     = vars_perm,
  pred_fun = predict_jarvis,
  B        = perm_B_jarvis,
  seed     = seed_perm
)

print(perm_jarvis)
write_csv(perm_jarvis, file.path(out_dir, "perm_Jarvis_fit.csv"))

# -----------------------
# VPD pathway decomposition (ET)
# -----------------------
recompute_K_ET_local <- function(d) {
  vVPD <- pmax(d$VPD, eps)
  (rhoAir * CpAir / gamma) * vVPD * (1/1000) *
    1 / ((2.501e6 - 2361 * d$Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

perm_importance_VPD_channels_ET_jarvis <- function(data, B = 200L, seed = 2024) {
  set.seed(seed)
  
  pred0 <- predict_jarvis(data)
  base_rmse_ET <- sqrt(mean((data$ET - pred0$ET)^2, na.rm = TRUE))
  
  K0    <- recompute_K_ET_local(data)
  geff0 <- pred0$geff
  
  d_full   <- numeric(B)
  d_stom   <- numeric(B)
  d_demand <- numeric(B)
  
  for (b in seq_len(B)) {
    dfp <- data
    dfp$VPD <- sample(dfp$VPD)
    
    dfp$inv_sqrt_VPD <- NULL
    dfp$K_ET <- NULL
    
    predp <- predict_jarvis(dfp)
    rmse_full <- sqrt(mean((data$ET - predp$ET)^2, na.rm = TRUE))
    d_full[b] <- rmse_full - base_rmse_ET
    
    geff_p <- predict_jarvis(dfp)$geff
    ET_stom <- K0 * geff_p
    rmse_stom <- sqrt(mean((data$ET - ET_stom)^2, na.rm = TRUE))
    d_stom[b] <- rmse_stom - base_rmse_ET
    
    K_p <- recompute_K_ET_local(dfp)
    ET_dem <- K_p * geff0
    rmse_dem <- sqrt(mean((data$ET - ET_dem)^2, na.rm = TRUE))
    d_demand[b] <- rmse_dem - base_rmse_ET
  }
  
  tibble(
    feature = "VPD",
    channel = c("Full (g + K_ET)", "Via g_eff (stomatal; K fixed)", "Via K_ET (demand; g fixed)"),
    mean_delta_ET = c(mean(d_full, na.rm = TRUE),
                      mean(d_stom, na.rm = TRUE),
                      mean(d_demand, na.rm = TRUE)),
    lo_ET = c(quantile(d_full, 0.05, na.rm = TRUE),
              quantile(d_stom, 0.05, na.rm = TRUE),
              quantile(d_demand, 0.05, na.rm = TRUE)),
    hi_ET = c(quantile(d_full, 0.95, na.rm = TRUE),
              quantile(d_stom, 0.95, na.rm = TRUE),
              quantile(d_demand, 0.95, na.rm = TRUE)),
    B = B
  )
}

perm_vpd_channels <- perm_importance_VPD_channels_ET_jarvis(
  df_opt, B = perm_B_jarvis, seed = seed_perm
)

print(perm_vpd_channels)
write_csv(perm_vpd_channels, file.path(out_dir, "perm_Jarvis_VPD_channels_ET.csv"))

# -----------------------
# Plots
# -----------------------
p_perm_ET <- ggplot(perm_jarvis, aes(x = reorder(feature, mean_delta_ET), y = mean_delta_ET)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = expression(Delta * " RMSE_ET"),
       title = "Permutation importance — Jarvis ET")

p_perm_geff <- ggplot(perm_jarvis, aes(x = reorder(feature, mean_delta_geff), y = mean_delta_geff)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_geff, ymax = hi_geff), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = expression(Delta * " RMSE_" * g[eff]),
       title = expression("Permutation importance — Jarvis " * g[eff]))

# VPD channels plot
p_perm_vpd_channels_ET <- ggplot(perm_vpd_channels,
                                 aes(x = reorder(channel, mean_delta_ET), y = mean_delta_ET)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo_ET, ymax = hi_ET), width = 0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = expression(Delta * " RMSE_ET (permute VPD)"),
       title = "VPD permutation importance decomposed by pathway (ET) — Jarvis")

print(p_perm_ET)
print(p_perm_geff)
print(p_perm_vpd_channels_ET)

ggsave(file.path(plot_dir, "perm_Jarvis_ET.png"), p_perm_ET, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(plot_dir, "perm_Jarvis_geff.png"), p_perm_geff, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(plot_dir, "perm_Jarvis_VPD_channels_ET.png"), p_perm_vpd_channels_ET, width = 7, height = 4.5, dpi = 300)
