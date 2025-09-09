# ---- Packages ----
library(DEoptim)
library(numDeriv)
library(tidyverse)

eps <- 1e-6  # guard

# ======================
# SWITCHBOARD (set FALSE to drop the effect)
# ======================
use_P   <- TRUE
use_Rg  <- TRUE
use_Ta  <- TRUE
use_VPD <- TRUE
use_CO2 <- TRUE

# ---- CO2 term ----
dln_CO2 <- log(CO2_2076_2100_RCP85 / CO2_1981_2005)
Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    CO2_term = if_else(grepl("CMIP", model) & PERIOD == "2076_2100", dln_CO2, 0)
  )

# ======================
# Build drivers ONCE
# ======================
build_drivers <- function(df) {

  # 1/sqrt(VPD)
  invVPD <- 1 / sqrt(pmax(df$VPD, eps))
  
  tibble(
    P      = df$P,
    Rg     = df$Rg,
    Ta     = df$Ta,
    invVPD = invVPD
  )
}

# ======================
# Jarvis g_eff (mm/s) with regression-style modifiers
# f(x) = a + b*x
# ======================
jarvis_g_eff <- function(drv, CO2_term,
                         a_P, b_P,
                         a_Rg, b_Rg,
                         a_Ta, b_Ta,
                         a_VPD, b_VPD,
                         k_CO2 = 0) {
  
  fP   <- if (use_P)   (a_P   + b_P   * drv$P)   else 1
  fRg  <- if (use_Rg)  (a_Rg  + b_Rg  * drv$Rg)     else 1
  fTa  <- if (use_Ta)  (a_Ta  + b_Ta  * drv$Ta)     else 1
  fVPD <- if (use_VPD) (a_VPD + b_VPD * drv$invVPD) else 1
  
  fCO2 <- if (use_CO2) exp(-pmax(k_CO2, 0) * CO2_term) else 1
  
  g_eff <- fP * fRg * fTa * fVPD * fCO2
  g_eff
}

# ======================
# Model wrapper -> ET
# par order: c(a_P,b_P, a_Rg,b_Rg, a_Ta,b_Ta, a_VPD,b_VPD, k_CO2)
# ======================
model_ET_components <- function(par, data, drv) {
  a_P   <- par[1];  b_P   <- par[2]
  a_Rg  <- par[3];  b_Rg  <- par[4]
  a_Ta  <- par[5];  b_Ta  <- par[6]
  a_VPD <- par[7];  b_VPD <- par[8]
  k_CO2 <- par[9]
  
  g_eff_pred <- jarvis_g_eff(
    drv      = drv,
    CO2_term = data$CO2_term,
    a_P=a_P, b_P=b_P,
    a_Rg=a_Rg, b_Rg=b_Rg,
    a_Ta=a_Ta, b_Ta=b_Ta,
    a_VPD=a_VPD, b_VPD=b_VPD,
    k_CO2=k_CO2
  )
  
  VPDg <- pmax(data$VPD, eps)
  
  # --- ET conversion ---
  # g_eff_pred is mm/s; convert to m/s inside resistance
  ET_pred <- (rhoAir * CpAir / gamma) * VPDg /
    (1 / (g_eff_pred / 1000)) *                    # s/m
    1 / ((2.501e6 - 2361 * data$Ta) / 1e6) *       # LE conversion (W m^-2)
    (3600 * 24) / 1e6 * 365.25                     # mm/yr
  
  list(ET_pred = ET_pred, g_eff_pred = g_eff_pred)
}

# ---- Objective ----
obj_fun_ET <- function(par, data, drv) {
  comps <- model_ET_components(par, data, drv)
  if (!all(is.finite(comps$ET_pred)) || any(!is.finite(comps$g_eff_pred))) return(1e12)
  if (any(comps$g_eff_pred <= 0)) return(1e12 + 1e6 * sum(comps$g_eff_pred <= 0))
  res <- data$ET - comps$ET_pred
  # res <- data$g_eff - comps$g_eff_pred
  if (!all(is.finite(res))) return(1e12)
  sum(res^2)
}

# ---- Fit subset ----
df_opt <- Data_to_plot_II |>
  filter(is.finite(ET),
         is.finite(P),
         is.finite(Rg),
         is.finite(VPD),
         is.finite(Ta))

drv_opt  <- build_drivers(df_opt)
drv_full <- build_drivers(Data_to_plot_II)

# ---- Bounds: c(a_P,b_P, a_Rg,b_Rg, a_Ta,b_Ta, a_VPD,b_VPD, k_CO2) ----
lower <- c(
  a_P=-2,   b_P=-2,
  a_Rg=-2,  b_Rg=-2,
  a_Ta=-2,  b_Ta=-2,
  a_VPD=-2, b_VPD=-2,
  k_CO2=0
)
upper <- c(
  a_P=  2,  b_P= 2,
  a_Rg= 2,  b_Rg= 2,
  a_Ta= 2,  b_Ta= 2,
  a_VPD= 2, b_VPD= 2,
  k_CO2= 5
)

set.seed(153)
de_fit <- DEoptim(
  fn      = obj_fun_ET,
  lower   = lower,
  upper   = upper,
  data    = df_opt,
  drv     = drv_opt,
  control = DEoptim.control(NP = 140, itermax = 1000, reltol = 1e-8, trace = TRUE)
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- c("a_P","b_P","a_Rg","b_Rg","a_Ta","b_Ta","a_VPD","b_VPD","k_CO2")
SSE_hat <- de_fit$optim$bestval

# ---- Optional polish ----
polished <- try(
  optim(par = par_hat, fn = obj_fun_ET, data = df_opt, drv = drv_opt,
        method = "BFGS", control = list(reltol = 1e-12, maxit = 1000)),
  silent = TRUE
)
if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < SSE_hat) {
  par_hat <- polished$par
  names(par_hat) <- c("a_P","b_P","a_Rg","b_Rg","a_Ta","b_Ta","a_VPD","b_VPD","k_CO2")
  SSE_hat <- polished$value
}

# ---- Predict on full data ----
comps_full <- model_ET_components(par_hat, Data_to_plot_II, drv_full)
Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    g_eff_predicted = comps_full$g_eff_pred,  # mm/s
    ET_predicted    = comps_full$ET_pred,
    ET_resids       = ET - ET_predicted
  )

# ---- Fit metrics ----
pred_fit <- model_ET_components(par_hat, df_opt, drv_opt)$ET_pred
obs_fit  <- df_opt$ET
RMSE     <- sqrt(mean((obs_fit - pred_fit)^2, na.rm = TRUE))
R2       <- suppressWarnings(cor(obs_fit, pred_fit, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# ---- Standard errors ----
residuals_fun <- function(th, data, drv) {
  comps <- model_ET_components(th, data, drv)
  r <- data$ET - comps$ET_pred
  r[!is.finite(r)] <- 0
  r
}
J_all <- jacobian(function(th) residuals_fun(th, df_opt, drv_opt), par_hat)

k_used <- if (isTRUE(use_CO2)) max(unname(par_hat["k_CO2"]), 0) else 0
active_flags <- c(
  a_P=use_P, b_P=use_P,
  a_Rg=use_Rg, b_Rg=use_Rg,
  a_Ta=use_Ta, b_Ta=use_Ta,
  a_VPD=use_VPD, b_VPD=use_VPD,
  k_CO2=(use_CO2 && k_used > 0)
)
active_idx <- which(active_flags)

SSE_hat <- sum(residuals_fun(par_hat, df_opt, drv_opt)^2)
n <- nrow(df_opt); p <- length(par_hat)
sigma2 <- SSE_hat / (n - length(active_idx))

se <- rep(NA_real_, p); names(se) <- names(par_hat)
if (length(active_idx) > 0) {
  J <- J_all[, active_idx, drop = FALSE]
  XtX <- crossprod(J)
  invXtX <- try(solve(XtX), silent = TRUE)
  if (inherits(invXtX, "try-error")) invXtX <- MASS::ginv(XtX)
  cov_theta_act <- sigma2 * invXtX
  se_act <- sqrt(pmax(diag(cov_theta_act), 0))
  se[active_idx] <- se_act
}

z  <- par_hat / se
pv <- 2 * pnorm(abs(z), lower.tail = FALSE)
param_table <- tibble(
  parameter = names(par_hat),
  estimate  = as.numeric(par_hat),
  std_error = as.numeric(se),
  z_value   = as.numeric(z),
  p_value   = as.numeric(pv),
  note      = case_when(
    !active_flags ~ "disabled by switch",
    parameter == "k_CO2" & k_used == 0 ~ "inactive (no CO2 effect at optimum)",
    TRUE ~ NA_character_
  )
)
print(param_table)

# ---- Diagnostic ----
with(Data_to_plot_II, {
  plot(ET, ET_predicted, xlab = "Observed ET", ylab = "Predicted ET")
  abline(0, 1, col = "red")
})
