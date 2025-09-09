# ---- Packages ----
library(DEoptim)
library(numDeriv)
library(tidyverse)

eps <- 1e-6  # guard for VPD

# ======================
# SWITCHBOARD (turn terms on/off)
# ======================
use_AI        <- TRUE
use_Rg        <- TRUE
use_intercept <- TRUE
use_VPD       <- TRUE
use_Ta        <- TRUE
use_CO2       <- TRUE   # you already had this (global CO2 toggle)

# ---- CO2 term prep (unchanged) ----
dln_CO2 <- log(CO2_2076_2100_RCP85 / CO2_1981_2005)  # ~1.003

Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    # Only CMIP*/CMIP6 AND future period get a CO2 term
    CO2_term = if_else(grepl("CMIP", model) & PERIOD == "2076_2100", dln_CO2, 0)
  )

# ---- Model (parameters = c(par_AI, par_Rg, par_intercept, par_VPD, par_Ta, par_CO2)) ----
model_ET_components <- function(par, data) {
  # Named parameters for clarity:
  par_AI        <- par[1]  # a: slope on AI_FAO56_alfalfa
  par_Rg        <- par[2]  # b: slope on Rg
  par_intercept <- par[3]  # c: intercept
  par_VPD       <- par[4]  # m: amplitude on 1/sqrt(VPD)
  par_Ta        <- par[5]  # d: slope on Ta
  par_CO2       <- par[6]  # k: CO2 elasticity (used via pmax)
  
  # Apply switches (turn param contributions off if switch is FALSE)
  par_AI_eff        <- if (use_AI)        par_AI        else 0
  par_Rg_eff        <- if (use_Rg)        par_Rg        else 0
  par_intercept_eff <- if (use_intercept) par_intercept else 0
  par_VPD_eff       <- if (use_VPD)       par_VPD       else 0
  par_Ta_eff        <- if (use_Ta)        par_Ta        else 0
  
  VPDg <- pmax(data$VPD, eps)
  
  # Base/reference conductance + VPD term
  g_eff_ref  <- par_AI_eff*data$AI_FAO56_alfalfa + par_Rg_eff*data$Rg + par_intercept_eff + par_Ta_eff*data$Ta
  g_eff_base <- g_eff_ref + par_VPD_eff / sqrt(VPDg)
  
  # CO2 decline only where flagged; allow "no-CO2 effect" via pmax(k,0) and switch
  k_eff <- if (isTRUE(use_CO2)) pmax(par_CO2, 0) else 0
  g_eff_pred <- g_eff_base * exp(-k_eff * data$CO2_term)  # stays positive
  
  # ET conversion (check your unit chain carefully)
  ET_pred <- (rhoAir * CpAir / gamma) * VPDg /
    (1 / (g_eff_pred / 1000)) *                    # g_eff in mmol m^-2 s^-1 -> mol m^-2 s^-1
    1 / ((2.501e6 - 2361 * data$Ta) / 1e6) *       # latent heat (J kg^-1) -> MJ kg^-1
    (3600 * 24) / 1e6 * 365.25                     # W m^-2 -> MJ m^-2 d^-1 -> per year
  
  list(ET_pred = ET_pred, g_eff_pred = g_eff_pred)
}

# ---- Objective: SSE in ET space, with penalties for nonphysical/NA ----
obj_fun_ET <- function(par, data) {
  comps <- model_ET_components(par, data)
  if (!all(is.finite(comps$ET_pred)) || any(!is.finite(comps$g_eff_pred))) return(1e12)
  if (any(comps$g_eff_pred <= 0)) return(1e12 + 1e6 * sum(comps$g_eff_pred <= 0))
  res <- data$ET - comps$ET_pred
  if (!all(is.finite(res))) return(1e12)
  sum(res^2)
}

# ---- Fit data ----
df_opt <- Data_to_plot_II |>
  filter(is.finite(ET),
         is.finite(AI_FAO56_alfalfa),
         is.finite(Rg),
         is.finite(VPD),
         is.finite(Ta))

# ---- DEoptim bounds (respect switches by fixing disabled params at 0) ----
# par = c(par_AI, par_Rg, par_intercept, par_VPD, par_Ta, par_CO2)
lower <- c(
  par_AI        = if (use_AI)        -10   else 0,
  par_Rg        = if (use_Rg)        -10   else 0,
  par_intercept = if (use_intercept) -1000 else 0,
  par_VPD       = if (use_VPD)       -100  else 0,
  par_Ta        = if (use_Ta)        -10   else 0,
  par_CO2       = if (use_CO2)       -5    else 0
)
upper <- c(
  par_AI        = if (use_AI)        10    else 0,
  par_Rg        = if (use_Rg)        10    else 0,
  par_intercept = if (use_intercept) 1000  else 0,
  par_VPD       = if (use_VPD)       100   else 0,
  par_Ta        = if (use_Ta)        10    else 0,
  par_CO2       = if (use_CO2)       5     else 0
)

set.seed(153)
de_fit <- DEoptim(
  fn     = obj_fun_ET,
  lower  = lower,
  upper  = upper,
  data   = df_opt,
  control = DEoptim.control(NP = 100, itermax = 600, reltol = 1e-8, trace = TRUE)
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- c("par_AI","par_Rg","par_intercept","par_VPD","par_Ta","par_CO2")
SSE_hat <- de_fit$optim$bestval

# ---- Optional local polish from DE solution ----
polished <- try(
  optim(par = par_hat, fn = obj_fun_ET, data = df_opt, method = "BFGS",
        control = list(reltol = 1e-12, maxit = 1000)),
  silent = TRUE
)
if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < SSE_hat) {
  par_hat <- polished$par
  names(par_hat) <- c("par_AI","par_Rg","par_intercept","par_VPD","par_Ta","par_CO2")
  SSE_hat <- polished$value
}

# ---- Prediction on the full data ----
comps_full <- model_ET_components(par_hat, Data_to_plot_II)
Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    g_eff_predicted = comps_full$g_eff_pred,
    ET_predicted    = comps_full$ET_pred,
    ET_resids       = ET - ET_predicted
  )

# ---- Quick fit metrics (optional) ----
pred_fit <- model_ET_components(par_hat, df_opt)$ET_pred
obs_fit  <- df_opt$ET
RMSE     <- sqrt(mean((obs_fit - pred_fit)^2, na.rm = TRUE))
R2       <- suppressWarnings(cor(obs_fit, pred_fit, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# ======================
#  Jacobian at the optimum (Gauss–Newton covariance)
#  Cov(theta) ≈ σ² * (J'J)^(-1) on ACTIVE parameters only
# ======================
p <- length(par_hat)
n <- nrow(df_opt)

residuals_fun <- function(th, data) {
  comps <- model_ET_components(th, data)
  r <- data$ET - comps$ET_pred
  r[!is.finite(r)] <- 0
  r
}

# Jacobian at θ̂ (all params)
J_all <- jacobian(function(th) residuals_fun(th, df_opt), par_hat)

# Decide active parameters for SEs:
k_used <- if (isTRUE(use_CO2)) max(unname(par_hat["par_CO2"]), 0) else 0
active_flags <- c(
  par_AI        = use_AI,
  par_Rg        = use_Rg,
  par_intercept = use_intercept,
  par_VPD       = use_VPD,
  par_Ta        = use_Ta,
  par_CO2       = (use_CO2 && k_used > 0)  # CO2 only if switch ON and effective
)
active_idx <- which(active_flags)

# Residual variance on ET scale
SSE_hat <- sum(residuals_fun(par_hat, df_opt)^2)
sigma2  <- SSE_hat / (n - length(active_idx))

# Invert only the active block (more stable), fill SEs for actives, NA for inactives
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
    parameter == "par_CO2" & k_used == 0 ~ "inactive (no CO2 effect at optimum)",
    TRUE ~ NA_character_
  )
)
print(param_table)

# ---- Quick model diagnostics plot ----
with(Data_to_plot_II, {
  plot(ET, ET_predicted, xlab = "Observed ET", ylab = "Predicted ET")
  abline(0, 1, col = "red")
})

# ---- CO2 effect summary (on flagged rows) ----
k_raw  <- unname(par_hat["par_CO2"])
k_eff  <- if (isTRUE(use_CO2)) max(k_raw, 0) else 0
decline_factor <- exp(-k_eff * dln_CO2)
pct_change <- (decline_factor - 1) * 100
message(sprintf("CO2: k_raw = %.4f, k_used = %.4f, effect on flagged rows = %.2f%% in g_eff",
                k_raw, k_eff, pct_change))
