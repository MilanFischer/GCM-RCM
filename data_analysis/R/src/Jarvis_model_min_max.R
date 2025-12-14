# ================
# Jarvis ET model
# ===============

suppressPackageStartupMessages({
  library(tidyverse)
  library(DEoptim)
  library(numDeriv)
  library(MASS)
})

source("./src/Jarvis_preprocessing.R")

# -----------------------
# Switchboard
# -----------------------
use_P   <- TRUE
use_Rg  <- TRUE
use_Ta  <- TRUE
use_CO2 <- TRUE
use_VPD <- TRUE

# fVPD handling
fVPD_mode  <- "optimize"   # "given" or "optimize"
b0_VPD_raw <- -12.613      # used if fVPD_mode == "given"
b1_VPD_raw <-  11.093

# Optimizer settings
trace_flag <- TRUE
rng_seed   <- 1990
polish     <- TRUE

clamp_Rg   <- TRUE
clamp_Ta   <- TRUE
clamp_P    <- TRUE

n_max_iter <- 25000

# Metadata / output
metadata <- "Jarvis parameter set and output tibble. Created 2025-12-14 by Milan Fischer."
out_file <- "./RData/20251214_jarvis_objects.RData"

# numeric guard
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# -----------------------
# Parameter helpers
# -----------------------
get_par_names <- function(fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
  if (fVPD_mode == "optimize") {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2","b0_VPD","b1_VPD")
  } else {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2")
  }
}

.require_pars <- function(par, needed) {
  miss <- setdiff(needed, names(par))
  if (length(miss)) stop("Missing parameters in 'par': ", paste(miss, collapse = ", "))
}

# -----------------------
# Jarvis: g_eff (mm/s)
# -----------------------
jarvis_g_eff <- function(par, input_df, fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
  
  base_needed <- if (fVPD_mode == "optimize") {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2","b0_VPD","b1_VPD")
  } else {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2")
  }
  .require_pars(par, base_needed)
  
  Rg_min <- par[["Rg_min"]]; Rg_max <- par[["Rg_max"]]
  Ta_min <- par[["Ta_min"]]; Ta_max <- par[["Ta_max"]]
  P_min  <- par[["P_min"]] ; P_max  <- par[["P_max"]]
  k_CO2  <- par[["k_CO2"]]
  
  # keep bounds ordered
  if (Rg_min > Rg_max) { tmp <- Rg_min; Rg_min <- Rg_max; Rg_max <- tmp }
  if (Ta_min > Ta_max) { tmp <- Ta_min; Ta_min <- Ta_max; Ta_max <- tmp }
  if (P_min  > P_max ) { tmp <- P_min ; P_min  <- P_max ; P_max  <- tmp }
  
  # VPD coefs
  if (fVPD_mode == "optimize") {
    b0_VPD <- par[["b0_VPD"]]; b1_VPD <- par[["b1_VPD"]]
  } else {
    b0_VPD <- b0_VPD_raw;      b1_VPD <- b1_VPD_raw
  }
  
  # pull
  Rg <- input_df$Rg
  Ta <- input_df$Ta
  P  <- input_df$P
  CO2_term     <- input_df$CO2_term
  inv_sqrt_VPD <- input_df$inv_sqrt_VPD
  
  # clamped affine responses
  spanRg <- pmax(Rg_max - Rg_min, eps)
  spanTa <- pmax(Ta_max - Ta_min, eps)
  spanP  <- pmax(P_max  - P_min , eps)
  
  # use the clamped versions to keep scale in [0,1]
  fRg <- if (use_Rg && clamp_Rg) pmin(1, pmax(0, (Rg - Rg_min) / spanRg)) else if (use_Rg) (Rg - Rg_min) / spanRg else 1
  fTa <- if (use_Ta && clamp_Ta) pmin(1, pmax(0, (Ta - Ta_min) / spanTa)) else if (use_Ta) (Ta - Ta_min) / spanTa else 1
  fP  <- if (use_P  && clamp_P)  pmin(1, pmax(0, (P  - P_min ) / spanP )) else if (use_P)  (P  - P_min ) / spanP  else 1
  
  fCO2 <- if (use_CO2) exp(-pmax(k_CO2, 0) * CO2_term) else 1
  fVPD <- if (use_VPD) (b0_VPD + b1_VPD * inv_sqrt_VPD) else 1
  
  fP * fRg * fTa * fVPD * fCO2
}

# -----------------------
# Model wrapper -> ET
# -----------------------
model_ET_components <- function(par, input_df, fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
  g_eff_pred <- jarvis_g_eff(par, input_df, fVPD_mode)
  ET_pred <- input_df$K_ET * g_eff_pred
  list(ET_pred = ET_pred, g_eff_pred = g_eff_pred)
}

# -----------------------
# Objective function
# -----------------------
obj_fun_ET <- function(par, input_df, fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
  names(par) <- get_par_names(fVPD_mode)
  ET <- input_df$ET
  
  g_eff_pred <- jarvis_g_eff(par, input_df, fVPD_mode)
  ET_pred  <- input_df$K_ET * g_eff_pred
  
  if (!all(is.finite(ET_pred))) return(1e12)
  sqrt(mean((ET - ET_pred)^2, na.rm = TRUE))
}

# -----------------------
# Bounds
# -----------------------
par_names <- get_par_names(fVPD_mode)
if (fVPD_mode == "optimize") {
  lower <- c(Rg_min=-500, Rg_max=50, Ta_min=-1000, Ta_max=5, P_min=-1000, P_max=400, k_CO2=0, b0_VPD=-50000, b1_VPD=-50000)
  upper <- c(Rg_min=190,  Rg_max=5000, Ta_min=20,  Ta_max=100, P_min=1000,  P_max=50000, k_CO2=5, b0_VPD=50000,  b1_VPD=50000)
} else {
  lower <- c(Rg_min=-500, Rg_max=50, Ta_min=-1000, Ta_max=5, P_min=-1000, P_max=400, k_CO2=0)
  upper <- c(Rg_min=190,  Rg_max=5000, Ta_min=20,  Ta_max=100, P_min=1000,  P_max=50000, k_CO2=5)
}
lower <- lower[par_names]; upper <- upper[par_names]

D  <- length(par_names)
NP <- max(10L * D, 60L)
set.seed(rng_seed)

# speed
compiler::enableJIT(3)
obj_fun_ET          <- compiler::cmpfun(obj_fun_ET)
model_ET_components <- compiler::cmpfun(model_ET_components)
jarvis_g_eff        <- compiler::cmpfun(jarvis_g_eff)

# -----------------------
# Global → local fit
# -----------------------
de_fit <- DEoptim::DEoptim(
  fn        = obj_fun_ET,
  lower     = lower,
  upper     = upper,
  input_df  = df_opt,
  fVPD_mode = fVPD_mode,
  control = DEoptim::DEoptim.control(
    NP = NP, itermax = n_max_iter, reltol = 1e-8, trace = trace_flag,
    parallelType = 0
  )
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- par_names
SSE_hat <- de_fit$optim$bestval

if (isTRUE(polish)) {
  polished <- try(
    optim(
      par       = par_hat,
      fn        = obj_fun_ET,
      input_df  = df_opt,
      fVPD_mode = fVPD_mode,
      method    = "L-BFGS-B",
      lower     = lower,
      upper     = upper,
      control   = list(factr = 1e4, pgtol = 0, maxit = 10000)
    ),
    silent = TRUE
  )
  if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < SSE_hat) {
    par_hat <- polished$par
    names(par_hat) <- par_names
    SSE_hat <- polished$value
  }
}

# ============================================================
# Unified prediction wrapper
#  - Always recomputes inv_sqrt_VPD and K_ET from the *current* VPD/Ta
# ============================================================
predict_jarvis <- function(newdata) {
  nd <- newdata |> 
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
# Predictions
# -----------------------
preds <- predict_jarvis(df)

output_df <- df |> 
  mutate(
    g_eff_predicted = preds$geff,
    ET_predicted    = preds$ET,
    ET_resids       = ET - ET_predicted
  )

RMSE <- sqrt(mean((output_df$ET - output_df$ET_predicted)^2, na.rm = TRUE))
R2   <- suppressWarnings(cor(output_df$ET, output_df$ET_predicted, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# -----------------------
# Parameter table (SEs via Jacobian)
# -----------------------
residuals_fun <- function(th, data) {
  comps <- model_ET_components(th, data, fVPD_mode = fVPD_mode)
  r <- data$ET - comps$ET_pred
  r[!is.finite(r)] <- 0
  r
}

J_all <- numDeriv::jacobian(function(th) residuals_fun(th, df_opt), par_hat)
SSE_hat <- sum(residuals_fun(par_hat, df_opt)^2)
n <- nrow(df_opt); p <- length(par_hat)
sigma2 <- SSE_hat / max(1, n - p)

XtX <- crossprod(J_all)
invXtX <- try(solve(XtX), silent = TRUE)
if (inherits(invXtX, "try-error")) invXtX <- MASS::ginv(XtX)
cov_theta <- sigma2 * invXtX

se <- sqrt(pmax(diag(cov_theta), 0))
se <- setNames(se, names(par_hat))

param_table <- tibble(
  parameter = names(par_hat),
  estimate  = as.numeric(par_hat),
  std_error = as.numeric(se[names(par_hat)]),
  z_value   = estimate / std_error,
  p_value   = 2 * pnorm(abs(z_value), lower.tail = FALSE)
)
print(param_table)

# -----------------------
# Save bundle
# -----------------------
jarvis_bundle <- list(
  metadata         = metadata,
  jarvis_fVPD_mode = fVPD_mode,
  jarvis_pars      = par_hat,
  jarvis_df_opt    = df_opt,
  jarvis_out       = output_df,
  
  jarvis_g_eff               = jarvis_g_eff,
  jarvis_model_ET_components = model_ET_components,
  .require_pars              = .require_pars,
  
  rhoAir = rhoAir, CpAir = CpAir, gamma = gamma,
  eps   = eps,
  
  use_P   = use_P,   use_Rg  = use_Rg,  use_Ta = use_Ta,
  use_CO2 = use_CO2, use_VPD = use_VPD,
  clamp_Rg = clamp_Rg, clamp_Ta = clamp_Ta, clamp_P = clamp_P,
  b0_VPD_raw = b0_VPD_raw, b1_VPD_raw = b1_VPD_raw
  
)

save(jarvis_bundle, file = out_file)
