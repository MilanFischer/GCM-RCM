# ================================
# Jarvis ET model (universal bundle)
# ================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DEoptim)
  library(numDeriv)
  library(MASS)
})

# ------------------------------------------------------------------------------
# REQUIRED INPUTS (must already exist in your session)
# ------------------------------------------------------------------------------
if (!exists("Data_to_plot", inherits = TRUE) || is.null(Data_to_plot[["abs"]])) {
  stop("Missing Data_to_plot[['abs']] in your environment.", call. = FALSE)
}
if (!all(c("rhoAir","CpAir","gamma") %in% ls(envir = .GlobalEnv))) {
  stop("Missing physical constants: 'rhoAir', 'CpAir', 'gamma'.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Load shared preprocessing (creates df, df_opt, eps, recompute_K_ET, CO2_term, etc.)
# ------------------------------------------------------------------------------
source("./src/Jarvis&RF_preprocessing.R")

stopifnot(exists("df", inherits = TRUE), exists("df_opt", inherits = TRUE))
stopifnot(exists("recompute_K_ET", inherits = TRUE))
stopifnot(all(c("ET","g_eff","Rg","Ta","P","VPD","CO2_term","inv_sqrt_VPD","K_ET") %in% names(df_opt)))

# ------------------------------------------------------------------------------
# Switchboard
# ------------------------------------------------------------------------------
use_P   <- TRUE
use_Rg  <- TRUE
use_Ta  <- TRUE
use_CO2 <- TRUE
use_VPD <- TRUE

# VPD response handling
fVPD_mode  <- "optimize"   # "given" or "optimize"
b0_VPD_raw <- -12.613
b1_VPD_raw <-  11.093

# Optimizer settings
trace_flag <- TRUE
rng_seed   <- 1990
polish     <- TRUE
n_max_iter <- 25000

# Clamping (min–max factors)
clamp_Rg <- TRUE
clamp_Ta <- TRUE
clamp_P  <- TRUE

# Metadata / output
metadata <- "Jarvis ET model (min–max factors + VPD + CO2). Universal bundle."
out_file <- "./RData/20251214_jarvis_objects.RData"

# numeric guard (already set in preprocessing, but keep safe fallback)
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

set.seed(rng_seed)

# ------------------------------------------------------------------------------
# Parameter helpers
# ------------------------------------------------------------------------------
get_par_names <- function(mode = c("given","optimize")) {
  mode <- match.arg(mode)
  if (mode == "optimize") {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2","b0_VPD","b1_VPD")
  } else {
    c("Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2")
  }
}

.require_pars <- function(par, needed) {
  miss <- setdiff(needed, names(par))
  if (length(miss)) stop("Missing parameters in 'par': ", paste(miss, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Jarvis g_eff
# ------------------------------------------------------------------------------
jarvis_g_eff <- function(par, d, mode = c("given","optimize")) {
  mode <- match.arg(mode)
  
  needed <- get_par_names(mode)
  .require_pars(par, needed)
  
  Rg_min <- par[["Rg_min"]]; Rg_max <- par[["Rg_max"]]
  Ta_min <- par[["Ta_min"]]; Ta_max <- par[["Ta_max"]]
  P_min  <- par[["P_min"]] ; P_max  <- par[["P_max"]]
  k_CO2  <- par[["k_CO2"]]
  
  # keep bounds ordered
  if (Rg_min > Rg_max) { tmp <- Rg_min; Rg_min <- Rg_max; Rg_max <- tmp }
  if (Ta_min > Ta_max) { tmp <- Ta_min; Ta_min <- Ta_max; Ta_max <- tmp }
  if (P_min  > P_max ) { tmp <- P_min ; P_min  <- P_max ; P_max  <- tmp }
  
  # VPD coefficients
  if (mode == "optimize") {
    b0 <- par[["b0_VPD"]]
    b1 <- par[["b1_VPD"]]
  } else {
    b0 <- b0_VPD_raw
    b1 <- b1_VPD_raw
  }
  
  # spans
  spanRg <- pmax(Rg_max - Rg_min, eps)
  spanTa <- pmax(Ta_max - Ta_min, eps)
  spanP  <- pmax(P_max  - P_min,  eps)
  
  # factors (clamped or unclamped)
  fRg <- if (!use_Rg) 1 else {
    x <- (d$Rg - Rg_min) / spanRg
    if (isTRUE(clamp_Rg)) pmin(1, pmax(0, x)) else x
  }
  
  fTa <- if (!use_Ta) 1 else {
    x <- (d$Ta - Ta_min) / spanTa
    if (isTRUE(clamp_Ta)) pmin(1, pmax(0, x)) else x
  }
  
  fP <- if (!use_P) 1 else {
    x <- (d$P - P_min) / spanP
    if (isTRUE(clamp_P)) pmin(1, pmax(0, x)) else x
  }
  
  fCO2 <- if (!use_CO2) 1 else exp(-pmax(k_CO2, 0) * d$CO2_term)
  fVPD <- if (!use_VPD) 1 else (b0 + b1 * d$inv_sqrt_VPD)
  
  fRg * fTa * fP * fVPD * fCO2
}

# ------------------------------------------------------------------------------
# Model wrapper -> ET
# ------------------------------------------------------------------------------
model_ET_components <- function(par, d, mode = c("given","optimize")) {
  mode <- match.arg(mode)
  g <- jarvis_g_eff(par, d, mode)
  list(
    geff = as.numeric(g),
    ET   = as.numeric(d$K_ET * g)
  )
}

# ------------------------------------------------------------------------------
# Objective
# ------------------------------------------------------------------------------
obj_fun_ET <- function(par, input_df, mode = c("given","optimize")) {
  mode <- match.arg(mode)
  names(par) <- get_par_names(mode)
  pred <- model_ET_components(par, input_df, mode)$ET
  if (!all(is.finite(pred))) return(1e12)
  sqrt(mean((input_df$ET - pred)^2, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Bounds
# ------------------------------------------------------------------------------
par_names <- get_par_names(fVPD_mode)

if (fVPD_mode == "optimize") {
  lower_all <- c(
    Rg_min=-500, Rg_max=50,
    Ta_min=-1000, Ta_max=5,
    P_min=-1000, P_max=400,
    k_CO2=0,
    b0_VPD=-5e4, b1_VPD=-5e4
  )
  upper_all <- c(
    Rg_min=190, Rg_max=5000,
    Ta_min=20,  Ta_max=100,
    P_min=1000, P_max=50000,
    k_CO2=5,
    b0_VPD=5e4, b1_VPD=5e4
  )
} else {
  lower_all <- c(
    Rg_min=-500, Rg_max=50,
    Ta_min=-1000, Ta_max=5,
    P_min=-1000, P_max=400,
    k_CO2=0
  )
  upper_all <- c(
    Rg_min=190, Rg_max=5000,
    Ta_min=20,  Ta_max=100,
    P_min=1000, P_max=50000,
    k_CO2=5
  )
}

lower <- lower_all[par_names]
upper <- upper_all[par_names]

# ------------------------------------------------------------------------------
# Fit (DEoptim + optional local polish)
# ------------------------------------------------------------------------------
D  <- length(par_names)
NP <- max(10L * D, 60L)

compiler::enableJIT(3)
obj_fun_ET          <- compiler::cmpfun(obj_fun_ET)
model_ET_components <- compiler::cmpfun(model_ET_components)
jarvis_g_eff        <- compiler::cmpfun(jarvis_g_eff)

de_fit <- DEoptim::DEoptim(
  fn        = obj_fun_ET,
  lower     = lower,
  upper     = upper,
  input_df  = df_opt,
  mode      = fVPD_mode,
  control   = DEoptim::DEoptim.control(
    NP = NP,
    itermax = n_max_iter,
    reltol = 1e-8,
    trace = trace_flag,
    parallelType = 0
  )
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- par_names
best_val <- de_fit$optim$bestval

if (isTRUE(polish)) {
  polished <- try(
    optim(
      par      = par_hat,
      fn       = obj_fun_ET,
      input_df = df_opt,
      mode     = fVPD_mode,
      method   = "L-BFGS-B",
      lower    = lower,
      upper    = upper,
      control  = list(factr = 1e4, pgtol = 0, maxit = 10000)
    ),
    silent = TRUE
  )
  if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < best_val) {
    par_hat <- polished$par
    names(par_hat) <- par_names
    best_val <- polished$value
  }
}

# ------------------------------------------------------------------------------
# Universal predictor (KEY)
#   - recomputes inv_sqrt_VPD + K_ET from the *current* VPD/Ta
#   - uses recompute_K_ET() from preprocessing as the single source of truth
# ------------------------------------------------------------------------------
predict_jarvis <- function(newdata) {
  d <- newdata |>
    mutate(
      VPD = pmax(VPD, eps),
      inv_sqrt_VPD = 1 / sqrt(VPD),
      K_ET = recompute_K_ET(cur_data())
    )
  
  model_ET_components(par_hat, d, mode = fVPD_mode)
}

# ------------------------------------------------------------------------------
# Predictions on full df
# ------------------------------------------------------------------------------
pred <- predict_jarvis(df)

output_df <- df |>
  mutate(
    g_eff_pred = pred$geff,
    ET_pred    = pred$ET,
    ET_resid   = ET - ET_pred
  )

RMSE <- sqrt(mean((output_df$ET - output_df$ET_pred)^2, na.rm = TRUE))
R2   <- suppressWarnings(cor(output_df$ET, output_df$ET_pred, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# ------------------------------------------------------------------------------
# Parameter table (SEs via Jacobian)
# ------------------------------------------------------------------------------
residuals_fun <- function(th, data) {
  names(th) <- names(par_hat)
  comps <- model_ET_components(th, data, mode = fVPD_mode)
  r <- data$ET - comps$ET
  r[!is.finite(r)] <- 0
  r
}

J_all   <- numDeriv::jacobian(function(th) residuals_fun(th, df_opt), par_hat)
SSE_hat <- sum(residuals_fun(par_hat, df_opt)^2)
n <- nrow(df_opt); p <- length(par_hat)
sigma2 <- SSE_hat / max(1, n - p)

XtX    <- crossprod(J_all)
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

# ------------------------------------------------------------------------------
# Save UNIVERSAL bundle
# ------------------------------------------------------------------------------
jarvis_bundle <- list(
  metadata = metadata,
  model_id = "jarvis",
  
  # universal interface (for the shared permutation script)
  df_opt   = df_opt,
  df       = df,
  pred_fun = predict_jarvis,
  recompute_K_ET = recompute_K_ET,
  eps = eps,
  
  # model content
  par_hat   = par_hat,
  fVPD_mode = fVPD_mode,
  output_df = output_df,
  RMSE = RMSE,
  R2   = R2,
  
  # switches
  use_P = use_P, use_Rg = use_Rg, use_Ta = use_Ta, use_CO2 = use_CO2, use_VPD = use_VPD,
  clamp_Rg = clamp_Rg, clamp_Ta = clamp_Ta, clamp_P = clamp_P,
  
  get_par_names = get_par_names,
  .require_pars = .require_pars,
  
  # optional internals (nice to keep)
  jarvis_g_eff = jarvis_g_eff,
  model_ET_components = model_ET_components
)

save(jarvis_bundle, file = out_file)
message("Jarvis model fitted and universal bundle saved: ", out_file)
