eps <- 1e-6  # guard

# ======================
# SWITCHBOARD
# ======================
use_P   <- TRUE
use_Rg  <- TRUE
use_Ta  <- TRUE
use_CO2 <- TRUE
use_VPD <- TRUE

input_df_raw <- Data_to_plot_II

# Select only the variables that are used for the purpose of g_eff modeling
input_df <- input_df_raw |> 
  select(g_eff, ET, Rg, Ta, P, VPD, PERIOD, model, label, 
         color, fill, border, shape, linetype)

g_eff_max <- input_df$g_eff |> max(na.rm=T) |> ceiling()

# ---- CHOOSE HOW TO HANDLE fVPD ----
# "given": fVPD is supplied (not optimized)
# "optimize": fVPD coefficients are optimized (a_VPD, b_VPD)
fVPD_mode <- "given"     # <- set to "given" or "optimize"

# Trace DEoptim?
trace_flag <- TRUE

# If fVPD_mode == "given", you must provide the coefficients
# RAW-scale coefficients for fVPD = b0 + b1 * invVPD_raw
b0_VPD_raw <- -6.909279   # from prior fit on RAW invVPD
b1_VPD_raw <-  8.282555

# ---- CO2 term ----
dln_CO2 <- log(CO2_2076_2100_RCP85 / CO2_1981_2005)
input_df <- input_df |>
  mutate(
    CO2_term = if_else(grepl("CMIP", model) & PERIOD == "2076_2100", dln_CO2, 0)
  )

# ==========================
# Initial set of parameters
# ==========================
make_par0 <- function(fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
    par0 <- c(
      g_eff_max  =   14,
      Rg_min     =    0, Rg_max = 200,
      Ta_min     =  -10, Ta_max = 50,
      P_min      =    0, P_max  = 2000,
      k_CO2      =    0.1,
      b0_VPD     =  if (fVPD_mode == "optimize") -7 else if (fVPD_mode == "given") b0_VPD_raw,
      b1_VPD     =  if (fVPD_mode == "optimize") 8  else if (fVPD_mode == "given") b1_VPD_raw 
    )
  return(par0)
}

# -----------------------
# Helper: validate names
# -----------------------
.require_pars <- function(par, needed) {
  miss <- setdiff(needed, names(par))
  if (length(miss)) stop("Missing parameters in 'par': ", paste(miss, collapse=", "))
}

# -----------------------------
# Helper: get parameters names
# -----------------------------
get_par_names <- function(fVPD_mode = c("given","optimize")) {
  fVPD_mode <- match.arg(fVPD_mode)
  if (fVPD_mode == "optimize") {
    c("g_eff_max","Rg_min","Rg_max","Ta_min","Ta_max",
      "P_min","P_max","k_CO2","b0_VPD","b1_VPD")
  } else {
    c("g_eff_max","Rg_min","Rg_max","Ta_min","Ta_max",
      "P_min","P_max","k_CO2")
  }
}

# ======================
# Jarvis g_eff (mm/s)
# ======================
jarvis_g_eff <- function(par, input_df, fVPD_mode  = c("given", "optimize")) {
  
  # ensure fVPD_mode is matched properly
  fVPD_mode <- match.arg(fVPD_mode)
  
  # base parameters always needed
  # only require VPD coefs if we're optimizing them
  base_needed <- if (fVPD_mode == "optimize") {
    c("g_eff_max","Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2","b0_VPD","b1_VPD")
  } else {
    c("g_eff_max","Rg_min","Rg_max","Ta_min","Ta_max","P_min","P_max","k_CO2")
  }

  .require_pars(par, base_needed)
  
  # unpack (explicit is clearer than attach)
  g_eff_max <- par[["g_eff_max"]]
  Rg_min    <- par[["Rg_min"]];  Rg_max <- par[["Rg_max"]]
  Ta_min    <- par[["Ta_min"]];  Ta_max <- par[["Ta_max"]]
  P_min     <- par[["P_min"]] ;  P_max  <- par[["P_max"]]
  k_CO2     <- par[["k_CO2"]]
  
  # VPD coefs: from par if optimizing, otherwise use the provided "given" ones
  if (fVPD_mode == "optimize") {
    b0_VPD <- par[["b0_VPD"]]
    b1_VPD <- par[["b1_VPD"]]
  } else {
    b0_VPD <- b0_VPD_raw
    b1_VPD <- b1_VPD_raw
  }
  
  # pull inputs
  Rg           <- input_df |> pull(Rg)
  Ta           <- input_df |> pull(Ta)
  P            <- input_df |> pull(P)
  CO2_term     <- input_df |> pull(CO2_term)
  VPD          <- input_df |> pull(VPD)
  inv_sqrt_VPD <- 1 / sqrt(VPD)
  g_eff        <- input_df |> pull(g_eff)
  ET           <- input_df |> pull(ET)
  
  # P, Rg, Ta terms (affine)
  # fRg  <- if (use_Rg)  pmax(0, pmin(1, (Rg - Rg_min) / (Rg_max - Rg_min))) else 1
  # fTa  <- if (use_Ta)  pmax(0, pmin(1, (Ta - Ta_min) / (Ta_max - Ta_min))) else 1
  # fP   <- if (use_P)   pmax(0, pmin(1, (P - P_min)   / (P_max  - P_min)))  else 1
  
  fRg  <- if (use_Rg)  (Rg - Rg_min) / (Rg_max - Rg_min) else 1
  fTa  <- if (use_Ta)  (Ta - Ta_min) / (Ta_max - Ta_min) else 1
  fP   <- if (use_P)   (P - P_min)   / (P_max  - P_min)  else 1
  
  # CO2 modifier
  fCO2 <- if (use_CO2) exp(-pmax(k_CO2, 0) * CO2_term) else 1
  
  # VPD modifier – need to decide if a1_VPD, perhaps only if g_eff_max is fixed
  # Or a1_VPD * (b0 + b1 / sqrt(VPD)) / max(a1_VPD * (b0 + b1 / sqrt(VPD)))
  # fVPD <- if (use_VPD) a1_VPD * (b0 + b1 * inv_sqrt_VPD) else 1

  # VPD term
  if (!use_VPD) {
    fVPD <- 1
  } else {
    fVPD <- (b0_VPD + b1_VPD * inv_sqrt_VPD) / max( (b0_VPD + b1_VPD * inv_sqrt_VPD), na.rm = TRUE)
  }

  # Final Jarvis multiplicative form – all function have maximum of 1 and minimum of 0
  g_eff_pred <- g_eff_max * fP * fRg * fTa * fVPD * fCO2
  
  return(g_eff_pred)
}

# ======================
# Model wrapper -> ET
# ======================
model_ET_components <- function(par, input_df, fVPD_mode  = c("given", "optimize")) {
  
  fVPD_mode <- match.arg(fVPD_mode)

  # predict g_eff
  g_eff_pred <- jarvis_g_eff(par, input_df, fVPD_mode )
  
  # pull inputs
  Ta           <- input_df |> pull(Ta)
  VPD          <- input_df |> pull(VPD)
  
  # --- ET conversion ---
  ET_pred <- (rhoAir * CpAir / gamma) * VPD /
    (1 / (g_eff_pred / 1000)) *                # s/m
    1 / ((2.501e6 - 2361 * Ta) / 1e6) *        # LE conversion (W m^-2)
    (3600 * 24) / 1e6 * 365.25                 # mm/yr
  
  list(ET_pred = ET_pred, g_eff_pred = g_eff_pred)
}

# ---- Objective ----
obj_fun_ET <- function(par, input_df, fVPD_mode  = c("given", "optimize")) {
  
  names(par) <- get_par_names(fVPD_mode) 
  fVPD_mode <- match.arg(fVPD_mode)
  
  # pull inputs
  ET <- input_df |> pull(ET)
  
  comps <- model_ET_components(par, input_df, fVPD_mode)
  if (!all(is.finite(comps$ET_pred)) || any(!is.finite(comps$g_eff_pred))) return(1e12)
  if (any(comps$g_eff_pred <= 0)) return(1e12 + 1e6 * sum(comps$g_eff_pred <= 0))
  res <- ET - comps$ET_pred
  if (!all(is.finite(res))) return(1e12)
  sum(res^2)
}

# ---- Fit subset ----
df_opt <- input_df |>
  filter(is.finite(ET),
         is.finite(P),
         is.finite(Rg),
         is.finite(VPD),
         is.finite(Ta))

par0 <- make_par0("optimize")

# ---- Bounds depend on mode ----
if (fVPD_mode == "optimize") {
  par_names <- get_par_names(fVPD_mode)
  lower <- c(g_eff_max=1, Rg_min=0, Rg_max=0, Ta_min=-10,Ta_max=-10, P_min=0,P_max=0, k_CO2=0, b0_VPD=-20, b1_VPD=-20)
  upper <- c(g_eff_max=100, Rg_min=500, Rg_max=500, Ta_min=50,Ta_max=50, P_min=3000,P_max=3000, k_CO2=5, b0_VPD=20, b1_VPD=20)
} else {
  par_names <- get_par_names(fVPD_mode)
  lower <- c(g_eff_max=1, Rg_min=0, Rg_max=0, Ta_min=-10,Ta_max=-10, P_min=0,P_max=0, k_CO2=0)
  upper <- c(g_eff_max=100, Rg_min=500, Rg_max=500, Ta_min=50,Ta_max=50, P_min=3000,P_max=3000, k_CO2=5)
}

set.seed(153)
de_fit <- DEoptim(
  fn       = obj_fun_ET,
  lower    = lower,
  upper    = upper,
  input_df = df_opt,
  # extra args forwarded to obj_fun_ET:
  fVPD_mode = fVPD_mode,
  control = DEoptim.control(NP = 140, itermax = 1000, reltol = 1e-8, trace = trace_flag)
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- par_names
SSE_hat <- de_fit$optim$bestval

# ---- Optional polish (bounded; keeps you inside the box) ----
# ---- Optional polish (bounded; keeps you inside the box) ----
polished <- try(
  optim(
    par     = par_hat,
    fn      = obj_fun_ET,
    input_df = df_opt,
    fVPD_mode = fVPD_mode,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(
      factr = 1e4,    # replaces reltol for L-BFGS-B
      pgtol = 0,      # gradient tolerance (0 = very strict)
      maxit = 1000
    )
  ),
  silent = TRUE
)
if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < SSE_hat) {
  par_hat <- polished$par
  names(par_hat) <- par_names
  SSE_hat <- polished$value
}

# ---- Predict on full data ----
comps_full <- model_ET_components(par_hat, input_df, fVPD_mode=fVPD_mode)

input_df <- input_df |>
  mutate(
    g_eff_predicted = comps_full$g_eff_pred,  # mm/s
    ET_predicted    = comps_full$ET_pred,     # mm/yr
    ET_resids       = ET - ET_predicted       # mm/yr
  )

# ---- Fit metrics ----
pred_fit <- model_ET_components(par_hat, input_df, fVPD_mode=fVPD_mode)$ET_pred

obs_fit  <- input_df$ET
RMSE     <- sqrt(mean((obs_fit - pred_fit)^2, na.rm = TRUE))
R2       <- suppressWarnings(cor(obs_fit, pred_fit, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# ---- Standard errors (active set depends on mode & switches) ----
residuals_fun <- function(th, data, drv) {
  comps <- model_ET_components(th, data, fVPD_mode=fVPD_mode)
  r <- data$ET - comps$ET_pred
  r[!is.finite(r)] <- 0
  r
}
J_all <- jacobian(function(th) residuals_fun(th, df_opt, drv_opt), par_hat)

k_used <- if (isTRUE(use_CO2)) max(unname(par_hat["k_CO2"]), 0) else 0
if (fVPD_mode == "optimize") {
  active_flags <- c(
    a_P=use_P, b_P=use_P,
    a_Rg=use_Rg, b_Rg=use_Rg,
    a_Ta=use_Ta, b_Ta=use_Ta,
    a_VPD=use_VPD, b_VPD=use_VPD,
    k_CO2=(use_CO2 && k_used > 0)
  )
} else {
  active_flags <- c(
    a_P=use_P, b_P=use_P,
    a_Rg=use_Rg, b_Rg=use_Rg,
    a_Ta=use_Ta, b_Ta=use_Ta,
    k_CO2=(use_CO2 && k_used > 0)
  )
}
active_idx <- which(active_flags)

SSE_hat <- sum(residuals_fun(par_hat, df_opt, drv_opt)^2)
n <- nrow(df_opt); p <- length(par_hat)
sigma2 <- SSE_hat / (n - length(active_idx))

se <- rep(NA_real_, length(par_hat)); names(se) <- names(par_hat)
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
    !(names(par_hat) %in% names(active_flags)) ~ "n/a (not estimated in this mode)",
    !active_flags[match(names(par_hat), names(active_flags))] ~ "disabled by switch",
    names(par_hat) == "k_CO2" & k_used == 0 ~ "inactive (no CO2 effect at optimum)",
    TRUE ~ NA_character_
  )
)
print(param_table)

# ---- Diagnostic ----
with(input_df, {
  plot(input_df$ET, input_df$ET_predicted, xlab = "Observed ET", ylab = "Predicted ET")
  abline(0, 1, col = "red")
})

jarvis_out <- input_df
