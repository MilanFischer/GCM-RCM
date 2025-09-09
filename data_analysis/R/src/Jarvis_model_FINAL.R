eps <- 1e-6  # guard

# ======================
# SWITCHBOARD
# ======================
use_P   <- TRUE
use_Rg  <- TRUE
use_Ta  <- TRUE
use_VPD <- TRUE
use_CO2 <- TRUE

# Standardize (z-score) drivers?
standardize_drivers  <- TRUE   # P, Rg, Ta will follow this flag
standardize_invVPD   <- TRUE   # invVPD standardization can be toggled separately

# ---- CHOOSE HOW TO HANDLE fVPD ----
# "given": fVPD is supplied (not optimized)
# "optimize": fVPD coefficients are optimized (a_VPD, b_VPD)
fVPD_mode <- "optimize"     # <- set to "given" or "optimize"

# Trace DEoptim?
trace_flag <- FALSE

# If fVPD_mode == "given", you can provide either:
#   (A) RAW-scale coefficients for fVPD = b0 + b1 * invVPD_raw
#   (B) a full precomputed per-row vector fVPD_given_vec (advanced; overrides A)
# IMPORTANT: Keep these RAW and immutable. We'll convert to z-scale on-the-fly if needed.
b0_VPD_raw <- -6.909279   # from prior fit on RAW invVPD
b1_VPD_raw <-  8.282555
fVPD_given_vec <- NULL    # or supply a numeric vector length nrow(Data_to_plot_II)

# ---- CO2 term ----
dln_CO2 <- log(CO2_2076_2100_RCP85 / CO2_1981_2005)
Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    CO2_term = if_else(grepl("CMIP", model) & PERIOD == "2076_2100", dln_CO2, 0)
  )

# ======================
# Helpers: scalers and drivers
# ======================

compute_scalers <- function(df) {
  invVPD_raw <- 1 / sqrt(pmax(df$VPD, eps))
  mu <- c(
    P      = mean(df$P,    na.rm = TRUE),
    Rg     = mean(df$Rg,   na.rm = TRUE),
    Ta     = mean(df$Ta,   na.rm = TRUE),
    invVPD = mean(invVPD_raw, na.rm = TRUE)
  )
  sd <- c(
    P      = sd(df$P,    na.rm = TRUE),
    Rg     = sd(df$Rg,   na.rm = TRUE),
    Ta     = sd(df$Ta,   na.rm = TRUE),
    invVPD = sd(invVPD_raw, na.rm = TRUE)
  )
  sd[!is.finite(sd) | sd < eps] <- 1
  list(mu = mu, sd = sd)
}

build_drivers <- function(df, scalers, standardize = FALSE, standardize_invVPD = TRUE) {
  invVPD_raw <- 1 / sqrt(pmax(df$VPD, eps))
  
  P_col  <- if (!standardize) df$P  else (df$P  - scalers$mu[["P"]])  / scalers$sd[["P"]]
  Rg_col <- if (!standardize) df$Rg else (df$Rg - scalers$mu[["Rg"]]) / scalers$sd[["Rg"]]
  Ta_col <- if (!standardize) df$Ta else (df$Ta - scalers$mu[["Ta"]]) / scalers$sd[["Ta"]]
  
  invVPD_col <- if (!standardize_invVPD) invVPD_raw
  else (invVPD_raw - scalers$mu[["invVPD"]]) / scalers$sd[["invVPD"]]
  
  tibble(P = P_col, Rg = Rg_col, Ta = Ta_col, invVPD = invVPD_col)
}

# Convert RAW (b0,b1) for fVPD = b0 + b1*invVPD_raw to the scale in use
convert_vpd_coefs <- function(b0_raw, b1_raw, scalers, standardize_invVPD) {
  if (isTRUE(standardize_invVPD)) {
    # invVPD_z = (invVPD_raw - mu)/sd -> f = (b0 + b1*mu) + (b1*sd)*invVPD_z
    b0_star <- b0_raw + b1_raw * scalers$mu[["invVPD"]]
    b1_star <- b1_raw * scalers$sd[["invVPD"]]
    c(b0 = b0_star, b1 = b1_star)
  } else {
    c(b0 = b0_raw, b1 = b1_raw)
  }
}

# ======================
# Jarvis g_eff (mm/s) with regression-style modifiers
# f(x) = a + b*x
# Supports both fVPD modes
# ======================
jarvis_g_eff <- function(drv, CO2_term,
                         a_P, b_P,
                         a_Rg, b_Rg,
                         a_Ta, b_Ta,
                         # VPD pieces:
                         fVPD_mode = c("given", "optimize"),
                         a_VPD = NA_real_, b_VPD = NA_real_,           # only used if optimize
                         b0_VPD = NA_real_, b1_VPD = NA_real_,         # only used if given & no vector
                         fVPD_given_vec = NULL,                        # only used if given & supplied
                         k_CO2 = 0) {
  
  fVPD_mode <- match.arg(fVPD_mode)
  
  # P, Rg, Ta terms (affine)
  fP   <- if (use_P)   (a_P   + b_P   * drv$P)    else 1
  fRg  <- if (use_Rg)  (a_Rg  + b_Rg  * drv$Rg)   else 1
  fTa  <- if (use_Ta)  (a_Ta  + b_Ta  * drv$Ta)   else 1
  
  # VPD term (affine)
  if (!use_VPD) {
    fVPD <- 1
  } else if (fVPD_mode == "optimize") {
    fVPD <- (a_VPD + b_VPD * drv$invVPD)
  } else {  # "given"
    if (!is.null(fVPD_given_vec)) {
      fVPD <- fVPD_given_vec
    } else {
      fVPD <- (b0_VPD + b1_VPD * drv$invVPD)
    }
  }
  
  # CO2 modifier
  fCO2 <- if (use_CO2) exp(-pmax(k_CO2, 0) * CO2_term) else 1
  
  # Final Jarvis multiplicative form
  g_eff <- fP * fRg * fTa * fVPD * fCO2
  g_eff
}

# ======================
# Model wrapper -> ET
# par order depends on fVPD_mode
#   If "optimize": c(a_P,b_P, a_Rg,b_Rg, a_Ta,b_Ta, a_VPD,b_VPD, k_CO2)  (9 params)
#   If "given":    c(a_P,b_P, a_Rg,b_Rg, a_Ta,b_Ta, k_CO2)               (7 params)
# ======================
model_ET_components <- function(par, data, drv,
                                fVPD_mode,
                                b0_VPD, b1_VPD, fVPD_given_vec) {
  
  # unpack depending on mode
  if (fVPD_mode == "optimize") {
    a_P   <- par[1];  b_P   <- par[2]
    a_Rg  <- par[3];  b_Rg  <- par[4]
    a_Ta  <- par[5];  b_Ta  <- par[6]
    a_VPD <- par[7];  b_VPD <- par[8]
    k_CO2 <- par[9]
  } else { # "given"
    a_P   <- par[1];  b_P   <- par[2]
    a_Rg  <- par[3];  b_Rg  <- par[4]
    a_Ta  <- par[5];  b_Ta  <- par[6]
    k_CO2 <- par[7]
    a_VPD <- NA_real_; b_VPD <- NA_real_
  }
  
  g_eff_pred <- jarvis_g_eff(
    drv      = drv,
    CO2_term = data$CO2_term,
    a_P=a_P, b_P=b_P,
    a_Rg=a_Rg, b_Rg=b_Rg,
    a_Ta=a_Ta, b_Ta=b_Ta,
    fVPD_mode = fVPD_mode,
    a_VPD=a_VPD, b_VPD=b_VPD,
    b0_VPD=b0_VPD, b1_VPD=b1_VPD,
    fVPD_given_vec = fVPD_given_vec,
    k_CO2=k_CO2
  )
  
  VPDg <- pmax(data$VPD, eps)
  
  # --- ET conversion ---
  ET_pred <- (rhoAir * CpAir / gamma) * VPDg /
    (1 / (g_eff_pred / 1000)) *                    # s/m
    1 / ((2.501e6 - 2361 * data$Ta) / 1e6) *       # LE conversion (W m^-2)
    (3600 * 24) / 1e6 * 365.25                     # mm/yr
  
  list(ET_pred = ET_pred, g_eff_pred = g_eff_pred)
}

# ---- Objective ----
obj_fun_ET <- function(par, data, drv,
                       fVPD_mode, b0_VPD, b1_VPD, fVPD_given_vec) {
  comps <- model_ET_components(par, data, drv,
                               fVPD_mode=fVPD_mode,
                               b0_VPD=b0_VPD, b1_VPD=b1_VPD,
                               fVPD_given_vec=fVPD_given_vec)
  if (!all(is.finite(comps$ET_pred)) || any(!is.finite(comps$g_eff_pred))) return(1e12)
  if (any(comps$g_eff_pred <= 0)) return(1e12 + 1e6 * sum(comps$g_eff_pred <= 0))
  res <- data$ET - comps$ET_pred
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

# Compute scalers on the FIT subset (df_opt) ONCE
scalers <- compute_scalers(df_opt)

# Build drivers (using the SAME scalers everywhere)
drv_opt  <- build_drivers(df_opt,          scalers,
                          standardize = standardize_drivers,
                          standardize_invVPD = standardize_invVPD)
drv_full <- build_drivers(Data_to_plot_II, scalers,
                          standardize = standardize_drivers,
                          standardize_invVPD = standardize_invVPD)

# If fVPD is "given" and no custom vector: derive the correct coefs for the current scaling
if (fVPD_mode == "given" && is.null(fVPD_given_vec)) {
  vpd_coefs_fit <- convert_vpd_coefs(b0_VPD_raw, b1_VPD_raw, scalers, standardize_invVPD)
} else {
  vpd_coefs_fit <- c(b0 = NA_real_, b1 = NA_real_) # ignored in optimize mode or custom vec
}

# ---- Bounds depend on mode ----
if (fVPD_mode == "optimize") {
  par_names <- c("a_P","b_P","a_Rg","b_Rg","a_Ta","b_Ta","a_VPD","b_VPD","k_CO2")
  lower <- c(a_P=-2, b_P=-2, a_Rg=-2, b_Rg=-2, a_Ta=-2, b_Ta=-2, a_VPD=-2, b_VPD=-2, k_CO2=0)
  upper <- c(a_P= 2, b_P= 2, a_Rg= 2, b_Rg= 2, a_Ta= 2, b_Ta= 2, a_VPD= 2, b_VPD= 2, k_CO2=5)
} else {
  par_names <- c("a_P","b_P","a_Rg","b_Rg","a_Ta","b_Ta","k_CO2")
  lower <- c(a_P=-2, b_P=-2, a_Rg=-2, b_Rg=-2, a_Ta=-2, b_Ta=-2, k_CO2=0)
  upper <- c(a_P= 2, b_P= 2, a_Rg= 2, b_Rg= 2, a_Ta= 2, b_Ta= 2, k_CO2=5)
}

set.seed(153)
de_fit <- DEoptim(
  fn      = obj_fun_ET,
  lower   = lower,
  upper   = upper,
  data    = df_opt,
  drv     = drv_opt,
  # extra args forwarded to obj_fun_ET:
  fVPD_mode = fVPD_mode,
  b0_VPD = unname(vpd_coefs_fit["b0"]),
  b1_VPD = unname(vpd_coefs_fit["b1"]),
  fVPD_given_vec = if (is.null(fVPD_given_vec)) NULL else fVPD_given_vec[seq_len(nrow(df_opt))],
  control = DEoptim.control(NP = 140, itermax = 1000, reltol = 1e-8, trace = trace_flag)
)

par_hat <- de_fit$optim$bestmem
names(par_hat) <- par_names
SSE_hat <- de_fit$optim$bestval

# ---- Optional polish (bounded; keeps you inside the box) ----
polished <- try(
  optim(par = par_hat, fn = obj_fun_ET, data = df_opt, drv = drv_opt,
        fVPD_mode = fVPD_mode,
        b0_VPD = unname(vpd_coefs_fit["b0"]),
        b1_VPD = unname(vpd_coefs_fit["b1"]),
        fVPD_given_vec = if (is.null(fVPD_given_vec)) NULL else fVPD_given_vec[seq_len(nrow(df_opt))],
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(reltol = 1e-12, maxit = 1000)),
  silent = TRUE
)
if (!inherits(polished, "try-error") && is.finite(polished$value) && polished$value < SSE_hat) {
  par_hat <- polished$par
  names(par_hat) <- par_names
  SSE_hat <- polished$value
}

# ---- Predict on full data ----
comps_full <- model_ET_components(par_hat, Data_to_plot_II, drv_full,
                                  fVPD_mode=fVPD_mode,
                                  b0_VPD=unname(vpd_coefs_fit["b0"]),
                                  b1_VPD=unname(vpd_coefs_fit["b1"]),
                                  fVPD_given_vec=fVPD_given_vec)

Data_to_plot_II <- Data_to_plot_II |>
  mutate(
    g_eff_predicted = comps_full$g_eff_pred,  # mm/s
    ET_predicted    = comps_full$ET_pred,
    ET_resids       = ET - ET_predicted
  )

# ---- Fit metrics ----
pred_fit <- model_ET_components(par_hat, df_opt, drv_opt,
                                fVPD_mode=fVPD_mode,
                                b0_VPD=unname(vpd_coefs_fit["b0"]),
                                b1_VPD=unname(vpd_coefs_fit["b1"]),
                                fVPD_given_vec = if (is.null(fVPD_given_vec)) NULL
                                else fVPD_given_vec[seq_len(nrow(df_opt))])$ET_pred

obs_fit  <- df_opt$ET
RMSE     <- sqrt(mean((obs_fit - pred_fit)^2, na.rm = TRUE))
R2       <- suppressWarnings(cor(obs_fit, pred_fit, use = "complete.obs")^2)
message(sprintf("In-sample fit: RMSE = %.3f, R^2 = %.3f", RMSE, R2))

# ---- Standard errors (active set depends on mode & switches) ----
residuals_fun <- function(th, data, drv) {
  comps <- model_ET_components(th, data, drv,
                               fVPD_mode=fVPD_mode,
                               b0_VPD=unname(vpd_coefs_fit["b0"]),
                               b1_VPD=unname(vpd_coefs_fit["b1"]),
                               fVPD_given_vec = if (is.null(fVPD_given_vec)) NULL
                               else fVPD_given_vec[seq_len(nrow(df_opt))])
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
with(Data_to_plot_II, {
  plot(ET, ET_predicted, xlab = "Observed ET", ylab = "Predicted ET")
  abline(0, 1, col = "red")
})
