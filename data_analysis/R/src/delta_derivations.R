suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(mgcv)        # for GAM fits (predicted mode)
  library(ggnewscale)  # only needed if you layer into an existing ggplot 'p9'
})

# ===========================
# Guards
# ===========================
if (!exists("Data_to_plot", inherits = TRUE) || is.null(Data_to_plot[["abs"]])) {
  stop("Missing Data_to_plot[['abs']] in your environment.", call. = FALSE)
}
stopifnot(all(c("rhoAir","CpAir","gamma") %in% ls(.GlobalEnv)))
if (!exists("jarvis_bundle", inherits = TRUE) || is.null(jarvis_bundle$par_hat)) {
  stop("Missing jarvis_bundle$par_hat in your environment.", call. = FALSE)
}

# ===========================
# Shared helpers
# ===========================
eps <- 1e-6

K_ET_fun <- function(VPD, Ta, rhoAir, CpAir, gamma){
  (rhoAir * CpAir / gamma) * VPD * (1/1000) *
    1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

# Fixed Jarvis parameters
par_fixed <- c(
  Rg_min    = as.numeric(jarvis_bundle$par_hat["Rg_min"]),
  Rg_max    = as.numeric(jarvis_bundle$par_hat["Rg_max"]),
  Ta_min    = as.numeric(jarvis_bundle$par_hat["Ta_min"]),
  Ta_max    = as.numeric(jarvis_bundle$par_hat["Ta_max"]),
  P_min     = as.numeric(jarvis_bundle$par_hat["P_min"]),
  P_max     = as.numeric(jarvis_bundle$par_hat["P_max"]),
  k_CO2     = as.numeric(jarvis_bundle$par_hat["k_CO2"]),
  b0_VPD    = as.numeric(jarvis_bundle$par_hat["b0_VPD"]),
  b1_VPD    = as.numeric(jarvis_bundle$par_hat["b1_VPD"])
)

use_P <- TRUE; use_Rg <- TRUE; use_Ta <- TRUE; use_CO2 <- TRUE; use_VPD <- TRUE
clamp_Rg <- TRUE; clamp_Ta <- TRUE; clamp_P <- TRUE

jarvis_g_eff_fixed <- function(par, input_df) {
  Rg_min <- par[["Rg_min"]];  Rg_max <- par[["Rg_max"]]
  Ta_min <- par[["Ta_min"]];  Ta_max <- par[["Ta_max"]]
  P_min  <- par[["P_min"]];   P_max  <- par[["P_max"]]
  k_CO2  <- par[["k_CO2"]]
  b0_VPD <- par[["b0_VPD"]];  b1_VPD <- par[["b1_VPD"]]
  
  Rg <- input_df$Rg; Ta <- input_df$Ta; P <- input_df$P
  
  spanRg <- pmax(Rg_max - Rg_min, eps)
  spanTa <- pmax(Ta_max - Ta_min, eps)
  spanP  <- pmax(P_max  - P_min , eps)
  
  fRg <- if (use_Rg && clamp_Rg) pmin(1, pmax(0, (Rg - Rg_min)/spanRg)) else 1
  fTa <- if (use_Ta && clamp_Ta) pmin(1, pmax(0, (Ta - Ta_min)/spanTa)) else 1
  fP  <- if (use_P  && clamp_P ) pmin(1, pmax(0, (P  - P_min )/spanP )) else 1
  
  # CO2 optional (neutral if no CO2_term provided)
  fCO2 <- if (use_CO2 && "CO2_term" %in% names(input_df)) exp(-pmax(k_CO2,0) * input_df$CO2_term) else 1
  
  if (!use_VPD) {
    fVPD <- 1
  } else {
    inv_sqrt_v <- if ("inv_sqrt_VPD" %in% names(input_df)) input_df$inv_sqrt_VPD else 1 / sqrt(pmax(input_df$VPD, eps))
    raw  <- b0_VPD + b1_VPD * inv_sqrt_v
    fVPD <- pmax(raw, 0)
  }
  
  fP * fRg * fTa * fVPD * fCO2
}

predict_ET_fixed <- function(par, input_df) {
  g_eff <- jarvis_g_eff_fixed(par, input_df)
  K_ET  <- if ("K_ET" %in% names(input_df)) input_df$K_ET else K_ET_fun(input_df$VPD, input_df$Ta, rhoAir, CpAir, gamma)
  ET    <- K_ET * g_eff
  list(g_eff = g_eff, ET = ET)
}

# Log-linear helper (a,b,C)
fit_log_abC <- function(y, vpd, name="Y") {
  d <- data.frame(y = y, vpd = vpd)
  fit <- lm(log(y) ~ log(vpd) + vpd, data = d)
  list(
    fit = fit,
    a   = unname(coef(fit)[["log(vpd)"]]),
    b   = unname(coef(fit)[["vpd"]]),
    C   = unname(coef(fit)[["(Intercept)"]]),
    pred_fun = function(v){ exp(coef(fit)[["log(vpd)"]]*log(v) + coef(fit)[["vpd"]]*v + coef(fit)[["(Intercept)"]]) },
    name = name
  )
}

# ===========================
# Base data
# ===========================
Data_to_plot$abs <- Data_to_plot$abs |> mutate(.row_id = row_number())

df_abs <- Data_to_plot$abs |> 
  transmute(
    .row_id,
    VPD = as.numeric(VPD),
    AI  = as.numeric(AI_FAO56_alfalfa),
    PET = as.numeric(ETo_FAO56_alfalfa),
    P   = as.numeric(P),
    Rg  = as.numeric(Rg),
    Ta  = as.numeric(Ta),
    ET  = as.numeric(ET)
  )

# For observed pipeline
df_obs_ok <- df_abs |>
  filter(is.finite(VPD), is.finite(P), is.finite(Rg), is.finite(Ta), is.finite(ET),
         VPD > 0, P > 0)

# For predicted pipeline
df_pred_ok <- df_abs |>
  filter(is.finite(VPD), VPD > 0,
         is.finite(AI),  AI  > 0,
         is.finite(PET), PET > 0,
         is.finite(P),   P   > 0,
         is.finite(Rg),  is.finite(Ta), is.finite(ET)) |>
  mutate(EI_obs = ET / P)

stopifnot(nrow(df_obs_ok) > 0, nrow(df_pred_ok) > 0)

# ==================
# OBSERVED pipeline
# ==================
df0 <- df_abs |> 
  transmute(.row_id, VPD, PET, P) |>
  filter(!is.na(VPD), !is.na(PET), !is.na(P), VPD > 0, PET > 0, P > 0) |>
  mutate(AI_obs = PET / P)

fit_PET <- lm(log(PET) ~ log(VPD) + VPD, data = df0)
a_PET   <- coef(fit_PET)[["log(VPD)"]]
b_PET   <- coef(fit_PET)[["VPD"]]
C_PET   <- coef(fit_PET)[["(Intercept)"]]

df0 <- df0 |>
  mutate(PET_hat = exp(a_PET*log(VPD) + b_PET*VPD + C_PET),
         AI_hat  = PET_hat / P)  # AI consistent with observed P

# Inputs with observed Rg/Ta/P for Jarvis
df_obs_in <- df_obs_ok |>
  transmute(
    .row_id,
    ET_obs = ET,
    VPD, P, Rg, Ta,
    CO2_term = 0,
    inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
    K_ET = K_ET_fun(VPD, Ta, rhoAir, CpAir, gamma)
  )

preds_obs <- predict_ET_fixed(par_fixed, df_obs_in)

df_out <- df_obs_in |>
  mutate(
    g_eff_pred = preds_obs$g_eff,
    ET_pred    = preds_obs$ET,
    resid      = ET_obs - ET_pred,
    EI_obs     = ET_obs / P,
    EI_pred    = ET_pred / P
  )

RMSE_obs <- sqrt(mean(df_out$resid^2, na.rm = TRUE))
R2_obs   <- suppressWarnings(cor(df_out$ET_obs, df_out$ET_pred, use = "complete.obs")^2)
cat(sprintf("\n[Observed] Jarvis with OBSERVED P/Rg/Ta: RMSE = %.3f, R^2 = %.3f\n", RMSE_obs, R2_obs))

# ===================
# PREDICTED pipeline
# ===================
m_AI <- fit_log_abC(df_pred_ok$AI, df_pred_ok$VPD, "AI")
m_P  <- fit_log_abC(df_pred_ok$P , df_pred_ok$VPD, "P")
m_EI_dir <- fit_log_abC(df_pred_ok$EI_obs, df_pred_ok$VPD, "EI_direct")
m_Rg <- gam(Rg ~ s(VPD, k = 20), data = df_pred_ok, method = "REML")
m_Ta <- gam(Ta ~ s(VPD, k = 20), data = df_pred_ok, method = "REML")

vpd_seq <- seq(max(min(df_pred_ok$VPD), .Machine$double.eps),
               max(df_pred_ok$VPD), length.out = 1000)

pred_df <- tibble(
  VPD    = vpd_seq,
  AI_hat = m_AI$pred_fun(vpd_seq),
  P_hat  = m_P$pred_fun(vpd_seq),
  EI_hat_direct = m_EI_dir$pred_fun(vpd_seq),
  Rg_hat = as.numeric(predict(m_Rg, tibble(VPD = vpd_seq))),
  Ta_hat = as.numeric(predict(m_Ta, tibble(VPD = vpd_seq)))
)

jarvis_in <- pred_df |> transmute(VPD, P = P_hat, Rg = Rg_hat, Ta = Ta_hat)
jarvis_out <- predict_ET_fixed(par_fixed, jarvis_in)

pred_df <- pred_df |>
  mutate(
    g_eff_hat   = jarvis_out$g_eff,
    ET_hat      = jarvis_out$ET,
    EI_hat_J    = ET_hat / P_hat,
    AI_hat_cons = NA_real_
  )

# ============================
# Plots from predicted drivers
# ============================
p_top <- list(
  ggplot(pred_df, aes(VPD, AI_hat)) + geom_line() +
    labs(title="AI(VPD) prediction", x="VPD (kPa)", y="AI") + theme_bw(),
  ggplot(pred_df, aes(VPD, P_hat)) + geom_line() +
    labs(title="P(VPD) prediction", x="VPD (kPa)", y="P") + theme_bw(),
  ggplot(pred_df, aes(VPD, Rg_hat)) + geom_line() +
    labs(title="Rg(VPD) prediction", x="VPD (kPa)", y="Rg") + theme_bw(),
  ggplot(pred_df, aes(VPD, Ta_hat)) + geom_line() +
    labs(title="Ta(VPD) prediction", x="VPD (kPa)", y="Ta") + theme_bw()
)
print(p_top[[1]]); print(p_top[[2]]); print(p_top[[3]]); print(p_top[[4]])

p_ei <- ggplot(pred_df, aes(VPD)) +
  geom_line(aes(y = EI_hat_direct), linetype = 2) +
  geom_line(aes(y = EI_hat_J)) +
  labs(title = "EI(VPD): direct fit (dashed) vs Jarvis-consistent (solid)",
       x = "VPD (kPa)", y = "EI") +
  theme_bw()
print(p_ei)

p_ai_ei_pred <- ggplot(pred_df, aes(x = AI_hat, y = EI_hat_J, color = VPD)) +
  geom_path(linewidth = 1.2) +
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(title = "AI vs EI (predicted from VPD; EI via Jarvis & P(VPD))",
       x = "AI(VPD) prediction", y = "EI(VPD) via Jarvis") +
  theme_bw() + theme(plot.title = element_text(face="bold"))
print(p_ai_ei_pred)

# =============================
# Plots (2×2): AI, PET, ET, EI
# =============================
op <- par(mfrow = c(2,2))

# (1) AI: observed vs AI_hat = PET_hat / P_obs (reference)
plot(df0$VPD, df0$AI_obs, pch=16, cex=.7,
     main="AI_obs vs VPD (and AI_hat ref.)",
     xlab="VPD (kPa)", ylab="AI")
points(df0$VPD, df0$AI_hat, col=2, pch=1)
legend("topleft", legend=c("AI_obs = PET/P", "AI_hat = PET_hat/P"),
       pch=c(16,1), col=c(1,2), bty="n")

# (2) PET fit
plot(df0$VPD, df0$PET, pch=16, cex=.7,
     main="PET: obs vs a,b,C fit",
     xlab="VPD (kPa)", ylab="PET (ETo)")
points(df0$VPD, df0$PET_hat, col=4, pch=1)
legend("topleft", legend=c("Observed","Fitted"),
       pch=c(16,1), col=c(1,4), bty="n")

# (3) ET vs VPD (Jarvis with observed P, Rg, Ta)
plot(df_out$VPD, df_out$ET_obs, pch=16, cex=.7,
     main="ET vs VPD (Jarvis, obs P/Rg/Ta)",
     xlab="VPD (kPa)", ylab="ET")
points(df_out$VPD, df_out$ET_pred, col=3, pch=1)
legend("topleft", legend=c("Observed","Jarvis Pred"),
       pch=c(16,1), col=c(1,3), bty="n")

# (4) EI vs VPD (EI = ET / P_obs)
plot(df_out$VPD, df_out$EI_obs, pch=16, cex=.7,
     main="EI = ET/P_obs vs VPD (obs Rg/Ta)",
     xlab="VPD (kPa)", ylab="Evaporative Index (EI)")
points(df_out$VPD, df_out$EI_pred, col=4, pch=1)
legend("topright", legend=c("Observed","Jarvis Pred"),
       pch=c(16,1), col=c(1,4), bty="n")

par(op)


plot(df0$AI_obs, df0$AI_hat)
plot(df0$PET, df0$PET_hat)
plot(df_out$ET_obs, df_out$ET_pred)
plot(df_out$EI_obs, df_out$EI_pred)


# ==========================================
# AI–EI scatter/curve for OBSERVED pipeline
# ==========================================
df_plot_obs <- df_out |>
  select(.row_id, VPD, EI_pred) |>
  left_join(df0 |> select(.row_id, AI_hat), by = ".row_id") |>
  filter(is.finite(AI_hat), is.finite(EI_pred), is.finite(VPD))

curve_df_obs <- df_out |>
  select(.row_id, VPD, P, Rg, Ta) |>
  left_join(df0 |> select(.row_id, PET_hat), by = ".row_id") |>
  mutate(
    AI_hat = PET_hat / P,
    inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
    K_ET = K_ET_fun(VPD, Ta, rhoAir, CpAir, gamma)
  )

pred_curve_obs <- predict_ET_fixed(
  par_fixed,
  curve_df_obs |> transmute(VPD, P, Rg, Ta, CO2_term = 0, inv_sqrt_VPD, K_ET)
)

curve_df_obs <- curve_df_obs |>
  mutate(EI_pred = pred_curve_obs$ET / P) |>
  arrange(VPD) |>
  select(AI_hat, EI_pred, VPD)

p_ai_ei_obs <- ggplot(df_plot_obs, aes(x = AI_hat, y = EI_pred, color = VPD)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_path(data = curve_df_obs, aes(x = AI_hat, y = EI_pred, color = VPD),
            linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(
    title = "AI (= PET_hat / P_obs) vs EI (= ET_pred / P_obs)",
    x = "AI (pred via PET_hat / P_obs)",
    y = "EI (ET_pred / P_obs)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
print(p_ai_ei_obs)

# =======================
# Efficiency diagnostics
# =======================
stopifnot(all(c("VPD","AI_hat","EI_hat_J") %in% names(pred_df)))
geo <- pred_df |>
  arrange(VPD) |>
  transmute(VPD, AI_hat, EI_hat = EI_hat_J) |>
  mutate(
    dV  = lead(VPD)   - VPD,
    dAI = lead(AI_hat) - AI_hat,
    dEI = lead(EI_hat) - EI_hat
  ) |>
  head(-1) |>
  mutate(
    dAI_dV = dAI / pmax(dV, .Machine$double.eps),
    dEI_dV = dEI / pmax(dV, .Machine$double.eps),
    ds_dV  = sqrt(dAI_dV^2 + dEI_dV^2),
    eff_h  = pmin(1, pmax(0, abs(dAI_dV) / pmax(ds_dV, .Machine$double.eps)))
  )

i_peak   <- which.max(pred_df$EI_hat_J)
vpd_peak <- pred_df$VPD[i_peak]
AI_peak  <- pred_df$AI_hat[i_peak]

p_eff_vpd <- ggplot(geo, aes(VPD, eff_h)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = vpd_peak, linetype = 2) +
  scale_y_continuous(limits = c(0.9,1)) +
  labs(title = "Horizontal efficiency vs VPD",
       subtitle = "eff = |dAI/dVPD| / sqrt((dAI/dVPD)^2 + (dEI/dVPD)^2)",
       x = "VPD (kPa)", y = "Horizontal efficiency [0–1]") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

p_eff_ai <- ggplot(geo, aes(AI_hat, eff_h)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = AI_peak, linetype = 2) +
  scale_y_continuous(limits = c(0.93,1)) +
  labs(title = "Horizontal efficiency vs AI",
       subtitle = "1 → purely horizontal motion (dEI ≈ 0)",
       x = "AI (predicted)", y = "Horizontal efficiency") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

print(p_eff_vpd); print(p_eff_ai)

# =======================================================
# Add layers into existing ggplot object p9 (if present)
# =======================================================
if (exists("p9", inherits = TRUE)) {
  p99 <- p9 + ggnewscale::new_scale_color() +
    # predicted path overlays
    geom_point(
      data = pred_df,
      aes(x = AI_hat, y = EI_hat_J, color = VPD),
      size = 2, alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_path(
      data = pred_df %>% arrange(VPD),
      aes(x = AI_hat, y = EI_hat_J, color = VPD),
      linewidth = 1.2,
      inherit.aes = FALSE
    ) +
    # observed path overlays
    geom_point(data = df_plot_obs, aes(AI_hat, EI_pred, color = VPD),
               size = 2, alpha = 0.85, inherit.aes = FALSE) +
    geom_point(data = curve_df_obs, aes(AI_hat, EI_pred, color = VPD),
               size = 1.2, inherit.aes = FALSE) +
    scale_color_viridis_c(name = "VPD (kPa)") +
    labs(
      title = "AI vs EI (observed & predicted)",
      x = "AI", y = "EI"
    ) +
    theme_bw() + theme(plot.title = element_text(face = "bold"))
  print(p99)
}

# ========================================
# Model bundle (predicted path) for reuse
# ========================================
model_bundle <- list(
  fits = list(
    AI_log = m_AI, P_log = m_P, EI_log_direct = m_EI_dir,
    Rg_gam = m_Rg, Ta_gam = m_Ta
  ),
  par_fixed = par_fixed,
  predict_all = function(v){
    Rg_hat <- as.numeric(predict(m_Rg, data.frame(VPD=v)))
    Ta_hat <- as.numeric(predict(m_Ta, data.frame(VPD=v)))
    P_hat  <- m_P$pred_fun(v)
    AI_hat <- m_AI$pred_fun(v)
    jarvis_out <- predict_ET_fixed(par_fixed, data.frame(VPD=v, P=P_hat, Rg=Rg_hat, Ta=Ta_hat))
    tibble(
      VPD = v,
      AI_hat = AI_hat,
      P_hat  = P_hat,
      Rg_hat = Rg_hat,
      Ta_hat = Ta_hat,
      ET_hat = jarvis_out$ET,
      EI_hat_J = jarvis_out$ET / P_hat,
      EI_hat_direct = m_EI_dir$pred_fun(v)
    )
  }
)
# Example:
# preds_at <- model_bundle$predict_all(seq(0.15, 0.9, by=0.01))
