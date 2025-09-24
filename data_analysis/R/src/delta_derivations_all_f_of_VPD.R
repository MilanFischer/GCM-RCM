suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(mgcv)
})

# ---- Guards (you must already have these constants defined) ----
stopifnot(exists("Data_to_plot_II"))
stopifnot(all(c("rhoAir","CpAir","gamma") %in% ls(.GlobalEnv)))

# ============================================================
# 0) Data prep
# ============================================================
df_raw <- Data_to_plot_II %>%
  transmute(
    VPD = as.numeric(VPD),
    AI  = as.numeric(AI_FAO56_alfalfa),
    PET = as.numeric(ETo_FAO56_alfalfa),
    P   = as.numeric(P),
    Rg  = as.numeric(Rg),
    Ta  = as.numeric(Ta),
    ET  = as.numeric(ET)
  ) %>%
  filter(
    is.finite(VPD), VPD > 0,
    is.finite(AI),  AI  > 0,
    is.finite(PET), PET > 0,
    is.finite(P),   P   > 0,
    is.finite(Rg),
    is.finite(Ta),
    is.finite(ET)
  ) %>%
  mutate(
    EI_obs = ET / P
  )

# ============================================================
# 1) Fits: everything as a function of VPD
#     - Log-normal style for positive variables (AI, P, EI)
#     - Smooth GAMs for Rg, Ta
# ============================================================
# Log-form helper
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

# AI(VPD) and P(VPD) with your preferred log-form
m_AI <- fit_log_abC(df_raw$AI, df_raw$VPD, "AI")
m_P  <- fit_log_abC(df_raw$P , df_raw$VPD, "P")

# EI(VPD) direct (diagnostic; may differ from ET/P consistency)
m_EI_dir <- fit_log_abC(df_raw$EI_obs, df_raw$VPD, "EI_direct")

# Rg(VPD) and Ta(VPD): smooth GAMs
m_Rg <- gam(Rg ~ s(VPD, k = 20), data = df_raw, method = "REML")
m_Ta <- gam(Ta ~ s(VPD, k = 20), data = df_raw, method = "REML")

# ============================================================
# 2) Build a VPD grid and predict all drivers
# ============================================================
vpd_seq <- seq(max(min(df_raw$VPD), .Machine$double.eps),
               max(df_raw$VPD), length.out = 1000)

# vpd_seq <- seq(0.18, 0.5, length.out = 1000)

pred_df <- tibble(
  VPD    = vpd_seq,
  AI_hat = m_AI$pred_fun(vpd_seq),
  P_hat  = m_P$pred_fun(vpd_seq),
  EI_hat_direct = m_EI_dir$pred_fun(vpd_seq),
  Rg_hat = as.numeric(predict(m_Rg, tibble(VPD = vpd_seq))),
  Ta_hat = as.numeric(predict(m_Ta, tibble(VPD = vpd_seq)))
)

# ============================================================
# 3) Jarvis ET using *predicted* P(VPD), Rg(VPD), Ta(VPD)
# ============================================================
# Fixed Jarvis parameters (use your best set here)
par_fixed <- c(
  Rg_min  = -227.49578,
  Rg_max  = 3919.44218,
  Ta_min  = -172.17573,
  Ta_max  = 81.30892,
  P_min   = -181.16213,
  P_max   = 43192.41972,
  k_CO2   = 0.00000,
  b0_VPD  = -3956.17778,
  b1_VPD  = 5148.86787
)

use_P <- TRUE; use_Rg <- TRUE; use_Ta <- TRUE; use_CO2 <- TRUE; use_VPD <- TRUE
clamp_Rg <- TRUE; clamp_Ta <- TRUE; clamp_P <- TRUE
eps <- 1e-6

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
  
  fCO2 <- if (use_CO2) 1 else 1  # k_CO2=0 → neutral, keep simple
  raw  <- if (!use_VPD) 1 else (b0_VPD + b1_VPD * (1 / sqrt(pmax(input_df$VPD, eps))))
  fVPD <- pmax(raw, 0)
  
  fP * fRg * fTa * fVPD * fCO2
}

predict_ET_fixed <- function(par, input_df) {
  g_eff <- jarvis_g_eff_fixed(par, input_df)
  K_ET  <- (rhoAir * CpAir / gamma) * input_df$VPD * (1/1000) *
    1 / ((2.501e6 - 2361 * input_df$Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
  ET    <- K_ET * g_eff
  list(g_eff = g_eff, ET = ET)
}

jarvis_in <- pred_df %>%
  transmute(
    VPD = VPD,
    P   = P_hat,
    Rg  = Rg_hat,
    Ta  = Ta_hat
  )

jarvis_out <- predict_ET_fixed(par_fixed, jarvis_in)
pred_df <- pred_df %>%
  mutate(
    g_eff_hat   = jarvis_out$g_eff,
    ET_hat      = jarvis_out$ET,
    EI_hat_J    = ET_hat / P_hat,           # EI predicted via Jarvis consistency
    AI_hat_cons = (NA_real_)                # (slot if you later want PET_hat/P_hat)
  )

# ============================================================
# 4) Quick diagnostics / plots
# ============================================================

# AI, P, Rg, Ta vs VPD (predictions)
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

# EI: direct vs Jarvis-consistent
p_ei <- ggplot(pred_df, aes(VPD)) +
  geom_line(aes(y = EI_hat_direct), linetype = 2) +
  geom_line(aes(y = EI_hat_J)) +
  labs(title = "EI(VPD): direct fit (dashed) vs Jarvis-consistent (solid)",
       x = "VPD (kPa)", y = "EI") +
  theme_bw()

# AI–EI parametric curve (color = VPD)
p_ai_ei <- ggplot(pred_df, aes(x = AI_hat, y = EI_hat_J, color = VPD)) +
  geom_path(linewidth = 1.2) +
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(title = "AI vs EI (predicted from VPD; EI via Jarvis & P(VPD))",
       x = "AI(VPD) prediction", y = "EI(VPD) via Jarvis") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

# Print the key plots
print(p_top[[1]]); print(p_top[[2]]); print(p_top[[3]]); print(p_top[[4]])
print(p_ei)
print(p_ai_ei)


# ============================================================
# 5) Export a compact model bundle (optional)
# ============================================================
model_bundle <- list(
  fits = list(
    AI_log = m_AI, P_log = m_P, EI_log_direct = m_EI_dir,
    Rg_gam = m_Rg, Ta_gam = m_Ta
  ),
  par_fixed = par_fixed,
  predict_all = function(v){
    # one-stop predictions for any VPD vector v
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

# Example: preds_at <- model_bundle$predict_all(seq(0.15, 0.9, by=0.01))

p99 <- p9 +
  ggnewscale::new_scale_color() +
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
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(
    title = "AI vs EI (all predicted from VPD)",
    x = "AI (predicted)",
    y = "EI (Jarvis-consistent, using predicted P, Rg, Ta)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

print(p99)

################################################################################
# pred_df must have: VPD, AI_hat, EI_hat_J
stopifnot(all(c("VPD","AI_hat","EI_hat_J") %in% names(pred_df)))

eps <- .Machine$double.eps

geo <- pred_df %>%
  arrange(VPD) %>%
  transmute(
    VPD,
    AI_hat,
    EI_hat = EI_hat_J
  ) %>%
  mutate(
    dV  = lead(VPD)   - VPD,
    dAI = lead(AI_hat) - AI_hat,
    dEI = lead(EI_hat) - EI_hat
  ) %>%
  head(-1) %>%  # drop last NA-diff row
  mutate(
    dAI_dV = dAI / pmax(dV, eps),
    dEI_dV = dEI / pmax(dV, eps),
    ds_dV  = sqrt(dAI_dV^2 + dEI_dV^2),
    
    # Horizontal efficiency (fraction of motion along AI)
    eff_h  = pmin(1, pmax(0, abs(dAI_dV) / pmax(ds_dV, eps))),
    # Equivalent geometric form (independent of dV):
    eff_h_geom = pmin(1, pmax(0, abs(dAI) / pmax(sqrt(dAI^2 + dEI^2), eps)))
  )

# Locate EI peak (where EI is maximal; often near dEI/dVPD ~ 0)
i_peak   <- which.max(pred_df$EI_hat_J)
vpd_peak <- pred_df$VPD[i_peak]
AI_peak  <- pred_df$AI_hat[i_peak]

# 1) Efficiency vs VPD
p_eff_vpd <- ggplot(geo, aes(VPD, eff_h)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = vpd_peak, linetype = 2) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "Horizontal efficiency vs VPD",
       subtitle = "eff = |dAI/dVPD| / sqrt((dAI/dVPD)^2 + (dEI/dVPD)^2)",
       x = "VPD (kPa)", y = "Horizontal efficiency [0–1]") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

# 2) Efficiency vs AI with ylim
p_eff_ai <- ggplot(geo, aes(AI_hat, eff_h)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = AI_peak, linetype = 2) +
  scale_y_continuous(limits = c(0.93,1)) +   # <-- set ylim here
  labs(title = "Horizontal efficiency vs AI",
       subtitle = "1 → purely horizontal motion (dEI ≈ 0)",
       x = "AI (predicted)", y = "Horizontal efficiency") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

print(p_eff_ai)

# 3) (Optional) show components to verify where EI flattens
p_components <- ggplot(geo, aes(VPD)) +
  geom_line(aes(y = dAI_dV), color = "steelblue", linewidth = 1.0) +
  geom_line(aes(y = dEI_dV), color = "darkorange", linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = vpd_peak, linetype = 2) +
  labs(title = "Derivatives vs VPD",
       x = "VPD (kPa)", y = "Derivative",
       subtitle = "steelblue: dAI/dVPD,  darkorange: dEI/dVPD") +
  theme_bw() + theme(plot.title = element_text(face="bold"))

print(p_eff_vpd)
print(p_eff_ai)
# print(p_components)  # uncomment if you want to inspect derivatives
