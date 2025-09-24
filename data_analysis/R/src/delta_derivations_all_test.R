suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggnewscale)  # only needed for optional p9 layering at the end
})

# ---- Guards ----
stopifnot(exists("Data_to_plot_II"))
stopifnot(all(c("rhoAir","CpAir","gamma") %in% ls(.GlobalEnv)))

# ============================================================
# 0) Prep: row id for safe alignment
# ============================================================
Data_to_plot_II <- Data_to_plot_II %>%
  mutate(.row_id = row_number())

# ============================================================
# 1) AI from observed P: AI_obs = PET / P; fit PET ~ log(VPD)+VPD
# ============================================================
df0 <- Data_to_plot_II %>%
  transmute(
    .row_id,
    VPD = as.numeric(VPD),
    PET = as.numeric(ETo_FAO56_alfalfa),
    P   = as.numeric(P)
  ) %>%
  filter(!is.na(VPD), !is.na(PET), !is.na(P),
         VPD > 0, PET > 0, P > 0) %>%
  mutate(AI_obs = PET / P)

# PET fit (a, b, C)
fit_PET <- lm(log(PET) ~ log(VPD) + VPD, data = df0)
a_PET   <- coef(fit_PET)[["log(VPD)"]]
b_PET   <- coef(fit_PET)[["VPD"]]
C_PET   <- coef(fit_PET)[["(Intercept)"]]
df0 <- df0 %>% mutate(PET_hat = exp(a_PET*log(VPD) + b_PET*VPD + C_PET),
                      AI_hat  = PET_hat / P)  # ensure AI = PET_hat / P_obs

# (Optional) AI_obs ~ VPD (diagnostic only)
fit_AI  <- lm(log(AI_obs) ~ log(VPD) + VPD, data = df0)
a_AI    <- coef(fit_AI)[["log(VPD)"]]
b_AI    <- coef(fit_AI)[["VPD"]]
C_AI    <- coef(fit_AI)[["(Intercept)"]]

cat("\nFitted relations (info):\n")
cat(sprintf("log(PET)    = %.4f*log(VPD) + %.4f*VPD + %.4f\n", a_PET, b_PET, C_PET))
cat(sprintf("log(AI_obs) = %.4f*log(VPD) + %.4f*VPD + %.4f  [diagnostic]\n", a_AI, b_AI, C_AI))

# ============================================================
# 2) Jarvis ET with FIXED params, OBSERVED P, OBSERVED Rg & Ta
# ============================================================
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

# Inputs with observed Rg/Ta/P
df <- Data_to_plot_II %>%
  transmute(
    .row_id,
    ET_obs = as.numeric(ET),
    VPD    = as.numeric(VPD),
    P      = as.numeric(P),      # observed P
    Rg     = as.numeric(Rg),     # observed Rg
    Ta     = as.numeric(Ta)      # observed Ta
  ) %>%
  mutate(
    CO2_term     = 0,
    inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
    K_ET = (rhoAir * CpAir / gamma) * VPD * (1/1000) *
      1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
  ) %>%
  filter(is.finite(ET_obs), is.finite(VPD), is.finite(P),
         is.finite(Rg), is.finite(Ta))

stopifnot(nrow(df) > 0)

# Jarvis functions
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
  
  fCO2 <- if (use_CO2) exp(-pmax(k_CO2, 0) * input_df$CO2_term) else 1
  
  if (!use_VPD) {
    fVPD <- 1
  } else {
    raw  <- b0_VPD + b1_VPD * input_df$inv_sqrt_VPD
    fVPD <- pmax(raw, 0)
  }
  
  fP * fRg * fTa * fVPD * fCO2
}

predict_ET_fixed <- function(par, input_df) {
  g_eff <- jarvis_g_eff_fixed(par, input_df)
  ET    <- input_df$K_ET * g_eff
  list(g_eff = g_eff, ET = ET)
}

# Predict ET and evaluate
preds <- predict_ET_fixed(par_fixed, df)

df_out <- df %>%
  mutate(
    g_eff_pred = preds$g_eff,
    ET_pred    = preds$ET,
    resid      = ET_obs - ET_pred
  )

RMSE <- sqrt(mean(df_out$resid^2, na.rm = TRUE))
R2   <- suppressWarnings(cor(df_out$ET_obs, df_out$ET_pred, use = "complete.obs")^2)
cat(sprintf("\nJarvis ET with OBSERVED P, Rg, Ta: RMSE = %.3f, R^2 = %.3f\n", RMSE, R2))

# EI = ET / P   (P already present in df_out)
df_out <- df_out %>%
  mutate(
    EI_obs  = ET_obs / P,
    EI_pred = ET_pred / P
  )

# ============================================================
# 3) Plots (2×2): AI, PET, ET, EI
# ============================================================
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

# ============================================================
# 4) 2D AI–EI plot colored by VPD (AI_hat = PET_hat/P_obs; EI_pred = ET_pred/P_obs)
# ============================================================
# Align df_out to df0 via .row_id
df_plot <- df_out %>%
  select(.row_id, VPD, EI_pred) %>%
  left_join(df0 %>% select(.row_id, AI_hat), by = ".row_id") %>%
  filter(is.finite(AI_hat), is.finite(EI_pred), is.finite(VPD))

# Parametric line on the SAME rows (observed P, Rg, Ta for ET; PET_hat for AI)
curve_df <- df_out %>%
  select(.row_id, VPD, P, Rg, Ta) %>%
  left_join(df0 %>% select(.row_id, PET_hat), by = ".row_id") %>%
  mutate(
    AI_hat = PET_hat / P,
    inv_sqrt_VPD = 1 / sqrt(pmax(VPD, 1e-6)),
    K_ET = (rhoAir * CpAir / gamma) * VPD * (1/1000) *
      1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
  )

pred_curve <- predict_ET_fixed(
  par_fixed,
  curve_df %>%
    transmute(VPD, P, Rg, Ta, CO2_term = 0, inv_sqrt_VPD, K_ET)
)

curve_df <- curve_df %>%
  mutate(EI_pred = pred_curve$ET / P) %>%
  arrange(VPD) %>%
  select(AI_hat, EI_pred, VPD)

# AI–EI scatter + parametric path colored by VPD
p_ai_ei <- ggplot(df_plot, aes(x = AI_hat, y = EI_pred, color = VPD)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_path(data = curve_df, aes(x = AI_hat, y = EI_pred, color = VPD),
            linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(
    title = "AI (= PET_hat / P_obs) vs EI (= ET_pred / P_obs)",
    x = "AI (pred via PET_hat / P_obs)",
    y = "EI (ET_pred / P_obs)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

print(p_ai_ei)

# ============================================================
# 5) Optional: add the AI–EI layers into an existing ggplot object p9
#    that already uses a discrete color scale
# ============================================================
p99 <- p9 +
  ggnewscale::new_scale_color() +
  geom_point(data = df_plot, aes(AI_hat, EI_pred, color = VPD),
             size = 2, alpha = 0.85, inherit.aes = FALSE) +
  geom_point(data = curve_df, aes(AI_hat, EI_pred, color = VPD),
            linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_viridis_c(name = "VPD (kPa)") +
  labs(title = "AI vs EI (obs P, obs Rg/Ta)",
       x = "AI (PET_hat/P_obs)", y = "EI (ET_pred/P_obs)") +
  theme_bw() + theme(plot.title = element_text(face = "bold"))
print(p99)
