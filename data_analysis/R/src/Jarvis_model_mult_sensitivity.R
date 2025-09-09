# ======================
# Response curves for each driver (g_eff and ET)
# ======================

# Choose a baseline (medians from the fit subset)
baseline <- list(
  P   = median(df_opt$P, na.rm = TRUE),
  Rg   = median(df_opt$Rg, na.rm = TRUE),
  Ta   = median(df_opt$Ta, na.rm = TRUE),
  VPD  = median(df_opt$VPD, na.rm = TRUE),
  CO2_term = 0  # set to 0 for historical; change to dln_CO2 to show future
)

# Helper to build 'data' + 'drv' for arbitrary P, Rg, Ta, VPD vectors
make_inputs <- function(P, Rg, Ta, VPD, CO2_term_val) {
  data <- tibble(
    P = P,
    Rg = Rg,
    Ta = Ta,
    VPD = VPD,
    CO2_term = CO2_term_val
  )
  drv  <- build_drivers(data)
  list(data = data, drv = drv)
}

# Compute response curve for a single variable while fixing others at baseline
# var_name ∈ {"P","Rg","Ta","VPD"}
response_curve <- function(var_name, seq_vals, par_hat, baseline, co2_term_label = "Historical (CO2_term=0)") {
  # Expand to vectors for all inputs
  P  <- rep(baseline$P,    length(seq_vals))
  Rg  <- rep(baseline$Rg,  length(seq_vals))
  Ta  <- rep(baseline$Ta,  length(seq_vals))
  VPD <- rep(baseline$VPD, length(seq_vals))
  
  if (var_name == "P")   P  <- seq_vals
  if (var_name == "Rg")  Rg  <- seq_vals
  if (var_name == "Ta")  Ta  <- seq_vals
  if (var_name == "VPD") VPD <- seq_vals
  
  inputs <- make_inputs(P, Rg, Ta, VPD, baseline$CO2_term)
  comps  <- model_ET_components(par_hat, inputs$data, inputs$drv)
  
  tibble(
    var      = var_name,
    x        = seq_vals,
    g_eff    = as.numeric(comps$g_eff_pred),   # mm/s
    ET       = as.numeric(comps$ET_pred),      # mm/yr
    scenario = co2_term_label
  )
}

# Build sequences for each variable over the 1st–99th percentile to avoid outliers
pseq <- function(x, n = 200) {
  qs <- quantile(x[is.finite(x)], probs = c(0.01, 0.99), na.rm = TRUE)
  seq(qs[1], qs[2], length.out = n)
}

seq_P  <-  pseq(df_opt$P)
seq_Rg  <- pseq(df_opt$Rg)
seq_Ta  <- pseq(df_opt$Ta)
seq_VPD <- pseq(df_opt$VPD)

# Historical CO2 (CO2_term = 0)
baseline$CO2_term <- 0
rc_P_h  <- response_curve("P",    seq_P,  par_hat, baseline, "Historical")
rc_Rg_h  <- response_curve("Rg",  seq_Rg,  par_hat, baseline, "Historical")
rc_Ta_h  <- response_curve("Ta",  seq_Ta,  par_hat, baseline, "Historical")
rc_VPD_h <- response_curve("VPD", seq_VPD, par_hat, baseline, "Historical")

# Optional: Future CO2 (uncomment to add)
# baseline$CO2_term <- as.numeric(dln_CO2)
# rc_P_f  <- response_curve("P",  seq_P,  par_hat, baseline, "Future")
# rc_Rg_f  <- response_curve("Rg",  seq_Rg,  par_hat, baseline, "Future")
# rc_Ta_f  <- response_curve("Ta",  seq_Ta,  par_hat, baseline, "Future")
# rc_VPD_f <- response_curve("VPD", seq_VPD, par_hat, baseline, "Future")

# Combine (only historical by default)
resp_all <- bind_rows(rc_P_h, rc_Rg_h, rc_Ta_h, rc_VPD_h)
# If you computed future too, use: resp_all <- bind_rows(rc_P_h, rc_Rg_h, rc_Ta_h, rc_VPD_h, rc_P_f, rc_Rg_f, rc_Ta_f, rc_VPD_f)

# Nice labels
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")
x_lab <- function(var) {
  switch(var,
         P  = "P",
         Rg  = "Rg",
         Ta  = "Ta (°C)",
         VPD = "VPD (kPa)"
  )
}

# -------- PLOT: g_eff responses --------
library(ggplot2)
plot_geff <- ggplot(resp_all, aes(x = x, y = g_eff, color = scenario)) +
  geom_line(linewidth = 1, color = "black") +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*")"),
       title = expression("Single-variable response of " * g[eff]),
       subtitle = "Each curve varies one driver across its range; others held at median") +
  theme_bw() +
  theme(legend.position = "top")

# print(plot_geff)

# -------- PLOT: ET responses --------
plot_et <- ggplot(resp_all, aes(x = x, y = ET, color = scenario)) +
  geom_line(linewidth = 1, color = "black") +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = "ET (mm/yr)",
       title = "Single-variable response of ET",
       subtitle = "Each curve varies one driver across its range; others held at median") +
  theme_bw() +
  theme(legend.position = "top")

# print(plot_et)

# Save g_eff plot
ggsave("../Plots/ggplot2/response_geff.png", plot_geff, width = 8, height = 6, dpi = 300)

# Save ET plot
ggsave("../Plots/ggplot2/response_et.png", plot_et, width = 8, height = 6, dpi = 300)

# -------- Optional: mark the medians used for the held variables on x-axes --------
# Example for VPD panel: add a vertical line at median VPD
# plot_geff + geom_vline(xintercept = median(df_opt$VPD, na.rm = TRUE), linetype = 2)
# plot_et   + geom_vline(xintercept = median(df_opt$VPD, na.rm = TRUE), linetype = 2)