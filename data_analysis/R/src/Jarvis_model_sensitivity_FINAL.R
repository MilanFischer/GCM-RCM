# ======================
# Response curves with percentile-held baselines
# ======================

# ======================
# Response curves with percentile-held baselines
# ======================

# Must match fit settings
fVPD_mode_resp <- fVPD_mode

# Derive response coefs from RAW every time (no global mutation)
if (identical(fVPD_mode_resp, "given") && is.null(fVPD_given_vec)) {
  vpd_coefs_resp <- convert_vpd_coefs(b0_VPD_raw, b1_VPD_raw, scalers, standardize_invVPD)
} else {
  vpd_coefs_resp <- c(b0 = NA_real_, b1 = NA_real_)
}

# Optional builder if you used a per-row fVPD vector in fitting (rare)
fVPD_given_builder <- NULL

# Percentiles for held variables
held_probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
held_lab   <- c("p10","p25","p50","p75","p90")

# Compute percentile table from df_opt
perc_tbl <- tibble(
  P   = quantile(df_opt$P,   probs = held_probs, na.rm = TRUE),
  Rg  = quantile(df_opt$Rg,  probs = held_probs, na.rm = TRUE),
  Ta  = quantile(df_opt$Ta,  probs = held_probs, na.rm = TRUE),
  VPD = quantile(df_opt$VPD, probs = held_probs, na.rm = TRUE),
  prob  = held_probs,
  label = held_lab
)

# Helper to build 'data' + 'drv' using same scalers/standardization as fitting
make_inputs <- function(P, Rg, Ta, VPD, CO2_term_val) {
  data <- tibble(P=P, Rg=Rg, Ta=Ta, VPD=VPD, CO2_term=CO2_term_val)
  drv  <- build_drivers(
    df = data, scalers = scalers,
    standardize = standardize_drivers,
    standardize_invVPD = standardize_invVPD
  )
  list(data = data, drv = drv)
}

# Wrapper calling your model components
predict_components <- function(par_hat, data, drv,
                               fVPD_mode_resp,
                               vpd_coefs_resp,
                               fVPD_given_builder = NULL) {
  fVPD_vec <- NULL
  if (identical(fVPD_mode_resp, "given") && !is.null(fVPD_given_builder)) {
    fVPD_vec <- fVPD_given_builder(drv)
    if (length(fVPD_vec) != nrow(drv)) stop("fVPD_given_builder returned wrong length.")
  }
  model_ET_components(
    par  = par_hat,
    data = data,
    drv  = drv,
    fVPD_mode      = fVPD_mode_resp,
    b0_VPD         = unname(vpd_coefs_resp["b0"]),
    b1_VPD         = unname(vpd_coefs_resp["b1"]),
    fVPD_given_vec = fVPD_vec
  )
}

# Build sequences for the x-variable (1–99% to avoid outliers)
pseq <- function(x, n = 200) {
  qs <- quantile(x[is.finite(x)], probs = c(0.01, 0.99), na.rm = TRUE)
  seq(qs[1], qs[2], length.out = n)
}
seq_P   <- pseq(df_opt$P)
seq_Rg  <- pseq(df_opt$Rg)
seq_Ta  <- pseq(df_opt$Ta)
seq_VPD <- pseq(df_opt$VPD)

# Curves for one focal variable, holding the other three at the SAME percentile q
make_curves_for_var <- function(var_name, x_seq, par_hat, co2_term_val = 0) {
  out_list <- vector("list", length(held_probs))
  for (i in seq_along(held_probs)) {
    P_h   <- perc_tbl$P[i]
    Rg_h  <- perc_tbl$Rg[i]
    Ta_h  <- perc_tbl$Ta[i]
    VPD_h <- perc_tbl$VPD[i]
    
    P_vec   <- rep(P_h,   length(x_seq))
    Rg_vec  <- rep(Rg_h,  length(x_seq))
    Ta_vec  <- rep(Ta_h,  length(x_seq))
    VPD_vec <- rep(VPD_h, length(x_seq))
    
    if (var_name == "P")   P_vec   <- x_seq
    if (var_name == "Rg")  Rg_vec  <- x_seq
    if (var_name == "Ta")  Ta_vec  <- x_seq
    if (var_name == "VPD") VPD_vec <- x_seq
    
    inp  <- make_inputs(P_vec, Rg_vec, Ta_vec, VPD_vec, co2_term_val)
    comp <- predict_components(par_hat, inp$data, inp$drv,
                               fVPD_mode_resp,
                               vpd_coefs_resp,
                               fVPD_given_builder)
    
    out_list[[i]] <- tibble(
      var      = var_name,
      x        = x_seq,
      g_eff    = as.numeric(comp$g_eff_pred),
      ET       = as.numeric(comp$ET_pred),
      scenario = paste0("held_", perc_tbl$label[i])  # e.g., held_p10, held_p50, ...
    )
  }
  bind_rows(out_list)
}

# Build all curves (Historical CO2_term=0; flip to dln_CO2 for future if desired)
rc_P   <- make_curves_for_var("P",   seq_P,   par_hat, co2_term_val = 0)
rc_Rg  <- make_curves_for_var("Rg",  seq_Rg,  par_hat, co2_term_val = 0)
rc_Ta  <- make_curves_for_var("Ta",  seq_Ta,  par_hat, co2_term_val = 0)
rc_VPD <- make_curves_for_var("VPD", seq_VPD, par_hat, co2_term_val = 0)

resp_all <- bind_rows(rc_P, rc_Rg, rc_Ta, rc_VPD)

# Labels
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")

# -------- PLOT: g_eff responses --------
plot_geff <- ggplot(resp_all, aes(x = x, y = g_eff, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*")"),
       title = expression("Single-variable response of " * g[eff]),
       subtitle = "Each curve varies one driver; others held at the SAME percentile (10/25/50/75/90)") +
  theme_bw() +
  theme(legend.position = "top")

# -------- PLOT: ET responses --------
plot_et <- ggplot(resp_all, aes(x = x, y = ET, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = "ET (mm/yr)",
       title = "Single-variable response of ET",
       subtitle = "Each curve varies one driver; others held at the SAME percentile (10/25/50/75/90)") +
  theme_bw() +
  theme(legend.position = "top")

# Save
ggsave("../Plots/ggplot2/response_geff.png", plot_geff, width = 9, height = 6.5, dpi = 300)
ggsave("../Plots/ggplot2/response_et.png",   plot_et,  width = 9, height = 6.5, dpi = 300)