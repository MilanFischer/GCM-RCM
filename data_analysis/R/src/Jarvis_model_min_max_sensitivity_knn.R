# ============================================================
# Response curves using KNN-conditional percentiles (Jarvis-compatible)
# Save as: ./src/Jarvis_response_curves_knn.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# -----------------------
# Requirements from the fit step (must exist)
# -----------------------
req_objs <- c("par_hat", "df_opt", "fVPD_mode", "rhoAir", "CpAir", "gamma", "model_ET_components")
missing <- setdiff(req_objs, ls(envir = .GlobalEnv))
if (length(missing)) stop("Missing in environment: ", paste(missing, collapse = ", "))

# Small numeric guard (use existing if you already defined one)
if (!exists("eps", inherits = TRUE)) eps <- 1e-6

# -----------------------
# Settings
# -----------------------
# Which fVPD handling to mirror
fVPD_mode_resp <- fVPD_mode

# K in the KNN neighborhood used to compute conditional percentiles
knn_k <- NULL  # default: ~5% of n (>=30). Set an integer to override, e.g., knn_k <- 50

# Percentiles for the "held" variables (used as conditional targets)
held_probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
held_lab   <- c("p10","p25","p50","p75","p90")

# For historical curves, hold CO2_term at 0 (set to dln_CO2 for future scenarios)
co2_term_val <- 0

# -----------------------
# Helpers
# -----------------------

# Build the exact inputs Jarvis expects (no scalers/driver objects)
make_inputs <- function(P, Rg, Ta, VPD, CO2_term_val) {
  tibble(P = P, Rg = Rg, Ta = Ta, VPD = VPD,
         CO2_term = CO2_term_val) |>
    mutate(
      inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
      K_ET = (rhoAir * CpAir / gamma) * VPD * (1/1000) *
        1 / ((2.501e6 - 2361 * Ta) / 1e6) * (3600 * 24) / 1e6 * 365.25
    )
}

# Predict components using your current Jarvis model
predict_components <- function(par_hat, input_df, fVPD_mode_resp) {
  model_ET_components(
    par       = par_hat,
    input_df  = input_df,
    fVPD_mode = fVPD_mode_resp
  )
}

# Create x-sequences inside the *active* fitted min/max range (avoids saturated flats)
active_seq <- function(name, n = 200, pad = 1e-6) {
  lo <- as.numeric(par_hat[paste0(name, "_min")]) + pad
  hi <- as.numeric(par_hat[paste0(name, "_max")]) - pad
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    # fallback to robust data range
    q <- stats::quantile(df_opt[[name]], c(0.01, 0.99), na.rm = TRUE)
    lo <- q[1]; hi <- q[2]
  }
  seq(lo, hi, length.out = n)
}

seq_P   <- active_seq("P")
seq_Rg  <- active_seq("Rg")
seq_Ta  <- active_seq("Ta")
seq_VPD <- active_seq("VPD")

# Conditional percentiles via KNN neighborhoods in the real data
cond_holds_knn <- function(var_name, x_seq, prob, data = df_opt, k = NULL) {
  if (is.null(k)) k <- max(30L, round(0.05 * nrow(data)))   # ~5% of data, at least 30
  X <- data[[var_name]]
  n <- length(x_seq)
  
  P_vec   <- numeric(n)
  Rg_vec  <- numeric(n)
  Ta_vec  <- numeric(n)
  VPD_vec <- numeric(n)
  
  for (j in seq_len(n)) {
    d   <- abs(X - x_seq[j])
    idx <- order(d)[1:k]
    P_vec[j]   <- as.numeric(stats::quantile(data$P[idx],   probs = prob, na.rm = TRUE))
    Rg_vec[j]  <- as.numeric(stats::quantile(data$Rg[idx],  probs = prob, na.rm = TRUE))
    Ta_vec[j]  <- as.numeric(stats::quantile(data$Ta[idx],  probs = prob, na.rm = TRUE))
    VPD_vec[j] <- as.numeric(stats::quantile(data$VPD[idx], probs = prob, na.rm = TRUE))
  }
  list(P = P_vec, Rg = Rg_vec, Ta = Ta_vec, VPD = VPD_vec)
}

# Build curves for one focal variable X, holding the others at *conditional* percentiles given X≈x
make_curves_for_var <- function(var_name, x_seq, par_hat, co2_term_val = 0) {
  out_list <- vector("list", length(held_probs))
  
  for (i in seq_along(held_probs)) {
    ch <- cond_holds_knn(var_name, x_seq, prob = held_probs[i], data = df_opt, k = knn_k)
    
    # start from conditional vectors, then set the focal variable explicitly
    P_vec   <- ch$P
    Rg_vec  <- ch$Rg
    Ta_vec  <- ch$Ta
    VPD_vec <- ch$VPD
    
    if (var_name == "P")   P_vec   <- x_seq
    if (var_name == "Rg")  Rg_vec  <- x_seq
    if (var_name == "Ta")  Ta_vec  <- x_seq
    if (var_name == "VPD") VPD_vec <- x_seq
    
    input_i <- make_inputs(P_vec, Rg_vec, Ta_vec, VPD_vec, co2_term_val)
    comp    <- predict_components(par_hat, input_i, fVPD_mode_resp)
    
    out_list[[i]] <- tibble(
      var      = var_name,
      x        = x_seq,
      g_eff    = as.numeric(comp$g_eff_pred),
      ET       = as.numeric(comp$ET_pred),
      scenario = paste0("held_", held_lab[i])
    )
  }
  
  dplyr::bind_rows(out_list)
}

# -----------------------
# Build all response curves
# -----------------------
rc_P   <- make_curves_for_var("P",   seq_P,   par_hat, co2_term_val = co2_term_val)
rc_Rg  <- make_curves_for_var("Rg",  seq_Rg,  par_hat, co2_term_val = co2_term_val)
rc_Ta  <- make_curves_for_var("Ta",  seq_Ta,  par_hat, co2_term_val = co2_term_val)
rc_VPD <- make_curves_for_var("VPD", seq_VPD, par_hat, co2_term_val = co2_term_val)

resp_all <- dplyr::bind_rows(rc_P, rc_Rg, rc_Ta, rc_VPD)

# -----------------------
# Plotting
# -----------------------
var_labels <- c(P = "P (precipitation)", Rg = "Rg (radiation)", Ta = "Ta (°C)", VPD = "VPD (kPa)")

plot_geff <- ggplot(resp_all, aes(x = x, y = g_eff, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*")"),
       title = expression("Single-variable response of " * g[eff]),
       subtitle = "KNN-conditional holds: others set to local percentiles of the real data") +
  theme_bw() +
  theme(legend.position = "top")

plot_et <- ggplot(resp_all, aes(x = x, y = ET, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
  labs(x = NULL, y = "ET (mm/yr)",
       title = "Single-variable response of ET",
       subtitle = "KNN-conditional holds: others set to local percentiles of the real data") +
  theme_bw() +
  theme(legend.position = "top")

# -----------------------
# Save (adjust paths as you like)
# -----------------------
ggsave("../Plots/ggplot2/response_geff_knn.png", plot_geff, width = 9, height = 6.5, dpi = 300)
ggsave("../Plots/ggplot2/response_et_knn.png",   plot_et,  width = 9, height = 6.5, dpi = 300)

# Expose objects if you want to use them interactively
jarvis_resp_curves <- resp_all
