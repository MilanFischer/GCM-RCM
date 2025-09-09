# --- Add after you have par_hat, df_opt, model_ET_components(), dln_CO2, eps defined ---

library(tidyverse)

# 1) Choose whether to plot both CO2 scenarios
show_both_CO2_scenarios <- TRUE

# 2) Baseline "hold-at" values (medians in the fit sample)
baseline <- df_opt |>
  summarise(
    AI_FAO56_alfalfa = median(AI_FAO56_alfalfa, na.rm = TRUE),
    Rg  = median(Rg,  na.rm = TRUE),
    Ta  = median(Ta,  na.rm = TRUE),
    VPD = median(VPD, na.rm = TRUE)
  )

# 3) Plot ranges (5th–95th percentiles to avoid outliers)
rng <- df_opt |>
  summarise(
    AI_lo  = quantile(AI_FAO56_alfalfa, 0.01, na.rm = TRUE),
    AI_hi  = quantile(AI_FAO56_alfalfa, 0.99, na.rm = TRUE),
    Rg_lo  = quantile(Rg,  0.01, na.rm = TRUE),
    Rg_hi  = quantile(Rg,  0.99, na.rm = TRUE),
    Ta_lo  = quantile(Ta,  0.01, na.rm = TRUE),
    Ta_hi  = quantile(Ta,  0.99, na.rm = TRUE),
    VPD_lo = quantile(VPD, 0.01, na.rm = TRUE),
    VPD_hi = quantile(VPD, 0.99, na.rm = TRUE)
  )

seq_AI  <- seq(rng$AI_lo,  rng$AI_hi,  length.out = 1000)
seq_Rg  <- seq(rng$Rg_lo,  rng$Rg_hi,  length.out = 1000)
seq_Ta  <- seq(rng$Ta_lo,  rng$Ta_hi,  length.out = 1000)
seq_VPD <- seq(0.1, rng$VPD_hi, length.out = 1000)

# 4) Helper to make one curve for a given variable and CO2 scenario
make_curve <- function(var_name, seq_vals, CO2_term_value) {
  ref <- baseline[rep(1, length(seq_vals)), ] |> as.data.frame()
  ref[[var_name]] <- seq_vals
  
  newdata <- tibble(
    AI_FAO56_alfalfa = ref$AI_FAO56_alfalfa,
    Rg  = ref$Rg,
    Ta  = ref$Ta,
    VPD = ref$VPD,
    CO2_term = CO2_term_value
  )
  comps <- model_ET_components(par_hat, newdata)
  tibble(
    var = var_name,
    x   = seq_vals,
    ET_predicted     = comps$ET_pred,
    g_eff_predicted  = comps$g_eff_pred
  )
}

# 5) Build curves for chosen CO2 scenario(s)
curves_baseline <- bind_rows(
  make_curve("AI_FAO56_alfalfa", seq_AI,  0),
  make_curve("Rg",               seq_Rg,  0),
  make_curve("Ta",               seq_Ta,  0),
  make_curve("VPD",              seq_VPD, 0)
) |> mutate(scenario = "No CO2 effect")

curves_flagged <- if (isTRUE(show_both_CO2_scenarios)) {
  bind_rows(
    make_curve("AI_FAO56_alfalfa", seq_AI,  dln_CO2),
    make_curve("Rg",               seq_Rg,  dln_CO2),
    make_curve("Ta",               seq_Ta,  dln_CO2),
    make_curve("VPD",              seq_VPD, dln_CO2)
  ) |> mutate(scenario = "CMIP future (CO2 active)")
} else {
  tibble()
}

curves <- bind_rows(curves_baseline, curves_flagged) |>
  mutate(
    # drop non-finite / non-physical
    ET_predicted    = ifelse(!is.finite(ET_predicted) | !is.finite(g_eff_predicted) | g_eff_predicted <= 0, NA, ET_predicted),
    g_eff_predicted = ifelse(!is.finite(g_eff_predicted), NA, g_eff_predicted)
  )

# Pretty facet labels
var_lab <- c(
  AI_FAO56_alfalfa = "Aridity Index (AI)",
  Rg               = "Global Radiation (Rg)",
  Ta               = "Air Temperature (°C)",
  VPD              = "VPD (kPa)"
)
curves$var <- factor(curves$var, levels = names(var_lab), labels = unname(var_lab))

# Data for vertical median lines in each facet
vline_df <- tibble(
  var = factor(names(var_lab), levels = names(var_lab), labels = unname(var_lab)),
  x   = c(baseline$AI_FAO56_alfalfa, baseline$Rg, baseline$Ta, baseline$VPD)
)

# 6) Plot: ET response
p_ET <- ggplot(curves, aes(x = x, y = ET_predicted, linetype = scenario)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_vline(data = vline_df, aes(xintercept = x), linewidth = 0.3, linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ var, scales = "free", ncol = 2) +   # <- was "free_x"
  labs(
    x = NULL,
    y = "ET (mm/yr)",
    title = "Sensitivity of ET to AI, Rg, Ta, and VPD"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = if (isTRUE(show_both_CO2_scenarios)) "bottom" else "none")

# 7) Plot: g_eff response
p_geff <- ggplot(curves, aes(x = x, y = g_eff_predicted, linetype = scenario)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_vline(data = vline_df, aes(xintercept = x), linewidth = 0.3, linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ var, scales = "free", ncol = 2) +  
  labs(
    x = NULL,
    y = expression(g[eff]~"(mm/s)"),
    title = expression("Sensitivity of " * g[eff] * " to AI, Rg, Ta, and VPD")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = if (isTRUE(show_both_CO2_scenarios)) "bottom" else "none")

print(p_ET)
print(p_geff)

# ======================
# BOOTSTRAP (with jittered starts) — place this AFTER par_hat/df_opt exist
# ======================

bootstrap_params <- function(B = 200, data = df_opt, par_init = par_hat, jitter_sd = 0.05) {
  # one bootstrap draw
  boot_one <- function() {
    idx <- sample.int(nrow(data), replace = TRUE)
    dat_b <- data[idx, , drop = FALSE]
    
    # Closure over dat_b so we don't need extra args in optim()
    obj_b <- function(par) obj_fun_ET(par, dat_b)
    
    # ---- jittered starting point (5% by default) ----
    par0 <- par_init * (1 + rnorm(length(par_init), sd = jitter_sd))
    
    # Try BFGS from jittered start; fallback to Nelder–Mead
    fit <- try(
      optim(par = par0, fn = obj_b, method = "BFGS",
            control = list(reltol = 1e-10, maxit = 800)),
      silent = TRUE
    )
    if (inherits(fit, "try-error") || !is.finite(fit$value)) {
      fit <- try(
        optim(par = par0, fn = obj_b, method = "Nelder-Mead",
              control = list(maxit = 1500)),
        silent = TRUE
      )
    }
    if (inherits(fit, "try-error") || !is.finite(fit$value)) {
      return(rep(NA_real_, length(par_init)))
    }
    as.numeric(fit$par)
  }
  
  mat <- replicate(B, boot_one())
  mat <- t(mat)
  colnames(mat) <- names(par_init)
  tibble::as_tibble(mat)
}

# ---- Run the bootstrap and summarise ----
set.seed(123)
boot_draws <- bootstrap_params(B = 200, data = df_opt, par_init = par_hat, jitter_sd = 0.05)

# Success rate per parameter (just for sanity)
succ_rate <- 1 - colMeans(is.na(boot_draws))
message(sprintf("Bootstrap success rate per parameter: %s",
                paste(sprintf("%s=%.0f%%", names(succ_rate), 100*succ_rate), collapse = ", ")))

boot_summary <- boot_draws |>
  tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "estimate") |>
  dplyr::group_by(parameter) |>
  dplyr::summarise(
    mean  = if (all(is.na(estimate))) NA_real_ else mean(estimate, na.rm = TRUE),
    sd    = if (all(is.na(estimate))) NA_real_ else sd(estimate, na.rm = TRUE),
    p2.5  = if (all(is.na(estimate))) NA_real_ else quantile(estimate, 0.025, na.rm = TRUE),
    p97.5 = if (all(is.na(estimate))) NA_real_ else quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

param_table_boot <- tibble::tibble(parameter = names(par_hat), point_est = as.numeric(par_hat)) |>
  dplyr::left_join(boot_summary, by = "parameter")

print(param_table_boot)

# ------------------------------------------------------------
# BOOTSTRAP-DRIVEN SENSITIVITY PLOTS WITH CONFIDENCE BELTS
# (place this AFTER you have par_hat/df_opt/model_ET_components/dln_CO2/eps)
# ------------------------------------------------------------
library(tidyverse)
library(numDeriv)

# --- 0) Bootstrap parameter draws (reuse if already available) ---
if (!exists("boot_draws")) {
  bootstrap_params <- function(B = 200, data = df_opt, par_init = par_hat, jitter_sd = 0.05) {
    boot_one <- function() {
      idx <- sample.int(nrow(data), replace = TRUE)
      dat_b <- data[idx, , drop = FALSE]
      obj_b <- function(par) obj_fun_ET(par, dat_b)
      par0  <- par_init * (1 + rnorm(length(par_init), sd = jitter_sd))   # jitter
      fit <- try(optim(par = par0, fn = obj_b, method = "BFGS",
                       control = list(reltol = 1e-10, maxit = 800)), silent = TRUE)
      if (inherits(fit, "try-error") || !is.finite(fit$value)) {
        fit <- try(optim(par = par0, fn = obj_b, method = "Nelder-Mead",
                         control = list(maxit = 1500)), silent = TRUE)
      }
      if (inherits(fit, "try-error") || !is.finite(fit$value)) return(rep(NA_real_, length(par_init)))
      as.numeric(fit$par)
    }
    mat <- replicate(B, boot_one())
    mat <- t(mat)
    colnames(mat) <- names(par_init)
    as_tibble(mat)
  }
  
  set.seed(123)
  boot_draws <- bootstrap_params(B = 200, data = df_opt, par_init = par_hat, jitter_sd = 0.05)
}

# Clean draws: keep rows with all finite parameters
boot_draws_clean <- boot_draws |> filter(if_all(everything(), is.finite))

# --- 1) CO2 scenario decision (auto-hide if no effect) ---
tol_pct <- 0.5  # treat |effect| < 0.5% in g_eff as "no effect"
k_used <- if (exists("use_CO2") && isTRUE(use_CO2)) max(unname(par_hat["par_CO2"]), 0) else 0
has_flagged <- any(Data_to_plot_II$CO2_term > 0, na.rm = TRUE)
decline_factor <- exp(-k_used * dln_CO2)
pct_change     <- (decline_factor - 1) * 100
show_CO2_line  <- has_flagged && (k_used > 0) && (abs(pct_change) >= tol_pct)
legend_pos     <- if (show_CO2_line) "bottom" else "none"

# --- 2) Baselines & ranges (hold others at medians, plot 1–99% range) ---
baseline <- df_opt |>
  summarise(
    AI_FAO56_alfalfa = median(AI_FAO56_alfalfa, na.rm = TRUE),
    Rg  = median(Rg,  na.rm = TRUE),
    Ta  = median(Ta,  na.rm = TRUE),
    VPD = median(VPD, na.rm = TRUE)
  )

rng <- df_opt |>
  summarise(
    AI_lo  = quantile(AI_FAO56_alfalfa, 0.01, na.rm = TRUE),
    AI_hi  = quantile(AI_FAO56_alfalfa, 0.99, na.rm = TRUE),
    Rg_lo  = quantile(Rg,  0.01, na.rm = TRUE),
    Rg_hi  = quantile(Rg,  0.99, na.rm = TRUE),
    Ta_lo  = quantile(Ta,  0.01, na.rm = TRUE),
    Ta_hi  = quantile(Ta,  0.99, na.rm = TRUE),
    VPD_lo = quantile(VPD, 0.01, na.rm = TRUE),
    VPD_hi = quantile(VPD, 0.99, na.rm = TRUE)
  )

seq_AI  <- seq(rng$AI_lo,  rng$AI_hi,  length.out = 1000)
seq_Rg  <- seq(rng$Rg_lo,  rng$Rg_hi,  length.out = 1000)
seq_Ta  <- seq(rng$Ta_lo,  rng$Ta_hi,  length.out = 1000)
seq_VPD <- seq(0.1,        rng$VPD_hi, length.out = 1000)  # guard very small VPD

# --- 3) Helpers to build curves for any parameter vector ---
make_curve_par <- function(par_vec, var_name, seq_vals, CO2_term_value) {
  ref <- baseline[rep(1, length(seq_vals)), ] |> as.data.frame()
  ref[[var_name]] <- seq_vals
  newdata <- tibble(
    AI_FAO56_alfalfa = ref$AI_FAO56_alfalfa,
    Rg  = ref$Rg,
    Ta  = ref$Ta,
    VPD = ref$VPD,
    CO2_term = CO2_term_value
  )
  comps <- model_ET_components(par_vec, newdata)
  tibble(
    var = var_name,
    x   = seq_vals,
    ET_predicted     = comps$ET_pred,
    g_eff_predicted  = comps$g_eff_pred
  )
}

# --- 4) Point-estimate curves (par_hat) for lines ---
curves_point <- bind_rows(
  make_curve_par(par_hat, "AI_FAO56_alfalfa", seq_AI,  0),
  make_curve_par(par_hat, "Rg",               seq_Rg,  0),
  make_curve_par(par_hat, "Ta",               seq_Ta,  0),
  make_curve_par(par_hat, "VPD",              seq_VPD, 0)
) |> mutate(scenario = "No CO2 effect")

if (show_CO2_line) {
  curves_point <- bind_rows(
    curves_point,
    bind_rows(
      make_curve_par(par_hat, "AI_FAO56_alfalfa", seq_AI,  dln_CO2),
      make_curve_par(par_hat, "Rg",               seq_Rg,  dln_CO2),
      make_curve_par(par_hat, "Ta",               seq_Ta,  dln_CO2),
      make_curve_par(par_hat, "VPD",              seq_VPD, dln_CO2)
    ) |> mutate(scenario = "CMIP future (CO2 active)")
  )
}

# --- 5) Bootstrap ribbons: per var & scenario, 5–95% bands across draws ---
compute_belt <- function(var_name, seq_vals, CO2_term_value, scen_label) {
  if (nrow(boot_draws_clean) == 0) return(
    list(ET = tibble(), g = tibble())
  )
  # build matrices of predictions across draws
  ET_mat <- matrix(NA_real_, nrow = length(seq_vals), ncol = nrow(boot_draws_clean))
  g_mat  <- matrix(NA_real_, nrow = length(seq_vals), ncol = nrow(boot_draws_clean))
  for (j in seq_len(nrow(boot_draws_clean))) {
    th <- as.numeric(boot_draws_clean[j, ])
    names(th) <- names(par_hat)
    crv <- make_curve_par(th, var_name, seq_vals, CO2_term_value)
    ET_mat[, j] <- crv$ET_predicted
    g_mat[, j]  <- crv$g_eff_predicted
  }
  tib_ET <- tibble(
    var = var_name, x = seq_vals,
    lo = apply(ET_mat, 1, quantile, probs = 0.05, na.rm = TRUE),
    hi = apply(ET_mat, 1, quantile, probs = 0.95, na.rm = TRUE),
    scenario = scen_label
  )
  tib_g <- tibble(
    var = var_name, x = seq_vals,
    lo = apply(g_mat, 1, quantile, probs = 0.05, na.rm = TRUE),
    hi = apply(g_mat, 1, quantile, probs = 0.95, na.rm = TRUE),
    scenario = scen_label
  )
  list(ET = tib_ET, g = tib_g)
}

# baseline ribbons
belts <- list()
belts$AI0  <- compute_belt("AI_FAO56_alfalfa", seq_AI,  0, "No CO2 effect")
belts$Rg0  <- compute_belt("Rg",               seq_Rg,  0, "No CO2 effect")
belts$Ta0  <- compute_belt("Ta",               seq_Ta,  0, "No CO2 effect")
belts$VPD0 <- compute_belt("VPD",              seq_VPD, 0, "No CO2 effect")

# flagged ribbons (only if CO2 active)
if (show_CO2_line) {
  belts$AI1  <- compute_belt("AI_FAO56_alfalfa", seq_AI,  dln_CO2, "CMIP future (CO2 active)")
  belts$Rg1  <- compute_belt("Rg",               seq_Rg,  dln_CO2, "CMIP future (CO2 active)")
  belts$Ta1  <- compute_belt("Ta",               seq_Ta,  dln_CO2, "CMIP future (CO2 active)")
  belts$VPD1 <- compute_belt("VPD",              seq_VPD, dln_CO2, "CMIP future (CO2 active)")
}

ribbons_ET <- bind_rows(purrr::map(belts, "ET"))
ribbons_g  <- bind_rows(purrr::map(belts, "g"))

# --- 6) Labels & facets ---
var_lab <- c(
  AI_FAO56_alfalfa = "Aridity Index (AI)",
  Rg               = "Global Radiation (Rg)",
  Ta               = "Air Temperature (°C)",
  VPD              = "VPD (kPa)"
)

curves_point$var <- factor(curves_point$var, levels = names(var_lab), labels = unname(var_lab))
ribbons_ET$var   <- factor(ribbons_ET$var,   levels = names(var_lab), labels = unname(var_lab))
ribbons_g$var    <- factor(ribbons_g$var,    levels = names(var_lab), labels = unname(var_lab))

# median vlines
vline_df <- tibble(
  var = factor(names(var_lab), levels = names(var_lab), labels = unname(var_lab)),
  x   = c(baseline$AI_FAO56_alfalfa, baseline$Rg, baseline$Ta, baseline$VPD)
)

# --- 7) Plot: ET with ribbons ---
p_ET <- ggplot() +
  geom_ribbon(data = ribbons_ET,
              aes(x = x, ymin = lo, ymax = hi),
              fill = "grey40", alpha = 0.20, show.legend = FALSE) +   # <- constant fill
  geom_line(data = curves_point,
            aes(x = x, y = ET_predicted, linetype = scenario),
            linewidth = 0.9) +
  geom_vline(data = vline_df, aes(xintercept = x),
             linewidth = 0.3, linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ var, scales = "free", ncol = 2) +
  labs(x = NULL, y = "ET (mm/yr)",
       title = "Sensitivity of ET to AI, Rg, Ta, and VPD") +
  theme_bw(base_size = 12) +
  theme(legend.position = legend_pos)

# --- 8) Plot: g_eff with ribbons ---
p_geff <- ggplot() +
  geom_ribbon(data = ribbons_g,
              aes(x = x, ymin = lo, ymax = hi),
              fill = "grey40", alpha = 0.20, show.legend = FALSE) +   # <- constant fill
  geom_line(data = curves_point,
            aes(x = x, y = g_eff_predicted, linetype = scenario),
            linewidth = 0.9) +
  geom_vline(data = vline_df, aes(xintercept = x),
             linewidth = 0.3, linetype = "dashed", alpha = 0.6) +
  facet_wrap(~ var, scales = "free", ncol = 2) +
  labs(x = NULL, y = expression(g[eff]~"(mm/s)"),
       title = expression("Sensitivity of " * g[eff] * " to AI, Rg, Ta, and VPD")) +
  theme_bw(base_size = 12) +
  theme(legend.position = legend_pos)

print(p_ET)
print(p_geff)

