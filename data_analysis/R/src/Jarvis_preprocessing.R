# ==============================================================================
# Helpers: data prep + physics helper + RF tuning helpers
# ==============================================================================

# ---- numeric stability constant ----
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# -----------------------
# Physics helper: K_ET per row (single source of truth)
# -----------------------
recompute_K_ET <- function(d) {
  stopifnot(all(c("VPD","Ta") %in% names(d)))
  vVPD <- pmax(d$VPD, eps)
  vTa  <- d$Ta
  (rhoAir * CpAir / gamma) * vVPD * (1/1000) *
    1 / ((2.501e6 - 2361 * vTa) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

# -----------------------
# Data prep (create needed columns FIRST)
# -----------------------
df <- Data_to_plot$abs |>
  dplyr::select(
    g_eff, ET, Rg, A, Ta, P, VPD, PERIOD, ensemble, model,
    color, fill, border, shape, linetype
  ) |>
  dplyr::mutate(
    # Clip VPD once, early
    VPD = pmax(VPD, eps),
    
    # CO2_term (safe, no use of '.')
    CO2_term = {
      if (exists("CO2_1981_2005", inherits = TRUE) &&
          exists("CO2_2076_2100_RCP85_CMIP5", inherits = TRUE) &&
          exists("CO2_2076_2100_RCP85_CMIP6", inherits = TRUE)) {
        
        dln_cmip5 <- log(CO2_2076_2100_RCP85_CMIP5 / CO2_1981_2005)
        dln_cmip6 <- log(CO2_2076_2100_RCP85_CMIP6 / CO2_1981_2005)
        
        dplyr::if_else(
          PERIOD == "2076_2100" & grepl("CMIP5", ensemble),
          dln_cmip5,
          dplyr::if_else(
            PERIOD == "2076_2100" & grepl("CMIP6", ensemble),
            dln_cmip6,
            0
          )
        )
      } else {
        rep(0, dplyr::n())
      }
    },
    
    inv_sqrt_VPD = 1 / sqrt(VPD),
    
    # K_ET via single source of truth
    K_ET = recompute_K_ET(dplyr::cur_data())
  )

# -----------------------
# Training subset
# -----------------------
# For the hybrid g_eff-residual RF you described, you typically need:
# g_eff, VPD, Ta, P, Rg, CO2_term, ET (ET not needed for RF fit but used for diagnostics)
df_opt <- df |>
  dplyr::filter(
    is.finite(g_eff), is.finite(ET),
    is.finite(VPD), is.finite(Ta),
    is.finite(P), is.finite(Rg),
    is.finite(CO2_term)
  )

# If you truly need A later, uncomment:
# df_opt <- df_opt |> dplyr::filter(is.finite(A))

stopifnot(nrow(df_opt) > 0)

# -----------------------
# RF tuning helpers (NA-robust)
# -----------------------
tune_rf_oob <- function(df_rf, target_col,
                        num.trees = 2000L, seed = 123,
                        mtry_vals = NULL,
                        min.node.size = c(3L, 5L, 10L, 20L),
                        sample.fraction = c(0.6, 0.8, 1.0),
                        splitrule = c("variance", "extratrees"),
                        max.depth = c(0L, 12L, 20L)) {
  set.seed(seed)
  stopifnot(target_col %in% names(df_rf))
  
  preds <- setdiff(names(df_rf), target_col)
  p <- length(preds)
  if (p < 1) stop("Need at least 1 predictor column for RF tuning.")
  
  if (is.null(mtry_vals)) {
    mtry_vals <- sort(unique(pmax(
      1L, pmin(p, c(1L, floor(sqrt(p)), max(1L, floor(0.33 * p)), max(1L, floor(0.5 * p)), p))
    )))
  }
  
  best <- list(rmse = Inf, fit = NULL, params = NULL)
  rows <- list()
  idx  <- 0L
  form <- as.formula(sprintf("%s ~ .", target_col))
  
  for (mtry in mtry_vals)
    for (minns in min.node.size)
      for (sf in sample.fraction)
        for (sr in splitrule)
          for (md in max.depth) {
            
            idx <- idx + 1L
            mtry_use <- max(1L, min(mtry, p))
            
            fit <- try(
              ranger::ranger(
                formula = form, data = df_rf,
                num.trees = num.trees,
                mtry = mtry_use,
                min.node.size = minns,
                sample.fraction = sf,
                replace = TRUE,
                splitrule = sr,
                max.depth = md,
                importance = "impurity",
                seed = seed + idx,
                oob.error = TRUE
              ),
              silent = TRUE
            )
            
            rmse <- NA_real_
            if (!inherits(fit, "try-error")) {
              pe <- suppressWarnings(fit$prediction.error)
              if (is.finite(pe)) rmse <- sqrt(pe)
            }
            
            rows[[idx]] <- tibble::tibble(
              mtry = mtry_use, min.node.size = minns, sample.fraction = sf,
              splitrule = sr, max.depth = md, oob_rmse = rmse
            )
            
            if (is.finite(rmse) && rmse < best$rmse) {
              best <- list(
                rmse = rmse,
                fit = fit,
                params = list(
                  mtry = mtry_use, min.node.size = minns,
                  sample.fraction = sf, splitrule = sr, max.depth = md
                )
              )
            }
          }
  
  grid_results <- dplyr::bind_rows(rows)
  
  if (is.null(best$fit)) {
    warning("RF tuner found no finite OOB RMSE; falling back to defaults.")
    fallback_mtry <- max(1L, floor(sqrt(p)))
    
    fit <- ranger::ranger(
      formula = form, data = df_rf,
      num.trees = num.trees,
      mtry = fallback_mtry,
      min.node.size = 5L,
      sample.fraction = 0.8,
      replace = TRUE,
      splitrule = "variance",
      max.depth = 0L,
      importance = "impurity",
      seed = seed,
      oob.error = TRUE
    )
    
    best <- list(
      rmse = sqrt(fit$prediction.error),
      fit = fit,
      params = list(
        mtry = fallback_mtry, min.node.size = 5L,
        sample.fraction = 0.8, splitrule = "variance", max.depth = 0L
      )
    )
  }
  
  list(
    best_fit = best$fit,
    best_rmse = best$rmse,
    best_params = best$params,
    grid_results = grid_results
  )
}

fit_rf_fixed <- function(df_rf, target_col, num.trees = 1200L, seed = 123) {
  stopifnot(target_col %in% names(df_rf))
  preds <- setdiff(names(df_rf), target_col)
  
  form <- as.formula(sprintf("%s ~ .", target_col))
  
  ranger::ranger(
    formula = form, data = df_rf,
    num.trees = num.trees,
    mtry = max(1L, min(3L, length(preds))),
    min.node.size = 5L,
    sample.fraction = 1.0,
    replace = TRUE,
    splitrule = "variance",
    max.depth = 0L,
    importance = "impurity",
    seed = seed,
    oob.error = TRUE
  )
}
