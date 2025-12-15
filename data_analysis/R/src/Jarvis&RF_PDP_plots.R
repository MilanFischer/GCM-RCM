# ============================================================
# Jarvis&RF_PDP_plots.R  (universal)
#   - Partial Dependence Plots (1-D) for ET and g_eff
#   - Others held fixed at mean/median (like ALE Option B)
# Saves ONLY TWO PLOTS (anchored + observed overlay), WITH model_id:
#   - PDP_ET_abs_<model_id>.png
#   - PDP_geff_abs_<model_id>.png
#
# Usage:
#   source("./src/Jarvis&RF_PDP_plots.R")
#   run_PDP_plots(bundle_path="./RData/20251214_jarvis_objects.RData")
#   run_PDP_plots(bundle_path="./RData/20251214_RF_objects.RData")
# ============================================================

run_PDP_plots <- function(
    bundle_path,
    vars = c("P","Rg","Ta","VPD"),
    hold_others = c("median","mean"),
    grid_points = 300,
    # observed overlay
    obs_max_points = 5000,
    # output
    plot_dir = "../Plots/ggplot2",
    width = 9,
    height = 6.5,
    dpi = 300
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(purrr)
    # For RF/ranger bundles: make sure predict.ranger is available
    suppressWarnings(suppressMessages(requireNamespace("ranger", quietly = TRUE)))
  })
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  hold_others <- match.arg(hold_others)
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -----------------------
  # Load bundle & detect object
  # -----------------------
  e0 <- new.env(parent = emptyenv())
  load(bundle_path, envir = e0)
  
  bundle <-
    if (exists("jarvis_bundle", envir = e0, inherits = FALSE)) get("jarvis_bundle", envir = e0) else
      if (exists("rf_hybrid_bundle", envir = e0, inherits = FALSE)) get("rf_hybrid_bundle", envir = e0) else
        stop("No recognized bundle object found in: ", bundle_path,
             "\nExpected: jarvis_bundle or rf_hybrid_bundle")
  
  stopifnot(is.list(bundle))
  stopifnot(is.function(bundle$pred_fun))
  stopifnot(!is.null(bundle$df_opt))
  
  model_id <- bundle$model_id %||% tools::file_path_sans_ext(basename(bundle_path))
  df_opt   <- bundle$df_opt
  
  # required response cols for observed overlay
  req_cols <- c("ET","g_eff")
  miss <- setdiff(req_cols, names(df_opt))
  if (length(miss)) stop("df_opt missing required columns: ", paste(miss, collapse = ", "))
  
  vars <- intersect(vars, names(df_opt))
  if (length(vars) < 1) stop("None of requested vars exist in df_opt.")
  
  message("Loaded bundle: ", model_id)
  message("Rows in df_opt: ", nrow(df_opt))
  message("PDP vars: ", paste(vars, collapse = ", "))
  
  # -----------------------
  # Unified predictor wrapper: pred_fun(newdata)-> list(ET, geff)
  # -----------------------
  pred_fun <- function(newdata) {
    out <- bundle$pred_fun(newdata)
    
    ET   <- out$ET   %||% out$ET_pred   %||% out$ET_hat
    geff <- out$geff %||% out$geff_pred %||% out$g_eff %||% out$g_eff_pred
    
    if (is.null(ET) || is.null(geff)) {
      stop("bundle$pred_fun(newdata) must return ET and geff (or compatible names).")
    }
    if (length(ET) != nrow(newdata) || length(geff) != nrow(newdata)) {
      stop("Prediction lengths do not match nrow(newdata).")
    }
    list(ET = as.numeric(ET), geff = as.numeric(geff))
  }
  
  # -----------------------
  # Reference row: others held fixed
  # -----------------------
  .mode <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  ref_row <- df_opt[1, , drop = FALSE]
  num_cols <- names(df_opt)[vapply(df_opt, is.numeric, logical(1))]
  
  if (hold_others == "mean") {
    ref_row[1, num_cols] <- lapply(df_opt[num_cols], function(col) mean(col, na.rm = TRUE))
  } else {
    ref_row[1, num_cols] <- lapply(df_opt[num_cols], function(col) median(col, na.rm = TRUE))
  }
  
  non_num_cols <- setdiff(names(df_opt), num_cols)
  if (length(non_num_cols)) {
    ref_modes <- lapply(df_opt[non_num_cols], .mode)
    for (nm in non_num_cols) ref_row[[nm]] <- ref_modes[[nm]]
  }
  
  # -----------------------
  # PDP for one variable (sweep var, others fixed to ref_row)
  # -----------------------
  make_pdp_one <- function(v) {
    x_min <- min(df_opt[[v]], na.rm = TRUE)
    x_max <- max(df_opt[[v]], na.rm = TRUE)
    x_grid <- seq(x_min, x_max, length.out = grid_points)
    
    d <- ref_row[rep(1, length(x_grid)), , drop = FALSE]
    d[[v]] <- x_grid
    
    preds <- pred_fun(d)
    
    tibble(
      var = v,
      x = x_grid,
      ET_pred   = preds$ET,
      geff_pred = preds$geff
    )
  }
  
  pdp <- map_dfr(vars, make_pdp_one)
  
  # -----------------------
  # Observed points (anchored overlay)
  # -----------------------
  make_obs_pts <- function(y_name) {
    out <- bind_rows(lapply(vars, function(v) {
      tibble(var = v, x = df_opt[[v]], y = df_opt[[y_name]])
    })) %>% filter(is.finite(x), is.finite(y))
    if (nrow(out) > obs_max_points) out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
    out
  }
  
  obs_pts_et   <- make_obs_pts("ET")
  obs_pts_geff <- make_obs_pts("g_eff")
  
  # Pretty facet labels
  var_labels <- c(
    P = "P (precipitation)",
    Rg = "Rg (radiation)",
    Ta = "Ta (°C)",
    VPD = "VPD (kPa)",
    CO2_term = "CO₂ term"
  )
  
  theme_base <- theme_bw()
  
  # -----------------------
  # TWO plots: ET anchored + g_eff anchored (absolute PDP)
  # -----------------------
  p_et_abs <- ggplot(pdp, aes(x = x, y = ET_pred)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(
      x = NULL, y = "ET (mm/yr)",
      title = paste0("PDP (ET) — others held fixed at ", hold_others, " | ", model_id)
    ) +
    theme_base +
    geom_point(
      data = obs_pts_et,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 3, size = 0.5, alpha = 0.25, color = "black"
    )
  
  p_geff_abs <- ggplot(pdp, aes(x = x, y = geff_pred)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ var, scales = "free_x", labeller = as_labeller(var_labels)) +
    labs(
      x = NULL, y = "g_eff (mm s^-1)",
      title = paste0("PDP (g_eff) — others held fixed at ", hold_others, " | ", model_id)
    ) +
    theme_base +
    geom_point(
      data = obs_pts_geff,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 3, size = 0.5, alpha = 0.25, color = "black"
    )
  
  print(p_et_abs)
  print(p_geff_abs)
  
  # filenames WITH model_id
  safe_id <- gsub("[^A-Za-z0-9_\\-]+", "_", model_id)
  f1 <- file.path(plot_dir, paste0("PDP_ET_abs_",   safe_id, ".png"))
  f2 <- file.path(plot_dir, paste0("PDP_geff_abs_", safe_id, ".png"))
  
  ggsave(f1, p_et_abs,   width = width, height = height, dpi = dpi)
  ggsave(f2, p_geff_abs, width = width, height = height, dpi = dpi)
  
  invisible(list(
    model_id = model_id,
    bundle_path = bundle_path,
    pdp = pdp,
    plots = list(ET_abs = p_et_abs, geff_abs = p_geff_abs),
    files = c(ET_abs = f1, geff_abs = f2)
  ))
}

# Optional CLI:
#   Rscript src/Jarvis&RF_PDP_plots.R path/to/bundle.RData
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    run_PDP_plots(bundle_path = args[[1]])
  } else {
    message("No bundle_path provided. Example:\n",
            "  Rscript src/Jarvis&RF_PDP_plots.R ./RData/20251214_jarvis_objects.RData")
  }
}
