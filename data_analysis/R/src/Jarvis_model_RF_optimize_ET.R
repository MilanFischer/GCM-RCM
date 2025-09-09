suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(ranger)
  library(DEoptim)
  library(patchwork)
})

# ===============================
# 0) Knobs, guards, prerequisites
# ===============================
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# ---- FLEXIBLE predictors for the RF-residual model ----------------------------
# Put ANY columns you want here (must exist in df_opt): e.g., c("Rn","RH","Ta","CO2_term")
rf_candidate_vars <- c("Rg","Ta","P","CO2_term")

# Optional simple interaction features among candidate predictors (products)
use_interactions <- FALSE   # set TRUE to add pairwise products of numeric candidates

# ---- RF hyperparam strategy ---------------------------------------------------
# "none"  -> fixed baseline settings
# "quick" -> small grid
# "full"  -> larger grid (slower)
rf_tune_mode   <- "full"
rf_num_trees   <- 1200L     # used for both tuned and fixed fits

# ---- GENERAL variable gate knobs (applies to ALL predictors) ------------------
# Permutation gate: keep if (CI lower bound > 0) AND (relative ΔOOB >= threshold)
var_gate_when      <- "insignificant"   # "insignificant" or "significant" (rare)
perm_B             <- 1000               # repeats for permutation importance
perm_ci_levels     <- c(0.05, 0.95)     # CI for permutation deltas
min_rel_keep_pct   <- 1.0               # practical effect threshold (% ΔOOB RMSE)

# Impurity gate: also keep if impurity importance is "large enough"
use_impurity_gate  <- TRUE
imp_keep_frac_max  <- 0.05              # keep if imp >= 5% of max impurity importance (pre-gate)
imp_keep_top_k     <- Inf               # or keep top-K by impurity (set to a number to cap)

# ---- Optional "always keep" & optional drop-column CV guard -------------------
protect_vars         <- character(0)    # always-keep, e.g. c("Ta")
use_dropcol_cv_guard <- FALSE           # TRUE to turn on adaptive CV guard (slow)
# (Guard code is present but off by default; enable if needed.)

# ---- CO2 physics SIGN gate ----------------------------------------------------
# Keep CO2 only if ALE_geff over observed CO2 range is overall negative
require_negative_CO2_effect <- TRUE
co2_ale_bins                <- 50L
co2_neg_majority_frac       <- 0.60     # at least 60% local effects negative
co2_neg_total_tol           <- -1e-6    # require total ALE change < this (overall negative)

# ---- Data frames --------------------------------------------------------------
# df_opt: training/optimization; df: evaluation/prediction
# If you have only one data set, do: df_opt <- df
stopifnot(exists("df"), exists("df_opt"))
stopifnot(all(c("VPD","Ta","ET","g_eff") %in% names(df_opt)))
stopifnot(all(c("VPD","Ta","ET","g_eff") %in% names(df)))

# ---- Physics helper: K_ET per row --------------------------------------------
recompute_K_ET <- function(d) {
  stopifnot(all(c("VPD","Ta") %in% names(d)))
  vVPD <- pmax(d$VPD, eps); vTa <- d$Ta
  (rhoAir * CpAir / gamma) * vVPD * (1/1000) *
    1 / ((2.501e6 - 2361 * vTa) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

# Derive inv_sqrt_VPD in both frames
df_opt <- df_opt %>% mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))
df     <- df     %>% mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)))

# Optional: simple interactions among candidate predictors (in df_opt, and later used in df)
if (use_interactions) {
  num_cands <- intersect(rf_candidate_vars, names(df_opt))
  num_cands <- num_cands[vapply(df_opt[num_cands], is.numeric, logical(1))]
  if (length(num_cands) >= 2) {
    pairs <- combn(num_cands, 2, simplify = FALSE)
    for (pr in pairs) {
      nm <- paste0(pr[1], "_x_", pr[2])
      if (!nm %in% names(df_opt)) df_opt[[nm]] <- df_opt[[pr[1]]] * df_opt[[pr[2]]]
      if (!nm %in% names(df))     df[[nm]]     <- df[[pr[1]]]     * df[[pr[2]]]
      rf_candidate_vars <- unique(c(rf_candidate_vars, nm))
    }
  }
}

# ============================================
# 1) Parametric VPD-only model (optimizing ET)
# ============================================
obj_fun_geff_vpd <- function(par, input_df) {
  b0 <- par[[1]]; b1 <- par[[2]]
  inv_sqrt_vpd <- 1 / sqrt(pmax(input_df$VPD, eps))
  g_hat        <- b0 + b1 * inv_sqrt_vpd
  pen <- 0  # set to 1e6*mean(pmax(-g_hat,0)^2) to enforce g_hat >= 0
  K_ET  <- recompute_K_ET(input_df)
  ET_hat<- K_ET * g_hat
  sqrt(mean((input_df$ET - ET_hat)^2, na.rm = TRUE)) + pen
}

set.seed(if (exists("rng_seed")) rng_seed else 123)
de_geff <- DEoptim::DEoptim(
  fn       = obj_fun_geff_vpd,
  lower    = c(-10, -1e3),
  upper    = c( 10,  1e3),
  input_df = df_opt,
  control  = DEoptim::DEoptim.control(
    NP = 60L, itermax = 2000L, reltol = 1e-8,
    trace = if (exists("trace_flag")) trace_flag else FALSE
  )
)
b_hat <- de_geff$optim$bestmem
names(b_hat) <- c("b0_VPD","b1_VPD")
message(sprintf("Parametric fit (minimizing ET): b0=%.6f, b1=%.6f", b_hat[1], b_hat[2]))

# Parametric predictions on full df
g_eff_vpd_hat <- with(df, b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD)
K_ET_full     <- recompute_K_ET(df)
ET_vpd_hat    <- K_ET_full * g_eff_vpd_hat

# ================================================
# 2) Random Forest residual model on ET residuals
#     Fit RF to: resid_ET = ET_obs - ET_vpd_hat
# ================================================
rf_vars_present <- intersect(rf_candidate_vars, names(df_opt))
if (!length(rf_vars_present)) stop("None of rf_candidate_vars found in df_opt.")
if (!all(rf_vars_present %in% names(df))) {
  message("Note: some candidate predictors not present in df; final prediction uses available columns only.")
}

K_ET_opt    <- recompute_K_ET(df_opt)
g_param_opt <- with(df_opt, b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD)
ET_vpd_opt  <- K_ET_opt * g_param_opt

df_rf <- df_opt %>%
  mutate(resid_ET = ET - ET_vpd_opt) %>%
  select(resid_ET, all_of(rf_vars_present)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

# ---------- RF tuning helpers (NA-robust) ----------
tune_rf_oob <- function(df_rf, target_col,
                        num.trees = 2000L, seed = 123,
                        mtry_vals = NULL,
                        min.node.size = c(3L, 5L, 10L, 20L),
                        sample.fraction = c(0.6, 0.8, 1.0),
                        splitrule = c("variance", "extratrees"),
                        max.depth = c(0L, 12L, 20L)) {
  set.seed(seed)
  preds <- setdiff(names(df_rf), target_col); p <- length(preds)
  if (is.null(mtry_vals)) {
    mtry_vals <- sort(unique(pmax(1L, pmin(p, c(1L, floor(sqrt(p)), max(1L, floor(0.33*p)), max(1L, floor(0.5*p)), p)))))
  }
  best <- list(rmse = Inf, fit = NULL, params = NULL)
  rows <- list(); idx <- 0L; form <- as.formula(sprintf("%s ~ .", target_col))
  for (mtry in mtry_vals)
    for (minns in min.node.size)
      for (sf in sample.fraction)
        for (sr in splitrule)
          for (md in max.depth) {
            idx <- idx + 1L
            fit <- try(ranger(formula = form, data = df_rf,
                              num.trees = num.trees,
                              mtry = mtry, min.node.size = minns,
                              sample.fraction = sf, replace = TRUE, # define OOB
                              splitrule = sr, max.depth = md,
                              importance = "impurity",
                              seed = seed + idx, oob.error = TRUE), silent = TRUE)
            rmse <- NA_real_
            if (!inherits(fit, "try-error")) {
              pe <- suppressWarnings(fit$prediction.error)
              if (is.finite(pe)) rmse <- sqrt(pe)
            }
            rows[[idx]] <- tibble(mtry = mtry, min.node.size = minns, sample.fraction = sf,
                                  splitrule = sr, max.depth = md, oob_rmse = rmse)
            if (is.finite(rmse) && rmse < best$rmse) {
              best <- list(rmse = rmse, fit = fit,
                           params = list(mtry = mtry, min.node.size = minns,
                                         sample.fraction = sf, splitrule = sr, max.depth = md))
            }
          }
  grid_results <- dplyr::bind_rows(rows)
  if (is.null(best$fit)) {
    warning("RF tuner found no finite OOB RMSE; falling back to defaults.")
    fallback_mtry <- max(1L, floor(sqrt(p)))
    fit <- ranger(formula = form, data = df_rf, num.trees = num.trees,
                  mtry = fallback_mtry, min.node.size = 5L,
                  sample.fraction = 0.8, replace = TRUE, splitrule = "variance", max.depth = 0L,
                  importance = "impurity", seed = seed, oob.error = TRUE)
    best <- list(rmse = sqrt(fit$prediction.error), fit = fit,
                 params = list(mtry = fallback_mtry, min.node.size = 5L,
                               sample.fraction = 0.8, splitrule = "variance", max.depth = 0L))
  }
  list(best_fit = best$fit, best_rmse = best$rmse,
       best_params = best$params, grid_results = grid_results)
}

fit_rf_fixed <- function(df_rf, target_col, num.trees = 1200L, seed = 123) {
  preds <- setdiff(names(df_rf), target_col)
  form <- as.formula(sprintf("%s ~ .", target_col))
  ranger(formula = form, data = df_rf, num.trees = num.trees,
         mtry = min(3L, length(preds)), min.node.size = 5L,
         sample.fraction = 1.0, replace = TRUE, splitrule = "variance", max.depth = 0L,
         importance = "impurity", seed = seed, oob.error = TRUE)
}

# ---------- train initial RF (or tuned) ----------
set.seed(if (exists("rng_seed")) rng_seed else 123)
if (rf_tune_mode == "none") {
  rf_fit <- fit_rf_fixed(df_rf, target_col = "resid_ET", num.trees = rf_num_trees,
                         seed = if (exists("rng_seed")) rng_seed else 123)
  best_params <- list(mtry = rf_fit$mtry, min.node.size = rf_fit$min.node.size,
                      sample.fraction = 1.0, splitrule = "variance", max.depth = 0L)
} else {
  if (rf_tune_mode == "quick") {
    tuned <- tune_rf_oob(df_rf, "resid_ET", num.trees = rf_num_trees,
                         min.node.size = c(3L,5L,10L), sample.fraction = c(0.8,1.0),
                         splitrule = c("variance","extratrees"), max.depth = c(0L,12L))
  } else {
    tuned <- tune_rf_oob(df_rf, "resid_ET", num.trees = rf_num_trees)  # full grid
  }
  rf_fit     <- tuned$best_fit
  best_params<- tuned$best_params
}
base_oob_rmse <- sqrt(rf_fit$prediction.error)
message(sprintf("Initial RF (all candidates) OOB RMSE = %.3f", base_oob_rmse))
message(sprintf("Initial RF params: mtry=%s, min.node.size=%s, sample.fraction=%.2f, splitrule=%s, max.depth=%s",
                best_params$mtry, best_params$min.node.size, best_params$sample.fraction,
                as.character(best_params$splitrule), best_params$max.depth))

# ---- VIP (initial, BEFORE gating) — shows all candidates (incl. CO2_term) ----
if (!is.null(rf_fit$variable.importance)) {
  vi_init <- tibble(
    feature    = names(rf_fit$variable.importance),
    importance = as.numeric(rf_fit$variable.importance)
  ) %>% arrange(desc(importance))
  p_vip_init <- ggplot(vi_init, aes(x = reorder(feature, importance), y = importance)) +
    geom_col() + coord_flip() + theme_bw() +
    labs(x = NULL, y = "RF importance (impurity)",
         title = "ET-residual RF: Variable Importance (initial, pre-gate)")
  print(p_vip_init)
}

# ====================================================================
# 2a) (Optional) drop-column CV guard to force-keep some vars
# ====================================================================
guard_keep_vars <- character(0)
if (isTRUE(use_dropcol_cv_guard) && length(protect_vars)) {
  dropcol_cv_for <- function(df_rf, drop_var, seed = 2024) {
    n <- nrow(df_rf); K <- max(3L, min(10L, as.integer(round(sqrt(n)))))
    R <- if (n >= 5000) 2L else 3L
    deltas <- numeric(K * R); set.seed(seed)
    for (r in seq_len(R)) {
      fold_ids <- sample(rep(1:K, length.out = n))
      for (k in seq_len(K)) {
        test_idx  <- which(fold_ids == k)
        train_idx <- which(fold_ids != k)
        df_tr_with <- df_rf[train_idx, , drop = FALSE]
        df_tr_wo   <- df_tr_with[, setdiff(names(df_tr_with), drop_var), drop = FALSE]
        df_te      <- df_rf[test_idx, , drop = FALSE]
        df_te_wo   <- df_te[, intersect(names(df_te), names(df_tr_wo)), drop = FALSE]
        rf_w <- ranger(resid_ET ~ ., data = df_tr_with, num.trees = rf_num_trees,
                       mtry = min(3, ncol(df_tr_with) - 1L),
                       min.node.size = 5, seed = seed + r*100 + k)
        rf_wo <- ranger(resid_ET ~ ., data = df_tr_wo, num.trees = rf_num_trees,
                        mtry = min(3, ncol(df_tr_wo) - 1L),
                        min.node.size = 5, seed = seed + r*100 + k + 1)
        y_true   <- df_te$resid_ET
        y_pred_w <- as.numeric(predict(rf_w,  data = df_te[, setdiff(names(df_te), "resid_ET"), drop = FALSE])$predictions)
        y_pred_wo<- as.numeric(predict(rf_wo, data = df_te_wo[, setdiff(names(df_te_wo), "resid_ET"), drop = FALSE])$predictions)
        rmse_w  <- sqrt(mean((y_true - y_pred_w )^2, na.rm = TRUE))
        rmse_wo <- sqrt(mean((y_true - y_pred_wo)^2, na.rm = TRUE))
        deltas[(r-1)*K + k] <- rmse_wo - rmse_w
      }
    }
    tibble(mean_delta = mean(deltas, na.rm = TRUE),
           median = median(deltas, na.rm = TRUE),
           q05 = quantile(deltas, 0.05, na.rm = TRUE),
           q95 = quantile(deltas, 0.95, na.rm = TRUE),
           K = K, R = R)
  }
  for (pv in intersect(protect_vars, colnames(df_rf))) {
    cvr <- dropcol_cv_for(df_rf, pv)
    message(sprintf("%s drop-column CV: mean ΔRMSE(without - with)=%.3f (5–95%% %.3f–%.3f) K=%d R=%d",
                    pv, cvr$mean_delta, cvr$q05, cvr$q95, cvr$K, cvr$R))
    if (is.finite(cvr$mean_delta) && cvr$mean_delta > 0) {
      guard_keep_vars <- union(guard_keep_vars, pv)
    }
  }
  if (length(guard_keep_vars)) message("Guard-keep variables: ", paste(guard_keep_vars, collapse = ", "))
}

# ====================================================================
# 2b) Permutation importance for ALL predictors (pre-gate, tuned params)
# ====================================================================
# ---- Permutation importance helper (robust; supports any target col) ----
perm_importance_all <- function(df_rf,
                                base_oob_rmse = NULL,
                                B = 100,
                                seed = 1000,
                                target_col = "resid_ET") {
  stopifnot(target_col %in% names(df_rf))
  preds <- setdiff(names(df_rf), target_col)
  if (!length(preds)) {
    return(tibble(feature = character(0), mean_delta = numeric(0),
                  lo = numeric(0), hi = numeric(0), rel = numeric(0), B = integer(0)))
  }
  
  # If baseline OOB RMSE not provided, estimate it with a fresh RF
  if (is.null(base_oob_rmse)) {
    rf0 <- ranger::ranger(
      formula = as.formula(paste(target_col, "~ .")),
      data = df_rf,
      num.trees = 800,
      mtry = min(3, length(preds)),
      min.node.size = 5,
      seed = seed,
      oob.error = TRUE
    )
    base_oob_rmse <- sqrt(rf0$prediction.error)
  }
  
  purrr::map_dfr(preds, function(v) {
    set.seed(seed)
    deltas <- numeric(B)
    for (b in seq_len(B)) {
      dfp <- df_rf
      dfp[[v]] <- sample(dfp[[v]]) # permute v
      rf_p <- ranger::ranger(
        formula = as.formula(paste(target_col, "~ .")),
        data = dfp,
        num.trees = 800,
        mtry = min(3, length(preds)),
        min.node.size = 5,
        seed = seed + b,
        oob.error = TRUE
      )
      deltas[b] <- sqrt(rf_p$prediction.error) - base_oob_rmse
    }
    md <- mean(deltas, na.rm = TRUE)
    tibble::tibble(
      feature    = v,
      mean_delta = md,
      lo         = stats::quantile(deltas, 0.05, na.rm = TRUE),
      hi         = stats::quantile(deltas, 0.95, na.rm = TRUE),
      rel        = 100 * md / base_oob_rmse,
      B          = B
    )
  }) %>% dplyr::arrange(dplyr::desc(mean_delta))
}

perm_all <- perm_importance_all(
  df_rf,
  base_oob_rmse = base_oob_rmse,
  B = perm_B,
  target_col = "resid_ET"
)
cat("\n=== Permutation importance (pre-gate) ===\n")
print(perm_all %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))))

# Plot with CI (pre-gate)
p_perm <- ggplot(perm_all, aes(x = reorder(feature, mean_delta), y = mean_delta)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  coord_flip() + theme_bw() +
  labs(x = NULL, y = expression(Delta*" OOB RMSE (permute)"),
       title = "Permutation importance (global, pre-gate)",
       subtitle = sprintf("Keep rule: lo>0 & rel≥%.1f%%", min_rel_keep_pct))
print(p_perm)

# ====================================================================
# 2c) Hybrid variable gate → choose predictors
#      Keep if (perm passes) OR (impurity passes); also keep protected/guard.
# ====================================================================
# Perm gate
if (tolower(var_gate_when) == "insignificant") {
  keep_perm <- perm_all %>% filter(lo > 0, rel >= min_rel_keep_pct) %>% pull(feature)
} else {
  drop_perm <- perm_all %>% filter(lo > 0, rel >= min_rel_keep_pct) %>% pull(feature)
  keep_perm <- setdiff(perm_all$feature, drop_perm)
}

# Impurity gate (pre-gate VIP)
keep_imp <- character(0)
if (exists("vi_init")) {
  thr <- imp_keep_frac_max * max(vi_init$importance, na.rm = TRUE)
  keep_imp <- vi_init %>% filter(importance >= thr) %>% arrange(desc(importance)) %>% pull(feature)
  if (is.finite(imp_keep_top_k)) keep_imp <- head(keep_imp, imp_keep_top_k)
}

vars_keep <- union(keep_perm, keep_imp)
vars_keep <- union(vars_keep, intersect(protect_vars, colnames(df_rf)))
vars_keep <- union(vars_keep, guard_keep_vars)

# ====================================================================
# 2d) CO2 physics SIGN gate (drop CO2 if it does NOT reduce g_eff)
# ====================================================================
if (require_negative_CO2_effect && "CO2_term" %in% vars_keep) {
  # local predictor using current rf_fit & rf_vars_present:
  predict_geff_ET_local <- function(newdata, rf_fit_current, rf_vars) {
    inv_sqrt_VPD <- 1 / sqrt(pmax(newdata$VPD, eps))
    g_param      <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD
    K_ET_new     <- recompute_K_ET(newdata)
    vars <- intersect(rf_vars, colnames(newdata))
    rf_resid_ET  <- if (length(vars)) as.numeric(predict(rf_fit_current, data = newdata[, vars, drop = FALSE])$predictions) else 0
    ET_hat       <- K_ET_new * g_param + rf_resid_ET
    g_eff_hat    <- ET_hat / K_ET_new
    list(geff = as.numeric(g_eff_hat), ET = as.numeric(ET_hat))
  }
  co2_ale_check <- function(var = "CO2_term", data = df_opt, n_bins = 50L) {
    z_all <- data[[var]]; z <- z_all[is.finite(z_all)]
    if (!length(z)) return(list(ok = FALSE, msg = "CO2 has no finite values"))
    edges <- as.numeric(quantile(z, probs = seq(0,1,length.out=n_bins+1), na.rm=TRUE))
    edges[1] <- min(z, na.rm=TRUE); edges[length(edges)] <- max(z, na.rm=TRUE)
    dG <- numeric(n_bins); cnt <- integer(n_bins)
    for (j in seq_len(n_bins)) {
      lo <- edges[j]; hi <- edges[j+1]
      idx <- which(z_all >= lo & (z_all < hi | (j == n_bins & z_all <= hi)))
      cnt[j] <- length(idx); if (cnt[j]==0) { dG[j] <- 0; next }
      df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
      df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
      pred_lo <- predict_geff_ET_local(df_lo, rf_fit, rf_vars_present)
      pred_hi <- predict_geff_ET_local(df_hi, rf_fit, rf_vars_present)
      dG[j]   <- mean(pred_hi$geff - pred_lo$geff, na.rm = TRUE)
    }
    ale_G <- cumsum(dG); N <- sum(cnt)
    if (N > 0) ale_G <- ale_G - sum(ale_G * cnt, na.rm = TRUE) / N
    total_change <- tail(ale_G, 1) - ale_G[1]
    frac_neg     <- mean(dG < 0, na.rm = TRUE)
    list(ok=TRUE, total_change=total_change, frac_neg=frac_neg)
  }
  rf_vars_present <- rf_vars_present
  chk <- co2_ale_check(n_bins = co2_ale_bins)
  if (isTRUE(chk$ok)) {
    message(sprintf("CO2 sign check (ALE of g_eff): total Δ=%.4g, frac(negative increments)=%.1f%%",
                    chk$total_change, 100*chk$frac_neg))
    if (!(chk$total_change <= co2_neg_total_tol && chk$frac_neg >= co2_neg_majority_frac)) {
      message("Decision: CO2_term does NOT reduce g_eff → dropping CO2_term by physics sign gate.")
      vars_keep <- setdiff(vars_keep, "CO2_term")
    } else {
      message("Decision: CO2_term reduces g_eff → keeping CO2_term (physics sign gate passed).")
    }
  } else {
    message("CO2 sign check skipped: ", chk$msg)
  }
}

# =========================
# 2e) Final RF (post gate)
# =========================
if (!length(vars_keep)) {
  message("Gate selected no predictors; keeping initial RF unchanged.")
  # rf_fit, best_params, rf_vars_present already defined
} else {
  rf_vars_present <- vars_keep
  df_rf <- df_opt %>%
    mutate(resid_ET = ET - ET_vpd_opt) %>%
    select(resid_ET, all_of(rf_vars_present)) %>%
    filter(if_all(where(is.numeric), ~ is.finite(.x)))
  
  # Refit with tuned (or fixed) params for consistency
  form <- resid_ET ~ .
  rf_fit <- ranger(
    formula = form, data = df_rf,
    num.trees = rf_num_trees,
    mtry = min(best_params$mtry, max(1L, length(rf_vars_present))),  # guard
    min.node.size = best_params$min.node.size,
    sample.fraction = best_params$sample.fraction,
    replace = TRUE,
    splitrule = as.character(best_params$splitrule),
    max.depth = best_params$max.depth,
    importance = "impurity",
    seed = if (exists("rng_seed")) rng_seed else 123,
    oob.error = TRUE
  )
  base_oob_rmse <- sqrt(rf_fit$prediction.error)
  message(sprintf("Final RF (gated predictors: %s) OOB RMSE = %.3f",
                  paste(rf_vars_present, collapse = ", "), base_oob_rmse))
  
  # VIP (final, after gate)
  if (!is.null(rf_fit$variable.importance)) {
    vi_df_final <- tibble(
      feature    = names(rf_fit$variable.importance),
      importance = as.numeric(rf_fit$variable.importance)
    ) %>% arrange(desc(importance))
    p_vip_final <- ggplot(vi_df_final,
                          aes(x = reorder(feature, importance), y = importance)) +
      geom_col() + coord_flip() + theme_bw() +
      labs(x = NULL, y = "RF importance (impurity)",
           title = "ET-residual RF: Variable Importance (final, post-gate)")
    print(p_vip_final)
  }
}

# ---- Safe RF prediction helper (ET residuals) using *current* rf_fit/vars ----
.predict_rf_resid_ET <- function(newdata) {
  vars <- intersect(rf_vars_present, colnames(newdata))
  if (length(vars) == 0) return(rep(0, nrow(newdata)))
  as.numeric(predict(rf_fit, data = newdata[, vars, drop = FALSE])$predictions)
}

# Predict ET residuals on full df
rf_pred_ET <- .predict_rf_resid_ET(df)

perm_final <- perm_importance_all(
  df_rf,
  base_oob_rmse = sqrt(rf_fit$prediction.error),
  B = perm_B,
  target_col = "resid_ET"
)
cat("\n=== Permutation importance (FINAL RF, post-gate) ===\n")
print(perm_final %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))))

# ================================
# 3) Final predictions & metrics
# ================================
ET_final    <- ET_vpd_hat + rf_pred_ET
g_eff_final <- ET_final / K_ET_full

output_df_geff <- df %>%
  mutate(
    g_eff_vpd_hat = g_eff_vpd_hat,
    ET_vpd_hat    = ET_vpd_hat,
    resid_ET      = ET - ET_vpd_hat,
    resid_rf_ET   = rf_pred_ET,
    ET_final      = ET_final,
    g_eff_final   = g_eff_final,
    g_eff_resids  = g_eff - g_eff_final
  )

RMSE_vpd_geff   <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_eff_vpd_hat)^2, na.rm = TRUE))
R2_vpd_geff     <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_eff_vpd_hat, use = "complete.obs")^2)
RMSE_final_geff <- sqrt(mean((output_df_geff$g_eff - output_df_geff$g_eff_final  )^2, na.rm = TRUE))
R2_final_geff   <- suppressWarnings(cor(output_df_geff$g_eff, output_df_geff$g_eff_final,   use = "complete.obs")^2)

RMSE_vpd_ET   <- sqrt(mean((output_df_geff$ET - output_df_geff$ET_vpd_hat)^2, na.rm = TRUE))
R2_vpd_ET     <- suppressWarnings(cor(output_df_geff$ET, output_df_geff$ET_vpd_hat, use = "complete.obs")^2)
RMSE_final_ET <- sqrt(mean((output_df_geff$ET - output_df_geff$ET_final  )^2, na.rm = TRUE))
R2_final_ET   <- suppressWarnings(cor(output_df_geff$ET, output_df_geff$ET_final,   use = "complete.obs")^2)

message(sprintf("g_eff — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Stage2(Final): RMSE=%.3f, R^2=%.3f",
                RMSE_vpd_geff, R2_vpd_geff, RMSE_final_geff, R2_final_geff))
message(sprintf("ET    — Stage1(VPD-only): RMSE=%.3f, R^2=%.3f  |  Stage2(Final): RMSE=%.3f, R^2=%.3f",
                RMSE_vpd_ET,   R2_vpd_ET,   RMSE_final_ET,   R2_final_ET))

# ================================
# 4) Plots: observed vs predicted
# ================================
p_geff <- output_df_geff %>%
  transmute(g_eff_obs = g_eff, geff_final = g_eff_final) %>%
  filter(is.finite(g_eff_obs), is.finite(geff_final)) %>%
  ggplot(aes(g_eff_obs, geff_final)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(
    x = expression(paste("Observed ", g[eff])),
    y = expression(paste("Predicted ", g[eff], " (Final)")),
    title = expression(paste("Observed vs Predicted ", g[eff])),
    subtitle = sprintf("Stage1: RMSE=%.3f R^2=%.3f   |   Stage2: RMSE=%.3f R^2=%.3f",
                       RMSE_vpd_geff, R2_vpd_geff, RMSE_final_geff, R2_final_geff)
  )

p_et <- output_df_geff %>%
  transmute(ET_obs = ET, ET_final = ET_final) %>%
  filter(is.finite(ET_obs), is.finite(ET_final)) %>%
  ggplot(aes(ET_obs, ET_final)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(
    x = "Observed ET",
    y = "Predicted ET (Final)",
    title = "Observed vs Predicted ET",
    subtitle = sprintf("Stage1: RMSE=%.3f R^2=%.3f   |   Stage2: RMSE=%.3f R^2=%.3f",
                       RMSE_vpd_ET, R2_vpd_ET, RMSE_final_ET, R2_final_ET)
  )

print(p_geff + p_et)

# ====================================================
# 5) ALE (Accumulated Local Effects), FINAL RF
#     (show only variables that actually PASSED the gate + VPD for context)
# ====================================================
.predict_rf_resid_ET <- function(newdata) {
  vars <- intersect(rf_vars_present, colnames(newdata))
  if (length(vars) == 0) return(rep(0, nrow(newdata)))
  as.numeric(predict(rf_fit, data = newdata[, vars, drop = FALSE])$predictions)
}
predict_geff_ET <- function(newdata) {
  inv_sqrt_VPD <- 1 / sqrt(pmax(newdata$VPD, eps))
  g_param      <- b_hat[["b0_VPD"]] + b_hat[["b1_VPD"]] * inv_sqrt_VPD
  K_ET_new     <- recompute_K_ET(newdata)
  rf_resid_ET  <- .predict_rf_resid_ET(newdata)
  ET_hat       <- K_ET_new * g_param + rf_resid_ET
  g_eff_hat    <- ET_hat / K_ET_new
  list(geff = as.numeric(g_eff_hat), ET = as.numeric(ET_hat))
}

ale_one <- function(var, data = df, n_bins = 50) {
  stopifnot(var %in% names(data))
  z_all <- data[[var]]; z <- z_all[is.finite(z_all)]
  if (!length(z)) stop("Variable ", var, " has no finite values.")
  edges <- as.numeric(quantile(z, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  edges[1] <- min(z, na.rm = TRUE); edges[length(edges)] <- max(z, na.rm = TRUE)
  dET <- numeric(n_bins); dG  <- numeric(n_bins); cnt <- integer(n_bins)
  for (j in seq_len(n_bins)) {
    lo <- edges[j]; hi <- edges[j + 1]
    idx <- which(z_all >= lo & (z_all < hi | (j == n_bins & z_all <= hi)))
    cnt[j] <- length(idx); if (cnt[j] == 0) { dET[j] <- 0; dG[j] <- 0; next }
    df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
    df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
    pred_lo <- predict_geff_ET(df_lo); pred_hi <- predict_geff_ET(df_hi)
    dET[j] <- mean(pred_hi$ET   - pred_lo$ET,   na.rm = TRUE)
    dG[j]  <- mean(pred_hi$geff - pred_lo$geff, na.rm = TRUE)
  }
  ale_ET <- cumsum(dET); ale_G <- cumsum(dG)
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  N <- sum(cnt)
  if (N > 0) {
    ale_ET <- ale_ET - sum(ale_ET * cnt, na.rm = TRUE) / N
    ale_G  <- ale_G  - sum(ale_G  * cnt, na.rm = TRUE) / N
  }
  tibble(var = var, x = mids, ALE_ET = ale_ET, ALE_geff = ale_G, n_in_bin = cnt)
}

# Only show VPD + variables that survived the gate in the final RF
vars_for_ALE <- unique(c("VPD", rf_vars_present))
vars_for_ALE <- intersect(vars_for_ALE, names(df))
stopifnot(length(vars_for_ALE) >= 1)

# Labels only for variables we will plot
var_labels <- setNames(vars_for_ALE, vars_for_ALE)
if ("VPD" %in% vars_for_ALE) var_labels["VPD"] <- "VPD (kPa)"
if ("Ta"  %in% vars_for_ALE) var_labels["Ta"]  <- "Ta (°C)"

n_bins <- 50L
ale_raw <- bind_rows(lapply(vars_for_ALE, ale_one, data = df, n_bins = n_bins))

# Anchor to absolute scale using model means
pred0 <- predict_geff_ET(df)
mu_ET_pred   <- mean(pred0$ET,   na.rm = TRUE)
mu_geff_pred <- mean(pred0$geff, na.rm = TRUE)

# Smooth & plot
.fit_loess_safe <- function(y, x, w, span) {
  try(loess(y ~ x, weights = w, span = span,
            control = loess.control(surface = "direct")), silent = TRUE)
}
.predict_at <- function(fit, x_train, y_train, x_new) {
  if (!inherits(fit, "try-error")) as.numeric(predict(fit, data.frame(x = x_new)))
  else as.numeric(approx(x = x_train, y = y_train, xout = x_new, rule = 2, ties = mean)$y)
}

smooth_span <- 0.5
smooth_grid_points <- 300L

smooth_one_var <- function(vname) {
  dv <- ale_raw %>% filter(var == vname)
  x_min <- min(df[[vname]], na.rm = TRUE); x_max <- max(df[[vname]], na.rm = TRUE)
  x_grid <- seq(x_min, x_max, length.out = smooth_grid_points)
  fit_ET <- .fit_loess_safe(dv$ALE_ET,   dv$x, dv$n_in_bin, smooth_span)
  fit_G  <- .fit_loess_safe(dv$ALE_geff, dv$x, dv$n_in_bin, smooth_span)
  y_ET_grid <- .predict_at(fit_ET, dv$x, dv$ALE_ET,   x_grid)
  y_G_grid  <- .predict_at(fit_G,  dv$x, dv$ALE_geff, x_grid)
  y_ET_mid <- .predict_at(fit_ET, dv$x, dv$ALE_ET,   dv$x)
  y_G_mid  <- .predict_at(fit_G,  dv$x, dv$ALE_geff, dv$x)
  off_ET <- weighted.mean(y_ET_mid, w = dv$n_in_bin, na.rm = TRUE)
  off_G  <- weighted.mean(y_G_mid,  w = dv$n_in_bin, na.rm = TRUE)
  tibble(
    var = vname, x = x_grid,
    ALE_ET_rel   = y_ET_grid - off_ET,
    ALE_geff_rel = y_G_grid  - off_G,
    ALE_ET_abs   = (y_ET_grid - off_ET) + mu_ET_pred,
    ALE_geff_abs = (y_G_grid  - off_G)  + mu_geff_pred
  )
}
ale_smoothed <- map_dfr(vars_for_ALE, smooth_one_var)

# Observed points for anchored panels
obs_max_points <- 5000L
make_obs_pts <- function(y_name) {
  out <- bind_rows(lapply(vars_for_ALE, function(v) {
    tibble(var = v, x = df[[v]], y = df[[y_name]])
  })) %>% filter(is.finite(x), is.finite(y))
  if (nrow(out) > obs_max_points) out <- out[sample.int(nrow(out), obs_max_points), , drop = FALSE]
  out
}
obs_pts_ET   <- make_obs_pts("ET")
obs_pts_geff <- make_obs_pts("g_eff")

plot_ale_et_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_ET_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(setNames(vars_for_ALE, vars_for_ALE))) +
  labs(x = NULL, y = "ET (mm/yr) — anchored to model mean",
       title = "ALE (ET): absolute scale") +
  theme_bw() +
  geom_point(data = obs_pts_ET, aes(x = x, y = y),
             inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25)

plot_ale_geff_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_geff_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(setNames(vars_for_ALE, vars_for_ALE))) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*") — anchored to model mean"),
       title = expression("ALE (" * g[eff] * "): absolute scale")) +
  theme_bw() +
  geom_point(data = obs_pts_geff, aes(x = x, y = y),
             inherit.aes = FALSE, shape = 3, size = 0.5, alpha = 0.25)

print(plot_ale_et_abs);  print(plot_ale_geff_abs)

# ==============================
# Save bundle (optional)
# ==============================
if (!exists("perm_all")) perm_all <- NULL
geff_vpd_rf_bundle <- list(
  metadata      = if (exists("metadata")) paste0(metadata, " | ET-optimized, RF on ET residuals + ALE (hybrid gate + optional tuning).") else
    "ET-optimized, RF on ET residuals + ALE (hybrid gate + optional tuning).",
  b0_VPD        = b_hat[["b0_VPD"]],
  b1_VPD        = b_hat[["b1_VPD"]],
  rf_fit        = rf_fit,
  rf_vars       = rf_vars_present,
  rf_params     = best_params,
  out_df        = output_df_geff,
  RMSE_vpd_geff = RMSE_vpd_geff,   R2_vpd_geff = R2_vpd_geff,
  RMSE_fin_geff = RMSE_final_geff, R2_fin_geff = R2_final_geff,
  RMSE_vpd_ET   = RMSE_vpd_ET,     R2_vpd_ET   = R2_vpd_ET,
  RMSE_fin_ET   = RMSE_final_ET,   R2_fin_ET   = R2_final_ET,
  rng_seed      = if (exists("rng_seed")) rng_seed else NA_integer_,
  ale_raw       = ale_raw,
  ale_smooth    = ale_smoothed,
  perm_all      = perm_all,
  perm_final    = perm_final,
  gate_mode     = var_gate_when,
  gate_rel_pct  = min_rel_keep_pct,
  impurity_gate = list(enabled = use_impurity_gate, frac_max = imp_keep_frac_max, top_k = imp_keep_top_k),
  protected     = protect_vars,
  guard_keep    = guard_keep_vars
)
if (exists("out_file")) {
  save(geff_vpd_rf_bundle, file = sub("\\.RData$", "_ETopt_RFresidET_ALE_VARgate_TUNED.RData", out_file))
}
