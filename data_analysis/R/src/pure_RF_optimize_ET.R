suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(ranger)
  library(patchwork)
})

# ============================================================
# 0) Knobs, guards, prerequisites
# ============================================================
eps <- if (exists("eps", inherits = TRUE)) eps else 1e-6

# ---- FLEXIBLE user predictors for the ET RF (besides physics core) ------------
# Add any columns that exist in df_opt/df: e.g., c("Rn","RH","Ta","CO2_term")
rf_candidate_vars <- c("Rg","Ta","P","CO2_term")

# Optional: add simple pairwise interaction features among candidates
use_interactions <- FALSE  # TRUE to create products of numeric candidates

# RF hyperparameters (robust defaults)
rf_num_trees   <- 1200L
rf_mtry        <- NULL       # NULL -> sqrt(p); else integer
rf_min_node    <- 5L
rf_splitrule   <- "variance"
rf_max_depth   <- 0L         # 0 = unlimited depth
rf_sample_frac <- 0.8        # OOB with subsampling

# ---- Variable gate knobs (applies to NON-core vars) ---------------------------
var_gate_when    <- "insignificant"  # keep if useful; or "significant" (invert)
perm_B           <- 200              # ↑ for more stable CIs
perm_ci_levels   <- c(0.05, 0.95)
min_rel_keep_pct <- 1.0              # practical effect threshold (% ΔOOB RMSE)

use_impurity_gate <- TRUE            # also keep if impurity is large enough
imp_keep_frac_max <- 0.05            # ≥5% of max impurity passes
imp_keep_top_k    <- Inf             # or cap to top-K

# ---- Optional "always-keep" & drop-column CV guard (off by default) ----------
protect_vars         <- character(0) # e.g. c("Ta")
use_dropcol_cv_guard <- FALSE        # TRUE to enable (slow)

# ---- CO2 physics SIGN gate ----------------------------------------------------
# Keep CO2 ONLY if it reduces g_eff. We test sign via ALE on g_eff where
# g_eff_hat = ET_hat / K_ET (even though RF trains on ET).
require_negative_CO2_effect <- TRUE
co2_ale_bins                <- 50L
co2_neg_majority_frac       <- 0.60   # ≥60% local effects negative
co2_neg_total_tol           <- -1e-6  # overall ALE change must be negative

# ============================================================
# 1) Data frames & physics helpers
# ============================================================
# df_opt: training/calibration; df: evaluation/prediction
# If you have one dataset, do: df_opt <- df
stopifnot(exists("df"), exists("df_opt"))
stopifnot(all(c("VPD","Ta","ET","g_eff") %in% names(df_opt)))
stopifnot(all(c("VPD","Ta","ET","g_eff") %in% names(df)))

# ---- physics: K_ET(VPD, Ta) row-wise -----------------------------------------
recompute_K_ET <- function(d) {
  stopifnot(all(c("VPD","Ta") %in% names(d)))
  vVPD <- pmax(d$VPD, eps); vTa <- d$Ta
  (rhoAir * CpAir / gamma) * vVPD * (1/1000) *
    1 / ((2.501e6 - 2361 * vTa) / 1e6) * (3600 * 24) / 1e6 * 365.25
}

# Core physics feature = K_ET * 1/sqrt(VPD)
df_opt <- df_opt %>%
  mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
         K_ET = recompute_K_ET(cur_data_all()),
         phys_core = K_ET * inv_sqrt_VPD)
df <- df %>%
  mutate(inv_sqrt_VPD = 1 / sqrt(pmax(VPD, eps)),
         K_ET = recompute_K_ET(cur_data_all()),
         phys_core = K_ET * inv_sqrt_VPD)

# Optional interaction features among candidate predictors
if (isTRUE(use_interactions)) {
  cand_in_opt <- intersect(rf_candidate_vars, names(df_opt))
  cand_in_opt <- cand_in_opt[vapply(df_opt[cand_in_opt], is.numeric, logical(1))]
  if (length(cand_in_opt) >= 2) {
    pairs <- combn(cand_in_opt, 2, simplify = FALSE)
    for (pr in pairs) {
      nm <- paste0(pr[1], "_x_", pr[2])
      if (!nm %in% names(df_opt)) df_opt[[nm]] <- df_opt[[pr[1]]] * df_opt[[pr[2]]]
      if (!nm %in% names(df))     df[[nm]]     <- df[[pr[1]]]     * df[[pr[2]]]
      rf_candidate_vars <- unique(c(rf_candidate_vars, nm))
    }
  }
}

# Physics core we ALWAYS keep in the ET model
rf_core <- "phys_core"

# ============================================================
# 2) Train RF to predict ET directly (single-stage)
# ============================================================
rf_vars_present_user <- intersect(rf_candidate_vars, names(df_opt))
rf_vars_all <- unique(c(rf_core, rf_vars_present_user))
if (!length(rf_vars_all)) stop("No predictors found.")
if (!all(rf_vars_all %in% names(df))) {
  message("Note: some predictors not present in df; prediction uses available columns only.")
}

df_rf <- df_opt %>%
  select(ET, all_of(rf_vars_all)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

p <- length(setdiff(names(df_rf), "ET"))
if (is.null(rf_mtry)) rf_mtry <- max(1L, floor(sqrt(p)))

set.seed(if (exists("rng_seed")) rng_seed else 123)
rf_fit <- ranger(
  formula = ET ~ .,
  data = df_rf,
  num.trees = rf_num_trees,
  mtry = min(rf_mtry, p),
  min.node.size = rf_min_node,
  sample.fraction = rf_sample_frac,
  replace = TRUE,
  splitrule = rf_splitrule,
  max.depth = rf_max_depth,
  importance = "impurity",
  seed = if (exists("rng_seed")) rng_seed else 123,
  oob.error = TRUE
)

base_oob_rmse <- sqrt(rf_fit$prediction.error)
rf_vars_initial <- setdiff(names(df_rf), "ET")  # used in the initial fit
message(sprintf("Initial RF(ET) OOB RMSE = %.3f", base_oob_rmse))

# ---- VIP (initial, shows ALL candidates incl. CO2_term & phys_core) ----------
if (!is.null(rf_fit$variable.importance)) {
  vi_init <- tibble(
    feature    = names(rf_fit$variable.importance),
    importance = as.numeric(rf_fit$variable.importance)
  ) %>% arrange(desc(importance))
  print(
    ggplot(vi_init, aes(x = reorder(feature, importance), y = importance)) +
      geom_col() + coord_flip() + theme_bw() +
      labs(x = NULL, y = "RF importance (impurity)",
           title = "ET RF: Variable Importance (initial, pre-gate)")
  )
}

# ============================================================
# 2a) (Optional) drop-column CV guard (off by default)
# ============================================================
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
        rf_w <- ranger(ET ~ ., data = df_tr_with, num.trees = rf_num_trees,
                       mtry = min(rf_mtry, ncol(df_tr_with) - 1L),
                       min.node.size = rf_min_node, seed = seed + r*100 + k)
        rf_wo <- ranger(ET ~ ., data = df_tr_wo, num.trees = rf_num_trees,
                        mtry = min(rf_mtry, ncol(df_tr_wo) - 1L),
                        min.node.size = rf_min_node, seed = seed + r*100 + k + 1)
        y_true   <- df_te$ET
        y_pred_w <- as.numeric(predict(rf_w,  data = df_te[, setdiff(names(df_te), "ET"), drop = FALSE])$predictions)
        y_pred_wo<- as.numeric(predict(rf_wo, data = df_te_wo[, setdiff(names(df_te_wo), "ET"), drop = FALSE])$predictions)
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

# ============================================================
# 2b) Permutation importance (pre-gate) on ET
# ============================================================
perm_importance_all <- function(df_rf, base_oob_rmse, B = 200, seed = 1000, target_col = "ET") {
  preds <- setdiff(names(df_rf), target_col)
  if (!length(preds)) return(tibble(feature=character(), mean_delta=numeric(),
                                    lo=numeric(), hi=numeric(), rel=numeric(), B=integer()))
  purrr::map_dfr(preds, function(v) {
    set.seed(seed)
    deltas <- numeric(B)
    for (b in seq_len(B)) {
      dfp <- df_rf
      dfp[[v]] <- sample(dfp[[v]])
      rf_p <- ranger(
        formula = as.formula(paste(target_col, "~ .")),
        data = dfp,
        num.trees = rf_num_trees,
        mtry = min(rf_mtry, length(preds)),
        min.node.size = rf_min_node,
        sample.fraction = rf_sample_frac,
        replace = TRUE,
        splitrule = rf_splitrule,
        max.depth = rf_max_depth,
        seed = seed + b,
        oob.error = TRUE
      )
      deltas[b] <- sqrt(rf_p$prediction.error) - base_oob_rmse
    }
    md <- mean(deltas, na.rm = TRUE)
    tibble(
      feature    = v,
      mean_delta = md,
      lo         = quantile(deltas, perm_ci_levels[1], na.rm = TRUE),
      hi         = quantile(deltas, perm_ci_levels[2], na.rm = TRUE),
      rel        = 100 * md / base_oob_rmse,
      B          = B
    )
  }) %>% arrange(desc(mean_delta))
}

perm_all <- perm_importance_all(df_rf, base_oob_rmse, B = perm_B, target_col = "ET")
cat("\n=== Permutation importance (pre-gate, ET target) ===\n")
print(perm_all %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

print(
  ggplot(perm_all, aes(x = reorder(feature, mean_delta), y = mean_delta)) +
    geom_col() +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    coord_flip() + theme_bw() +
    labs(x = NULL, y = expression(Delta*" OOB RMSE (permute)"),
         title = "Permutation importance (ET, pre-gate)",
         subtitle = sprintf("Keep rule: lo>0 & rel≥%.1f%%", min_rel_keep_pct))
)

# ============================================================
# 2c) Hybrid gate (perm + impurity) for NON-core vars
# ============================================================
# Perm gate
if (tolower(var_gate_when) == "insignificant") {
  keep_perm <- perm_all %>% filter(lo > 0, rel >= min_rel_keep_pct) %>% pull(feature)
} else {
  drop_perm <- perm_all %>% filter(lo > 0, rel >= min_rel_keep_pct) %>% pull(feature)
  keep_perm <- setdiff(perm_all$feature, drop_perm)
}

# Impurity gate (from initial VIP)
keep_imp <- character(0)
if (exists("vi_init")) {
  thr <- imp_keep_frac_max * max(vi_init$importance, na.rm = TRUE)
  keep_imp <- vi_init %>% filter(importance >= thr) %>% arrange(desc(importance)) %>% pull(feature)
  if (is.finite(imp_keep_top_k)) keep_imp <- head(keep_imp, imp_keep_top_k)
}

# Final keep set = physics core + (perm OR impurity) + protected/guard
vars_keep <- union(keep_perm, keep_imp)
vars_keep <- union(vars_keep, intersect(protect_vars, rf_vars_initial))
vars_keep <- union(vars_keep, guard_keep_vars)
vars_keep <- unique(c(rf_core, intersect(vars_keep, rf_vars_initial)))  # ALWAYS keep physics core

# ============================================================
# 2d) CO2 physics SIGN gate using ALE on g_eff = ET_hat / K_ET
# ============================================================
predict_ET_geff_local <- function(newdata, rf_fit_current, rf_vars_current) {
  X <- newdata[, intersect(rf_vars_current, colnames(newdata)), drop = FALSE]
  ET_hat <- as.numeric(predict(rf_fit_current, data = X)$predictions)
  g_hat  <- ET_hat / newdata$K_ET
  list(ET = ET_hat, geff = g_hat)
}

if (require_negative_CO2_effect && "CO2_term" %in% vars_keep) {
  co2_ale_check <- function(var = "CO2_term", data = df_opt, n_bins = 50L) {
    z_all <- data[[var]]; z <- z_all[is.finite(z_all)]
    if (!length(z)) return(list(ok = FALSE, msg = "CO2 has no finite values"))
    edges <- as.numeric(quantile(z, probs = seq(0,1,length.out = n_bins + 1), na.rm = TRUE))
    edges[1] <- min(z, na.rm = TRUE); edges[length(edges)] <- max(z, na.rm = TRUE)
    dG <- numeric(n_bins); cnt <- integer(n_bins)
    for (j in seq_len(n_bins)) {
      lo <- edges[j]; hi <- edges[j + 1]
      idx <- which(z_all >= lo & (z_all < hi | (j == n_bins & z_all <= hi)))
      cnt[j] <- length(idx); if (cnt[j] == 0) { dG[j] <- 0; next }
      df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
      df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
      pr_lo <- predict_ET_geff_local(df_lo, rf_fit, rf_vars_initial)
      pr_hi <- predict_ET_geff_local(df_hi, rf_fit, rf_vars_initial)
      dG[j] <- mean(pr_hi$geff - pr_lo$geff, na.rm = TRUE)
    }
    ale_G <- cumsum(dG); N <- sum(cnt)
    if (N > 0) ale_G <- ale_G - sum(ale_G * cnt, na.rm = TRUE) / N
    total_change <- tail(ale_G, 1) - ale_G[1]
    frac_neg     <- mean(dG < 0, na.rm = TRUE)
    list(ok=TRUE, total_change=total_change, frac_neg=frac_neg)
  }
  chk <- co2_ale_check(n_bins = co2_ale_bins)
  if (isTRUE(chk$ok)) {
    message(sprintf("CO2 sign check (ALE on g_eff): total Δ=%.4g, frac(neg)=%.1f%%",
                    chk$total_change, 100*chk$frac_neg))
    if (!(chk$total_change <= co2_neg_total_tol && chk$frac_neg >= co2_neg_majority_frac)) {
      message("Decision: CO2_term does NOT reduce g_eff → dropping CO2_term.")
      vars_keep <- setdiff(vars_keep, "CO2_term")
    } else {
      message("Decision: CO2_term reduces g_eff → keeping CO2_term.")
    }
  } else {
    message("CO2 sign check skipped: ", chk$msg)
  }
}

# ============================================================
# 2e) Final RF (post gate)
# ============================================================
rf_vars_present <- vars_keep
df_rf <- df_opt %>%
  select(ET, all_of(rf_vars_present)) %>%
  filter(if_all(where(is.numeric), ~ is.finite(.x)))

set.seed(if (exists("rng_seed")) rng_seed else 123)
rf_fit <- ranger(
  formula = ET ~ .,
  data = df_rf,
  num.trees = rf_num_trees,
  mtry = min(rf_mtry, max(1L, ncol(df_rf) - 1L)),
  min.node.size = rf_min_node,
  sample.fraction = rf_sample_frac,
  replace = TRUE,
  splitrule = rf_splitrule,
  max.depth = rf_max_depth,
  importance = "impurity",
  seed = if (exists("rng_seed")) rng_seed else 123,
  oob.error = TRUE
)
base_oob_rmse <- sqrt(rf_fit$prediction.error)
message(sprintf("Final RF(ET) predictors: %s", paste(rf_vars_present, collapse = ", ")))
message(sprintf("Final RF(ET) OOB RMSE = %.3f", base_oob_rmse))

# VIP (final)
if (!is.null(rf_fit$variable.importance)) {
  vi_final <- tibble(
    feature    = names(rf_fit$variable.importance),
    importance = as.numeric(rf_fit$variable.importance)
  ) %>% arrange(desc(importance))
  print(
    ggplot(vi_final, aes(x = reorder(feature, importance), y = importance)) +
      geom_col() + coord_flip() + theme_bw() +
      labs(x = NULL, y = "RF importance (impurity)",
           title = "ET RF: Variable Importance (final, post-gate)")
  )
}

# (Optional) permutation importance for FINAL RF
perm_final <- perm_importance_all(df_rf, base_oob_rmse, B = perm_B, target_col = "ET")
cat("\n=== Permutation importance (FINAL RF, ET target) ===\n")
print(perm_final %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# ============================================================
# 3) Predictions & metrics on df
# ============================================================
predict_ET_geff <- function(newdata) {
  X <- newdata[, intersect(rf_vars_present, colnames(newdata)), drop = FALSE]
  ET_hat <- as.numeric(predict(rf_fit, data = X)$predictions)
  g_hat  <- ET_hat / newdata$K_ET
  list(ET = ET_hat, geff = g_hat)
}

pred <- predict_ET_geff(df)
output_df <- df %>%
  mutate(ET_hat = pred$ET,
         g_eff_hat = pred$geff)

RMSE_ET <- sqrt(mean((output_df$ET - output_df$ET_hat)^2, na.rm = TRUE))
R2_ET   <- suppressWarnings(cor(output_df$ET, output_df$ET_hat, use = "complete.obs")^2)
RMSE_geff <- sqrt(mean((output_df$g_eff - output_df$g_eff_hat)^2, na.rm = TRUE))
R2_geff   <- suppressWarnings(cor(output_df$g_eff, output_df$g_eff_hat, use = "complete.obs")^2)

message(sprintf("ET   : RMSE = %.3f, R^2 = %.3f", RMSE_ET,   R2_ET))
message(sprintf("g_eff: RMSE = %.3f, R^2 = %.3f", RMSE_geff, R2_geff))

# Scatter plots
p_et <- output_df %>%
  filter(is.finite(ET), is.finite(ET_hat)) %>%
  ggplot(aes(ET, ET_hat)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(x = "Observed ET", y = "Predicted ET",
       title = "Observed vs Predicted ET",
       subtitle = sprintf("RMSE=%.3f, R^2=%.3f (OOB RMSE=%.3f)", RMSE_ET, R2_ET, base_oob_rmse))

p_geff <- output_df %>%
  filter(is.finite(g_eff), is.finite(g_eff_hat)) %>%
  ggplot(aes(g_eff, g_eff_hat)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_bw() +
  labs(x = expression(paste("Observed ", g[eff])),
       y = expression(paste("Predicted ", g[eff])),
       title = expression(paste("Observed vs Predicted ", g[eff])),
       subtitle = sprintf("RMSE=%.3f, R^2=%.3f", RMSE_geff, R2_geff))

print(p_et + p_geff)

# ============================================================
# 4) ALE (Accumulated Local Effects), FINAL RF(ET)
#     Show ONLY variables actually used (including physics core).
# ============================================================
ale_one <- function(var, data = df, n_bins = 50L) {
  stopifnot(var %in% names(data))
  z_all <- data[[var]]; z <- z_all[is.finite(z_all)]
  if (!length(z)) stop("Variable ", var, " has no finite values.")
  edges <- as.numeric(quantile(z, probs = seq(0,1,length.out = n_bins + 1), na.rm = TRUE))
  edges[1] <- min(z, na.rm = TRUE); edges[length(edges)] <- max(z, na.rm = TRUE)
  dET <- numeric(n_bins); dG  <- numeric(n_bins); cnt <- integer(n_bins)
  for (j in seq_len(n_bins)) {
    lo <- edges[j]; hi <- edges[j + 1]
    idx <- which(z_all >= lo & (z_all < hi | (j == n_bins & z_all <= hi)))
    cnt[j] <- length(idx); if (cnt[j] == 0) { dET[j] <- 0; dG[j] <- 0; next }
    df_lo <- data[idx, , drop = FALSE]; df_lo[[var]] <- lo
    df_hi <- data[idx, , drop = FALSE]; df_hi[[var]] <- hi
    pr_lo <- predict_ET_geff(df_lo); pr_hi <- predict_ET_geff(df_hi)
    dET[j] <- mean(pr_hi$ET   - pr_lo$ET,   na.rm = TRUE)
    dG[j]  <- mean(pr_hi$geff - pr_lo$geff, na.rm = TRUE)
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

vars_for_ALE <- intersect(rf_vars_present, names(df))   # only used vars (includes phys_core)
stopifnot(length(vars_for_ALE) >= 1)

n_bins <- 50L
ale_raw <- bind_rows(lapply(vars_for_ALE, ale_one, data = df, n_bins = n_bins))

# Anchor to absolute scales with model means
mu_ET_pred   <- mean(output_df$ET_hat,   na.rm = TRUE)
mu_geff_pred <- mean(output_df$g_eff_hat, na.rm = TRUE)

# Smooth & anchor
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
    ALE_ET_abs   = (y_ET_grid - off_ET) + mu_ET_pred,
    ALE_geff_abs = (y_G_grid  - off_G)  + mu_geff_pred
  )
}
ale_smoothed <- purrr::map_dfr(vars_for_ALE, smooth_one_var)

# Nice labels
lab_map <- setNames(vars_for_ALE, vars_for_ALE)
if ("phys_core" %in% vars_for_ALE) lab_map["phys_core"] <- "K_ET × 1/√VPD"
if ("CO2_term" %in% vars_for_ALE) lab_map["CO2_term"] <- "CO2_term"

plot_ale_et_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_ET_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(lab_map)) +
  labs(x = NULL, y = "ET (mm/yr) — ALE, anchored to model mean",
       title = "ALE (ET) — final RF(ET)") +
  theme_bw()

plot_ale_geff_abs <- ggplot(ale_smoothed, aes(x = x, y = ALE_geff_abs)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var, scales = "free_x", labeller = as_labeller(lab_map)) +
  labs(x = NULL, y = expression(g[eff]~"(mm s"^{-1}*") — ALE, anchored to model mean"),
       title = expression("ALE ("*g[eff]*") — derived from ET_hat/K_ET")) +
  theme_bw()

print(plot_ale_et_abs); print(plot_ale_geff_abs)

# ============================================================
# 5) Save bundle (optional)
# ============================================================
geff_et_rf_bundle <- list(
  metadata      = if (exists("metadata")) paste0(metadata, " | Single-stage RF(ET) with phys_core = K_ET*inv_sqrt_VPD.") else
    "Single-stage RF(ET) with phys_core = K_ET*inv_sqrt_VPD",
  rf_fit        = rf_fit,
  rf_vars       = rf_vars_present,
  out_df        = output_df,
  RMSE_ET       = RMSE_ET,   R2_ET = R2_ET,
  RMSE_geff     = RMSE_geff, R2_geff = R2_geff,
  rng_seed      = if (exists("rng_seed")) rng_seed else NA_integer_,
  vip_init      = if (exists("vi_init")) vi_init else NULL,
  perm_all      = perm_all,
  perm_final    = perm_final,
  ale_raw       = ale_raw,
  ale_smooth    = ale_smoothed,
  gates         = list(mode = var_gate_when, rel_pct = min_rel_keep_pct,
                       impurity = list(enabled = use_impurity_gate, frac_max = imp_keep_frac_max,
                                       top_k = imp_keep_top_k),
                       co2_sign = list(required = require_negative_CO2_effect,
                                       bins = co2_ale_bins,
                                       frac_req = co2_neg_majority_frac,
                                       tol = co2_neg_total_tol))
)
if (exists("out_file")) {
  save(geff_et_rf_bundle, file = sub("\\.RData$", "_RF_ET_singleStage_physCORE.RData", out_file))
}
