# ============================================================
# Var_Model_MI_Rubin.R
#
# Purpose:
#   1. Fit single-marker and combined logistic models in all
#      completed MI datasets and pool coefficients using Rubin's rules.
#   2. Assess model performance with patient-level bootstrap
#      optimism correction across the 20 imputations.
#   3. Generate PCA/t-SNE exploratory plots, pooled OR forest plots,
#      calibration, ROC, DCA, confusion matrices, and a pooled-coefficient
#      nomogram.
#
# Input:
#   Z_Norm_Markers.rds       (all 20 completed datasets, stacked)
#   Screening_Results.rds    (selected marker names from Var_Screen.R)
#
# Important:
#   - Screening was previously performed on imputation_id = 1.
#   - Existing z_Age and z_marker columns are used as-is.
#   - Sex and CCI are taken from their original columns because the
#     current Z_Norm_Markers.rds does not contain z_Sex or z_CCI.
#   - No additional standardization is performed.
#   - Bootstrap resampling uses the same patient indices in all imputations.
#   - This is an MI-before-bootstrap (MI-prior-IMV) practical strategy using
#     the already completed datasets. It does not repeat the original imputation
#     model or the LASSO selection inside each bootstrap replicate.
# ============================================================

rm(list = ls())
gc()
graphics.off()

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

# ------------------------------
# 1. Parameters and input files
# ------------------------------

base_dir <- "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS"
z_file <- file.path(base_dir, "Z_Norm_Markers.rds")
screening_file <- file.path(base_dir, "Screening_Results.rds")
output_dir <- file.path(base_dir, "Final_Model_MI_Rubin")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

B <- 1000
random_seed <- 20261112
classification_threshold <- 0.50
dca_thresholds <- seq(0.01, 0.99, by = 0.01)

# Use the current Z_Norm_Markers structure:
# Age is already z-scored; Sex and CCI remain in their original coding.
covariates <- c("z_Age", "Sex", "CCI")
id_var <- "No"
group_var <- "Tag_No"
imp_var <- "imputation_id"

# Existing binary comparisons from Var_Screen.R
comparison_groups <- list(
  "0_vs_1plus2plus3" = c("0", "1", "2", "3"),
  "0_vs_3plus4" = c("0", "3", "4"),
  "0_vs_1plus2plus3plus4" = c("0", "1", "2", "3", "4")
)

# NULL: retrieve the LASSO-selected markers from Screening_Results.rds.
# To force the same four markers in all comparisons, uncomment:
# final_markers_override <- c("z_tau217", "z_Ab40", "z_tau181", "z_miR26b")
final_markers_override <- NULL

set.seed(random_seed)

save_checkpoint <- function(object, filename) {
  saveRDS(object, file.path(output_dir, filename))
  cat("RDS saved:", filename, "\n")
}

write_csv <- function(object, filename) {
  utils::write.csv(
    object,
    file.path(output_dir, filename),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# ------------------------------
# 2. Package checks and data input
# ------------------------------

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("缺少 mice 程序包。请先运行 install.packages('mice')。")
}
if (!requireNamespace("Rtsne", quietly = TRUE)) {
  message("未检测到 Rtsne：PCA 将正常运行，t-SNE 将跳过。需要 t-SNE 时请安装 install.packages('Rtsne')。")
}

if (!file.exists(z_file)) stop("未找到：", z_file)
if (!file.exists(screening_file)) stop("未找到：", screening_file)

mi_long <- readRDS(z_file)
screening_objects <- readRDS(screening_file)

required_screening_objects <- c("binary_results", "markers")
missing_screening_objects <- setdiff(required_screening_objects, names(screening_objects))
if (length(missing_screening_objects) > 0) {
  stop("Screening_Results.rds 缺少对象：",
       paste(missing_screening_objects, collapse = ", "))
}

binary_results <- screening_objects$binary_results
screening_covariates <- screening_objects$covariates
if (!is.null(screening_covariates) && !identical(screening_covariates, covariates)) {
  message(
    "注意：Screening_Results.rds 记录的协变量为 [",
    paste(screening_covariates, collapse = ", "),
    "]；本脚本按当前 Z_Norm_Markers.rds 使用 [",
    paste(covariates, collapse = ", "),
    "]。筛选结果仅用于读取标志物名称，不直接读取旧模型数据或旧协变量系数。"
  )
}

required_data_cols <- c(id_var, group_var, imp_var, covariates)
missing_data_cols <- setdiff(required_data_cols, names(mi_long))
if (length(missing_data_cols) > 0) {
  stop("Z_Norm_Markers.rds 缺少必要列：",
       paste(missing_data_cols, collapse = ", "))
}

imputation_ids <- sort(unique(mi_long[[imp_var]]))
M <- length(imputation_ids)
if (M != 20L) {
  warning("当前检测到 ", M, " 个插补数据集，而预期为 20 个。脚本将按实际数量处理。")
}

# Check that each imputation contains one row per patient and identical IDs/groups.
mi_data_list <- lapply(imputation_ids, function(k) {
  d <- mi_long[mi_long[[imp_var]] == k, , drop = FALSE]
  d <- d[order(as.character(d[[id_var]])), , drop = FALSE]
  rownames(d) <- NULL
  d
})
names(mi_data_list) <- paste0("Imp", imputation_ids)

ref_ids <- as.character(mi_data_list[[1]][[id_var]])
ref_groups <- as.character(mi_data_list[[1]][[group_var]])
if (anyDuplicated(ref_ids)) stop("参考插补集存在重复病例 ID。")

for (i in seq_along(mi_data_list)) {
  d <- mi_data_list[[i]]
  if (anyDuplicated(as.character(d[[id_var]]))) {
    stop("imputation_id = ", imputation_ids[i], " 存在重复病例 ID。")
  }
  if (!identical(as.character(d[[id_var]]), ref_ids)) {
    stop("不同插补集的病例 ID 不一致：imputation_id = ", imputation_ids[i])
  }
  if (!identical(as.character(d[[group_var]]), ref_groups)) {
    stop("不同插补集的 Tag_No 不一致：imputation_id = ", imputation_ids[i])
  }
}

n_all <- nrow(mi_data_list[[1]])
cat("\n插补数据检查通过：M = ", M, "；每个插补集 N = ", n_all, "\n", sep = "")

# Validate covariate types. Sex is retained in its original binary numeric coding.
for (i in seq_along(mi_data_list)) {
  d <- mi_data_list[[i]]
  if (!is.numeric(d$z_Age)) stop("z_Age 必须为数值型。")
  if (!is.numeric(d$Sex)) stop("当前模型要求 Sex 为原始数值型二分类编码；请检查 MI_data.rds。")
  if (length(unique(d$Sex)) != 2) stop("Sex 在插补集 ", imputation_ids[i], " 中不是两个取值。")
  if (!is.numeric(d$CCI)) stop("CCI 必须为数值型原始评分。")
}

# ------------------------------
# 3. Helper functions
# ------------------------------

get_final_markers <- function(comparison_name) {
  if (!is.null(final_markers_override)) {
    selected <- final_markers_override
  } else {
    br <- binary_results[[comparison_name]]
    if (is.null(br)) stop("Screening_Results.rds 中缺少比较：", comparison_name)
    selected <- br$selected_markers

    # Preserve the fallback logic used in the previous Var_Model.R.
    if (length(selected) == 0) {
      warning(comparison_name, " 的 LASSO 未选出标志物，尝试使用单变量 P 值最小的前 3 个标志物。")
      uni <- br$univariate
      if (!is.null(uni) && nrow(uni) > 0 && all(c("Marker", "P_value") %in% names(uni))) {
        uni <- uni[is.finite(uni$P_value), , drop = FALSE]
        uni <- uni[order(uni$P_value), , drop = FALSE]
        selected <- head(uni$Marker, 3)
      }
      if (length(selected) == 0 && !is.null(br$candidates)) selected <- head(br$candidates, 3)
      if (length(selected) == 0) stop("无法为比较 ", comparison_name, " 确定最终标志物。")
    }
  }

  selected <- unique(as.character(selected))
  selected <- selected[!is.na(selected) & nzchar(selected)]
  selected <- setdiff(selected, covariates)
  if (length(selected) == 0) stop("比较 ", comparison_name, " 没有有效标志物。")
  selected
}

make_comparison_data <- function(data_list, selected_groups) {
  lapply(data_list, function(d) {
    group_chr <- as.character(d[[group_var]])
    keep <- !is.na(group_chr) & group_chr %in% selected_groups
    out <- d[keep, , drop = FALSE]
    out$Y_binary <- as.integer(as.character(out[[group_var]]) != "0")
    out$Y_binary <- as.integer(out$Y_binary)
    out
  })
}

check_model_data <- function(data_list, predictors, comparison_name) {
  required <- c("Y_binary", predictors)
  for (i in seq_along(data_list)) {
    d <- data_list[[i]]
    missing <- setdiff(required, names(d))
    if (length(missing) > 0) {
      stop(comparison_name, " / 插补集 ", imputation_ids[i], " 缺少变量：",
           paste(missing, collapse = ", "))
    }
    if (any(!stats::complete.cases(d[, required, drop = FALSE]))) {
      stop(comparison_name, " / 插补集 ", imputation_ids[i], " 的模型变量仍有缺失值。")
    }
    for (v in predictors) {
      if (is.numeric(d[[v]]) && any(!is.finite(d[[v]]))) {
        stop(comparison_name, " / 插补集 ", imputation_ids[i], " 的变量 ", v, " 存在非有限值。")
      }
      if (length(unique(d[[v]])) < 2) {
        stop(comparison_name, " / 插补集 ", imputation_ids[i], " 的变量 ", v, " 为零方差。")
      }
    }
    if (length(unique(d$Y_binary)) != 2) {
      stop(comparison_name, " / 插补集 ", imputation_ids[i], " 未同时包含 0 和 1。")
    }
  }
  invisible(TRUE)
}

fit_logistic <- function(data, predictors) {
  fit <- suppressWarnings(
    stats::glm(
      stats::reformulate(predictors, response = "Y_binary"),
      data = data,
      family = stats::binomial(link = "logit"),
      control = stats::glm.control(maxit = 100)
    )
  )
  if (!isTRUE(fit$converged)) stop("Logistic 模型未收敛。")
  if (anyNA(stats::coef(fit)) || any(!is.finite(stats::coef(fit)))) {
    stop("Logistic 模型系数存在 NA 或非有限值，可能存在完全分离或共线性。")
  }
  fit
}

pool_logistic <- function(data_list, predictors) {
  fits <- lapply(data_list, fit_logistic, predictors = predictors)
  mira <- mice::as.mira(fits)
  pooled <- mice::pool(mira)
  tab <- as.data.frame(summary(pooled, conf.int = TRUE, exponentiate = FALSE))

  # Explicit, stable column names for downstream plotting/export.
  names(tab)[names(tab) == "2.5 %"] <- "CI_Lower_Beta"
  names(tab)[names(tab) == "97.5 %"] <- "CI_Upper_Beta"
  tab$OR <- exp(tab$estimate)
  tab$OR_Lower95 <- exp(tab$CI_Lower_Beta)
  tab$OR_Upper95 <- exp(tab$CI_Upper_Beta)
  list(fits = fits, pooled = pooled, table = tab)
}

calculate_auc <- function(y, probability) {
  tryCatch({
    r <- pROC::roc(y, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
    as.numeric(pROC::auc(r))
  }, error = function(e) NA_real_)
}

safe_logit <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

calculate_calibration <- function(y, probability) {
  lp <- safe_logit(probability)
  ans <- tryCatch({
    fit <- suppressWarnings(stats::glm(y ~ lp, family = stats::binomial()))
    b <- stats::coef(fit)
    if (length(b) < 2 || any(!is.finite(b))) stop("无效校准系数")
    c(intercept = unname(b[1]), slope = unname(b[2]))
  }, error = function(e) c(intercept = NA_real_, slope = NA_real_))
  ans
}

calculate_net_benefit <- function(y, probability, thresholds) {
  vapply(thresholds, function(pt) {
    pred <- probability >= pt
    tp <- sum(pred & y == 1)
    fp <- sum(pred & y == 0)
    tp / length(y) - fp / length(y) * pt / (1 - pt)
  }, numeric(1))
}

calculate_classification <- function(y, probability, threshold) {
  pred <- as.integer(probability >= threshold)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  tp <- sum(pred == 1 & y == 1)
  div <- function(a, b) if (b == 0) NA_real_ else a / b
  metrics <- data.frame(
    Threshold = threshold, TN = tn, FP = fp, FN = fn, TP = tp,
    Accuracy = div(tn + tp, length(y)),
    Sensitivity = div(tp, tp + fn),
    Specificity = div(tn, tn + fp),
    PPV = div(tp, tp + fp),
    NPV = div(tn, tn + fn),
    F1 = div(2 * tp, 2 * tp + fp + fn)
  )
  confusion <- data.frame(
    Actual = c("Actual_0", "Actual_1"),
    Predicted_0 = c(tn, fn),
    Predicted_1 = c(fp, tp)
  )
  list(metrics = metrics, confusion = confusion)
}

make_calibration_curve <- function(y, probability, comparison_name, model_id, model_label, groups = 5) {
  n <- length(y)
  group_id <- pmin(groups, ceiling(rank(probability, ties.method = "first") / n * groups))
  idx_list <- split(seq_len(n), group_id)
  do.call(rbind, lapply(seq_along(idx_list), function(i) {
    idx <- idx_list[[i]]
    events <- sum(y[idx] == 1)
    ci <- stats::binom.test(events, length(idx))$conf.int
    data.frame(
      Comparison = comparison_name, Model = model_id, Model_Label = model_label,
      Group = i, N = length(idx), Mean_Predicted = mean(probability[idx]),
      Observed_Rate = events / length(idx),
      Observed_Lower95 = ci[1], Observed_Upper95 = ci[2]
    )
  }))
}

# ------------------------------
# 4. Bootstrap validation for MI-pooled predictions
# ------------------------------

evaluate_mi_model <- function(data_list, predictors, comparison_name,
                              model_id, model_label, B = 1000) {
  cat("\n模型：", model_label, "\n变量：", paste(predictors, collapse = ", "), "\n")
  pooled_obj <- pool_logistic(data_list, predictors)
  fits <- pooled_obj$fits
  pooled_table <- pooled_obj$table

  y <- data_list[[1]]$Y_binary
  n <- length(y)
  m <- length(data_list)

  # Average individual predicted probabilities across completed datasets.
  p_full_by_imp <- lapply(seq_len(m), function(i) {
    as.numeric(stats::predict(fits[[i]], newdata = data_list[[i]], type = "response"))
  })
  p_full <- rowMeans(do.call(cbind, p_full_by_imp))
  auc_apparent <- calculate_auc(y, p_full)
  cal_apparent <- calculate_calibration(y, p_full)
  if (!is.finite(auc_apparent)) stop("无法计算表观 AUC：", model_id)
  if (any(!is.finite(cal_apparent))) stop("无法计算表观校准指标：", model_id)

  # Precompute model matrices for faster bootstrap fitting.
  form <- stats::reformulate(predictors, response = "Y_binary")
  x_list <- lapply(data_list, function(d) stats::model.matrix(form, data = d))

  auc_optimism <- numeric(0)
  auc_corrected_replicates <- numeric(0)
  cal_int_optimism <- numeric(0)
  cal_slope_optimism <- numeric(0)
  bootstrap_test_predictions <- list()
  bootstrap_train_predictions <- list()
  bootstrap_train_y <- list()
  valid_boot <- 0L

  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    y_boot <- y[idx]
    if (length(unique(y_boot)) < 2) next

    p_train_by_imp <- vector("list", m)
    p_test_by_imp <- vector("list", m)
    successful <- TRUE

    for (j in seq_len(m)) {
      xj <- x_list[[j]]
      fit_b <- tryCatch(
        suppressWarnings(stats::glm.fit(
          x = xj[idx, , drop = FALSE],
          y = y_boot,
          family = stats::binomial(),
          control = stats::glm.control(maxit = 100)
        )),
        error = function(e) NULL
      )
      if (is.null(fit_b) || !isTRUE(fit_b$converged) ||
          anyNA(fit_b$coefficients) || any(!is.finite(fit_b$coefficients))) {
        successful <- FALSE
        break
      }
      beta_b <- fit_b$coefficients
      p_train_by_imp[[j]] <- stats::plogis(as.numeric(xj[idx, , drop = FALSE] %*% beta_b))
      p_test_by_imp[[j]] <- stats::plogis(as.numeric(xj %*% beta_b))
    }
    if (!successful) next

    p_train <- rowMeans(do.call(cbind, p_train_by_imp))
    p_test <- rowMeans(do.call(cbind, p_test_by_imp))
    auc_train <- calculate_auc(y_boot, p_train)
    auc_test <- calculate_auc(y, p_test)
    if (!is.finite(auc_train) || !is.finite(auc_test)) next

    auc_opt <- auc_train - auc_test
    auc_optimism <- c(auc_optimism, auc_opt)
    auc_corrected_replicates <- c(auc_corrected_replicates, auc_apparent - auc_opt)

    cal_train <- calculate_calibration(y_boot, p_train)
    cal_test <- calculate_calibration(y, p_test)
    if (all(is.finite(cal_train)) && all(is.finite(cal_test))) {
      cal_int_optimism <- c(cal_int_optimism, cal_train["intercept"] - cal_test["intercept"])
      cal_slope_optimism <- c(cal_slope_optimism, cal_train["slope"] - cal_test["slope"])
    }

    bootstrap_test_predictions[[length(bootstrap_test_predictions) + 1L]] <- p_test
    bootstrap_train_predictions[[length(bootstrap_train_predictions) + 1L]] <- p_train
    bootstrap_train_y[[length(bootstrap_train_y) + 1L]] <- y_boot
    valid_boot <- valid_boot + 1L
    if (b %% 100 == 0) cat("  Bootstrap：", b, "/", B, "；有效：", valid_boot, "\n")
  }

  if (length(auc_optimism) < 100 || length(cal_int_optimism) < 100 || valid_boot < 100) {
    warning(model_id, " 有效 Bootstrap 次数不足 100；请谨慎解释内部验证结果。")
  }
  if (length(auc_optimism) == 0 || length(cal_int_optimism) == 0) {
    stop("没有足够有效的 Bootstrap 结果：", model_id)
  }

  auc_corrected <- auc_apparent - mean(auc_optimism)
  auc_ci <- as.numeric(stats::quantile(
    auc_corrected_replicates, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE
  ))

  int_corrected <- unname(cal_apparent["intercept"] - mean(cal_int_optimism))
  slope_corrected <- unname(cal_apparent["slope"] - mean(cal_slope_optimism))
  if (!is.finite(int_corrected) || !is.finite(slope_corrected)) {
    stop("校准乐观偏倚校正失败：", model_id)
  }

  p_corrected <- stats::plogis(int_corrected + slope_corrected * safe_logit(p_full))
  class_result <- calculate_classification(y, p_corrected, classification_threshold)
  class_result$metrics$Comparison <- comparison_name
  class_result$metrics$Model <- model_id
  class_result$metrics$Model_Label <- model_label
  class_result$metrics$Estimate_Type <- "Full-sample classification; not optimism-corrected"
  class_result$confusion$Comparison <- comparison_name
  class_result$confusion$Model <- model_id
  class_result$confusion$Model_Label <- model_label

  calibration_curve <- make_calibration_curve(
    y, p_corrected, comparison_name, model_id, model_label,
    groups = min(5, n)
  )

  # DCA optimism correction using paired bootstrap training/test predictions.
  nb_apparent <- calculate_net_benefit(y, p_corrected, dca_thresholds)
  # Training outcomes are retained together with the paired bootstrap predictions.
  dca_optimism_matrix <- do.call(rbind, lapply(seq_along(bootstrap_test_predictions), function(k) {
    p_train_b <- stats::plogis(int_corrected + slope_corrected * safe_logit(bootstrap_train_predictions[[k]]))
    p_test_b <- stats::plogis(int_corrected + slope_corrected * safe_logit(bootstrap_test_predictions[[k]]))
    nb_train_b <- calculate_net_benefit(bootstrap_train_y[[k]], p_train_b, dca_thresholds)
    nb_test_b <- calculate_net_benefit(y, p_test_b, dca_thresholds)
    nb_train_b - nb_test_b
  }))
  mean_dca_optimism <- colMeans(dca_optimism_matrix, na.rm = TRUE)
  nb_corrected <- nb_apparent - mean_dca_optimism
  nb_test_matrix <- do.call(rbind, lapply(bootstrap_test_predictions, function(p) {
    calculate_net_benefit(y, p, dca_thresholds)
  }))
  nb_test_mean <- colMeans(nb_test_matrix, na.rm = TRUE)
  dca_data <- data.frame(
    Comparison = comparison_name, Model = model_id, Model_Label = model_label,
    Threshold = dca_thresholds,
    NB_Apparent_Calibrated = nb_apparent,
    NB_Optimism_Corrected = nb_corrected,
    NB_Bootstrap_Test_Mean = nb_test_mean
  )

  # ROC curve from the MI-averaged full-sample probabilities.
  roc_full <- pROC::roc(y, p_corrected, levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_main <- data.frame(
    FPR = 1 - as.numeric(roc_full$specificities),
    TPR = as.numeric(roc_full$sensitivities)
  )
  roc_main <- roc_main[order(roc_main$FPR), , drop = FALSE]

  fpr_grid <- seq(0, 1, by = 0.01)
  tpr_mat <- matrix(NA_real_, nrow = length(fpr_grid), ncol = length(bootstrap_test_predictions))
  for (k in seq_along(bootstrap_test_predictions)) {
    p_test <- bootstrap_test_predictions[[k]]
    rb <- tryCatch(pROC::roc(y, p_test, levels = c(0, 1), direction = "<", quiet = TRUE),
                   error = function(e) NULL)
    if (is.null(rb)) next
    xx <- 1 - as.numeric(rb$specificities)
    yy <- as.numeric(rb$sensitivities)
    ord <- order(xx)
    tpr_mat[, k] <- tryCatch(
      stats::approx(xx[ord], yy[ord], xout = fpr_grid, ties = mean, rule = 2)$y,
      error = function(e) rep(NA_real_, length(fpr_grid))
    )
  }
  tpr_lower <- apply(tpr_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  tpr_upper <- apply(tpr_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  roc_ci <- list(
    main_curve = roc_main,
    auc_apparent = auc_apparent,
    auc_optimism_corrected = auc_corrected,
    auc_corrected_ci = auc_ci,
    ci_band = data.frame(FPR = fpr_grid, TPR_Lower = tpr_lower, TPR_Upper = tpr_upper)
  )

  prediction_data <- data.frame(
    No = data_list[[1]][[id_var]], Y_binary = y,
    Probability_MI_Averaged = p_full,
    Probability_Calibration_Adjusted = p_corrected,
    LinearPredictor_MI_Averaged = safe_logit(p_full)
  )

  summary_data <- data.frame(
    Comparison = comparison_name, Model = model_id, Model_Label = model_label,
    N = n, Events = sum(y == 1), NonEvents = sum(y == 0),
    Predictors = paste(predictors, collapse = "; "),
    AUC_Apparent = auc_apparent,
    AUC_Optimism_Corrected = auc_corrected,
    AUC_Corrected_Lower95 = auc_ci[1],
    AUC_Corrected_Upper95 = auc_ci[2],
    Calibration_Intercept_Apparent = unname(cal_apparent["intercept"]),
    Calibration_Intercept_Corrected = int_corrected,
    Calibration_Slope_Apparent = unname(cal_apparent["slope"]),
    Calibration_Slope_Corrected = slope_corrected,
    Bootstrap_Valid = valid_boot
  )

  list(
    fits = fits, pooled_model = pooled_obj$pooled, pooled_table = pooled_table,
    summary = summary_data, predictions = prediction_data,
    calibration_curve = calibration_curve,
    classification_metrics = class_result$metrics,
    confusion_matrix = class_result$confusion,
    dca = dca_data, roc_ci = roc_ci,
    auc_corrected_replicates = auc_corrected_replicates
  )
}

# ------------------------------
# 5. PCA and t-SNE for selected markers
# ------------------------------

make_subject_mean_data <- function(data_list, vars) {
  base <- data_list[[1]][, c(id_var, group_var, "Y_binary"), drop = FALSE]
  for (v in vars) {
    mat <- do.call(cbind, lapply(data_list, function(d) as.numeric(d[[v]])))
    base[[v]] <- rowMeans(mat)
  }
  base
}

plot_dimension_reduction <- function(subject_data, selected_markers, comparison_name) {
  X <- as.matrix(subject_data[, selected_markers, drop = FALSE])
  if (any(!is.finite(X))) stop("PCA/t-SNE 输入矩阵存在非有限值：", comparison_name)
  if (ncol(X) < 2) {
    warning(comparison_name, " 的标志物少于 2 个，跳过 PCA/t-SNE。")
    return(NULL)
  }

  # Markers are already z-scored. Center only; do not scale again.
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  pca_var <- pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100
  pca_data <- data.frame(
    PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2],
    Y_binary = factor(subject_data$Y_binary, levels = c(0, 1),
                      labels = c("Control", "AD"))
  )
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Y_binary)) +
    geom_point(size = 2.7, alpha = 0.85) +
    stat_ellipse(aes(group = Y_binary), type = "norm", level = 0.68,
                 linewidth = 0.7, show.legend = FALSE, na.rm = TRUE) +
    labs(
      title = paste0("PCA: ", comparison_name),
      subtitle = paste0("Markers: ", paste(selected_markers, collapse = ", ")),
      x = sprintf("PC1 (%.1f%%)", pca_var[1]),
      y = sprintf("PC2 (%.1f%%)", pca_var[2]), color = NULL
    ) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                           legend.position = "bottom")
  print(p_pca)
  ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_PCA.png")),
                  p_pca, width = 6.5, height = 5.2, dpi = 300)

  tsne_data <- NULL
  p_tsne <- NULL
  if (requireNamespace("Rtsne", quietly = TRUE)) {
    perplexity_use <- min(10, floor((nrow(X) - 1) / 3) - 1)
    if (perplexity_use >= 2) {
      set.seed(random_seed)
      tsne_fit <- Rtsne::Rtsne(
        X, dims = 2, perplexity = perplexity_use,
        theta = 0.5, max_iter = 1000,
        pca = FALSE, check_duplicates = FALSE,
        verbose = FALSE
      )
      tsne_data <- data.frame(
        tSNE1 = tsne_fit$Y[, 1], tSNE2 = tsne_fit$Y[, 2],
        Y_binary = factor(subject_data$Y_binary, levels = c(0, 1),
                          labels = c("Control", "AD"))
      )
      p_tsne <- ggplot(tsne_data, aes(tSNE1, tSNE2, color = Y_binary)) +
        geom_point(size = 2.7, alpha = 0.85) +
        labs(
          title = paste0("t-SNE: ", comparison_name),
          subtitle = paste0("Perplexity = ", perplexity_use,
                            "; markers: ", paste(selected_markers, collapse = ", ")),
          x = "t-SNE 1", y = "t-SNE 2", color = NULL
        ) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                               legend.position = "bottom")
      print(p_tsne)
      ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_tSNE.png")),
                      p_tsne, width = 6.5, height = 5.2, dpi = 300)
    } else {
      warning("样本量不足以设置有效 perplexity，跳过 t-SNE：", comparison_name)
    }
  }

  list(
    input_data = subject_data,
    pca = pca_fit,
    pca_data = pca_data,
    pca_variance_percent = pca_var,
    pca_plot = p_pca,
    tsne_data = tsne_data,
    tsne_plot = p_tsne
  )
}

# ------------------------------
# 6. Pooled OR forest plot and nomogram
# ------------------------------

plot_pooled_forest <- function(pooled_table, comparison_name, model_id) {
  d <- pooled_table[pooled_table$term != "(Intercept)", , drop = FALSE]
  d <- d[is.finite(d$OR) & is.finite(d$OR_Lower95) & is.finite(d$OR_Upper95) &
           d$OR > 0 & d$OR_Lower95 > 0 & d$OR_Upper95 > 0, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$term <- factor(d$term, levels = rev(d$term))
  p <- ggplot(d, aes(x = OR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = OR_Lower95, xend = OR_Upper95,
                     y = term, yend = term), linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_x_log10() +
    labs(title = paste0("Rubin-pooled OR: ", comparison_name, " / ", model_id),
         x = "Odds ratio (95% CI; log scale)", y = NULL) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_", model_id, "_Pooled_OR_Forest.png")),
                  p, width = 7, height = max(4, 0.45 * nrow(d) + 1.5), dpi = 300)
  p
}

make_nomogram <- function(pooled_table, data_ref, predictors, comparison_name) {
  beta <- stats::setNames(pooled_table$estimate, pooled_table$term)
  if (!"(Intercept)" %in% names(beta)) stop("Rubin 合并结果中缺少截距，无法绘制列线图。")
  if (any(!predictors %in% names(beta))) {
    stop("列线图变量未能与合并系数一一对应；请检查是否存在分类变量展开或共线性。")
  }
  intercept <- beta[["(Intercept)"]]
  ranges <- lapply(predictors, function(v) range(data_ref[[v]], na.rm = TRUE))
  names(ranges) <- predictors
  effect_ranges <- vapply(seq_along(predictors), function(i) {
    abs(beta[[predictors[i]]]) * diff(ranges[[i]])
  }, numeric(1))
  max_effect <- max(effect_ranges)
  if (!is.finite(max_effect) || max_effect <= 0) stop("列线图变量贡献范围无效。")

  min_contrib <- vapply(seq_along(predictors), function(i) {
    b <- beta[[predictors[i]]]
    rr <- ranges[[i]]
    min(b * rr)
  }, numeric(1))
  max_points_by_var <- effect_ranges / max_effect * 100
  total_max <- sum(max_points_by_var)
  baseline_logit <- intercept + sum(min_contrib)

  # Axis rows: Points, each predictor, Total Points, Predicted probability.
  row_names <- c("Points", predictors, "Total Points", "Predicted probability")
  y_positions <- rev(seq_along(row_names))
  names(y_positions) <- row_names
  axis_lines <- list()
  axis_ticks <- list()
  axis_text <- list()

  # Common points axis.
  axis_lines[[1]] <- data.frame(x = 0, xend = 100, y = y_positions["Points"], yend = y_positions["Points"])
  pts <- seq(0, 100, by = 20)
  axis_ticks[[1]] <- data.frame(x = pts, xend = pts,
                                 y = y_positions["Points"] - 0.08,
                                 yend = y_positions["Points"] + 0.08)
  axis_text[[1]] <- data.frame(x = pts, y = y_positions["Points"] + 0.22,
                                label = as.character(pts))

  for (i in seq_along(predictors)) {
    v <- predictors[i]
    b <- beta[[v]]
    rr <- ranges[[i]]
    vals <- unique(as.numeric(round(seq(rr[1], rr[2], length.out = 5), 3)))
    if (isTRUE(all.equal(rr[1], rr[2]))) next
    contribution <- b * vals
    pts_v <- (contribution - min(b * rr)) / max_effect * 100
    yy <- y_positions[v]
    axis_lines[[length(axis_lines) + 1]] <- data.frame(
      x = min(pts_v), xend = max(pts_v), y = yy, yend = yy
    )
    axis_ticks[[length(axis_ticks) + 1]] <- data.frame(
      x = pts_v, xend = pts_v, y = yy - 0.08, yend = yy + 0.08
    )
    axis_text[[length(axis_text) + 1]] <- data.frame(
      x = pts_v, y = yy + 0.22, label = format(vals, trim = TRUE, scientific = FALSE)
    )
  }

  # Total points axis.
  total_ticks <- pretty(c(0, total_max), n = 6)
  total_ticks <- total_ticks[total_ticks >= 0 & total_ticks <= total_max]
  total_x <- total_ticks / total_max * 100
  y_total <- y_positions["Total Points"]
  axis_lines[[length(axis_lines) + 1]] <- data.frame(x = 0, xend = 100, y = y_total, yend = y_total)
  axis_ticks[[length(axis_ticks) + 1]] <- data.frame(x = total_x, xend = total_x,
                                                      y = y_total - 0.08, yend = y_total + 0.08)
  axis_text[[length(axis_text) + 1]] <- data.frame(x = total_x, y = y_total + 0.22,
                                                    label = format(round(total_ticks, 1), trim = TRUE))

  # Probability axis derived from the pooled intercept and slopes.
  probs <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.975, 0.99)
  total_for_prob <- (stats::qlogis(probs) - baseline_logit) / (max_effect / 100)
  keep_prob <- is.finite(total_for_prob) & total_for_prob >= 0 & total_for_prob <= total_max
  probs <- probs[keep_prob]
  total_for_prob <- total_for_prob[keep_prob]
  prob_x <- total_for_prob / total_max * 100
  y_prob <- y_positions["Predicted probability"]
  axis_lines[[length(axis_lines) + 1]] <- data.frame(x = 0, xend = 100, y = y_prob, yend = y_prob)
  axis_ticks[[length(axis_ticks) + 1]] <- data.frame(x = prob_x, xend = prob_x,
                                                      y = y_prob - 0.08, yend = y_prob + 0.08)
  axis_text[[length(axis_text) + 1]] <- data.frame(x = prob_x, y = y_prob + 0.22,
                                                    label = format(probs, trim = TRUE))

  lines_df <- do.call(rbind, axis_lines)
  ticks_df <- do.call(rbind, axis_ticks)
  text_df <- do.call(rbind, axis_text)
  label_df <- data.frame(
    y = unname(y_positions),
    label = names(y_positions)
  )

  p <- ggplot() +
    geom_segment(data = lines_df, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.55) +
    geom_segment(data = ticks_df, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.45) +
    geom_text(data = text_df, aes(x = x, y = y, label = label), size = 3, vjust = 0) +
    scale_x_continuous(limits = c(0, 100), breaks = NULL, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = label_df$y, labels = label_df$label,
                       limits = c(0.5, max(label_df$y) + 0.8), expand = c(0, 0)) +
    labs(
      title = paste0("Nomogram based on Rubin-pooled coefficients: ", comparison_name),
      x = NULL, y = NULL,
      caption = "Predictor values are shown on their corresponding axes. Points are calculated from pooled logistic coefficients."
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  print(p)
  ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_Combined_Nomogram.png")),
                  p, width = 9, height = max(5, 0.65 * length(row_names) + 2), dpi = 300)

  # Save the exact coefficients and scaling constants used by the nomogram.
  list(
    plot = p, coefficients = beta, predictor_ranges = ranges,
    max_effect = max_effect, max_points_by_variable = max_points_by_var,
    total_points_max = total_max, baseline_logit = baseline_logit
  )
}

# ------------------------------
# 7. Run all comparisons and models
# ------------------------------

all_results <- list()
all_summaries <- list()
all_pooled_coefficients <- list()
all_predictions <- list()
all_calibration <- list()
all_classification <- list()
all_confusion <- list()
all_dca <- list()
all_pca_tsne <- list()
all_forest_plots <- list()
all_nomograms <- list()

comparison_names <- names(binary_results)
comparison_names <- comparison_names[comparison_names %in% names(comparison_groups)]
if (length(comparison_names) == 0) {
  stop("Screening_Results.rds 中没有与当前三种比较名称匹配的 binary_results。")
}

for (comparison_name in comparison_names) {
  cat("\n\n========================================\n")
  cat("开始分析：", comparison_name, "\n")
  cat("========================================\n")

  selected_markers <- get_final_markers(comparison_name)
  required_markers <- paste0("z_", sub("^z_", "", selected_markers))
  missing_markers <- setdiff(selected_markers, names(mi_long))
  if (length(missing_markers) > 0) {
    stop(comparison_name, " 的筛选标志物在 Z_Norm_Markers.rds 中不存在：",
         paste(missing_markers, collapse = ", "))
  }
  cat("最终纳入标志物：", paste(selected_markers, collapse = ", "), "\n")
  cat("模型协变量：", paste(covariates, collapse = ", "), "\n")

  data_list <- make_comparison_data(
    mi_data_list,
    selected_groups = comparison_groups[[comparison_name]]
  )
  y <- data_list[[1]]$Y_binary
  cat("N = ", length(y), "；对照 = ", sum(y == 0), "；病例 = ", sum(y == 1), "\n", sep = "")

  # PCA/t-SNE are exploratory displays, not pooled inferential tests.
  # The subject-level marker values are averaged across imputations to obtain
  # one coordinate per subject; this can smooth between-imputation variation.
  viz_vars <- selected_markers
  subject_mean_data <- make_subject_mean_data(data_list, viz_vars)
  dim_results <- plot_dimension_reduction(subject_mean_data, viz_vars, comparison_name)
  all_pca_tsne[[comparison_name]] <- dim_results

  model_specs <- list()
  model_labels <- character(0)
  model_specs[["AgeSex_Base"]] <- covariates
  model_labels["AgeSex_Base"] <- "Age + Sex + CCI"

  single_ids <- make.unique(paste0("Single_", make.names(selected_markers)))
  for (i in seq_along(selected_markers)) {
    model_specs[[single_ids[i]]] <- c(covariates, selected_markers[i])
    model_labels[single_ids[i]] <- paste("Covariates +", selected_markers[i])
  }
  model_specs[["Combined"]] <- c(covariates, selected_markers)
  model_labels["Combined"] <- "Covariates + selected markers"

  comparison_results <- list()
  for (model_id in names(model_specs)) {
    predictors <- unique(model_specs[[model_id]])
    check_model_data(data_list, predictors, comparison_name)
    res <- evaluate_mi_model(
      data_list = data_list,
      predictors = predictors,
      comparison_name = comparison_name,
      model_id = model_id,
      model_label = unname(model_labels[model_id]),
      B = B
    )
    comparison_results[[model_id]] <- res

    # Save Rubin-pooled coefficient table for every fitted model.
    pooled_tab <- res$pooled_table
    pooled_tab$Comparison <- comparison_name
    pooled_tab$Model <- model_id
    pooled_tab$Model_Label <- unname(model_labels[model_id])
    all_pooled_coefficients[[paste(comparison_name, model_id, sep = "__")]] <- pooled_tab

    # ROC plot and bootstrap test-set envelope.
    roc_df <- res$roc_ci$main_curve
    band_df <- res$roc_ci$ci_band
    p_roc <- ggplot() +
      geom_ribbon(data = band_df, aes(FPR, ymin = TPR_Lower, ymax = TPR_Upper),
                  fill = "lightblue", alpha = 0.45) +
      geom_line(data = roc_df, aes(FPR, TPR), color = "#D9534F", linewidth = 1.1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      annotate("text", x = 0.60, y = 0.18,
               label = sprintf("AUC = %.3f\n95%% interval: %.3f–%.3f",
                               res$summary$AUC_Optimism_Corrected,
                               res$summary$AUC_Corrected_Lower95,
                               res$summary$AUC_Corrected_Upper95),
               size = 3.5, fontface = "bold") +
      labs(title = paste0(comparison_name, " / ", model_labels[model_id]),
           x = "1 - Specificity", y = "Sensitivity") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    print(p_roc)
    ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_", model_id, "_ROC.png")),
                    p_roc, width = 6.2, height = 5.5, dpi = 300)
    res$roc_plot <- p_roc

    all_summaries[[length(all_summaries) + 1L]] <- res$summary
    all_predictions[[length(all_predictions) + 1L]] <- transform(res$predictions,
                                                                  Comparison = comparison_name,
                                                                  Model = model_id,
                                                                  Model_Label = unname(model_labels[model_id]))
    all_calibration[[length(all_calibration) + 1L]] <- res$calibration_curve
    all_classification[[length(all_classification) + 1L]] <- res$classification_metrics
    all_confusion[[length(all_confusion) + 1L]] <- res$confusion_matrix
    all_dca[[length(all_dca) + 1L]] <- res$dca
  }

  # Calibration plot for all models in this comparison.
  cal_df <- do.call(rbind, lapply(comparison_results, `[[`, "calibration_curve"))
  p_cal <- ggplot(cal_df, aes(Mean_Predicted, Observed_Rate, color = Model_Label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = Observed_Lower95, ymax = Observed_Upper95), width = 0.015, alpha = 0.6) +
    geom_line(aes(group = Model_Label), linewidth = 0.7) +
    geom_point(size = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = paste0("Calibration: ", comparison_name),
         x = "Mean predicted probability", y = "Observed event rate", color = NULL) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  print(p_cal)
  ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_Calibration.png")),
                  p_cal, width = 7.5, height = 5.5, dpi = 300)

  # DCA plot: bootstrap optimism-corrected net benefit and reference strategies.
  dca_df <- do.call(rbind, lapply(comparison_results, `[[`, "dca"))
  prevalence <- mean(y == 1)
  ref_all <- data.frame(Threshold = dca_thresholds,
                        NetBenefit = prevalence - (1 - prevalence) * dca_thresholds / (1 - dca_thresholds),
                        Strategy = "Treat all")
  ref_none <- data.frame(Threshold = dca_thresholds,
                         NetBenefit = 0,
                         Strategy = "Treat none")
  dca_plot_df <- data.frame(
    Threshold = dca_df$Threshold,
    NetBenefit = dca_df$NB_Optimism_Corrected,
    Strategy = dca_df$Model_Label
  )
  dca_plot_df <- rbind(dca_plot_df, ref_all, ref_none)
  p_dca <- ggplot(dca_plot_df, aes(Threshold, NetBenefit, color = Strategy)) +
    geom_line(linewidth = 0.85) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
    coord_cartesian(xlim = c(0.01, 0.99), ylim = c(-0.05, max(dca_df$NB_Optimism_Corrected) * 1.2)) +
    labs(title = paste0("Decision curve analysis: ", comparison_name),
         x = "Threshold probability", y = "Net benefit", color = NULL) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  print(p_dca)
  ggplot2::ggsave(file.path(output_dir, paste0(comparison_name, "_DCA.png")),
                  p_dca, width = 7.5, height = 5.5, dpi = 300)

  all_results[[comparison_name]] <- comparison_results

  # Pooled OR forest plot and nomogram for the final combined model.
  combined <- comparison_results[["Combined"]]
  forest_plot <- plot_pooled_forest(combined$pooled_table, comparison_name, "Combined")
  all_forest_plots[[comparison_name]] <- forest_plot

  # Average the predictor values across imputations for nomogram display ranges only.
  # The regression coefficients and intercept are Rubin-pooled estimates.
  nomogram_ref <- data_list[[1]][, c(id_var, covariates, selected_markers), drop = FALSE]
  for (v in c(covariates, selected_markers)) {
    vv <- do.call(cbind, lapply(data_list, function(d) as.numeric(d[[v]])))
    nomogram_ref[[v]] <- rowMeans(vv)
  }
  nomogram_obj <- make_nomogram(combined$pooled_table, nomogram_ref,
                                c(covariates, selected_markers), comparison_name)
  all_nomograms[[comparison_name]] <- nomogram_obj

  comparison_checkpoint <- list(
    comparison = comparison_name,
    selected_markers = selected_markers,
    covariates = covariates,
    comparison_groups = comparison_groups[[comparison_name]],
    data_by_imputation = data_list,
    dimension_reduction = dim_results,
    model_specs = model_specs,
    model_results = comparison_results,
    calibration_plot = p_cal,
    dca_plot = p_dca,
    pooled_or_forest_plot = forest_plot,
    nomogram = nomogram_obj
  )
  safe_name <- gsub("[^A-Za-z0-9_-]+", "_", comparison_name)
  save_checkpoint(comparison_checkpoint, paste0(safe_name, "_MI_Model_Results.rds"))

  cat("\n模型性能汇总：\n")
  print(do.call(rbind, lapply(comparison_results, `[[`, "summary")), row.names = FALSE)
  cat("\n分类性能：\n")
  print(do.call(rbind, lapply(comparison_results, `[[`, "classification_metrics")), row.names = FALSE)
}

# ------------------------------
# 8. Export pooled results
# ------------------------------

final_model_summary <- do.call(rbind, all_summaries)
final_pooled_coefficients <- do.call(rbind, all_pooled_coefficients)
final_predictions <- do.call(rbind, all_predictions)
final_calibration <- do.call(rbind, all_calibration)
final_classification <- do.call(rbind, all_classification)
final_confusion <- do.call(rbind, all_confusion)
final_dca <- do.call(rbind, all_dca)

write_csv(final_model_summary, "All_Comparisons_MI_Model_Performance.csv")
write_csv(final_pooled_coefficients, "All_Comparisons_Rubin_Pooled_Coefficients.csv")
write_csv(final_predictions, "All_Comparisons_MI_Averaged_Predictions.csv")
write_csv(final_calibration, "All_Comparisons_Calibration_Curve_Data.csv")
write_csv(final_classification, "All_Comparisons_Classification_Metrics.csv")
write_csv(final_confusion, "All_Comparisons_Confusion_Matrices.csv")
write_csv(final_dca, "All_Comparisons_DCA_Data.csv")

final_results <- list(
  input_z_file = z_file,
  input_screening_file = screening_file,
  output_dir = output_dir,
  imputation_ids = imputation_ids,
  M = M,
  N = n_all,
  bootstrap_repetitions = B,
  random_seed = random_seed,
  classification_threshold = classification_threshold,
  covariates = covariates,
  comparison_groups = comparison_groups,
  model_summary = final_model_summary,
  pooled_coefficients = final_pooled_coefficients,
  predictions = final_predictions,
  calibration_curve = final_calibration,
  classification_metrics = final_classification,
  confusion_matrix = final_confusion,
  dca = final_dca,
  results_by_comparison = all_results,
  pca_tsne = all_pca_tsne,
  pooled_or_forest_plots = all_forest_plots,
  nomograms = all_nomograms
)
save_checkpoint(final_results, "Var_Model_MI_Rubin_Results.rds")

cat("\n========================================\n")
cat("MI Rubin 合并及模型评估完成。\n")
cat("插补集数量：", M, "\n")
cat("每个插补集病例数：", n_all, "\n")
cat("Bootstrap 次数：", B, "\n")
cat("输出目录：", output_dir, "\n")
cat("========================================\n")
