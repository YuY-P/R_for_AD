# ============================================================
# Var_Screen.R
#
# Purpose:
#   Initial variable screening, LASSO-based variable selection,
#   ROC evaluation, and forest plot generation for odds ratios (ORs).
#     1-2. Input and initiation
#     3. Functions
#       select_candidates
#       make_stratified_foldid
#       extract_selected_variables
#       plot_lasso_path
#       run_lasso
#       build_model_list
#       plot_roc_overlay
#       evaluate_multiclass_models
#       run_binary_univariate
#       evaluate_binary_models
#     4. 4-Class variable LASSO Screening
#     5. ROC assessment
#     6. 2-Class variable LASSO Screening
#     7. ROC assessment
#     8. Forest plot
#     9. Output
#
# Input:
#   Z-score-standardized variables.
#     Z_Norm_Markers.rds
#
# Analysis:
#   Uses imputation_id = 1 only.
#   Three binary comparisons are performed independently:
#   1) 0 vs 1+2+3 (group 4 excluded)
#   2) 0 vs 3+4 (groups 1 and 2 excluded)
#   3) 0 vs 1+2+3+4 (all groups included)
#
# Output:
#   ...
#     Screening_Results.rds
#
# Note:
#   This script uses only the first completed imputed dataset
#   (imputation_id = 1) for exploratory analyses. It does not
#   pool estimates across imputations.
# ============================================================

# 1. 初始化环境与加载程序包 ------------------------------------
rm(list = ls())
gc()
graphics.off()

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(nnet)
})

input_file <- "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/Z_Norm_Markers.rds"
output_dir <- "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_result <- function(object, filename) {
  saveRDS(object, file.path(output_dir, filename))
}
write_result <- function(object, filename) {
  write.csv(object, file.path(output_dir, filename), row.names = FALSE,
            fileEncoding = "UTF-8")
}

# 2. 数据读取与基础设置 ----------------------------------------
df_final <- readRDS(input_file)
id_var <- "No"
group_var <- "Tag_No"
imp_var <- "imputation_id"
covariates <- c("z_Age", "Sex", "CCI")

# 四分类结局保持为：0、1+2、3、4
group_levels <- c("0", "1+2", "3", "4")

# 三种独立二分类比较。每次只纳入列表中指定的原始组别。
binary_comparisons <- list(
  "0_vs_1plus2plus3" = c("0", "1", "2", "3"),
  "0_vs_3plus4" = c("0", "3", "4"),
  "0_vs_1plus2plus3plus4" = c("0", "1", "2", "3", "4")
)

required_vars <- c(imp_var, group_var, covariates)
missing_vars <- setdiff(required_vars, names(df_final))
if (length(missing_vars) > 0) {
  stop("数据中缺少必要变量：", paste(missing_vars, collapse = ", "))
}

# 仅提取第1个插补数据集
df_1 <- df_final %>%
  dplyr::filter(.data[[imp_var]] == 1) %>%
  droplevels()
if (nrow(df_1) == 0) stop("未找到 imputation_id = 1 的数据。")

non_numeric_covariates <- covariates[
  !vapply(df_1[covariates], is.numeric, logical(1))
]
if (length(non_numeric_covariates) > 0) {
  stop("以下协变量不是数值型：", paste(non_numeric_covariates, collapse = ", "))
}

# 自动识别数值型Z-score标志物，排除协变量
exclude_z_vars <- covariates
markers <- names(df_1)[
  grepl("^z_", names(df_1)) &
    vapply(df_1, is.numeric, logical(1))
]
markers <- setdiff(markers, exclude_z_vars)
if (length(markers) == 0) stop("未识别到Z-score标准化的标志物变量。")

cat("\n================ 标志物变量 ================\n")
print(markers)

df_1$Group_original <- as.character(df_1[[group_var]])
cat("\n================ 原始分组情况 ================\n")
print(table(df_1$Group_original, useNA = "ifany"))

# 构建四分类结局
df_1$Y4 <- NA_character_
df_1$Y4[df_1$Group_original == "0"] <- "0"
df_1$Y4[df_1$Group_original %in% c("1", "2")] <- "1+2"
df_1$Y4[df_1$Group_original == "3"] <- "3"
df_1$Y4[df_1$Group_original == "4"] <- "4"
df_1$Y4 <- factor(df_1$Y4, levels = group_levels)

cat("\n================ 四分类分组 ================\n")
print(table(df_1$Y4, useNA = "ifany"))
if (any(table(df_1$Y4) == 0)) stop("四分类中至少有一个类别没有样本。")

# 四分类及各二分类分析均基于同一套完整病例数据，保证模型间样本一致
model_vars <- c("Group_original", "Y4", covariates, markers)
df_model <- df_1 %>% dplyr::select(dplyr::all_of(model_vars))
valid_rows <- complete.cases(df_model)
for (v in c(covariates, markers)) {
  valid_rows <- valid_rows & is.finite(df_model[[v]])
}
df_model <- df_model[valid_rows, , drop = FALSE]
df_model$Y4 <- factor(df_model$Y4, levels = group_levels)

# 排除共同建模数据中的零方差变量
zero_var <- vapply(
  df_model[, c(covariates, markers), drop = FALSE],
  function(x) length(unique(x)) < 2,
  logical(1)
)
if (any(zero_var)) {
  excluded_vars <- names(zero_var)[zero_var]
  cat("\n以下变量无变异，将从后续模型中排除：\n")
  print(excluded_vars)
  markers <- setdiff(markers, excluded_vars)
  covariates <- setdiff(covariates, excluded_vars)
}
if (length(markers) == 0) stop("排除无变异变量后，没有可用的标志物。")
if (length(unique(df_model$Y4)) != 4) stop("删除缺失值后，四分类结局不足4个类别。")

cat("\n================ 共同建模数据 ================\n")
cat("有效样本量：", nrow(df_model), "\n")
cat("标志物数量：", length(markers), "\n")
cat("协变量：", paste(covariates, collapse = ", "), "\n")
print(table(df_model$Y4))

# miRNA名称兼容 miRNA、miR 等常见写法，如 z_miR34a
miRNA_markers <- markers[grepl("miRNA|miR", markers, ignore.case = TRUE)]
other_markers <- setdiff(markers, miRNA_markers)
cat("\n识别到的miRNA标志物：\n")
print(miRNA_markers)

# 检查每种二分类比较是否同时包含对照组与病例组
for (nm in names(binary_comparisons)) {
  included <- binary_comparisons[[nm]]
  absent <- setdiff(included, unique(df_model$Group_original))
  if (length(absent) > 0) {
    warning(nm, "中以下原始组别在共同建模数据中不存在：",
            paste(absent, collapse = ", "))
  }
}

data_preparation <- list(
  df_model = df_model,
  markers = markers,
  covariates = covariates,
  miRNA_markers = miRNA_markers,
  other_markers = other_markers,
  group_levels = group_levels,
  binary_comparisons = binary_comparisons,
  imputation_id = 1
)
save_result(data_preparation, "Screening_Data.rds")

# 3. 通用函数 --------------------------------------------------

# 3.1 根据单变量P值确定候选标志物
select_candidates <- function(uni_results, analysis_name) {
  valid_results <- uni_results %>%
    dplyr::filter(is.finite(P_value)) %>%
    dplyr::arrange(P_value)
  candidate <- valid_results %>%
    dplyr::filter(P_value < 0.10) %>%
    dplyr::pull(Marker)
  if (length(candidate) < 4) {
    cat("\n", analysis_name,
        "：P < 0.10的变量不足4个，按P值顺序补足至最多4个。\n", sep = "")
    candidate <- head(valid_results$Marker, 4)
  }
  if (length(candidate) == 0) stop(analysis_name, "：未获得有效候选标志物。")
  candidate
}

# 3.2 构建分层交叉验证fold
make_stratified_foldid <- function(y, nfolds, seed) {
  set.seed(seed)
  foldid <- integer(length(y))
  y_char <- as.character(y)
  for (g in unique(y_char)) {
    idx <- which(y_char == g)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  foldid
}

# 3.3 提取LASSO非零系数变量
extract_selected_variables <- function(coef_object) {
  if (is.list(coef_object)) {
    selected <- unlist(lapply(coef_object, function(x) {
      x <- as.matrix(x)
      rownames(x)[as.numeric(x) != 0]
    }))
  } else {
    coef_object <- as.matrix(coef_object)
    selected <- rownames(coef_object)[as.numeric(coef_object) != 0]
  }
  unique(setdiff(selected, "(Intercept)"))
}

# 3.4 绘制LASSO系数路径
plot_lasso_path <- function(fit, family, main_title) {
  
  lambda_values <- fit$lambda
  
  if (is.null(lambda_values) || length(lambda_values) == 0) {
    warning("LASSO模型中没有可用的lambda值：", main_title)
    return(invisible(NULL))
  }
  
  log_lambda <- log(lambda_values)
  path_data_list <- list()
  
  # multinomial：fit$beta为按结局类别分别存储的系数矩阵
  if (family == "multinomial") {
    
    if (!is.list(fit$beta)) {
      warning("多分类LASSO模型的系数结构异常：", main_title)
      return(invisible(NULL))
    }
    
    class_names <- names(fit$beta)
    
    if (is.null(class_names) || any(class_names == "")) {
      class_names <- paste0("Class_", seq_along(fit$beta))
    }
    
    for (i in seq_along(fit$beta)) {
      
      coef_matrix <- as.matrix(fit$beta[[i]])
      
      if (ncol(coef_matrix) != length(lambda_values)) {
        warning("系数矩阵与lambda数量不一致：", class_names[i])
        next
      }
      
      variable_names <- rownames(coef_matrix)
      
      if (is.null(variable_names)) next
      
      for (j in seq_along(variable_names)) {
        
        path_data_list[[length(path_data_list) + 1]] <- data.frame(
          Log_Lambda = log_lambda,
          Coefficient = as.numeric(coef_matrix[j, ]),
          Path_Label = variable_names[j],
          Class = class_names[i],
          stringsAsFactors = FALSE
        )
      }
    }
    
  } else {
    
    # binomial：fit$beta为变量×lambda的系数矩阵
    coef_matrix <- as.matrix(fit$beta)
    
    if (ncol(coef_matrix) != length(lambda_values)) {
      warning("系数矩阵与lambda数量不一致：", main_title)
      return(invisible(NULL))
    }
    
    variable_names <- rownames(coef_matrix)
    
    if (is.null(variable_names)) {
      warning("无法识别LASSO系数对应的变量名：", main_title)
      return(invisible(NULL))
    }
    
    for (j in seq_along(variable_names)) {
      
      path_data_list[[length(path_data_list) + 1]] <- data.frame(
        Log_Lambda = log_lambda,
        Coefficient = as.numeric(coef_matrix[j, ]),
        Path_Label = variable_names[j],
        Class = "Binary",
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(path_data_list) == 0) {
    warning("未提取到可用于绘制的LASSO系数路径：", main_title)
    return(invisible(NULL))
  }
  
  path_data <- dplyr::bind_rows(path_data_list)
  
  # 去除非有限值，避免ggplot绘图异常
  path_data <- path_data %>%
    dplyr::filter(
      is.finite(Log_Lambda),
      is.finite(Coefficient)
    )
  
  if (nrow(path_data) == 0) {
    warning("LASSO系数路径数据为空：", main_title)
    return(invisible(NULL))
  }
  
  if (family == "multinomial") {
    
    p <- ggplot(
      path_data,
      aes(
        x = Log_Lambda,
        y = Coefficient,
        group = interaction(Path_Label, Class),
        color = Path_Label,
        linetype = Class
      )
    ) +
      geom_line(linewidth = 0.7, show.legend = TRUE) +
      scale_color_discrete(name = "Variables") +
      scale_linetype_discrete(name = "Outcome\nclass") +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey60"
      ) +
      labs(
        title = main_title,
        x = "Log(lambda)",
        y = "Coefficient"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0)
      )
    
  } else {
    
    p <- ggplot(
      path_data,
      aes(
        x = Log_Lambda,
        y = Coefficient,
        group = Path_Label,
        color = Path_Label
      )
    ) +
      geom_line(linewidth = 0.7, show.legend = TRUE) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey60"
      ) +
      labs(
        title = main_title,
        x = "Log(lambda)",
        y = "Coefficient"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"), # 缩小图例符号
        legend.margin = margin(0, 0, 0, 0)
      )
  }
  
  print("===== path_data 的列名 =====")
  print(names(path_data))
  print("===== path_data 的结构 =====")
  str(path_data)
  print("===== path_data 的前几行 =====")
  print(head(path_data))
  
  print(p)
  invisible(p)
}

# 3.5 LASSO分析
run_lasso <- function(df, outcome, candidate_markers, covariates,
                      family, analysis_name,
                      seed_list = c(42, 123, 456, 789, 999)) {
  
  cat("\n================ ", analysis_name,
      "：LASSO分析 ================\n", sep = "")
  
  if (length(candidate_markers) == 0) {
    stop(analysis_name, "：没有候选标志物可用于LASSO。")
  }
  
  # 检查候选变量是否存在
  final_vars <- c(covariates, candidate_markers)
  missing_vars <- setdiff(final_vars, names(df))
  
  if (length(missing_vars) > 0) {
    stop(
      analysis_name, "：以下变量在数据中不存在：",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  # 构建自变量矩阵和结局变量
  X <- as.matrix(df[, final_vars, drop = FALSE])
  Y <- df[[outcome]]
  
  storage.mode(X) <- "double"
  
  if (any(!is.finite(X))) {
    stop(analysis_name, "：LASSO输入矩阵中存在非有限值。")
  }
  
  if (anyNA(Y)) {
    stop(analysis_name, "：结局变量中存在缺失值。")
  }
  
  # 协变量不惩罚，候选标志物进行LASSO惩罚
  penalty_factors <- c(
    rep(0, length(covariates)),
    rep(1, length(candidate_markers))
  )
  
  # 根据最小类别样本量确定交叉验证折数
  class_counts <- table(Y)
  min_class_n <- min(class_counts)
  nfolds_use <- min(5, min_class_n)
  
  if (nfolds_use < 3) {
    stop(
      analysis_name,
      "：最小类别样本量不足，无法进行至少3折交叉验证。"
    )
  }
  
  cat("交叉验证折数：", nfolds_use, "\n")
  cat("各类别样本量：\n")
  print(class_counts)
  
  # 统一设置模型参数
  # 输入变量已完成Z-score标准化，因此不再进行内部标准化
  # 增大maxit以改善部分lambda下的收敛情况
  fit_args <- list(
    x = X,
    y = Y,
    alpha = 1,
    penalty.factor = penalty_factors,
    standardize = FALSE,
    maxit = 1000000
  )
  
  if (family == "multinomial") {
    
    fit_args$family <- "multinomial"
    fit_args$type.multinomial <- "grouped"
    
  } else if (family == "binomial") {
    
    fit_args$family <- "binomial"
    
  } else {
    
    stop("不支持的LASSO family：", family)
  }
  
  # 交叉验证
  foldid <- make_stratified_foldid(
    Y,
    nfolds = nfolds_use,
    seed = 123
  )
  
  if (family == "multinomial") {
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = X, y = Y, family = "multinomial", type.multinomial = "grouped",
        alpha = 1, penalty.factor = penalty_factors,
        standardize = FALSE, maxit = 1000000,
        foldid = foldid, type.measure = "deviance"
      ),
      error = function(e) stop(analysis_name, "：LASSO交叉验证拟合失败：", e$message)
    )
  } else {
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = X, y = Y, family = "binomial",
        alpha = 1, penalty.factor = penalty_factors,
        standardize = FALSE, maxit = 1000000,
        foldid = foldid, type.measure = "deviance"
      ),
      error = function(e) stop(analysis_name, "：LASSO交叉验证拟合失败：", e$message)
    )
  }
  
  # 使用全部数据拟合完整lambda路径
  if (family == "multinomial") {
    fit <- tryCatch(
      glmnet::glmnet(
        x = X, y = Y, family = "multinomial", type.multinomial = "grouped",
        alpha = 1, penalty.factor = penalty_factors,
        standardize = FALSE, maxit = 1000000
      ),
      error = function(e) stop(analysis_name, "：LASSO完整路径拟合失败：", e$message)
    )
  } else {
    fit <- tryCatch(
      glmnet::glmnet(
        x = X, y = Y, family = "binomial",
        alpha = 1, penalty.factor = penalty_factors,
        standardize = FALSE, maxit = 1000000
      ),
      error = function(e) stop(analysis_name, "：LASSO完整路径拟合失败：", e$message)
    )
  }
  
  # 绘制交叉验证曲线和系数路径
  old_par <- par(no.readonly = TRUE)
  par(mar = c(3.5, 3.5, 4, 2) + 0.1,
      mgp = c(2.2, 0.3, 0))
  plot(
    cv_fit,
    main = ""
  )
  title(
    main = paste0(analysis_name, "\n-LASSO CV Curve"),
    cex.main = 1,
    font.main = 2,
    line = 1.8
  )
  par(old_par)
  
  plot_lasso_path(
    fit,
    family,
    paste0(analysis_name, "\n - LASSO Coefficient Path")
  )
  
  # 提取lambda.min和lambda.1se对应的非零系数
  coef_min <- coef(cv_fit, s = "lambda.min")
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  
  selected_min <- extract_selected_variables(coef_min)
  selected_1se <- extract_selected_variables(coef_1se)
  
  # 仅保留候选标志物，不将未惩罚协变量作为标志物输出
  selected_markers <- intersect(
    selected_min,
    candidate_markers
  )
  
  cat("\nLambda.min选中的变量：\n")
  print(selected_min)
  
  cat("\nLambda.min选中的标志物：\n")
  print(selected_markers)
  
  cat("\nLambda.1se选中的变量：\n")
  print(selected_1se)
  
  # 整理lambda.min下的系数
  if (is.list(coef_min)) {
    
    coef_df <- dplyr::bind_rows(
      lapply(names(coef_min), function(class_name) {
        
        x <- as.matrix(coef_min[[class_name]])
        
        data.frame(
          Class = class_name,
          Variable = rownames(x),
          Coefficient = as.numeric(x),
          stringsAsFactors = FALSE
        )
      })
    )
    
  } else {
    
    x <- as.matrix(coef_min)
    
    coef_df <- data.frame(
      Variable = rownames(x),
      Coefficient = as.numeric(x),
      stringsAsFactors = FALSE
    )
  }
  
  # LASSO稳定性检验：不同随机种子下重复交叉验证
  cat("\n--- LASSO稳定性检验 ---\n")
  
  stable_selected <- list()
  
  for (s in seed_list) {
    
    cat("当前随机种子：", s, "\n")
    
    fold_tmp <- make_stratified_foldid(
      Y,
      nfolds = nfolds_use,
      seed = s
    )
    
    if (family == "multinomial") {
      cv_tmp <- tryCatch(
        glmnet::cv.glmnet(
          x = X, y = Y, family = "multinomial", type.multinomial = "grouped",
          alpha = 1, penalty.factor = penalty_factors,
          standardize = FALSE, maxit = 1000000,
          foldid = fold_tmp, type.measure = "deviance"
        ),
        error = function(e) NULL
      )
    } else {
      cv_tmp <- tryCatch(
        glmnet::cv.glmnet(
          x = X, y = Y, family = "binomial",
          alpha = 1, penalty.factor = penalty_factors,
          standardize = FALSE, maxit = 1000000,
          foldid = fold_tmp, type.measure = "deviance"
        ),
        error = function(e) NULL
      )
    }
    
    if (is.null(cv_tmp)) {
      stable_selected[[as.character(s)]] <- character(0)
    } else {
      stable_selected[[as.character(s)]] <-
        extract_selected_variables(
          coef(cv_tmp, s = "lambda.min")
        )
    }
  }
  
  # 统计不同随机种子下的变量选择频率
  all_variables <- unique(unlist(stable_selected))
  
  if (length(all_variables) > 0) {
    
    selection_frequency_df <- data.frame(
      Variable = all_variables,
      Selection_Count = vapply(
        all_variables,
        function(v) {
          sum(
            vapply(
              stable_selected,
              function(x) v %in% x,
              logical(1)
            )
          )
        },
        numeric(1)
      )
    )
    
    selection_frequency_df$Selection_Frequency <-
      selection_frequency_df$Selection_Count /
      length(seed_list)
    
    selection_frequency_df <- selection_frequency_df %>%
      dplyr::arrange(
        dplyr::desc(Selection_Frequency)
      )
    
  } else {
    
    selection_frequency_df <- data.frame(
      Variable = character(),
      Selection_Count = integer(),
      Selection_Frequency = numeric()
    )
  }
  
  cat("\n不同随机种子下的变量选择频率：\n")
  print(selection_frequency_df)
  
  # 保存结果
  lasso_result <- list(
    cv_fit = cv_fit,
    fit = fit,
    coef_min = coef_min,
    coef_1se = coef_1se,
    coef_df = coef_df,
    selected_lasso = selected_min,
    selected_markers = selected_markers,
    stable_selected = stable_selected,
    selection_frequency = selection_frequency_df,
    nfolds = nfolds_use
  )
  
  save_result(
    lasso_result,
    paste0(analysis_name, "_LASSO_Results.rds")
  )
  
  write_result(
    coef_df,
    paste0(analysis_name, "_LASSO_Coefficients.csv")
  )
  
  write_result(
    selection_frequency_df,
    paste0(analysis_name, "_LASSO_Selection_Frequency.csv")
  )
  
  lasso_result
}

# 3.6 构建ROC模型组合
build_model_list <- function(markers, miRNA_markers, other_markers, selected_markers) {
  model_list <- list(Model_Combined = markers)
  if (length(miRNA_markers) > 0) model_list$Model_miRNA <- miRNA_markers
  if (length(other_markers) > 0 && !setequal(other_markers, markers)) {
    model_list$Model_Other_Markers <- other_markers
  }
  if (length(selected_markers) > 0) model_list$Model_LASSO_Selected <- selected_markers
  signatures <- vapply(model_list, function(x) paste(sort(unique(x)), collapse = ";"),
                       character(1))
  model_list[!duplicated(signatures)]
}

# 3.7 ROC曲线叠加绘制
plot_roc_overlay <- function(roc_list, auc_table, title) {
  if (length(roc_list) == 0) {
    warning("没有可绘制的ROC曲线：", title)
    return(invisible(NULL))
  }
  model_names <- names(roc_list)
  colors <- grDevices::rainbow(length(model_names))
  plot(roc_list[[1]], col = colors[1], lwd = 2, legacy.axes = TRUE, main = title)
  if (length(roc_list) > 1) {
    for (i in 2:length(roc_list)) plot(roc_list[[i]], add = TRUE,
                                       col = colors[i], lwd = 2)
  }
  abline(a = 0, b = 1, lty = 2, col = "grey60")
  # 使用简短模型名称，避免图例过长
  short_names <- c(
    "Model_Combined" = "Combined",
    "Model_miRNA" = "miRNA",
    "Model_Other_Markers" = "Other markers",
    "Model_LASSO_Selected" = "LASSO"
  )
  
  auc_labels <- vapply(model_names, function(model_name) {
    auc_val <- auc_table$AUC[match(model_name, auc_table$Model)]
    
    display_name <- if (model_name %in% names(short_names)) {
      short_names[[model_name]]
    } else {
      model_name
    }
    
    paste0(display_name, " (AUC=", sprintf("%.3f", auc_val), ")")
  }, character(1))
  
  legend(
    "bottomright",
    legend = auc_labels,
    col = colors,
    lwd = 2,
    cex = 0.6,
    bty = "n",
    ncol = 1,
    y.intersp = 0.7,
    x.intersp = 0.5,
    inset = c(0.02, 0.05)
  )
}

# 3.8 多分类模型拟合与ROC分析
evaluate_multiclass_models <- function(df, model_list, outcome,
                                        class_levels, analysis_name) {
  predictions <- list(); roc_results <- list(); fitted_models <- list()
  auc_results <- data.frame(Model = character(), Class = character(), AUC = numeric())
  model_variables <- data.frame(
    Model = names(model_list),
    Markers = vapply(model_list, function(x) paste(x, collapse = "; "), character(1))
  )

  for (model_name in names(model_list)) {
    formula_model <- reformulate(c(covariates, model_list[[model_name]]), response = outcome)
    cat("\n拟合四分类模型：", model_name, "\n"); print(formula_model)
    fit_model <- tryCatch(
      suppressWarnings(nnet::multinom(formula_model, data = df,
                                      trace = FALSE, maxit = 1000)),
      error = function(e) { warning(model_name, "模型拟合失败：", e$message); NULL }
    )
    if (is.null(fit_model)) {
      cat("\n模型拟合失败，跳过：", model_name, "\n")
      next
    } else {
      cat("\n模型拟合成功：", model_name, "\n")
    }
    fitted_models[[model_name]] <- fit_model
    pred_prob <- as.matrix(predict(fit_model, newdata = df, type = "probs"))
    predictions[[model_name]] <- pred_prob
    roc_results[[model_name]] <- list()

    for (class_name in class_levels) {
      if (!class_name %in% colnames(pred_prob)) next
      binary_y <- as.numeric(df[[outcome]] == class_name)
      roc_obj <- tryCatch(
        pROC::roc(response = binary_y, predictor = pred_prob[, class_name],
                  levels = c(0, 1), direction = "<", quiet = TRUE),
        error = function(e) {
          warning(model_name, " / ", class_name, " ROC计算失败：", e$message); NULL
        }
      )
      if (is.null(roc_obj)) next
      roc_results[[model_name]][[class_name]] <- roc_obj
      auc_results <- dplyr::bind_rows(
        auc_results,
        data.frame(Model = model_name, Class = class_name,
                   AUC = as.numeric(pROC::auc(roc_obj)))
      )
    }
  }

  for (class_name in class_levels) {
    valid_models <- names(roc_results)[vapply(
      roc_results, function(x) !is.null(x[[class_name]]), logical(1)
    )]
    if (length(valid_models) == 0) next
    roc_list_class <- lapply(valid_models, function(x) roc_results[[x]][[class_name]])
    names(roc_list_class) <- valid_models
    auc_table_class <- auc_results %>% dplyr::filter(Class == class_name)
    plot_roc_overlay(roc_list_class, auc_table_class,
                     paste0(analysis_name, ": Class ", class_name, " vs Rest"))
  }
  
  cat("\n========== 四分类模型拟合情况 ==========\n")
  cat("成功拟合的模型：\n")
  print(names(fitted_models))
  
  cat("\n各模型成功计算ROC的类别：\n")
  for (nm in names(roc_results)) {
    cat(nm, "：", paste(names(roc_results[[nm]]), collapse = ", "), "\n")
  }

  result <- list(fitted_models = fitted_models, predictions = predictions,
                 roc_results = roc_results, auc_results = auc_results,
                 model_variables = model_variables)
  save_result(result, paste0(analysis_name, "_ROC_Results.rds"))
  write_result(model_variables, paste0(analysis_name, "_Model_Variables.csv"))
  write_result(auc_results, paste0(analysis_name, "_ROC_AUC.csv"))
  result
}

# 3.9 二分类单变量Logistic回归
run_binary_univariate <- function(df, markers, outcome) {
  results <- lapply(markers, function(m) {
    fit <- tryCatch(
      suppressWarnings(glm(reformulate(m, response = outcome), data = df,
                           family = binomial)),
      error = function(e) { warning(m, "单变量Logistic模型拟合失败：", e$message); NULL }
    )
    p_value <- NA_real_
    if (!is.null(fit)) {
      coef_table <- summary(fit)$coefficients
      if (m %in% rownames(coef_table)) p_value <- coef_table[m, "Pr(>|z|)"]
    }
    data.frame(Marker = m, P_value = p_value)
  })
  dplyr::bind_rows(results) %>% dplyr::arrange(P_value)
}

# 3.10 二分类模型拟合与ROC分析
evaluate_binary_models <- function(df, model_list, outcome, analysis_name) {
  fitted_models <- list(); predictions <- list(); roc_results <- list()
  auc_results <- data.frame(Model = character(), AUC = numeric())
  model_variables <- data.frame(
    Model = names(model_list),
    Markers = vapply(model_list, function(x) paste(x, collapse = "; "), character(1))
  )

  for (model_name in names(model_list)) {
    formula_model <- reformulate(c(covariates, model_list[[model_name]]), response = outcome)
    cat("\n拟合模型：", model_name, "\n"); print(formula_model)
    fit_model <- tryCatch(
      suppressWarnings(glm(formula_model, data = df, family = binomial)),
      error = function(e) { warning(model_name, "模型拟合失败：", e$message); NULL }
    )
    if (is.null(fit_model)) next
    fitted_models[[model_name]] <- fit_model
    pred_prob <- predict(fit_model, type = "response")
    predictions[[model_name]] <- pred_prob
    roc_obj <- tryCatch(
      pROC::roc(response = df[[outcome]], predictor = pred_prob,
                levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) { warning(model_name, " ROC计算失败：", e$message); NULL }
    )
    if (is.null(roc_obj)) next
    roc_results[[model_name]] <- roc_obj
    auc_results <- dplyr::bind_rows(
      auc_results,
      data.frame(Model = model_name, AUC = as.numeric(pROC::auc(roc_obj)))
    )
  }

  # 比较各模型在同一批样本上的预测概率
  if (length(predictions) >= 2) {
    model_pairs <- combn(names(predictions), 2, simplify = FALSE)
    cat("\n================ 预测概率比较 ================\n")
    for (pair in model_pairs) {
      pred_a <- predictions[[pair[1]]]; pred_b <- predictions[[pair[2]]]
      cat("\n", pair[1], " vs ", pair[2], "\n", sep = "")
      cat("预测概率是否相同：", isTRUE(all.equal(pred_a, pred_b)), "\n")
      cat("预测概率最大绝对差值：", max(abs(pred_a - pred_b)), "\n")
      cat("预测概率排序是否相同：",
          identical(rank(pred_a, ties.method = "average"),
                    rank(pred_b, ties.method = "average")), "\n")
      auc_a <- auc_results$AUC[match(pair[1], auc_results$Model)]
      auc_b <- auc_results$AUC[match(pair[2], auc_results$Model)]
      cat("AUC绝对差值：", abs(auc_a - auc_b), "\n")
    }
  }
  if (length(roc_results) > 0) {
    plot_roc_overlay(roc_results, auc_results, paste0(analysis_name, " ROC Curves"))
  }
  result <- list(fitted_models = fitted_models, predictions = predictions,
                 roc_results = roc_results, auc_results = auc_results,
                 model_variables = model_variables)
  save_result(result, paste0(analysis_name, "_ROC_Results.rds"))
  write_result(model_variables, paste0(analysis_name, "_Model_Variables.csv"))
  write_result(auc_results, paste0(analysis_name, "_ROC_AUC.csv"))
  result
}

# 4. 四分类分析：单变量筛选与LASSO ------------------------------
cat("\n\n##################################################\n")
cat("########## 四分类分析：0、1+2、3、4 ##########\n")
cat("##################################################\n")

uni_results_4 <- dplyr::bind_rows(lapply(markers, function(m) {
  fit_full <- tryCatch(
    suppressWarnings(nnet::multinom(reformulate(m, response = "Y4"),
                                    data = df_model, trace = FALSE, maxit = 1000)),
    error = function(e) { warning(m, "单变量多项Logistic模型拟合失败：", e$message); NULL }
  )
  fit_null <- tryCatch(
    suppressWarnings(nnet::multinom(Y4 ~ 1, data = df_model,
                                    trace = FALSE, maxit = 1000)),
    error = function(e) NULL
  )
  p_value <- NA_real_
  if (!is.null(fit_full) && !is.null(fit_null)) {
    ll_full <- logLik(fit_full); ll_null <- logLik(fit_null)
    LR <- 2 * (as.numeric(ll_full) - as.numeric(ll_null))
    df_diff <- attr(ll_full, "df") - attr(ll_null, "df")
    if (is.finite(LR) && LR >= 0 && df_diff > 0) {
      p_value <- pchisq(LR, df = df_diff, lower.tail = FALSE)
    }
  }
  data.frame(Marker = m, P_value = p_value)
})) %>% dplyr::arrange(P_value)

cat("\n四分类单变量多项Logistic回归结果：\n"); print(uni_results_4)
candidate_markers_4 <- select_candidates(uni_results_4, "四分类")
cat("\n四分类候选标志物：\n"); print(candidate_markers_4)

lasso_4 <- run_lasso(df = df_model, outcome = "Y4",
                     candidate_markers = candidate_markers_4,
                     covariates = covariates, family = "multinomial",
                     analysis_name = "FourClass")
selected_markers_4 <- lasso_4$selected_markers
fourclass_screening <- list(uni_results = uni_results_4,
                            candidate_markers = candidate_markers_4,
                            lasso = lasso_4,
                            selected_markers = selected_markers_4)
save_result(fourclass_screening, "FourClass_Screening_Results.rds")
write_result(uni_results_4, "FourClass_Univariate_Screening.csv")

# 5. 四分类分析：模型组合与ROC评价 ------------------------------
model_list_4 <- build_model_list(markers, miRNA_markers,
                                 other_markers, selected_markers_4)
cat("\n四分类模型组合：\n"); print(model_list_4)

cat("\n========== 四分类 ROC 模型变量组成 ==========\n")
for (nm in names(model_list_4)) {
  cat("\n模型：", nm, "\n")
  cat("标志物数量：", length(model_list_4[[nm]]), "\n")
  print(model_list_4[[nm]])
}

result_4 <- evaluate_multiclass_models(
  df = df_model, model_list = model_list_4, outcome = "Y4",
  class_levels = group_levels, analysis_name = "FourClass"
)
predictions_4 <- result_4$predictions
roc_results_4 <- result_4$roc_results
auc_results_4 <- result_4$auc_results
cat("\n四分类AUC结果：\n"); print(auc_results_4)
fourclass_results <- list(model_list = model_list_4, result = result_4,
                          predictions = predictions_4,
                          roc_results = roc_results_4,
                          auc_results = auc_results_4)
save_result(fourclass_results, "FourClass_Model_Results.rds")

# 6. 二分类分析：三种独立比较 ----------------------------------
cat("\n\n##################################################\n")
cat("########## 二分类分析：三种独立比较 ##########\n")
cat("##################################################\n")

binary_results <- list()
for (comparison_name in names(binary_comparisons)) {
  cat("\n\n============================================\n")
  cat("当前二分类比较：", comparison_name, "\n")
  cat("============================================\n")

  selected_groups <- binary_comparisons[[comparison_name]]
  df_binary <- df_model %>%
    dplyr::filter(Group_original %in% selected_groups) %>%
    droplevels()

  # 每个比较中：原始0组编码为0，其余纳入组别统一编码为1
  df_binary$Y_binary <- as.integer(df_binary$Group_original != "0")

  cat("\n纳入的原始组别：\n"); print(selected_groups)
  cat("\n原始分组分布：\n"); print(table(df_binary$Group_original))
  cat("\n二分类结局分布（0=对照，1=病例）：\n"); print(table(df_binary$Y_binary))
  cat("样本量：", nrow(df_binary), "\n")

  if (length(unique(df_binary$Y_binary)) != 2 ||
      any(table(factor(df_binary$Y_binary, levels = c(0, 1))) == 0)) {
    warning(comparison_name, "未同时包含对照组和病例组，跳过。")
    next
  }

  # 6.1 单变量Logistic回归
  uni_results_2 <- run_binary_univariate(df_binary, markers, "Y_binary")
  cat("\n单变量Logistic回归结果：\n"); print(uni_results_2)
  candidate_markers_2 <- select_candidates(uni_results_2, comparison_name)
  cat("\n候选标志物：\n"); print(candidate_markers_2)

  save_result(list(data = df_binary, univariate = uni_results_2,
                   candidates = candidate_markers_2),
              paste0(comparison_name, "_Screening_Results.rds"))
  write_result(uni_results_2, paste0(comparison_name, "_Univariate_Screening.csv"))

  # 6.2 Binomial LASSO
  lasso_2 <- run_lasso(df = df_binary, outcome = "Y_binary",
                       candidate_markers = candidate_markers_2,
                       covariates = covariates, family = "binomial",
                       analysis_name = paste0("LASSO_", comparison_name))
  selected_markers_2 <- lasso_2$selected_markers
  cat("\nLASSO选中的标志物：\n"); print(selected_markers_2)

  # 6.3 模型组合
  model_list_2 <- build_model_list(markers, miRNA_markers,
                                   other_markers, selected_markers_2)
  cat("\n模型组合：\n"); print(model_list_2)

  # 6.4 模型拟合及ROC评价
  result_2 <- evaluate_binary_models(df = df_binary,
                                     model_list = model_list_2,
                                     outcome = "Y_binary",
                                     analysis_name = comparison_name)
  cat("\nAUC结果：\n"); print(result_2$auc_results)

  binary_results[[comparison_name]] <- list(
    data = df_binary,
    univariate = uni_results_2,
    candidates = candidate_markers_2,
    lasso = lasso_2,
    selected_markers = selected_markers_2,
    model_list = model_list_2,
    roc = result_2
  )
  save_result(binary_results, "Binary_Screening_Results.rds")
}

# 7. 汇总二分类AUC --------------------------------------------
binary_auc_summary <- dplyr::bind_rows(lapply(names(binary_results), function(nm) {
  auc_tmp <- binary_results[[nm]]$roc$auc_results
  if (nrow(auc_tmp) == 0) return(NULL)
  auc_tmp$Comparison <- nm
  auc_tmp[, c("Comparison", "Model", "AUC")]
}))
cat("\n================ 二分类AUC汇总 ================\n")
print(binary_auc_summary)
save_result(binary_auc_summary, "Binary_All_Comparisons_ROC_AUC.rds")
write_result(binary_auc_summary, "Binary_All_Comparisons_ROC_AUC.csv")

# 8. 森林图：单变量与多变量OR及95% CI --------------------------
cat("\n\n##################################################\n")
cat("########## 森林图：OR及95% CI ##########\n")
cat("##################################################\n")

# 8.1 通用森林图绘制函数
plot_forest <- function(forest_data, title, facet_variable = NULL) {
  forest_data <- forest_data %>% dplyr::filter(
    is.finite(OR), is.finite(Lower_CI), is.finite(Upper_CI),
    OR > 0, Lower_CI > 0, Upper_CI > 0
  )
  if (nrow(forest_data) == 0) {
    warning("没有可绘制的森林图数据：", title)
    return(invisible(NULL))
  }
  
  marker_order <- rev(unique(forest_data$Marker))
  forest_data$Marker <- factor(forest_data$Marker, levels = marker_order)
  
  # 兼容所有模型名称的因子排序
  model_levels <- c("Univariate", "Model_Combined", "Model_miRNA", 
                    "Model_Other_Markers", "Model_LASSO_Selected")
  forest_data$Model <- factor(forest_data$Model, levels = intersect(model_levels, unique(forest_data$Model)))
  
  old_par <- par(no.readonly = TRUE)
  par(mar = c(3.5, 3.5, 4, 2) + 0.1,
      mgp = c(2.2, 0.3, 0))
  
  p <- ggplot(forest_data, aes(x = OR, y = Marker, color = Model, shape = Model)) +
    # 将 dodge width 调大至 0.8，让点和线错开
    geom_segment(aes(x = Lower_CI, xend = Upper_CI, y = Marker, yend = Marker),
                 position = position_dodge(width = 0.8), linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.8), size = 2.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    scale_x_log10(
      labels = scales::label_number()
    ) +
    coord_cartesian(xlim = c(0.01, 1000)) +
    # 【修改点】调整图例的标签名称，避免过长导致显示不全
    scale_color_manual(
      values = c(
        "Univariate" = "grey50",
        "Model_Combined" = "#D9534F",
        "Model_miRNA" = "#2878B5",
        "Model_Other_Markers" = "#28A745",
        "Model_LASSO_Selected" = "#8E44AD"
      ),
      labels = c(
        "Univariate" = "Univariate",
        "Model_Combined" = "Combined",
        "Model_miRNA" = "miRNA",
        "Model_Other_Markers" = "Others",
        "Model_LASSO_Selected" = "LASSO"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Univariate" = 16,
        "Model_Combined" = 17,
        "Model_miRNA" = 15,
        "Model_Other_Markers" = 18,
        "Model_LASSO_Selected" = 8
      ),
      labels = c(
        "Univariate" = "Univariate",
        "Model_Combined" = "Combined",
        "Model_miRNA" = "miRNA",
        "Model_Other_Markers" = "Others",
        "Model_LASSO_Selected" = "LASSO"
      )
    ) +
    labs(title = title, x = "Odds Ratio (95% CI, log scale)",
         y = NULL, color = NULL, shape = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.position = "bottom", 
      legend.text = element_text(size = 6),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    # 【修改点】强制图例排成一行，防止因空间不足被截断
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    )
  
  if (!is.null(facet_variable)) {
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_variable)),
                        scales = "free_y")
  }
  par(old_par)
  
  print(p)
  invisible(p)
}

# 8.2 提取多项Logistic回归OR及95% CI
extract_multinom_forest <- function(fit, markers, model_type) {
  coef_matrix <- as.matrix(summary(fit)$coefficients)
  se_matrix <- as.matrix(summary(fit)$standard.errors)
  results <- list()
  for (class_name in rownames(coef_matrix)) {
    for (m in markers) {
      if (!m %in% colnames(coef_matrix) || !m %in% colnames(se_matrix)) next
      beta <- coef_matrix[class_name, m]
      se <- se_matrix[class_name, m]
      if (!is.finite(beta) || !is.finite(se) || se <= 0) next
      z_value <- beta / se
      results[[length(results) + 1]] <- data.frame(
        Contrast = paste0(class_name, " vs 0"), Marker = m,
        Model = model_type, Beta = beta, SE = se, OR = exp(beta),
        Lower_CI = exp(beta - 1.96 * se),
        Upper_CI = exp(beta + 1.96 * se),
        P_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
      )
    }
  }
  dplyr::bind_rows(results)
}

# 8.3 四分类单变量模型
forest_uni_4 <- dplyr::bind_rows(lapply(markers, function(m) {
  fit_uni <- tryCatch(
    suppressWarnings(nnet::multinom(reformulate(m, response = "Y4"),
                                    data = df_model, trace = FALSE, maxit = 1000)),
    error = function(e) { NULL }
  )
  if (is.null(fit_uni)) return(NULL)
  extract_multinom_forest(fit_uni, m, "Univariate")
}))

# 8.4 四分类多变量模型（按 Model 组合提取）
forest_multi_4_list <- list()
for (model_name in names(model_list_4)) {
  vars_in_model <- model_list_4[[model_name]]
  formula_m <- reformulate(c(covariates, vars_in_model), response = "Y4")
  
  fit_m <- tryCatch(
    suppressWarnings(nnet::multinom(formula_m, data = df_model, trace = FALSE, maxit = 1000)),
    error = function(e) NULL
  )
  
  if (!is.null(fit_m)) {
    res <- extract_multinom_forest(fit_m, vars_in_model, model_name)
    forest_multi_4_list[[model_name]] <- res
  }
}

# 合并单变量与多变量数据（使用 do.call 防止 bind_rows 报错）
forest_results_4 <- do.call(rbind, c(list(forest_uni_4), forest_multi_4_list))

if (nrow(forest_results_4) > 0) {
  plot_forest(forest_results_4, "Four-Class Multinomial (Model Groups)", facet_variable = "Contrast")
}

# 8.5 提取二分类Logistic回归OR及95% CI
extract_binary_forest <- function(fit, markers, model_type) {
  coef_table <- summary(fit)$coefficients
  results <- list()
  for (m in markers) {
    if (!m %in% rownames(coef_table)) next
    beta <- coef_table[m, "Estimate"]
    se <- coef_table[m, "Std. Error"]
    if (!is.finite(beta) || !is.finite(se) || se <= 0) next
    z_value <- beta / se
    results[[length(results) + 1]] <- data.frame(
      Marker = m, Model = model_type, Beta = beta, SE = se,
      OR = exp(beta), Lower_CI = exp(beta - 1.96 * se),
      Upper_CI = exp(beta + 1.96 * se),
      P_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
    )
  }
  dplyr::bind_rows(results)
}

# 8.6 三种二分类比较的森林图
binary_forest_results <- list()
for (comparison_name in names(binary_results)) {
  cat("\n================ 森林图：", comparison_name, " ================\n", sep = "")
  df_binary <- binary_results[[comparison_name]]$data
  
  # 8.6.1 单变量模型
  forest_uni_2 <- dplyr::bind_rows(lapply(markers, function(m) {
    fit_uni <- tryCatch(
      suppressWarnings(glm(reformulate(m, response = "Y_binary"),
                           data = df_binary, family = binomial)),
      error = function(e) NULL
    )
    if (is.null(fit_uni)) return(NULL)
    extract_binary_forest(fit_uni, m, "Univariate")
  }))
  
  # 8.6.2 多变量模型组合（按 Model 提取）
  model_list_2 <- build_model_list(markers, miRNA_markers, other_markers, 
                                   binary_results[[comparison_name]]$selected_markers)
  
  forest_multi_2_list <- list()
  for (model_name in names(model_list_2)) {
    vars_in_model <- model_list_2[[model_name]]
    formula_m2 <- reformulate(c(covariates, vars_in_model), response = "Y_binary")
    
    fit_m2 <- tryCatch(
      suppressWarnings(glm(formula_m2, data = df_binary, family = binomial)),
      error = function(e) NULL
    )
    
    if (!is.null(fit_m2)) {
      forest_multi_2_list[[model_name]] <- extract_binary_forest(fit_m2, vars_in_model, model_name)
    }
  }
  
  # 8.6.3 合并数据并画图
  forest_combined_2 <- do.call(rbind, c(list(forest_uni_2), forest_multi_2_list))
  
  binary_forest_results[[comparison_name]] <- forest_combined_2
  cat("\n森林图数据：\n"); print(forest_combined_2)
  save_result(binary_forest_results, "Binary_Forest_Plot_Results.rds")
  write_result(forest_combined_2, paste0(comparison_name, "_Forest_Plot_Data.csv"))
  
  if (nrow(forest_combined_2) > 0) {
    plot_forest(forest_combined_2,
                paste0(comparison_name, " Logistic Regression"))
  }
}

# 8.7 汇总二分类森林图结果
binary_forest_summary <- dplyr::bind_rows(lapply(names(binary_forest_results), function(nm) {
  temp <- binary_forest_results[[nm]]
  if (nrow(temp) == 0) return(NULL)
  temp$Comparison <- nm
  temp
}))
cat("\n================ 二分类森林图汇总 ================\n")
print(binary_forest_summary)
save_result(binary_forest_summary, "Binary_All_Comparisons_Forest_Plot_Data.rds")
write_result(binary_forest_summary, "Binary_All_Comparisons_Forest_Plot_Data.csv")

# 9. 保存最终分析结果 ------------------------------------------
screening_objects <- list(
  df_model = df_model,
  markers = markers,
  covariates = covariates,
  miRNA_markers = miRNA_markers,
  other_markers = other_markers,
  group_levels = group_levels,
  binary_comparisons = binary_comparisons,
  imputation_id = 1,
  uni_results_4 = uni_results_4,
  candidate_markers_4 = candidate_markers_4,
  lasso_4 = lasso_4,
  selected_markers_4 = selected_markers_4,
  model_list_4 = model_list_4,
  result_4 = result_4,
  auc_results_4 = auc_results_4,
  binary_results = binary_results,
  binary_auc_summary = binary_auc_summary,
  forest_results_4 = forest_results_4,
  binary_forest_results = binary_forest_results,
  binary_forest_summary = binary_forest_summary
)
save_result(screening_objects, "Screening_Results.rds")

cat("\n============================================\n")
cat("变量筛选、初步建模及森林图分析全部完成。\n")
cat("============================================\n")
cat("\n四分类样本分布：\n"); print(table(df_model$Y4))
cat("\n四分类LASSO选中的标志物：\n"); print(selected_markers_4)
cat("\n四分类AUC：\n"); print(auc_results_4)
cat("\n二分类比较结果：\n")
for (comparison_name in names(binary_results)) {
  cat("\n", comparison_name, "\n", sep = "")
  cat("样本量：", nrow(binary_results[[comparison_name]]$data), "\n")
  print(table(binary_results[[comparison_name]]$data$Y_binary))
  cat("LASSO选中的标志物：\n"); print(binary_results[[comparison_name]]$selected_markers)
  cat("AUC：\n"); print(binary_results[[comparison_name]]$roc$auc_results)
}
cat("\n最终结果文件：\n", file.path(output_dir, "Screening_Results.rds"), "\n")
