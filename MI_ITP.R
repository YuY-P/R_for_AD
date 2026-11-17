# ============================================================
# MI_ITP.R
#
# Purpose:
#   Interpolate for any left-censored and missing digits.
#
# Input:
#   Sheet "Merge" from Cohort_Data_Collect.xlsx
#
# Output:
#   Intermediate RDS files for downstream analysis
#     MI_data.rds
#   CSV files containing screening and model evaluation results
#     Merged_MI_Imputed.csv
# 
# Note:
#   ...
# ============================================================


rm(list = ls())
gc()

library(VGAM)
library(readxl)
library(truncnorm)
library(MASS)

# ==========================================
# 1. 读取数据
# ==========================================

mydata <- read_excel(
  "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/Cohort_Data_Collect.xlsx",
  sheet = "Merge"
)

id_col <- "No"
group_col <- "Tag_No"

mydata[[group_col]] <- as.factor(mydata[[group_col]])

if (anyNA(mydata[[group_col]])) {
  stop("分组变量存在缺失值，请先处理。")
}

# ==========================================
# 2. 自动识别需要插补的指标列
# ==========================================

# 2.1 识别数值型二分类变量
# 根据非缺失值判断：恰好存在两个不同取值的数值型变量
binary_cols <- names(mydata)[vapply(mydata, function(x) {
  
  if (!is.numeric(x)) {
    return(FALSE)
  }
  
  unique_values <- unique(x[!is.na(x)])
  
  length(unique_values) == 2
  
}, logical(1))]

# 排除ID、分组变量和二分类变量、离散变量
exclude_cols <- c(
  id_col,
  group_col,
  binary_cols,
  "Age",
  "CCI"
)

# 2.2 自动识别需要插补的连续数值型指标
cols_to_impute <- names(mydata)[vapply(mydata, function(x) {
  
  is.numeric(x) && any(is.na(x) | x == 0, na.rm = TRUE)
  
}, logical(1))]

cols_to_impute <- setdiff(
  cols_to_impute,
  exclude_cols
)

# 2.3 输出识别结果
cat("\n识别出的数值型二分类变量（已排除）：\n")
print(binary_cols)

cat("\n需要进行连续变量插补的指标：\n")
print(cols_to_impute)

cat("\n未纳入本次插补的其他变量：\n")
print(setdiff(
  names(mydata),
  c(cols_to_impute, binary_cols, id_col, group_col)
))

if (length(cols_to_impute) == 0) {
  stop("未检测到需要插补的连续数值型指标。")
}

# ==========================================
# 3. 多重插补
# ==========================================

M <- 20
imputed_datasets <- vector("list", M)

set.seed(2026)

for (m in seq_len(M)) {
  
  cat("\n正在生成第", m, "个插补数据集...\n")
  
  df_m <- mydata
  df_m$imputation_id <- m
  
  for (col in cols_to_impute) {
    
    cat("  处理指标：", col, "\n")
    
    val <- mydata[[col]]
    
    # 识别左删失和普通缺失
    is_censored <- !is.na(val) & val == 0
    is_missing  <- is.na(val)
    
    # 检查数据
    if (any(val < 0, na.rm = TRUE)) {
      stop(paste("指标", col, "仍存在负值，请检查数据。"))
    }
    
    # 确定经验性左删失上限
    pos_vals <- val[val > 0 & !is.na(val)]
    
    if (length(pos_vals) == 0) {
      stop(paste("指标", col, "没有正浓度值，无法确定插补上限。"))
    }
    
    # 略低于最小正浓度，避免最小正值被VGAM误识别为删失值
    censor_val <- min(pos_vals) * (1 - 1e-8)
    log_censor <- log(censor_val)
    
    # --------------------------------------
    # A. 构建Tobit模型数据
    # --------------------------------------
    
    df_temp <- data.frame(
      group = mydata[[group_col]],
      log_value = rep(NA_real_, nrow(mydata))
    )
    
    # 检出值：保留实际浓度的对数
    det_idx <- which(!is_censored & !is_missing)
    
    df_temp$log_value[det_idx] <- log(val[det_idx])
    
    # 左删失值：赋予统一的删失边界
    cens_idx <- which(is_censored)
    
    df_temp$log_value[cens_idx] <- log_censor
    
    # 普通缺失值保持NA，不参与Tobit参数估计
    
    # --------------------------------------
    # B. 拟合
    # --------------------------------------
    
    # 根据是否存在左删失选择模型
    family_use <- if (any(is_censored)) {
      tobit(
        Lower = log_censor,
        zero = 2
      )
    } else {
      uninormal(zero = 2)
    }
    
    fit <- tryCatch(
      vglm(
        log_value ~ group,
        family = family_use,
        data = df_temp,
        na.action = na.omit,
        trace = FALSE
      ),
      error = function(e) {
        stop(paste(
          "指标", col, "模型拟合失败：", e$message
        ))
      }
    )
    
    # --------------------------------------
    # C. 抽取模型参数并引入参数不确定性
    # --------------------------------------
    
    # colon=TRUE：明确区分两个线性预测量的系数
    coef_est <- coef(fit, colon = TRUE)
    
    # 默认系数顺序与上面的系数数值顺序一致
    coef_default <- coef(fit)
    
    if (!isTRUE(all.equal(
      as.numeric(coef_est),
      as.numeric(coef_default)
    ))) {
      stop(paste("指标", col, "的模型系数顺序不一致。"))
    }
    
    vcov_est <- as.matrix(vcov(fit))
    
    if (any(!is.finite(vcov_est))) {
      stop(paste("指标", col, "的协方差矩阵存在非有限值。"))
    }
    
    # 协方差矩阵对称化
    vcov_est <- (vcov_est + t(vcov_est)) / 2
    
    # 检查协方差矩阵是否存在明显的非正定问题
    eig_values <- eigen(
      vcov_est,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    if (min(eig_values) < -1e-8 * max(1, max(abs(eig_values)))) {
      stop(paste("指标", col, "的协方差矩阵非正定，请检查模型。"))
    }
    
    # 按照系数向量的顺序进行多元正态抽样
    coef_draw <- MASS::mvrnorm(
      n = 1,
      mu = as.numeric(coef_default),
      Sigma = vcov_est
    )
    
    names(coef_draw) <- names(coef_est)
    
    # 均值模型系数：预测量1
    mu_idx <- grepl(":1$", names(coef_draw))
    
    # 标准差模型系数：预测量2
    sd_idx <- grepl(":2$", names(coef_draw))
    
    if (!any(mu_idx) || sum(sd_idx) != 1) {
      stop(paste("指标", col, "的Tobit模型系数提取异常。"))
    }
    
    beta_mu <- coef_draw[mu_idx]
    
    names(beta_mu) <- sub(
      ":1$", "",
      names(beta_mu)
    )
    
    # 构建分组设计矩阵
    X_new <- model.matrix(
      ~ group,
      data = data.frame(group = mydata[[group_col]])
    )
    
    # 按设计矩阵顺序排列均值模型系数
    if (!setequal(names(beta_mu), colnames(X_new))) {
      stop(paste("指标", col, "的均值模型系数与设计矩阵不匹配。"))
    }
    
    beta_mu <- beta_mu[colnames(X_new)]
    
    mu_log <- as.numeric(X_new %*% beta_mu)
    
    # 标准差：第二个线性预测量为log(sigma)
    sigma_m <- exp(as.numeric(coef_draw[sd_idx]))
    
    if (!is.finite(sigma_m) || sigma_m <= 0) {
      stop(paste("指标", col, "抽样得到的标准差无效。"))
    }
    
    # --------------------------------------
    # D. 插补左删失值（0）
    # --------------------------------------
    
    if (length(cens_idx) > 0) {
      
      df_m[[col]][cens_idx] <- exp(
        truncnorm::rtruncnorm(
          n = length(cens_idx),
          a = -Inf,
          b = log_censor,
          mean = mu_log[cens_idx],
          sd = sigma_m
        )
      )
    }
    
    # --------------------------------------
    # E. 插补普通缺失值（NA）—— 已移除上限限制
    # --------------------------------------
    
    miss_idx <- which(is_missing)
    
    if (length(miss_idx) > 0) {
      
      # 直接从无界正态分布抽样，不做任何上限截断
      df_m[[col]][miss_idx] <- exp(
        rnorm(
          n = length(miss_idx),
          mean = mu_log[miss_idx],
          sd = sigma_m
        )
      )
      
    }
    
  }
  
  imputed_datasets[[m]] <- df_m
}

# ==========================================
# 4. 合并与导出
# ==========================================

df_final <- do.call(rbind, imputed_datasets)

cat("\n多重插补完成！\n")
cat("原始样本数：", nrow(mydata), "\n")
cat("插补数据集数量：", M, "\n")
cat("合并后总行数：", nrow(df_final), "\n")

# 检查是否仍存在NA或0
cat("\n各指标剩余NA数量：\n")
print(colSums(is.na(df_final[, cols_to_impute])))

cat("\n各指标剩余0值数量：\n")
print(colSums(df_final[, cols_to_impute] == 0, na.rm = TRUE))


saveRDS(
  df_final, 
  file = "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/MI_data.rds"
)

write.csv(
  df_final,
  "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/Merged_MI_Imputed.csv",
  row.names = FALSE
)

cat("\n数据已保存至 Merged_MI_Imputed.csv\n")
