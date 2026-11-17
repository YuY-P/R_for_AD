# ============================================================
# Norm_Test_Zscore.R
#
# Purpose:
#   ...
#     1. Z-score standardization
#     2. Statistic analysis
#     3. Q-Q plot
#
# Input:
#   MI_data.rds (stacked completed datasets; one row per subject
#   and imputation_id)
#
# Output:
#   ...
#
# Note:
#   This script uses imputation_id = 1 for exploratory analyses.
#   It does not pool inferential estimates across imputations.
# ============================================================

rm(list = ls())
gc()

library(e1071)
library(ggplot2)

df_final <- readRDS("D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/MI_data.rds")

markers <- setdiff(names(df_final)[sapply(df_final, is.numeric)], 
                   c("imputation_id", "Sex", "CCI"))

# ================= 1. 生成 Z-score 标准化列 =================
for (m in markers) {
  # 按照 imputation_id 分组进行 Z-score 标准化
  # ave() 函数确保每个插补数据集内部独立计算均值和标准差
  df_final[[paste0("z_", m)]] <- ave(
    df_final[[m]], 
    df_final$imputation_id, 
    FUN = function(x) as.numeric(scale(x))
  )
}

# 提取刚刚生成的 z-score 列
z_markers <- paste0("z_", markers)

# ================= 2. 统计检验部分 =================
results_z <- data.frame()
for (i in 1:20) {
  df_m <- df_final[df_final$imputation_id == i, ]
  for (marker in z_markers) {
    x <- df_m[[marker]]
    # 过滤掉因原本全为同一个值导致的 NaN，或者原有数据中可能包含的 NA
    x <- x[is.finite(x)] 
    
    if(length(x) > 3) { # 确保样本量足够进行 Shapiro 检验
      sw_test <- shapiro.test(x)
      skew_val <- e1071::skewness(x)
      
      results_z <- rbind(results_z, data.frame(
        Imputation_ID = i,
        Marker = marker,
        P_value = sw_test$p.value,
        Skewness = skew_val,
        Normal = ifelse(sw_test$p.value > 0.05, "Yes", "No")
      ))
    }
  }
}

print(results_z)

saveRDS(
  df_final, 
  file = "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/Z_Norm_Markers.rds"
)

# 导出为 Excel
write.csv(results_z, 
          "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/Normality_Test_Results_Z.csv", 
          row.names = FALSE)


# ================= 3. Q-Q plot 绘图部分 (基于 Z-score 数据) =================

qq_z_main  <- list()  

# 设置主图的横坐标范围 (依然保持 -3.5 到 3.5)
main_xlim <- c(-3.5, 3.5)

for (marker in z_markers) {
  
  cat("正在绘制 Z-score 转换后的 Q-Q 图：", marker, "\n")
  
  # ---------- 动态计算 Y 轴范围 ----------
  raw_data <- df_final[[marker]]
  raw_data <- raw_data[is.finite(raw_data)] # 剔除 NA 和 NaN
  
  # 计算理论分位数（X轴）
  theoretical_q <- qnorm(ppoints(length(raw_data)))
  
  # 找到理论分位数落在 [-3.5, 3.5] 区间内的样本真实值，用来确定 Y 轴范围
  y_in_range <- raw_data[theoretical_q >= main_xlim[1] & theoretical_q <= main_xlim[2]]
  
  # 设定 Y 轴范围，并留出 5% 的上下边距
  y_min <- min(y_in_range, na.rm = TRUE)
  y_max <- max(y_in_range, na.rm = TRUE)
  y_padding <- (y_max - y_min) * 0.05
  main_ylim <- c(y_min - y_padding, y_max + y_padding)
  
  # ---------- A. 绘制主图 ----------
  p_main <- ggplot(
    df_final,
    aes(sample = .data[[marker]])
  ) +
    stat_qq(
      aes(color = factor(imputation_id),
          group = imputation_id),
      alpha = 0.6,
      size = 1.2
    ) +
    stat_qq_line(
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # 同时限制 X 和 Y 轴，使点铺满画布，恒等线对角
    coord_cartesian(
      xlim = main_xlim,
      ylim = main_ylim
    ) +
    theme_bw() +
    labs(
      title = paste("Q-Q Plot for Z-scored", sub("^z_", "", marker)), # 自动去掉 z_ 前缀
      x = "Theoretical Quantiles (Normal)",
      y = "Sample Quantiles (Z-score)", # Y 轴标签改为 Z-score
      color = "Imputation ID"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(
        ncol = 2,
        override.aes = list(alpha = 1, size = 2)
      )
    )
  
  # ---------- B. 保存并展示主图 ----------
  qq_z_main[[marker]] <- p_main
  
  # 显示当前 Marker 的最终图
  print(p_main)
}

# 查看保存的 Marker 名称
names(qq_z_main)

# 修改坐标轴范围示例
# qq_z_main[["z_Abeta42C_LODcut"]]$coordinates$limits$y <- c(-3, 3)
# print(qq_z_main[["z_Abeta42C_LODcut"]])