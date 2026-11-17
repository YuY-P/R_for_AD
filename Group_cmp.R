# ============================================================
# Group_cmp.R
#
# Purpose:
#   Exploratory comparison of biomarker concentrations across
#   diagnostic groups using one completed imputed dataset.
#     1-2. Input and initiation
#     3. Get biomarkers
#     4. log transformation
#     5. Baseline table
#     6. Min-max normalized heat-map
#     7. Boxplot
#       5-Class: Kruskal-Wallis + Dunn
#       2-Class: Wilcoxon 秩和检验
#     8. 效应量汇总
#       5-Class: Kruskal-Wallis epsilon-squared
#       2-Class: 二组比较报告 Wilcoxon r
#     9. Spearman correlation heat-map
#     10. Output
# 
# Input:
#   MI_data.rds (stacked completed datasets; one row per subject
#   and imputation_id)
#
# Output:
#   ...
#     Table1_Baseline.csv
#     Group_Cmp_Heatmap
#
# Note:
#   This script uses imputation_id = 1 for exploratory analyses.
#   It does not pool inferential estimates across imputations.
# ============================================================

rm(list = ls())
gc()
graphics.off()

# 1. 加载程序包 ------------------------------------------------

library(dplyr)
library(ggplot2)
library(tableone)
library(ggpubr)
library(rstatix)
library(tidyr)
library(Hmisc)
library(corrplot)
library(scales)
library(ggdendro)
library(patchwork)


# 2. 数据导入与参数设置 ----------------------------------------

input_file <- paste0(
  "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS/",
  "MI_data.rds"
)
output_dir <- "D:/0 BME/0006 PEC+Neuro/EXP_DATA/ANALYSIS"

id_var <- "No"
group_var <- "Tag_No"
age_var <- "Age"
sex_var <- "Sex"
cci_var <- "CCI"
imp_var <- "imputation_id"
selected_imputation <- 1

# 原始五分类及三种二分类比较
original_group_levels <- c("0", "1", "2", "3", "4")
binary_comparisons <- list(
  "0_vs_1plus2plus3" = c("0", "1", "2", "3"),
  "0_vs_3plus4" = c("0", "3", "4"),
  "0_vs_1plus2plus3plus4" = c("0", "1", "2", "3", "4")
)

# 显著性水平及多重比较校正方法
alpha <- 0.05
p_adjust_method <- "BH"

# 统一绘图边距：上、右、下、左，单位为 pt
common_plot_theme <- theme(
  plot.margin = margin(t = 20, r = 15, b = 10, l = 20, unit = "pt")
)

if (!file.exists(input_file)) {
  stop("未找到输入文件：", input_file)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

save_checkpoint <- function(object, filename) {
  saveRDS(object, file.path(output_dir, filename))
  cat("RDS saved:", filename, "\n")
}

export_csv <- function(object, filename, row_names = FALSE) {
  write.csv(
    object,
    file.path(output_dir, filename),
    row.names = row_names,
    fileEncoding = "UTF-8"
  )
  cat("CSV saved:", filename, "\n")
}

# 3. 读取数据并识别标志物 --------------------------------------

df_all <- readRDS(input_file)

fixed_vars <- c(id_var, group_var, age_var, sex_var, cci_var, imp_var)
missing_vars <- setdiff(fixed_vars, names(df_all))
if (length(missing_vars) > 0) {
  stop("数据中缺少必要变量：", paste(missing_vars, collapse = ", "))
}

if (!selected_imputation %in% df_all[[imp_var]]) {
  stop("数据中不存在 imputation_id = ", selected_imputation)
}

df_raw <- df_all %>%
  filter(.data[[imp_var]] == selected_imputation) %>%
  droplevels()

# 统一分组标签，确保绘图顺序为 0、1、2、3、4
group_values <- as.character(df_raw[[group_var]])
unknown_groups <- setdiff(unique(na.omit(group_values)), original_group_levels)
if (length(unknown_groups) > 0) {
  stop("发现未定义的疾病分组：", paste(unknown_groups, collapse = ", "))
}
df_raw[[group_var]] <- factor(group_values, levels = original_group_levels)

observed_groups <- unique(as.character(na.omit(df_raw[[group_var]])))
if (length(observed_groups) < 2) {
  stop("有效疾病分组不足两个。")
}

markers <- setdiff(
  names(df_raw)[vapply(df_raw, is.numeric, logical(1))],
  fixed_vars
)
if (length(markers) == 0) {
  stop("未识别到数值型标志物，请检查数据列。")
}

cat("\n本次分析使用 imputation_id =", selected_imputation, "\n")
cat("识别到的标志物：\n")
print(markers)
cat("各组样本量：\n")
print(table(df_raw[[group_var]], useNA = "ifany"))

# 标志物浓度应为非负值；负值不进行静默修正
negative_counts <- vapply(
  markers,
  function(v) sum(df_raw[[v]] < 0, na.rm = TRUE),
  integer(1)
)
if (any(negative_counts > 0)) {
  stop(
    "检测到负浓度值。请先核对浓度换算及数据处理；本脚本不会自动修改原始数据。\n",
    paste(names(negative_counts)[negative_counts > 0], collapse = ", ")
  )
}

# 4. 构建对数转换数据 ------------------------------------------

# 对含 0 值的标志物，使用该标志物最小正值的一半作为 log 转换偏移量。
# 该偏移量仅用于 log 转换，不会改写原始浓度数据。
log_offsets <- data.frame(
  Marker = markers,
  Offset = NA_real_,
  stringsAsFactors = FALSE
)

df_log <- df_raw
for (i in seq_along(markers)) {
  marker <- markers[i]
  x <- df_raw[[marker]]
  positive_values <- x[is.finite(x) & x > 0]

  if (length(positive_values) == 0) {
    log_offsets$Offset[i] <- NA_real_
    df_log[[marker]] <- NA_real_
  } else {
    offset <- min(positive_values) / 2
    log_offsets$Offset[i] <- offset
    df_log[[marker]] <- log10(x + offset)
  }
}

save_checkpoint(
  list(
    imputation_id = selected_imputation,
    df_raw = df_raw,
    df_log = df_log,
    markers = markers,
    log_offsets = log_offsets,
    group_levels = original_group_levels,
    binary_comparisons = binary_comparisons
  ),
  "Group_Cmp_Data.rds"
)
export_csv(log_offsets, "Marker_Log_Offsets.csv")

# 5. 基线特征表（Table 1） -------------------------------------

continuous_vars <- age_var
categorical_vars <- c(sex_var, cci_var)
table_vars <- c(markers, continuous_vars, categorical_vars)

df_table <- df_raw
for (v in categorical_vars) {
  df_table[[v]] <- as.factor(df_table[[v]])
}

table1 <- CreateTableOne(
  vars = table_vars,
  strata = group_var,
  data = df_table,
  factorVars = categorical_vars,
  addOverall = TRUE,
  test = TRUE
)

table1_print <- print(
  table1,
  nonnormal = markers,
  showAllLevels = TRUE,
  smd = TRUE,
  catDigits = 2,
  contDigits = 2,
  pDigits = 3
)

table1_csv <- print(
  table1,
  nonnormal = markers,
  showAllLevels = TRUE,
  smd = TRUE,
  catDigits = 2,
  contDigits = 2,
  pDigits = 3,
  printToggle = FALSE,
  quote = FALSE,
  noSpaces = TRUE
)

save_checkpoint(
  list(table1 = table1, table1_print = table1_print, table1_csv = table1_csv),
  "Group_Cmp_Table1.rds"
)
export_csv(table1_csv, "Table1_Baseline.csv", row_names = TRUE)

# 6. 病例-标志物热图：log 转换后进行逐标志物归一化 ------------

heatmap_data <- df_log %>%
  dplyr::select(dplyr::all_of(c(id_var, group_var, markers))) %>%
  mutate(
    CaseID = as.character(.data[[id_var]]),
    Group = factor(as.character(.data[[group_var]]), levels = original_group_levels)
  ) %>%
  arrange(Group, CaseID)

if (anyDuplicated(heatmap_data$CaseID)) {
  stop("当前插补数据集中病例编号不唯一，无法绘制病例热图。")
}

heatmap_long <- heatmap_data %>%
  dplyr::select(CaseID, Group, dplyr::all_of(markers)) %>%
  pivot_longer(
    cols = all_of(markers),
    names_to = "Marker",
    values_to = "Log_Value"
  ) %>%
  mutate(Marker = factor(Marker, levels = markers)) %>%
  group_by(Marker) %>%
  mutate(
    Value_norm = {
      valid <- is.finite(Log_Value)
      if (!any(valid)) {
        rep(NA_real_, length(Log_Value))
      } else {
        value_min <- min(Log_Value[valid])
        value_max <- max(Log_Value[valid])
        if (value_max == value_min) {
          ifelse(valid, 0.5, NA_real_)
        } else {
          ifelse(valid, (Log_Value - value_min) / (value_max - value_min), NA_real_)
        }
      }
    }
  ) %>%
  ungroup()

heatmap_long$CaseID <- factor(
  heatmap_long$CaseID,
  levels = rev(unique(heatmap_data$CaseID))
)

heatmap_colors <- c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026")

p_heatmap <- ggplot(
  heatmap_long,
  aes(x = Marker, y = CaseID, fill = Value_norm)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(
    Group ~ ., scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_fill_gradientn(
    colours = heatmap_colors,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    na.value = "grey90",
    name = "Scaled\nlog(Conc.)",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth = unit(0.5, "cm"),
      title.position = "top"
    )
  ) +
  labs(
    x = "Biomarker",
    y = "Case"
  ) +
  theme_bw() +
  common_plot_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0, size = 8),
    strip.placement = "outside",
    panel.spacing.y = grid::unit(0.1, "lines"),
    legend.position = "right"
  )

print(p_heatmap)

save_checkpoint(
  list(heatmap_data = heatmap_data, heatmap_long = heatmap_long, plot = p_heatmap),
  "Group_Cmp_Heatmap.rds"
)
export_csv(heatmap_long, "All_Cases_Marker_Heatmap_Data.csv")

# 7. 组间比较与箱线图 ------------------------------------------

# 箱线图使用原始浓度；pseudo-log 坐标可保留 0 值。
# 统计检验及效应量使用 log 转换浓度，秩检验结果对严格单调变换不变。

make_pseudo_log_sigma <- function(x) {
  positive <- x[is.finite(x) & x > 0]
  if (length(positive) == 0) return(1)
  max(min(positive) / 10, .Machine$double.eps)
}

make_concentration_scale <- function(x) {
  if (any(x == 0, na.rm = TRUE)) {
    scale_y_continuous(
      trans = scales::pseudo_log_trans(
        base = 10,
        sigma = make_pseudo_log_sigma(x)
      ),
      labels = scales::label_number()
    )
  } else {
    scale_y_log10(labels = scales::label_number())
  }
}

make_y_positions <- function(data, n) {
  ymax <- max(data$Value[is.finite(data$Value)], na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) {
    return(seq_len(n) * 0.1)
  }
  ymax * (1.25 ^ seq_len(n))
}

# 7.1 五组总体比较：Kruskal-Wallis + Dunn 两两比较
multi_group_results <- vector("list", length(markers))
names(multi_group_results) <- markers
multi_group_dunn <- vector("list", length(markers))
names(multi_group_dunn) <- markers
multi_group_plots <- vector("list", length(markers))
names(multi_group_plots) <- markers

for (marker in markers) {
  plot_data <- data.frame(
    Group = factor(as.character(df_raw[[group_var]]), levels = original_group_levels),
    Value = df_raw[[marker]]
  ) %>%
    filter(!is.na(Group), is.finite(Value)) %>%
    droplevels()

  test_data <- data.frame(
    Group = factor(as.character(df_log[[group_var]]), levels = original_group_levels),
    Log_Value = df_log[[marker]]
  ) %>%
    filter(!is.na(Group), is.finite(Log_Value)) %>%
    droplevels()

  if (n_distinct(test_data$Group) < 2) {
    multi_group_results[[marker]] <- data.frame(
      Marker = marker, N = nrow(test_data), Test = "Kruskal-Wallis",
      Effect_Size_Type = "Kruskal-Wallis epsilon-squared",
      Effect_Size = NA_real_, P_value = NA_real_
    )
  } else {

  kw <- tryCatch(
    kruskal.test(Log_Value ~ Group, data = test_data),
    error = function(e) NULL
  )

  eff <- tryCatch(
    rstatix::kruskal_effsize(test_data, Log_Value ~ Group),
    error = function(e) NULL
  )

  p_value <- if (is.null(kw)) NA_real_ else as.numeric(kw$p.value)
  effect_value <- if (!is.null(eff) && nrow(eff) > 0) {
    as.numeric(eff$effsize[1])
  } else {
    NA_real_
  }

  multi_group_results[[marker]] <- data.frame(
    Marker = marker,
    N = nrow(test_data),
    Test = "Kruskal-Wallis",
    Effect_Size_Type = "Kruskal-Wallis epsilon-squared",
    Effect_Size = effect_value,
    P_value = p_value,
    stringsAsFactors = FALSE
  )

  dunn <- tryCatch(
    rstatix::dunn_test(
      test_data,
      Log_Value ~ Group,
      p.adjust.method = p_adjust_method
    ),
    error = function(e) NULL
  )

  if (!is.null(dunn) && nrow(dunn) > 0) {
    dunn$Marker <- marker
    multi_group_dunn[[marker]] <- dunn
  } else {
    multi_group_dunn[[marker]] <- data.frame()
  }

  }

  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.62) +
    geom_jitter(width = 0.14, alpha = 1, size = 1, aes(color = Group)) +
    scale_x_discrete(drop = FALSE) +
    make_concentration_scale(plot_data$Value) +
    labs(
      title = marker,
      x = "Diagnostic group",
      y = "Conc."
    ) +
    theme_bw() +
    common_plot_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  multi_group_plots[[marker]] <- p
}

multi_group_summary <- bind_rows(multi_group_results) %>%
  mutate(P_adjusted_across_markers = p.adjust(P_value, method = p_adjust_method))
multi_group_dunn_summary <- bind_rows(multi_group_dunn)
if (ncol(multi_group_dunn_summary) == 0) {
  multi_group_dunn_summary <- data.frame(
    Marker = character(),
    group1 = character(),
    group2 = character(),
    p = numeric(),
    p.adj = numeric(),
    p.adj.signif = character(),
    P_adjusted_global = numeric()
  )
}

# Dunn 检验先在每个标志物内部进行两两比较校正，再对所有标志物的
# 两两比较 P 值进行一次 BH 校正，用于显著性标记。
# 没必要让所有标志物再进行一次校正
# if (nrow(multi_group_dunn_summary) > 0) {
#   multi_group_dunn_summary$P_adjusted_global <- p.adjust(
#     multi_group_dunn_summary$p,
#     method = p_adjust_method
#   )
# }

for (marker in markers) {
  p <- multi_group_plots[[marker]]
  sig_dunn <- multi_group_dunn_summary %>%
    filter(Marker == marker, p.adj < alpha)
  if (nrow(sig_dunn) > 0) {
    raw_values <- df_raw[[marker]]
    plot_data_for_position <- data.frame(Value = raw_values)
    sig_dunn <- sig_dunn %>%
      mutate(
        y.position = make_y_positions(plot_data_for_position, n()),
        p.adj.signif = case_when(
          p.adj <= 0.001 ~ "***",
          p.adj <= 0.01 ~ "**",
          p.adj < alpha ~ "*",
          TRUE ~ "ns"
        )
      )
    p <- p + ggpubr::stat_pvalue_manual(
      sig_dunn,
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = TRUE,
      bracket.size = 0.35,
      size = 3,
      y.position = max(plot_data$Value) * 1.5
    )
  }
  multi_group_plots[[marker]] <- p
  print(p)
}

cat("\n五组总体比较结果：\n")
print(multi_group_summary)

save_checkpoint(
  list(
    summary = multi_group_summary,
    dunn = multi_group_dunn_summary,
    plots = multi_group_plots
  ),
  "Group_Cmp_Five_Group_Results.rds"
)
export_csv(multi_group_summary, "Five_Group_Comparison.csv")
export_csv(multi_group_dunn_summary, "Five_Group_Dunn_Posthoc.csv")

# 7.2 三种二分类比较：Wilcoxon 秩和检验
binary_group_results <- list()
binary_group_plots <- list()

for (comparison_name in names(binary_comparisons)) {
  selected_groups <- binary_comparisons[[comparison_name]]
  combined_label <- paste(selected_groups[-1], collapse = "+")

  df_binary_raw <- df_raw %>%
    filter(as.character(.data[[group_var]]) %in% selected_groups)
  df_binary_log <- df_log %>%
    filter(as.character(.data[[group_var]]) %in% selected_groups)

  binary_group <- ifelse(
    as.character(df_binary_raw[[group_var]]) == "0",
    "0",
    combined_label
  )
  binary_levels <- c("0", combined_label)

  marker_stats <- vector("list", length(markers))
  names(marker_stats) <- markers
  marker_plots <- vector("list", length(markers))
  names(marker_plots) <- markers

  for (marker in markers) {
    plot_data <- data.frame(
      Group = factor(binary_group, levels = binary_levels),
      Value = df_binary_raw[[marker]]
    ) %>%
      filter(!is.na(Group), is.finite(Value)) %>%
      droplevels()

    test_data <- data.frame(
      Group = factor(binary_group, levels = binary_levels),
      Log_Value = df_binary_log[[marker]]
    ) %>%
      filter(!is.na(Group), is.finite(Log_Value)) %>%
      droplevels()

    # 即使某标志物全为 0，仍显示原始浓度箱线图；该情形不进行秩检验。
      p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.62) +
      geom_jitter(width = 0.14, alpha = 1, size = 1, aes(color = Group)) +
      scale_x_discrete(drop = FALSE) +
      make_concentration_scale(plot_data$Value) +
      labs(
        title = marker,
        subtitle = paste0("Comparison: 0 vs ", combined_label),
        x = "Diagnostic group",
        y = "Conc."
      ) +
      theme_bw() +
      common_plot_theme +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    marker_plots[[marker]] <- p

    if (n_distinct(test_data$Group) < 2) {
      marker_stats[[marker]] <- data.frame(
        Comparison = comparison_name,
        Marker = marker,
        N = nrow(test_data),
        Test = "Wilcoxon rank-sum",
        Effect_Size_Type = "Wilcoxon effect size r",
        Effect_Size = NA_real_,
        P_value = NA_real_
      )
      next
    }

    wt <- tryCatch(
      wilcox.test(Log_Value ~ Group, data = test_data, exact = FALSE),
      error = function(e) NULL
    )
    eff <- tryCatch(
      rstatix::wilcox_effsize(test_data, Log_Value ~ Group),
      error = function(e) NULL
    )

    p_value <- if (is.null(wt)) NA_real_ else as.numeric(wt$p.value)
    effect_value <- if (!is.null(eff) && nrow(eff) > 0) {
      as.numeric(eff$effsize[1])
    } else {
      NA_real_
    }

    marker_stats[[marker]] <- data.frame(
      Comparison = comparison_name,
      Marker = marker,
      N = nrow(test_data),
      Test = "Wilcoxon rank-sum",
      Effect_Size_Type = "Wilcoxon effect size r",
      Effect_Size = effect_value,
      P_value = p_value,
      stringsAsFactors = FALSE
    )
  }

  binary_summary <- bind_rows(marker_stats) %>%
    mutate(P_adjusted = P_value)

  for (marker in markers) {
    p <- marker_plots[[marker]]
    if (is.null(p)) next
    stat_row <- binary_summary %>% filter(Marker == marker)
    if (nrow(stat_row) == 1 && is.finite(stat_row$P_adjusted) &&
        stat_row$P_adjusted < alpha) {
      plot_data <- data.frame(
        Value = df_binary_raw[[marker]],
        Group = factor(binary_group, levels = binary_levels)
      ) %>%
        filter(is.finite(Value), !is.na(Group))
      y_pos <- make_y_positions(plot_data, 1)
      annotation <- data.frame(
        group1 = binary_levels[1],
        group2 = binary_levels[2],
        y.position = y_pos,
        p.adj = stat_row$P_adjusted,
        p.adj.signif = ifelse(
          stat_row$P_adjusted <= 0.001, "***",
          ifelse(stat_row$P_adjusted <= 0.01, "**", "*")
        )
      )
      p <- p + ggpubr::stat_pvalue_manual(
        annotation,
        label = "p.adj.signif",
        tip.length = 0.01,
        hide.ns = TRUE,
        bracket.size = 0.35,
        size = 3
      )
    }
    marker_plots[[marker]] <- p
    print(p)
  }

  binary_group_results[[comparison_name]] <- list(
    summary = binary_summary,
    plots = marker_plots,
    selected_groups = selected_groups,
    combined_label = combined_label
  )

  cat("\n二分类比较结果：", comparison_name, "\n")
  print(binary_summary)
}

binary_group_summary <- bind_rows(
  lapply(binary_group_results, function(x) x$summary)
)

save_checkpoint(
  list(results = binary_group_results, summary = binary_group_summary),
  "Group_Cmp_Binary_Results.rds"
)
export_csv(binary_group_summary, "Binary_Group_Comparison.csv")

# 8. 效应量汇总 ------------------------------------------------

# 五组比较报告 Kruskal-Wallis epsilon-squared；二组比较报告 Wilcoxon r。
effect_size_results <- bind_rows(
  multi_group_summary %>%
    transmute(
      Comparison = "Five_groups_0_1_2_3_4",
      Marker,
      N,
      Test,
      Effect_Size_Type,
      Effect_Size,
      P_value,
      P_adjusted = P_adjusted_across_markers
    ),
  binary_group_summary %>%
    transmute(
      Comparison,
      Marker,
      N,
      Test,
      Effect_Size_Type,
      Effect_Size,
      P_value,
      P_adjusted
    )
)

save_checkpoint(effect_size_results, "Group_Cmp_Effect_Size_Results.rds")
export_csv(effect_size_results, "Marker_Effect_Size.csv")

# 9. 标志物 Spearman 相关性热图 --------------------------------

# 相关性使用 log 转换后的浓度；Spearman 相关对严格单调变换不变。
cor_data <- df_log %>%
  dplyr::select(dplyr::all_of(markers)) %>%
  as.data.frame()

valid_markers <- vapply(
  cor_data,
  function(x) {
    x <- x[is.finite(x)]
    length(x) >= 5 && length(unique(x)) >= 2
  },
  logical(1)
)

if (any(!valid_markers)) {
  cat("\n以下标志物因有效样本不足或无变异而不参与相关性分析：\n")
  print(names(cor_data)[!valid_markers])
}

cor_data <- cor_data[, valid_markers, drop = FALSE]
if (ncol(cor_data) < 2) {
  warning("有效标志物不足2个，跳过相关性分析。")
  correlation_results <- NULL
} else {
  res_cor <- Hmisc::rcorr(as.matrix(cor_data), type = "spearman")

  correlation_results <- list(
    correlation = res_cor$r,
    p_value = res_cor$P,
    n = res_cor$n,
    markers = colnames(cor_data)
  )

  correlation_colors <- colorRampPalette(
    c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")
  )(200)

  corrplot::corrplot(
    res_cor$r,
    method = "color",
    type = "upper",
    order = "hclust",
    col = correlation_colors,
    tl.col = "black",
    tl.srt = 45,
    tl.cex = 0.8,
    tl.pos = "td",
    addCoef.col = "black",
    number.cex = 0.6,
    number.font = 1,
    title = "Spearman correlation matrix\n(log Conc.)",
    mar = c(0, 0, 2, 0)
  )

  save_checkpoint(correlation_results, "Group_Cmp_Correlation_Results.rds")
  export_csv(res_cor$r, "Spearman_Correlation.csv", row_names = TRUE)
  export_csv(res_cor$P, "Spearman_Pvalues.csv", row_names = TRUE)
  export_csv(res_cor$n, "Spearman_Sample_Size.csv", row_names = TRUE)
}

# 10. 保存完整分析结果 -----------------------------------------

final_results <- list(
  imputation_id = selected_imputation,
  df_raw = df_raw,
  df_log = df_log,
  markers = markers,
  log_offsets = log_offsets,
  group_levels = original_group_levels,
  binary_comparisons = binary_comparisons,
  table1 = table1,
  table1_csv = table1_csv,
  heatmap_long = heatmap_long,
  multi_group_summary = multi_group_summary,
  multi_group_dunn = multi_group_dunn_summary,
  binary_group_results = binary_group_results,
  binary_group_summary = binary_group_summary,
  effect_size_results = effect_size_results,
  correlation_results = correlation_results,
  plots = list(
    heatmap = p_heatmap,
    five_group = multi_group_plots,
    binary_group = lapply(binary_group_results, `[[`, "plots")
  )
)

save_checkpoint(final_results, "Group_Cmp_Results.rds")

cat("\n========================================\n")
cat("Exploratory group comparisons completed.\n")
cat("Imputation dataset:", selected_imputation, "\n")
cat("Markers:", length(markers), "\n")
cat("Output directory:", output_dir, "\n")
cat("Figures were displayed in RStudio and were not automatically saved.\n")
cat("========================================\n")
