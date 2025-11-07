# ===================================================================
# 06_plot_combined.R (全局统一配色版)
# 综合分析图绘制
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 绘制综合分析图（箱线图 + 趋势图）- 全局统一配色版
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表（必须包含 CONFIG$colors）
#'
#' @return patchwork组合图
#'
#' @examples
#' p <- plot_combined_analysis(combined_data, CONFIG)
#'
plot_combined_analysis <- function(combined_data, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # ========================================
  # 1. 验证全局颜色方案
  # ========================================
  
  if (is.null(CONFIG$colors) || is.null(CONFIG$colors$celltype)) {
    stop("❌ 全局颜色方案未初始化，请先调用 create_global_color_scheme()")
  }
  
  if (is.null(CONFIG$colors$density_zone)) {
    stop("❌ 密度区域颜色未初始化")
  }
  
  # ========================================
  # 2. 获取全局颜色方案
  # ========================================
  
  # ✅ 使用全局颜色
  celltype_colors_global <- CONFIG$colors$celltype
  zone_colors_global <- CONFIG$colors$density_zone
  
  # 全局所有细胞类型（保持顺序）
  all_celltypes <- names(celltype_colors_global)
  
  n_zones <- length(unique(combined_data$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf("\n📊 绘制综合分析图...\n"))
  cat(sprintf("   样本数: %d\n", length(unique(combined_data$sample))))
  cat(sprintf("   细胞类型: %d (全局)\n", length(all_celltypes)))
  cat(sprintf("   密度区域: %d\n", n_zones))
  
  # ========================================
  # 3. 准备数据
  # ========================================
  
  # 确保zone按顺序排列
  combined_data <- combined_data %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      # ✅ 使用全局细胞类型顺序
      celltype_clean = factor(celltype_clean, levels = all_celltypes),
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    )
  
  # ========================================
  # 4. 图1：箱线图（分面）
  # ========================================
  
  p1 <- ggplot(combined_data, aes(x = density_zone, y = percentage, fill = density_zone)) +
    geom_boxplot(
      alpha = 0.8, 
      outlier.shape = 16, 
      outlier.size = 1.5, 
      color = "gray30", 
      linewidth = 0.5
    ) +
    scale_fill_manual(
      values = zone_colors_global,  # ✅ 使用全局区域颜色
      guide = "none"
    ) +
    # ✅ 使用全局细胞类型分面
    facet_wrap(
      ~celltype_clean, 
      scales = "free_y", 
      ncol = 4,
      drop = FALSE  # ✅ 显示所有细胞类型
    ) +
    labs(
      title = "Cell Type Percentage Distribution by Density Zone",
      subtitle = sprintf("Data from %d samples | %d cell types (global colors)", 
                        length(unique(combined_data$sample)),
                        length(all_celltypes)),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30", margin = margin(b = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 10, 15)
    )
  
  # ========================================
  # 5. 图2：趋势图
  # ========================================
  
  # 计算趋势数据
  trend_data <- combined_data %>%
    dplyr::group_by(celltype_clean, zone_numeric, density_zone) %>%
    dplyr::summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      se_pct = sd(percentage, na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = "drop"
    )
  
  # ✅ 补全缺失的组合（某些细胞类型在某些区域可能没有数据）
  complete_trend_grid <- expand.grid(
    celltype_clean = all_celltypes,
    zone_numeric = 0:(n_zones - 1),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      density_zone = factor(sprintf("Zone_%d", zone_numeric), levels = zone_levels)
    )
  
  trend_data <- complete_trend_grid %>%
    dplyr::left_join(trend_data, by = c("celltype_clean", "zone_numeric", "density_zone")) %>%
    dplyr::mutate(
      mean_pct = ifelse(is.na(mean_pct), 0, mean_pct),
      se_pct = ifelse(is.na(se_pct), 0, se_pct),
      n_obs = ifelse(is.na(n_obs), 0, n_obs),
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  p2 <- ggplot(trend_data, aes(x = zone_numeric, y = mean_pct, color = celltype_clean, group = celltype_clean)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = pmax(mean_pct - se_pct, 0), ymax = mean_pct + se_pct),  # 不低于0
      width = 0.2, 
      linewidth = 0.8,
      alpha = 0.7
    ) +
    scale_color_manual(
      values = celltype_colors_global,  # ✅ 使用全局细胞类型颜色
      name = "Cell Type",
      breaks = all_celltypes,  # ✅ 显示所有类型
      drop = FALSE  # ✅ 不丢弃未使用的级别
    ) +
    scale_x_continuous(
      breaks = 0:(n_zones - 1),
      labels = zone_levels
    ) +
    labs(
      title = "Cell Type Enrichment Trend Across Density Zones",
      subtitle = "Mean ± SE across all samples",
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Mean Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30", margin = margin(b = 10)),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.8, "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 15, 15, 15)
    )
  
  # ========================================
  # 6. 添加zone颜色参考条（底部）
  # ========================================
  
  zone_ref_data <- data.frame(
    zone_numeric = 0:(n_zones - 1),
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 0
  )
  
  max_y <- max(trend_data$mean_pct + trend_data$se_pct, na.rm = TRUE)
  
  p2 <- p2 +
    geom_tile(
      data = zone_ref_data,
      aes(x = zone_numeric, y = y_position, fill = density_zone),
      height = max_y * 0.05,  # 动态高度
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = zone_colors_global,  # ✅ 使用全局区域颜色
      guide = "none"
    )
  
  # ========================================
  # 7. 合并两个图
  # ========================================
  
  p_combined <- p1 / p2 + 
    plot_layout(heights = c(2, 1.2)) +
    plot_annotation(
      caption = "Global color scheme applied: colors are consistent across all samples and plots",
      theme = theme(
        plot.caption = element_text(size = 9, color = "gray40", hjust = 1, margin = margin(t = 10))
      )
    )
  
  cat("   ✅ 综合分析图绘制完成\n")
  
  return(p_combined)
}


# ========================================
# 向后兼容函数（已弃用）
# ========================================

#' @deprecated 请使用 get_zone_colors() from 01_color_schemes.R
get_zone_colors <- function(n_zones) {
  warning("get_zone_colors() 已弃用，请使用 01_color_schemes.R 中的版本")
  
  colorRampPalette(c(
    "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
  ))(n_zones)
}

#' @deprecated 请使用 get_celltype_colors() from 01_color_schemes.R
get_celltype_colors <- function(celltypes) {
  warning("get_celltype_colors() 已弃用，请使用 01_color_schemes.R 中的版本")
  
  require(RColorBrewer)
  n <- length(celltypes)
  
  if (n <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n), "Set2")
  } else if (n <= 12) {
    colors <- RColorBrewer::brewer.pal(n, "Set3")
  } else {
    colors <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3")
    )[1:n]
  }
  
  names(colors) <- celltypes
  return(colors)
}

cat("✅ 06_plot_combined.R 已加载（全局统一配色版）\n")