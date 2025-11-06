# ===================================================================
# 06_plot_combined.R
# 综合分析图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制综合分析图（箱线图 + 趋势图）
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表
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
  
  # 获取统一的颜色方案
  n_zones <- length(unique(combined_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  celltype_colors <- get_celltype_colors(unique(combined_data$celltype_clean))
  
  # 确保zone按顺序排列
  combined_data <- combined_data %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    )
  
  # 1. 箱线图
  p1 <- ggplot(combined_data, aes(x = density_zone, y = percentage, fill = density_zone)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    facet_wrap(~celltype_clean, scales = "free_y", ncol = 4) +
    labs(
      title = "Cell Type Percentage Distribution by Density Zone",
      subtitle = sprintf("Data from %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
  
  # 2. 趋势图
  trend_data <- combined_data %>%
    dplyr::group_by(celltype_clean, zone_numeric, density_zone) %>%
    dplyr::summarise(
      mean_pct = mean(percentage),
      se_pct = sd(percentage) / sqrt(n()),
      .groups = "drop"
    )
  
  p2 <- ggplot(trend_data, aes(x = zone_numeric, y = mean_pct, color = celltype_clean, group = celltype_clean)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct), 
      width = 0.2, 
      linewidth = 0.8,
      alpha = 0.7
    ) +
    scale_color_manual(values = celltype_colors, name = "Cell Type") +
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
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  # 3. 添加zone颜色参考条
  zone_ref_data <- data.frame(
    zone_numeric = 0:(n_zones - 1),
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 0
  )
  
  p2 <- p2 +
    geom_tile(
      data = zone_ref_data,
      aes(x = zone_numeric, y = y_position, fill = density_zone),
      height = max(trend_data$mean_pct) * 0.05,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = zone_colors, guide = "none")
  
  # 合并
  p_combined <- p1 / p2 + plot_layout(heights = c(2, 1.2))
  
  return(p_combined)
}