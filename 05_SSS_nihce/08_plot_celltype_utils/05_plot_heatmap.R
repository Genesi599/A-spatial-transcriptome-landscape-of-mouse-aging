# ===================================================================
# 05_plot_heatmap.R
# 合并热图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制合并热图（所有样本）
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表
#'
#' @return patchwork组合图（zone颜色条 + 热图）
#'
#' @examples
#' p <- plot_combined_heatmap(combined_data, CONFIG)
#'
plot_combined_heatmap <- function(combined_data, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # 计算平均百分比
  heatmap_data <- combined_data %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(
      mean_pct = mean(percentage),
      sd_pct = sd(percentage),
      n_samples = n(),
      .groups = "drop"
    )
  
  # 确保zone按顺序排列
  n_zones <- length(unique(heatmap_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  heatmap_data <- heatmap_data %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 热图主体
  p <- ggplot(heatmap_data, aes(x = density_zone, y = celltype_clean, fill = mean_pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_pct)), size = 3.5, color = "black", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white", "#fee090", "#fc8d59", "#d73027"),
      name = "Mean %",
      limits = c(0, NA),
      guide = guide_colorbar(
        barwidth = 1.5,
        barheight = 15,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = "Cell Type Composition Across Density Zones (All Samples)",
      subtitle = sprintf("Averaged across %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30", margin = margin(b = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1)
    )
  
  # 添加zone颜色参考条（顶部）
  zone_bar_data <- data.frame(
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 1
  )
  
  p_zone_bar <- ggplot(zone_bar_data, aes(x = density_zone, y = y_position, fill = density_zone)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # 合并图形
  p_final <- p_zone_bar / p + plot_layout(heights = c(0.05, 1))
  
  return(p_final)
}