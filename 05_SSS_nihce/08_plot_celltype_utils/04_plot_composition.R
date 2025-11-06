# ===================================================================
# 04_plot_composition.R
# 区域组成柱状图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制zone组成柱状图
#'
#' @param zone_composition zone组成数据框
#' @param sample_id 样本ID
#' @param CONFIG 配置列表
#'
#' @return patchwork组合图（细胞类型组成 + spot数量）
#'
#' @examples
#' p <- plot_zone_composition(zone_comp, "Sample_01", CONFIG)
#'
plot_zone_composition <- function(zone_composition, sample_id, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # 使用统一的颜色方案
  n_zones <- length(unique(zone_composition$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  celltype_colors <- get_celltype_colors(unique(zone_composition$celltype_clean))
  
  # 确保zone按 Zone_0, Zone_1, ... 排序
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  zone_composition <- zone_composition %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 图1：细胞类型组成堆叠柱状图
  p1 <- ggplot(zone_composition, aes(x = density_zone, y = percentage, fill = celltype_clean)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = sprintf("Cell Type Composition by Density Zone - %s", sample_id),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  # 图2：Zone的spot数量
  zone_totals <- zone_composition %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(total = sum(count), .groups = "drop")
  
  p2 <- ggplot(zone_totals, aes(x = density_zone, y = total, fill = density_zone)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.5) +
    geom_text(aes(label = total), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Total Spots per Density Zone",
      x = "Density Zone (Zone_0=Core → Higher=Outer)",
      y = "Count"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # 合并
  p_combined <- p1 / p2 + plot_layout(heights = c(2, 1))
  
  return(p_combined)
}