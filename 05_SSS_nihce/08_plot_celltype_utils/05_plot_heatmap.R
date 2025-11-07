# ===================================================================
# 05_plot_heatmap.R (全局统一配色版)
# 合并热图绘制
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 绘制合并热图（所有样本）- 全局统一配色版
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表（必须包含 CONFIG$colors）
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
  
  # 确保zone按顺序排列
  n_zones <- length(unique(combined_data$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf("\n📊 绘制跨样本热图...\n"))
  cat(sprintf("   样本数: %d\n", length(unique(combined_data$sample))))
  cat(sprintf("   细胞类型: %d (全局)\n", length(all_celltypes)))
  cat(sprintf("   密度区域: %d\n", n_zones))
  
  # ========================================
  # 3. 计算平均百分比
  # ========================================
  
  heatmap_data <- combined_data %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      sd_pct = sd(percentage, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    )
  
  # ✅ 补全缺失的细胞类型（某些类型可能在所有样本中都不存在某些区域）
  complete_grid <- expand.grid(
    density_zone = zone_levels,
    celltype_clean = all_celltypes,
    stringsAsFactors = FALSE
  )
  
  heatmap_data <- complete_grid %>%
    dplyr::left_join(heatmap_data, by = c("density_zone", "celltype_clean")) %>%
    dplyr::mutate(
      mean_pct = ifelse(is.na(mean_pct), 0, mean_pct),
      sd_pct = ifelse(is.na(sd_pct), 0, sd_pct),
      n_samples = ifelse(is.na(n_samples), 0, n_samples),
      # ✅ 使用全局顺序
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  # ========================================
  # 4. 热图主体
  # ========================================
  
  p <- ggplot(heatmap_data, aes(x = density_zone, y = celltype_clean, fill = mean_pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = ifelse(mean_pct > 0.5, sprintf("%.1f", mean_pct), "")),  # 只显示 >0.5% 的值
      size = 3.5, 
      color = "black", 
      fontface = "bold"
    ) +
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
    # ✅ Y轴使用全局细胞类型顺序
    scale_y_discrete(
      breaks = all_celltypes,
      drop = FALSE  # ✅ 显示所有细胞类型
    ) +
    labs(
      title = "Cell Type Composition Across Density Zones (All Samples)",
      subtitle = sprintf("Averaged across %d samples | %d cell types (global)", 
                        length(unique(combined_data$sample)), 
                        length(all_celltypes)),
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
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # ========================================
  # 5. Zone颜色参考条（顶部）
  # ========================================
  
  zone_bar_data <- data.frame(
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 1
  )
  
  p_zone_bar <- ggplot(zone_bar_data, aes(x = density_zone, y = y_position, fill = density_zone)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(
      values = zone_colors_global,  # ✅ 使用全局区域颜色
      guide = "none"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # ========================================
  # 6. 细胞类型颜色参考条（左侧）
  # ========================================
  
  celltype_bar_data <- data.frame(
    celltype_clean = factor(all_celltypes, levels = all_celltypes),
    x_position = 1
  )
  
  p_celltype_bar <- ggplot(celltype_bar_data, aes(x = x_position, y = celltype_clean, fill = celltype_clean)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_manual(
      values = celltype_colors_global,  # ✅ 使用全局细胞类型颜色
      guide = "none"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(
      breaks = all_celltypes,
      drop = FALSE
    ) +
    theme_void() +
    theme(
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # ========================================
  # 7. 合并图形（添加细胞类型颜色条）
  # ========================================
  
  # 空白占位符（左上角）
  p_blank <- ggplot() + theme_void()
  
  # 使用 patchwork 布局
  p_final <- (p_blank | p_zone_bar) / 
             (p_celltype_bar | p) + 
    plot_layout(
      widths = c(0.03, 1),   # 左侧颜色条宽度
      heights = c(0.05, 1)   # 顶部颜色条高度
    ) +
    plot_annotation(
      caption = "Colors are consistent across all samples (global color scheme)",
      theme = theme(
        plot.caption = element_text(size = 9, color = "gray40", hjust = 1, margin = margin(t = 10))
      )
    )
  
  cat("   ✅ 热图绘制完成\n")
  
  return(p_final)
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

cat("✅ 05_plot_heatmap.R 已加载（全局统一配色版）\n")