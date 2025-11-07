# ===================================================================
# 03_plot_overlay.R
# 细胞类型+密度叠加图（使用raster，无网格线）
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制细胞类型和密度区域叠加图
#'
#' @param df 数据框，包含细胞类型和坐标信息
#' @param density_data 密度计算结果（来自 calculate_density_zones）
#' @param sample_id 样本ID
#' @param CONFIG 配置列表
#'
#' @return ggplot对象
#'
#' @details
#' 图层从下到上：
#' 1. 细胞类型（geom_tile，color+fill）
#' 2. 密度区域（geom_raster，fill，透明度0.3）
#' 3. 等高线边界（geom_contour，细线）
#'
#' @examples
#' p <- plot_celltype_density_overlay(df, density_data, "Sample_01", CONFIG)
#'
plot_celltype_density_overlay <- function(df, density_data, sample_id, CONFIG) {
  
  # 加载必要的包
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  require(RANN)
  
  # 获取zone信息
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  # 使用统一的颜色方案
  zone_colors <- get_zone_colors(n_zones)
  celltype_colors <- get_celltype_colors(unique(df$celltype_clean))
  
  # 使用切片的实际范围，并添加边距
  col_range_raw <- density_data$col_range
  row_range_raw <- density_data$row_range
  
  # 计算边距（范围的5%）
  col_margin <- diff(col_range_raw) * 0.05
  row_margin <- diff(row_range_raw) * 0.05
  
  # 应用边距
  col_limits <- c(col_range_raw[1] - col_margin, col_range_raw[2] + col_margin)
  row_limits <- c(row_range_raw[1] - row_margin, row_range_raw[2] + row_margin)
  
  cat(sprintf("   ✅ 原始切片范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range_raw[1], col_range_raw[2], row_range_raw[1], row_range_raw[2]))
  cat(sprintf("   ✅ 添加边距后范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_limits[1], col_limits[2], row_limits[1], row_limits[2]))
  
  # 计算每个zone的密度范围
  zone_density_ranges <- density_data$grid %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(
      min_density = min(density_norm, na.rm = TRUE),
      max_density = max(density_norm, na.rm = TRUE),
      mean_density = mean(density_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels)) %>%
    dplyr::arrange(density_zone)
  
  # 为zone创建带密度信息的标签
  zone_labels <- zone_density_ranges %>%
    dplyr::mutate(
      zone_label = sprintf("%s (%.2f-%.2f)", 
                          density_zone, 
                          min_density, 
                          max_density)
    ) %>%
    dplyr::pull(zone_label)
  
  names(zone_labels) <- as.character(zone_density_ranges$density_zone)
  
  # 使用等距边界
  contour_breaks <- density_data$equal_breaks
  
  cat(sprintf("   ✅ 等高线边界（共%d条）:\n", length(contour_breaks)))
  for (i in seq_along(contour_breaks)) {
    cat(sprintf("      %.2f\n", contour_breaks[i]))
  }
  
  # 准备数据
  contour_data <- density_data$grid %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  df_filtered <- df %>% 
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 自动计算细胞正方形大小
  if (nrow(df_filtered) > 10000) {
    sample_idx <- sample(nrow(df_filtered), 10000)
    coords_sample <- df_filtered[sample_idx, c("col", "row")]
  } else {
    coords_sample <- df_filtered[, c("col", "row")]
  }
  
  nn_dist <- RANN::nn2(coords_sample, k = 2)$nn.dists[, 2]
  median_dist <- median(nn_dist, na.rm = TRUE)
  square_size <- median_dist * 1.0
  
  cat(sprintf("   📏 细胞正方形大小: %.3f\n", square_size))
  
  # 等高线颜色
  contour_colors <- get_contour_colors(length(contour_breaks))
  
  # 获取实际存在的细胞类型
  celltypes_present <- unique(df_filtered$celltype_clean)
  
  # 绘图
  p <- ggplot() +
    # 1. 细胞类型（color + fill）
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, color = celltype_clean, fill = celltype_clean),
      width = square_size,
      height = square_size,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = celltype_colors,
      name = "Cell Type",
      breaks = celltypes_present,
      guide = guide_legend(
        order = 2,
        override.aes = list(alpha = 1),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        keywidth = unit(1.2, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    scale_color_manual(
      values = celltype_colors,
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    
    # 2. Zone填充
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.3,
      interpolate = FALSE,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones\n(Zone_0=Core Red → Higher=Outer Blue)",
      breaks = zone_levels,
      na.value = "transparent",
      drop = FALSE,
      guide = guide_legend(
        order = 1,
        override.aes = list(alpha = 0.7),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        keywidth = unit(1.0, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    
    # 在等高线之前重置 color 通道
    ggnewscale::new_scale_color() +

    # 3. 等高线边界
    geom_contour(
      data = contour_data,
      aes(x = col, y = row, z = density_norm, color = after_stat(level)),
      breaks = contour_breaks,
      linewidth = 0.5,
      alpha = 0.7
    ) +
    scale_color_gradientn(
      colors = contour_colors,
      limits = c(min(contour_breaks), max(contour_breaks)),
      guide = "none"
    ) +
    
    # 坐标和主题
    scale_x_continuous(
      limits = col_limits,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = rev(row_limits),
      expand = c(0, 0)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = col_limits,
      ylim = rev(row_limits),
      clip = "on"
    ) +
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = sprintf("Bottom = Cell types | Middle = Density zones (raster) | Top = %d contour lines", 
                        length(contour_breaks))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = unit(0.8, "cm"),
      legend.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 6)),
      legend.text = element_text(size = 9.5, lineheight = 1.3, hjust = 0, 
                                margin = margin(l = 2, r = 5, t = 2, b = 2)),
      legend.key = element_rect(color = "gray70", fill = NA, linewidth = 0.3),
      legend.key.width = unit(1.0, "cm"),
      legend.key.height = unit(0.7, "cm"),
      legend.key.spacing.y = unit(0.25, "cm"),
      legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.5),
      legend.margin = margin(12, 12, 12, 12),
      plot.margin = margin(15, 25, 15, 15)
    )
  
  return(p)
}