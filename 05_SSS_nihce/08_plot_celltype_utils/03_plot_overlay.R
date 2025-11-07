# ===================================================================
# 03_plot_overlay.R (修复版)
# 细胞类型+密度叠加图（使用raster，无网格线）
# Author: Assistant (Fixed Version)
# Date: 2025-11-07
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
plot_celltype_density_overlay <- function(df, density_data, sample_id, CONFIG) {
  
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  require(RANN)
  
  # ========================================
  # 1. 准备数据
  # ========================================
  
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  zone_colors <- CONFIG$colors$zone_colors %||% get_zone_colors(n_zones)
  
  # 清理 celltype
  df$celltype_clean <- as.character(df$celltype_clean)
  df$celltype_clean[is.na(df$celltype_clean) | df$celltype_clean == ""] <- "Unknown"
  all_celltypes <- sort(unique(df$celltype_clean))
  
  # 获取配置的颜色
  celltype_colors <- CONFIG$colors$celltype_colors
  
  # 确保所有类型都有颜色
  missing_types <- setdiff(all_celltypes, names(celltype_colors))
  if (length(missing_types) > 0) {
    extra_colors <- rainbow(length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors <- c(celltype_colors, extra_colors)
  }
  
  # 只保留实际存在的类型
  celltype_colors <- celltype_colors[all_celltypes]
  
  cat("   📊 Celltype 颜色映射:\n")
  for (ct in all_celltypes) {
    cat(sprintf("      %s → %s\n", ct, celltype_colors[ct]))
  }
  
  # ========================================
  # 2. 坐标范围
  # ========================================
  
  col_range_raw <- density_data$col_range
  row_range_raw <- density_data$row_range
  
  if (!is.null(density_data$col_range_expanded)) {
    col_limits <- density_data$col_range_expanded
    row_limits <- density_data$row_range_expanded
  } else {
    expand_margin <- 0.1
    col_margin <- diff(col_range_raw) * expand_margin
    row_margin <- diff(row_range_raw) * expand_margin
    col_limits <- c(col_range_raw[1] - col_margin, col_range_raw[2] + col_margin)
    row_limits <- c(row_range_raw[1] - row_margin, row_range_raw[2] + row_margin)
  }
  
  cat(sprintf("   ✅ 绘图范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_limits[1], col_limits[2], row_limits[1], row_limits[2]))
  
  # ========================================
  # 3. 准备等高线数据
  # ========================================
  
  zone_density_ranges <- density_data$grid %>%
    group_by(density_zone) %>%
    summarise(
      min_density = min(density_norm, na.rm = TRUE),
      max_density = max(density_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(density_zone = factor(density_zone, levels = zone_levels)) %>%
    arrange(density_zone)
  
  zone_labels <- zone_density_ranges %>%
    mutate(zone_label = sprintf("%s (%.2f-%.2f)", density_zone, min_density, max_density)) %>%
    pull(zone_label)
  names(zone_labels) <- as.character(zone_density_ranges$density_zone)
  
  contour_breaks <- density_data$equal_breaks
  contour_data <- density_data$grid %>%
    mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # ========================================
  # 4. 准备细胞数据
  # ========================================
  
  df_filtered <- df %>% 
    filter(!is.na(density_zone)) %>%
    mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  # 计算细胞大小
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
  
  # ========================================
  # 5. 绘图（关键修复）
  # ========================================
  
  p <- ggplot() +
    # Layer 1: 细胞类型（底层）
    # ✅ 关键：使用 geom_point 代替 geom_tile 测试
    geom_point(
      data = df_filtered,
      aes(x = col, y = row, color = celltype_clean),
      size = square_size * 2,  # point 的 size 需要调整
      shape = 15,  # 正方形
      alpha = 1
    ) +
    scale_color_manual(
      values = celltype_colors,
      name = "Cell Type",
      breaks = all_celltypes,
      guide = guide_legend(
        order = 2,
        override.aes = list(size = 5, alpha = 1),
        title.position = "top",
        ncol = 1
      )
    ) +
    
    # 新的scale用于density zones
    ggnewscale::new_scale_fill() +
    
    # Layer 2: Zone填充
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.25,
      interpolate = TRUE
    ) +
    scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones\n(Zone_0=Core Red → Zone_9=Outer Blue)",
      breaks = zone_levels,
      na.value = "transparent",
      guide = guide_legend(
        order = 1,
        override.aes = list(alpha = 0.7),
        title.position = "top",
        ncol = 1
      )
    ) +
    
    # 为等高线准备新的color scale
    ggnewscale::new_scale_color() +
    
    # Layer 3: 等高线
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
    
    # 坐标系统
    scale_x_continuous(limits = col_limits, expand = c(0, 0)) +
    scale_y_reverse(limits = rev(row_limits), expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = col_limits, ylim = rev(row_limits), clip = "off") +
    
    # 标题
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = sprintf("Bottom = Cell types | Middle = Density zones | Top = %d contour lines", 
                        length(contour_breaks)),
      x = NULL, y = NULL
    ) +
    
    # 主题
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9.5),
      legend.key.size = unit(0.8, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      plot.margin = margin(15, 25, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

