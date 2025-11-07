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
  
  # 加载必要的包
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  require(RANN)
  
  # ========================================
  # 1. 数据准备和验证
  # ========================================
  
  # 获取zone信息
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  # 获取颜色方案
  zone_colors <- CONFIG$colors$zone_colors %||% get_zone_colors(n_zones)
  celltype_colors <- CONFIG$colors$celltype_colors %||% get_celltype_colors(unique(df$celltype_clean))
  
  # ✅ 关键修复1：确保celltype_clean无NA且都有对应颜色
  df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
  
  # 获取所有实际存在的细胞类型
  actual_celltypes <- unique(as.character(df$celltype_clean))
  
  # 检查是否所有celltype都有颜色，如果没有则补充
  missing_types <- setdiff(actual_celltypes, names(celltype_colors))
  if (length(missing_types) > 0) {
    cat(sprintf("   ⚠️  发现未配色的细胞类型: %s\n", paste(missing_types, collapse=", ")))
    # 为缺失的类型生成颜色
    extra_colors <- rainbow(length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors <- c(celltype_colors, extra_colors)
  }
  
  # ✅ 关键修复2：确保factor水平与颜色表完全一致
  df_filtered <- df %>% 
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = names(celltype_colors))
    )
  
  # 调试输出
  cat(sprintf("   📊 细胞类型统计:\n"))
  celltype_table <- table(df_filtered$celltype_clean)
  for(i in 1:length(celltype_table)) {
    cat(sprintf("      %s: %d (颜色: %s)\n", 
                names(celltype_table)[i], 
                celltype_table[i],
                celltype_colors[names(celltype_table)[i]]))
  }
  
  # ========================================
  # 2. 坐标范围设置
  # ========================================
  
  # 使用切片的实际范围
  col_range_raw <- density_data$col_range
  row_range_raw <- density_data$row_range
  
  # ✅ 关键修复3：增大边距以确保zone明显溢出
  expand_margin <- CONFIG$plot$expand_margin %||% 0.1  # 使用配置的边距，默认10%
  col_margin <- diff(col_range_raw) * expand_margin
  row_margin <- diff(row_range_raw) * expand_margin
  
  # 应用边距
  col_limits <- c(col_range_raw[1] - col_margin, col_range_raw[2] + col_margin)
  row_limits <- c(row_range_raw[1] - row_margin, row_range_raw[2] + row_margin)
  
  cat(sprintf("   ✅ 原始切片范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range_raw[1], col_range_raw[2], row_range_raw[1], row_range_raw[2]))
  cat(sprintf("   ✅ 扩展边距: %.0f%%\n", expand_margin * 100))
  cat(sprintf("   ✅ 添加边距后范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_limits[1], col_limits[2], row_limits[1], row_limits[2]))
  
  # ========================================
  # 3. 准备等高线数据
  # ========================================
  
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
  
  # 准备contour数据
  contour_data <- density_data$grid %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # ========================================
  # 4. 计算细胞大小
  # ========================================
  
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
  
  # 获取实际存在的细胞类型（用于图例）
  celltypes_present <- levels(droplevels(df_filtered$celltype_clean))
  
  # ========================================
  # 5. 绘图
  # ========================================
  
  p <- ggplot() +
    # Layer 1: 细胞类型（底层）
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, fill = celltype_clean),
      width = square_size,
      height = square_size,
      alpha = 1  # ✅ 修复：提高不透明度使颜色更清晰
    ) +
    scale_fill_manual(
      values = celltype_colors,
      name = "Cell Type",
      breaks = celltypes_present,  # 只显示实际存在的类型
      drop = FALSE,  # 保持所有水平
      guide = guide_legend(
        order = 2,
        override.aes = list(alpha = 1, size = 5),  # ✅ 修复：确保图例清晰可见
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        keywidth = unit(0.8, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    
    # 新的scale用于density zones
    ggnewscale::new_scale_fill() +
    
    # Layer 2: Zone填充（半透明覆盖）
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.25,  # ✅ 修复：降低透明度，让底层celltype更清晰
      interpolate = TRUE,  # 平滑插值
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones\n(Zone_0=Core Red → Zone_9=Outer Blue)",
      breaks = zone_levels,
      na.value = "transparent",
      drop = FALSE,
      guide = guide_legend(
        order = 1,
        override.aes = list(alpha = 0.7),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        keywidth = unit(0.8, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    
    # 为等高线准备新的color scale
    ggnewscale::new_scale_color() +
    
    # Layer 3: 等高线边界（保留颜色渐变）
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
      guide = "none"  # 不显示等高线的图例
    ) +
    
    # ✅ 关键修复4：明确设置坐标系统和限制
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
      clip = "off"  # 允许图形元素超出边界
    ) +
    
    # 标题和标签
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = sprintf("Bottom = Cell types | Middle = Density zones (raster) | Top = %d contour lines", 
                        length(contour_breaks)),
      x = NULL,
      y = NULL
    ) +
    
    # 主题设置
    theme_void() +
    theme(
      # 标题
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      
      # 图例
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
      
      # 图边距
      plot.margin = margin(15, 25, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  
  return(p)
}

