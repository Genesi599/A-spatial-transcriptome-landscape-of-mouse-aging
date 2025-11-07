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
  
  # ========================================
  # ✅ 强制修复 celltype 颜色映射
  # ========================================
  df$celltype_clean <- as.character(df$celltype_clean)
  df$celltype_clean[is.na(df$celltype_clean) | df$celltype_clean == ""] <- "Unknown"
  
  # 确保所有 celltype 都有颜色
  all_celltypes <- unique(df$celltype_clean)
  missing_types <- setdiff(all_celltypes, names(celltype_colors))
  if (length(missing_types) > 0) {
    cat(sprintf("   ⚠️  补充颜色: %s\n", paste(missing_types, collapse=", ")))
    extra_colors <- rainbow(length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors <- c(celltype_colors, extra_colors)
  }
  
  # ✅ 关键修复：确保颜色向量包含所有实际存在的类型，并且顺序一致
  celltype_colors <- celltype_colors[names(celltype_colors) %in% all_celltypes]
  
  # 打印验证（带颜色值）
  cat("   📊 Celltype 颜色映射:\n")
  for (ct in sort(all_celltypes)) {
    cat(sprintf("      %s → %s\n", ct, celltype_colors[ct]))
  }
  
  # ✅ 创建统一的 celltype levels（按字母顺序）
  celltype_levels <- sort(names(celltype_colors))
  
  # ========================================
  # 2. 坐标范围设置
  # ========================================
  
  # 使用切片的实际范围
  col_range_raw <- density_data$col_range
  row_range_raw <- density_data$row_range
  
  # 使用已经扩展好的范围
  if (!is.null(density_data$col_range_expanded)) {
    col_limits <- density_data$col_range_expanded
    row_limits <- density_data$row_range_expanded
    expand_margin <- 0  # 已经扩展过了
  } else {
    # 如果没有扩展范围，则手动扩展
    expand_margin <- CONFIG$plot$expand_margin %||% 0.1
    col_margin <- diff(col_range_raw) * expand_margin
    row_margin <- diff(row_range_raw) * expand_margin
    col_limits <- c(col_range_raw[1] - col_margin, col_range_raw[2] + col_margin)
    row_limits <- c(row_range_raw[1] - row_margin, row_range_raw[2] + row_margin)
  }
  
  cat(sprintf("   ✅ 原始切片范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range_raw[1], col_range_raw[2], row_range_raw[1], row_range_raw[2]))
  cat(sprintf("   ✅ 绘图范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
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
  
  cat(sprintf("   ✅ 等高线边界（共%d条）\n", length(contour_breaks)))
  
  # 准备contour数据
  contour_data <- density_data$grid %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # ========================================
  # 4. 准备细胞数据（关键！）
  # ========================================
  
  # ✅ 强制转换为 factor，使用统一的 levels
  df_filtered <- df %>% 
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = celltype_levels)  # 使用统一的 levels
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
  
  # 获取实际存在的细胞类型（用于图例显示）
  celltypes_present <- levels(droplevels(df_filtered$celltype_clean))
  
  # ✅ 调试输出：检查数据
  cat(sprintf("   🔍 调试信息:\n"))
  cat(sprintf("      df_filtered 行数: %d\n", nrow(df_filtered)))
  cat(sprintf("      celltype_clean 类型: %s\n", class(df_filtered$celltype_clean)))
  cat(sprintf("      celltype_clean levels: %s\n", paste(levels(df_filtered$celltype_clean), collapse=", ")))
  cat(sprintf("      实际存在的类型: %s\n", paste(celltypes_present, collapse=", ")))
  cat(sprintf("      颜色向量长度: %d\n", length(celltype_colors)))
  cat(sprintf("      颜色向量名字: %s\n", paste(names(celltype_colors), collapse=", ")))
  
  # ========================================
  # 5. 绘图
  # ========================================
  
  p <- ggplot() +
    # Layer 1: 细胞类型（底层）
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, fill = celltype_clean),  # 只用 fill，不用 color
      width = square_size,
      height = square_size,
      alpha = 1,
      color = NA  # ✅ 明确设置边框为 NA
    ) +
    # ✅ 关键修复：强制指定 values，使用命名向量
    scale_fill_manual(
      values = celltype_colors,  # 这必须是命名向量
      breaks = celltypes_present,  # 只显示实际存在的
      labels = celltypes_present,  # 使用相同的标签
      name = "Cell Type",
      drop = TRUE,  # 删除未使用的 levels
      na.value = "gray50",  # NA 值用灰色
      guide = guide_legend(
        order = 2,
        override.aes = list(
          alpha = 1,
          color = NA,  # ✅ 图例中也不要边框
          size = 1
        ),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        keywidth = unit(1.2, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    
    # 新的scale用于density zones
    ggnewscale::new_scale_fill() +
    
    # Layer 2: Zone填充（半透明覆盖）
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.25,
      interpolate = TRUE,
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
      guide = "none"
    ) +
    
    # 坐标系统
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
      clip = "off"
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
      plot.margin = margin(15, 25, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 添加原始组织边界框
  if (TRUE) {
    p <- p + 
      annotate(
        "rect",
        xmin = col_range_raw[1], xmax = col_range_raw[2],
        ymin = row_range_raw[1], ymax = row_range_raw[2],
        fill = NA,
        color = "black",
        linewidth = 0.8,
        linetype = "dashed",
        alpha = 0.3
      )
  }
  
  return(p)
}

