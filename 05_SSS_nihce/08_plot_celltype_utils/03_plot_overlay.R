# ===================================================================
# 03_plot_overlay.R (全局统一配色版)
# 细胞类型+密度叠加图（使用raster，无网格线）
# Author: Assistant (Global Color Scheme)
# Date: 2025-11-07
# ===================================================================

#' 绘制细胞类型和密度区域叠加图（全局统一配色版）
#'
#' @param df 数据框，包含细胞类型和坐标信息
#' @param density_data 密度计算结果（来自 calculate_density_zones）
#' @param sample_id 样本ID
#' @param CONFIG 配置列表（必须包含 CONFIG$colors）
#'
#' @return ggplot对象
#'
plot_celltype_density_overlay <- function(df, density_data, sample_id, CONFIG) {
  
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  require(RANN)
  
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
  # 2. 准备数据
  # ========================================
  
  # 获取全局颜色方案
  celltype_colors_global <- CONFIG$colors$celltype
  zone_colors_global <- CONFIG$colors$density_zone
  
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  # 清理 celltype
  df$celltype_clean <- as.character(df$celltype_clean)
  df$celltype_clean[is.na(df$celltype_clean) | df$celltype_clean == ""] <- "Unknown"
  
  # ✅ 获取当前样本的细胞类型（可能是全局的子集）
  sample_celltypes <- sort(unique(df$celltype_clean))
  
  # ✅ 检查是否有未知细胞类型
  missing_types <- setdiff(sample_celltypes, names(celltype_colors_global))
  if (length(missing_types) > 0) {
    warning(sprintf("样本 %s 包含未在全局颜色方案中的细胞类型: %s",
                   sample_id, paste(missing_types, collapse = ", ")))
    # 为缺失类型分配灰色
    extra_colors <- rep("#CCCCCC", length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors_global <- c(celltype_colors_global, extra_colors)
  }
  
  # ✅ 使用全局颜色（保留所有细胞类型，即使当前样本没有）
  celltype_colors <- celltype_colors_global
  all_celltypes <- names(celltype_colors)  # 全局所有细胞类型
  
  cat(sprintf("   📊 当前样本细胞类型: %d 个 (全局: %d 个)\n", 
              length(sample_celltypes), length(all_celltypes)))
  
  # 只打印当前样本的颜色映射
  if (length(sample_celltypes) <= 10) {
    cat("   📊 当前样本细胞类型颜色:\n")
    for (ct in sample_celltypes) {
      cat(sprintf("      • %-25s → %s\n", ct, celltype_colors[ct]))
    }
  } else {
    cat("   📊 当前样本细胞类型颜色（前10个）:\n")
    for (i in 1:min(10, length(sample_celltypes))) {
      ct <- sample_celltypes[i]
      cat(sprintf("      • %-25s → %s\n", ct, celltype_colors[ct]))
    }
    cat(sprintf("      ... 还有 %d 个细胞类型\n", length(sample_celltypes) - 10))
  }
  
  # ========================================
  # 3. 坐标范围
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
  # 4. 准备等高线数据
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
  # 5. 准备细胞数据
  # ========================================
  
  df_filtered <- df %>% 
    filter(!is.na(density_zone)) %>%
    mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      # ✅ 使用全局所有细胞类型作为 factor levels
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
  
  # 统一的图例样式参数
  legend_key_width <- 1.0
  legend_key_height <- 0.7
  legend_text_size <- 10
  legend_title_size <- 11
  
  # ========================================
  # 6. 绘图
  # ========================================
  
  p <- ggplot() +
    # ========================================
    # Layer 1: 细胞类型（底层）- 使用 FILL
    # ========================================
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, fill = celltype_clean),
      width = square_size,
      height = square_size,
      color = NA,
      alpha = 1
    ) +
    scale_fill_manual(
      values = celltype_colors,  # ✅ 使用全局颜色
      name = "Cell Type",
      breaks = all_celltypes,  # ✅ 显示所有细胞类型
      drop = FALSE,  # ✅ 不丢弃未使用的级别
      na.value = "gray50",
      guide = guide_legend(
        order = 2,
        override.aes = list(
          alpha = 1,
          color = NA
        ),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keywidth = unit(legend_key_width, "cm"),
        keyheight = unit(legend_key_height, "cm")
      )
    ) +
    
    # ========================================
    # 新的scale用于density zones
    # ========================================
    ggnewscale::new_scale_fill() +
    
    # Layer 2: Zone填充（半透明覆盖层）
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.3,
      interpolate = TRUE
    ) +
    scale_fill_manual(
      values = zone_colors_global,  # ✅ 使用全局区域颜色
      labels = zone_labels,
      name = "Density Zones\n(Zone_0 = Core Red → Zone_N = Outer Blue)",
      breaks = zone_levels,
      na.value = "transparent",
      drop = FALSE,
      guide = guide_legend(
        order = 1,
        override.aes = list(
          alpha = 0.8,
          color = "gray40",
          linewidth = 0.3
        ),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keywidth = unit(legend_key_width, "cm"),
        keyheight = unit(legend_key_height, "cm")
      )
    ) +
    
    # ========================================
    # 为等高线准备新的color scale
    # ========================================
    ggnewscale::new_scale_color() +
    
    # Layer 3: 等高线边界
    geom_contour(
      data = contour_data,
      aes(x = col, y = row, z = density_norm, color = after_stat(level)),
      breaks = contour_breaks,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    scale_color_gradientn(
      colors = contour_colors,
      limits = c(min(contour_breaks), max(contour_breaks)),
      guide = "none"
    ) +
    
    # ========================================
    # 坐标系统
    # ========================================
    scale_x_continuous(limits = col_limits, expand = c(0, 0)) +
    scale_y_reverse(limits = rev(row_limits), expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = col_limits, ylim = rev(row_limits), clip = "off") +
    
    # ========================================
    # 标题和主题
    # ========================================
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = sprintf("Bottom = Cell types | Middle = Density zones (α=0.3) | Top = %d contour lines", 
                        length(contour_breaks)),
      x = NULL, y = NULL
    ) +
    
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = unit(0.6, "cm"),
      legend.title = element_text(size = legend_title_size, face = "bold", hjust = 0, margin = margin(b = 8)),
      legend.text = element_text(size = legend_text_size, lineheight = 1.2, margin = margin(l = 3, r = 5, t = 2, b = 2)),
      legend.key = element_rect(color = "gray60", fill = NA, linewidth = 0.3),
      legend.key.spacing.y = unit(0.2, "cm"),
      legend.background = element_rect(fill = "white", color = "gray50", linewidth = 0.5),
      legend.margin = margin(10, 10, 10, 10),
      plot.margin = margin(15, 20, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# ========================================
# 辅助函数（向后兼容）
# ========================================

#' 生成等高线颜色（紫色渐变）
#' 
#' @param n_contours 等高线数量
#' @return 颜色向量
#'
get_contour_colors <- function(n_contours) {
  colorRampPalette(c(
    "#542788",  # 深紫
    "#8073AC",  # 中紫
    "#B2ABD2",  # 浅紫
    "#D8DAEB"   # 淡紫
  ))(n_contours)
}


# ========================================
# ⚠️ 以下函数已弃用，请使用 01_color_schemes.R 中的函数
# ========================================

#' @deprecated 请使用 get_zone_colors() from 01_color_schemes.R
get_zone_colors <- function(n_zones) {
  warning("get_zone_colors() 已弃用，请使用 01_color_schemes.R 中的版本")
  
  colorRampPalette(c(
    "#B2182B",  # 深红（核心高密度）
    "#EF8A62",  # 浅红
    "#FDDBC7",  # 粉色
    "#F7F7F7",  # 白色（中间）
    "#D1E5F0",  # 浅蓝
    "#67A9CF",  # 蓝色
    "#2166AC"   # 深蓝（外围低密度）
  ))(n_zones)
}

#' @deprecated 请使用 get_celltype_colors() from 01_color_schemes.R
get_celltype_colors <- function(celltypes) {
  warning("get_celltype_colors() 已弃用，请使用 01_color_schemes.R 中的版本")
  
  require(RColorBrewer)
  n <- length(celltypes)
  
  if (n <= 3) {
    colors <- brewer.pal(3, "Set2")[1:n]
  } else if (n <= 12) {
    colors <- brewer.pal(n, "Set3")
  } else {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n)
  }
  
  names(colors) <- celltypes
  return(colors)
}

cat("✅ 03_plot_overlay.R 已加载（全局统一配色版）\n")