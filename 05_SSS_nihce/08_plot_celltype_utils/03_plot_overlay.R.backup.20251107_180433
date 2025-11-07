# ==================================================================
# 03_plot_overlay.R - 细胞类型+密度叠加图
# ==================================================================

#' 绘制细胞类型和密度区域叠加图
plot_celltype_density_overlay <- function(
    df, 
    density_data, 
    sample_id, 
    CONFIG) {
  
  # 验证配置
  if (is.null(CONFIG$colors$celltype) || 
      is.null(CONFIG$colors$density_zone)) {
    stop("全局颜色方案未初始化")
  }
  
  # 准备颜色和分级
  celltype_colors <- CONFIG$colors$celltype
  zone_colors <- CONFIG$colors$density_zone
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  # 清理细胞类型
  df$celltype_clean <- as.character(df$celltype_clean)
  df$celltype_clean[
    is.na(df$celltype_clean) | df$celltype_clean == ""
  ] <- "Unknown"
  
  # 处理缺失颜色
  sample_celltypes <- unique(df$celltype_clean)
  missing_types <- setdiff(
    sample_celltypes, 
    names(celltype_colors)
  )
  
  if (length(missing_types) > 0) {
    extra_colors <- rep("#CCCCCC", length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors <- c(celltype_colors, extra_colors)
  }
  
  all_celltypes <- names(celltype_colors)
  
  # 坐标范围
  if (!is.null(density_data$col_range_expanded)) {
    col_limits <- density_data$col_range_expanded
    row_limits <- density_data$row_range_expanded
  } else {
    col_range <- density_data$col_range
    row_range <- density_data$row_range
    margin <- 0.1
    col_limits <- c(
      col_range[1] - diff(col_range) * margin,
      col_range[2] + diff(col_range) * margin
    )
    row_limits <- c(
      row_range[1] - diff(row_range) * margin,
      row_range[2] + diff(row_range) * margin
    )
  }
  
  # 准备等高线
  zone_density_ranges <- density_data$grid %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(
      min_density = min(density_norm, na.rm = TRUE),
      max_density = max(density_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels)
    )
  
  zone_labels <- zone_density_ranges %>%
    dplyr::mutate(
      zone_label = sprintf(
        "%s (%.2f-%.2f)", 
        density_zone, 
        min_density, 
        max_density
      )
    ) %>%
    dplyr::pull(zone_label)
  
  names(zone_labels) <- as.character(
    zone_density_ranges$density_zone
  )
  
  contour_data <- density_data$grid %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels)
    )
  
  # 过滤细胞数据
  df_filtered <- df %>% 
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  # 计算细胞大小
  coords_sample <- if (nrow(df_filtered) > 10000) {
    df_filtered[sample(nrow(df_filtered), 10000), c("col", "row")]
  } else {
    df_filtered[, c("col", "row")]
  }
  
  square_size <- median(
    RANN::nn2(coords_sample, k = 2)$nn.dists[, 2], 
    na.rm = TRUE
  )
  
  # 等高线颜色
  contour_colors <- grDevices::colorRampPalette(
    c("#542788", "#8073AC", "#B2ABD2", "#D8DAEB")
  )(length(density_data$equal_breaks))
  
  # 绘图
  p <- ggplot2::ggplot() +
    # 细胞类型层
    ggplot2::geom_tile(
      data = df_filtered,
      ggplot2::aes(x = col, y = row, fill = celltype_clean),
      width = square_size,
      height = square_size,
      color = NA
    ) +
    ggplot2::scale_fill_manual(
      values = celltype_colors,
      name = "Cell Type",
      breaks = all_celltypes,
      drop = FALSE,
      na.value = "gray50",
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(alpha = 1, color = NA),
        title.position = "top",
        ncol = 1,
        keywidth = ggplot2::unit(1, "cm"),
        keyheight = ggplot2::unit(0.7, "cm")
      )
    ) +
    ggnewscale::new_scale_fill() +
    
    # 密度区域层
    ggplot2::geom_raster(
      data = contour_data,
      ggplot2::aes(x = col, y = row, fill = density_zone),
      alpha = 0.3,
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones\n(Zone_0 = Core → Zone_N = Outer)",
      breaks = zone_levels,
      na.value = "transparent",
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = 0.8,
          color = "gray40",
          linewidth = 0.3
        ),
        title.position = "top",
        ncol = 1,
        keywidth = ggplot2::unit(1, "cm"),
        keyheight = ggplot2::unit(0.7, "cm")
      )
    ) +
    ggnewscale::new_scale_color() +
    
    # 等高线层
    ggplot2::geom_contour(
      data = contour_data,
      ggplot2::aes(
        x = col, 
        y = row, 
        z = density_norm, 
        color = ggplot2::after_stat(level)
      ),
      breaks = density_data$equal_breaks,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    ggplot2::scale_color_gradientn(
      colors = contour_colors,
      guide = "none"
    ) +
    
    # 坐标系
    ggplot2::scale_x_continuous(
      limits = col_limits, 
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = rev(row_limits), 
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(
      ratio = 1, 
      xlim = col_limits, 
      ylim = rev(row_limits), 
      clip = "off"
    ) +
    
    # 主题
    ggplot2::labs(
      title = sprintf(
        "Cell Type Distribution in Density Zones - %s", 
        sample_id
      ),
      subtitle = sprintf(
        "Cell types + Density zones (α=0.3) + %d contours", 
        length(density_data$equal_breaks)
      ),
      x = NULL, 
      y = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 16, face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = 9, color = "gray30"
      ),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(0.6, "cm"),
      legend.title = ggplot2::element_text(
        size = 11, face = "bold"
      ),
      legend.text = ggplot2::element_text(size = 10),
      legend.key = ggplot2::element_rect(
        color = "gray60", fill = NA
      ),
      legend.background = ggplot2::element_rect(
        fill = "white", color = "gray50"
      ),
      plot.background = ggplot2::element_rect(
        fill = "white", color = NA
      )
    )
  
  return(p)
}

cat("✅ 03_plot_overlay.R 已加载\n")