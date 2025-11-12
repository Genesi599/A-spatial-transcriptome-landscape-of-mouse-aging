# ===================================================================
# 06_plot_combined.R (全局统一配色版 - 叠加柱状图)
# 综合分析图绘制
# Author: Assistant
# Date: 2025-11-12
# ===================================================================

#' 绘制综合分析图（箱线图 + 叠加柱状图）- 全局统一配色版
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表（必须包含 CONFIG$colors）
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
  
  if (is.null(CONFIG$colors) || is.null(CONFIG$colors$celltype)) {
    stop("❌ 全局颜色方案未初始化")
  }
  
  if (is.null(CONFIG$colors$density_zone)) {
    stop("❌ 密度区域颜色未初始化")
  }
  
  celltype_colors_global <- CONFIG$colors$celltype
  zone_colors_global <- CONFIG$colors$density_zone
  
  all_celltypes <- names(celltype_colors_global)
  
  n_zones <- length(unique(combined_data$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf("\n📊 绘制综合分析图...\n"))
  cat(sprintf("   样本数: %d\n", 
              length(unique(combined_data$sample))))
  cat(sprintf("   细胞类型: %d (全局)\n", 
              length(all_celltypes)))
  cat(sprintf("   密度区域: %d\n", n_zones))
  
  combined_data <- combined_data %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, 
                              levels = all_celltypes),
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    )
  
  p1 <- ggplot(combined_data, 
               aes(x = density_zone, y = percentage, 
                   fill = density_zone)) +
    geom_boxplot(
      alpha = 0.8, 
      outlier.shape = 16, 
      outlier.size = 2, 
      color = "gray20", 
      linewidth = 0.6
    ) +
    scale_fill_manual(
      values = zone_colors_global,
      guide = "none"
    ) +
    facet_wrap(
      ~celltype_clean, 
      scales = "free_y", 
      ncol = 4,
      drop = FALSE
    ) +
    labs(
      x = "Density Zone",
      y = "Percentage (%)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "gray95", 
                                      color = "gray50", 
                                      linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_rect(color = "gray50", linewidth = 0.5)
    )
  
  # 修复后的叠加柱状图数据准备
  stacked_data <- combined_data %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(
      mean_pct = mean(percentage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(density_zone) %>%
    dplyr::mutate(
      zone_total = sum(mean_pct, na.rm = TRUE),
      normalized_pct = if_else(
        zone_total > 0,
        mean_pct / zone_total * 100,  # 强制归一化到100%
        0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-zone_total, -mean_pct) %>%
    dplyr::rename(mean_pct = normalized_pct)

  complete_stacked_grid <- expand.grid(
    celltype_clean = all_celltypes,
    density_zone = zone_levels,
    stringsAsFactors = FALSE
  )

  stacked_data <- complete_stacked_grid %>%
    dplyr::left_join(stacked_data, 
                    by = c("celltype_clean", "density_zone")) %>%
    dplyr::mutate(
      mean_pct = ifelse(is.na(mean_pct), 0, mean_pct),
      celltype_clean = factor(celltype_clean, levels = all_celltypes),
      density_zone = factor(density_zone, levels = zone_levels)
    )

  # 验证归一化结果
  verify_totals <- stacked_data %>%
    group_by(density_zone) %>%
    summarise(total = sum(mean_pct))

  cat("\n✅ 归一化后的Zone总和:\n")
  print(verify_totals)
  
  p2 <- ggplot(stacked_data, 
               aes(x = density_zone, y = mean_pct, 
                   fill = celltype_clean)) +
    geom_col(position = "stack", width = 0.75, alpha = 0.95) +
    scale_fill_manual(
      values = celltype_colors_global,
      name = "Cell Type",
      breaks = all_celltypes,
      drop = FALSE
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      x = "Density Zone",
      y = "Composition (%)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(linewidth = 0.6, color = "gray20"),
      axis.ticks = element_line(linewidth = 0.6, color = "gray20"),
      legend.position = "right",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.height = unit(0.7, "cm"),
      legend.key.width = unit(0.6, "cm"),
      panel.grid.major.y = element_line(color = "gray85", 
                                       linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  p_combined <- p1 / p2 + 
    plot_layout(heights = c(2.5, 1))
  
  cat("   ✅ 综合分析图绘制完成\n")
  
  return(list(
    combined = p_combined,
    boxplot = p1,
    barplot = p2
  ))
}

get_zone_colors <- function(n_zones) {
  warning(paste0(
    "get_zone_colors() 已弃用，",
    "请使用 01_color_schemes.R 中的版本"
  ))
  
  colorRampPalette(c(
    "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
  ))(n_zones)
}

get_celltype_colors <- function(celltypes) {
  warning(paste0(
    "get_celltype_colors() 已弃用，",
    "请使用 01_color_schemes.R 中的版本"
  ))
  
  require(RColorBrewer)
  n <- length(celltypes)
  
  if (n <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n), "Set2")
  } else if (n <= 12) {
    colors <- RColorBrewer::brewer.pal(n, "Set3")
  } else {
    colors <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3")
    )[1:n]
  }
  
  names(colors) <- celltypes
  return(colors)
}

cat("✅ 06_plot_combined.R 已加载（叠加柱状图版）\n")