# ===================================================================
# 04_plot_composition.R (全局统一配色版)
# 区域组成柱状图绘制
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 绘制zone组成柱状图（全局统一配色版）
#'
#' @param zone_composition zone组成数据框
#' @param sample_id 样本ID
#' @param CONFIG 配置列表（必须包含 CONFIG$colors）
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
  
  # ✅ 使用全局颜色（所有样本一致）
  celltype_colors_global <- CONFIG$colors$celltype
  zone_colors_global <- CONFIG$colors$density_zone
  
  # 全局所有细胞类型
  all_celltypes <- names(celltype_colors_global)
  
  # 当前样本的细胞类型
  sample_celltypes <- unique(zone_composition$celltype_clean)
  
  # 检查未知类型
  missing_types <- setdiff(sample_celltypes, all_celltypes)
  if (length(missing_types) > 0) {
    warning(sprintf("样本 %s 包含未知细胞类型: %s", 
                   sample_id, paste(missing_types, collapse = ", ")))
  }
  
  # ========================================
  # 3. 准备数据
  # ========================================
  
  # 确保zone按 Zone_0, Zone_1, ... 排序
  n_zones <- length(unique(zone_composition$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  zone_composition <- zone_composition %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      # ✅ 使用全局所有细胞类型作为 factor levels
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  cat(sprintf("   📊 绘制组成图: %d zones, %d celltypes (全局: %d)\n",
              n_zones, length(sample_celltypes), length(all_celltypes)))
  
  # ========================================
  # 4. 图1：细胞类型组成堆叠柱状图
  # ========================================
  
  p1 <- ggplot(zone_composition, aes(x = density_zone, y = percentage, fill = celltype_clean)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = celltype_colors_global,  # ✅ 使用全局颜色
      name = "Cell Type",
      breaks = all_celltypes,  # ✅ 显示所有细胞类型
      drop = FALSE  # ✅ 不丢弃未使用的级别
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = sprintf("Cell Type Composition by Density Zone - %s", sample_id),
      subtitle = sprintf("%d cell types (global colors applied)", length(all_celltypes)),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "cm"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    )
  
  # ========================================
  # 5. 图2：Zone的spot数量
  # ========================================
  
  zone_totals <- zone_composition %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(total = sum(count), .groups = "drop") %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  p2 <- ggplot(zone_totals, aes(x = density_zone, y = total, fill = density_zone)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.5) +
    geom_text(aes(label = total), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(
      values = zone_colors_global,  # ✅ 使用全局区域颜色
      guide = "none"
    ) +
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
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10, face = "bold"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    )
  
  # ========================================
  # 6. 合并两个图
  # ========================================
  
  p_combined <- p1 / p2 + 
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      caption = sprintf("Global color scheme applied across all samples | Sample: %s", sample_id),
      theme = theme(
        plot.caption = element_text(size = 8, color = "gray40", hjust = 1, margin = margin(t = 10))
      )
    )
  
  return(p_combined)
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

#' @deprecated 请使用 get_celltype_colors() from 01_color_schemes.R
get_celltype_colors <- function(celltypes) {
  warning("get_celltype_colors() 已弃用，请使用 01_color_schemes.R 中的版本")
  
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

cat("✅ 04_plot_composition.R 已加载（全局统一配色版）\n")