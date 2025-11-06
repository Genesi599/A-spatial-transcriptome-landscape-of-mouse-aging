# ===================================================================
# 08_plot_celltype.R
# 细胞类型 + 等高线分析完整工作流（修复版）
# Author: Assistant
# Date: 2025-11-05
# ===================================================================

# ===================================================================
# 辅助函数 0：统一的颜色方案
# ===================================================================

#' 生成统一的zone颜色方案（密度从高到低）
#'
#' @param n_zones zone数量
#' @return 命名的颜色向量
get_zone_colors <- function(n_zones = 5) {
  # Zone_0 = 核心（深红），Zone_4 = 外围（浅黄）
  zone_colors <- colorRampPalette(c(
    "#800026",  # 深红（核心 Zone_0）
    "#bd0026",
    "#e31a1c",
    "#fc4e2a",
    "#fd8d3c",
    "#feb24c",
    "#fed976",
    "#ffeda0",
    "#ffffcc"   # 浅黄（外围 Zone_n-1）
  ))(n_zones)
  
  zone_names <- sprintf("Zone_%d", 0:(n_zones - 1))
  names(zone_colors) <- zone_names
  
  return(zone_colors)
}


#' 生成统一的细胞类型颜色
#'
#' @param celltypes 细胞类型向量
#' @return 命名的颜色向量
get_celltype_colors <- function(celltypes) {
  n_celltypes <- length(celltypes)
  
  if (n_celltypes <= 8) {
    colors <- brewer.pal(max(3, n_celltypes), "Set2")
  } else if (n_celltypes <= 12) {
    colors <- brewer.pal(n_celltypes, "Set3")
  } else {
    colors <- c(
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Set2"),
      brewer.pal(12, "Set3")
    )[1:n_celltypes]
  }
  
  names(colors) <- celltypes
  return(colors)
}


# ===================================================================
# 主函数：细胞类型等高线分析
# ===================================================================

#' 细胞类型 + Clock Gene Niche 等高线综合分析
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 要分析的样本列表
#' @param CONFIG 配置列表
#' @param density_bins 等高线分级数量，默认 5（对应5个区域）
#' @param celltype_col 细胞类型列名，默认 "celltype"
#' @param plot_overlay 是否绘制叠加图，默认 TRUE
#' @param plot_composition 是否绘制组成图，默认 TRUE
#' @param plot_heatmap 是否绘制热图，默认 TRUE
#' @param plot_combined 是否绘制合并分析图，默认 TRUE
#'
#' @return 返回统计数据列表
#'
#' @examples
#' result <- analyze_celltype_niche(seurat_obj, samples_to_plot, CONFIG)
#'
analyze_celltype_niche <- function(
    seurat_obj,
    samples_to_plot,
    CONFIG,
    density_bins = 5,
    celltype_col = "celltype",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE
) {
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("🧬 细胞类型 + Clock Gene Niche 等高线分析\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  required_cols <- c("ClockGene_High", "orig.ident", celltype_col)
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 创建输出目录
  output_dirs <- c(
    CONFIG$dirs$overlay,
    CONFIG$dirs$celltype,
    CONFIG$dirs$composition,
    CONFIG$dirs$heatmaps,
    CONFIG$dirs$combined
  )
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 验证样本
  available_samples <- unique(seurat_obj$orig.ident)
  samples_to_plot <- intersect(samples_to_plot, available_samples)
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本可分析")
  }
  
  cat(sprintf("✅ 将分析 %d 个样本\n", length(samples_to_plot)))
  cat(sprintf("✅ 等高线分为 %d 个区域 (Zone_0=核心高密度, Zone_%d=外围低密度)\n", 
              density_bins, density_bins - 1))
  
  # ========================================
  # 2. 初始化结果容器
  # ========================================
  all_sample_stats <- list()
  combined_data <- data.frame()
  
  # ========================================
  # 3. 逐样本分析
  # ========================================
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    cat(sprintf("\n[%d/%d] 📊 分析样本: %s\n", i, length(samples_to_plot), sample_id))
    cat(rep("-", 80), "\n", sep = "")
    
    tryCatch({
      # -------------------------------
      # 3.1 提取样本数据
      # -------------------------------
      seurat_subset <- subset(seurat_obj, subset = orig.ident == sample_id)
      
      if (ncol(seurat_subset) == 0) {
        warning(sprintf("样本 %s 无数据，跳过", sample_id))
        next
      }
      
      # 获取坐标
      coords <- GetTissueCoordinates(
        seurat_subset,
        cols = c("row", "col"),
        scale = NULL
      )
      
      # 合并数据
      df <- seurat_subset@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(coords %>% rownames_to_column("barcode"), by = "barcode") %>%
        filter(!is.na(col), !is.na(row))
      
      # 检查细胞类型
      df$celltype_clean <- as.character(df[[celltype_col]])
      df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
      
      cat(sprintf("   ✅ 有效spots: %d\n", nrow(df)))
      cat(sprintf("   ✅ 高表达spots: %d (%.2f%%)\n", 
                  sum(df$ClockGene_High), 
                  100 * mean(df$ClockGene_High)))
      
      # -------------------------------
      # 3.2 计算密度并分级
      # -------------------------------
      density_data <- calculate_density_zones(
        df = df,
        density_bins = density_bins,
        expand_margin = CONFIG$plot$expand_margin %||% 0.05
      )
      
      if (is.null(density_data)) {
        warning(sprintf("样本 %s 密度计算失败，跳过", sample_id))
        next
      }
      
      # 合并密度信息到df
      df <- df %>%
        left_join(
          density_data$spot_zones %>% select(col, row, density_zone, density_value),
          by = c("col", "row")
        )
      
      # 检查NA情况
      n_na <- sum(is.na(df$density_zone))
      if (n_na > 0) {
        cat(sprintf("   ⚠️  警告: %d 个spots未分配到zone (%.2f%%)\n", 
                    n_na, 100 * n_na / nrow(df)))
      }
      
      # 统计每个区域的细胞类型组成
      zone_composition <- df %>%
        filter(!is.na(density_zone)) %>%
        group_by(density_zone, celltype_clean) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(density_zone) %>%
        mutate(
          total = sum(count),
          percentage = 100 * count / total
        ) %>%
        ungroup() %>%
        mutate(sample = sample_id)
      
      cat(sprintf("   ✅ 密度分区完成，共 %d 个区域\n", 
                  length(unique(zone_composition$density_zone))))
      
      # 打印每个zone的统计（修复：使用已合并density_value的df）
      zone_stats <- df %>%
        filter(!is.na(density_zone)) %>%
        group_by(density_zone) %>%
        summarise(
          n_spots = n(),
          mean_density = mean(density_value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(density_zone)
      
      cat("   Zone统计:\n")
      for (j in 1:nrow(zone_stats)) {
        cat(sprintf("     %s: %d spots (mean density: %.3f)\n",
                    zone_stats$density_zone[j],
                    zone_stats$n_spots[j],
                    zone_stats$mean_density[j]))
      }
      
      # 保存到总体数据
      all_sample_stats[[sample_id]] <- zone_composition
      combined_data <- bind_rows(combined_data, zone_composition)
      
      # -------------------------------
      # 3.3 绘制叠加图
      # -------------------------------
      if (plot_overlay) {
        p_overlay <- plot_celltype_density_overlay(
          df = df,
          density_data = density_data,
          sample_id = sample_id,
          CONFIG = CONFIG
        )
        
        safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
        ggsave(
          file.path(CONFIG$dirs$overlay, sprintf("celltype_overlay_%s.pdf", safe_name)),
          plot = p_overlay,
          width = 12, height = 10, dpi = CONFIG$plot$dpi %||% 300
        )
        cat("   ✅ 保存叠加图\n")
      }
      
      # -------------------------------
      # 3.4 绘制组成图
      # -------------------------------
      if (plot_composition) {
        p_comp <- plot_zone_composition(
          zone_composition = zone_composition,
          sample_id = sample_id,
          CONFIG = CONFIG
        )
        
        safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
        ggsave(
          file.path(CONFIG$dirs$composition, sprintf("composition_%s.pdf", safe_name)),
          plot = p_comp,
          width = 12, height = 6, dpi = CONFIG$plot$dpi %||% 300
        )
        cat("   ✅ 保存组成图\n")
      }
      
    }, error = function(e) {
      cat(sprintf("   ❌ 错误: %s\n", e$message))
    })
  }
  
  # ========================================
  # 4. 合并所有样本的统计分析
  # ========================================
  if (nrow(combined_data) > 0) {
    cat("\n")
    cat(rep("=", 80), "\n", sep = "")
    cat("📈 合并所有样本进行统计分析\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    # -------------------------------
    # 4.1 绘制热图
    # -------------------------------
    if (plot_heatmap) {
      p_heatmap <- plot_combined_heatmap(
        combined_data = combined_data,
        CONFIG = CONFIG
      )
      
      ggsave(
        file.path(CONFIG$dirs$heatmaps, "celltype_heatmap_all_samples.pdf"),
        plot = p_heatmap,
        width = 14, height = 10, dpi = CONFIG$plot$dpi %||% 300
      )
      cat("✅ 保存合并热图\n")
    }
    
    # -------------------------------
    # 4.2 绘制综合分析图
    # -------------------------------
    if (plot_combined) {
      p_combined <- plot_combined_analysis(
        combined_data = combined_data,
        CONFIG = CONFIG
      )
      
      ggsave(
        file.path(CONFIG$dirs$combined, "combined_analysis.pdf"),
        plot = p_combined,
        width = 16, height = 12, dpi = CONFIG$plot$dpi %||% 300
      )
      cat("✅ 保存综合分析图\n")
    }
    
    # -------------------------------
    # 4.3 保存统计数据
    # -------------------------------
    write.csv(
      combined_data,
      file.path(CONFIG$dirs$composition, "celltype_composition_all_samples.csv"),
      row.names = FALSE
    )
    cat("✅ 保存统计数据 CSV\n")
    
    # -------------------------------
    # 4.4 统计摘要
    # -------------------------------
    summary_stats <- generate_summary_statistics(combined_data)
    write.csv(
      summary_stats,
      file.path(CONFIG$dirs$composition, "summary_statistics.csv"),
      row.names = FALSE
    )
    cat("✅ 保存统计摘要\n")
  }
  
  # ========================================
  # 5. 返回结果
  # ========================================
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("✅ 分析完成！\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  invisible(list(
    sample_stats = all_sample_stats,
    combined_data = combined_data,
    summary_stats = if(exists("summary_stats")) summary_stats else NULL
  ))
}


# ===================================================================
# 辅助函数 1：计算密度分区（修改版：反转zone序号）
# ===================================================================

calculate_density_zones <- function(df, density_bins = 5, expand_margin = 0.05) {
  
  # 只使用高表达点计算密度
  high_points <- df %>% filter(ClockGene_High)
  
  if (nrow(high_points) < 10) {
    warning("高表达点数量不足（< 10），无法计算密度")
    return(NULL)
  }
  
  # 计算坐标范围
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  
  col_expand <- diff(col_range) * expand_margin
  row_expand <- diff(row_range) * expand_margin
  
  col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
  row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
  
  # 使用 MASS::kde2d 计算密度
  kde_result <- tryCatch({
    MASS::kde2d(
      x = high_points$col,
      y = high_points$row,
      n = 200,
      lims = c(col_limits, row_limits)
    )
  }, error = function(e) {
    warning(sprintf("密度估计失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(kde_result)) return(NULL)
  
  # 转换为 data frame
  density_df <- expand.grid(
    col = kde_result$x,
    row = kde_result$y
  )
  density_df$density <- as.vector(kde_result$z)
  
  # 归一化密度
  density_df$density_norm <- density_df$density / max(density_df$density, na.rm = TRUE)
  
  # 分级（反转：Zone_0 = 最高密度核心，Zone_n-1 = 最低密度外围）
  density_df$density_zone <- cut(
    density_df$density_norm,
    breaks = seq(0, 1, length.out = density_bins + 1),
    labels = sprintf("Zone_%d", (density_bins - 1):0),  # 反转标签顺序
    include.lowest = TRUE,
    right = TRUE
  )
  
  # 为每个spot分配最近的密度区域（使用更可靠的方法）
  # 创建一个密度网格的索引
  spot_zones <- df %>%
    select(col, row) %>%
    mutate(
      # 找到最近的网格点
      col_idx = sapply(col, function(x) which.min(abs(kde_result$x - x))),
      row_idx = sapply(row, function(y) which.min(abs(kde_result$y - y)))
    ) %>%
    rowwise() %>%
    mutate(
      # 根据索引获取密度zone
      grid_idx = (col_idx - 1) * length(kde_result$y) + row_idx,
      density_zone = density_df$density_zone[grid_idx],
      density_value = density_df$density_norm[grid_idx]
    ) %>%
    ungroup() %>%
    select(col, row, density_zone, density_value)
  
  # 检查是否有NA，如果有，使用最近邻方法填充
  if (any(is.na(spot_zones$density_zone))) {
    # 对有NA的点，使用KNN方法分配
    na_spots <- which(is.na(spot_zones$density_zone))
    
    for (idx in na_spots) {
      spot_col <- spot_zones$col[idx]
      spot_row <- spot_zones$row[idx]
      
      # 计算到所有网格点的距离
      distances <- sqrt((density_df$col - spot_col)^2 + (density_df$row - spot_row)^2)
      
      # 找到最近的非NA点
      valid_idx <- which(!is.na(density_df$density_zone))
      nearest_valid <- valid_idx[which.min(distances[valid_idx])]
      
      spot_zones$density_zone[idx] <- density_df$density_zone[nearest_valid]
      spot_zones$density_value[idx] <- density_df$density_norm[nearest_valid]
    }
  }
  
  return(list(
    grid = density_df,
    spot_zones = spot_zones,
    kde_result = kde_result
  ))
}


# ===================================================================
# 辅助函数 2：绘制细胞类型+密度叠加图（等高线和Zone都在最上层）
# ===================================================================

plot_celltype_density_overlay <- function(df, density_data, sample_id, CONFIG) {
  
  # 获取zone信息
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sort(unique(as.character(density_data$grid$density_zone)))
  
  # 使用统一的颜色方案
  zone_colors <- get_zone_colors(n_zones)
  celltype_colors <- get_celltype_colors(unique(df$celltype_clean))
  
  # 坐标范围
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  expand <- CONFIG$plot$expand_margin %||% 0.05
  
  col_limits <- col_range + c(-1, 1) * diff(col_range) * expand
  row_limits <- row_range + c(-1, 1) * diff(row_range) * expand
  
  # 计算每个zone的密度范围
  zone_density_ranges <- density_data$grid %>%
    group_by(density_zone) %>%
    summarise(
      min_density = min(density_norm, na.rm = TRUE),
      max_density = max(density_norm, na.rm = TRUE),
      mean_density = mean(density_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_density))
  
  # 为zone创建带密度信息的标签
  zone_labels <- zone_density_ranges %>%
    mutate(
      zone_label = sprintf("%s\n(%.2f-%.2f)", 
                          density_zone, 
                          min_density, 
                          max_density)
    ) %>%
    pull(zone_label)
  
  names(zone_labels) <- zone_density_ranges$density_zone
  
  # 计算等高线的具体数值（对应zone边界）
  contour_breaks <- seq(0, 1, length.out = n_zones + 1)
  
  # 为等高线创建颜色映射数据
  contour_data <- density_data$grid
  
  # =============================================
  # 自动计算正方形大小（每个点一个正方形）
  # =============================================
  df_filtered <- df %>% filter(!is.na(density_zone))
  
  # 计算最近邻距离来确定合适的正方形大小
  # 使用采样以加快计算（如果细胞太多）
  if (nrow(df_filtered) > 10000) {
    sample_idx <- sample(nrow(df_filtered), 10000)
    coords_sample <- df_filtered[sample_idx, c("col", "row")]
  } else {
    coords_sample <- df_filtered[, c("col", "row")]
  }
  
  # 计算最近邻距离
  nn_dist <- RANN::nn2(coords_sample, k = 2)$nn.dists[, 2]  # 第2列是最近邻距离
  median_dist <- median(nn_dist, na.rm = TRUE)
  
  # 正方形的边长应该等于最近邻距离，使其刚好无缝铺满
  square_size <- median_dist * 1.0  # 可以微调这个系数（0.95-1.05）
  
  # 绘图
  p <- ggplot() +
    # =============================================
    # 1. 细胞类型正方形（底层，每个点一个正方形）
    # =============================================
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, fill = celltype_clean),
      width = square_size,
      height = square_size,
      alpha = 0.85,
      color = NA  # 无描边
    ) +
    scale_fill_manual(
      values = celltype_colors,
      name = "Cell Type",
      guide = guide_legend(
        order = 2,
        override.aes = list(size = 4, alpha = 1),
        title.position = "top",
        title.hjust = 0.5,
        ncol = 1
      )
    ) +
    new_scale_fill() +
    
    # =============================================
    # 2. Zone区域填充（上层，半透明覆盖）
    # =============================================
    geom_tile(
      data = density_data$grid,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.25  # 更透明一些，可以看到下面的细胞类型
    ) +
    scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones & Contour Lines\n(Normalized Range, Zone_0=Core)",
      breaks = zone_levels,
      guide = guide_legend(
        order = 1,
        override.aes = list(alpha = 0.7, size = 0),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        label.hjust = 0,
        keywidth = unit(1.2, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    
    # =============================================
    # 3. 等高线（最上层，无描边）
    # =============================================
    {
      contour_layers <- list()
      for (i in 1:length(contour_breaks)) {
        # 计算该等高线对应的zone
        if (i == 1) {
          zone_idx <- n_zones - 1
        } else {
          zone_idx <- n_zones - i + 1
        }
        
        zone_name <- sprintf("Zone_%d", zone_idx)
        zone_color <- zone_colors[zone_name]
        
        contour_layers[[i]] <- geom_contour(
          data = contour_data,
          aes(x = col, y = row, z = density_norm),
          breaks = contour_breaks[i],
          color = zone_color,
          linewidth = 1.5,
          alpha = 0.9
        )
      }
      contour_layers
    } +
    
    # =============================================
    # 坐标和主题
    # =============================================
    scale_x_continuous(limits = col_limits, expand = expansion(mult = 0.02)) +
    scale_y_reverse(limits = rev(row_limits), expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1) +
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = "Bottom layer = Cell types (squares) | Middle layer = Density zones (transparent) | Top layer = Contour lines (Zone_0=Core/High)"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 8.5, lineheight = 1.2),
      legend.key = element_rect(color = "gray70", linewidth = 0.3),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p)
}


# ===================================================================
# 辅助函数 3：绘制区域组成柱状图（统一配色版）
# ===================================================================

plot_zone_composition <- function(zone_composition, sample_id, CONFIG) {
  
  # 使用统一的颜色方案
  n_zones <- length(unique(zone_composition$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  celltype_colors <- get_celltype_colors(unique(zone_composition$celltype_clean))
  
  # 确保zone按顺序排列
  zone_composition <- zone_composition %>%
    mutate(density_zone = factor(density_zone, levels = names(zone_colors)))
  
  # 图1：细胞类型组成堆叠柱状图
  p1 <- ggplot(zone_composition, aes(x = density_zone, y = percentage, fill = celltype_clean)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = sprintf("Cell Type Composition by Density Zone - %s", sample_id),
      x = "Density Zone (0=Core/High, Higher=Outer/Low)",
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
  
  # 图2：Zone的spot数量（使用统一的zone颜色）
  zone_totals <- zone_composition %>%
    group_by(density_zone) %>%
    summarise(total = sum(count), .groups = "drop")
  
  p2 <- ggplot(zone_totals, aes(x = density_zone, y = total, fill = density_zone)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.5) +
    geom_text(aes(label = total), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Total Spots per Density Zone",
      x = "Density Zone (0=Core → Higher=Outer)",
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


# ===================================================================
# 辅助函数 4：绘制合并热图（统一配色版）
# ===================================================================

plot_combined_heatmap <- function(combined_data, CONFIG) {
  
  # 计算平均百分比
  heatmap_data <- combined_data %>%
    group_by(density_zone, celltype_clean) %>%
    summarise(
      mean_pct = mean(percentage),
      sd_pct = sd(percentage),
      n_samples = n(),
      .groups = "drop"
    )
  
  # 确保zone按顺序排列
  n_zones <- length(unique(heatmap_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- names(zone_colors)
  
  heatmap_data <- heatmap_data %>%
    mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 热图主体
  p <- ggplot(heatmap_data, aes(x = density_zone, y = celltype_clean, fill = mean_pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_pct)), size = 3.5, color = "black", fontface = "bold") +
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
    labs(
      title = "Cell Type Composition Across Density Zones (All Samples)",
      subtitle = sprintf("Averaged across %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (0=Core/High → Higher=Outer/Low)",
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
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1)
    )
  
  # 添加zone颜色条（顶部）
  zone_bar_data <- data.frame(
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 1
  )
  
  p_zone_bar <- ggplot(zone_bar_data, aes(x = density_zone, y = y_position, fill = density_zone)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # 合并图形
  p_final <- p_zone_bar / p + plot_layout(heights = c(0.05, 1))
  
  return(p_final)
}


# ===================================================================
# 辅助函数 5：绘制综合分析图（统一配色版）
# ===================================================================

plot_combined_analysis <- function(combined_data, CONFIG) {
  
  # 获取统一的颜色方案
  n_zones <- length(unique(combined_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- names(zone_colors)
  celltype_colors <- get_celltype_colors(unique(combined_data$celltype_clean))
  
  # 确保zone按顺序排列
  combined_data <- combined_data %>%
    mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    )
  
  # 1. 箱线图：每个区域的细胞类型比例分布
  p1 <- ggplot(combined_data, aes(x = density_zone, y = percentage, fill = density_zone)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    facet_wrap(~celltype_clean, scales = "free_y", ncol = 4) +
    labs(
      title = "Cell Type Percentage Distribution by Density Zone",
      subtitle = sprintf("Data from %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
  
  # 2. 趋势图：核心到外围的变化
  trend_data <- combined_data %>%
    group_by(celltype_clean, zone_numeric, density_zone) %>%
    summarise(
      mean_pct = mean(percentage),
      se_pct = sd(percentage) / sqrt(n()),
      .groups = "drop"
    )
  
  p2 <- ggplot(trend_data, aes(x = zone_numeric, y = mean_pct, color = celltype_clean, group = celltype_clean)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct), 
      width = 0.2, 
      linewidth = 0.8,
      alpha = 0.7
    ) +
    scale_color_manual(values = celltype_colors, name = "Cell Type") +
    scale_x_continuous(
      breaks = 0:(n_zones - 1),
      labels = zone_levels
    ) +
    labs(
      title = "Cell Type Enrichment Trend Across Density Zones",
      subtitle = "Mean ± SE across all samples",
      x = "Density Zone (0=Core/High → Higher=Outer/Low)",
      y = "Mean Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  # 3. 添加zone颜色参考条
  zone_ref_data <- data.frame(
    zone_numeric = 0:(n_zones - 1),
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 0
  )
  
  p2 <- p2 +
    geom_tile(
      data = zone_ref_data,
      aes(x = zone_numeric, y = y_position, fill = density_zone),
      height = max(trend_data$mean_pct) * 0.05,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = zone_colors, guide = "none")
  
  # 合并
  p_combined <- p1 / p2 + plot_layout(heights = c(2, 1.2))
  
  return(p_combined)
}


# ===================================================================
# 辅助函数 6：生成统计摘要
# ===================================================================

generate_summary_statistics <- function(combined_data) {
  
  # 计算每种细胞类型在不同区域的富集情况
  summary <- combined_data %>%
    mutate(zone_numeric = as.numeric(gsub("Zone_", "", density_zone))) %>%
    group_by(celltype_clean) %>%
    summarise(
      mean_pct_all = mean(percentage),
      sd_pct_all = sd(percentage),
      max_zone = density_zone[which.max(percentage)],
      max_pct = max(percentage),
      min_zone = density_zone[which.min(percentage)],
      min_pct = min(percentage),
      # Zone_0和Zone_1是核心区，其他是外围
      core_enrichment = mean(percentage[zone_numeric <= 1]) - mean(percentage[zone_numeric > 1]),
      n_samples = length(unique(sample)),
      .groups = "drop"
    ) %>%
    arrange(desc(core_enrichment))
  
  return(summary)
}


# ===================================================================
# 辅助函数：%||% 操作符
# ===================================================================
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}