# ===================================================================
# 08_plot_celltype.R
# 细胞类型 + 等高线分析完整工作流（模块化版本）
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

# ===================================================================
# 加载工具函数
# ===================================================================

utils_dir <- "08_plot_celltype_utils"

source(file.path(utils_dir, "00_operators.R"))
source(file.path(utils_dir, "01_color_schemes.R"))
source(file.path(utils_dir, "02_density_zones.R"))
source(file.path(utils_dir, "03_plot_overlay.R"))
source(file.path(utils_dir, "04_plot_composition.R"))
source(file.path(utils_dir, "05_plot_heatmap.R"))
source(file.path(utils_dir, "06_plot_combined.R"))
source(file.path(utils_dir, "07_statistics.R"))

cat("✅ 已加载所有工具函数\n")


# ===================================================================
# 主函数：细胞类型等高线分析
# ===================================================================

#' 细胞类型 + Clock Gene Niche 等高线综合分析
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 要分析的样本列表
#' @param CONFIG 配置列表
#' @param density_bins 等高线分级数量，默认 10（对应10个区域）
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
    density_bins = 10,
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
        tibble::rownames_to_column("barcode") %>%
        dplyr::left_join(coords %>% tibble::rownames_to_column("barcode"), by = "barcode") %>%
        dplyr::filter(!is.na(col), !is.na(row))
      
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
        dplyr::left_join(
          density_data$spot_zones %>% dplyr::select(col, row, density_zone, density_value),
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
        dplyr::filter(!is.na(density_zone)) %>%
        dplyr::group_by(density_zone, celltype_clean) %>%
        dplyr::summarise(count = n(), .groups = "drop") %>%
        dplyr::group_by(density_zone) %>%
        dplyr::mutate(
          total = sum(count),
          percentage = 100 * count / total
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sample = sample_id)
      
      cat(sprintf("   ✅ 密度分区完成，共 %d 个区域\n", 
                  length(unique(zone_composition$density_zone))))
      
      # 打印每个zone的统计
      zone_stats <- df %>%
        dplyr::filter(!is.na(density_zone)) %>%
        dplyr::group_by(density_zone) %>%
        dplyr::summarise(
          n_spots = n(),
          mean_density = mean(density_value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(density_zone)
      
      cat("   Zone统计:\n")
      for (j in 1:nrow(zone_stats)) {
        cat(sprintf("     %s: %d spots (mean density: %.3f)\n",
                    zone_stats$density_zone[j],
                    zone_stats$n_spots[j],
                    zone_stats$mean_density[j]))
      }
      
      # 保存到总体数据
      all_sample_stats[[sample_id]] <- zone_composition
      combined_data <- dplyr::bind_rows(combined_data, zone_composition)
      
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