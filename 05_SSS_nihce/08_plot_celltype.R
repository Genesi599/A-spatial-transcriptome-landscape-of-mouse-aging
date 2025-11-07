#!/usr/bin/env Rscript
# ===================================================================
# 细胞类型 Niche 分析模块（串联版 - 无并行依赖）
# 功能：分析不同密度区域的细胞类型分布和富集
# ===================================================================

library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

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
source(file.path(utils_dir, "08_validation.R"))
source(file.path(utils_dir, "09_save_plots.R"))
source(file.path(utils_dir, "10_summary.R"))

cat("✅ 已加载所有工具函数\n")


# ===================================================================
# 主函数：细胞类型 Niche 分析
# ===================================================================

#' 细胞类型 Niche 分析
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param density_bins 密度分区数量
#' @param celltype_col 细胞类型列名
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
#' @param seurat_basename 文件基础名
#' 
#' @return 处理结果列表
#'
analyze_celltype_niche <- function(
    sample_list,
    CONFIG,
    density_bins = 10,
    celltype_col = "celltype",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE,
    seurat_basename = NULL
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   细胞类型 Niche 分析\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  validate_inputs(sample_list, CONFIG)
  validate_required_functions()
  
  # ========================================
  # 2. 初始化配置
  # ========================================
  
  setup_colors(sample_list[[1]], CONFIG, celltype_col, density_bins)
  
  cat(sprintf("📊 将分析 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 密度分区: %d 个区域 (Zone_0=核心, Zone_%d=外围)\n", 
              density_bins, density_bins - 1))
  cat(sprintf("🔧 使用串联模式（稳定可靠）\n\n"))
  
  start_time <- Sys.time()
  
  # ========================================
  # 3. 串联处理样本
  # ========================================
  
  cat("🔬 开始分析样本...\n\n")
  
  results <- list()
  total_samples <- length(sample_list)
  
  for (i in seq_along(sample_list)) {
    
    sample_id <- names(sample_list)[i]
    
    cat(sprintf("[%2d/%2d] ", i, total_samples))
    
    tryCatch({
      
      # 调用单样本处理函数
      result <- process_single_sample(
        sample_id = sample_id,
        sample_list = sample_list,
        CONFIG = CONFIG,
        celltype_col = celltype_col,
        density_bins = density_bins,
        plot_overlay = plot_overlay,
        plot_composition = plot_composition,
        progressor = NULL  # 串联模式不需要进度对象
      )
      
      results[[sample_id]] <- result
      
      # 输出成功信息
      if (result$success) {
        cat(sprintf("✅ %s", sample_id))
        
        # 添加额外统计信息
        if (!is.null(result$n_spots)) {
          cat(sprintf(" (%d spots", result$n_spots))
          
          if (!is.null(result$n_high)) {
            cat(sprintf(", %d high", result$n_high))
          }
          
          if (!is.null(result$n_celltypes)) {
            cat(sprintf(", %d celltypes", result$n_celltypes))
          }
          
          cat(")")
        }
        
        cat("\n")
      } else {
        cat(sprintf("❌ %s - %s\n", sample_id, result$error %||% "Unknown error"))
      }
      
      # 清理内存
      if (i %% 3 == 0) gc(verbose = FALSE)
      
    }, error = function(e) {
      cat(sprintf("❌ %s - %s\n", sample_id, e$message))
      results[[sample_id]] <- list(
        sample = sample_id,
        success = FALSE,
        error = as.character(e$message)
      )
    })
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  cat(sprintf("\n⏱️  分析耗时: %.2f 分钟\n", elapsed / 60))
  
  # ========================================
  # 4. 统计样本处理结果
  # ========================================
  
  print_sample_summary(results, sample_list, elapsed)
  
  # ========================================
  # 5. 生成综合分析
  # ========================================
  
  combined_data <- collect_combined_data(results)
  
  if (nrow(combined_data) > 0) {
    generate_combined_analysis(
      combined_data = combined_data,
      CONFIG = CONFIG,
      seurat_basename = seurat_basename,
      plot_heatmap = plot_heatmap,
      plot_combined = plot_combined
    )
  }
  
  # ========================================
  # 6. 最终总结
  # ========================================
  
  print_final_summary(results, sample_list, start_time, combined_data,
                     plot_overlay, plot_composition, plot_heatmap, plot_combined,
                     CONFIG)
  
  # ========================================
  # 7. 返回结果
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    elapsed_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    combined_data = combined_data,
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 08_plot_celltype.R 已加载\n")