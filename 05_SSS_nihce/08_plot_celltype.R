#!/usr/bin/env Rscript
# ===================================================================
# 细胞类型 Niche 分析模块（简化版 + 进度条）
# 功能：分析不同密度区域的细胞类型分布和富集
# ===================================================================

library(future)
library(future.apply)
library(progressr)
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
    celltype_col = "predicted.id",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE,
    seurat_basename = NULL
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   细胞类型 Niche 分析（多线程并行）\n")
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
  
  n_workers <- CONFIG$n_workers %||% 4
  
  cat(sprintf("📊 将分析 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 密度分区: %d 个区域 (Zone_0=核心, Zone_%d=外围)\n", 
              density_bins, density_bins - 1))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 3. 设置并行和进度条
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  # 设置进度条
  has_handlers <- !is.null(progressr::handlers(NULL))
  
  if (!has_handlers) {
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "[:bar] :percent | 已完成: :current/:total | 预计剩余: :eta | :message",
        width    = 80,
        complete = "=",
        clear    = FALSE
      )
    ))
  }
  
  start_time <- Sys.time()
  
  # ========================================
  # 4. 并行处理样本
  # ========================================
  
  cat("🔬 开始分析样本...\n\n")
  
  progressr::with_progress({
    
    p <- progressr::progressor(
      steps = length(sample_list),
      message = "分析细胞类型 Niche"
    )
    
    results <- future.apply::future_lapply(
      
      names(sample_list),
      
      function(sample_id) {
        
        process_single_sample(
          sample_id = sample_id,
          sample_list = sample_list,
          CONFIG = CONFIG,
          celltype_col = celltype_col,
          density_bins = density_bins,
          plot_overlay = plot_overlay,
          plot_composition = plot_composition,
          progressor = p
        )
      },
      
      future.seed = TRUE,
      future.chunk.size = 1
    )
  })
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  future::plan(future::sequential)
  
  # ========================================
  # 5. 统计样本处理结果
  # ========================================
  
  print_sample_summary(results, sample_list, elapsed)
  
  # ========================================
  # 6. 生成综合分析
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
  # 7. 最终总结
  # ========================================
  
  print_final_summary(results, sample_list, start_time, combined_data,
                     plot_overlay, plot_composition, plot_heatmap, plot_combined,
                     CONFIG)
  
  # ========================================
  # 8. 返回结果
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