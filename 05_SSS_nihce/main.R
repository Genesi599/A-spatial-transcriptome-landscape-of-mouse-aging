#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis - Main Script (Optimized)
# Author: Zhangbin
# Date: 2024-11-04
# Optimized: 2024-11-07
#   - 模块化拆分
#   - 统一环境初始化
#   - 内存管理优化
# ===================================================================

# ===================================================================
# 加载配置和模块
# ===================================================================

source("00_config.R")
source("01_setup.R")
source("02_utils.R")
source("03_load_data.R")
source("04_module_score.R")
source("05_niche_analysis.R")
source("06_plot_isoheight.R")
source("07_plot_spatial.R")
source("08_plot_celltype.R")
source("09_save_results.R")
source("10_batch_processing.R")       # 批量处理
source("11_sample_preprocessing.R")   # 样本预处理
source("12_file_utils.R")             # 文件工具
source("13_reporting.R")              # 报告生成


# ===================================================================
# 单文件处理函数（核心流程）
# ===================================================================

#' 处理单个 Seurat 文件
#'
#' @param seurat_path Seurat 文件路径
#' @param gene_list 基因列表
#' @param base_config 基础配置
#' 
#' @return 处理结果列表
#'
process_seurat_file <- function(seurat_path, gene_list, base_config) {
  
  # 1. 更新配置
  config <- update_config_for_file(seurat_path, base_config)
  seurat_basename <- tools::file_path_sans_ext(basename(seurat_path))
  
  # 2. 打印处理信息
  print_file_header(seurat_basename)
  
  file_start_time <- Sys.time()
  
  tryCatch({
    
    # ----------------------------------------
    # 步骤 1-5: 数据准备和分析
    # ----------------------------------------
    cat("\n【步骤 1/9】环境设置\n")
    setup_environment(config)
    
    cat("\n【步骤 2/9】加载 Seurat 对象\n")
    seurat_obj <- load_seurat_object(config)
    genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
    
    cat("\n【步骤 3/9】计算 Clock Gene Score\n")
    seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, config)
    
    cat("\n【步骤 4/9】识别高表达区域\n")
    result <- define_high_expression(seurat_obj, config)
    seurat_obj <- result$seurat_obj
    threshold <- result$threshold
    
    cat("\n【步骤 5/9】Niche 分析\n")
    seurat_obj <- perform_niche_analysis(seurat_obj, threshold, config)
    
    # ----------------------------------------
    # 步骤 6: 样本预处理（统一切分）
    # ----------------------------------------
    cat("\n【步骤 6/9】样本预处理\n")
    
    samples <- unique(seurat_obj$orig.ident)
    samples_to_plot <- if (config$debug_mode) {
      head(samples, config$debug_sample_limit %||% 3)
    } else {
      samples
    }
    
    # 一次性切分所有样本
    sample_list <- preprocess_samples(seurat_obj, samples_to_plot, config)
    
    # 更新配置中的线程数（基于内存估算）
    recommended_workers <- attr(sample_list, "recommended_workers")
    if (!is.null(recommended_workers)) {
      config$n_workers <- recommended_workers
    }
    
    # ----------------------------------------
    # 步骤 7-9: 可视化分析
    # ----------------------------------------
    cat("\n【步骤 7/9】绘制等高线密度图\n")
    iso_results <- plot_isoheight(
      sample_list = sample_list,
      CONFIG = config
    )
    
    cat("\n【步骤 8/9】绘制空间梯度图\n")
    spatial_results <- plot_spatial_gradient(
      sample_list = sample_list,
      CONFIG = config
    )
    
    cat("\n【步骤 9/9】细胞类型 Niche 分析\n")
    celltype_results <- analyze_celltype_niche(
      sample_list = sample_list,
      CONFIG = config,
      seurat_basename = seurat_basename
    )
    
    # ----------------------------------------
    # 保存结果
    # ----------------------------------------
    save_results(seurat_obj, config)
    
    # ----------------------------------------
    # 完成
    # ----------------------------------------
    file_end_time <- Sys.time()
    file_elapsed <- difftime(file_end_time, file_start_time, units = "mins")
    
    print_file_success(seurat_basename, length(sample_list), file_elapsed, config)
    
    # 清理内存
    rm(seurat_obj, sample_list)
    gc(verbose = FALSE)
    
    return(list(
      success = TRUE,
      file = seurat_basename,
      processing_time = as.numeric(file_elapsed),
      n_samples = length(samples_to_plot),
      error = NULL
    ))
    
  }, error = function(e) {
    
    file_end_time <- Sys.time()
    file_elapsed <- difftime(file_end_time, file_start_time, units = "mins")
    
    print_file_failure(seurat_basename, e$message, file_elapsed)
    
    # 清理内存
    gc(verbose = FALSE)
    
    return(list(
      success = FALSE,
      file = seurat_basename,
      processing_time = as.numeric(file_elapsed),
      n_samples = 0,
      error = e$message
    ))
  })
}


# ===================================================================
# 批量处理主函数（简化版）
# ===================================================================

#' 批量处理主函数
#'
#' @return 批量处理结果
#'
main_batch <- function() {
  
  batch_start_time <- Sys.time()
  
  print_batch_header()
  
  # ----------------------------------------
  # 1. 统一初始化环境
  # ----------------------------------------
  cat("\n【初始化】环境设置\n")
  
  init_result <- initialize_environment(
    config = CONFIG,
    custom_scripts = c("niche_marker.R", "SSS_isoheight_plot.R")
  )
  
  if (length(init_result$packages$failed) > 0) {
    warning("⚠️  部分包加载失败，可能影响分析")
  }
  
  # ----------------------------------------
  # 2. 验证输出目录
  # ----------------------------------------
  validate_output_directory(CONFIG)
  
  # ----------------------------------------
  # 3. 扫描输入文件
  # ----------------------------------------
  seurat_files <- scan_seurat_files(CONFIG)
  
  if (length(seurat_files) == 0) {
    stop("❌ 未找到可处理的文件")
  }
  
  print_file_list(seurat_files)
  
  # 确认处理
  if (!confirm_batch_processing(seurat_files, CONFIG)) {
    cat("❌ 已取消处理\n")
    return(invisible(NULL))
  }
  
  # ----------------------------------------
  # 4. 加载基因列表（只加载一次）
  # ----------------------------------------
  gene_list <- load_gene_list_once(CONFIG)
  
  # ----------------------------------------
  # 5. 批量处理文件
  # ----------------------------------------
  results <- process_all_files(seurat_files, gene_list, CONFIG)
  
  # ----------------------------------------
  # 6. 生成总结报告
  # ----------------------------------------
  batch_end_time <- Sys.time()
  total_elapsed <- difftime(batch_end_time, batch_start_time, units = "mins")
  
  print_batch_summary(results, total_elapsed, CONFIG)
  
  log_files <- save_batch_logs(results, batch_start_time, batch_end_time, CONFIG)
  
  cat("\n🎉 批量处理完成！\n\n")
  
  return(invisible(list(
    results = results,
    summary = create_summary_object(results, total_elapsed, log_files)
  )))
}


# ===================================================================
# 辅助操作符
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}


# ===================================================================
# 运行主流程
# ===================================================================

if (!interactive()) {
  main_batch()
}

cat("✅ main.R 已加载\n")
cat("📚 使用 main_batch() 开始批量处理\n\n")