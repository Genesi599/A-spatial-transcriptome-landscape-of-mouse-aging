# main.R

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
# 🔧 强制使用 dplyr 函数（全局修复）
# ===================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# 全局覆盖（确保在任何地方都使用 dplyr 版本）
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
group_by <- dplyr::group_by
summarize <- dplyr::summarize
summarise <- dplyr::summarise
left_join <- dplyr::left_join
right_join <- dplyr::right_join
inner_join <- dplyr::inner_join
full_join <- dplyr::full_join
rownames_to_column <- tibble::rownames_to_column
column_to_rownames <- tibble::column_to_rownames

cat("✅ 已全局设置 dplyr/tibble 函数\n\n")

# ===================================================================
# 错误追踪设置
# ===================================================================
options(error = function() {
  traceback(max.lines = 30)
  cat("\n📝 错误位置已显示\n")
})

# ===================================================================
# 加载配置和模块
# ===================================================================
if (!exists("CONFIG")) {
  source("00_config.R")
}
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
source("14_gene_list_utils.R")       # 基因列表工具



# ===================================================================
# 单文件处理函数（核心流程 - 已移除 tryCatch）
# ===================================================================

#' 处理单个 Seurat 文件（无错误捕获，直接暴露错误）
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
  
  # ========================================
  # 步骤 1-5: 数据准备和分析
  # ========================================
  
  cat("\n【步骤 1/9】环境设置\n")
  setup_environment(config)
  
  # ----------------------------------------
  cat("\n【步骤 2/9】加载 Seurat 对象\n")
  seurat_obj <- load_seurat_object(config)
  
  # 🔍 类型检查
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf(
      "步骤2失败: load_seurat_object 返回了 %s 类型（应该是 Seurat）",
      class(seurat_obj)[1]
    ))
  }
  cat(sprintf("   ✅ 对象类型: %s | 细胞数: %d\n", 
              class(seurat_obj)[1], ncol(seurat_obj)))
  
  genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
  
  # ----------------------------------------
  cat("\n【步骤 3/9】计算 Clock Gene Score\n")
  seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, config)
  
  # 🔍 类型检查
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf(
      "步骤3失败: calculate_module_score 返回了 %s 类型（应该是 Seurat）",
      class(seurat_obj)[1]
    ))
  }
  cat(sprintf("   ✅ 对象类型: %s\n", class(seurat_obj)[1]))
  
  # ----------------------------------------
  cat("\n【步骤 4/9】识别高表达区域\n")
  result <- define_high_expression(seurat_obj, config)
  
  # 🔍 类型检查
  if (!is.list(result)) {
    stop(sprintf(
      "步骤4失败: define_high_expression 返回了 %s 类型（应该是 list）",
      class(result)[1]
    ))
  }
  
  if (!"seurat_obj" %in% names(result)) {
    stop("步骤4失败: 返回结果中缺少 seurat_obj 元素")
  }
  
  seurat_obj <- result$seurat_obj
  threshold <- result$threshold
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf(
      "步骤4失败: result$seurat_obj 是 %s 类型（应该是 Seurat）",
      class(seurat_obj)[1]
    ))
  }
  cat(sprintf("   ✅ 对象类型: %s | 阈值: %.3f\n", 
              class(seurat_obj)[1], threshold))
  
  # ----------------------------------------
  cat("\n【步骤 5/9】Niche 分析\n")
  cat(sprintf("   🔍 输入对象类型: %s\n", class(seurat_obj)[1]))
  
  seurat_obj <- perform_niche_analysis(seurat_obj, threshold, config)
  
  # 🔍 类型检查
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf(
      "步骤5失败: perform_niche_analysis 返回了 %s 类型（应该是 Seurat）",
      class(seurat_obj)[1]
    ))
  }
  cat(sprintf("   ✅ 返回对象类型: %s\n", class(seurat_obj)[1]))
  
  # ========================================
  # 步骤 6: 样本预处理（统一切分）
  # ========================================
  
  cat("\n【步骤 6/9】样本预处理\n")
  
  samples <- unique(seurat_obj$orig.ident)
  samples_to_plot <- if (config$debug_mode) {
    head(samples, config$debug_sample_limit %||% 3)
  } else {
    samples
  }
  
  cat(sprintf("   🔬 将处理 %d 个样本\n", length(samples_to_plot)))
  
  # 一次性切分所有样本
  sample_list <- preprocess_samples(seurat_obj, samples_to_plot, config)
  
  # 🔍 验证 sample_list
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop(sprintf(
      "步骤6失败: preprocess_samples 返回了无效的 sample_list (类型: %s, 长度: %d)",
      class(sample_list)[1], length(sample_list)
    ))
  }
  cat(sprintf("   ✅ 样本列表: %d 个样本\n", length(sample_list)))
  
  # 更新配置中的线程数（基于内存估算）
  recommended_workers <- attr(sample_list, "recommended_workers")
  if (!is.null(recommended_workers)) {
    config$n_workers <- recommended_workers
  }
  
  # ========================================
  # 步骤 7-9: 可视化分析
  # ========================================
  
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
  
  # ========================================
  # 保存结果
  # ========================================
  
  save_results(seurat_obj, config)
  
  # ========================================
  # 完成
  # ========================================
  
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
}


# ===================================================================
# 批量处理主函数（简化版）
# ===================================================================

#' 批量处理主函数
#'
#' @return 批量处理结果
#'
main_batch <- function() {

  ## —— 瞬间快照：谁踩了 CONFIG$gene_list_path ——
  cat("=== entering main_batch ===\n")
  cat(sprintf("CONFIG$gene_list_path: '%s'\n", CONFIG$gene_list_path))
  cat(sprintf("class: %s ; length: %d\n", class(CONFIG$gene_list_path), length(CONFIG$gene_list_path)))

  
  batch_start_time <- Sys.time()
  
  print_batch_header()
  
  cat("\n【初始化】环境设置\n")
  init_result <- initialize_environment(
    config = CONFIG,
    custom_scripts = c("niche_marker.R", "SSS_isoheight_plot.R")
  )
  CONFIG <- init_result$config

  cat(sprintf("post-init CONFIG$gene_list_path: '%s' cls=%s len=%d\n",
            CONFIG$gene_list_path %||% "<NULL>",
            class(CONFIG$gene_list_path)[1],
            length(CONFIG$gene_list_path)))
  
  # 验证基因列表文件（使用当前 CONFIG，而非默认值）
  if (!file.exists(CONFIG$gene_list_path)) {
    stop(sprintf("❌ 基因列表文件不存在: %s", CONFIG$gene_list_path))
  }
  
  validate_output_directory(CONFIG)
  
  seurat_files <- scan_seurat_files(CONFIG)
  if (length(seurat_files) == 0) {
    stop("❌ 未找到可处理的 Seurat 文件")
  }
  print_file_list(seurat_files)
  
  cat("\n【加载】基因列表\n")
  gene_list <- load_gene_list_once(CONFIG)
  gene_name <- tools::file_path_sans_ext(
    basename(CONFIG$gene_list_path)
  )
  
  cat(sprintf("   📋 %s (%d genes)\n", 
              gene_name, length(gene_list)))
  
  total_tasks <- length(seurat_files)
  cat(sprintf(
    "\n⚠️  将处理 %d 个 Seurat 文件\n",
    total_tasks
  ))
  
  if (!confirm_batch_processing(seurat_files, CONFIG)) {
    cat("❌ 已取消处理\n")
    return(invisible(NULL))
  }
  
  results <- lapply(seurat_files, function(sf) {
    process_seurat_file(sf, gene_list, CONFIG)
  })
  
  batch_end_time <- Sys.time()
  total_elapsed <- difftime(
    batch_end_time, batch_start_time, units = "mins"
  )
  
  print_batch_summary(results, total_elapsed, CONFIG)
  
  log_files <- save_batch_logs(
    results, batch_start_time, batch_end_time, CONFIG
  )

  aggregate_score_statistics(CONFIG$output_base_dir)

  cat("\n🎉 批量处理完成！\n\n")
  
  return(invisible(list(
    results = results,
    gene_list = gene_name,
    summary = create_summary_object(
      results, total_elapsed, log_files
    )
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