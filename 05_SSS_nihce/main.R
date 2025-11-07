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
# 加载配置和模块
# ===================================================================
options(error = function() {
  traceback(max.lines = 20)
  cat("\n📝 错误位置已显示\n")
})

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
#' 处理单个 Seurat 文件（带详细错误追踪）
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
    
  }, error = function(e) {
    
    # ========================================
    # 🔴 详细错误追踪
    # ========================================
    
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════╗\n")
    cat("║            🔴 详细错误信息                                 ║\n")
    cat("╚═══════════════════════════════════════════════════════════╝\n\n")
    
    cat(sprintf("📂 文件: %s\n", seurat_basename))
    cat(sprintf("❌ 错误: %s\n\n", e$message))
    
    # 打印错误调用位置
    if (!is.null(e$call)) {
      cat("📍 错误调用:\n")
      call_str <- deparse(e$call)
      for (line in call_str) {
        cat(sprintf("   %s\n", line))
      }
      cat("\n")
    }
    
    # 🔑 打印完整调用堆栈
    cat("📚 调用堆栈（最近20层）:\n")
    cat(paste(rep("─", 70), collapse = ""), "\n")
    
    all_calls <- sys.calls()
    n_calls <- length(all_calls)
    
    # 从后往前打印最近的调用
    start_idx <- max(1, n_calls - 19)
    for (i in start_idx:n_calls) {
      call_text <- deparse(all_calls[[i]])[1]
      
      # 截断过长的调用
      if (nchar(call_text) > 100) {
        call_text <- paste0(substr(call_text, 1, 97), "...")
      }
      
      # 高亮包含 filter 的调用
      marker <- if (grepl("filter", call_text, ignore.case = TRUE)) {
        " ⚠️ "
      } else {
        "    "
      }
      
      cat(sprintf("%s%2d: %s\n", marker, i, call_text))
    }
    
    cat(paste(rep("─", 70), collapse = ""), "\n\n")
    
    # 保存详细日志到文件
    log_file <- sprintf(
      "error_detail_%s_%s.txt",
      seurat_basename,
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
    
    tryCatch({
      sink(log_file)
      cat("╔═══════════════════════════════════════════════════════════╗\n")
      cat("║            错误详细信息                                    ║\n")
      cat("╚═══════════════════════════════════════════════════════════╝\n\n")
      
      cat("文件:", seurat_basename, "\n")
      cat("时间:", format(Sys.time()), "\n")
      cat("错误消息:", e$message, "\n\n")
      
      cat("错误调用:\n")
      print(e$call)
      cat("\n")
      
      cat("完整堆栈:\n")
      cat(paste(rep("─", 70), collapse = ""), "\n")
      for (i in seq_along(all_calls)) {
        cat(sprintf("%d: %s\n", i, deparse(all_calls[[i]])[1]))
      }
      cat(paste(rep("─", 70), collapse = ""), "\n\n")
      
      cat("会话信息:\n")
      print(sessionInfo())
      
      sink()
      
      cat(sprintf("📝 详细日志已保存: %s\n\n", log_file))
      
    }, error = function(log_err) {
      sink()  # 确保 sink 被关闭
      cat(sprintf("⚠️  无法保存日志: %s\n\n", log_err$message))
    })
    
    # ========================================
    # 原有的错误处理
    # ========================================
    
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
      error = e$message,
      error_call = if (!is.null(e$call)) deparse(e$call) else NULL,
      log_file = if (exists("log_file")) log_file else NULL
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

  # ✅ 新增：接收更新后的 CONFIG
  CONFIG <- init_result$config

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