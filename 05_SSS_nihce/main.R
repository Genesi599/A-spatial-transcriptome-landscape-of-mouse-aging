#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis - Main Script (Optimized)
# Author: Zhangbin
# Date: 2024-11-04
# Optimized: 2024-11-07
#   - 统一环境初始化
#   - 样本预处理优化
#   - 内存管理改进
#   - 批量处理增强
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


# ===================================================================
# 样本预处理模块（统一切分）
# ===================================================================

#' 预处理样本（一次性切分所有样本）
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 要处理的样本列表
#' @param config 配置对象
#' @return 切分后的样本列表
#'
preprocess_samples <- function(seurat_obj, samples_to_plot, config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本预处理（统一切分）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 验证样本
  available_samples <- unique(seurat_obj$orig.ident)
  invalid_samples <- setdiff(samples_to_plot, available_samples)
  
  if (length(invalid_samples) > 0) {
    warning(sprintf("⚠️ 以下样本不存在，将跳过: %s", 
                    paste(invalid_samples, collapse = ", ")))
    samples_to_plot <- intersect(samples_to_plot, available_samples)
  }
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本")
  }
  
  cat(sprintf("📊 原始数据: %d spots, %d 个样本\n", 
              ncol(seurat_obj), length(available_samples)))
  cat(sprintf("📊 将处理: %d 个样本\n\n", length(samples_to_plot)))
  
  # 开始切分
  cat("🔧 切分样本...\n")
  start_time <- Sys.time()
  
  sample_list <- list()
  sample_stats <- data.frame(
    sample_id = character(),
    n_spots = integer(),
    n_high = integer(),
    high_pct = numeric(),
    size_mb = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    
    seurat_subset <- tryCatch({
      subset(seurat_obj, subset = orig.ident == sample_id)
    }, error = function(e) {
      seurat_obj[, seurat_obj$orig.ident == sample_id]
    })
    
    if (ncol(seurat_subset) > 0) {
      sample_list[[sample_id]] <- seurat_subset
      
      # 统计信息
      n_spots <- ncol(seurat_subset)
      n_high <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
      high_pct <- 100 * mean(seurat_subset$ClockGene_High, na.rm = TRUE)
      size_mb <- as.numeric(object.size(seurat_subset)) / 1024^2
      
      sample_stats <- rbind(sample_stats, data.frame(
        sample_id = sample_id,
        n_spots = n_spots,
        n_high = n_high,
        high_pct = high_pct,
        size_mb = size_mb
      ))
      
      cat(sprintf("  [%2d/%2d] ✅ %-30s | %6d spots | %4d high (%.2f%%) | %.2f MB\n",
                  i, length(samples_to_plot),
                  sample_id,
                  n_spots, n_high, high_pct, size_mb))
    } else {
      warning(sprintf("⚠️ 样本 %s 无数据，已跳过", sample_id))
    }
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 汇总统计
  total_spots <- sum(sample_stats$n_spots)
  total_size_mb <- sum(sample_stats$size_mb)
  avg_size_mb <- mean(sample_stats$size_mb)
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("✅ 成功切分 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 总计: %d spots (%.2f MB)\n", total_spots, total_size_mb))
  cat(sprintf("📊 平均: %.0f spots/样本 (%.2f MB/样本)\n", 
              total_spots / length(sample_list), avg_size_mb))
  cat(sprintf("⏱️  耗时: %.2f 秒\n", as.numeric(elapsed)))
  
  # 动态调整线程数
  max_memory_gb <- config$max_memory_gb %||% 100
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  recommended_workers <- min(safe_workers, length(sample_list), config$n_workers)
  
  cat(sprintf("\n💡 推荐线程数: %d (基于内存 %.0f GB)\n", 
              recommended_workers, max_memory_gb))
  
  if (recommended_workers < config$n_workers) {
    cat(sprintf("⚠️  原配置 %d 线程可能导致内存不足，已自动调整\n", 
                config$n_workers))
    config$n_workers <- recommended_workers
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 保存统计到配置
  attr(sample_list, "stats") <- sample_stats
  attr(sample_list, "recommended_workers") <- recommended_workers
  
  return(sample_list)
}


# ===================================================================
# 单文件处理函数（优化版）
# ===================================================================

process_seurat_file <- function(seurat_path, gene_list, base_config) {
  
  # 1. 更新配置
  config <- base_config
  config$seurat_path <- seurat_path
  
  # 提取文件名
  seurat_basename <- tools::file_path_sans_ext(basename(seurat_path))
  config$output_dir <- file.path(config$output_base_dir, seurat_basename)
  
  # 更新所有目录路径
  config <- update_config_paths(config)
  
  # 2. 打印处理信息
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  处理文件: %-46s ║\n", seurat_basename))
  cat("╚═══════════════════════════════════════════════════════════╝\n")
  
  file_start_time <- Sys.time()
  
  tryCatch({
    
    # ----------------------------------------
    # 步骤 1: 环境设置
    # ----------------------------------------
    cat("\n【步骤 1/9】环境设置\n")
    setup_environment(config)
    
    # ----------------------------------------
    # 步骤 2: 加载数据
    # ----------------------------------------
    cat("\n【步骤 2/9】加载 Seurat 对象\n")
    seurat_obj <- load_seurat_object(config)
    genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
    
    # ----------------------------------------
    # 步骤 3: 计算评分
    # ----------------------------------------
    cat("\n【步骤 3/9】计算 Clock Gene Score\n")
    seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, config)
    
    # ----------------------------------------
    # 步骤 4: 识别高表达区域
    # ----------------------------------------
    cat("\n【步骤 4/9】识别高表达区域\n")
    result <- define_high_expression(seurat_obj, config)
    seurat_obj <- result$seurat_obj
    threshold <- result$threshold
    
    # ----------------------------------------
    # 步骤 5: Niche 分析
    # ----------------------------------------
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
    
    # 【关键优化】一次性切分所有样本
    sample_list <- preprocess_samples(seurat_obj, samples_to_plot, config)
    
    # 更新配置中的线程数
    recommended_workers <- attr(sample_list, "recommended_workers")
    if (!is.null(recommended_workers)) {
      config$n_workers <- recommended_workers
    }
    
    # ----------------------------------------
    # 步骤 7: 绘制等高线图
    # ----------------------------------------
    cat("\n【步骤 7/9】绘制等高线密度图\n")
    iso_results <- plot_isoheight(
      sample_list = sample_list,
      CONFIG = config
    )
    
    # ----------------------------------------
    # 步骤 8: 绘制空间梯度图
    # ----------------------------------------
    cat("\n【步骤 8/9】绘制空间梯度图\n")
    spatial_results <- plot_spatial_gradient(
      sample_list = sample_list,
      CONFIG = config
    )
    
    # ----------------------------------------
    # 步骤 9: 细胞类型分析
    # ----------------------------------------
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
    
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════╗\n")
    cat("║                    处理完成                                ║\n")
    cat("╚═══════════════════════════════════════════════════════════╝\n\n")
    
    cat(sprintf("✅ 文件: %s\n", seurat_basename))
    cat(sprintf("📊 处理样本: %d\n", length(sample_list)))
    cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(file_elapsed)))
    cat(sprintf("📁 输出: %s\n", config$output_dir))
    
    print_summary(config)
    
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
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
    
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════╗\n")
    cat("║                    处理失败                                ║\n")
    cat("╚═══════════════════════════════════════════════════════════╝\n\n")
    
    cat(sprintf("❌ 文件: %s\n", seurat_basename))
    cat(sprintf("❌ 错误: %s\n", e$message))
    cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(file_elapsed)))
    
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
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
# 辅助函数：更新配置路径
# ===================================================================

update_config_paths <- function(config) {
  # 更新基础目录
  config$cache_dir <- file.path(config$output_dir, "cache")
  config$figure_dir <- file.path(config$output_dir, "figure")
  config$metadata_dir <- file.path(config$output_dir, "metadata")
  
  # 更新详细目录
  config$dirs <- list(
    cache = config$cache_dir,
    figure = config$figure_dir,
    metadata = config$metadata_dir,
    isoheight = file.path(config$figure_dir, "isoheight"),
    spatial = file.path(config$figure_dir, "spatial"),
    overlay = file.path(config$figure_dir, "isoheight", "01_overlay_plots"),
    celltype = file.path(config$figure_dir, "isoheight", "02_celltype_only"),
    composition = file.path(config$figure_dir, "isoheight", "03_composition_stats"),
    heatmaps = file.path(config$figure_dir, "isoheight", "04_heatmaps"),
    combined = file.path(config$figure_dir, "isoheight", "05_combined_analysis")
  )
  
  return(config)
}


# ===================================================================
# 批量处理主函数（优化版）
# ===================================================================

main_batch <- function() {
  
  batch_start_time <- Sys.time()
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║        Clock Gene Niche Analysis - Batch Processing       ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n")
  
  # ----------------------------------------
  # 0. 统一初始化环境
  # ----------------------------------------
  cat("\n【初始化】环境设置\n")
  
  init_result <- initialize_environment(
    config = CONFIG,
    custom_scripts = c("niche_marker.R", "SSS_isoheight_plot.R")
  )
  
  # 检查初始化结果
  if (length(init_result$packages$failed) > 0) {
    warning("部分包加载失败，可能影响分析")
  }
  
  # ----------------------------------------
  # 1. 验证输出目录
  # ----------------------------------------
  if (is.null(CONFIG$output_base_dir) || CONFIG$output_base_dir == "") {
    stop("❌ 未配置 output_base_dir")
  }
  
  if (!dir.exists(CONFIG$output_base_dir)) {
    cat(sprintf("📁 创建输出基础目录: %s\n", CONFIG$output_base_dir))
    dir.create(CONFIG$output_base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ----------------------------------------
  # 2. 扫描输入文件
  # ----------------------------------------
  seurat_files <- scan_seurat_files(CONFIG)
  
  if (length(seurat_files) == 0) {
    stop("❌ 未找到可处理的文件")
  }
  
  # 打印文件列表
  print_file_list(seurat_files)
  
  # 确认处理
  if (CONFIG$batch_mode && length(seurat_files) > 1 && interactive()) {
    response <- readline(prompt = sprintf(
      "即将处理 %d 个文件，是否继续? (y/n): ", length(seurat_files)))
    if (tolower(response) != "y") {
      cat("❌ 已取消处理\n")
      return(invisible(NULL))
    }
  }
  
  # ----------------------------------------
  # 3. 加载基因列表（只加载一次）
  # ----------------------------------------
  cat("\n【准备】加载基因列表\n")
  gene_list <- load_gene_list(CONFIG)
  cat(sprintf("✅ 加载了 %d 个基因\n\n", length(gene_list)))
  
  # ----------------------------------------
  # 4. 批量处理
  # ----------------------------------------
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   开始批量处理\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  results <- list()
  
  for (i in seq_along(seurat_files)) {
    
    cat(sprintf("\n╔═══════════════════════════════════════════════════════════╗\n"))
    cat(sprintf("║  进度: [%d/%d] (%.1f%%)                                     \n", 
                i, length(seurat_files), (i/length(seurat_files))*100))
    cat(sprintf("╚═══════════════════════════════════════════════════════════╝\n"))
    
    # 处理单个文件
    result <- process_seurat_file(
      seurat_path = seurat_files[i],
      gene_list = gene_list,
      base_config = CONFIG
    )
    
    results[[i]] <- result
    
    # 估计剩余时间
    if (i < length(seurat_files)) {
      avg_time <- mean(sapply(results, function(x) x$processing_time), na.rm = TRUE)
      remaining_time <- avg_time * (length(seurat_files) - i)
      cat(sprintf("\n📊 预计剩余时间: %.2f 分钟 (%.2f 小时)\n", 
                  remaining_time, remaining_time/60))
    }
    
    # 强制内存清理
    gc(verbose = FALSE)
  }
  
  # ----------------------------------------
  # 5. 生成总结报告
  # ----------------------------------------
  batch_end_time <- Sys.time()
  total_elapsed <- difftime(batch_end_time, batch_start_time, units = "mins")
  
  # 打印总结
  print_batch_summary(results, total_elapsed, CONFIG)
  
  # 保存日志
  log_files <- save_batch_logs(results, batch_start_time, batch_end_time, CONFIG)
  
  cat("\n🎉 批量处理完成！\n\n")
  
  return(invisible(list(
    results = results,
    summary = list(
      total = length(results),
      success = sum(sapply(results, function(x) x$success)),
      failed = sum(sapply(results, function(x) !x$success)),
      total_time = as.numeric(total_elapsed),
      log_file = log_files$log,
      csv_file = log_files$csv
    )
  )))
}


# ===================================================================
# 辅助函数：扫描 Seurat 文件
# ===================================================================

scan_seurat_files <- function(config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   扫描输入文件\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  if (config$batch_mode) {
    # 批量模式
    cat(sprintf("📁 扫描目录: %s\n", config$seurat_dir))
    cat(sprintf("🔍 文件模式: %s\n", config$seurat_pattern))
    cat(sprintf("🔍 递归搜索: %s\n\n", config$recursive_search))
    
    if (!dir.exists(config$seurat_dir)) {
      stop(sprintf("❌ 目录不存在: %s", config$seurat_dir))
    }
    
    seurat_files <- list.files(
      path = config$seurat_dir,
      pattern = config$seurat_pattern,
      full.names = TRUE,
      recursive = config$recursive_search
    )
    
    if (length(seurat_files) == 0) {
      stop(sprintf("❌ 未找到匹配文件 (模式: %s)", config$seurat_pattern))
    }
    
    cat(sprintf("✅ 找到 %d 个文件\n", length(seurat_files)))
    
    # 过滤文件
    if (!is.null(config$specific_files) || !is.null(config$exclude_files)) {
      original_count <- length(seurat_files)
      seurat_files <- filter_seurat_files(seurat_files, config)
      cat(sprintf("📋 过滤后剩余 %d 个文件 (原始: %d)\n", 
                  length(seurat_files), original_count))
    }
    
  } else {
    # 单文件模式
    if (!file.exists(config$seurat_path)) {
      stop(sprintf("❌ 文件不存在: %s", config$seurat_path))
    }
    seurat_files <- config$seurat_path
    cat(sprintf("📄 单文件模式: %s\n", basename(seurat_files)))
  }
  
  cat("\n")
  
  return(seurat_files)
}


# ===================================================================
# 辅助函数：打印文件列表
# ===================================================================

print_file_list <- function(seurat_files) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   待处理文件列表\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("%-4s %-50s %10s\n", "No.", "文件名", "大小"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (i in seq_along(seurat_files)) {
    file_size_gb <- file.size(seurat_files[i]) / (1024^3)
    cat(sprintf("%3d. %-50s %8.2f GB\n", 
                i, 
                basename(seurat_files[i]), 
                file_size_gb))
  }
  
  total_size_gb <- sum(file.size(seurat_files)) / (1024^3)
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-55s %8.2f GB\n", "总计:", total_size_gb))
  
  cat("\n")
}


# ===================================================================
# 辅助函数：打印批量处理总结
# ===================================================================

print_batch_summary <- function(results, total_elapsed, config) {
  
  success_count <- sum(sapply(results, function(x) x$success))
  fail_count <- length(results) - success_count
  
  cat("\n\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    批量处理总结                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("📊 总文件数: %d\n", length(results)))
  cat(sprintf("✅ 成功: %d (%.1f%%)\n", 
              success_count, (success_count/length(results))*100))
  cat(sprintf("❌ 失败: %d (%.1f%%)\n", 
              fail_count, (fail_count/length(results))*100))
  cat(sprintf("⏱️  总耗时: %.2f 分钟 (%.2f 小时)\n", 
              as.numeric(total_elapsed), as.numeric(total_elapsed)/60))
  
  if (success_count > 0) {
    successful_results <- results[sapply(results, function(x) x$success)]
    avg_time <- mean(sapply(successful_results, function(x) x$processing_time))
    total_samples <- sum(sapply(successful_results, function(x) x$n_samples))
    
    cat(sprintf("📈 平均耗时: %.2f 分钟/文件\n", avg_time))
    cat(sprintf("📊 总样本数: %d\n", total_samples))
  }
  
  cat(sprintf("📁 输出目录: %s\n\n", config$output_base_dir))
  
  # 成功列表
  if (success_count > 0) {
    cat("✅ 成功处理的文件:\n")
    cat(sprintf("%-4s %-40s %10s %10s\n", "No.", "文件名", "耗时(分)", "样本数"))
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    j <- 1
    for (i in seq_along(results)) {
      if (results[[i]]$success) {
        cat(sprintf("%3d. %-40s %10.2f %10d\n", 
                    j,
                    results[[i]]$file,
                    results[[i]]$processing_time,
                    results[[i]]$n_samples))
        j <- j + 1
      }
    }
    cat("\n")
  }
  
  # 失败列表
  if (fail_count > 0) {
    cat("❌ 失败的文件:\n")
    cat(sprintf("%-4s %-40s %s\n", "No.", "文件名", "错误信息"))
    cat(paste(rep("-", 100), collapse = ""), "\n")
    
    j <- 1
    for (i in seq_along(results)) {
      if (!results[[i]]$success) {
        cat(sprintf("%3d. %-40s %s\n", 
                    j,
                    results[[i]]$file,
                    substr(results[[i]]$error, 1, 50)))
        j <- j + 1
      }
    }
    cat("\n")
  }
}


# ===================================================================
# 辅助函数：保存批量处理日志
# ===================================================================

save_batch_logs <- function(results, start_time, end_time, config) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # 文本日志
  log_file <- file.path(config$output_base_dir, 
                        sprintf("batch_processing_log_%s.txt", timestamp))
  
  # CSV 汇总
  csv_file <- file.path(config$output_base_dir, 
                        sprintf("batch_summary_%s.csv", timestamp))
  
  # 保存文本日志
  tryCatch({
    sink(log_file)
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("           Batch Processing Log\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    total_time <- difftime(end_time, start_time, units = "mins")
    
    cat(sprintf("Start time: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("End time:   %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Total time: %.2f minutes (%.2f hours)\n\n", 
                as.numeric(total_time), as.numeric(total_time)/60))
    
    # 详细结果
    for (i in seq_along(results)) {
      result <- results[[i]]
      status <- if(result$success) "SUCCESS" else "FAILED"
      
      cat(sprintf("[%s] File %2d/%d: %s\n", 
                  status, i, length(results), result$file))
      
      if (result$success) {
        cat(sprintf("           Time: %.2f min, Samples: %d\n", 
                    result$processing_time, result$n_samples))
      } else {
        cat(sprintf("           Error: %s\n", result$error))
      }
      cat("\n")
    }
    
    sink()
    
    cat(sprintf("📝 日志已保存:\n   %s\n", log_file))
    
  }, error = function(e) {
    sink()
    warning(sprintf("无法保存日志: %s", e$message))
  })
  
  # 保存 CSV
  tryCatch({
    summary_df <- data.frame(
      File_Number = seq_along(results),
      File_Name = sapply(results, function(x) x$file),
      Status = sapply(results, function(x) ifelse(x$success, "Success", "Failed")),
      Processing_Time_Minutes = sapply(results, function(x) round(x$processing_time, 2)),
      Number_of_Samples = sapply(results, function(x) x$n_samples),
      Error_Message = sapply(results, function(x) ifelse(!x$success, x$error, "")),
      stringsAsFactors = FALSE
    )
    
    write.csv(summary_df, csv_file, row.names = FALSE, quote = TRUE)
    cat(sprintf("📊 CSV已保存:\n   %s\n\n", csv_file))
    
  }, error = function(e) {
    warning(sprintf("无法保存CSV: %s", e$message))
  })
  
  return(list(log = log_file, csv = csv_file))
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