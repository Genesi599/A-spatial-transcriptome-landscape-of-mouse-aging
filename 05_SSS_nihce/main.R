# main.R

#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis - Main Script (Batch Processing)
# Author: Zhangbin (optimized)
# Date: 2024-11-04
# Modified: 2024-11-06 - Added batch processing capability
# ===================================================================

# 加载配置
source("00_config.R")

# 加载模块
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
# 单个文件处理函数
# ===================================================================

process_single_file <- function(seurat_path, gene_list, base_config) {
  # 更新配置
  config <- base_config
  config$seurat_path <- seurat_path
  
  # 提取文件名并生成输出目录
  seurat_basename <- tools::file_path_sans_ext(basename(seurat_path))
  config$output_dir <- file.path(config$output_base_dir, seurat_basename)
  
  # 更新所有目录路径
  config$cache_dir <- file.path(config$output_dir, "cache")
  config$figure_dir <- file.path(config$output_dir, "figure")
  config$metadata_dir <- file.path(config$output_dir, "metadata")
  
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
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("📂 开始处理: %s\n", seurat_basename))
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  tryCatch({
    # ✅ 1. 环境设置（创建所有必要的目录）
    cat("🔧 设置环境...\n")
    setup_environment(config)
    cat("✓ 环境设置完成\n\n")
    
    # 2. 加载数据
    cat(sprintf("📥 加载 Seurat 对象...\n"))
    seurat_obj <- load_seurat_object(config)
    genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
    
    # 3. 计算评分
    seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, config)
    result <- define_high_expression(seurat_obj, config)
    seurat_obj <- result$seurat_obj
    threshold <- result$threshold
    
    # 4. Niche分析
    seurat_obj <- perform_niche_analysis(seurat_obj, threshold, config)
    
    # 5. 确定要处理的样本
    samples <- unique(seurat_obj$orig.ident)
    samples_to_plot <- if (config$debug_mode) {
      head(samples, config$debug_sample_limit)
    } else {
      samples
    }
    cat(sprintf("📋 将处理 %d 个样本\n\n", length(samples_to_plot)))
    
    # 6. 绘图
    plot_isoheight(seurat_obj, samples_to_plot, config)
    plot_spatial_gradient(seurat_obj, samples_to_plot, config)
    analyze_celltype_niche(seurat_obj, samples_to_plot, config)
    
    # 7. 保存结果
    save_results(seurat_obj, config)
    print_summary(config)
    
    cat("\n")
    cat(sprintf("✅ 成功完成: %s\n", seurat_basename))
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    return(list(success = TRUE, file = seurat_basename, error = NULL))
    
  }, error = function(e) {
    cat("\n")
    cat(sprintf("❌ 处理失败: %s\n", seurat_basename))
    cat(sprintf("错误信息: %s\n", e$message))
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    return(list(success = FALSE, file = seurat_basename, error = e$message))
  })
}

# ===================================================================
# 批量处理主函数（改进版）
# ===================================================================

main_batch <- function() {
  start_time <- Sys.time()
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║        Clock Gene Niche Analysis - Batch Processing       ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  # 0. 加载必要的包和函数
  cat("🔧 加载包和函数...\n")
  load_packages()
  load_custom_functions()
  cat("✓ 完成\n\n")
  

    # ✅ 0.5 确保输出基础目录存在
  if (!is.null(CONFIG$output_base_dir) && CONFIG$output_base_dir != "") {
    if (!dir.exists(CONFIG$output_base_dir)) {
      cat(sprintf("📁 创建输出基础目录: %s\n", CONFIG$output_base_dir))
      dir.create(CONFIG$output_base_dir, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    stop("❌ 未配置 output_base_dir")
  }
  
  # 1. 扫描输入目录
  if (CONFIG$batch_mode) {
    # 批量模式：扫描目录
    cat(sprintf("📁 扫描目录: %s\n", CONFIG$seurat_dir))
    
    # 检查目录是否存在
    if (!dir.exists(CONFIG$seurat_dir)) {
      stop(sprintf("❌ 目录不存在: %s", CONFIG$seurat_dir))
    }
    
    seurat_files <- list.files(
      path = CONFIG$seurat_dir,
      pattern = CONFIG$seurat_pattern,
      full.names = TRUE,
      recursive = CONFIG$recursive_search
    )
    
    if (length(seurat_files) == 0) {
      stop(sprintf("❌ 在目录 %s 中未找到匹配的文件 (模式: %s)", 
                   CONFIG$seurat_dir, CONFIG$seurat_pattern))
    }
    
    cat(sprintf("✓ 找到 %d 个 Seurat 文件\n\n", length(seurat_files)))
    
    # 应用文件过滤
    if (!is.null(CONFIG$specific_files) || !is.null(CONFIG$exclude_files)) {
      original_count <- length(seurat_files)
      seurat_files <- filter_seurat_files(seurat_files, CONFIG)
      cat(sprintf("📋 过滤后剩余 %d 个文件 (原始: %d)\n\n", 
                  length(seurat_files), original_count))
    }
    
    # 打印文件列表
    cat("将处理以下文件:\n")
    for (i in seq_along(seurat_files)) {
      file_size <- file.size(seurat_files[i]) / (1024^3)  # GB
      cat(sprintf("  %2d. %-50s (%.2f GB)\n", 
                  i, basename(seurat_files[i]), file_size))
    }
    cat("\n")
    
  } else {
    # 单文件模式
    if (!file.exists(CONFIG$seurat_path)) {
      stop(sprintf("❌ 文件不存在: %s", CONFIG$seurat_path))
    }
    seurat_files <- CONFIG$seurat_path
    file_size <- file.size(seurat_files) / (1024^3)
    cat(sprintf("📄 单文件模式: %s (%.2f GB)\n\n", 
                basename(seurat_files), file_size))
  }
  
  # 确认是否继续
  if (CONFIG$batch_mode && length(seurat_files) > 1 && interactive()) {
    response <- readline(prompt = sprintf(
      "即将处理 %d 个文件，是否继续? (y/n): ", length(seurat_files)))
    if (tolower(response) != "y") {
      cat("已取消处理\n")
      return(invisible(NULL))
    }
  }
  
  # 2. 加载基因列表（只加载一次）
  cat("📋 加载基因列表...\n")
  gene_list <- load_gene_list(CONFIG)
  cat(sprintf("✓ 加载了 %d 个基因\n\n", length(gene_list)))
  
  # 3. 批量处理
  results <- list()
  successful_files <- c()
  failed_files <- c()
  
  for (i in seq_along(seurat_files)) {
    file_start_time <- Sys.time()
    
    cat(sprintf("\n═══════════════════════════════════════════════════════════\n"))
    cat(sprintf("进度: [%d/%d] (%.1f%%)\n", 
                i, length(seurat_files), (i/length(seurat_files))*100))
    cat(sprintf("═══════════════════════════════════════════════════════════\n"))
    
    result <- process_single_file(
      seurat_path = seurat_files[i],
      gene_list = gene_list,
      base_config = CONFIG
    )
    
    file_end_time <- Sys.time()
    file_time <- difftime(file_end_time, file_start_time, units = "mins")
    
    result$processing_time <- as.numeric(file_time)
    results[[i]] <- result
    
    # 记录成功/失败的文件
    if (result$success) {
      successful_files <- c(successful_files, basename(seurat_files[i]))
      cat(sprintf("\n⏱️  文件处理耗时: %.2f 分钟\n", as.numeric(file_time)))
    } else {
      failed_files <- c(failed_files, basename(seurat_files[i]))
    }
    
    # 估计剩余时间
    if (i < length(seurat_files)) {
      avg_time <- mean(sapply(results, function(x) x$processing_time), na.rm = TRUE)
      remaining_time <- avg_time * (length(seurat_files) - i)
      cat(sprintf("📊 预计剩余时间: %.2f 分钟\n", remaining_time))
    }
    
    # 清理内存
    gc(verbose = FALSE)
  }
  
  # 4. 生成总结报告
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  cat("\n\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    处理完成总结                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  success_count <- sum(sapply(results, function(x) x$success))
  fail_count <- length(results) - success_count
  
  cat(sprintf("📊 总文件数: %d\n", length(results)))
  cat(sprintf("✅ 成功: %d (%.1f%%)\n", 
              success_count, (success_count/length(results))*100))
  cat(sprintf("❌ 失败: %d (%.1f%%)\n", 
              fail_count, (fail_count/length(results))*100))
  cat(sprintf("⏱️  总耗时: %.2f 分钟 (%.2f 小时)\n", 
              as.numeric(total_time), as.numeric(total_time)/60))
  
  if (success_count > 0) {
    avg_time <- mean(sapply(results[sapply(results, function(x) x$success)], 
                           function(x) x$processing_time), na.rm = TRUE)
    cat(sprintf("📈 平均每个文件: %.2f 分钟\n", avg_time))
  }
  
  cat(sprintf("📁 输出目录: %s\n\n", CONFIG$output_base_dir))
  
  # 显示成功的文件
  if (success_count > 0) {
    cat("✅ 成功处理的文件:\n")
    for (i in seq_along(results)) {
      if (results[[i]]$success) {
        cat(sprintf("  %2d. %-50s (%.2f 分钟)\n", 
                    i, results[[i]]$file, results[[i]]$processing_time))
      }
    }
    cat("\n")
  }
  
  # 显示失败的文件
  if (fail_count > 0) {
    cat("❌ 失败的文件:\n")
    for (i in seq_along(results)) {
      if (!results[[i]]$success) {
        cat(sprintf("  %2d. %s\n", i, results[[i]]$file))
        cat(sprintf("      错误: %s\n", results[[i]]$error))
      }
    }
    cat("\n")
  }
  
  # 5. 保存处理日志
  log_file <- file.path(CONFIG$output_base_dir, 
                        sprintf("batch_processing_log_%s.txt", 
                                format(Sys.time(), "%Y%m%d_%H%M%S")))
  
  # 确保输出目录存在
  if (!dir.exists(CONFIG$output_base_dir)) {
    dir.create(CONFIG$output_base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tryCatch({
    sink(log_file)
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("           Batch Processing Log\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    cat(sprintf("Start time: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("End time:   %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Total time: %.2f minutes (%.2f hours)\n\n", 
                as.numeric(total_time), as.numeric(total_time)/60))
    
    cat("───────────────────────────────────────────────────────────\n")
    cat("Summary Statistics\n")
    cat("───────────────────────────────────────────────────────────\n\n")
    
    cat(sprintf("Total files:     %d\n", length(results)))
    cat(sprintf("Successful:      %d (%.1f%%)\n", 
                success_count, (success_count/length(results))*100))
    cat(sprintf("Failed:          %d (%.1f%%)\n", 
                fail_count, (fail_count/length(results))*100))
    
    if (success_count > 0) {
      avg_time <- mean(sapply(results[sapply(results, function(x) x$success)], 
                             function(x) x$processing_time), na.rm = TRUE)
      cat(sprintf("Average time:    %.2f minutes per file\n", avg_time))
    }
    
    cat("\n───────────────────────────────────────────────────────────\n")
    cat("Configuration\n")
    cat("───────────────────────────────────────────────────────────\n\n")
    
    cat(sprintf("Input directory:     %s\n", CONFIG$seurat_dir))
    cat(sprintf("Pattern:             %s\n", CONFIG$seurat_pattern))
    cat(sprintf("Recursive search:    %s\n", CONFIG$recursive_search))
    cat(sprintf("Threshold quantile:  %.2f\n", CONFIG$threshold_quantile))
    cat(sprintf("Number of workers:   %d\n", CONFIG$n_workers))
    cat(sprintf("Debug mode:          %s\n", CONFIG$debug_mode))
    
    cat("\n───────────────────────────────────────────────────────────\n")
    cat("Detailed Results\n")
    cat("───────────────────────────────────────────────────────────\n\n")
    
    for (i in seq_along(results)) {
      result <- results[[i]]
      status <- if(result$success) "SUCCESS" else "FAILED "
      
      cat(sprintf("[%s] File %2d/%d: %s\n", 
                  status, i, length(results), result$file))
      
      if (result$success) {
        cat(sprintf("           Processing time: %.2f minutes\n", 
                    result$processing_time))
      } else {
        cat(sprintf("           Error: %s\n", result$error))
      }
      cat("\n")
    }
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("End of Log\n")
    cat("═══════════════════════════════════════════════════════════\n")
    
    sink()
    
    cat(sprintf("📝 处理日志已保存至:\n   %s\n\n", log_file))
    
  }, error = function(e) {
    sink()
    warning(sprintf("无法保存日志文件: %s", e$message))
  })
  
  # 6. 生成 CSV 汇总报告
  csv_file <- file.path(CONFIG$output_base_dir, 
                        sprintf("batch_summary_%s.csv", 
                                format(Sys.time(), "%Y%m%d_%H%M%S")))
  
  tryCatch({
    summary_df <- data.frame(
      File_Number = seq_along(results),
      File_Name = sapply(results, function(x) x$file),
      Status = sapply(results, function(x) ifelse(x$success, "Success", "Failed")),
      Processing_Time_Minutes = sapply(results, function(x) {
        if (!is.null(x$processing_time)) round(x$processing_time, 2) else NA
      }),
      Error_Message = sapply(results, function(x) {
        if (!x$success) x$error else ""
      }),
      stringsAsFactors = FALSE
    )
    
    write.csv(summary_df, csv_file, row.names = FALSE, quote = TRUE)
    cat(sprintf("📊 CSV汇总已保存至:\n   %s\n\n", csv_file))
    
  }, error = function(e) {
    warning(sprintf("无法保存CSV文件: %s", e$message))
  })
  
  # 7. 返回结果
  cat("🎉 批量处理完成！\n\n")
  
  return(invisible(list(
    results = results,
    summary = list(
      total = length(results),
      success = success_count,
      failed = fail_count,
      total_time = as.numeric(total_time),
      log_file = log_file,
      csv_file = csv_file
    )
  )))
}

# ===================================================================
# 运行主流程
# ===================================================================

if (!interactive()) {
  main_batch()
}