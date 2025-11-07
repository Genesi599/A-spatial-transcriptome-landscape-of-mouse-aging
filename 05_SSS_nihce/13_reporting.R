#!/usr/bin/env Rscript
# ===================================================================
# 报告生成模块
# ===================================================================

#' 打印批量处理头部
#'
print_batch_header <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║        Clock Gene Niche Analysis - Batch Processing       ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n")
}


#' 打印文件处理头部
#'
#' @param seurat_basename 文件基础名
#'
print_file_header <- function(seurat_basename) {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  处理文件: %-46s ║\n", seurat_basename))
  cat("╚═══════════════════════════════════════════════════════════╝\n")
}


#' 打印文件处理成功
#'
#' @param seurat_basename 文件基础名
#' @param n_samples 样本数
#' @param elapsed 耗时
#' @param config 配置对象
#'
print_file_success <- function(seurat_basename, n_samples, elapsed, config) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    处理完成                                ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("✅ 文件: %s\n", seurat_basename))
  cat(sprintf("📊 处理样本: %d\n", n_samples))
  cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(elapsed)))
  cat(sprintf("📁 输出: %s\n", config$output_dir))
  
  print_summary(config)
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 打印文件处理失败
#'
#' @param seurat_basename 文件基础名
#' @param error_msg 错误信息
#' @param elapsed 耗时
#'
print_file_failure <- function(seurat_basename, error_msg, elapsed) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    处理失败                                ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("❌ 文件: %s\n", seurat_basename))
  cat(sprintf("❌ 错误: %s\n", error_msg))
  cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(elapsed)))
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 打印批量处理总结
#'
#' @param results 结果列表
#' @param total_elapsed 总耗时
#' @param config 配置对象
#'
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
    print_success_statistics(results)
  }
  
  cat(sprintf("📁 输出目录: %s\n\n", config$output_base_dir))
  
  if (success_count > 0) {
    print_success_table(results)
  }
  
  if (fail_count > 0) {
    print_failure_table(results)
  }
}


#' 打印成功统计
#'
#' @param results 结果列表
#'
print_success_statistics <- function(results) {
  
  successful_results <- results[sapply(results, function(x) x$success)]
  avg_time <- mean(sapply(successful_results, function(x) x$processing_time))
  total_samples <- sum(sapply(successful_results, function(x) x$n_samples))
  
  cat(sprintf("📈 平均耗时: %.2f 分钟/文件\n", avg_time))
  cat(sprintf("📊 总样本数: %d\n", total_samples))
}


#' 打印成功文件表格
#'
#' @param results 结果列表
#'
print_success_table <- function(results) {
  
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


#' 打印失败文件表格
#'
#' @param results 结果列表
#'
print_failure_table <- function(results) {
  
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


#' 保存批量处理日志
#'
#' @param results 结果列表
#' @param start_time 开始时间
#' @param end_time 结束时间
#' @param config 配置对象
#' 
#' @return 日志文件路径
#'
save_batch_logs <- function(results, start_time, end_time, config) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  log_file <- file.path(config$output_base_dir, 
                       sprintf("batch_processing_log_%s.txt", timestamp))
  
  csv_file <- file.path(config$output_base_dir, 
                       sprintf("batch_summary_%s.csv", timestamp))
  
  # 保存文本日志
  save_text_log(results, start_time, end_time, log_file)
  
  # 保存 CSV
  save_csv_summary(results, csv_file)
  
  return(list(log = log_file, csv = csv_file))
}


#' 保存文本日志
#'
#' @param results 结果列表
#' @param start_time 开始时间
#' @param end_time 结束时间
#' @param log_file 日志文件路径
#'
save_text_log <- function(results, start_time, end_time, log_file) {
  
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
    warning(sprintf("⚠️  无法保存日志: %s", e$message))
  })
}


#' 保存 CSV 汇总
#'
#' @param results 结果列表
#' @param csv_file CSV 文件路径
#'
save_csv_summary <- function(results, csv_file) {
  
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
    warning(sprintf("⚠️  无法保存CSV: %s", e$message))
  })
}

cat("✅ 13_reporting.R 已加载\n")