#!/usr/bin/env Rscript
# ===================================================================
# 批量处理核心模块
# ===================================================================

#' 处理所有文件
#'
#' @param seurat_files 文件列表
#' @param gene_list 基因列表
#' @param CONFIG 配置对象
#' 
#' @return 处理结果列表
#'
process_all_files <- function(seurat_files, gene_list, CONFIG) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   开始批量处理\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  results <- list()
  
  for (i in seq_along(seurat_files)) {
    
    # 打印进度
    print_progress_header(i, length(seurat_files))
    
    # 处理单个文件
    result <- process_seurat_file(
      seurat_path = seurat_files[i],
      gene_list = gene_list,
      base_config = CONFIG
    )
    
    results[[i]] <- result
    
    # 估计剩余时间
    if (i < length(seurat_files)) {
      estimate_remaining_time(results, seurat_files, i)
    }
    
    # 强制内存清理
    gc(verbose = FALSE)
  }
  
  return(results)
}


#' 打印进度头部
#'
#' @param current 当前序号
#' @param total 总数
#'
print_progress_header <- function(current, total) {
  cat(sprintf("\n╔═══════════════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  进度: [%d/%d] (%.1f%%)%*s║\n", 
              current, total, (current/total)*100,
              60 - nchar(sprintf("  进度: [%d/%d] (%.1f%%)", current, total, (current/total)*100)), ""))
  cat(sprintf("╚═══════════════════════════════════════════════════════════╝\n"))
}


#' 估算剩余时间
#'
#' @param results 已完成的结果
#' @param seurat_files 所有文件
#' @param current_idx 当前索引
#'
estimate_remaining_time <- function(results, seurat_files, current_idx) {
  avg_time <- mean(sapply(results, function(x) x$processing_time), na.rm = TRUE)
  remaining_time <- avg_time * (length(seurat_files) - current_idx)
  
  cat(sprintf("\n📊 预计剩余时间: %.2f 分钟 (%.2f 小时)\n", 
              remaining_time, remaining_time/60))
}


#' 确认批量处理
#'
#' @param seurat_files 文件列表
#' @param CONFIG 配置对象
#' 
#' @return 逻辑值
#'
confirm_batch_processing <- function(seurat_files, CONFIG) {
  
  if (!CONFIG$batch_mode) return(TRUE)
  if (length(seurat_files) <= 1) return(TRUE)
  if (!interactive()) return(TRUE)
  
  response <- readline(prompt = sprintf(
    "即将处理 %d 个文件，是否继续? (y/n): ", length(seurat_files)))
  
  return(tolower(response) == "y")
}


#' 创建汇总对象
#'
#' @param results 结果列表
#' @param total_elapsed 总耗时
#' @param log_files 日志文件
#' 
#' @return 汇总对象
#'
create_summary_object <- function(results, total_elapsed, log_files) {
  list(
    total = length(results),
    success = sum(sapply(results, function(x) x$success)),
    failed = sum(sapply(results, function(x) !x$success)),
    total_time = as.numeric(total_elapsed),
    log_file = log_files$log,
    csv_file = log_files$csv
  )
}

cat("✅ 10_batch_processing.R 已加载\n")