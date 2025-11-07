#!/usr/bin/env Rscript
# ===================================================================
# 样本预处理模块
# ===================================================================

#' 预处理样本（一次性切分所有样本）
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 要处理的样本列表
#' @param config 配置对象
#' 
#' @return 切分后的样本列表
#'
preprocess_samples <- function(seurat_obj, samples_to_plot, config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本预处理（统一切分）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 验证样本
  validation <- validate_samples(seurat_obj, samples_to_plot)
  samples_to_plot <- validation$valid_samples
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本")
  }
  
  print_preprocessing_info(seurat_obj, samples_to_plot)
  
  # 开始切分
  cat("🔧 切分样本...\n")
  start_time <- Sys.time()
  
  sample_list <- list()
  sample_stats <- initialize_sample_stats()
  
  # 切分每个样本
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    
    result <- split_single_sample(seurat_obj, sample_id, i, length(samples_to_plot))
    
    if (!is.null(result$seurat_subset)) {
      sample_list[[sample_id]] <- result$seurat_subset
      sample_stats <- rbind(sample_stats, result$stats)
    }
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 打印汇总
  print_preprocessing_summary(sample_list, sample_stats, elapsed, config)
  
  # 动态调整线程数
  recommended_workers <- calculate_safe_workers(sample_stats, config)
  
  # 保存统计到属性
  attr(sample_list, "stats") <- sample_stats
  attr(sample_list, "recommended_workers") <- recommended_workers
  
  return(sample_list)
}


#' 验证样本
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 样本列表
#' 
#' @return 验证结果
#'
validate_samples <- function(seurat_obj, samples_to_plot) {
  
  available_samples <- unique(seurat_obj$orig.ident)
  invalid_samples <- setdiff(samples_to_plot, available_samples)
  
  if (length(invalid_samples) > 0) {
    warning(sprintf("⚠️  以下样本不存在，将跳过: %s", 
                    paste(invalid_samples, collapse = ", ")))
  }
  
  valid_samples <- intersect(samples_to_plot, available_samples)
  
  return(list(
    valid_samples = valid_samples,
    invalid_samples = invalid_samples
  ))
}


#' 初始化样本统计表
#'
#' @return 空数据框
#'
initialize_sample_stats <- function() {
  data.frame(
    sample_id = character(),
    n_spots = integer(),
    n_high = integer(),
    high_pct = numeric(),
    size_mb = numeric(),
    stringsAsFactors = FALSE
  )
}


#' 切分单个样本
#'
#' @param seurat_obj Seurat 对象
#' @param sample_id 样本 ID
#' @param idx 当前索引
#' @param total 总数
#' 
#' @return 切分结果
#'
split_single_sample <- function(seurat_obj, sample_id, idx, total) {
  
  seurat_subset <- tryCatch({
    subset(seurat_obj, subset = orig.ident == sample_id)
  }, error = function(e) {
    seurat_obj[, seurat_obj$orig.ident == sample_id]
  })
  
  if (ncol(seurat_subset) == 0) {
    warning(sprintf("⚠️  样本 %s 无数据，已跳过", sample_id))
    return(list(seurat_subset = NULL, stats = NULL))
  }
  
  # 统计信息
  n_spots <- ncol(seurat_subset)
  n_high <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
  high_pct <- 100 * mean(seurat_subset$ClockGene_High, na.rm = TRUE)
  size_mb <- as.numeric(object.size(seurat_subset)) / 1024^2
  
  stats <- data.frame(
    sample_id = sample_id,
    n_spots = n_spots,
    n_high = n_high,
    high_pct = high_pct,
    size_mb = size_mb
  )
  
  cat(sprintf("  [%2d/%2d] ✅ %-30s | %6d spots | %4d high (%.2f%%) | %.2f MB\n",
              idx, total, sample_id, n_spots, n_high, high_pct, size_mb))
  
  return(list(seurat_subset = seurat_subset, stats = stats))
}


#' 打印预处理信息
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 样本列表
#'
print_preprocessing_info <- function(seurat_obj, samples_to_plot) {
  available_samples <- unique(seurat_obj$orig.ident)
  
  cat(sprintf("📊 原始数据: %d spots, %d 个样本\n", 
              ncol(seurat_obj), length(available_samples)))
  cat(sprintf("📊 将处理: %d 个样本\n\n", length(samples_to_plot)))
}


#' 打印预处理汇总
#'
#' @param sample_list 样本列表
#' @param sample_stats 统计数据
#' @param elapsed 耗时
#' @param config 配置对象
#'
print_preprocessing_summary <- function(sample_list, sample_stats, elapsed, config) {
  
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
  
  # 动态调整线程数建议
  max_memory_gb <- config$max_memory_gb %||% 100
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  recommended_workers <- min(safe_workers, length(sample_list), config$n_workers)
  
  cat(sprintf("\n💡 推荐线程数: %d (基于内存 %.0f GB)\n", 
              recommended_workers, max_memory_gb))
  
  if (recommended_workers < config$n_workers) {
    cat(sprintf("⚠️  原配置 %d 线程可能导致内存不足，已自动调整\n", 
                config$n_workers))
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 计算安全线程数
#'
#' @param sample_stats 样本统计
#' @param config 配置对象
#' 
#' @return 推荐线程数
#'
calculate_safe_workers <- function(sample_stats, config) {
  
  max_memory_gb <- config$max_memory_gb %||% 100
  avg_size_mb <- mean(sample_stats$size_mb)
  
  # 每个线程需要约 1.5 倍样本大小
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  
  # 不超过样本数和配置的线程数
  recommended_workers <- min(
    safe_workers, 
    nrow(sample_stats), 
    config$n_workers
  )
  
  return(max(1, recommended_workers))
}

cat("✅ 11_sample_preprocessing.R 已加载\n")