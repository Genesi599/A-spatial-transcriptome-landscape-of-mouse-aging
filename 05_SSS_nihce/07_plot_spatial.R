#!/usr/bin/env Rscript
# ===================================================================
# 空间梯度图绘制模块（串联版 - 无并行依赖）
# 功能：绘制 Clock Gene 空间距离梯度图
# ===================================================================

library(Seurat)
library(ggplot2)


#' 绘制空间梯度图（接收预切分样本）
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param pt_size_factor 点大小因子
#' @param alpha 透明度
#' @param color_option viridis 色谱选项
#' @param color_direction 色谱方向
#' @param plot_width 图宽
#' @param plot_height 图高
#' 
#' @return 处理结果列表
#'
plot_spatial <- function(sample_list,
                        CONFIG,
                        pt_size_factor = 1.6,
                        alpha = 0.8,
                        color_option = "magma",
                        color_direction = -1,
                        plot_width = 10,
                        plot_height = 8) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   空间梯度图绘制\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  if (is.null(CONFIG$dirs$spatial)) {
    stop("❌ CONFIG$dirs$spatial 未定义")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$spatial)) {
    dir.create(CONFIG$dirs$spatial, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 提取参数
  dpi <- CONFIG$plot$dpi %||% 300
  
  cat(sprintf("📊 将绘制 %d 个样本\n\n", length(sample_list)))
  
  start_time <- Sys.time()
  
  # ========================================
  # 2. 串联绘图
  # ========================================
  
  cat("🗺️  开始绘图...\n\n")
  
  success_list <- list()
  failed_list <- list()
  total_samples <- length(sample_list)
  
  for (i in seq_along(sample_list)) {
    
    sample_id <- names(sample_list)[i]
    
    cat(sprintf("[%2d/%2d] ", i, total_samples))
    
    tryCatch({
      
      # 获取样本数据
      seurat_subset <- sample_list[[sample_id]]
      
      # 验证数据
      if (ncol(seurat_subset) == 0) {
        cat(sprintf("⚠️  %s - 无数据\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "No data"
        )
        next
      }
      
      if (!"ClockGene_Distance" %in% colnames(seurat_subset@meta.data)) {
        cat(sprintf("⚠️  %s - 缺少距离数据\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "Missing ClockGene_Distance column"
        )
        next
      }
      
      # 检查空间数据
      if (length(Seurat::Images(seurat_subset)) == 0) {
        cat(sprintf("⚠️  %s - 无空间图像数据\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "No spatial image data"
        )
        next
      }
      
      # 统计距离数据
      distance_values <- seurat_subset$ClockGene_Distance
      distance_stats <- list(
        min = min(distance_values, na.rm = TRUE),
        max = max(distance_values, na.rm = TRUE),
        mean = mean(distance_values, na.rm = TRUE),
        median = median(distance_values, na.rm = TRUE),
        sd = sd(distance_values, na.rm = TRUE)
      )
      
      # 绘制空间分布图
      p_spatial <- Seurat::SpatialFeaturePlot(
        seurat_subset,
        features = "ClockGene_Distance",
        pt.size.factor = pt_size_factor,
        alpha = alpha,
        stroke = 0
      ) + 
        ggplot2::scale_fill_viridis_c(
          option = color_option,
          direction = color_direction,
          name = "Distance\nto High",
          limits = c(0, NA)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "right",
          legend.title = ggplot2::element_text(size = 12, face = "bold"),
          legend.text = ggplot2::element_text(size = 10),
          plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5)
        ) +
        ggplot2::ggtitle(
          sprintf("Clock Gene Distance Field - %s", sample_id),
          subtitle = sprintf(
            "Mean: %.2f | Median: %.2f | Range: [%.2f, %.2f]",
            distance_stats$mean,
            distance_stats$median,
            distance_stats$min,
            distance_stats$max
          )
        )
      
      # 保存图形
      safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
      output_file <- sprintf("ClockGene_spatial_%s.pdf", safe_name)
      output_path <- file.path(CONFIG$dirs$spatial, output_file)
      
      ggplot2::ggsave(
        filename = output_path,
        plot = p_spatial,
        width = plot_width,
        height = plot_height,
        dpi = dpi
      )
      
      # 统计信息
      file_size_mb <- file.size(output_path) / 1024^2
      n_spots <- ncol(seurat_subset)
      
      # 输出成功信息
      cat(sprintf("✅ %s (%.2f MB, %d spots, dist: %.2f±%.2f)\n", 
                 sample_id, file_size_mb, n_spots, 
                 distance_stats$mean, distance_stats$sd))
      
      success_list[[sample_id]] <- list(
        sample = sample_id,
        success = TRUE,
        file = output_path,
        file_size_mb = file_size_mb,
        n_spots = n_spots,
        distance_stats = distance_stats
      )
      
      # 清理内存
      rm(seurat_subset, p_spatial)
      if (i %% 3 == 0) gc(verbose = FALSE)
      
    }, error = function(e) {
      cat(sprintf("❌ %s - %s\n", sample_id, e$message))
      failed_list[[sample_id]] <- list(
        sample = sample_id,
        success = FALSE,
        error = as.character(e$message)
      )
    })
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 合并结果
  results <- c(success_list, failed_list)
  
  # ========================================
  # 3. 统计输出
  # ========================================
  
  n_success <- length(success_list)
  n_failed <- length(failed_list)
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, 
              total_samples,
              100 * n_success / total_samples))
  
  if (n_failed > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", n_failed, total_samples))
    cat("失败样本:\n")
    for (sample_id in names(failed_list)) {
      res <- failed_list[[sample_id]]
      cat(sprintf("  • %s: %s\n", res$sample, res$error))
    }
    cat("\n")
  }
  
  if (n_success > 0) {
    cat("成功样本:\n")
    cat(sprintf("%-30s %10s %12s %12s %10s\n", 
                "样本", "Spots", "Mean Dist", "SD Dist", "文件大小"))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    total_file_size <- 0
    
    for (sample_id in names(success_list)) {
      res <- success_list[[sample_id]]
      cat(sprintf("%-30s %10d %12.2f %12.2f %8.2f MB\n",
                  res$sample,
                  res$n_spots,
                  res$distance_stats$mean,
                  res$distance_stats$sd,
                  res$file_size_mb))
      total_file_size <- total_file_size + res$file_size_mb
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("%-30s %10s %12s %12s %8.2f MB\n",
                "总计", "", "", "", total_file_size))
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / total_samples))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 4. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = total_samples,
    output_dir = CONFIG$dirs$spatial,
    elapsed_time = as.numeric(elapsed),
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 07_plot_spatial.R 已加载\n")