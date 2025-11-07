#!/usr/bin/env Rscript
# ===================================================================
# 空间梯度图绘制模块
# 功能：绘制 Clock Gene 距离场的空间分布图
# ===================================================================

#' 绘制空间梯度图
#'
#' @param sample_list 预切分的样本列表
#' @param CONFIG 配置对象
#' @param plot_width 图宽（英寸）
#' @param plot_height 图高（英寸）
#' @param pt_size_factor 点大小因子
#' @param alpha 透明度范围
#' @param color_option viridis 色板选项
#' @param color_direction 色板方向
#' 
#' @return 处理结果列表
#'
plot_spatial_gradient <- function(sample_list, 
                                  CONFIG,
                                  plot_width = 10,
                                  plot_height = 8,
                                  pt_size_factor = 1.5,
                                  alpha = c(0.3, 1),
                                  color_option = "magma",
                                  color_direction = -1) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   空间梯度图绘制（多线程并行）\n")
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
  n_workers <- CONFIG$n_workers %||% 4
  dpi <- CONFIG$plot$dpi %||% 300
  
  cat(sprintf("📊 将绘制 %d 个样本\n", length(sample_list)))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 2. 设置并行环境
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 3. 并行绘图
  # ========================================
  
  progressr::with_progress({
    
    p <- progressr::progressor(steps = length(sample_list))
    
    results <- future.apply::future_lapply(
      
      names(sample_list),
      
      function(sample_id) {
        
        tryCatch({
          
          # 获取样本数据
          seurat_subset <- sample_list[[sample_id]]
          
          # 验证数据
          if (ncol(seurat_subset) == 0) {
            p(message = sprintf("⚠️  %s - 无数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No data"
            ))
          }
          
          if (!"ClockGene_Distance" %in% colnames(seurat_subset@meta.data)) {
            p(message = sprintf("⚠️  %s - 缺少距离数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "Missing ClockGene_Distance column"
            ))
          }
          
          # 检查空间数据
          if (length(Seurat::Images(seurat_subset)) == 0) {
            p(message = sprintf("⚠️  %s - 无空间图像数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No spatial image data"
            ))
          }
          
          # 统计距离数据
          distance_values <- seurat_subset$ClockGene_Distance
          distance_stats <- list(
            min = min(distance_values, na.rm = TRUE),
            max = max(distance_values, na.rm = TRUE),
            mean = mean(distance_values, na.rm = TRUE),
            median = median(distance_values, na.rm = TRUE)
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
          
          p(message = sprintf("✅ %s", sample_id))
          
          return(list(
            sample = sample_id,
            success = TRUE,
            file = output_path,
            file_size_mb = file_size_mb,
            n_spots = n_spots,
            distance_stats = distance_stats
          ))
          
        }, error = function(e) {
          p(message = sprintf("❌ %s - %s", sample_id, e$message))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = as.character(e$message)
          ))
        })
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
  # 4. 统计输出
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d\n", n_success, length(sample_list)))
  
  if (n_failed > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", n_failed, length(sample_list)))
    cat("失败样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  • %s: %s\n", res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  if (n_success > 0) {
    cat("成功样本:\n")
    cat(sprintf("%-30s %10s %15s %15s %10s\n", 
                "样本", "Spots", "平均距离", "中位距离", "文件大小"))
    cat(paste(rep("-", 85), collapse = ""), "\n")
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %10d %15.2f %15.2f %8.2f MB\n",
                    res$sample,
                    res$n_spots,
                    res$distance_stats$mean,
                    res$distance_stats$median,
                    res$file_size_mb))
      }
    }
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 5. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
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