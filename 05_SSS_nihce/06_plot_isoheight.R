#!/usr/bin/env Rscript
# ===================================================================
# 等高线密度图绘制模块（优化版 + 进度条）
# 功能：多线程并行绘制 Clock Gene 等高线密度图
# ===================================================================

library(future)
library(future.apply)
library(progressr)
library(ggplot2)


#' 绘制等高线密度图（接收预切分样本）
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param col_bg 背景点颜色
#' @param col_top 高表达点颜色
#' @param col_isoheight 等高线颜色
#' @param col_white_ratio 白色占比
#' @param cols_fill_isoheight 填充色谱
#' @param plot_width 图宽
#' @param plot_height 图高
#' @param nrow 子图排列行数
#' 
#' @return 处理结果列表
#'
plot_isoheight <- function(sample_list, 
                          CONFIG,
                          col_bg = "gray92",
                          col_top = "#d62728",
                          col_isoheight = "white",
                          col_white_ratio = 0.25,
                          cols_fill_isoheight = NULL,
                          plot_width = 8,
                          plot_height = 8,
                          nrow = 1) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   等高线密度图绘制（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  if (is.null(CONFIG$dirs$isoheight)) {
    stop("❌ CONFIG$dirs$isoheight 未定义")
  }
  
  if (!exists("celltype_isoheight_plot")) {
    stop("❌ 未找到 celltype_isoheight_plot 函数，请先加载 SSS_isoheight_plot.R")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$isoheight)) {
    dir.create(CONFIG$dirs$isoheight, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 提取参数
  size_bg <- CONFIG$plot$point_size_bg %||% 0.3
  size_top <- CONFIG$plot$point_size_top %||% 1.2
  dpi <- CONFIG$plot$dpi %||% 300
  n_workers <- CONFIG$n_workers %||% 4
  
  # 默认色谱
  if (is.null(cols_fill_isoheight)) {
    cols_fill_isoheight <- c(
      rep("white", 25),
      colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])(75)
    )
  }
  
  cat(sprintf("📊 将绘制 %d 个样本\n", length(sample_list)))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 2. 设置并行环境
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  
  # ========================================
  # 4. 并行绘图
  # ========================================
  
  cat("🎨 开始绘图...\n\n")
  
  progressr::with_progress({
    
    p <- progressr::progressor(
      steps = length(sample_list),
      message = "绘制等高线图"
    )
    
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
          
          if (!"ClockGene_High" %in% colnames(seurat_subset@meta.data)) {
            p(message = sprintf("⚠️  %s - 缺少 ClockGene_High 列", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "Missing ClockGene_High column"
            ))
          }
          
          n_high <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
          
          if (n_high == 0) {
            p(message = sprintf("⚠️  %s - 无高表达点", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No high expression spots"
            ))
          }
          
          # 绘图
          p_iso <- celltype_isoheight_plot(
            .data = seurat_subset,
            density_top = ClockGene_High,
            col_bg = col_bg,
            col_top = col_top,
            col_isoheight = col_isoheight,
            col_white_ratio = col_white_ratio,
            cols_fill_isoheight = cols_fill_isoheight,
            size_bg = size_bg,
            size_top = size_top,
            nrow = nrow
          )
          
          # 保存
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          output_file <- sprintf("ClockGene_isoheight_%s.pdf", safe_name)
          output_path <- file.path(CONFIG$dirs$isoheight, output_file)
          
          ggplot2::ggsave(
            filename = output_path,
            plot = p_iso, 
            width = plot_width, 
            height = plot_height, 
            dpi = dpi
          )
          
          # 统计
          file_size_mb <- file.size(output_path) / 1024^2
          n_spots <- ncol(seurat_subset)
          high_pct <- 100 * n_high / n_spots
          
          # 更新进度（显示样本名）
          p(message = sprintf("✅ %s (%.2f MB)", sample_id, file_size_mb))
          
          return(list(
            sample = sample_id,
            success = TRUE,
            file = output_path,
            file_size_mb = file_size_mb,
            n_spots = n_spots,
            n_high = n_high,
            high_pct = high_pct
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
  # 5. 统计输出
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, 
              length(sample_list),
              100 * n_success / length(sample_list)))
  
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
    cat(sprintf("%-30s %10s %10s %10s %10s\n", 
                "样本", "Spots", "High", "High%", "文件大小"))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    total_file_size <- 0
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %10d %10d %9.1f%% %8.2f MB\n",
                    res$sample,
                    res$n_spots,
                    res$n_high,
                    res$high_pct,
                    res$file_size_mb))
        total_file_size <- total_file_size + res$file_size_mb
      }
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("%-30s %10s %10s %10s %8.2f MB\n",
                "总计", "", "", "", total_file_size))
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$isoheight))
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 6. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    output_dir = CONFIG$dirs$isoheight,
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

cat("✅ 06_plot_isoheight.R 已加载\n")