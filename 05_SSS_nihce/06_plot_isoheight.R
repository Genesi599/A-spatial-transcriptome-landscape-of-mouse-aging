# 06_plot_isoheight.R (多线程并行版)

# ===================================================================
# 函数：绘制 Isoheight 图（多线程并行版）
# 作者：Assistant (优化版)
# 日期：2025-11-05
# 更新：2025-11-06 - 添加多线程并行支持
# ===================================================================

library(future)
library(future.apply)

#' 绘制 Clock Gene High 表达的 Isoheight 密度图（并行版）
#'
#' @param seurat_obj Seurat 对象，必须包含以下列：
#'   - ClockGene_High: 布尔值，标记高表达点
#'   - orig.ident: 样本ID
#' @param samples_to_plot 字符向量，要绘制的样本ID列表
#' @param CONFIG 配置列表，必须包含：
#'   - dirs$isoheight: 输出目录路径
#'   - plot$point_size_bg: 背景点大小（可选，默认0.3）
#'   - plot$point_size_top: 高表达点大小（可选，默认1.2）
#'   - plot$dpi: 图形分辨率（可选，默认300）
#'   - plot$contour_bins: 等高线数量（可选，默认8）
#'   - n_workers: 并行线程数（可选，默认4）
#' @param col_bg 背景点颜色，默认 "gray92"
#' @param col_top 高表达点颜色，默认 "#d62728"（红色）
#' @param col_isoheight 等高线颜色，默认 "white"
#' @param col_white_ratio 白色占比，默认 0.25
#' @param cols_fill_isoheight 等高线填充颜色向量，默认为黄-橙-红渐变
#' @param plot_width 图形宽度，默认 8
#' @param plot_height 图形高度，默认 8
#' @param nrow 图形布局行数，默认 1
#'
#' @return 返回结果列表（隐式），包含成功/失败统计
#'
#' @examples
#' # 基础调用（使用默认颜色）
#' plot_isoheight(seurat_obj, samples_to_plot, CONFIG)
#' 
#' # 自定义颜色和线程数
#' plot_isoheight(
#'   seurat_obj, samples_to_plot, CONFIG,
#'   col_bg = "lightgray",
#'   col_top = "darkred"
#' )
#'
plot_isoheight <- function(seurat_obj, 
                          samples_to_plot, 
                          CONFIG,
                          col_bg = "gray92",
                          col_top = "#d62728",
                          col_isoheight = "white",
                          col_white_ratio = 0.25,
                          cols_fill_isoheight = NULL,
                          plot_width = 8,
                          plot_height = 8,
                          nrow = 1) {
  
  # ========================================
  # 1. 参数验证
  # ========================================
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   等高线图绘制（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 检查必需的列
  required_cols <- c("ClockGene_High", "orig.ident")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 检查 ClockGene_High 是否为布尔值
  if (!is.logical(seurat_obj$ClockGene_High)) {
    stop("❌ ClockGene_High 必须是逻辑值（TRUE/FALSE）")
  }
  
  # 检查必需的配置
  if (is.null(CONFIG$dirs$isoheight)) {
    stop("❌ CONFIG$dirs$isoheight 未定义")
  }
  
  # 检查 celltype_isoheight_plot 函数是否可用
  if (!exists("celltype_isoheight_plot")) {
    stop("❌ 未找到 celltype_isoheight_plot 函数，请先加载 SSS_isoheight_plot.R")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$isoheight)) {
    dir.create(CONFIG$dirs$isoheight, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("✅ 创建输出目录: %s\n", CONFIG$dirs$isoheight))
  }
  
  # 提取配置参数（设置默认值）
  size_bg <- CONFIG$plot$point_size_bg %||% 0.3
  size_top <- CONFIG$plot$point_size_top %||% 1.2
  dpi <- CONFIG$plot$dpi %||% 300
  n_workers <- CONFIG$n_workers %||% 4
  
  # 默认颜色方案
  if (is.null(cols_fill_isoheight)) {
    cols_fill_isoheight <- c(
      rep("white", 25),
      colorRampPalette(brewer.pal(9, "YlOrRd")[3:9])(75)
    )
  }
  
  # ========================================
  # 2. 样本验证
  # ========================================
  available_samples <- unique(seurat_obj$orig.ident)
  invalid_samples <- setdiff(samples_to_plot, available_samples)
  
  if (length(invalid_samples) > 0) {
    warning(sprintf("⚠️ 以下样本不存在，将跳过: %s", 
                    paste(invalid_samples, collapse = ", ")))
    samples_to_plot <- intersect(samples_to_plot, available_samples)
  }
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本可绘制")
  }
  
  cat(sprintf("📊 将绘制 %d 个样本\n", length(samples_to_plot)))
  
  # 统计高表达点
  high_count <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
  high_pct <- 100 * high_count / ncol(seurat_obj)
  cat(sprintf("📊 高表达点: %d / %d (%.2f%%)\n", 
              high_count, ncol(seurat_obj), high_pct))
  
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 3. 设置并行计划
  # ========================================
  plan(multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 4. 并行处理每个样本
  # ========================================
  results <- future_lapply(seq_along(samples_to_plot), function(i) {
    
    sample_id <- samples_to_plot[i]
    
    result <- tryCatch({
      
      # --------------------------------
      # 4.1 提取子集
      # --------------------------------
      seurat_subset <- tryCatch(
        subset(seurat_obj, subset = orig.ident == sample_id),
        error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
      )
      
      if (ncol(seurat_subset) == 0) {
        return(list(
          sample = sample_id,
          index = i,
          success = FALSE,
          error = "No data for this sample"
        ))
      }
      
      # 统计该样本的高表达点
      sample_high_count <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
      sample_high_pct <- 100 * sample_high_count / ncol(seurat_subset)
      
      if (sample_high_count == 0) {
        return(list(
          sample = sample_id,
          index = i,
          success = FALSE,
          error = "No high expression spots"
        ))
      }
      
      # --------------------------------
      # 4.2 调用 celltype_isoheight_plot
      # --------------------------------
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
      
      # --------------------------------
      # 4.3 保存图形
      # --------------------------------
      safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
      output_path <- file.path(
        CONFIG$dirs$isoheight, 
        sprintf("ClockGene_isoheight_%s.pdf", safe_name)
      )
      
      ggsave(
        filename = output_path,
        plot = p_iso, 
        width = plot_width, 
        height = plot_height, 
        dpi = dpi
      )
      
      file_size_mb <- file.size(output_path) / 1024^2
      
      return(list(
        sample = sample_id,
        index = i,
        success = TRUE,
        file = output_path,
        file_size_mb = file_size_mb,
        n_spots = ncol(seurat_subset),
        n_high = sample_high_count,
        high_pct = sample_high_pct
      ))
      
    }, error = function(e) {
      return(list(
        sample = sample_id,
        index = i,
        success = FALSE,
        error = e$message
      ))
    })
    
    return(result)
    
  }, future.seed = TRUE, future.chunk.size = 1)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  plan(sequential)
  
  # ========================================
  # 5. 统计和输出结果
  # ========================================
  success_count <- sum(sapply(results, function(x) x$success))
  error_count <- length(results) - success_count
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              success_count, 
              length(samples_to_plot),
              100 * success_count / length(samples_to_plot)))
  
  if (error_count > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", error_count, length(samples_to_plot)))
    
    cat("失败的样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  [%d] %s: %s\n", res$index, res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  # 成功样本详情
  if (success_count > 0) {
    cat("成功绘制的样本:\n")
    for (res in results) {
      if (res$success) {
        cat(sprintf("  [%d] %-30s | %5d spots | %4d high (%.1f%%) | %.2f MB\n",
                    res$index,
                    res$sample,
                    res$n_spots,
                    res$n_high,
                    res$high_pct,
                    res$file_size_mb))
      }
    }
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(samples_to_plot)))
  
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$isoheight))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 6. 返回统计信息
  # ========================================
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(samples_to_plot),
    output_dir = CONFIG$dirs$isoheight,
    high_expr_total = high_count,
    high_expr_pct = high_pct,
    elapsed_time = as.numeric(elapsed),
    results = results
  ))
}


# ===================================================================
# 辅助函数：%||% 操作符（如果左侧为NULL则返回右侧）
# ===================================================================
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}