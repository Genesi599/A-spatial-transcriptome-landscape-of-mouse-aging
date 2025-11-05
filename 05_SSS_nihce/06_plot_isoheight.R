# ===================================================================
# 函数：绘制 Isoheight 图
# 作者：Assistant
# 日期：2025-11-05
# ===================================================================

#' 绘制 Clock Gene High 表达的 Isoheight 密度图
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
#' @param col_bg 背景点颜色，默认 "gray92"
#' @param col_top 高表达点颜色，默认 "#d62728"（红色）
#' @param col_isoheight 等高线颜色，默认 "white"
#' @param col_white_ratio 白色占比，默认 0.25
#' @param cols_fill_isoheight 等高线填充颜色向量，默认为黄-橙-红渐变
#' @param plot_width 图形宽度，默认 8
#' @param plot_height 图形高度，默认 8
#' @param nrow 图形布局行数，默认 1
#'
#' @return 无返回值，直接保存PDF文件到 CONFIG$dirs$isoheight
#'
#' @examples
#' # 基础调用（使用默认颜色）
#' plot_isoheight(seurat_obj, samples_to_plot, CONFIG)
#' 
#' # 自定义颜色
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
  cat("\n🎨 绘制 Isoheight 图...\n")
  
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
  
  cat(sprintf("✅ 将绘制 %d 个样本\n", length(samples_to_plot)))
  
  # 统计高表达点
  high_count <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
  high_pct <- 100 * high_count / ncol(seurat_obj)
  cat(sprintf("✅ 高表达点: %d / %d (%.2f%%)\n", 
              high_count, ncol(seurat_obj), high_pct))
  
  # ========================================
  # 3. 循环绘制每个样本
  # ========================================
  success_count <- 0
  error_count <- 0
  
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    cat(sprintf("\n[%d/%d] 正在处理: %s\n", i, length(samples_to_plot), sample_id))
    
    tryCatch({
      # --------------------------------
      # 3.1 提取子集
      # --------------------------------
      seurat_subset <- tryCatch(
        subset(seurat_obj, subset = orig.ident == sample_id),
        error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
      )
      
      if (ncol(seurat_subset) == 0) {
        stop(sprintf("样本 %s 无数据", sample_id))
      }
      
      # 统计该样本的高表达点
      sample_high_count <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
      sample_high_pct <- 100 * sample_high_count / ncol(seurat_subset)
      
      cat(sprintf("   📊 样本包含 %d 个spots，其中 %d 个高表达点 (%.2f%%)\n", 
                  ncol(seurat_subset), sample_high_count, sample_high_pct))
      
      if (sample_high_count == 0) {
        warning(sprintf("   ⚠️ 样本 %s 没有高表达点，跳过", sample_id))
        next
      }
      
      # --------------------------------
      # 3.2 调用 celltype_isoheight_plot
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
      # 3.3 保存图形
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
      
      cat(sprintf("   ✅ 已保存: %s (%.2f MB)\n", 
                  basename(output_path), 
                  file.size(output_path) / 1024^2))
      
      success_count <- success_count + 1
      
    }, error = function(e) {
      cat(sprintf("   ❌ 错误: %s\n", e$message))
      error_count <- error_count + 1
    })
  }
  
  # ========================================
  # 4. 总结
  # ========================================
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("✅ Isoheight 图绘制完成！\n")
  cat(sprintf("   成功: %d/%d\n", success_count, length(samples_to_plot)))
  if (error_count > 0) {
    cat(sprintf("   失败: %d/%d\n", error_count, length(samples_to_plot)))
  }
  cat(sprintf("   输出目录: %s\n", CONFIG$dirs$isoheight))
  cat(rep("=", 80), "\n\n", sep = "")
  
  # 返回统计信息（隐式）
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(samples_to_plot),
    output_dir = CONFIG$dirs$isoheight,
    high_expr_total = high_count,
    high_expr_pct = high_pct
  ))
}


# ===================================================================
# 辅助函数：%||% 操作符（如果左侧为NULL则返回右侧）
# ===================================================================
# 如果之前没定义过，添加这个
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}