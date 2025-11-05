# ===================================================================
# 函数：绘制空间梯度图
# 作者：Assistant
# 日期：2025-11-05
# ===================================================================

#' 绘制 Clock Gene Score 和 Distance 的空间梯度图
#'
#' @param seurat_obj Seurat 对象，必须包含以下列：
#'   - ClockGene_Score1: Clock基因评分
#'   - ClockGene_Distance: 到高表达区域的距离
#'   - orig.ident: 样本ID
#' @param samples_to_plot 字符向量，要绘制的样本ID列表
#' @param CONFIG 配置列表，必须包含：
#'   - dirs$spatial: 输出目录路径
#'   - plot$point_size_scatter: 点的大小
#'   - plot$expand_margin: 边界扩展比例
#'   - plot$dpi: 图形分辨率
#'
#' @return 无返回值，直接保存PDF文件到 CONFIG$dirs$spatial
#'
#' @examples
#' plot_spatial_gradient(seurat_obj, samples_to_plot, CONFIG)
#'
plot_spatial_gradient <- function(seurat_obj, samples_to_plot, CONFIG) {
  
  # ========================================
  # 1. 参数验证
  # ========================================
  cat("\n🔥 绘制空间梯度图（匹配 Isoheight 坐标）...\n")
  
  # 检查必需的列
  required_cols <- c("ClockGene_Score1", "ClockGene_Distance", "orig.ident")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 检查必需的配置
  if (is.null(CONFIG$dirs$spatial)) {
    stop("❌ CONFIG$dirs$spatial 未定义")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$spatial)) {
    dir.create(CONFIG$dirs$spatial, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("✅ 创建输出目录: %s\n", CONFIG$dirs$spatial))
  }
  
  # 提取配置参数（设置默认值）
  point_size <- CONFIG$plot$point_size_scatter %||% 2.5
  expand_margin <- CONFIG$plot$expand_margin %||% 0.05
  dpi <- CONFIG$plot$dpi %||% 300
  
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
      
      cat(sprintf("   📊 样本包含 %d 个spots\n", ncol(seurat_subset)))
      
      # --------------------------------
      # 3.2 获取坐标
      # --------------------------------
      coords <- GetTissueCoordinates(
        seurat_subset,
        cols = c("row", "col"),
        scale = NULL
      )
      
      # 检查坐标列名
      coord_cols <- colnames(coords)
      
      if (!all(c("row", "col") %in% coord_cols)) {
        stop(sprintf("坐标列不完整，可用列: %s", paste(coord_cols, collapse = ", ")))
      }
      
      cat(sprintf("   ✅ 坐标列: %s\n", paste(coord_cols, collapse = ", ")))
      
      # --------------------------------
      # 3.3 合并数据
      # --------------------------------
      plot_data <- seurat_subset@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
      
      # 检查缺失值
      na_coords <- sum(is.na(plot_data$col) | is.na(plot_data$row))
      if (na_coords > 0) {
        warning(sprintf("   ⚠️ %d 个spots缺少坐标，已过滤", na_coords))
        plot_data <- plot_data %>% filter(!is.na(col), !is.na(row))
      }
      
      # --------------------------------
      # 3.4 计算坐标范围
      # --------------------------------
      col_range <- range(plot_data$col, na.rm = TRUE)
      row_range <- range(plot_data$row, na.rm = TRUE)
      
      col_expand <- diff(col_range) * expand_margin
      row_expand <- diff(row_range) * expand_margin
      
      col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
      row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
      
      cat(sprintf("   📐 坐标范围: col[%.1f, %.1f], row[%.1f, %.1f]\n",
                  col_limits[1], col_limits[2], row_limits[1], row_limits[2]))
      
      # --------------------------------
      # 3.5 绘制左图：Clock Gene Score
      # --------------------------------
      p_score <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_point(
          aes(fill = ClockGene_Score1), 
          shape = 21, 
          size = point_size, 
          color = "white", 
          stroke = 0.1
        ) +
        scale_fill_gradientn(
          colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
          name = "Clock Gene\nScore",
          na.value = "gray90"
        ) +
        scale_x_continuous(
          limits = col_limits,
          expand = expansion(mult = 0.02)
        ) +
        scale_y_reverse(
          limits = rev(row_limits),
          expand = expansion(mult = 0.02)
        ) +
        coord_fixed(ratio = 1) +
        ggtitle("Clock Gene Score") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          aspect.ratio = 1,
          plot.margin = margin(5, 5, 5, 5)
        )
      
      # --------------------------------
      # 3.6 绘制右图：Distance
      # --------------------------------
      p_distance <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_point(
          aes(fill = ClockGene_Distance), 
          shape = 21, 
          size = point_size, 
          color = "white", 
          stroke = 0.1
        ) +
        scale_fill_gradientn(
          colors = rev(c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027")),
          name = "Distance to\nHigh Score\nRegion",
          na.value = "gray90"
        ) +
        scale_x_continuous(
          limits = col_limits,
          expand = expansion(mult = 0.02)
        ) +
        scale_y_reverse(
          limits = rev(row_limits),
          expand = expansion(mult = 0.02)
        ) +
        coord_fixed(ratio = 1) +
        ggtitle("Distance to High Score Region") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          aspect.ratio = 1,
          plot.margin = margin(5, 5, 5, 5)
        )
      
      # --------------------------------
      # 3.7 合并图形
      # --------------------------------
      p_combined <- (p_score | p_distance) +
        plot_annotation(
          title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )
      
      # --------------------------------
      # 3.8 保存图形
      # --------------------------------
      safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
      output_path <- file.path(
        CONFIG$dirs$spatial, 
        sprintf("ClockGene_spatial_%s.pdf", safe_name)
      )
      
      ggsave(
        filename = output_path,
        plot = p_combined, 
        width = 16, 
        height = 8, 
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
  cat("✅ 空间梯度图绘制完成！\n")
  cat(sprintf("   成功: %d/%d\n", success_count, length(samples_to_plot)))
  if (error_count > 0) {
    cat(sprintf("   失败: %d/%d\n", error_count, length(samples_to_plot)))
  }
  cat(sprintf("   输出目录: %s\n", CONFIG$dirs$spatial))
  cat("   Y 轴已反转以匹配 Isoheight 图\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # 返回统计信息（隐式）
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(samples_to_plot),
    output_dir = CONFIG$dirs$spatial
  ))
}


# ===================================================================
# 辅助函数：%||% 操作符（如果左侧为NULL则返回右侧）
# ===================================================================
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}