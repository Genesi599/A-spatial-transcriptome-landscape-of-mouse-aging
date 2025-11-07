# 07_plot_spatial.R (多线程并行版)

# ===================================================================
# 函数：绘制空间梯度图（正方形平铺版 + 多线程并行）
# 作者：Assistant (优化版)
# 日期：2025-11-05
# 更新：2025-11-06 - 添加多线程并行支持
# ===================================================================

library(future)
library(future.apply)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tibble)
library(RANN)

#' 绘制 Clock Gene Score 和 Distance 的空间梯度图（正方形平铺 + 并行版）
#'
#' @param seurat_obj Seurat 对象，必须包含以下列：
#'   - ClockGene_Score1: Clock基因评分
#'   - ClockGene_Distance: 到高表达区域的距离
#'   - orig.ident: 样本ID
#' @param samples_to_plot 字符向量，要绘制的样本ID列表
#' @param CONFIG 配置列表，必须包含：
#'   - dirs$spatial: 输出目录路径
#'   - plot$expand_margin: 边界扩展比例
#'   - plot$dpi: 图形分辨率
#'   - n_workers: 并行线程数
#'
#' @return 返回结果列表（隐式），包含成功/失败统计
#'
#' @examples
#' plot_spatial_gradient(seurat_obj, samples_to_plot, CONFIG)
#'
plot_spatial_gradient <- function(seurat_obj, samples_to_plot, CONFIG) {
  
  # ========================================
  # 1. 参数验证
  # ========================================
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   空间梯度图绘制（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
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
  expand_margin <- CONFIG$plot$expand_margin %||% 0.05
  dpi <- CONFIG$plot$dpi %||% 300
  n_workers <- CONFIG$n_workers %||% 4
  
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
      
      n_spots <- ncol(seurat_subset)
      
      # --------------------------------
      # 4.2 获取坐标
      # --------------------------------
      coords <- tryCatch({
        GetTissueCoordinates(
          seurat_subset,
          cols = c("row", "col"),
          scale = NULL
        )
      }, error = function(e) {
        # 尝试从 @images 直接提取
        if (sample_id %in% names(seurat_subset@images)) {
          coords_df <- seurat_subset@images[[sample_id]]@coordinates
          
          row_col <- intersect(
            colnames(coords_df),
            c("row", "imagerow", "array_row", "tissue_row")
          )[1]
          col_col <- intersect(
            colnames(coords_df),
            c("col", "imagecol", "array_col", "tissue_col")
          )[1]
          
          if (!is.na(row_col) && !is.na(col_col)) {
            data.frame(
              row = coords_df[[row_col]],
              col = coords_df[[col_col]],
              row.names = rownames(coords_df)
            )
          } else {
            stop("No valid coordinate columns")
          }
        } else {
          stop("No spatial coordinates available")
        }
      })
      
      # 检查坐标列
      if (!all(c("row", "col") %in% colnames(coords))) {
        return(list(
          sample = sample_id,
          index = i,
          success = FALSE,
          error = sprintf("Missing coordinate columns: %s", 
                         paste(colnames(coords), collapse = ", "))
        ))
      }
      
      # --------------------------------
      # 4.3 合并数据
      # --------------------------------
      plot_data <- seurat_subset@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
      
      # 检查缺失值
      na_coords <- sum(is.na(plot_data$col) | is.na(plot_data$row))
      if (na_coords > 0) {
        plot_data <- plot_data %>% filter(!is.na(col), !is.na(row))
      }
      
      if (nrow(plot_data) == 0) {
        return(list(
          sample = sample_id,
          index = i,
          success = FALSE,
          error = "No valid coordinates after filtering"
        ))
      }
      
      # --------------------------------
      # 4.4 自动计算正方形大小
      # --------------------------------
      if (nrow(plot_data) > 10000) {
        sample_idx <- sample(nrow(plot_data), 10000)
        coords_sample <- plot_data[sample_idx, c("col", "row")]
      } else {
        coords_sample <- plot_data[, c("col", "row")]
      }
      
      nn_dist <- RANN::nn2(coords_sample, k = 2)$nn.dists[, 2]
      median_dist <- median(nn_dist, na.rm = TRUE)
      square_size <- median_dist * 1.0
      
      # --------------------------------
      # 4.5 计算坐标范围
      # --------------------------------
      col_range <- range(plot_data$col, na.rm = TRUE)
      row_range <- range(plot_data$row, na.rm = TRUE)
      
      col_limits <- col_range
      row_limits <- row_range
      
      # --------------------------------
      # 4.6 绘制左图：Clock Gene Score
      # --------------------------------
      p_score <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_tile(
          aes(fill = ClockGene_Score1),
          width = square_size,
          height = square_size,
          color = NA
        ) +
        scale_fill_gradientn(
          colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
          name = "Clock Gene\nScore",
          na.value = "gray90"
        ) +
        scale_x_continuous(
          limits = col_limits,
          expand = c(0, 0)
        ) +
        scale_y_reverse(
          limits = rev(row_limits),
          expand = c(0, 0)
        ) +
        coord_fixed(
          ratio = 1,
          xlim = col_limits,
          ylim = rev(row_limits),
          clip = "on"
        ) +
        ggtitle("Clock Gene Score") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8),
          plot.margin = margin(10, 10, 10, 10)
        )
      
      # --------------------------------
      # 4.7 绘制右图：Distance
      # --------------------------------
      p_distance <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_tile(
          aes(fill = ClockGene_Distance),
          width = square_size,
          height = square_size,
          color = NA
        ) +
        scale_fill_gradientn(
          colors = rev(c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027")),
          name = "Distance to\nHigh Score\nRegion",
          na.value = "gray90"
        ) +
        scale_x_continuous(
          limits = col_limits,
          expand = c(0, 0)
        ) +
        scale_y_reverse(
          limits = rev(row_limits),
          expand = c(0, 0)
        ) +
        coord_fixed(
          ratio = 1,
          xlim = col_limits,
          ylim = rev(row_limits),
          clip = "on"
        ) +
        ggtitle("Distance to High Score Region") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8),
          plot.margin = margin(10, 10, 10, 10)
        )
      
      # --------------------------------
      # 4.8 合并图形
      # --------------------------------
      p_combined <- (p_score | p_distance) +
        plot_annotation(
          title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.margin = margin(10, 10, 10, 10)
          )
        )
      
      # --------------------------------
      # 4.9 保存图形
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
      
      file_size_mb <- file.size(output_path) / 1024^2
      
      # 统计信息
      score_stats <- list(
        min = min(plot_data$ClockGene_Score1, na.rm = TRUE),
        max = max(plot_data$ClockGene_Score1, na.rm = TRUE),
        mean = mean(plot_data$ClockGene_Score1, na.rm = TRUE)
      )
      
      dist_stats <- list(
        min = min(plot_data$ClockGene_Distance, na.rm = TRUE),
        max = max(plot_data$ClockGene_Distance, na.rm = TRUE),
        mean = mean(plot_data$ClockGene_Distance, na.rm = TRUE)
      )
      
      return(list(
        sample = sample_id,
        index = i,
        success = TRUE,
        file = output_path,
        file_size_mb = file_size_mb,
        n_spots = n_spots,
        n_valid_coords = nrow(plot_data),
        square_size = square_size,
        score_range = sprintf("[%.3f, %.3f]", score_stats$min, score_stats$max),
        dist_range = sprintf("[%.1f, %.1f]", dist_stats$min, dist_stats$max)
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
    cat(sprintf("%-4s %-30s | %6s | %6s | %8s | %20s | %20s\n",
                "No.", "Sample", "Spots", "Valid", "Size(MB)", 
                "Score Range", "Dist Range"))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("[%2d] %-30s | %6d | %6d | %8.2f | %20s | %20s\n",
                    res$index,
                    res$sample,
                    res$n_spots,
                    res$n_valid_coords,
                    res$file_size_mb,
                    res$score_range,
                    res$dist_range))
      }
    }
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(samples_to_plot)))
  
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  cat("📐 使用正方形平铺 (geom_tile)\n")
  cat("🔄 Y 轴已反转以匹配 Isoheight 图\n")
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 6. 返回统计信息
  # ========================================
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(samples_to_plot),
    output_dir = CONFIG$dirs$spatial,
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