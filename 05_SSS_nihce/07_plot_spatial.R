# 07_plot_spatial.R (多线程并行版 + 进度条 + 内存优化)

library(future)
library(future.apply)
library(progressr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tibble)
library(RANN)

plot_spatial_gradient <- function(seurat_obj, samples_to_plot, CONFIG) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   空间梯度图绘制（多线程并行 - 优化内存版）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  required_cols <- c("ClockGene_Score1", "ClockGene_Distance", "orig.ident")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  if (is.null(CONFIG$dirs$spatial)) {
    stop("❌ CONFIG$dirs$spatial 未定义")
  }
  
  if (!dir.exists(CONFIG$dirs$spatial)) {
    dir.create(CONFIG$dirs$spatial, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("✅ 创建输出目录: %s\n", CONFIG$dirs$spatial))
  }
  
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
  
  # ========================================
  # 3. 【关键改进】预先切分样本
  # ========================================
  cat(sprintf("\n🔧 预处理: 切分 %d 个样本...\n", length(samples_to_plot)))
  
  sample_list <- list()
  for (sample_id in samples_to_plot) {
    seurat_subset <- tryCatch(
      subset(seurat_obj, subset = orig.ident == sample_id),
      error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
    )
    
    if (ncol(seurat_subset) > 0) {
      sample_list[[sample_id]] <- seurat_subset
    }
  }
  
  cat(sprintf("✅ 已切分 %d 个样本\n", length(sample_list)))
  
  # 计算内存使用
  if (length(sample_list) > 0) {
    avg_size_mb <- object.size(sample_list[[1]]) / 1024^2
    total_size_mb <- avg_size_mb * length(sample_list)
    cat(sprintf("💾 单样本大小: %.2f MB, 总计: %.2f MB\n", avg_size_mb, total_size_mb))
    
    # 动态调整线程数
    max_memory_gb <- 100
    safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
    n_workers <- min(n_workers, safe_workers, length(sample_list))
  }
  
  cat(sprintf("🔧 使用 %d 个线程 (根据内存自动调整)\n\n", n_workers))
  
  # ========================================
  # 4. 设置并行
  # ========================================
  plan(multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 5. 并行处理（只传递 sample_list）
  # ========================================
  
  # 【修改】条件设置 handlers
  has_handlers <- !is.null(progressr::handlers(NULL))

  if (!has_handlers) {
    handlers(global = TRUE)
    handlers("txtprogressbar")
  }
  
  with_progress({
    p <- progressor(steps = length(sample_list))
    
    # 【关键】只传递 sample_list，不传递 seurat_obj
    results <- future_lapply(names(sample_list), function(sample_id) {
      
      result <- tryCatch({
        
        # 从预切分列表获取
        seurat_subset <- sample_list[[sample_id]]
        
        if (ncol(seurat_subset) == 0) {
          p(message = sprintf("❌ %s - 无数据", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "No data for this sample"
          ))
        }
        
        n_spots <- ncol(seurat_subset)
        
        # --------------------------------
        # 5.1 获取坐标
        # --------------------------------
        coords <- tryCatch({
          GetTissueCoordinates(
            seurat_subset,
            cols = c("row", "col"),
            scale = NULL
          )
        }, error = function(e) {
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
        
        if (!all(c("row", "col") %in% colnames(coords))) {
          p(message = sprintf("❌ %s - 坐标列不完整", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = sprintf("Missing coordinate columns: %s", 
                           paste(colnames(coords), collapse = ", "))
          ))
        }
        
        # --------------------------------
        # 5.2 合并数据
        # --------------------------------
        plot_data <- seurat_subset@meta.data %>%
          rownames_to_column("barcode") %>%
          left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
        
        na_coords <- sum(is.na(plot_data$col) | is.na(plot_data$row))
        if (na_coords > 0) {
          plot_data <- plot_data %>% filter(!is.na(col), !is.na(row))
        }
        
        if (nrow(plot_data) == 0) {
          p(message = sprintf("❌ %s - 无有效坐标", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "No valid coordinates after filtering"
          ))
        }
        
        # --------------------------------
        # 5.3 计算正方形大小
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
        # 5.4 计算坐标范围
        # --------------------------------
        col_range <- range(plot_data$col, na.rm = TRUE)
        row_range <- range(plot_data$row, na.rm = TRUE)
        
        col_limits <- col_range
        row_limits <- row_range
        
        # --------------------------------
        # 5.5 绘制左图：Score
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
        # 5.6 绘制右图：Distance
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
        # 5.7 合并图形
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
        # 5.8 保存图形
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
        
        # 更新进度条
        p(message = sprintf("✅ %s (%.2f MB)", sample_id, file_size_mb))
        
        return(list(
          sample = sample_id,
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
        p(message = sprintf("❌ %s - %s", sample_id, e$message))
        return(list(
          sample = sample_id,
          success = FALSE,
          error = e$message
        ))
      })
      
      return(result)
      
    }, future.seed = TRUE, future.chunk.size = 1)
  })
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  plan(sequential)
  
  # ========================================
  # 6. 统计和输出
  # ========================================
  success_count <- sum(sapply(results, function(x) x$success))
  error_count <- length(results) - success_count
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              success_count, 
              length(sample_list),
              100 * success_count / length(sample_list)))
  
  if (error_count > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", error_count, length(sample_list)))
    
    cat("失败的样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  %s: %s\n", res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  if (success_count > 0) {
    cat("成功绘制的样本:\n")
    cat(sprintf("%-30s | %6s | %6s | %8s | %20s | %20s\n",
                "Sample", "Spots", "Valid", "Size(MB)", 
                "Score Range", "Dist Range"))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s | %6d | %6d | %8.2f | %20s | %20s\n",
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
              as.numeric(elapsed) / length(sample_list)))
  
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(sample_list),
    output_dir = CONFIG$dirs$spatial,
    elapsed_time = as.numeric(elapsed),
    results = results
  ))
}

# 辅助函数
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}