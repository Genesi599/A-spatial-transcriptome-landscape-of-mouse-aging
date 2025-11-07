#!/usr/bin/env Rscript
# ==================================================================
# 空间梯度图绘制模块（串联版 - 正方形平铺）
# 功能：绘制 Clock Gene Score 和 Distance 的空间梯度图
# ==================================================================

#' 绘制空间梯度图（接收预切分样本，正方形平铺）
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param plot_width 图宽
#' @param plot_height 图高
#' 
#' @return 处理结果列表
#'
plot_spatial_gradient <- function(sample_list,
                                  CONFIG,
                                  plot_width = 16,
                                  plot_height = 8) {
  
  cat("\n")
  cat(paste0(
    "═══════════════════════════════════════════════",
    "══════════════\n"
  ))
  cat("   空间梯度图绘制（正方形平铺）\n")
  cat(paste0(
    "═══════════════════════════════════════════════",
    "══════════════\n\n"
  ))
  
  # ======================================
  # 1. 参数验证
  # ======================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  if (is.null(CONFIG$dirs$spatial)) {
    stop("❌ CONFIG$dirs$spatial 未定义")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$spatial)) {
    dir.create(
      CONFIG$dirs$spatial, 
      recursive = TRUE, 
      showWarnings = FALSE
    )
  }
  
  # 提取参数
  `%||%` <- function(a, b) if (is.null(a)) b else a
  expand_margin <- CONFIG$plot$expand_margin %||% 0.05
  dpi <- CONFIG$plot$dpi %||% 300
  
  cat(sprintf("📊 将绘制 %d 个样本\n\n", length(sample_list)))
  
  start_time <- Sys.time()
  
  # ======================================
  # 2. 串联绘图
  # ======================================
  
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
      
      # ----------------------------
      # 验证数据
      # ----------------------------
      if (ncol(seurat_subset) == 0) {
        cat(sprintf("⚠️  %s - 无数据\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "No data"
        )
        next
      }
      
      required_cols <- c("ClockGene_Score1", "ClockGene_Distance")
      missing_cols <- setdiff(
        required_cols, 
        colnames(seurat_subset@meta.data)
      )
      
      if (length(missing_cols) > 0) {
        cat(sprintf(
          "⚠️  %s - 缺少列: %s\n", 
          sample_id, 
          paste(missing_cols, collapse = ", ")
        ))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = sprintf(
            "Missing columns: %s", 
            paste(missing_cols, collapse = ", ")
          )
        )
        next
      }
      
      # ----------------------------
      # 获取坐标
      # ----------------------------
      coords <- Seurat::GetTissueCoordinates(
        seurat_subset,
        cols = c("row", "col"),
        scale = NULL
      )
      
      if (!all(c("row", "col") %in% colnames(coords))) {
        cat(sprintf("⚠️  %s - 坐标列不完整\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "Incomplete coordinate columns"
        )
        next
      }
      
      # ----------------------------
      # ✅ 修复：合并数据（添加命名空间前缀）
      # ----------------------------
      plot_data <- seurat_subset@meta.data %>%
        tibble::rownames_to_column("barcode") %>%
        dplyr::left_join(
          coords %>% tibble::rownames_to_column("barcode"), 
          by = "barcode"
        ) %>%
        dplyr::filter(!is.na(col), !is.na(row))
      
      if (nrow(plot_data) == 0) {
        cat(sprintf("⚠️  %s - 无有效坐标\n", sample_id))
        failed_list[[sample_id]] <- list(
          sample = sample_id,
          success = FALSE,
          error = "No valid coordinates"
        )
        next
      }
      
      # ----------------------------
      # 自动计算正方形大小
      # ----------------------------
      if (nrow(plot_data) > 10000) {
        sample_idx <- sample(nrow(plot_data), 10000)
        coords_sample <- plot_data[sample_idx, c("col", "row")]
      } else {
        coords_sample <- plot_data[, c("col", "row")]
      }
      
      nn_dist <- RANN::nn2(coords_sample, k = 2)$nn.dists[, 2]
      median_dist <- median(nn_dist, na.rm = TRUE)
      square_size <- median_dist * 1.0
      
      # ----------------------------
      # 计算坐标范围
      # ----------------------------
      col_range <- range(plot_data$col, na.rm = TRUE)
      row_range <- range(plot_data$row, na.rm = TRUE)
      
      col_limits <- col_range
      row_limits <- row_range
      
      # ----------------------------
      # 统计数据
      # ----------------------------
      score_stats <- list(
        min = min(plot_data$ClockGene_Score1, na.rm = TRUE),
        max = max(plot_data$ClockGene_Score1, na.rm = TRUE),
        mean = mean(plot_data$ClockGene_Score1, na.rm = TRUE),
        median = median(plot_data$ClockGene_Score1, na.rm = TRUE)
      )
      
      distance_stats <- list(
        min = min(plot_data$ClockGene_Distance, na.rm = TRUE),
        max = max(plot_data$ClockGene_Distance, na.rm = TRUE),
        mean = mean(plot_data$ClockGene_Distance, na.rm = TRUE),
        median = median(plot_data$ClockGene_Distance, na.rm = TRUE)
      )
      
      # ----------------------------
      # 绘制左图：Clock Gene Score
      # ----------------------------
      p_score <- ggplot2::ggplot(
        plot_data, 
        ggplot2::aes(x = col, y = row)
      ) +
        ggplot2::geom_tile(
          ggplot2::aes(fill = ClockGene_Score1),
          width = square_size,
          height = square_size,
          color = NA
        ) +
        ggplot2::scale_fill_gradientn(
          colors = c(
            "#313695", "#4575b4", "#abd9e9", 
            "#fee090", "#f46d43", "#d73027"
          ),
          name = "Clock Gene\nScore",
          na.value = "gray90"
        ) +
        ggplot2::scale_x_continuous(
          limits = col_limits,
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_reverse(
          limits = rev(row_limits),
          expand = c(0, 0)
        ) +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = col_limits,
          ylim = rev(row_limits),
          clip = "on"
        ) +
        ggplot2::ggtitle(
          "Clock Gene Score",
          subtitle = sprintf(
            "Mean: %.3f | Median: %.3f", 
            score_stats$mean, 
            score_stats$median
          )
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5, size = 14, face = "bold"
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5, size = 10
          ),
          legend.position = "right",
          legend.title = ggplot2::element_text(
            size = 10, face = "bold"
          ),
          legend.text = ggplot2::element_text(size = 8),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
      
      # ----------------------------
      # 绘制右图：Distance
      # ----------------------------
      p_distance <- ggplot2::ggplot(
        plot_data, 
        ggplot2::aes(x = col, y = row)
      ) +
        ggplot2::geom_tile(
          ggplot2::aes(fill = ClockGene_Distance),
          width = square_size,
          height = square_size,
          color = NA
        ) +
        ggplot2::scale_fill_gradientn(
          colors = rev(c(
            "#313695", "#4575b4", "#abd9e9", 
            "#fee090", "#f46d43", "#d73027"
          )),
          name = "Distance to\nHigh Score\nRegion",
          na.value = "gray90"
        ) +
        ggplot2::scale_x_continuous(
          limits = col_limits,
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_reverse(
          limits = rev(row_limits),
          expand = c(0, 0)
        ) +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = col_limits,
          ylim = rev(row_limits),
          clip = "on"
        ) +
        ggplot2::ggtitle(
          "Distance to High Score Region",
          subtitle = sprintf(
            "Mean: %.2f | Median: %.2f", 
            distance_stats$mean, 
            distance_stats$median
          )
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5, size = 14, face = "bold"
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5, size = 10
          ),
          legend.position = "right",
          legend.title = ggplot2::element_text(
            size = 10, face = "bold"
          ),
          legend.text = ggplot2::element_text(size = 8),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
      
      # ----------------------------
      # 合并图形
      # ----------------------------
      p_combined <- (p_score | p_distance) +
        patchwork::plot_annotation(
          title = sprintf(
            "Clock Gene Niche Analysis - %s", 
            sample_id
          ),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5, size = 16, face = "bold"
            ),
            plot.margin = ggplot2::margin(10, 10, 10, 10)
          )
        )
      
      # ----------------------------
      # 保存图形
      # ----------------------------
      safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
      output_file <- sprintf("ClockGene_spatial_%s.pdf", safe_name)
      output_path <- file.path(CONFIG$dirs$spatial, output_file)
      
      ggplot2::ggsave(
        filename = output_path,
        plot = p_combined,
        width = plot_width,
        height = plot_height,
        dpi = dpi
      )
      
      # 统计信息
      file_size_mb <- file.size(output_path) / 1024^2
      n_spots <- nrow(plot_data)
      
      # 输出成功信息
      cat(sprintf(
        paste0(
          "✅ %s (%d spots, score: %.3f±%.3f, ",
          "dist: %.2f±%.2f, %.2f MB)\n"
        ), 
        sample_id, n_spots, 
        score_stats$mean, 
        sd(plot_data$ClockGene_Score1, na.rm = TRUE),
        distance_stats$mean, 
        sd(plot_data$ClockGene_Distance, na.rm = TRUE),
        file_size_mb
      ))
      
      success_list[[sample_id]] <- list(
        sample = sample_id,
        success = TRUE,
        file = output_path,
        file_size_mb = file_size_mb,
        n_spots = n_spots,
        score_stats = score_stats,
        distance_stats = distance_stats
      )
      
      # 清理内存
      rm(seurat_subset, plot_data, p_score, p_distance, p_combined)
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
  
  # ======================================
  # 3. 统计输出
  # ======================================
  
  n_success <- length(success_list)
  n_failed <- length(failed_list)
  
  cat("\n")
  cat(paste0(
    "═══════════════════════════════════════════════",
    "══════════════\n"
  ))
  cat("   绘图完成\n")
  cat(paste0(
    "═══════════════════════════════════════════════",
    "══════════════\n\n"
  ))
  
  cat(sprintf(
    "✅ 成功: %d/%d (%.1f%%)\n", 
    n_success, 
    total_samples,
    100 * n_success / total_samples
  ))
  
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
    cat(sprintf(
      "%-30s %10s %12s %12s %10s\n", 
      "样本", "Spots", "Mean Score", "Mean Dist", "文件大小"
    ))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    total_file_size <- 0
    
    for (sample_id in names(success_list)) {
      res <- success_list[[sample_id]]
      cat(sprintf(
        "%-30s %10d %12.3f %12.2f %8.2f MB\n",
        res$sample,
        res$n_spots,
        res$score_stats$mean,
        res$distance_stats$mean,
        res$file_size_mb
      ))
      total_file_size <- total_file_size + res$file_size_mb
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf(
      "%-30s %10s %12s %12s %8.2f MB\n",
      "总计", "", "", "", total_file_size
    ))
    cat("\n")
  }
  
  cat(sprintf(
    "⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
    as.numeric(elapsed),
    as.numeric(elapsed) / total_samples
  ))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  cat("📐 使用正方形平铺 (geom_tile)\n")
  cat("🔄 Y轴已反转以匹配 Isoheight 图\n")
  cat(paste0(
    "\n═══════════════════════════════════════════════",
    "══════════════\n\n"
  ))
  
  # ======================================
  # 4. 返回结果
  # ======================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = total_samples,
    output_dir = CONFIG$dirs$spatial,
    elapsed_time = as.numeric(elapsed),
    results = results
  )))
}

cat("✅ 07_plot_spatial.R 已加载\n")