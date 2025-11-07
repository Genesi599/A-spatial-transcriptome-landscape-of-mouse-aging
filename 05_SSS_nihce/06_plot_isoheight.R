# 06_plot_isoheight.R (修复内存问题版)

library(future)
library(future.apply)
library(progressr)

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
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   等高线图绘制（多线程并行 - 优化内存版）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证（同之前）
  # ========================================
  required_cols <- c("ClockGene_High", "orig.ident")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.logical(seurat_obj$ClockGene_High)) {
    stop("❌ ClockGene_High 必须是逻辑值（TRUE/FALSE）")
  }
  
  if (is.null(CONFIG$dirs$isoheight)) {
    stop("❌ CONFIG$dirs$isoheight 未定义")
  }
  
  if (!exists("celltype_isoheight_plot")) {
    stop("❌ 未找到 celltype_isoheight_plot 函数，请先加载 SSS_isoheight_plot.R")
  }
  
  if (!dir.exists(CONFIG$dirs$isoheight)) {
    dir.create(CONFIG$dirs$isoheight, recursive = TRUE, showWarnings = FALSE)
  }
  
  size_bg <- CONFIG$plot$point_size_bg %||% 0.3
  size_top <- CONFIG$plot$point_size_top %||% 1.2
  dpi <- CONFIG$plot$dpi %||% 300
  n_workers <- CONFIG$n_workers %||% 4
  
  if (is.null(cols_fill_isoheight)) {
    cols_fill_isoheight <- c(
      rep("white", 25),
      colorRampPalette(brewer.pal(9, "YlOrRd")[3:9])(75)
    )
  }
  
  # 样本验证
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
  
  high_count <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
  high_pct <- 100 * high_count / ncol(seurat_obj)
  cat(sprintf("📊 高表达点: %d / %d (%.2f%%)\n", 
              high_count, ncol(seurat_obj), high_pct))
  
  # ========================================
  # 关键改进：预先切分样本，减少内存传输
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
  
  # 计算单个样本平均大小
  avg_size_mb <- object.size(sample_list[[1]]) / 1024^2
  total_size_mb <- avg_size_mb * length(sample_list)
  cat(sprintf("💾 单样本大小: %.2f MB, 总计: %.2f MB\n", avg_size_mb, total_size_mb))
  
  # 动态调整线程数（避免内存溢出）
  max_memory_gb <- 100  # 假设最大可用内存 100GB
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  n_workers <- min(n_workers, safe_workers, length(sample_list))
  
  cat(sprintf("🔧 使用 %d 个线程 (根据内存自动调整)\n\n", n_workers))
  
  # ========================================
  # 2. 设置并行
  # ========================================
  plan(multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 3. 并行处理（传递小对象）
  # ========================================
  
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  with_progress({
    p <- progressor(steps = length(sample_list))
    
    # 关键：只传递样本子集，不传递整个 seurat_obj
    results <- future_lapply(names(sample_list), function(sample_id) {
      
      result <- tryCatch({
        
        # 从预切分的列表中获取
        seurat_subset <- sample_list[[sample_id]]
        
        if (ncol(seurat_subset) == 0) {
          p(message = sprintf("❌ %s - 无数据", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "No data for this sample"
          ))
        }
        
        sample_high_count <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
        sample_high_pct <- 100 * sample_high_count / ncol(seurat_subset)
        
        if (sample_high_count == 0) {
          p(message = sprintf("⚠️  %s - 无高表达点", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "No high expression spots"
          ))
        }
        
        # 调用绘图函数
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
        
        # 保存图形
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
        
        # 更新进度条
        p(message = sprintf("✅ %s (%.2f MB)", sample_id, file_size_mb))
        
        return(list(
          sample = sample_id,
          success = TRUE,
          file = output_path,
          file_size_mb = file_size_mb,
          n_spots = ncol(seurat_subset),
          n_high = sample_high_count,
          high_pct = sample_high_pct
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
  # 4. 统计和输出（同之前）
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
    for (res in results) {
      if (res$success) {
        cat(sprintf("  %-30s | %5d spots | %4d high (%.1f%%) | %.2f MB\n",
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
              as.numeric(elapsed) / length(sample_list)))
  
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$isoheight))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(sample_list),
    output_dir = CONFIG$dirs$isoheight,
    high_expr_total = high_count,
    high_expr_pct = high_pct,
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