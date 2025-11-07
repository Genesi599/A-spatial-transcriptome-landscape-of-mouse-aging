# 08_plot_celltype.R (多线程并行版 + 进度条 + 内存优化)

library(future)
library(future.apply)
library(progressr)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

# ===================================================================
# 加载工具函数
# ===================================================================

utils_dir <- "08_plot_celltype_utils"

source(file.path(utils_dir, "00_operators.R"))
source(file.path(utils_dir, "01_color_schemes.R"))
source(file.path(utils_dir, "02_density_zones.R"))
source(file.path(utils_dir, "03_plot_overlay.R"))
source(file.path(utils_dir, "04_plot_composition.R"))
source(file.path(utils_dir, "05_plot_heatmap.R"))
source(file.path(utils_dir, "06_plot_combined.R"))
source(file.path(utils_dir, "07_statistics.R"))

cat("✅ 已加载所有工具函数\n")


# ===================================================================
# 主函数：细胞类型等高线分析（优化内存版）
# ===================================================================

analyze_celltype_niche <- function(
    seurat_obj,
    samples_to_plot,
    CONFIG,
    density_bins = 10,
    celltype_col = "celltype",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE,
    seurat_basename = NULL
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   细胞类型 + Clock Gene Niche 等高线分析（优化内存版）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 初始化颜色配置
  # ========================================
  all_celltypes <- sort(unique(as.character(seurat_obj[[celltype_col]][,1])))
  
  if (is.null(CONFIG$colors$celltype_colors)) {
    CONFIG$colors$celltype_colors <- get_celltype_colors(all_celltypes)
    cat(sprintf("🎨 已生成 %d 种细胞类型颜色方案\n", length(CONFIG$colors$celltype_colors)))
  }
  
  if (is.null(CONFIG$colors$zone_colors)) {
    CONFIG$colors$zone_colors <- get_zone_colors(density_bins)
  }
  
  # ========================================
  # 2. 参数验证
  # ========================================
  required_cols <- c("ClockGene_High", "orig.ident", celltype_col)
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ Seurat对象缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 创建输出目录
  output_dirs <- list(
    overlay = CONFIG$dirs$overlay,
    celltype = CONFIG$dirs$celltype,
    composition = CONFIG$dirs$composition,
    heatmaps = CONFIG$dirs$heatmaps,
    combined = CONFIG$dirs$combined
  )
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 验证样本
  available_samples <- unique(seurat_obj$orig.ident)
  samples_to_plot <- intersect(samples_to_plot, available_samples)
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本可分析")
  }
  
  cat(sprintf("📊 将分析 %d 个样本\n", length(samples_to_plot)))
  cat(sprintf("📊 等高线分为 %d 个区域 (Zone_0=核心, Zone_%d=外围)\n", 
              density_bins, density_bins - 1))
  
  # ========================================
  # 3. 【关键改进】预先切分样本
  # ========================================
  cat(sprintf("\n🔧 预处理: 切分 %d 个样本...\n", length(samples_to_plot)))
  
  sample_list <- list()
  for (sample_id in samples_to_plot) {
    seurat_subset <- subset(seurat_obj, subset = orig.ident == sample_id)
    
    if (ncol(seurat_subset) > 0) {
      sample_list[[sample_id]] <- seurat_subset
    }
  }
  
  cat(sprintf("✅ 已切分 %d 个样本\n", length(sample_list)))
  
  # 计算内存
  if (length(sample_list) > 0) {
    avg_size_mb <- object.size(sample_list[[1]]) / 1024^2
    total_size_mb <- avg_size_mb * length(sample_list)
    cat(sprintf("💾 单样本大小: %.2f MB, 总计: %.2f MB\n", avg_size_mb, total_size_mb))
    
    # 动态调整线程数
    max_memory_gb <- 100
    safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
    n_workers <- min(CONFIG$n_workers %||% 4, safe_workers, length(sample_list))
  } else {
    n_workers <- CONFIG$n_workers %||% 4
  }
  
  cat(sprintf("🔧 使用 %d 个线程 (根据内存自动调整)\n\n", n_workers))
  
  # ========================================
  # 4. 准备并行环境
  # ========================================
  required_functions <- c(
    "calculate_density_zones",
    "plot_celltype_density_overlay",
    "plot_zone_composition",
    "get_celltype_colors",
    "get_zone_colors"
  )
  
  missing_funcs <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_funcs) > 0) {
    stop(sprintf("❌ 缺少必需函数: %s", paste(missing_funcs, collapse = ", ")))
  }
  
  # 设置并行
  plan(multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 5. 并行处理（只传递 sample_list）
  # ========================================
  
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  with_progress({
    p <- progressor(steps = length(sample_list))
    
    # 【关键】只传递 sample_list
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
        
        # -------------------------------
        # 5.1 获取坐标
        # -------------------------------
        coords <- GetTissueCoordinates(
          seurat_subset,
          cols = c("row", "col"),
          scale = NULL
        )
        
        # 合并数据
        df <- seurat_subset@meta.data %>%
          tibble::rownames_to_column("barcode") %>%
          dplyr::left_join(coords %>% tibble::rownames_to_column("barcode"), by = "barcode") %>%
          dplyr::filter(!is.na(col), !is.na(row))
        
        if (nrow(df) == 0) {
          p(message = sprintf("❌ %s - 无有效坐标", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "No valid coordinates"
          ))
        }
        
        # 检查细胞类型
        df$celltype_clean <- as.character(df[[celltype_col]])
        df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
        
        n_spots <- nrow(df)
        n_high <- sum(df$ClockGene_High)
        high_pct <- 100 * mean(df$ClockGene_High)
        
        # -------------------------------
        # 5.2 计算密度并分级
        # -------------------------------
        density_data <- calculate_density_zones(
          df = df,
          density_bins = density_bins,
          expand_margin = CONFIG$plot$expand_margin %||% 0.1
        )
        
        if (is.null(density_data)) {
          p(message = sprintf("❌ %s - 密度计算失败", sample_id))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = "Density calculation failed"
          ))
        }
        
        df <- df %>%
          dplyr::left_join(
            density_data$spot_zones %>% dplyr::select(col, row, density_zone, density_value),
            by = c("col", "row")
          )
        
        n_na_zones <- sum(is.na(df$density_zone))
        
        # -------------------------------
        # 5.3 组成计算
        # -------------------------------
        zone_composition <- df %>%
          dplyr::filter(!is.na(density_zone)) %>%
          dplyr::group_by(density_zone, celltype_clean) %>%
          dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(density_zone) %>%
          dplyr::mutate(
            total = sum(count),
            percentage = 100 * count / total
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(sample = sample_id)
        
        n_zones <- length(unique(zone_composition$density_zone))
        n_celltypes <- length(unique(zone_composition$celltype_clean))
        
        # -------------------------------
        # 5.4 绘制叠加图
        # -------------------------------
        overlay_file <- NULL
        if (plot_overlay) {
          p_overlay <- plot_celltype_density_overlay(
            df = df,
            density_data = density_data,
            sample_id = sample_id,
            CONFIG = CONFIG
          )
          
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          overlay_file <- file.path(
            CONFIG$dirs$overlay, 
            sprintf("celltype_overlay_%s.pdf", safe_name)
          )
          
          ggsave(
            overlay_file,
            plot = p_overlay,
            width = 12, height = 10,
            dpi = CONFIG$plot$dpi %||% 300,
            bg = "white"
          )
        }
        
        # -------------------------------
        # 5.5 绘制组成图
        # -------------------------------
        composition_file <- NULL
        if (plot_composition) {
          p_comp <- plot_zone_composition(
            zone_composition = zone_composition,
            sample_id = sample_id,
            CONFIG = CONFIG
          )
          
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          composition_file <- file.path(
            CONFIG$dirs$composition, 
            sprintf("composition_%s.pdf", safe_name)
          )
          
          ggsave(
            composition_file,
            plot = p_comp,
            width = 12, height = 6,
            dpi = CONFIG$plot$dpi %||% 300,
            bg = "white"
          )
        }
        
        # 计算文件大小
        total_size <- 0
        if (!is.null(overlay_file) && file.exists(overlay_file)) {
          total_size <- total_size + file.size(overlay_file)
        }
        if (!is.null(composition_file) && file.exists(composition_file)) {
          total_size <- total_size + file.size(composition_file)
        }
        total_size_mb <- total_size / 1024^2
        
        # 更新进度条
        p(message = sprintf("✅ %s - %d zones, %d types (%.2f MB)", 
                           sample_id, n_zones, n_celltypes, total_size_mb))
        
        # -------------------------------
        # 5.6 返回结果
        # -------------------------------
        return(list(
          sample = sample_id,
          success = TRUE,
          zone_composition = zone_composition,
          n_spots = n_spots,
          n_high = n_high,
          high_pct = high_pct,
          n_zones = n_zones,
          n_celltypes = n_celltypes,
          n_na_zones = n_na_zones,
          overlay_file = overlay_file,
          composition_file = composition_file,
          total_size_mb = total_size_mb
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
  # 6. 收集和汇总结果
  # ========================================
  
  success_count <- sum(sapply(results, function(x) x$success))
  error_count <- length(results) - success_count
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本处理完成\n")
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
    cat("成功分析的样本:\n")
    cat(sprintf("%-30s | %6s | %5s | %8s | %6s | %6s | %8s | %10s\n",
                "Sample", "Spots", "High", "High%", "Zones", "Types", "NA_Zones", "Size(MB)"))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s | %6d | %5d | %7.2f%% | %6d | %6d | %8d | %10.2f\n",
                    res$sample,
                    res$n_spots,
                    res$n_high,
                    res$high_pct,
                    res$n_zones,
                    res$n_celltypes,
                    res$n_na_zones,
                    res$total_size_mb))
      }
    }
    cat("\n")
  }
  
  cat(sprintf("⏱️  样本处理耗时: %.2f 秒 (平均 %.2f 秒/样本)\n\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  
  # ========================================
  # 7. 合并数据并生成综合图
  # ========================================
  
  all_sample_stats <- list()
  combined_data <- data.frame()
  
  for (res in results) {
    if (res$success) {
      all_sample_stats[[res$sample]] <- res$zone_composition
      combined_data <- dplyr::bind_rows(combined_data, res$zone_composition)
    }
  }
  
  if (nrow(combined_data) > 0) {
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("   生成综合统计图\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    combined_start_time <- Sys.time()
    
    main_title <- if (!is.null(seurat_basename)) seurat_basename else "Seurat Object"
    
    # 绘制热图
    if (plot_heatmap) {
      cat("📊 生成细胞类型热图...\n")
      
      p_heatmap <- plot_combined_heatmap(
        combined_data = combined_data, 
        CONFIG = CONFIG
      ) + ggtitle(main_title)
      
      heatmap_file <- file.path(
        CONFIG$dirs$heatmaps, 
        "celltype_heatmap_all_samples.pdf"
      )
      
      ggsave(
        heatmap_file,
        plot = p_heatmap, 
        width = 14, 
        height = 10, 
        dpi = CONFIG$plot$dpi %||% 300, 
        bg = "white"
      )
      
      cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                  basename(heatmap_file),
                  file.size(heatmap_file) / 1024^2))
    }
    
    # 绘制综合分析图
    if (plot_combined) {
      cat("📊 生成综合分析图...\n")
      
      p_combined <- plot_combined_analysis(
        combined_data = combined_data, 
        CONFIG = CONFIG
      ) + ggtitle(main_title)
      
      combined_file <- file.path(
        CONFIG$dirs$combined, 
        "combined_analysis.pdf"
      )
      
      ggsave(
        combined_file,
        plot = p_combined, 
        width = 16, 
        height = 12, 
        dpi = CONFIG$plot$dpi %||% 300, 
        bg = "white"
      )
      
      cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                  basename(combined_file),
                  file.size(combined_file) / 1024^2))
    }
    
    # 保存数据
    cat("💾 保存统计数据...\n")
    
    composition_csv <- file.path(
      CONFIG$dirs$composition, 
      "celltype_composition_all_samples.csv"
    )
    write.csv(combined_data, composition_csv, row.names = FALSE)
    cat(sprintf("   ✅ 保存: %s\n", basename(composition_csv)))
    
    summary_stats <- generate_summary_statistics(combined_data)
    summary_csv <- file.path(
      CONFIG$dirs$composition, 
      "summary_statistics.csv"
    )
    write.csv(summary_stats, summary_csv, row.names = FALSE)
    cat(sprintf("   ✅ 保存: %s\n", basename(summary_csv)))
    
    combined_end_time <- Sys.time()
    combined_elapsed <- difftime(combined_end_time, combined_start_time, units = "secs")
    
    cat(sprintf("\n⏱️  综合图生成耗时: %.2f 秒\n", as.numeric(combined_elapsed)))
  }
  
  # ========================================
  # 8. 最终总结
  # ========================================
  
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功分析样本: %d/%d\n", success_count, length(sample_list)))
  cat(sprintf("⏱️  总耗时: %.2f 秒 (%.2f 分钟)\n", 
              as.numeric(total_time),
              as.numeric(total_time) / 60))
  
  if (success_count > 0) {
    cat(sprintf("📊 生成图表:\n"))
    if (plot_overlay) cat(sprintf("   - %d 个叠加图\n", success_count))
    if (plot_composition) cat(sprintf("   - %d 个组成图\n", success_count))
    if (plot_heatmap) cat("   - 1 个综合热图\n")
    if (plot_combined) cat("   - 1 个综合分析图\n")
  }
  
  cat(sprintf("\n📁 输出目录:\n"))
  cat(sprintf("   - 叠加图: %s\n", CONFIG$dirs$overlay))
  cat(sprintf("   - 组成图: %s\n", CONFIG$dirs$composition))
  cat(sprintf("   - 热图: %s\n", CONFIG$dirs$heatmaps))
  cat(sprintf("   - 综合图: %s\n", CONFIG$dirs$combined))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(list(
    success = success_count,
    failed = error_count,
    total = length(sample_list),
    elapsed_time = as.numeric(total_time),
    sample_stats = all_sample_stats,
    combined_data = combined_data,
    results = results
  ))
}