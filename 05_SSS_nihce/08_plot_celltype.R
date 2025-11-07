#!/usr/bin/env Rscript
# ===================================================================
# 细胞类型 Niche 分析模块（优化版）
# 功能：分析不同密度区域的细胞类型分布和富集
# ===================================================================

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
# 主函数：细胞类型等高线分析（接收预切分样本）
# ===================================================================

#' 细胞类型 Niche 分析
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param density_bins 密度分区数量
#' @param celltype_col 细胞类型列名
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
#' @param seurat_basename 文件基础名
#' 
#' @return 处理结果列表
#'
analyze_celltype_niche <- function(
    sample_list,
    CONFIG,
    density_bins = 10,
    celltype_col = "predicted.id",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE,
    seurat_basename = NULL
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   细胞类型 Niche 分析（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  # 验证必需目录
  required_dirs <- c("overlay", "celltype", "composition", "heatmaps", "combined")
  
  for (dir_name in required_dirs) {
    if (is.null(CONFIG$dirs[[dir_name]])) {
      stop(sprintf("❌ CONFIG$dirs$%s 未定义", dir_name))
    }
    if (!dir.exists(CONFIG$dirs[[dir_name]])) {
      dir.create(CONFIG$dirs[[dir_name]], recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 检查必需函数
  required_functions <- c(
    "calculate_density_zones",
    "plot_celltype_density_overlay",
    "plot_zone_composition",
    "plot_combined_heatmap",
    "plot_combined_analysis",
    "generate_summary_statistics"
  )
  
  missing_funcs <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_funcs) > 0) {
    stop(sprintf("❌ 缺少必需函数: %s", paste(missing_funcs, collapse = ", ")))
  }
  
  # ========================================
  # 2. 初始化颜色配置
  # ========================================
  
  # 从第一个样本获取所有细胞类型
  first_sample <- sample_list[[1]]
  all_celltypes <- sort(unique(as.character(first_sample[[celltype_col]][,1])))
  
  if (is.null(CONFIG$colors$celltype_colors)) {
    CONFIG$colors$celltype_colors <- get_celltype_colors(all_celltypes)
    cat(sprintf("🎨 已生成 %d 种细胞类型颜色方案\n", length(CONFIG$colors$celltype_colors)))
  }
  
  if (is.null(CONFIG$colors$zone_colors)) {
    CONFIG$colors$zone_colors <- get_zone_colors(density_bins)
  }
  
  # ========================================
  # 3. 准备并行环境
  # ========================================
  
  n_workers <- CONFIG$n_workers %||% 4
  
  cat(sprintf("📊 将分析 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 密度分区: %d 个区域 (Zone_0=核心, Zone_%d=外围)\n", 
              density_bins, density_bins - 1))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 4. 并行处理样本
  # ========================================
  
  progressr::with_progress({
    
    p <- progressr::progressor(steps = length(sample_list))
    
    results <- future.apply::future_lapply(
      
      names(sample_list),
      
      function(sample_id) {
        
        tryCatch({
          
          # -------------------------------
          # 4.1 准备数据
          # -------------------------------
          seurat_subset <- sample_list[[sample_id]]
          
          if (ncol(seurat_subset) == 0) {
            p(message = sprintf("⚠️  %s - 无数据", sample_id))
            return(list(sample = sample_id, success = FALSE, error = "No data"))
          }
          
          # 获取坐标
          coords <- Seurat::GetTissueCoordinates(
            seurat_subset,
            cols = c("row", "col"),
            scale = NULL
          )
          
          # 合并元数据
          df <- seurat_subset@meta.data %>%
            tibble::rownames_to_column("barcode") %>%
            dplyr::left_join(
              coords %>% tibble::rownames_to_column("barcode"), 
              by = "barcode"
            ) %>%
            dplyr::filter(!is.na(col), !is.na(row))
          
          if (nrow(df) == 0) {
            p(message = sprintf("⚠️  %s - 无有效坐标", sample_id))
            return(list(sample = sample_id, success = FALSE, error = "No valid coordinates"))
          }
          
          # 检查必需列
          if (!celltype_col %in% colnames(df)) {
            p(message = sprintf("⚠️  %s - 缺少细胞类型列", sample_id))
            return(list(sample = sample_id, success = FALSE, error = "Missing celltype column"))
          }
          
          if (!"ClockGene_High" %in% colnames(df)) {
            p(message = sprintf("⚠️  %s - 缺少 ClockGene_High 列", sample_id))
            return(list(sample = sample_id, success = FALSE, error = "Missing ClockGene_High column"))
          }
          
          # 清理细胞类型
          df$celltype_clean <- as.character(df[[celltype_col]])
          df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
          
          # 统计基本信息
          n_spots <- nrow(df)
          n_high <- sum(df$ClockGene_High, na.rm = TRUE)
          high_pct <- 100 * mean(df$ClockGene_High, na.rm = TRUE)
          
          # -------------------------------
          # 4.2 计算密度区域
          # -------------------------------
          density_data <- calculate_density_zones(
            df = df,
            density_bins = density_bins,
            expand_margin = CONFIG$plot$expand_margin %||% 0.1
          )
          
          if (is.null(density_data)) {
            p(message = sprintf("⚠️  %s - 密度计算失败", sample_id))
            return(list(sample = sample_id, success = FALSE, error = "Density calculation failed"))
          }
          
          # 合并密度信息
          df <- df %>%
            dplyr::left_join(
              density_data$spot_zones %>% 
                dplyr::select(col, row, density_zone, density_value),
              by = c("col", "row")
            )
          
          n_na_zones <- sum(is.na(df$density_zone))
          
          # -------------------------------
          # 4.3 计算区域组成
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
          # 4.4 绘制图形
          # -------------------------------
          output_files <- list()
          total_size <- 0
          
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          
          # 叠加图
          if (plot_overlay) {
            p_overlay <- plot_celltype_density_overlay(
              df = df,
              density_data = density_data,
              sample_id = sample_id,
              CONFIG = CONFIG
            )
            
            overlay_file <- file.path(
              CONFIG$dirs$overlay, 
              sprintf("celltype_overlay_%s.pdf", safe_name)
            )
            
            ggplot2::ggsave(
              overlay_file,
              plot = p_overlay,
              width = 12, 
              height = 10,
              dpi = CONFIG$plot$dpi %||% 300,
              bg = "white"
            )
            
            output_files$overlay <- overlay_file
            total_size <- total_size + file.size(overlay_file)
          }
          
          # 组成图
          if (plot_composition) {
            p_comp <- plot_zone_composition(
              zone_composition = zone_composition,
              sample_id = sample_id,
              CONFIG = CONFIG
            )
            
            composition_file <- file.path(
              CONFIG$dirs$composition, 
              sprintf("composition_%s.pdf", safe_name)
            )
            
            ggplot2::ggsave(
              composition_file,
              plot = p_comp,
              width = 12, 
              height = 6,
              dpi = CONFIG$plot$dpi %||% 300,
              bg = "white"
            )
            
            output_files$composition <- composition_file
            total_size <- total_size + file.size(composition_file)
          }
          
          total_size_mb <- total_size / 1024^2
          
          # 更新进度
          p(message = sprintf("✅ %s", sample_id))
          
          # -------------------------------
          # 4.5 返回结果
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
            output_files = output_files,
            total_size_mb = total_size_mb
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
  # 5. 统计样本处理结果
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本处理完成\n")
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
    cat(sprintf("%-30s %8s %7s %8s %7s %8s %7s %10s\n",
                "样本", "Spots", "High", "High%", "Zones", "Types", "NA", "大小(MB)"))
    cat(paste(rep("-", 100), collapse = ""), "\n")
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %8d %7d %7.2f%% %7d %8d %7d %10.2f\n",
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
  # 6. 生成综合分析
  # ========================================
  
  combined_data <- data.frame()
  
  for (res in results) {
    if (res$success) {
      combined_data <- dplyr::bind_rows(combined_data, res$zone_composition)
    }
  }
  
  if (nrow(combined_data) > 0 && n_success > 0) {
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("   生成综合统计图\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    combined_start <- Sys.time()
    
    main_title <- seurat_basename %||% "Seurat Object"
    combined_files <- list()
    
    # 热图
    if (plot_heatmap) {
      cat("📊 生成细胞类型热图...\n")
      
      tryCatch({
        p_heatmap <- plot_combined_heatmap(
          combined_data = combined_data, 
          CONFIG = CONFIG
        ) + ggplot2::ggtitle(main_title)
        
        heatmap_file <- file.path(
          CONFIG$dirs$heatmaps, 
          "celltype_heatmap_all_samples.pdf"
        )
        
        ggplot2::ggsave(
          heatmap_file,
          plot = p_heatmap, 
          width = 14, 
          height = 10, 
          dpi = CONFIG$plot$dpi %||% 300, 
          bg = "white"
        )
        
        combined_files$heatmap <- heatmap_file
        
        cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                    basename(heatmap_file),
                    file.size(heatmap_file) / 1024^2))
      }, error = function(e) {
        cat(sprintf("   ⚠️  热图生成失败: %s\n", e$message))
      })
    }
    
    # 综合分析图
    if (plot_combined) {
      cat("📊 生成综合分析图...\n")
      
      tryCatch({
        p_combined <- plot_combined_analysis(
          combined_data = combined_data, 
          CONFIG = CONFIG
        ) + ggplot2::ggtitle(main_title)
        
        combined_file <- file.path(
          CONFIG$dirs$combined, 
          "combined_analysis.pdf"
        )
        
        ggplot2::ggsave(
          combined_file,
          plot = p_combined, 
          width = 16, 
          height = 12, 
          dpi = CONFIG$plot$dpi %||% 300, 
          bg = "white"
        )
        
        combined_files$combined <- combined_file
        
        cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                    basename(combined_file),
                    file.size(combined_file) / 1024^2))
      }, error = function(e) {
        cat(sprintf("   ⚠️  综合图生成失败: %s\n", e$message))
      })
    }
    
    # 保存数据
    cat("💾 保存统计数据...\n")
    
    # 组成数据
    composition_csv <- file.path(
      CONFIG$dirs$composition, 
      "celltype_composition_all_samples.csv"
    )
    write.csv(combined_data, composition_csv, row.names = FALSE)
    cat(sprintf("   ✅ 组成数据: %s\n", basename(composition_csv)))
    
    # 汇总统计
    tryCatch({
      summary_stats <- generate_summary_statistics(combined_data)
      summary_csv <- file.path(
        CONFIG$dirs$composition, 
        "summary_statistics.csv"
      )
      write.csv(summary_stats, summary_csv, row.names = FALSE)
      cat(sprintf("   ✅ 汇总统计: %s\n", basename(summary_csv)))
    }, error = function(e) {
      cat(sprintf("   ⚠️  统计计算失败: %s\n", e$message))
    })
    
    combined_end <- Sys.time()
    combined_elapsed <- difftime(combined_end, combined_start, units = "secs")
    
    cat(sprintf("\n⏱️  综合图生成耗时: %.2f 秒\n", as.numeric(combined_elapsed)))
  }
  
  # ========================================
  # 7. 最终总结
  # ========================================
  
  total_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d\n", n_success, length(sample_list)))
  cat(sprintf("⏱️  总耗时: %.2f 秒 (%.2f 分钟)\n", 
              as.numeric(total_elapsed),
              as.numeric(total_elapsed) / 60))
  
  if (n_success > 0) {
    cat("\n📊 生成内容:\n")
    if (plot_overlay) 
      cat(sprintf("   • 叠加图: %d 个\n", n_success))
    if (plot_composition) 
      cat(sprintf("   • 组成图: %d 个\n", n_success))
    if (plot_heatmap && nrow(combined_data) > 0) 
      cat("   • 热图: 1 个\n")
    if (plot_combined && nrow(combined_data) > 0) 
      cat("   • 综合图: 1 个\n")
  }
  
  cat("\n📁 输出目录:\n")
  cat(sprintf("   • Overlay:     %s\n", CONFIG$dirs$overlay))
  cat(sprintf("   • Composition: %s\n", CONFIG$dirs$composition))
  cat(sprintf("   • Heatmaps:    %s\n", CONFIG$dirs$heatmaps))
  cat(sprintf("   • Combined:    %s\n", CONFIG$dirs$combined))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 8. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    elapsed_time = as.numeric(total_elapsed),
    combined_data = combined_data,
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 08_plot_celltype.R 已加载\n")