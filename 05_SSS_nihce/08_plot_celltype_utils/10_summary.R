#!/usr/bin/env Rscript
# ===================================================================
# 汇总统计模块
# ===================================================================

#' 打印样本汇总
#' 
#' @param results 结果列表
#' @param sample_list 样本列表
#' @param elapsed 耗时
print_sample_summary <- function(results, sample_list, elapsed) {
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本处理完成\n")
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
    cat(sprintf("%-30s %8s %7s %8s %7s %8s %7s %10s\n",
                "样本", "Spots", "High", "High%", "Zones", "Types", "NA", "大小(MB)"))
    cat(paste(rep("-", 100), collapse = ""), "\n")
    
    total_size <- 0
    total_spots <- 0
    
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
        
        total_size <- total_size + res$total_size_mb
        total_spots <- total_spots + res$n_spots
      }
    }
    
    if (n_success > 1) {
      cat(paste(rep("-", 100), collapse = ""), "\n")
      cat(sprintf("%-30s %8d %7s %8s %7s %8s %7s %10.2f\n",
                  "总计",
                  total_spots,
                  "-", "-", "-", "-", "-",
                  total_size))
    }
    
    cat("\n")
  }
  
  cat(sprintf("⏱️  样本处理耗时: %.2f 秒 (平均 %.2f 秒/样本)\n\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  
  invisible(NULL)
}


#' 收集合并数据
#' 
#' @param results 结果列表
#' @return 合并的数据框
collect_combined_data <- function(results) {
  
  combined_data <- data.frame()
  
  for (res in results) {
    if (res$success) {
      combined_data <- dplyr::bind_rows(combined_data, res$zone_composition)
    }
  }
  
  return(combined_data)
}


#' 生成综合分析
#' 
#' @param combined_data 合并数据
#' @param CONFIG 配置对象
#' @param seurat_basename 基础名
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
generate_combined_analysis <- function(combined_data, CONFIG, seurat_basename,
                                       plot_heatmap, plot_combined) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   生成综合统计图\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  combined_start <- Sys.time()
  
  main_title <- seurat_basename %||% "Seurat Object"
  
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
      
      cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                  basename(combined_file),
                  file.size(combined_file) / 1024^2))
    }, error = function(e) {
      cat(sprintf("   ⚠️  综合图生成失败: %s\n", e$message))
    })
  }
  
  # 保存数据
  cat("💾 保存统计数据...\n")
  
  composition_csv <- file.path(
    CONFIG$dirs$composition, 
    "celltype_composition_all_samples.csv"
  )
  write.csv(combined_data, composition_csv, row.names = FALSE)
  cat(sprintf("   ✅ 组成数据: %s\n", basename(composition_csv)))
  
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
  
  invisible(NULL)
}


#' 打印最终汇总
#' 
#' @param results 结果列表
#' @param sample_list 样本列表
#' @param start_time 开始时间
#' @param combined_data 合并数据
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
#' @param CONFIG 配置对象
print_final_summary <- function(results, sample_list, start_time, combined_data,
                               plot_overlay, plot_composition, plot_heatmap, 
                               plot_combined, CONFIG) {
  
  total_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  n_success <- sum(sapply(results, function(x) x$success))
  
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
  
  invisible(NULL)
}

cat("✅ 10_summary.R 已加载\n")