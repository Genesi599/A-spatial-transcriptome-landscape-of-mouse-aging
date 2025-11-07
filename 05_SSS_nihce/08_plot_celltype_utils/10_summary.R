# ===================================================================
# 10_summary.R (简化版)
# 结果汇总打印模块
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 打印分析完成汇总
#' 
#' @param results_list 结果列表（来自 run_celltype_analysis）
#' @param CONFIG 配置对象
#' @param total_elapsed 总耗时（秒）
#' @param combined_data 合并数据（可选）
#'
#' @details
#' 打印内容：
#' - 成功/失败样本统计
#' - 每个样本的详细信息（spots数、high密度、zones数、细胞类型数）
#' - 总耗时和平均耗时
#' - 输出文件位置
#' - 合并数据统计
#'
print_analysis_summary <- function(results_list, CONFIG, total_elapsed, combined_data = NULL) {
  
  n_samples <- length(results_list)
  n_success <- sum(sapply(results_list, function(x) !is.null(x$stats)))
  n_failed <- n_samples - n_success
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                    📊 分析汇总报告                            ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  # ========================================
  # 1. 基本统计
  # ========================================
  
  cat("📈 样本处理统计:\n")
  cat(sprintf("   ✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, n_samples, 100 * n_success / n_samples))
  
  if (n_failed > 0) {
    cat(sprintf("   ❌ 失败: %d/%d (%.1f%%)\n", 
                n_failed, n_samples, 100 * n_failed / n_samples))
  }
  
  cat(sprintf("   ⏱️  总耗时: %.2f 秒 (%.2f 分钟)\n", 
              total_elapsed, total_elapsed / 60))
  
  if (n_success > 0) {
    cat(sprintf("   📊 平均耗时: %.2f 秒/样本\n", total_elapsed / n_samples))
  }
  
  cat("\n")
  
  # ========================================
  # 2. 样本详情表格
  # ========================================
  
  if (n_success > 0) {
    cat("📋 样本详情:\n")
    cat(sprintf("%-35s %10s %12s %10s %12s\n",
                "样本ID", "Spots", "High密度", "Zones", "细胞类型"))
    cat(paste(rep("─", 80), collapse = ""), "\n")
    
    total_spots <- 0
    total_high <- 0
    
    for (sid in names(results_list)) {
      res <- results_list[[sid]]
      
      if (!is.null(res$stats)) {
        cat(sprintf("%-35s %10d %12d %10d %12d\n",
                    substr(sid, 1, 35),  # 限制样本ID长度
                    res$stats$n_spots,
                    res$stats$n_high_density,
                    length(unique(res$zone_composition$density_zone)),
                    res$stats$n_celltypes))
        
        total_spots <- total_spots + res$stats$n_spots
        total_high <- total_high + res$stats$n_high_density
      }
    }
    
    if (n_success > 1) {
      cat(paste(rep("─", 80), collapse = ""), "\n")
      cat(sprintf("%-35s %10d %12d %10s %12s\n",
                  "总计", total_spots, total_high, "-", "-"))
    }
    
    cat("\n")
  }
  
  # ========================================
  # 3. 失败样本列表（如果有）
  # ========================================
  
  if (n_failed > 0) {
    cat("❌ 失败样本:\n")
    
    for (sid in names(results_list)) {
      res <- results_list[[sid]]
      if (is.null(res$stats)) {
        cat(sprintf("   • %s: 处理失败\n", sid))
      }
    }
    
    cat("\n")
  }
  
  # ========================================
  # 4. 输出文件位置
  # ========================================
  
  cat("📁 输出位置:\n")
  cat(sprintf("   • 图表目录: %s\n", CONFIG$output$plot_dir))
  cat(sprintf("   • 数据目录: %s\n", CONFIG$output$data_dir))
  
  # 统计文件大小
  if (dir.exists(CONFIG$output$plot_dir)) {
    plot_files <- list.files(CONFIG$output$plot_dir, full.names = TRUE, recursive = TRUE)
    total_plot_size <- sum(file.size(plot_files), na.rm = TRUE) / 1024^2
    cat(sprintf("   • 图表文件: %d 个 (%.2f MB)\n", 
                length(plot_files), total_plot_size))
  }
  
  if (dir.exists(CONFIG$output$data_dir)) {
    data_files <- list.files(CONFIG$output$data_dir, full.names = TRUE, recursive = TRUE)
    total_data_size <- sum(file.size(data_files), na.rm = TRUE) / 1024^2
    cat(sprintf("   • 数据文件: %d 个 (%.2f MB)\n", 
                length(data_files), total_data_size))
  }
  
  cat("\n")
  
  # ========================================
  # 5. 合并数据统计
  # ========================================
  
  if (!is.null(combined_data) && nrow(combined_data) > 0) {
    cat("📊 合并数据统计:\n")
    cat(sprintf("   • 总记录数: %d\n", nrow(combined_data)))
    cat(sprintf("   • 细胞类型: %d 种\n", length(unique(combined_data$celltype_clean))))
    cat(sprintf("   • 密度区域: %d 个\n", length(unique(combined_data$density_zone))))
    cat(sprintf("   • 样本数: %d\n", length(unique(combined_data$sample))))
    
    # 数据完整性
    completeness <- 100 * (1 - sum(is.na(combined_data)) / (nrow(combined_data) * ncol(combined_data)))
    cat(sprintf("   • 数据完整性: %.2f%%\n", completeness))
    
    cat("\n")
  }
  
  # ========================================
  # 6. 结束标志
  # ========================================
  
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                   ✅ 分析完成！                               ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  invisible(NULL)
}


#' 打印简短进度信息
#' 
#' @param i 当前索引
#' @param total 总数
#' @param sample_id 样本ID
#' @param status 状态（"处理中"/"完成"/"失败"）
#' @param extra_info 额外信息（可选）
#'
print_progress <- function(i, total, sample_id, status = "处理中", extra_info = NULL) {
  
  # 状态图标
  status_icon <- switch(status,
    "处理中" = "⏳",
    "完成" = "✅",
    "失败" = "❌",
    "跳过" = "⚠️",
    "⚪"  # 默认
  )
  
  # 构建输出
  output <- sprintf("[%2d/%2d] %s %s", i, total, status_icon, sample_id)
  
  if (!is.null(extra_info)) {
    output <- paste0(output, " - ", extra_info)
  }
  
  cat(output, "\n")
}


#' 打印细胞类型颜色映射
#' 
#' @param CONFIG 配置对象（必须包含 CONFIG$colors$celltype）
#' @param max_display 最多显示的细胞类型数量
#'
print_color_mapping <- function(CONFIG, max_display = 10) {
  
  if (is.null(CONFIG$colors) || is.null(CONFIG$colors$celltype)) {
    cat("⚠️  全局颜色方案未初始化\n")
    return(invisible(NULL))
  }
  
  celltype_colors <- CONFIG$colors$celltype
  n_celltypes <- length(celltype_colors)
  
  cat("\n🎨 全局细胞类型颜色映射:\n")
  
  if (n_celltypes <= max_display) {
    # 显示所有
    for (ct in names(celltype_colors)) {
      cat(sprintf("   • %-30s → %s\n", ct, celltype_colors[ct]))
    }
  } else {
    # 只显示前 max_display 个
    celltypes_to_show <- names(celltype_colors)[1:max_display]
    for (ct in celltypes_to_show) {
      cat(sprintf("   • %-30s → %s\n", ct, celltype_colors[ct]))
    }
    cat(sprintf("   ... 还有 %d 个细胞类型\n", n_celltypes - max_display))
  }
  
  cat("\n")
  
  invisible(NULL)
}


#' 打印zone颜色映射
#' 
#' @param CONFIG 配置对象（必须包含 CONFIG$colors$density_zone）
#'
print_zone_colors <- function(CONFIG) {
  
  if (is.null(CONFIG$colors) || is.null(CONFIG$colors$density_zone)) {
    cat("⚠️  密度区域颜色方案未初始化\n")
    return(invisible(NULL))
  }
  
  zone_colors <- CONFIG$colors$density_zone
  
  cat("\n🎨 密度区域颜色映射:\n")
  cat("   (Zone_0=核心/高密度 → Zone_N=外围/低密度)\n\n")
  
  for (zone in names(zone_colors)) {
    cat(sprintf("   • %-10s → %s\n", zone, zone_colors[zone]))
  }
  
  cat("\n")
  
  invisible(NULL)
}


#' 生成Markdown格式的报告摘要
#' 
#' @param results_list 结果列表
#' @param CONFIG 配置对象
#' @param total_elapsed 总耗时
#' @param output_file 输出文件路径（默认在data_dir中）
#'
generate_markdown_summary <- function(results_list, CONFIG, total_elapsed, 
                                     output_file = NULL) {
  
  if (is.null(output_file)) {
    output_file <- file.path(CONFIG$output$data_dir, "analysis_summary.md")
  }
  
  n_samples <- length(results_list)
  n_success <- sum(sapply(results_list, function(x) !is.null(x$stats)))
  
  # 构建Markdown内容
  md_content <- c(
    "# 细胞类型密度分布分析报告",
    "",
    sprintf("**生成时间**: %s", Sys.time()),
    sprintf("**分析样本**: %d", n_samples),
    sprintf("**成功样本**: %d (%.1f%%)", n_success, 100 * n_success / n_samples),
    sprintf("**总耗时**: %.2f 秒 (%.2f 分钟)", total_elapsed, total_elapsed / 60),
    "",
    "## 样本详情",
    "",
    "| 样本ID | Spots | High密度 | Zones | 细胞类型 |",
    "|--------|-------|----------|-------|----------|"
  )
  
  for (sid in names(results_list)) {
    res <- results_list[[sid]]
    if (!is.null(res$stats)) {
      md_content <- c(md_content,
        sprintf("| %s | %d | %d | %d | %d |",
                sid,
                res$stats$n_spots,
                res$stats$n_high_density,
                length(unique(res$zone_composition$density_zone)),
                res$stats$n_celltypes)
      )
    }
  }
  
  md_content <- c(md_content,
    "",
    "## 输出文件",
    "",
    sprintf("- **图表目录**: `%s`", CONFIG$output$plot_dir),
    sprintf("- **数据目录**: `%s`", CONFIG$output$data_dir),
    "",
    "## 参数配置",
    "",
    sprintf("- **密度阈值**: %.2f", CONFIG$params$density_threshold_percentile),
    sprintf("- **区域数量**: %d", CONFIG$params$n_zones),
    "",
    "---",
    "",
    "*Report generated by 08_plot_celltype.R*"
  )
  
  # 写入文件
  writeLines(md_content, output_file)
  
  cat(sprintf("📝 Markdown报告已保存: %s\n", output_file))
  
  invisible(output_file)
}

cat("✅ 10_summary.R 已加载（简化版）\n")