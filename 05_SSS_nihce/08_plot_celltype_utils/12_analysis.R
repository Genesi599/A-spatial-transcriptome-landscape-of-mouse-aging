# ==================================================
# 12_analysis.R
# 综合分析函数
# ==================================================

create_global_color_scheme <- function(
    data_list, 
    celltype_col, 
    n_zones = 10) {
  
  cat("\n🎨 生成全局颜色方案...\n")
  
  all_celltypes <- unique(unlist(lapply(data_list, function(obj) {
    
    if (inherits(obj, "Seurat")) {
      ct <- obj@meta.data[[celltype_col]]
    } else if (is.data.frame(obj)) {
      ct <- obj[[celltype_col]]
    } else {
      warning(sprintf("未知对象类型: %s", class(obj)[1]))
      return(character(0))
    }
    
    ct <- as.character(ct)
    unique(ct[!is.na(ct) & ct != "" & ct != "Unknown"])
  })))
  
  all_celltypes <- sort(all_celltypes)
  n_celltypes <- length(all_celltypes)
  
  if (n_celltypes == 0) {
    stop("❌ 未找到任何细胞类型")
  }
  
  cat(sprintf("  📊 发现 %d 个细胞类型\n", n_celltypes))
  
  if (n_celltypes <= 10) {
    for (ct in all_celltypes) {
      cat(sprintf("     • %s\n", ct))
    }
  } else {
    for (i in 1:10) {
      cat(sprintf("     • %s\n", all_celltypes[i]))
    }
    cat(sprintf("     ... 还有 %d 个\n", n_celltypes - 10))
  }
  
  celltype_colors <- get_celltype_colors(all_celltypes)
  zone_colors <- get_zone_colors(n_zones)
  names(zone_colors) <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf(
    "\n  ✅ 完成 (%d 类型 + %d 区域)\n", 
    n_celltypes, n_zones
  ))
  
  list(celltype = celltype_colors, density_zone = zone_colors)
}

#' 保存全局颜色方案到 CSV 文件
#'
#' @param color_scheme 颜色方案对象（来自 create_global_color_scheme）
#' @param output_dir 输出目录（用于保存 CSV 文件的目录）
#' 
#' @return 无返回值
#'
save_global_color_scheme_to_csv <- function(color_scheme, output_dir) {
  # 创建数据框，只包含细胞类型和颜色
  color_data <- data.frame(
    Cell_Type = names(color_scheme$celltype),
    Color = color_scheme$celltype,
    stringsAsFactors = FALSE
  )
  
  # 定义 CSV 文件名和路径
  csv_file <- file.path(output_dir, "global_celltype_colors.csv")
  
  # 写入 CSV 文件
  write.csv(color_data, csv_file, row.names = FALSE)
  
  cat(sprintf("✅ 全局细胞类型颜色方案已保存到: %s\n", csv_file))
}

run_celltype_analysis <- function(
    data_list, 
    sample_ids, 
    CONFIG) {
  
  cat(paste0(
    "\n╔═══════════════════════════════════════════════",
    "═══════════╗\n"
  ))
  cat("║  细胞类型分布分析                            ")
  cat("          ║\n")
  cat(paste0(
    "╚═══════════════════════════════════════════════",
    "═══════════╝\n"
  ))
  
  if (CONFIG$debug_mode && 
      !is.null(CONFIG$cache_max_age_hours)) {
    clean_expired_cache(CONFIG)
  }
  
  if (CONFIG$debug_mode && !is.null(CONFIG$cache_dir)) {
    cat(sprintf(
      "\n🔧 调试模式开启\n📂 缓存: %s\n", 
      CONFIG$cache_dir
    ))
  }
  
  CONFIG$colors <- create_global_color_scheme(
    data_list, 
    CONFIG$params$celltype_col, 
    CONFIG$params$n_zones
  )

  # 在创建颜色方案后保存到 CSV
  save_global_color_scheme_to_csv(CONFIG$colors, CONFIG$output$data_dir)
  
  cat("\n🔬 处理样本...\n")
  
  results_list <- list()
  n_from_cache <- 0
  
  for (i in seq_along(data_list)) {
    cat(sprintf("\n[%2d/%2d]", i, length(data_list)))
    
    result <- process_single_sample(
      data_list[[i]], 
      sample_ids[i], 
      CONFIG
    )
    
    if (!is.null(result$from_cache) && result$from_cache) {
      n_from_cache <- n_from_cache + 1
    }
    
    results_list[[sample_ids[i]]] <- result
  }
  
  if (CONFIG$debug_mode && n_from_cache > 0) {
    cat(sprintf(
      "\n💾 缓存命中: %d/%d (%.1f%%)\n", 
      n_from_cache, 
      length(data_list), 
      100 * n_from_cache / length(data_list)
    ))
  }
  
  cat("\n\n📊 合并数据...\n")
  
  combined_data <- do.call(rbind, lapply(
    names(results_list), 
    function(sid) {
      comp <- results_list[[sid]]$zone_composition
      comp$sample <- sid
      comp
    }
  ))
  
  combined_file <- file.path(
    CONFIG$output$data_dir, 
    "combined_zone_composition.csv"
  )
  write.csv(combined_data, combined_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(combined_file)))
  
  cat("\n📊 生成热图...\n")
  p_heatmap <- plot_combined_heatmap(combined_data, CONFIG)
  heatmap_png <- file.path(
    CONFIG$output$plot_dir, 
    "combined_heatmap.png"
  )
  heatmap_pdf <- file.path(
    CONFIG$output$plot_dir, 
    "combined_heatmap.pdf"
  )
  ggsave(heatmap_png, p_heatmap, 
         width = 18, height = 14, dpi = 300, bg = "white")
  ggsave(heatmap_pdf, p_heatmap, 
         width = 18, height = 14, bg = "white")
  cat(sprintf("  ✅ %s (PNG + PDF)\n", "combined_heatmap"))
  
  cat("\n📊 生成综合图...\n")
  plots_combined <- plot_combined_analysis(combined_data, CONFIG)
  
  combined_png <- file.path(
    CONFIG$output$plot_dir, "combined_analysis.png"
  )
  combined_pdf <- file.path(
    CONFIG$output$plot_dir, "combined_analysis.pdf"
  )
  boxplot_png <- file.path(
    CONFIG$output$plot_dir, "combined_boxplot.png"
  )
  boxplot_pdf <- file.path(
    CONFIG$output$plot_dir, "combined_boxplot.pdf"
  )
  barplot_png <- file.path(
    CONFIG$output$plot_dir, "combined_barplot.png"
  )
  barplot_pdf <- file.path(
    CONFIG$output$plot_dir, "combined_barplot.pdf"
  )
  
  ggsave(combined_png, plots_combined$combined, 
         width = 20, height = 16, dpi = 300, bg = "white")
  ggsave(combined_pdf, plots_combined$combined, 
         width = 20, height = 16, bg = "white")
  ggsave(boxplot_png, plots_combined$boxplot, 
         width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(boxplot_pdf, plots_combined$boxplot, 
         width = 16, height = 12, bg = "white")
  ggsave(barplot_png, plots_combined$barplot, 
         width = 12, height = 7, dpi = 300, bg = "white")
  ggsave(barplot_pdf, plots_combined$barplot, 
         width = 12, height = 7, bg = "white")
  
  cat(sprintf("  ✅ %s (组合 + 单独)\n", "combined_analysis"))
  
  cat("\n📊 生成统计摘要...\n")
  summary_stats <- generate_summary_statistics(
    combined_data, CONFIG
  )
  summary_file <- file.path(
    CONFIG$output$data_dir, 
    "summary_statistics.csv"
  )
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(summary_file)))
  
  cat("\n📊 生成Zone-Celltype详细统计...\n")
  
  if (length(data_list) > 0) {
    first_obj <- data_list[[1]]
    cat(sprintf("   🔍 输入对象类型: %s\n", class(first_obj)[1]))
    if (inherits(first_obj, "Seurat")) {
      cat(sprintf("   ✅ meta.data列: %s\n", 
                  paste(head(names(first_obj@meta.data), 5), 
                        collapse = ", ")))
    }
  }
  
  zone_stats <- generate_zone_celltype_statistics(
    sample_ids, 
    results_list,
    data_list,
    CONFIG
  )
  
  saved_files <- save_zone_celltype_statistics(
    zone_stats, 
    CONFIG$output$data_dir
  )
  
  print_statistics_summary(zone_stats)
  
  gene_list_name <- if (!is.null(CONFIG$gene_list_path)) {
    tools::file_path_sans_ext(basename(CONFIG$gene_list_path))
  } else {
    NULL
  }
  
  summary_files <- save_aggregated_statistics(
    zone_stats, 
    CONFIG$output$data_dir,
    gene_list_name
  )
  
  cat("\n📦 汇总统计已生成\n")
  
  cat(paste0(
    "\n╔═══════════════════════════════════════════════",
    "═══════════╗\n"
  ))
  cat("║  ✅ 分析完成！                                ")
  cat("          ║\n")
  cat(paste0(
    "╚═══════════════════════════════════════════════",
    "═══════════╝\n\n"
  ))
  
  list(
    individual_results = results_list,
    combined_data = combined_data,
    summary_statistics = summary_stats,
    zone_celltype_stats = zone_stats,
    combined_plots = list(
      heatmap = p_heatmap, 
      analysis = plots_combined$combined,
      boxplot = plots_combined$boxplot,
      barplot = plots_combined$barplot
    ),
    config = CONFIG,
    cache_stats = list(
      n_from_cache = n_from_cache,
      total = length(data_list),
      cache_hit_rate = 100 * n_from_cache / length(data_list)
    ),
    saved_files = saved_files,
    summary_files = summary_files
  )
}