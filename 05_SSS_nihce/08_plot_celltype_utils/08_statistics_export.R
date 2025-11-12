# ==================================================================
# 08_plot_celltype_utils\08_statistics_export.R
# ==================================================================

generate_zone_celltype_statistics <- function(
    sample_ids, 
    results_list,
    seurat_list,
    CONFIG) {
  
  cat("\n📊 生成Zone-Celltype统计数据...\n")
  
  zone_celltype_stats <- do.call(rbind, lapply(
    sample_ids, 
    function(sid) {
      result <- results_list[[sid]]
      if (is.null(result)) return(NULL)
      
      comp <- result$zone_composition
      meta <- seurat_list[[sid]]@meta.data
      
      data.frame(
        sample_id = sid,
        tissue = unique(meta$tissue)[1],
        age = unique(meta$age)[1],
        zone = comp$density_zone,
        celltype = comp$celltype_clean,
        count = comp$count,
        zone_total = comp$total,
        percentage = round(comp$percentage, 2),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  zone_celltype_stats %>%
    dplyr::arrange(sample_id, zone, dplyr::desc(count))
}

save_zone_celltype_statistics <- function(
    stats_table, 
    output_dir, 
    prefix = "") {
  
  cat("\n💾 保存统计数据...\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  filename <- if (prefix != "") {
    paste0(prefix, "_zone_celltype_counts.csv")
  } else {
    "zone_celltype_counts.csv"
  }
  
  filepath <- file.path(output_dir, filename)
  write.csv(stats_table, filepath, row.names = FALSE)
  cat(sprintf("  ✅ %s (%d 行)\n", basename(filepath), 
              nrow(stats_table)))
  
  invisible(filepath)
}

print_statistics_summary <- function(stats_table) {
  cat("\n", strrep("=", 60), "\n")
  cat("📊 统计摘要\n")
  cat(strrep("=", 60), "\n\n")
  
  cat(sprintf("🧬 样本数: %d\n", 
              length(unique(stats_table$sample_id))))
  cat(sprintf("🔬 细胞类型: %d\n", 
              length(unique(stats_table$celltype))))
  cat(sprintf("📊 密度区域: %d\n", 
              length(unique(stats_table$zone))))
  cat(sprintf("📈 总记录数: %s\n", 
              format(nrow(stats_table), big.mark = ",")))
  
  cat("\n", strrep("=", 60), "\n")
  invisible(NULL)
}


#' 生成所有样本的Zone-Celltype汇总统计
#'
#' @param stats_table 统计数据表
#' @param output_dir 输出目录
#' @param gene_list_name 基因列表名称
#'
#' @return 保存的文件路径列表
#'
save_aggregated_statistics <- function(
    stats_table, 
    output_dir, 
    gene_list_name = NULL) {
  
  cat("\n📊 生成汇总统计...\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  prefix <- if (!is.null(gene_list_name)) {
    paste0(gene_list_name, "_")
  } else {
    ""
  }
  
  # 1. 按样本汇总
  sample_summary <- stats_table %>%
    dplyr::group_by(sample_id, tissue, age) %>%
    dplyr::summarise(
      total_cells = sum(count),
      n_celltypes = dplyr::n_distinct(celltype),
      n_zones = dplyr::n_distinct(zone),
      .groups = "drop"
    ) %>%
    dplyr::arrange(tissue, age, sample_id)
  
  sample_file <- file.path(
    output_dir, 
    paste0(prefix, "summary_by_sample.csv")
  )
  write.csv(sample_summary, sample_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(sample_file)))
  
  # 2. 按细胞类型汇总
  celltype_summary <- stats_table %>%
    dplyr::group_by(celltype, zone) %>%
    dplyr::summarise(
      total_count = sum(count),
      n_samples = dplyr::n_distinct(sample_id),
      mean_pct = mean(percentage),
      sd_pct = sd(percentage),
      .groups = "drop"
    ) %>%
    dplyr::arrange(zone, dplyr::desc(total_count))
  
  celltype_file <- file.path(
    output_dir, 
    paste0(prefix, "summary_by_celltype.csv")
  )
  write.csv(celltype_summary, celltype_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(celltype_file)))
  
  # 3. 按Zone汇总
  zone_summary <- stats_table %>%
    dplyr::group_by(zone) %>%
    dplyr::summarise(
      total_cells = sum(count),
      n_samples = dplyr::n_distinct(sample_id),
      n_celltypes = dplyr::n_distinct(celltype),
      mean_cells_per_sample = mean(zone_total),
      .groups = "drop"
    ) %>%
    dplyr::arrange(zone)
  
  zone_file <- file.path(
    output_dir, 
    paste0(prefix, "summary_by_zone.csv")
  )
  write.csv(zone_summary, zone_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(zone_file)))
  
  # 4. Top细胞类型(按Zone)
  top_celltypes <- stats_table %>%
    dplyr::group_by(zone, celltype) %>%
    dplyr::summarise(
      total_count = sum(count),
      avg_percentage = mean(percentage),
      .groups = "drop"
    ) %>%
    dplyr::group_by(zone) %>%
    dplyr::slice_max(order_by = total_count, n = 10) %>%
    dplyr::arrange(zone, dplyr::desc(total_count))
  
  top_file <- file.path(
    output_dir, 
    paste0(prefix, "top10_celltypes_by_zone.csv")
  )
  write.csv(top_celltypes, top_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(top_file)))
  
  # 5. 全局统计
  global_stats <- list(
    gene_list = gene_list_name %||% "unknown",
    total_samples = length(unique(stats_table$sample_id)),
    total_celltypes = length(unique(stats_table$celltype)),
    total_zones = length(unique(stats_table$zone)),
    total_cells = sum(stats_table$count),
    generated_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  global_file <- file.path(
    output_dir, 
    paste0(prefix, "global_stats.txt")
  )
  
  sink(global_file)
  cat("==============================================\n")
  cat("  Zone-Celltype Global Statistics\n")
  cat("==============================================\n\n")
  for (key in names(global_stats)) {
    cat(sprintf("%-20s: %s\n", key, global_stats[[key]]))
  }
  sink()
  cat(sprintf("  ✅ %s\n", basename(global_file)))
  
  return(list(
    sample = sample_file,
    celltype = celltype_file,
    zone = zone_file,
    top = top_file,
    global = global_file
  ))
}