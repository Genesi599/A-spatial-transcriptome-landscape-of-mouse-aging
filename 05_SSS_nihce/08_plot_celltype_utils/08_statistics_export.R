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