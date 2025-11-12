# ==================================================================
# 08_statistics_export.R
# 生成zone-celltype统计表格
# ==================================================================

generate_zone_celltype_statistics <- function(
    sample_ids, 
    results_list, 
    CONFIG) {
  
  cat("\n📊 生成Zone-Celltype统计数据...\n")
  
  zone_celltype_stats <- do.call(rbind, lapply(
    sample_ids, 
    function(sid) {
      result <- results_list[[sid]]
      if (is.null(result)) return(NULL)
      
      comp <- result$zone_composition
      
      sample_df <- result$sample_metadata
      tissue_val <- if (!is.null(sample_df$tissue)) {
        unique(sample_df$tissue)[1]
      } else NA_character_
      age_val <- if (!is.null(sample_df$age)) {
        unique(sample_df$age)[1]
      } else NA_character_
      
      data.frame(
        sample_id = sid,
        tissue = tissue_val,
        age = age_val,
        zone = comp$density_zone,
        celltype = comp$celltype_clean,
        count = comp$count,
        zone_total = comp$total,
        percentage = round(comp$percentage, 2),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  zone_celltype_stats <- zone_celltype_stats %>%
    dplyr::arrange(sample_id, zone, dplyr::desc(count))
  
  zone_summary <- zone_celltype_stats %>%
    dplyr::group_by(zone, celltype) %>%
    dplyr::summarise(
      total_count = sum(count),
      n_samples = dplyr::n(),
      mean_count = round(mean(count), 1),
      mean_pct = round(mean(percentage), 2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(zone, dplyr::desc(total_count))
  
  celltype_summary <- zone_celltype_stats %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      total_count = sum(count),
      n_samples = length(unique(sample_id)),
      n_zones = length(unique(zone)),
      mean_pct_per_zone = round(mean(percentage), 2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_count))
  
  list(
    main_table = zone_celltype_stats,
    zone_summary = zone_summary,
    celltype_summary = celltype_summary
  )
}

save_zone_celltype_statistics <- function(
    stats_list, 
    output_dir, 
    prefix = "") {
  
  cat("\n💾 保存统计数据...\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  add_prefix <- function(name) {
    if (prefix != "") paste0(prefix, "_", name) else name
  }
  
  saved_files <- list()
  
  file1 <- file.path(
    output_dir, 
    add_prefix("zone_celltype_counts.csv")
  )
  write.csv(stats_list$main_table, file1, row.names = FALSE)
  cat(sprintf("  ✅ %s (%d 行)\n", 
              basename(file1), nrow(stats_list$main_table)))
  saved_files$main_table <- file1
  
  file2 <- file.path(
    output_dir, 
    add_prefix("zone_summary.csv")
  )
  write.csv(stats_list$zone_summary, file2, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(file2)))
  saved_files$zone_summary <- file2
  
  file3 <- file.path(
    output_dir, 
    add_prefix("celltype_summary.csv")
  )
  write.csv(stats_list$celltype_summary, file3, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(file3)))
  saved_files$celltype_summary <- file3
  
  invisible(saved_files)
}

print_statistics_summary <- function(stats_list) {
  cat("\n" %+% strrep("=", 60) %+% "\n")
  cat("📊 统计摘要\n")
  cat(strrep("=", 60) %+% "\n\n")
  
  n_samples <- length(unique(stats_list$main_table$sample_id))
  cat(sprintf("🧬 样本数: %d\n", n_samples))
  
  cat(sprintf("🔬 细胞类型: %d\n", 
              nrow(stats_list$celltype_summary)))
  
  cat(sprintf("📊 密度区域: %d\n", 
              length(unique(stats_list$main_table$zone))))
  
  cat(sprintf("📈 总记录数: %s\n", 
              format(nrow(stats_list$main_table), big.mark = ",")))
  
  cat("\n" %+% "前5个细胞类型:\n")
  top5 <- head(stats_list$celltype_summary, 5)
  for (i in seq_len(nrow(top5))) {
    cat(sprintf(
      "  %d. %s: %s 个细胞\n",
      i, top5$celltype[i],
      format(top5$total_count[i], big.mark = ",")
    ))
  }
  
  cat("\n" %+% strrep("=", 60) %+% "\n")
  
  invisible(NULL)
}