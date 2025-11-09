
# ==================================================================
# 08_statistics_export.R
# 生成zone-celltype统计表格
# ==================================================================

generate_zone_celltype_statistics <- function(
    sample_ids, 
    results_list, 
    CONFIG) {
  
  cat("\n📊 生成Zone-Celltype统计数据...\n")
  
  # ==============================================================
  # 主统计表：sample × zone × celltype 计数
  # ==============================================================
  
  zone_celltype_stats <- do.call(rbind, lapply(
    sample_ids, 
    function(sid) {
      result <- results_list[[sid]]
      if (is.null(result)) return(NULL)
      
      comp <- result$zone_composition
      data.frame(
        sample_id = sid,
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
  
  # ==============================================================
  # 样本级汇总统计
  # ==============================================================
  
  sample_summary <- do.call(rbind, lapply(
    sample_ids, 
    function(sid) {
      result <- results_list[[sid]]
      if (is.null(result)) return(NULL)
      
      data.frame(
        sample_id = sid,
        total_spots = result$stats$n_spots,
        high_density_spots = result$stats$n_high_density,
        high_density_pct = round(
          100 * result$stats$n_high_density / result$stats$n_spots, 
          2
        ),
        n_celltypes = result$stats$n_celltypes,
        n_zones = CONFIG$params$n_zones,
        stringsAsFactors = FALSE
      )
    }
  ))
  
  # ==============================================================
  # Zone汇总（跨样本）
  # ==============================================================
  
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
  
  # ==============================================================
  # 细胞类型汇总（跨样本和zone）
  # ==============================================================
  
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
  
  # ==============================================================
  # 宽格式矩阵：样本×zone（每个细胞类型一个矩阵）
  # ==============================================================
  
  wide_matrices <- list()
  
  unique_celltypes <- sort(unique(zone_celltype_stats$celltype))
  
  for (ct in unique_celltypes) {
    wide_matrices[[ct]] <- zone_celltype_stats %>%
      dplyr::filter(celltype == ct) %>%
      dplyr::select(sample_id, zone, count) %>%
      tidyr::pivot_wider(
        names_from = zone,
        names_prefix = "Zone_",
        values_from = count,
        values_fill = 0
      )
  }
  
  list(
    main_table = zone_celltype_stats,
    sample_summary = sample_summary,
    zone_summary = zone_summary,
    celltype_summary = celltype_summary,
    wide_matrices = wide_matrices
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
    if (prefix != "") {
      paste0(prefix, "_", name)
    } else {
      name
    }
  }
  
  saved_files <- list()
  
  # 1. 主表：sample × zone × celltype
  file1 <- file.path(
    output_dir, 
    add_prefix("zone_celltype_counts.csv")
  )
  write.csv(
    stats_list$main_table, 
    file1, 
    row.names = FALSE
  )
  cat(sprintf("  ✅ %s (%d 行)\n", 
              basename(file1), 
              nrow(stats_list$main_table)))
  saved_files$main_table <- file1
  
  # 2. 样本汇总
  file2 <- file.path(
    output_dir, 
    add_prefix("sample_summary.csv")
  )
  write.csv(
    stats_list$sample_summary, 
    file2, 
    row.names = FALSE
  )
  cat(sprintf("  ✅ %s\n", basename(file2)))
  saved_files$sample_summary <- file2
  
  # 3. Zone汇总
  file3 <- file.path(
    output_dir, 
    add_prefix("zone_summary.csv")
  )
  write.csv(
    stats_list$zone_summary, 
    file3, 
    row.names = FALSE
  )
  cat(sprintf("  ✅ %s\n", basename(file3)))
  saved_files$zone_summary <- file3
  
  # 4. 细胞类型汇总
  file4 <- file.path(
    output_dir, 
    add_prefix("celltype_summary.csv")
  )
  write.csv(
    stats_list$celltype_summary, 
    file4, 
    row.names = FALSE
  )
  cat(sprintf("  ✅ %s\n", basename(file4)))
  saved_files$celltype_summary <- file4
  
  # 5. 宽格式矩阵（每个细胞类型单独文件）
  if (length(stats_list$wide_matrices) > 0) {
    matrix_dir <- file.path(output_dir, "wide_matrices")
    if (!dir.exists(matrix_dir)) {
      dir.create(matrix_dir, showWarnings = FALSE)
    }
    
    saved_files$wide_matrices <- list()
    
    for (ct in names(stats_list$wide_matrices)) {
      safe_name <- gsub("[^A-Za-z0-9_-]", "_", ct)
      file_ct <- file.path(
        matrix_dir, 
        sprintf("%s.csv", safe_name)
      )
      write.csv(
        stats_list$wide_matrices[[ct]], 
        file_ct, 
        row.names = FALSE
      )
      saved_files$wide_matrices[[ct]] <- file_ct
    }
    
    cat(sprintf(
      "  ✅ wide_matrices/ (%d 个细胞类型)\n", 
      length(stats_list$wide_matrices)
    ))
  }
  
  invisible(saved_files)
}

print_statistics_summary <- function(stats_list) {
  cat("\n" %+% strrep("=", 60) %+% "\n")
  cat("📊 统计摘要\n")
  cat(strrep("=", 60) %+% "\n\n")
  
  cat(sprintf(
    "🧬 样本数: %d\n", 
    nrow(stats_list$sample_summary)
  ))
  
  cat(sprintf(
    "🔬 细胞类型: %d\n", 
    nrow(stats_list$celltype_summary)
  ))
  
  cat(sprintf(
    "📊 密度区域: %d\n", 
    length(unique(stats_list$main_table$zone))
  ))
  
  cat(sprintf(
    "📈 总记录数: %s\n", 
    format(nrow(stats_list$main_table), big.mark = ",")
  ))
  
  high_pct <- mean(stats_list$sample_summary$high_density_pct)
  cat(sprintf(
    "🔥 平均高密度占比: %.1f%%\n", 
    high_pct
  ))
  
  cat("\n" %+% "前5个细胞类型:\n")
  top5 <- head(stats_list$celltype_summary, 5)
  for (i in seq_len(nrow(top5))) {
    cat(sprintf(
      "  %d. %s: %s 个细胞\n",
      i,
      top5$celltype[i],
      format(top5$total_count[i], big.mark = ",")
    ))
  }
  
  cat("\n" %+% strrep("=", 60) %+% "\n")
  
  invisible(NULL)
}
