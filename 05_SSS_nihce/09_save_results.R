#!/usr/bin/env Rscript
# ===================================================================
# 保存结果
# ===================================================================
## —— 快照：谁踩了 CONFIG$gene_list_path ——
cat(sprintf("%s: '%s'  class=%s  len=%d\n",
            basename(getSrcDirectory(function() NULL)),
            CONFIG$gene_list_path,
            class(CONFIG$gene_list_path),
            length(CONFIG$gene_list_path)))

GetAllCoordinates <- function(.data) {
    .data@images %>%
        names() %>%
        unique() %>%
        map_dfr(~{
            GetTissueCoordinates(
                    .data,
                    image = .x,
                    cols = c("row", "col"),
                    scale = NULL
                ) %>%
            tibble::rownames_to_column(var = "cellid")
        })
}

save_results <- function(seurat_obj, config) {
  cat("💾 保存结果...\n")
  
  # 0. 检查输出目录
  if (is.null(config$metadata_dir) || length(config$metadata_dir) == 0) {
    stop("config$metadata_dir 为空，请在 CONFIG 或 update_config_for_file() 里设置")
  }
  if (!dir.exists(config$metadata_dir)) {
    dir.create(config$metadata_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 1. 检查 seurat_basename
  if (is.null(config$seurat_basename) || length(config$seurat_basename) == 0 || config$seurat_basename == "") {
    stop("config$seurat_basename 为空，请在 process_seurat_file() 里设置 config$seurat_basename")
  }
  
  metadata_file <- file.path(
    config$metadata_dir, 
    sprintf("%s_metadata.csv", config$seurat_basename)
  )
  
  cat("   metadata_dir   :", config$metadata_dir, "\n")
  cat("   seurat_basename:", config$seurat_basename, "\n")
  cat("   metadata_file  :", metadata_file, "\n")
  
  # 2. 取 meta + 坐标（注意这里已经有 cellid, row, col）
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cellid")
  
  coords_df <- GetAllCoordinates(seurat_obj)   # 你之前的统一坐标函数
  
  meta_with_coords <- dplyr::left_join(
    meta_df,
    coords_df,
    by = "cellid"
  )
  
  # 2.1 检查 ClockGene_High 是否存在（密度计算需要）
  if (!("ClockGene_High" %in% colnames(meta_with_coords))) {
    warning("meta_with_coords 缺少 'ClockGene_High' 列，无法计算密度 zone，仅导出 meta+coords")
    
    write.csv(meta_with_coords, metadata_file, row.names = FALSE)
    cat(sprintf("   ✅ Metadata+coords: %s（无 zone）\n", basename(metadata_file)))
    
  } else {
    # 2.2 计算密度 zone（使用你写的 calculate_density_zones）
    density_bins  <- config$params$n_zones        %||% 10
    expand_margin <- config$params$expand_margin  %||% 0.1
    
    density_res <- calculate_density_zones(
      df           = meta_with_coords,
      density_bins = density_bins,
      expand_margin = expand_margin,
      # 列名和你 meta_with_coords 中一致
      col_col      = "col",
      row_col      = "row",
      flag_col     = "ClockGene_High",
      quiet        = TRUE
    )
    
    if (is.null(density_res)) {
      warning("calculate_density_zones 返回 NULL（高表达点不足或 KDE 失败），仅导出 meta+coords")
      meta_with_zone <- meta_with_coords
    } else {
      # 2.3 把 zone 映射回每个 cell（按 col,row 对齐）
      meta_with_zone <- meta_with_coords %>%
        dplyr::left_join(
          density_res$spot_zones,
          by = c("col", "row")
        )
      # meta_with_zone 新增两列：
      #   - density_zone
      #   - density_value
    }
    
    # 2.4 写出包含 zone 的汇总表
    write.csv(meta_with_zone, metadata_file, row.names = FALSE)
    cat(sprintf("   ✅ Metadata+coords+zone: %s\n", basename(metadata_file)))
  }
  
  # 3. 导出打分统计（注意要传 seurat_basename）
  export_score_statistics(seurat_obj, config, config$seurat_basename)
  
  # 4. 可选保存 RDS
  if (isTRUE(config$save_full_object)) {
    rds_file <- file.path(
      config$metadata_dir, 
      sprintf("%s_with_niche.rds", config$seurat_basename)
    )
    saveRDS(seurat_obj, rds_file)
    cat(sprintf("   ✅ RDS: %s\n", basename(rds_file)))
  }
  
  cat("✅ 结果保存完成\n\n")
}

#' 导出模块评分统计
#'
#' @param seurat_obj Seurat对象
#' @param config 配置
#'
export_score_statistics <- function(seurat_obj, config) {
  
  score_col <- config$score_column_name %||% "ClockGene_Score1"
  
  if (is.null(score_col) || length(score_col) == 0 || score_col == "") {
    warning("评分列名未设置，跳过统计导出")
    return(invisible(NULL))
  }
  
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    warning(sprintf("评分列 %s 不存在于 meta.data", score_col))
    return(invisible(NULL))
  }
  
  scores <- seurat_obj@meta.data[[score_col]]
  
  stats_df <- data.frame(
    metric = c("mean", "median", "sd", "min", "max", 
               "q25", "q75", "n_cells"),
    value = c(
      mean(scores, na.rm = TRUE),
      median(scores, na.rm = TRUE),
      sd(scores, na.rm = TRUE),
      min(scores, na.rm = TRUE),
      max(scores, na.rm = TRUE),
      quantile(scores, 0.25, na.rm = TRUE),
      quantile(scores, 0.75, na.rm = TRUE),
      length(scores)
    )
  )
  
  output_file <- file.path(
    config$output_dir,
    sprintf("%s_score_statistics.csv", config$seurat_basename)
  )
  
  write.csv(stats_df, output_file, row.names = FALSE)
  cat(sprintf("   💾 评分统计: %s\n", basename(output_file)))
  
  return(invisible(stats_df))
}

print_summary <- function(config) {
  cat("📊 文件统计:\n")
  cat(sprintf("   图形文件夹: %s\n", config$figure_dir))
  cat(sprintf("   - Isoheight: %d 个文件\n", 
              length(list.files(config$dirs$isoheight))))
  cat(sprintf("   - Spatial: %d 个文件\n", 
              length(list.files(config$dirs$spatial))))
  cat("\n✅ 全部完成！\n")
}