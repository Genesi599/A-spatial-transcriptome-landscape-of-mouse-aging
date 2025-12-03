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
  
  # 确保输出目录存在
  if (!dir.exists(config$metadata_dir)) {
    dir.create(config$metadata_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 1. 构建 metadata 输出路径
  metadata_file <- file.path(
    config$metadata_dir, 
    sprintf("%s_metadata.csv", config$seurat_basename)
  )
  
  # 2. meta.data 加上 cellid
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cellid")
  
  # 3. 用统一的 GetAllCoordinates() 提取所有 cellid 的 row/col
  coords_df <- GetAllCoordinates(seurat_obj)
  # 此时 coords_df 至少有: cellid, row, col
  
  # 4. 按 cellid 合并 meta 和 坐标
  # 用 left_join 确保所有 meta 里的 cellid 都保留，缺坐标的行 row/col 为 NA
  meta_with_coords <- dplyr::left_join(
    meta_df,
    coords_df,
    by = "cellid"
  )
  
  # 5. 写出 CSV：包含 cellid + meta 列 + row + col
  # 已经有 cellid 列，就没必要再用 row.names 了
  write.csv(meta_with_coords, metadata_file, row.names = FALSE)
  cat(sprintf("   ✅ Metadata+coords: %s\n", basename(metadata_file)))
  
  # 6. 导出你的统计信息（保持原逻辑）
  export_score_statistics(seurat_obj, config)
  
  # 7. 选择性保存完整 Seurat 对象
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