#!/usr/bin/env Rscript
# ===================================================================
# 保存结果
# ===================================================================

save_results <- function(seurat_obj, config) {
  cat("💾 保存结果...\n")
  
  metadata_file <- file.path(
    config$metadata_dir, 
    sprintf("%s_metadata.csv", config$seurat_basename)
  )
  write.csv(seurat_obj@meta.data, metadata_file, row.names = TRUE)
  cat(sprintf("   ✅ Metadata: %s\n", basename(metadata_file)))
  
  export_score_statistics(seurat_obj, config)
  
  if (config$save_full_object) {
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
  
  score_col <- config$score_column_name %||% "Clock_Gene_Score"
  
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