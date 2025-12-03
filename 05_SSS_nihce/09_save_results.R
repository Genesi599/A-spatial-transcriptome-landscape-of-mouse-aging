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
  
  # 0. 检查并修正 metadata_dir
  if (is.null(config$metadata_dir) || length(config$metadata_dir) == 0) {
    stop("config$metadata_dir 为空，请在 config 里指定输出目录")
  }
  if (!dir.exists(config$metadata_dir)) {
    dir.create(config$metadata_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 1. 生成 / 推断 seurat_basename
  basename <- config$seurat_basename
  if (is.null(basename) || length(basename) == 0 || basename == "") {
    # 尝试从其他字段推断
    if (!is.null(config$sample_id) && nzchar(config$sample_id)) {
      basename <- config$sample_id
      cat("   ℹ️ seurat_basename 未设置，使用 config$sample_id 作为 basename: ",
          basename, "\n")
    } else if (!is.null(config$sample_name) && nzchar(config$sample_name)) {
      basename <- config$sample_name
      cat("   ℹ️ seurat_basename 未设置，使用 config$sample_name 作为 basename: ",
          basename, "\n")
    } else if (!is.null(seurat_obj@project.name) && nzchar(seurat_obj@project.name)) {
      basename <- seurat_obj@project.name
      cat("   ℹ️ seurat_basename 未设置，使用 Seurat@project.name 作为 basename: ",
          basename, "\n")
    } else {
      stop("config$seurat_basename 为空，且无法从 sample_id / sample_name / project.name 推断，请在 config 里显式设置 seurat_basename")
    }
  }
  
  # 2. 组合 metadata 文件路径
  metadata_file <- file.path(
    config$metadata_dir, 
    sprintf("%s_metadata.csv", basename)
  )
  
  cat("   metadata_dir   :", config$metadata_dir, "\n")
  cat("   seurat_basename:", basename, "\n")
  cat("   metadata_file  :", metadata_file, "\n")
  
  # 3. 取 meta.data + cellid
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cellid")
  
  # 4. 提取空间坐标（用你统一的 GetAllCoordinates）
  coords_df <- GetAllCoordinates(seurat_obj)
  
  # 5. 合并 meta + 坐标
  meta_with_coords <- dplyr::left_join(
    meta_df,
    coords_df,
    by = "cellid"
  )
  
  # 6. 写出 CSV（包含 cellid, meta 列, row, col）
  write.csv(meta_with_coords, metadata_file, row.names = FALSE)
  cat(sprintf("   ✅ Metadata+coords: %s\n", basename(metadata_file)))
  
  # 7. 导出统计信息
  export_score_statistics(seurat_obj, config)
  
  # 8. 可选保存 RDS
  if (isTRUE(config$save_full_object)) {
    rds_file <- file.path(
      config$metadata_dir, 
      sprintf("%s_with_niche.rds", basename)
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