#!/usr/bin/env Rscript
# ===================================================================
# 12_file_utils.R
# ===================================================================
## —— 快照：谁踩了 CONFIG$gene_list_path ——
cat(sprintf("%s: '%s'  class=%s  len=%d\n",
            basename(getSrcDirectory(function() NULL)),
            CONFIG$gene_list_path,
            class(CONFIG$gene_list_path),
            length(CONFIG$gene_list_path)))

scan_seurat_files <- function(config) {
  cat("\n", strrep("=", 60), "\n")
  cat("扫描输入文件\n")
  cat(strrep("=", 60), "\n\n")
  
  seurat_files <- if (config$batch_mode) {
    scan_batch_files(config)
  } else {
    scan_single_file(config)
  }
  
  cat("\n")
  return(seurat_files)
}

scan_batch_files <- function(config) {
  cat(sprintf("📁 扫描目录: %s\n", config$seurat_dir))
  cat(sprintf("🔍 文件模式: %s\n", config$seurat_pattern))
  cat(sprintf("🔍 递归搜索: %s\n\n", config$recursive_search))
  
  seurat_dir <- config$seurat_dir %||% ""
  
  if (seurat_dir == "" || !dir.exists(seurat_dir)) {
    stop(sprintf("❌ 目录不存在: %s", seurat_dir))
  }
  
  seurat_files <- list.files(
    path = seurat_dir,
    pattern = config$seurat_pattern,
    full.names = TRUE,
    recursive = config$recursive_search
  )
  
  if (length(seurat_files) == 0) {
    stop(sprintf("❌ 未找到匹配文件 (模式: %s)", 
                 config$seurat_pattern))
  }
  
  cat(sprintf("✅ 找到 %d 个文件\n", length(seurat_files)))
  
  if (!is.null(config$specific_files) || 
      !is.null(config$exclude_files)) {
    original_count <- length(seurat_files)
    seurat_files <- filter_seurat_files(seurat_files, config)
    cat(sprintf("📋 过滤后剩余 %d 个文件 (原始: %d)\n", 
                length(seurat_files), original_count))
  }
  
  return(seurat_files)
}

scan_single_file <- function(config) {
  if (!file.exists(config$seurat_path)) {
    stop(sprintf("❌ 文件不存在: %s", config$seurat_path))
  }
  
  seurat_files <- config$seurat_path
  cat(sprintf("📄 单文件模式: %s\n", basename(seurat_files)))
  
  return(seurat_files)
}

filter_seurat_files <- function(seurat_files, config) {
  if (!is.null(config$specific_files)) {
    basenames <- basename(seurat_files)
    seurat_files <- seurat_files[basenames %in% 
                                   config$specific_files]
  }
  
  if (!is.null(config$exclude_files)) {
    basenames <- basename(seurat_files)
    seurat_files <- seurat_files[!basenames %in% 
                                    config$exclude_files]
  }
  
  return(seurat_files)
}

print_file_list <- function(seurat_files) {
  cat(strrep("=", 60), "\n")
  cat("待处理文件列表\n")
  cat(strrep("=", 60), "\n\n")
  
  cat(sprintf("%-4s %-40s %10s\n", "No.", "文件名", "大小"))
  cat(strrep("-", 60), "\n")
  
  for (i in seq_along(seurat_files)) {
    file_size_gb <- file.size(seurat_files[i]) / (1024^3)
    cat(sprintf("%3d. %-40s %8.2f GB\n", 
                i, basename(seurat_files[i]), file_size_gb))
  }
  
  total_size_gb <- sum(file.size(seurat_files)) / (1024^3)
  cat(strrep("-", 60), "\n")
  cat(sprintf("%-45s %8.2f GB\n", "总计:", total_size_gb))
  cat("\n")
}

update_config_for_file <- function(seurat_path, base_config) {
  config <- base_config
  config$seurat_path <- seurat_path
  config$seurat_file <- seurat_path
  
  seurat_basename <- tools::file_path_sans_ext(
    basename(seurat_path))
  config$output_dir <- file.path(config$output_base_dir, 
                                 seurat_basename)
  
  config <- update_config_paths(config)
  return(config)
}

update_config_paths <- function(config) {
  config$cache_dir <- file.path(config$output_dir, "cache")
  config$figure_dir <- file.path(config$output_dir, "figure")
  config$metadata_dir <- file.path(config$output_dir, "metadata")
  
  config$dirs <- list(
    cache = config$cache_dir,
    figure = config$figure_dir,
    metadata = config$metadata_dir,
    isoheight = file.path(config$figure_dir, "isoheight"),
    spatial = file.path(config$figure_dir, "spatial")
  )
  
  return(config)
}

validate_output_directory <- function(CONFIG) {
  if (is.null(CONFIG$output_base_dir) || 
      CONFIG$output_base_dir == "") {
    stop("❌ 未配置 output_base_dir")
  }
  
  if (!dir.exists(CONFIG$output_base_dir)) {
    cat(sprintf("📁 创建输出基础目录: %s\n", 
                CONFIG$output_base_dir))
    dir.create(CONFIG$output_base_dir, recursive = TRUE, 
               showWarnings = FALSE)
  }
}

load_gene_list_once <- function(CONFIG) {
  cat("\n【准备】加载基因列表\n")
  gene_list <- load_gene_list(CONFIG)
  cat(sprintf("✅ 加载了 %d 个基因\n\n", length(gene_list)))
  return(gene_list)
}


export_score_statistics <- function(seurat_obj, config, seurat_basename) {
  
  score_col <- config$score_column
  
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    warning(sprintf("评分列 %s 不存在", score_col))
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
    sprintf("%s_score_statistics.csv", seurat_basename)
  )
  
  write.csv(stats_df, output_file, row.names = FALSE)
  cat(sprintf("   💾 评分统计: %s\n", basename(output_file)))
  
  return(invisible(stats_df))
}

cat("✅ 12_file_utils.R 已加载\n")