#!/usr/bin/env Rscript
# ===================================================================
# 文件操作工具模块
# ===================================================================

#' 扫描 Seurat 文件
#'
#' @param config 配置对象
#' @return 文件路径列表
#'
scan_seurat_files <- function(config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   扫描输入文件\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  if (config$batch_mode) {
    seurat_files <- scan_batch_files(config)
  } else {
    seurat_files <- scan_single_file(config)
  }
  
  cat("\n")
  
  return(seurat_files)
}


#' 批量模式扫描文件
#'
#' @param config 配置对象
#' @return 文件列表
#'
scan_batch_files <- function(config) {
  
  cat(sprintf("📁 扫描目录: %s\n", config$seurat_dir))
  cat(sprintf("🔍 文件模式: %s\n", config$seurat_pattern))
  cat(sprintf("🔍 递归搜索: %s\n\n", config$recursive_search))
  
  if (!dir.exists(config$seurat_dir)) {
    stop(sprintf("❌ 目录不存在: %s", config$seurat_dir))
  }
  
  seurat_files <- list.files(
    path = config$seurat_dir,
    pattern = config$seurat_pattern,
    full.names = TRUE,
    recursive = config$recursive_search
  )
  
  if (length(seurat_files) == 0) {
    stop(sprintf("❌ 未找到匹配文件 (模式: %s)", config$seurat_pattern))
  }
  
  cat(sprintf("✅ 找到 %d 个文件\n", length(seurat_files)))
  
  # 过滤文件
  if (!is.null(config$specific_files) || !is.null(config$exclude_files)) {
    original_count <- length(seurat_files)
    seurat_files <- filter_seurat_files(seurat_files, config)
    cat(sprintf("📋 过滤后剩余 %d 个文件 (原始: %d)\n", 
                length(seurat_files), original_count))
  }
  
  return(seurat_files)
}


#' 单文件模式扫描
#'
#' @param config 配置对象
#' @return 文件路径
#'
scan_single_file <- function(config) {
  
  if (!file.exists(config$seurat_path)) {
    stop(sprintf("❌ 文件不存在: %s", config$seurat_path))
  }
  
  seurat_files <- config$seurat_path
  cat(sprintf("📄 单文件模式: %s\n", basename(seurat_files)))
  
  return(seurat_files)
}


#' 过滤文件列表
#'
#' @param seurat_files 原始文件列表
#' @param config 配置对象
#' @return 过滤后的文件列表
#'
filter_seurat_files <- function(seurat_files, config) {
  
  # 特定文件过滤
  if (!is.null(config$specific_files)) {
    basenames <- basename(seurat_files)
    seurat_files <- seurat_files[basenames %in% config$specific_files]
  }
  
  # 排除文件过滤
  if (!is.null(config$exclude_files)) {
    basenames <- basename(seurat_files)
    seurat_files <- seurat_files[!basenames %in% config$exclude_files]
  }
  
  return(seurat_files)
}


#' 打印文件列表
#'
#' @param seurat_files 文件列表
#'
print_file_list <- function(seurat_files) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   待处理文件列表\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("%-4s %-50s %10s\n", "No.", "文件名", "大小"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (i in seq_along(seurat_files)) {
    file_size_gb <- file.size(seurat_files[i]) / (1024^3)
    cat(sprintf("%3d. %-50s %8.2f GB\n", 
                i, 
                basename(seurat_files[i]), 
                file_size_gb))
  }
  
  total_size_gb <- sum(file.size(seurat_files)) / (1024^3)
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-55s %8.2f GB\n", "总计:", total_size_gb))
  
  cat("\n")
}


#' 更新文件配置
#'
#' @param seurat_path Seurat 文件路径
#' @param base_config 基础配置
#' @return 更新后的配置
#'
update_config_paths <- function(config) {
  
  # 更新基础目录
  config$cache_dir <- file.path(config$output_dir, "cache")
  config$figure_dir <- file.path(config$output_dir, "figure")
  config$metadata_dir <- file.path(config$output_dir, "metadata")
  
  # 更新详细目录
  config$dirs <- list(
    cache = config$cache_dir,
    figure = config$figure_dir,
    metadata = config$metadata_dir,
    isoheight = file.path(config$figure_dir, "isoheight"),
    spatial = file.path(config$figure_dir, "spatial"),
    overlay = file.path(config$figure_dir, "isoheight", "01_overlay_plots"),
    celltype = file.path(config$figure_dir, "isoheight", "02_celltype_only"),
    composition = file.path(config$figure_dir, "isoheight", "03_composition_stats"),
    heatmaps = file.path(config$figure_dir, "isoheight", "04_heatmaps"),
    combined = file.path(config$figure_dir, "isoheight", "05_combined_analysis")
  )
  
  # ✅ 添加 output 结构（用于 08_plot_celltype.R）
  config$output <- list(
    base_dir = config$output_dir,
    plot_dir = file.path(config$figure_dir, "celltype"),
    data_dir = file.path(config$metadata_dir, "celltype")
  )
  
  return(config)
}


#' 更新配置路径
#'
#' @param config 配置对象
#' @return 更新后的配置
#'
update_config_paths <- function(config) {
  
  # 更新基础目录
  config$cache_dir <- file.path(config$output_dir, "cache")
  config$figure_dir <- file.path(config$output_dir, "figure")
  config$metadata_dir <- file.path(config$output_dir, "metadata")
  
  # 更新详细目录
  config$dirs <- list(
    cache = config$cache_dir,
    figure = config$figure_dir,
    metadata = config$metadata_dir,
    isoheight = file.path(config$figure_dir, "isoheight"),
    spatial = file.path(config$figure_dir, "spatial"),
    overlay = file.path(config$figure_dir, "isoheight", "01_overlay_plots"),
    celltype = file.path(config$figure_dir, "isoheight", "02_celltype_only"),
    composition = file.path(config$figure_dir, "isoheight", "03_composition_stats"),
    heatmaps = file.path(config$figure_dir, "isoheight", "04_heatmaps"),
    combined = file.path(config$figure_dir, "isoheight", "05_combined_analysis")
  )
  
  return(config)
}


#' 验证输出目录
#'
#' @param CONFIG 配置对象
#'
validate_output_directory <- function(CONFIG) {
  
  if (is.null(CONFIG$output_base_dir) || CONFIG$output_base_dir == "") {
    stop("❌ 未配置 output_base_dir")
  }
  
  if (!dir.exists(CONFIG$output_base_dir)) {
    cat(sprintf("📁 创建输出基础目录: %s\n", CONFIG$output_base_dir))
    dir.create(CONFIG$output_base_dir, recursive = TRUE, showWarnings = FALSE)
  }
}


#' 加载基因列表（仅一次）
#'
#' @param CONFIG 配置对象
#' @return 基因列表
#'
load_gene_list_once <- function(CONFIG) {
  
  cat("\n【准备】加载基因列表\n")
  gene_list <- load_gene_list(CONFIG)
  cat(sprintf("✅ 加载了 %d 个基因\n\n", length(gene_list)))
  
  return(gene_list)
}

cat("✅ 12_file_utils.R 已加载\n")