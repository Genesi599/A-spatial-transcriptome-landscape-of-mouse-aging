#!/usr/bin/env Rscript
# ===================================================================
# 环境设置和包加载
# ===================================================================

setup_environment <- function(config) {
  # 设置工作目录
  if (!is.null(config$work_dir) && config$work_dir != "") {
    if (!dir.exists(config$work_dir)) {
      dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
    }
    setwd(config$work_dir)
    cat(sprintf("✓ 工作目录: %s\n", config$work_dir))
  }
  
  # 创建所有必要的目录
  if (!is.null(config$dirs)) {
    cat("📁 创建输出目录...\n")
    
    for (dir_name in names(config$dirs)) {
      dir_path <- config$dirs[[dir_name]]
      
      if (is.null(dir_path) || dir_path == "") {
        warning(sprintf("⚠️  跳过空路径: %s", dir_name))
        next
      }
      
      if (!dir.exists(dir_path)) {
        tryCatch({
          dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
          cat(sprintf("  ✓ %s: %s\n", dir_name, dir_path))
        }, error = function(e) {
          warning(sprintf("⚠️  无法创建目录 %s: %s", dir_path, e$message))
        })
      }
    }
    cat("\n")
  }
  
  # 同时创建基础目录（防止遗漏）
  base_dirs <- c(
    config$output_dir,
    config$cache_dir,
    config$figure_dir,
    config$metadata_dir
  )
  
  for (dir_path in base_dirs) {
    if (!is.null(dir_path) && dir_path != "" && !dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(invisible(NULL))
}

load_packages <- function() {
  cat("📦 加载所需R包...\n")
  
  packages <- c(
    "Seurat", "tidyverse", "Matrix", 
    "future", "future.apply", "ggnewscale",
    "RColorBrewer", "patchwork", "digest", 
    "akima", "pheatmap", "tidyr", "tibble"
  )
  
  suppressPackageStartupMessages({
    for (pkg in packages) {
      library(pkg, character.only = TRUE)
    }
  })
  
  cat("✅ 所有包加载完成\n\n")
}

load_custom_functions <- function() {
  cat("📚 加载自定义函数...\n")
  source("niche_marker.R")
  source("SSS_isoheight_plot.R")
  cat("✅ 函数加载完成\n\n")
}