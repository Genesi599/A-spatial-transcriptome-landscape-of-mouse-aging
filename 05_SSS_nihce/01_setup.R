#!/usr/bin/env Rscript
# ===================================================================
# 环境设置和包加载
# ===================================================================

setup_environment <- function(config) {
  # 设置工作目录
  setwd(config$work_dir)
  
  # 创建目录
  for (dir_path in config$dirs) {
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  cat("✅ 工作目录:", getwd(), "\n")
  cat("✅ 输出目录:", config$output_dir, "\n")
  cat("✅ 图形目录:", config$figure_dir, "\n")
  cat("✅ 缓存目录:", config$cache_dir, "\n\n")
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