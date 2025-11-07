#!/usr/bin/env Rscript
# ===================================================================
# 环境设置、包加载和全局配置
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 设置环境和创建输出目录
#'
#' @param config 配置列表
#'
setup_environment <- function(config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   环境初始化\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ---------------------------
  # 1. 设置工作目录
  # ---------------------------
  if (!is.null(config$work_dir) && config$work_dir != "") {
    if (!dir.exists(config$work_dir)) {
      dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
    }
    setwd(config$work_dir)
    cat(sprintf("✓ 工作目录: %s\n\n", config$work_dir))
  }
  
  # ---------------------------
  # 2. 创建所有输出目录
  # ---------------------------
  cat("📁 创建输出目录...\n")
  
  # 主要目录列表
  all_dirs <- c(
    config$output_dir,
    config$cache_dir,
    config$figure_dir,
    config$metadata_dir
  )
  
  # 添加 dirs 配置中的所有目录
  if (!is.null(config$dirs)) {
    all_dirs <- c(all_dirs, unlist(config$dirs))
  }
  
  # 去重并创建
  all_dirs <- unique(all_dirs[!is.na(all_dirs) & all_dirs != ""])
  
  for (dir_path in all_dirs) {
    if (!dir.exists(dir_path)) {
      tryCatch({
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
        cat(sprintf("  ✓ %s\n", dir_path))
      }, error = function(e) {
        warning(sprintf("⚠️  无法创建目录 %s: %s", dir_path, e$message))
      })
    }
  }
  
  cat("\n")
  
  return(invisible(NULL))
}


#' 加载所有必需的 R 包
#'
#' @param verbose 是否显示详细信息
#'
load_packages <- function(verbose = TRUE) {
  
  if (verbose) {
    cat("═══════════════════════════════════════════════════════════\n")
    cat("   加载 R 包\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
  }
  
  # ---------------------------
  # 核心包列表
  # ---------------------------
  core_packages <- c(
    # Seurat 生态
    "Seurat",
    
    # 数据处理
    "dplyr",
    "tidyr",
    "tibble",
    "Matrix",
    
    # 可视化
    "ggplot2",
    "patchwork",
    "RColorBrewer",
    "ggnewscale",
    "pheatmap",
    
    # 并行计算
    "future",
    "future.apply",
    "progressr",
    
    # 空间分析
    "RANN",
    "akima",
    
    # 工具包
    "digest"
  )
  
  # ---------------------------
  # 加载包
  # ---------------------------
  if (verbose) {
    cat("📦 加载核心包:\n")
  }
  
  loaded_packages <- character(0)
  failed_packages <- character(0)
  
  suppressPackageStartupMessages({
    for (pkg in core_packages) {
      tryCatch({
        library(pkg, character.only = TRUE, quietly = !verbose)
        loaded_packages <- c(loaded_packages, pkg)
        if (verbose) {
          cat(sprintf("  ✓ %s\n", pkg))
        }
      }, error = function(e) {
        failed_packages <- c(failed_packages, pkg)
        warning(sprintf("⚠️  无法加载 %s: %s", pkg, e$message))
      })
    }
  })
  
  if (verbose) {
    cat("\n")
    cat(sprintf("✅ 成功加载 %d/%d 个包\n", 
                length(loaded_packages), length(core_packages)))
    
    if (length(failed_packages) > 0) {
      cat(sprintf("❌ 加载失败: %s\n", paste(failed_packages, collapse = ", ")))
    }
    cat("\n")
  }
  
  return(invisible(list(
    loaded = loaded_packages,
    failed = failed_packages
  )))
}


#' 配置并行计算环境
#'
#' @param n_workers 并行线程数
#' @param memory_limit 内存限制（GB）
#'
setup_parallel <- function(n_workers = 4, memory_limit = 100) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   并行计算配置\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ---------------------------
  # 1. 设置 future 并行策略
  # ---------------------------
  future::plan(future::sequential)  # 先重置为串行
  
  cat(sprintf("🔧 并行线程数: %d\n", n_workers))
  cat(sprintf("💾 内存限制: %d GB\n", memory_limit))
  
  # ---------------------------
  # 2. 设置全局选项
  # ---------------------------
  options(
    future.globals.maxSize = Inf,  # 取消对象大小限制
    future.rng.onMisuse = "ignore"  # 忽略随机数警告
  )
  
  cat("✓ future 全局选项已设置\n")
  
  # ---------------------------
  # 3. 设置 progressr handlers（全局唯一设置）
  # ---------------------------
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  
  cat("✓ progressr 进度条已启用\n\n")
  
  return(invisible(NULL))
}


#' 加载自定义函数
#'
#' @param script_paths 脚本文件路径向量
#'
load_custom_functions <- function(script_paths = c("niche_marker.R", 
                                                   "SSS_isoheight_plot.R")) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   加载自定义函数\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  loaded_scripts <- character(0)
  failed_scripts <- character(0)
  
  for (script in script_paths) {
    if (file.exists(script)) {
      tryCatch({
        source(script)
        loaded_scripts <- c(loaded_scripts, script)
        cat(sprintf("  ✓ %s\n", script))
      }, error = function(e) {
        failed_scripts <- c(failed_scripts, script)
        warning(sprintf("⚠️  加载失败 %s: %s", script, e$message))
      })
    } else {
      warning(sprintf("⚠️  文件不存在: %s", script))
      failed_scripts <- c(failed_scripts, script)
    }
  }
  
  cat("\n")
  cat(sprintf("✅ 成功加载 %d/%d 个脚本\n\n", 
              length(loaded_scripts), length(script_paths)))
  
  if (length(failed_scripts) > 0) {
    cat(sprintf("❌ 加载失败: %s\n\n", paste(failed_scripts, collapse = ", ")))
  }
  
  return(invisible(list(
    loaded = loaded_scripts,
    failed = failed_scripts
  )))
}


#' 完整初始化流程（推荐使用）
#'
#' @param config 配置列表
#' @param custom_scripts 自定义脚本路径
#'
initialize_environment <- function(config, 
                                  custom_scripts = c("niche_marker.R", 
                                                    "SSS_isoheight_plot.R")) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║          Clock Gene Niche Analysis Pipeline              ║\n")
  cat("║                  环境初始化                               ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n")
  
  start_time <- Sys.time()
  
  # 1. 设置环境
  setup_environment(config)
  
  # 2. 加载包
  pkg_result <- load_packages(verbose = TRUE)
  
  # 3. 配置并行
  setup_parallel(
    n_workers = config$n_workers %||% 4,
    memory_limit = 100
  )
  
  # 4. 加载自定义函数
  script_result <- load_custom_functions(custom_scripts)
  
  # 汇总
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   初始化完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ R 包: %d 个已加载\n", length(pkg_result$loaded)))
  cat(sprintf("✅ 脚本: %d 个已加载\n", length(script_result$loaded)))
  cat(sprintf("⏱️  耗时: %.2f 秒\n\n", as.numeric(elapsed)))
  
  # 显示关键配置
  cat("📋 关键配置:\n")
  cat(sprintf("  - 工作目录: %s\n", getwd()))
  cat(sprintf("  - 输出目录: %s\n", config$output_dir))
  cat(sprintf("  - 并行线程: %d\n", config$n_workers %||% 4))
  cat(sprintf("  - 图形 DPI: %d\n", config$plot$dpi %||% 300))
  cat("\n")
  
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  return(invisible(list(
    packages = pkg_result,
    scripts = script_result,
    elapsed_time = as.numeric(elapsed)
  )))
}


# ===================================================================
# 辅助函数：%||% 操作符（如果左侧为NULL则返回右侧）
# ===================================================================
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}


# ===================================================================
# 导出函数列表（方便检查）
# ===================================================================
cat("✅ 01_setup.R 已加载\n")
cat("📚 可用函数:\n")
cat("  - setup_environment(config)\n")
cat("  - load_packages(verbose = TRUE)\n")
cat("  - setup_parallel(n_workers = 4, memory_limit = 100)\n")
cat("  - load_custom_functions(script_paths)\n")
cat("  - initialize_environment(config, custom_scripts)  [推荐]\n\n")