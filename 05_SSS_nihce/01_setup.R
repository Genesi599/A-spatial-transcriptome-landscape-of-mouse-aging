#!/usr/bin/env Rscript
# ===================================================================
# 01_setup.R
# ===================================================================

setup_environment <- function(config) {
  cat("\n", strrep("=", 60), "\n")
  cat("环境初始化\n")
  cat(strrep("=", 60), "\n\n")
  
  if (!is.null(config$work_dir) && config$work_dir != "") {
    if (!dir.exists(config$work_dir)) {
      dir.create(config$work_dir, recursive = TRUE, 
                 showWarnings = FALSE)
    }
    setwd(config$work_dir)
    cat(sprintf("✓ 工作目录: %s\n\n", config$work_dir))
  }
  
  cat("📁 创建输出目录...\n")
  
  all_dirs <- c(
    config$output_dir,
    config$cache_dir,
    config$figure_dir,
    config$metadata_dir
  )
  
  if (!is.null(config$dirs)) {
    all_dirs <- c(all_dirs, unlist(config$dirs))
  }
  
  all_dirs <- unique(all_dirs[!is.na(all_dirs) & all_dirs != ""])
  
  for (dir_path in all_dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, 
                 showWarnings = FALSE)
      cat(sprintf("  ✓ %s\n", dir_path))
    }
  }
  
  cat("\n🔧 初始化分析参数...\n")
  
  if (!exists("%||%", mode = "function")) {
    `%||%` <- function(a, b) if (is.null(a)) b else a
  }
  
  if (is.null(config$params)) {
    config$params <- list()
  }
  
  config$params$celltype_col <- 
    config$params$celltype_col %||% "celltype"
  config$params$col_col <- config$params$col_col %||% "x"
  config$params$row_col <- config$params$row_col %||% "y"
  config$params$density_threshold_percentile <- 
    config$params$density_threshold_percentile %||% 0.95
  config$params$n_zones <- config$params$n_zones %||% 10
  config$params$grid_resolution <- 
    config$params$grid_resolution %||% 200
  
  cat(sprintf("  ✓ 细胞类型列: %s\n", 
              config$params$celltype_col))
  cat(sprintf("  ✓ 坐标列: %s, %s\n", 
              config$params$col_col, config$params$row_col))
  cat(sprintf("  ✓ 密度阈值: %.2f (分位数)\n", 
              config$params$density_threshold_percentile))
  cat(sprintf("  ✓ 区域数量: %d\n", config$params$n_zones))
  cat(sprintf("  ✓ 网格分辨率: %d\n", 
              config$params$grid_resolution))
  cat("\n")
  
  return(config)
}

load_packages <- function(verbose = TRUE) {
  if (verbose) {
    cat(strrep("=", 60), "\n")
    cat("加载 R 包\n")
    cat(strrep("=", 60), "\n\n")
  }
  
  core_packages <- c(
    "Seurat", "dplyr", "tidyr", "tibble", "Matrix",
    "ggplot2", "patchwork", "RColorBrewer", "ggnewscale", 
    "pheatmap", "future", "future.apply", "progressr",
    "RANN", "akima", "digest"
  )
  
  if (verbose) cat("📦 加载核心包:\n")
  
  loaded_packages <- character(0)
  failed_packages <- character(0)
  
  suppressPackageStartupMessages({
    for (pkg in core_packages) {
      result <- try({
        library(pkg, character.only = TRUE, quietly = !verbose)
        loaded_packages <- c(loaded_packages, pkg)
        if (verbose) cat(sprintf("  ✓ %s\n", pkg))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        failed_packages <- c(failed_packages, pkg)
        warning(sprintf("⚠️  无法加载 %s", pkg))
      }
    }
  })
  
  if (verbose) {
    cat("\n", sprintf("✅ 成功加载 %d/%d 个包\n", 
                      length(loaded_packages), 
                      length(core_packages)))
    if (length(failed_packages) > 0) {
      cat(sprintf("❌ 加载失败: %s\n", 
                  paste(failed_packages, collapse = ", ")))
    }
    cat("\n")
  }
  
  return(invisible(list(
    loaded = loaded_packages,
    failed = failed_packages
  )))
}

setup_parallel <- function(n_workers = 4, memory_limit = 100) {
  cat(strrep("=", 60), "\n")
  cat("并行计算配置\n")
  cat(strrep("=", 60), "\n\n")
  
  Sys.setenv(R_FUTURE_PLAN = "multisession")
  Sys.setenv(R_FUTURE_FORK_ENABLE = "false")
  future::plan(future::sequential)
  
  cat(sprintf("🔧 并行线程数: %d\n", n_workers))
  cat(sprintf("💾 内存限制: %d GB\n", memory_limit))
  
  options(
    future.globals.maxSize = Inf,
    future.rng.onMisuse = "ignore"
  )
  
  cat("✓ future 全局选项已设置\n")
  cat("✓ SLURM 检测已禁用\n")
  
  progressr::handlers(
    progressr::handler_progress(
      format = paste0("[:bar] :percent | 已完成: :current/:total",
                      " | 预计剩余: :eta | :message"),
      width = 80,
      complete = "=",
      clear = FALSE
    )
  )
  progressr::handlers(global = TRUE)
  
  cat("✓ progressr 进度条已启用（全局模式）\n\n")
  invisible(NULL)
}

load_custom_functions <- function(
    script_paths = c("niche_marker.R", "SSS_isoheight_plot.R")) {
  
  cat(strrep("=", 60), "\n")
  cat("加载自定义函数\n")
  cat(strrep("=", 60), "\n\n")
  
  loaded_scripts <- character(0)
  failed_scripts <- character(0)
  
  for (script in script_paths) {
    script_full <- file.path(getwd(), script)
    
    if (file.exists(script_full)) {
      result <- try({
        source(script_full)
        loaded_scripts <- c(loaded_scripts, script)
        cat(sprintf("  ✓ %s\n", script))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        failed_scripts <- c(failed_scripts, script)
        warning(sprintf("⚠️  加载失败 %s", script))
      }
    } else {
      warning(sprintf("⚠️  文件不存在: %s", script_full))
      failed_scripts <- c(failed_scripts, script)
    }
  }
  
  cat("\n", sprintf("✅ 成功加载 %d/%d 个脚本\n\n", 
                    length(loaded_scripts), 
                    length(script_paths)))
  
  if (length(failed_scripts) > 0) {
    cat(sprintf("❌ 加载失败: %s\n\n", 
                paste(failed_scripts, collapse = ", ")))
  }
  
  return(invisible(list(
    loaded = loaded_scripts,
    failed = failed_scripts
  )))
}

initialize_environment <- function(
    config, 
    custom_scripts = c("niche_marker.R", "SSS_isoheight_plot.R")) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("Clock Gene Niche Analysis Pipeline\n")
  cat("环境初始化\n")
  cat(strrep("=", 60), "\n")
  
  start_time <- Sys.time()
  
  config <- setup_environment(config)
  pkg_result <- load_packages(verbose = TRUE)
  setup_parallel(
    n_workers = config$n_workers %||% 4,
    memory_limit = 100
  )
  script_result <- load_custom_functions(custom_scripts)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  cat(strrep("=", 60), "\n")
  cat("初始化完成\n")
  cat(strrep("=", 60), "\n\n")
  
  cat(sprintf("✅ R 包: %d 个已加载\n", 
              length(pkg_result$loaded)))
  cat(sprintf("✅ 脚本: %d 个已加载\n", 
              length(script_result$loaded)))
  cat(sprintf("⏱️  耗时: %.2f 秒\n\n", as.numeric(elapsed)))
  
  cat("📋 关键配置:\n")
  cat(sprintf("  - 工作目录: %s\n", getwd()))
  cat(sprintf("  - 输出目录: %s\n", config$output_dir))
  cat(sprintf("  - 并行线程: %d\n", config$n_workers %||% 4))
  cat(sprintf("  - 图形 DPI: %d\n", config$plot$dpi %||% 300))
  cat(sprintf("  - 细胞类型列: %s\n", 
              config$params$celltype_col))
  cat("\n", strrep("=", 60), "\n\n")
  
  return(list(
    config = config,
    packages = pkg_result,
    scripts = script_result,
    elapsed_time = as.numeric(elapsed)
  ))
}

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 01_setup.R 已加载\n")