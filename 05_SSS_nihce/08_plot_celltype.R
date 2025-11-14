# ==================================================
# 08_plot_celltype.R (主入口)
# Version: 2.7.1
# ==================================================
## —— 快照：谁踩了 CONFIG$gene_list_path ——
cat(sprintf("%s: '%s'  class=%s  len=%d\n",
            basename(getSrcDirectory(function() NULL)),
            CONFIG$gene_list_path,
            class(CONFIG$gene_list_path),
            length(CONFIG$gene_list_path)))

cat("🔧 加载 08_plot_celltype.R (v2.7.1)...\n")

# 初始化脚本目录
if (!exists("script_dir")) {
  current_script <- tryCatch({
    normalizePath(sys.frame(1)$ofile, winslash = "/")
  }, error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      sub("^--file=", "", file_arg)
    } else {
      file.path(getwd(), "08_plot_celltype.R")
    }
  })
  script_dir <- dirname(current_script)
  cat(sprintf("   📂 脚本目录: %s\n", script_dir))
}

# 加载工具函数
utils_dir <- file.path(script_dir, "08_plot_celltype_utils")
if (!dir.exists(utils_dir)) {
  stop(sprintf("❌ 工具函数目录不存在: %s", utils_dir))
}

cat("   📦 加载工具函数...\n")
tryCatch({
  source(file.path(utils_dir, "00_operators.R"))
  source(file.path(utils_dir, "01_color_schemes.R"))
  source(file.path(utils_dir, "02_density_zones.R"))
  source(file.path(utils_dir, "03_plot_overlay.R"))
  source(file.path(utils_dir, "04_plot_composition.R"))
  source(file.path(utils_dir, "05_plot_heatmap.R"))
  source(file.path(utils_dir, "06_plot_combined.R"))
  source(file.path(utils_dir, "07_statistics.R"))
  source(file.path(utils_dir, "08_statistics_export.R"))
  source(file.path(utils_dir, "09_cache.R"))
  source(file.path(utils_dir, "10_summary.R"))
  source(file.path(utils_dir, "11_sample_process.R"))
  source(file.path(utils_dir, "12_analysis.R"))
}, error = function(e) {
  stop(sprintf("❌ 工具函数加载失败: %s", e$message))
})

# 主接口函数
analyze_celltype_niche <- function(
    sample_list, 
    CONFIG, 
    seurat_basename = NULL) {
  
  sample_ids <- names(sample_list)
  
  needs_init <- is.null(CONFIG$output) || 
    is.null(CONFIG$output$plot_dir) || 
    is.null(CONFIG$output$data_dir)
  
  if (needs_init) {
    cat("   🔧 初始化输出目录...\n")
    CONFIG <- init_output_dirs(CONFIG, seurat_basename)
  }
  
  for (dir_path in CONFIG$output) {
    if (!is.null(dir_path) && !dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, 
                 showWarnings = FALSE)
    }
  }
  
  run_celltype_analysis(
    data_list = sample_list, 
    sample_ids = sample_ids, 
    CONFIG = CONFIG
  )
}

# 初始化输出目录
init_output_dirs <- function(CONFIG, seurat_basename) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  base_figure_dir <- NULL
  base_metadata_dir <- NULL
  
  if (!is.null(CONFIG$dirs$figure)) {
    base_figure_dir <- CONFIG$dirs$figure
    base_metadata_dir <- CONFIG$metadata_dir %||% 
      file.path(dirname(base_figure_dir), "metadata")
  } else if (!is.null(CONFIG$figure_dir)) {
    base_figure_dir <- CONFIG$figure_dir
    base_metadata_dir <- CONFIG$metadata_dir %||% 
      file.path(dirname(base_figure_dir), "metadata")
  } else if (!is.null(CONFIG$output_dir)) {
    base_figure_dir <- file.path(CONFIG$output_dir, "figure")
    base_metadata_dir <- file.path(CONFIG$output_dir, "metadata")
  } else if (!is.null(CONFIG$output_base_dir)) {
    base_dir <- if (!is.null(seurat_basename)) {
      file.path(CONFIG$output_base_dir, seurat_basename)
    } else {
      CONFIG$output_base_dir
    }
    base_figure_dir <- file.path(base_dir, "figure")
    base_metadata_dir <- file.path(base_dir, "metadata")
  } else {
    stop("❌ 无法推断输出目录")
  }
  
  CONFIG$output <- list(
    base_dir = dirname(base_figure_dir),
    plot_dir = file.path(base_figure_dir, "celltype"),
    data_dir = file.path(base_metadata_dir, "celltype")
  )
  
  cat(sprintf("      📊 图形: %s\n", CONFIG$output$plot_dir))
  cat(sprintf("      📁 数据: %s\n", CONFIG$output$data_dir))
  
  CONFIG
}

cat("✅ 08_plot_celltype.R 加载完成 (v2.7.1)\n")
cat("📚 主要函数:\n")
cat("  - analyze_celltype_niche(sample_list, CONFIG)\n")
cat("  - list_cache_info(CONFIG)\n")
cat("  - clear_all_cache(CONFIG)\n\n")