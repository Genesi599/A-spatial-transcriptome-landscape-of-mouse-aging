# 项目导出: 05_SSS_nihce

**导出时间**: 2025-11-07 14:32:32
**项目路径**: C:/Users/yh109/Documents/GitHub/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce

---

## 目录结构

```
05_SSS_nihce
├── 00_config.R
├── 00_config_examples.R
├── 01_setup.R
├── 02_utils.R
├── 03_load_data.R
├── 04_module_score.R
├── 05_niche_analysis.R
├── 06_plot_isoheight.R
├── 07_plot_spatial.R
├── 08_plot_celltype.R
├── 08_plot_celltype_utils
│   ├── 00_operators.R
│   ├── 01_color_schemes.R
│   ├── 02_density_zones.R
│   ├── 03_plot_overlay.R
│   ├── 04_plot_composition.R
│   ├── 05_plot_heatmap.R
│   ├── 06_plot_combined.R
│   ├── 07_statistics.R
│   ├── 08_validation.R
│   ├── 09_save_plots.R
│   └── 10_summary.R
├── 09_save_results.R
├── 10_batch_processing.R
├── 11_sample_preprocessing.R
├── 12_file_utils.R
├── 13_reporting.R
├── AI_trans.R
├── main.R
├── niche_grade_entropy.R
├── niche_marker.R
└── SSS_isoheight_plot.R
```

## 文件统计

- **总文件数**: 31
- **文件类型分布**:
  - .R: 31 个
- **项目总大小**: 188.41 KB

## 文件内容

### 00_config.R

- **大小**: 3.91 KB
- **修改时间**: 2025-11-07 12:35:33

```r
# 00_config.R (添加缓存设置)

CONFIG <- list(
  # ===== 路径设置 =====
  work_dir = "/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce",
  output_base_dir = "/dellstorage09/quj_lab/yanghang/spatial",
  
  # 数据路径
  gene_list_path = "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt",
  
  # ✅ 缓存路径（新增）
  cache_dir = "/dellstorage09/quj_lab/yanghang/spatial/cache",
  
  # ===== 批量处理设置 =====
  batch_mode = FALSE,
  seurat_path = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Intestine_2-25M.rds",
  seurat_dir = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds",
  seurat_pattern = "\\.rds$",
  recursive_search = FALSE,
  specific_files = NULL,
  exclude_files = NULL,

  # ===== 分析参数 =====
  threshold_quantile = 0.95,
  niche_dist_method = "Euclidean",
  n_workers = 12,
  
  # ===== 绘图参数 =====
  plot = list(
    contour_bins = 8,
    point_size_bg = 0.3,
    point_size_top = 1.2,
    point_size_scatter = 2.5,
    contour_alpha = 0.25,
    interp_resolution = 200,
    expand_margin = 0.05,
    dpi = 300
  ),
  
  # ===== 调试参数 =====
  debug_mode = FALSE,
  debug_sample_limit = 3,
  save_full_object = FALSE,
  
  # ===== 缓存参数 =====
  cache_max_age_hours = NULL  # ✅ 缓存有效期（小时），NULL = 永久有效
)

# ===================================================================
# 初始化缓存目录
# ===================================================================

if (!is.null(CONFIG$cache_dir)) {
  if (!dir.exists(CONFIG$cache_dir)) {
    dir.create(CONFIG$cache_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("✅ 创建缓存目录: %s\n", CONFIG$cache_dir))
  } else {
    cat(sprintf("✅ 缓存目录: %s\n", CONFIG$cache_dir))
  }
}

# ===================================================================
# 打印配置信息
# ===================================================================

cat(sprintf("\n╔═══════════════════════════════════════════════════════════╗\n"))
cat(sprintf("║                    配置信息                                ║\n"))
cat(sprintf("╚═══════════════════════════════════════════════════════════╝\n\n"))

if (CONFIG$batch_mode) {
  cat(sprintf("运行模式: 批量处理\n"))
  cat(sprintf("输入目录: %s\n", CONFIG$seurat_dir))
  cat(sprintf("文件模式: %s\n", CONFIG$seurat_pattern))
  cat(sprintf("递归搜索: %s\n", ifelse(CONFIG$recursive_search, "是", "否")))
  
  if (!is.null(CONFIG$specific_files)) {
    cat(sprintf("指定文件: %s\n", paste(CONFIG$specific_files, collapse = ", ")))
  }
  
  if (!is.null(CONFIG$exclude_files)) {
    cat(sprintf("排除文件: %s\n", paste(CONFIG$exclude_files, collapse = ", ")))
  }
} else {
  cat(sprintf("运行模式: 单文件处理\n"))
  cat(sprintf("Seurat 文件: %s\n", basename(CONFIG$seurat_path)))
}

cat(sprintf("输出目录: %s\n", CONFIG$output_base_dir))
cat(sprintf("缓存目录: %s\n", ifelse(is.null(CONFIG$cache_dir), "未设置", CONFIG$cache_dir)))  # ✅ 新增
cat(sprintf("缓存有效期: %s\n", ifelse(is.null(CONFIG$cache_max_age_hours), "永久", sprintf("%d小时", CONFIG$cache_max_age_hours))))  # ✅ 新增
cat(sprintf("阈值分位数: %.2f (Top %.0f%%)\n", 
            CONFIG$threshold_quantile, 
            (1 - CONFIG$threshold_quantile) * 100))
cat(sprintf("并行工作数: %d\n", CONFIG$n_workers))
cat(sprintf("调试模式: %s\n", ifelse(CONFIG$debug_mode, "开启", "关闭")))
cat(sprintf("\n"))
```

---

### 00_config_examples.R

- **大小**: 2.76 KB
- **修改时间**: 2025-11-06 21:40:23

```r
#!/usr/bin/env Rscript
# ===================================================================
# 配置示例
# ===================================================================

# ===================================================================
# 示例 1: 批量处理目录中的所有 .rds 文件
# ===================================================================
CONFIG_EXAMPLE_1 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  recursive_search = FALSE,
  specific_files = NULL,
  exclude_files = NULL
)

# ===================================================================
# 示例 2: 批量处理指定的几个文件
# ===================================================================
CONFIG_EXAMPLE_2 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  specific_files = c("Lung_2-25M.rds", "Heart_2-25M.rds", "Brain_2-25M.rds"),
  exclude_files = NULL
)

# ===================================================================
# 示例 3: 批量处理，但排除某些文件
# ===================================================================
CONFIG_EXAMPLE_3 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  exclude_files = c("test.rds", "backup.rds", "old_version.rds")
)

# ===================================================================
# 示例 4: 递归搜索子目录
# ===================================================================
CONFIG_EXAMPLE_4 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  recursive_search = TRUE  # 会搜索所有子目录
)

# ===================================================================
# 示例 5: 单文件模式
# ===================================================================
CONFIG_EXAMPLE_5 <- list(
  batch_mode = FALSE,
  seurat_path = "/path/to/single_file.rds"
)

# ===================================================================
# 示例 6: 匹配特定命名模式的文件
# ===================================================================
CONFIG_EXAMPLE_6 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "^Lung.*\\.rds$",  # 只处理以 Lung 开头的文件
  recursive_search = FALSE
)

# ===================================================================
# 示例 7: 调试模式（只处理每个文件的前几个样本）
# ===================================================================
CONFIG_EXAMPLE_7 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  debug_mode = TRUE,
  debug_sample_limit = 2  # 每个文件只处理前2个样本
)
```

---

### 01_setup.R

- **大小**: 10.61 KB
- **修改时间**: 2025-11-07 13:41:55

```r
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
```

---

### 02_utils.R

- **大小**: 7 KB
- **修改时间**: 2025-11-07 12:32:28

```r


#!/usr/bin/env Rscript
# ===================================================================
# 工具函数
# ===================================================================

# -----------------------------
# 缓存管理
# -----------------------------
generate_cache_key <- function(...) {
  digest::digest(list(...), algo = "md5")
}

save_cache <- function(obj, file, desc = "") {
  tryCatch({
    saveRDS(obj, file)
    cat(sprintf("💾 缓存已保存: %s (%.2f MB) %s\n", 
                basename(file), file.size(file)/1024^2, 
                ifelse(desc != "", paste0("- ", desc), "")))
  }, error = function(e) {
    warning(sprintf("⚠️ 缓存保存失败: %s\n", e$message))
  })
}

load_cache <- function(file, desc = "") {
  if (file.exists(file)) {
    cat(sprintf("📂 加载缓存: %s (%.2f MB) %s\n", 
                basename(file), file.size(file)/1024^2,
                ifelse(desc != "", paste0("- ", desc), "")))
    return(readRDS(file))
  }
  NULL
}

is_cache_valid <- function(cache_file, max_age_hours = NULL) {
  # ✅ 检查 cache_file 是否为空或 NULL
  if (is.null(cache_file) || length(cache_file) == 0 || cache_file == "") {
    return(FALSE)
  }
  
  # 检查缓存文件是否存在
  if (!file.exists(cache_file)) {
    return(FALSE)
  }
  
  # 如果不限制缓存时间，只要文件存在就有效
  if (is.null(max_age_hours)) {
    return(TRUE)
  }
  
  # 检查缓存年龄
  cache_time <- file.info(cache_file)$mtime
  age_hours <- as.numeric(difftime(Sys.time(), cache_time, units = "hours"))
  
  return(age_hours < max_age_hours)
}

# -----------------------------
# 坐标获取（兼容多种格式）
# -----------------------------
get_coordinates_safely <- function(seurat_obj) {
  coord_attempts <- list(
    c("row", "col"),
    c("imagerow", "imagecol"),
    c("x", "y")
  )
  
  for (cols in coord_attempts) {
    coords <- tryCatch({
      GetTissueCoordinates(seurat_obj, cols = cols, scale = NULL)
    }, error = function(e) NULL)
    
    if (!is.null(coords) && all(cols %in% colnames(coords))) {
      colnames(coords)[match(cols, colnames(coords))] <- c("row", "col")
      return(coords)
    }
  }
  
  stop("❌ 无法获取坐标信息，请检查 Seurat 对象")
}

# -----------------------------
# 安全文件名
# -----------------------------
safe_filename <- function(name) {
  gsub("[^[:alnum:]]", "_", name)
}

# -----------------------------
# 计算坐标范围
# -----------------------------
calculate_coord_limits <- function(plot_data, expand = 0.05) {
  col_range <- range(plot_data$col, na.rm = TRUE)
  row_range <- range(plot_data$row, na.rm = TRUE)
  
  col_expand <- diff(col_range) * expand
  row_expand <- diff(row_range) * expand
  
  list(
    col = c(col_range[1] - col_expand, col_range[2] + col_expand),
    row = c(row_range[1] - row_expand, row_range[2] + row_expand)
  )
}

# ===================================================================
# 文件过滤函数
# ===================================================================

filter_seurat_files <- function(files, config) {
  original_count <- length(files)
  
  # 只保留指定的文件
  if (!is.null(config$specific_files)) {
    basenames <- basename(files)
    files <- files[basenames %in% config$specific_files]
    
    if (length(files) == 0) {
      stop("❌ 未找到任何指定的文件")
    }
    
    # 检查是否有未找到的文件
    found_files <- basename(files)
    missing <- setdiff(config$specific_files, found_files)
    if (length(missing) > 0) {
      warning(sprintf("⚠️  以下指定文件未找到:\n  %s", 
                      paste(missing, collapse = "\n  ")))
    }
    
    cat(sprintf("✓ 匹配到 %d/%d 个指定文件\n", 
                length(files), length(config$specific_files)))
  }
  
  # 排除指定的文件
  if (!is.null(config$exclude_files)) {
    basenames <- basename(files)
    excluded_count <- sum(basenames %in% config$exclude_files)
    files <- files[!basenames %in% config$exclude_files]
    
    if (length(files) == 0) {
      stop("❌ 过滤后没有剩余文件")
    }
    
    if (excluded_count > 0) {
      cat(sprintf("✓ 排除了 %d 个文件\n", excluded_count))
    }
  }
  
  if (length(files) != original_count) {
    cat(sprintf("📋 文件过滤: %d -> %d\n", original_count, length(files)))
  }
  
  return(files)
}

# ========== 新增函数：坐标标准化 ==========
standardize_spatial_coordinates <- function(seurat_obj) {
  # 检查是否是 Seurat 对象
  if (!inherits(seurat_obj, "Seurat")) {
    stop("输入必须是 Seurat 对象")
  }
  
  # 获取所有图像名称
  image_names <- names(seurat_obj@images)
  
  if (length(image_names) == 0) {
    warning("未找到空间图像数据")
    return(seurat_obj)
  }
  
  # 定义可能的坐标列名
  possible_row_names <- c("row", "imagerow", "array_row", "tissue_row", "pxl_row_in_fullres")
  possible_col_names <- c("col", "imagecol", "array_col", "tissue_col", "pxl_col_in_fullres")
  
  cat(sprintf(">> 检查 %d 个图像的坐标系统...\n", length(image_names)))
  
  coord_issues <- 0
  
  for (img_name in image_names) {
    img_obj <- seurat_obj@images[[img_name]]
    
    # 检查是否有 coordinates 槽
    if (!"coordinates" %in% slotNames(img_obj)) {
      warning(sprintf("图像 '%s' 没有 coordinates 槽", img_name))
      coord_issues <- coord_issues + 1
      next
    }
    
    coords <- img_obj@coordinates
    coord_cols <- colnames(coords)
    
    # 检查是否已经有标准的 row/col 列
    has_row <- "row" %in% coord_cols
    has_col <- "col" %in% coord_cols
    
    if (has_row && has_col) {
      # 已经有标准列名，跳过
      next
    }
    
    # 查找实际的行列名
    actual_row_name <- intersect(coord_cols, possible_row_names)[1]
    actual_col_name <- intersect(coord_cols, possible_col_names)[1]
    
    if (is.na(actual_row_name) || is.na(actual_col_name)) {
      warning(sprintf(
        "图像 '%s' 未找到有效的坐标列。\n可用列: %s", 
        img_name, 
        paste(coord_cols, collapse=", ")
      ))
      coord_issues <- coord_issues + 1
      next
    }
    
    # 标准化列名
    if (!has_row) {
      coords$row <- coords[[actual_row_name]]
      cat(sprintf("   %s: %s → row\n", img_name, actual_row_name))
    }
    
    if (!has_col) {
      coords$col <- coords[[actual_col_name]]
      cat(sprintf("   %s: %s → col\n", img_name, actual_col_name))
    }
    
    # 验证坐标值
    if (any(is.na(coords$row)) || any(is.na(coords$col))) {
      warning(sprintf("图像 '%s' 包含 NA 坐标值", img_name))
      coord_issues <- coord_issues + 1
    }
    
    # 更新坐标
    seurat_obj@images[[img_name]]@coordinates <- coords
  }
  
  if (coord_issues > 0) {
    warning(sprintf("发现 %d 个图像存在坐标问题", coord_issues))
  }
  
  return(seurat_obj)
}


```

---

### 03_load_data.R

- **大小**: 5.76 KB
- **修改时间**: 2025-11-06 22:00:52

```r
#!/usr/bin/env Rscript
# ===================================================================
# 数据加载模块
# ===================================================================

load_gene_list <- function(config) {
  cat("📄 读取基因列表...\n")
  
  # ✅ 确保目录存在
  if (!is.null(config$cache_dir)) {
    if (!dir.exists(config$cache_dir)) {
      dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # 生成缓存文件路径
    cache_file <- file.path(config$cache_dir, "gene_list.rds")
  } else {
    # 如果没有缓存目录，设置为 NULL（不使用缓存）
    cache_file <- NULL
    cat("⚠️  未设置缓存目录，将不使用缓存\n")
  }
  
  # 检查缓存
  if (!is.null(cache_file) && is_cache_valid(cache_file, config$cache_max_age_hours)) {
    cat("✓ 从缓存加载基因列表\n")
    gene_list <- readRDS(cache_file)
    cat(sprintf("✓ 加载了 %d 个基因\n", length(gene_list)))
    return(gene_list)
  }
  
  # 读取基因列表
  if (!file.exists(config$gene_list_path)) {
    stop(sprintf("❌ 基因列表文件不存在: %s", config$gene_list_path))
  }
  
  gene_list <- readLines(config$gene_list_path)
  gene_list <- gene_list[gene_list != ""]  # 移除空行
  
  cat(sprintf("✓ 从文件读取了 %d 个基因\n", length(gene_list)))
  
  # 保存缓存
  if (!is.null(cache_file)) {
    tryCatch({
      saveRDS(gene_list, cache_file)
      cat(sprintf("✓ 缓存已保存至: %s\n", cache_file))
    }, error = function(e) {
      warning(sprintf("⚠️  无法保存缓存: %s", e$message))
    })
  }
  
  return(gene_list)
}


load_seurat_object <- function(config) {
  cat("📥 加载 Seurat 对象...\n")
  
  # 确保目录存在
  if (!is.null(config$cache_dir)) {
    if (!dir.exists(config$cache_dir)) {
      dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # 生成缓存文件名
    seurat_basename <- tools::file_path_sans_ext(basename(config$seurat_path))
    cache_file <- file.path(config$cache_dir, sprintf("%s_processed.rds", seurat_basename))
  } else {
    cache_file <- NULL
    cat("⚠️  未设置缓存目录，将不使用缓存\n")
  }
  
  # 检查缓存
  if (!is.null(cache_file) && is_cache_valid(cache_file, config$cache_max_age_hours)) {
    cat("✓ 从缓存加载 Seurat 对象\n")
    seurat_obj <- readRDS(cache_file)
    cat(sprintf("✓ 加载完成: %d 个细胞\n", ncol(seurat_obj)))
    
    # ✅ 修复缓存的对象
    seurat_obj <- fix_seurat_object(seurat_obj)
    
    return(seurat_obj)
  }
  
  # 加载原始 Seurat 对象
  if (!file.exists(config$seurat_path)) {
    stop(sprintf("❌ Seurat 文件不存在: %s", config$seurat_path))
  }
  
  cat(sprintf("📂 文件: %s\n", basename(config$seurat_path)))
  cat(sprintf("📏 大小: %.2f GB\n", file.size(config$seurat_path) / (1024^3)))
  
  seurat_obj <- readRDS(config$seurat_path)
  
  cat(sprintf("✓ 加载完成: %d 个细胞, %d 个基因\n", 
              ncol(seurat_obj), nrow(seurat_obj)))
  
  # ✅ 修复对象
  seurat_obj <- fix_seurat_object(seurat_obj)
  
  # 保存缓存（可选）
  if (!is.null(cache_file) && config$save_full_object) {
    tryCatch({
      saveRDS(seurat_obj, cache_file)
      cat(sprintf("✓ 缓存已保存至: %s\n", cache_file))
    }, error = function(e) {
      warning(sprintf("⚠️  无法保存缓存: %s", e$message))
    })
  }
  
  return(seurat_obj)
}

check_gene_overlap <- function(gene_list, seurat_obj) {
  cat("🔍 检查基因匹配情况...\n")
  
  genes_in_data <- intersect(gene_list, rownames(seurat_obj))
  genes_missing <- setdiff(gene_list, rownames(seurat_obj))
  
  cat(sprintf("✅ 匹配基因: %d / %d (%.1f%%)\n",
              length(genes_in_data), length(gene_list),
              100 * length(genes_in_data) / length(gene_list)))
  
  if (length(genes_missing) > 0) {
    n_show <- min(10, length(genes_missing))
    cat(sprintf("⚠️ 缺失 %d 个基因 (前%d个): %s%s\n", 
                length(genes_missing), n_show,
                paste(head(genes_missing, n_show), collapse = ", "),
                ifelse(length(genes_missing) > n_show, " ...", "")))
  }
  
  cat("\n")
  return(genes_in_data)
}

fix_seurat_object <- function(seurat_obj) {
  cat("🔧 检查并修复 Seurat 对象...\n")
  
  # 检查是否是 Seurat 对象
  if (!inherits(seurat_obj, "Seurat")) {
    warning("⚠️  对象不是 Seurat 类")
    return(seurat_obj)
  }
  
  # 修复 VisiumV1 对象
  if (length(seurat_obj@images) > 0) {
    for (img_name in names(seurat_obj@images)) {
      img <- seurat_obj@images[[img_name]]
      
      # 检查是否是 VisiumV1 类
      if (inherits(img, "VisiumV1")) {
        cat(sprintf("   🔧 修复图像: %s\n", img_name))
        
        # 添加缺失的 misc 插槽
        if (!.hasSlot(img, "misc")) {
          tryCatch({
            img@misc <- list()
            seurat_obj@images[[img_name]] <- img
            cat(sprintf("   ✓ 已添加 misc 插槽\n"))
          }, error = function(e) {
            cat(sprintf("   ⚠️  无法添加 misc 插槽: %s\n", e$message))
          })
        }
        
        # 验证修复
        tryCatch({
          validObject(img)
          cat(sprintf("   ✓ 对象验证通过\n"))
        }, error = function(e) {
          cat(sprintf("   ⚠️  对象验证失败: %s\n", e$message))
          
          # 尝试更激进的修复：移除 images
          cat(sprintf("   🔧 尝试移除有问题的图像对象...\n"))
          seurat_obj@images[[img_name]] <- NULL
        })
      }
    }
  }
  
  cat("✓ Seurat 对象检查完成\n\n")
  return(seurat_obj)
}
```

---

### 04_module_score.R

- **大小**: 1.81 KB
- **修改时间**: 2025-11-05 17:09:13

```r
#!/usr/bin/env Rscript
# ===================================================================
# Module Score 计算
# ===================================================================

calculate_module_score <- function(seurat_obj, genes, config) {
  cat("🧮 计算 Clock Gene Module Score...\n")
  
  # 生成缓存key
  cache_key <- generate_cache_key(genes, dim(seurat_obj), "AddModuleScore")
  cache_file <- file.path(config$cache_dir, sprintf("module_score_%s.rds", cache_key))
  
  if (file.exists(cache_file)) {
    score_data <- load_cache(cache_file, "Module Score")
    seurat_obj$ClockGene_Score1 <- score_data$ClockGene_Score1
  } else {
    cat("🔄 正在计算 Module Score...\n")
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(clock_gene_set = genes),
      name = "ClockGene_Score"
    )
    score_data <- data.frame(ClockGene_Score1 = seurat_obj$ClockGene_Score1)
    save_cache(score_data, cache_file, "Module Score")
  }
  
  cat(sprintf("✅ 评分范围: %.3f ~ %.3f\n\n", 
              min(seurat_obj$ClockGene_Score1, na.rm = TRUE),
              max(seurat_obj$ClockGene_Score1, na.rm = TRUE)))
  
  return(seurat_obj)
}

define_high_expression <- function(seurat_obj, config) {
  cat("🎯 设置高表达阈值...\n")
  
  threshold <- quantile(seurat_obj$ClockGene_Score1, 
                       config$threshold_quantile, 
                       na.rm = TRUE)
  
  threshold_pct <- (1 - config$threshold_quantile) * 100
  cat(sprintf("✅ 高表达阈值: %.3f (Top %.1f%%)\n", threshold, threshold_pct))
  
  seurat_obj$ClockGene_High <- seurat_obj$ClockGene_Score1 > threshold
  
  cat("✅ 高/低表达分组:\n")
  print(table(seurat_obj$ClockGene_High))
  cat("\n")
  
  return(list(seurat_obj = seurat_obj, threshold = threshold))
}
```

---

### 05_niche_analysis.R

- **大小**: 10.56 KB
- **修改时间**: 2025-11-07 12:35:16

```r
# 05_niche_analysis.R

#!/usr/bin/env Rscript
# ===================================================================
# Niche 距离计算
# ===================================================================

perform_niche_analysis <- function(seurat_obj, threshold, config) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   Niche 距离分析\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 生成缓存key
  cache_key <- generate_cache_key(
    threshold, 
    sum(seurat_obj$ClockGene_High, na.rm = TRUE), 
    ncol(seurat_obj), 
    config$niche_dist_method
  )
  cache_file <- file.path(config$cache_dir, sprintf("niche_analysis_%s.rds", cache_key))
  
  # 检查缓存
  if (file.exists(cache_file)) {
    cat("📦 从缓存加载 Niche 距离数据...\n")
    niche_data <- load_cache(cache_file, "Niche 距离")
    
    # 验证缓存数据
    if (length(niche_data$ClockGene_Distance) != ncol(seurat_obj)) {
      warning("⚠️ 缓存数据大小不匹配，将重新计算")
      file.remove(cache_file)
    } else {
      seurat_obj$ClockGene_Distance <- niche_data$ClockGene_Distance
      
      cat(sprintf("✅ 距离范围: %.2f ~ %.2f\n",
                  min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
                  max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))
      cat("✅ Niche 分析完成（从缓存加载）\n\n")
      
      return(seurat_obj)
    }
  }
  
  # 如果没有缓存或缓存无效，进行计算
  cat("🔄 开始计算 Niche 距离...\n\n")
  
  # 数据统计
  n_total <- ncol(seurat_obj)
  n_marker <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
  cat(sprintf("数据概况:\n"))
  cat(sprintf("  总细胞数: %d\n", n_total))
  cat(sprintf("  标记细胞: %d (%.1f%%)\n", n_marker, 100 * n_marker / n_total))
  cat(sprintf("  使用核心数: %d\n", config$n_workers))
  cat(sprintf("  距离方法: %s\n\n", config$niche_dist_method))
  
  # 验证必需的列
  if (!"ClockGene_High" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Seurat 对象中缺少 'ClockGene_High' 列，请先运行 define_high_expression()")
  }
  
  if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Seurat 对象中缺少 'orig.ident' 列")
  }
  
  # 验证空间数据
  if (length(names(seurat_obj@images)) == 0) {
    stop("❌ Seurat 对象中没有空间图像数据")
  }
  
  # 调用 niche_marker 函数
  result <- tryCatch({
    
    niche_marker(
      .data = seurat_obj,
      marker = ClockGene_High,
      spot_type = ClockGene_Distance,
      slide = orig.ident,
      dist_method = config$niche_dist_method,
      FUN = NA,
      n_work = config$n_workers
    )
    
  }, error = function(e) {
    # 详细的错误诊断
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("   ❌ Niche 分析失败\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    cat(sprintf("错误信息: %s\n\n", e$message))
    
    # 诊断信息
    cat("🔍 诊断信息:\n\n")
    
    # 样本信息
    sample_names <- unique(seurat_obj$orig.ident)
    cat(sprintf("1. 样本数量: %d\n", length(sample_names)))
    cat(sprintf("   样本列表（前5个）: %s\n\n", 
                paste(head(sample_names, 5), collapse=", ")))
    
    # 空间数据信息
    image_names <- names(seurat_obj@images)
    cat(sprintf("2. 空间图像数: %d\n", length(image_names)))
    
    if (length(image_names) > 0) {
      cat("   图像列表（前5个）:\n")
      for (img in head(image_names, 5)) {
        cat(sprintf("     - %s\n", img))
        if (img %in% names(seurat_obj@images)) {
          coords <- seurat_obj@images[[img]]@coordinates
          cat(sprintf("       细胞数: %d\n", nrow(coords)))
          cat(sprintf("       坐标列: %s\n", paste(colnames(coords), collapse=", ")))
        }
      }
    }
    cat("\n")
    
    # 标记细胞信息
    cat(sprintf("3. 标记细胞统计:\n"))
    marker_table <- table(seurat_obj$ClockGene_High)
    print(marker_table)
    cat("\n")
    
    # 按样本统计标记细胞
    cat("4. 各样本标记细胞数:\n")
    marker_by_sample <- table(
      seurat_obj$orig.ident[seurat_obj$ClockGene_High]
    )
    print(head(marker_by_sample, 10))
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    # 抛出原始错误
    stop(sprintf("Niche 分析失败: %s", e$message))
  })
  
  seurat_obj <- result
  
  # 验证结果
  if (!"ClockGene_Distance" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Niche 分析未能生成 'ClockGene_Distance' 列")
  }
  
  if (any(is.na(seurat_obj$ClockGene_Distance))) {
    n_na <- sum(is.na(seurat_obj$ClockGene_Distance))
    warning(sprintf("⚠️ 警告：%d 个细胞的距离值为 NA", n_na))
  }
  
  # 保存缓存
  cat("\n💾 保存结果到缓存...\n")
  niche_data <- data.frame(
    ClockGene_Distance = seurat_obj$ClockGene_Distance,
    stringsAsFactors = FALSE
  )
  save_cache(niche_data, cache_file, "Niche 距离")
  
  # 输出结果摘要
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   结果摘要\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  dist_vals <- seurat_obj$ClockGene_Distance
  cat(sprintf("Distance 统计:\n"))
  cat(sprintf("  最小值: %.2f\n", min(dist_vals, na.rm = TRUE)))
  cat(sprintf("  第25百分位: %.2f\n", quantile(dist_vals, 0.25, na.rm = TRUE)))
  cat(sprintf("  中位数: %.2f\n", median(dist_vals, na.rm = TRUE)))
  cat(sprintf("  第75百分位: %.2f\n", quantile(dist_vals, 0.75, na.rm = TRUE)))
  cat(sprintf("  最大值: %.2f\n", max(dist_vals, na.rm = TRUE)))
  cat(sprintf("  平均值: %.2f\n", mean(dist_vals, na.rm = TRUE)))
  cat(sprintf("  标准差: %.2f\n", sd(dist_vals, na.rm = TRUE)))
  
  # 标记细胞的距离验证
  marker_dist <- dist_vals[seurat_obj$ClockGene_High]
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- sum(!is.na(marker_dist))
  
  cat(sprintf("\n标记细胞验证:\n"))
  cat(sprintf("  标记细胞数: %d\n", n_marker_total))
  cat(sprintf("  Distance=0: %d (%.1f%%)\n", 
              n_marker_zero, 
              100 * n_marker_zero / n_marker_total))
  
  if (n_marker_zero / n_marker_total < 0.95) {
    cat("\n⚠️ 警告：少于95%的标记细胞距离为0，可能存在问题\n")
  } else {
    cat("\n✅ 验证通过：标记细胞距离计算正确\n")
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("   Niche 分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  return(seurat_obj)
}


# ===================================================================
# 辅助函数：快速诊断（可选）
# ===================================================================

quick_diagnose_niche <- function(seurat_obj) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   快速 Niche 诊断\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 基本信息
  cat("1. 基本信息:\n")
  cat(sprintf("   总细胞数: %d\n", ncol(seurat_obj)))
  cat(sprintf("   总基因数: %d\n", nrow(seurat_obj)))
  
  # 检查必需的列
  cat("\n2. 必需列检查:\n")
  required_cols <- c("ClockGene_High", "orig.ident")
  for (col in required_cols) {
    exists <- col %in% colnames(seurat_obj@meta.data)
    cat(sprintf("   %s: %s\n", col, ifelse(exists, "✓", "✗")))
  }
  
  # 空间数据
  cat("\n3. 空间数据:\n")
  n_images <- length(names(seurat_obj@images))
  cat(sprintf("   图像数: %d\n", n_images))
  
  if (n_images > 0) {
    cat("   样本列表:\n")
    for (img in names(seurat_obj@images)) {
      coords <- seurat_obj@images[[img]]@coordinates
      cat(sprintf("     - %s: %d 个细胞, 坐标列: [%s]\n", 
                  img, 
                  nrow(coords),
                  paste(colnames(coords), collapse=", ")))
    }
  }
  
  # 标记细胞
  if ("ClockGene_High" %in% colnames(seurat_obj@meta.data)) {
    cat("\n4. 标记细胞:\n")
    n_marked <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
    cat(sprintf("   总数: %d (%.1f%%)\n", 
                n_marked, 
                100 * n_marked / ncol(seurat_obj)))
    
    # 按样本统计
    if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      cat("\n   各样本标记细胞数:\n")
      marker_by_sample <- table(
        seurat_obj$orig.ident[seurat_obj$ClockGene_High]
      )
      for (i in seq_along(marker_by_sample)) {
        sample_name <- names(marker_by_sample)[i]
        n_marker <- marker_by_sample[i]
        n_total <- sum(seurat_obj$orig.ident == sample_name)
        cat(sprintf("     - %s: %d/%d (%.1f%%)\n", 
                    sample_name, n_marker, n_total,
                    100 * n_marker / n_total))
      }
    }
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(NULL)
}
```

---

### 06_plot_isoheight.R

- **大小**: 10.14 KB
- **修改时间**: 2025-11-07 13:59:07

```r
#!/usr/bin/env Rscript
# ===================================================================
# 等高线密度图绘制模块（优化版 + 进度条）
# 功能：多线程并行绘制 Clock Gene 等高线密度图
# ===================================================================

library(future)
library(future.apply)
library(progressr)
library(ggplot2)


#' 绘制等高线密度图（接收预切分样本）
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param col_bg 背景点颜色
#' @param col_top 高表达点颜色
#' @param col_isoheight 等高线颜色
#' @param col_white_ratio 白色占比
#' @param cols_fill_isoheight 填充色谱
#' @param plot_width 图宽
#' @param plot_height 图高
#' @param nrow 子图排列行数
#' 
#' @return 处理结果列表
#'
plot_isoheight <- function(sample_list, 
                          CONFIG,
                          col_bg = "gray92",
                          col_top = "#d62728",
                          col_isoheight = "white",
                          col_white_ratio = 0.25,
                          cols_fill_isoheight = NULL,
                          plot_width = 8,
                          plot_height = 8,
                          nrow = 1) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   等高线密度图绘制（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  if (is.null(CONFIG$dirs$isoheight)) {
    stop("❌ CONFIG$dirs$isoheight 未定义")
  }
  
  if (!exists("celltype_isoheight_plot")) {
    stop("❌ 未找到 celltype_isoheight_plot 函数，请先加载 SSS_isoheight_plot.R")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$isoheight)) {
    dir.create(CONFIG$dirs$isoheight, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 提取参数
  size_bg <- CONFIG$plot$point_size_bg %||% 0.3
  size_top <- CONFIG$plot$point_size_top %||% 1.2
  dpi <- CONFIG$plot$dpi %||% 300
  n_workers <- CONFIG$n_workers %||% 4
  
  # 默认色谱
  if (is.null(cols_fill_isoheight)) {
    cols_fill_isoheight <- c(
      rep("white", 25),
      colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])(75)
    )
  }
  
  cat(sprintf("📊 将绘制 %d 个样本\n", length(sample_list)))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 2. 设置并行环境
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 3. 设置进度条
  # ========================================
  
  # 检查是否已经设置了 handlers
  has_handlers <- !is.null(progressr::handlers(NULL))
  
  if (!has_handlers) {
    # 如果没有设置，使用详细的进度条
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "[:bar] :percent | 已完成: :current/:total | 预计剩余: :eta | :message",
        width    = 80,
        complete = "=",
        clear    = FALSE
      )
    ))
  }
  
  # ========================================
  # 4. 并行绘图
  # ========================================
  
  cat("🎨 开始绘图...\n\n")
  
  progressr::with_progress({
    
    p <- progressr::progressor(
      steps = length(sample_list),
      message = "绘制等高线图"
    )
    
    results <- future.apply::future_lapply(
      
      names(sample_list), 
      
      function(sample_id) {
        
        tryCatch({
          
          # 获取样本数据
          seurat_subset <- sample_list[[sample_id]]
          
          # 验证数据
          if (ncol(seurat_subset) == 0) {
            p(message = sprintf("⚠️  %s - 无数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No data"
            ))
          }
          
          if (!"ClockGene_High" %in% colnames(seurat_subset@meta.data)) {
            p(message = sprintf("⚠️  %s - 缺少 ClockGene_High 列", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "Missing ClockGene_High column"
            ))
          }
          
          n_high <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
          
          if (n_high == 0) {
            p(message = sprintf("⚠️  %s - 无高表达点", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No high expression spots"
            ))
          }
          
          # 绘图
          p_iso <- celltype_isoheight_plot(
            .data = seurat_subset,
            density_top = ClockGene_High,
            col_bg = col_bg,
            col_top = col_top,
            col_isoheight = col_isoheight,
            col_white_ratio = col_white_ratio,
            cols_fill_isoheight = cols_fill_isoheight,
            size_bg = size_bg,
            size_top = size_top,
            nrow = nrow
          )
          
          # 保存
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          output_file <- sprintf("ClockGene_isoheight_%s.pdf", safe_name)
          output_path <- file.path(CONFIG$dirs$isoheight, output_file)
          
          ggplot2::ggsave(
            filename = output_path,
            plot = p_iso, 
            width = plot_width, 
            height = plot_height, 
            dpi = dpi
          )
          
          # 统计
          file_size_mb <- file.size(output_path) / 1024^2
          n_spots <- ncol(seurat_subset)
          high_pct <- 100 * n_high / n_spots
          
          # 更新进度（显示样本名）
          p(message = sprintf("✅ %s (%.2f MB)", sample_id, file_size_mb))
          
          return(list(
            sample = sample_id,
            success = TRUE,
            file = output_path,
            file_size_mb = file_size_mb,
            n_spots = n_spots,
            n_high = n_high,
            high_pct = high_pct
          ))
          
        }, error = function(e) {
          p(message = sprintf("❌ %s - %s", sample_id, e$message))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = as.character(e$message)
          ))
        })
      },
      
      future.seed = TRUE,
      future.chunk.size = 1
    )
  })
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  future::plan(future::sequential)
  
  # ========================================
  # 5. 统计输出
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, 
              length(sample_list),
              100 * n_success / length(sample_list)))
  
  if (n_failed > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", n_failed, length(sample_list)))
    cat("失败样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  • %s: %s\n", res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  if (n_success > 0) {
    cat("成功样本:\n")
    cat(sprintf("%-30s %10s %10s %10s %10s\n", 
                "样本", "Spots", "High", "High%", "文件大小"))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    total_file_size <- 0
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %10d %10d %9.1f%% %8.2f MB\n",
                    res$sample,
                    res$n_spots,
                    res$n_high,
                    res$high_pct,
                    res$file_size_mb))
        total_file_size <- total_file_size + res$file_size_mb
      }
    }
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("%-30s %10s %10s %10s %8.2f MB\n",
                "总计", "", "", "", total_file_size))
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$isoheight))
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 6. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    output_dir = CONFIG$dirs$isoheight,
    elapsed_time = as.numeric(elapsed),
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 06_plot_isoheight.R 已加载\n")
```

---

### 07_plot_spatial.R

- **大小**: 11.8 KB
- **修改时间**: 2025-11-07 14:00:36

```r
#!/usr/bin/env Rscript
# ===================================================================
# 空间梯度图绘制模块（优化版 + 进度条）
# 功能：绘制 Clock Gene 距离场的空间分布图
# ===================================================================

library(future)
library(future.apply)
library(progressr)
library(ggplot2)


#' 绘制空间梯度图
#'
#' @param sample_list 预切分的样本列表
#' @param CONFIG 配置对象
#' @param plot_width 图宽（英寸）
#' @param plot_height 图高（英寸）
#' @param pt_size_factor 点大小因子
#' @param alpha 透明度范围
#' @param color_option viridis 色板选项
#' @param color_direction 色板方向
#' 
#' @return 处理结果列表
#'
plot_spatial_gradient <- function(sample_list, 
                                  CONFIG,
                                  plot_width = 10,
                                  plot_height = 8,
                                  pt_size_factor = 1.5,
                                  alpha = c(0.3, 1),
                                  color_option = "magma",
                                  color_direction = -1) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   空间梯度图绘制（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  if (is.null(CONFIG$dirs$spatial)) {
    stop("❌ CONFIG$dirs$spatial 未定义")
  }
  
  # 创建输出目录
  if (!dir.exists(CONFIG$dirs$spatial)) {
    dir.create(CONFIG$dirs$spatial, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 提取参数
  n_workers <- CONFIG$n_workers %||% 4
  dpi <- CONFIG$plot$dpi %||% 300
  
  cat(sprintf("📊 将绘制 %d 个样本\n", length(sample_list)))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 2. 设置并行环境
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  start_time <- Sys.time()
  
  # ========================================
  # 3. 设置进度条
  # ========================================
  
  # 检查是否已经设置了 handlers
  has_handlers <- !is.null(progressr::handlers(NULL))
  
  if (!has_handlers) {
    # 如果没有设置，使用详细的进度条
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "[:bar] :percent | 已完成: :current/:total | 预计剩余: :eta | :message",
        width    = 80,
        complete = "=",
        clear    = FALSE
      )
    ))
  }
  
  # ========================================
  # 4. 并行绘图
  # ========================================
  
  cat("🗺️  开始绘图...\n\n")
  
  progressr::with_progress({
    
    p <- progressr::progressor(
      steps = length(sample_list),
      message = "绘制空间梯度图"
    )
    
    results <- future.apply::future_lapply(
      
      names(sample_list),
      
      function(sample_id) {
        
        tryCatch({
          
          # 获取样本数据
          seurat_subset <- sample_list[[sample_id]]
          
          # 验证数据
          if (ncol(seurat_subset) == 0) {
            p(message = sprintf("⚠️  %s - 无数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No data"
            ))
          }
          
          if (!"ClockGene_Distance" %in% colnames(seurat_subset@meta.data)) {
            p(message = sprintf("⚠️  %s - 缺少距离数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "Missing ClockGene_Distance column"
            ))
          }
          
          # 检查空间数据
          if (length(Seurat::Images(seurat_subset)) == 0) {
            p(message = sprintf("⚠️  %s - 无空间图像数据", sample_id))
            return(list(
              sample = sample_id,
              success = FALSE,
              error = "No spatial image data"
            ))
          }
          
          # 统计距离数据
          distance_values <- seurat_subset$ClockGene_Distance
          distance_stats <- list(
            min = min(distance_values, na.rm = TRUE),
            max = max(distance_values, na.rm = TRUE),
            mean = mean(distance_values, na.rm = TRUE),
            median = median(distance_values, na.rm = TRUE),
            sd = sd(distance_values, na.rm = TRUE)
          )
          
          # 绘制空间分布图
          p_spatial <- Seurat::SpatialFeaturePlot(
            seurat_subset,
            features = "ClockGene_Distance",
            pt.size.factor = pt_size_factor,
            alpha = alpha,
            stroke = 0
          ) + 
            ggplot2::scale_fill_viridis_c(
              option = color_option,
              direction = color_direction,
              name = "Distance\nto High",
              limits = c(0, NA)
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              legend.position = "right",
              legend.title = ggplot2::element_text(size = 12, face = "bold"),
              legend.text = ggplot2::element_text(size = 10),
              plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5)
            ) +
            ggplot2::ggtitle(
              sprintf("Clock Gene Distance Field - %s", sample_id),
              subtitle = sprintf(
                "Mean: %.2f | Median: %.2f | Range: [%.2f, %.2f]",
                distance_stats$mean,
                distance_stats$median,
                distance_stats$min,
                distance_stats$max
              )
            )
          
          # 保存图形
          safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
          output_file <- sprintf("ClockGene_spatial_%s.pdf", safe_name)
          output_path <- file.path(CONFIG$dirs$spatial, output_file)
          
          ggplot2::ggsave(
            filename = output_path,
            plot = p_spatial,
            width = plot_width,
            height = plot_height,
            dpi = dpi
          )
          
          # 统计信息
          file_size_mb <- file.size(output_path) / 1024^2
          n_spots <- ncol(seurat_subset)
          
          # 更新进度（显示样本名和文件大小）
          p(message = sprintf("✅ %s (%.2f MB)", sample_id, file_size_mb))
          
          return(list(
            sample = sample_id,
            success = TRUE,
            file = output_path,
            file_size_mb = file_size_mb,
            n_spots = n_spots,
            distance_stats = distance_stats
          ))
          
        }, error = function(e) {
          p(message = sprintf("❌ %s - %s", sample_id, e$message))
          return(list(
            sample = sample_id,
            success = FALSE,
            error = as.character(e$message)
          ))
        })
      },
      
      future.seed = TRUE,
      future.chunk.size = 1
    )
  })
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  future::plan(future::sequential)
  
  # ========================================
  # 5. 统计输出
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   绘图完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, 
              length(sample_list),
              100 * n_success / length(sample_list)))
  
  if (n_failed > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", n_failed, length(sample_list)))
    cat("失败样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  • %s: %s\n", res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  if (n_success > 0) {
    cat("成功样本:\n")
    cat(sprintf("%-30s %10s %12s %12s %12s %10s\n", 
                "样本", "Spots", "平均距离", "中位距离", "标准差", "文件大小"))
    cat(paste(rep("-", 95), collapse = ""), "\n")
    
    total_file_size <- 0
    total_spots <- 0
    all_means <- c()
    all_medians <- c()
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %10d %12.2f %12.2f %12.2f %8.2f MB\n",
                    res$sample,
                    res$n_spots,
                    res$distance_stats$mean,
                    res$distance_stats$median,
                    res$distance_stats$sd,
                    res$file_size_mb))
        
        total_file_size <- total_file_size + res$file_size_mb
        total_spots <- total_spots + res$n_spots
        all_means <- c(all_means, res$distance_stats$mean)
        all_medians <- c(all_medians, res$distance_stats$median)
      }
    }
    
    cat(paste(rep("-", 95), collapse = ""), "\n")
    
    if (n_success > 1) {
      cat(sprintf("%-30s %10d %12.2f %12.2f %12s %8.2f MB\n",
                  "平均",
                  as.integer(total_spots / n_success),
                  mean(all_means),
                  mean(all_medians),
                  "-",
                  total_file_size / n_success))
      cat(sprintf("%-30s %10d %12s %12s %12s %8.2f MB\n",
                  "总计",
                  total_spots,
                  "-",
                  "-",
                  "-",
                  total_file_size))
    } else {
      cat(sprintf("%-30s %10s %12s %12s %12s %8.2f MB\n",
                  "总计",
                  "",
                  "",
                  "",
                  "",
                  total_file_size))
    }
    
    cat("\n")
  }
  
  cat(sprintf("⏱️  总耗时: %.2f 秒 (平均 %.2f 秒/样本)\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  cat(sprintf("📁 输出目录: %s\n", CONFIG$dirs$spatial))
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 6. 返回结果
  # ========================================
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    output_dir = CONFIG$dirs$spatial,
    elapsed_time = as.numeric(elapsed),
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 07_plot_spatial.R 已加载\n")
```

---

### 08_plot_celltype.R

- **大小**: 6.56 KB
- **修改时间**: 2025-11-07 14:03:13

```r
#!/usr/bin/env Rscript
# ===================================================================
# 细胞类型 Niche 分析模块（简化版 + 进度条）
# 功能：分析不同密度区域的细胞类型分布和富集
# ===================================================================

library(future)
library(future.apply)
library(progressr)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

# ===================================================================
# 加载工具函数
# ===================================================================

utils_dir <- "08_plot_celltype_utils"

source(file.path(utils_dir, "00_operators.R"))
source(file.path(utils_dir, "01_color_schemes.R"))
source(file.path(utils_dir, "02_density_zones.R"))
source(file.path(utils_dir, "03_plot_overlay.R"))
source(file.path(utils_dir, "04_plot_composition.R"))
source(file.path(utils_dir, "05_plot_heatmap.R"))
source(file.path(utils_dir, "06_plot_combined.R"))
source(file.path(utils_dir, "07_statistics.R"))
source(file.path(utils_dir, "08_validation.R"))
source(file.path(utils_dir, "09_save_plots.R"))
source(file.path(utils_dir, "10_summary.R"))

cat("✅ 已加载所有工具函数\n")


# ===================================================================
# 主函数：细胞类型 Niche 分析
# ===================================================================

#' 细胞类型 Niche 分析
#'
#' @param sample_list 预切分的样本列表（来自 main.R）
#' @param CONFIG 配置对象
#' @param density_bins 密度分区数量
#' @param celltype_col 细胞类型列名
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
#' @param seurat_basename 文件基础名
#' 
#' @return 处理结果列表
#'
analyze_celltype_niche <- function(
    sample_list,
    CONFIG,
    density_bins = 10,
    celltype_col = "predicted.id",
    plot_overlay = TRUE,
    plot_composition = TRUE,
    plot_heatmap = TRUE,
    plot_combined = TRUE,
    seurat_basename = NULL
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   细胞类型 Niche 分析（多线程并行）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # ========================================
  # 1. 参数验证
  # ========================================
  
  validate_inputs(sample_list, CONFIG)
  validate_required_functions()
  
  # ========================================
  # 2. 初始化配置
  # ========================================
  
  setup_colors(sample_list[[1]], CONFIG, celltype_col, density_bins)
  
  n_workers <- CONFIG$n_workers %||% 4
  
  cat(sprintf("📊 将分析 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 密度分区: %d 个区域 (Zone_0=核心, Zone_%d=外围)\n", 
              density_bins, density_bins - 1))
  cat(sprintf("🔧 使用 %d 个线程\n\n", n_workers))
  
  # ========================================
  # 3. 设置并行和进度条
  # ========================================
  
  future::plan(future::multisession, workers = n_workers)
  options(future.globals.maxSize = Inf)
  
  # 设置进度条
  has_handlers <- !is.null(progressr::handlers(NULL))
  
  if (!has_handlers) {
    progressr::handlers(list(
      progressr::handler_progress(
        format   = "[:bar] :percent | 已完成: :current/:total | 预计剩余: :eta | :message",
        width    = 80,
        complete = "=",
        clear    = FALSE
      )
    ))
  }
  
  start_time <- Sys.time()
  
  # ========================================
  # 4. 并行处理样本
  # ========================================
  
  cat("🔬 开始分析样本...\n\n")
  
  progressr::with_progress({
    
    p <- progressr::progressor(
      steps = length(sample_list),
      message = "分析细胞类型 Niche"
    )
    
    results <- future.apply::future_lapply(
      
      names(sample_list),
      
      function(sample_id) {
        
        process_single_sample(
          sample_id = sample_id,
          sample_list = sample_list,
          CONFIG = CONFIG,
          celltype_col = celltype_col,
          density_bins = density_bins,
          plot_overlay = plot_overlay,
          plot_composition = plot_composition,
          progressor = p
        )
      },
      
      future.seed = TRUE,
      future.chunk.size = 1
    )
  })
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 关闭并行
  future::plan(future::sequential)
  
  # ========================================
  # 5. 统计样本处理结果
  # ========================================
  
  print_sample_summary(results, sample_list, elapsed)
  
  # ========================================
  # 6. 生成综合分析
  # ========================================
  
  combined_data <- collect_combined_data(results)
  
  if (nrow(combined_data) > 0) {
    generate_combined_analysis(
      combined_data = combined_data,
      CONFIG = CONFIG,
      seurat_basename = seurat_basename,
      plot_heatmap = plot_heatmap,
      plot_combined = plot_combined
    )
  }
  
  # ========================================
  # 7. 最终总结
  # ========================================
  
  print_final_summary(results, sample_list, start_time, combined_data,
                     plot_overlay, plot_composition, plot_heatmap, plot_combined,
                     CONFIG)
  
  # ========================================
  # 8. 返回结果
  # ========================================
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  return(invisible(list(
    success = n_success,
    failed = n_failed,
    total = length(sample_list),
    elapsed_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    combined_data = combined_data,
    results = results
  )))
}


# ===================================================================
# 辅助函数
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

cat("✅ 08_plot_celltype.R 已加载\n")
```

---

### 08_plot_celltype_utils/00_operators.R

- **大小**: 552 B
- **修改时间**: 2025-11-06 22:24:04

```r
# ===================================================================
# 00_operators.R
# 基础操作符定义
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' %||% 操作符（空值默认值）
#'
#' @param a 主值
#' @param b 默认值（当a为NULL时使用）
#' @return 返回a（如果非NULL）或b
#'
#' @examples
#' NULL %||% "default"  # 返回 "default"
#' "value" %||% "default"  # 返回 "value"
#'
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
```

---

### 08_plot_celltype_utils/01_color_schemes.R

- **大小**: 2.56 KB
- **修改时间**: 2025-11-06 22:24:23

```r
# ===================================================================
# 01_color_schemes.R
# 统一的颜色方案管理
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 生成统一的zone颜色方案（支持任意数量的区域）
#'
#' @param n_zones zone数量，默认10
#' @return 命名的颜色向量
#'
#' @details
#' Zone_0 = 核心（深红色，高密度）
#' Zone_N = 外围（深蓝色，低密度）
#' 使用深红到深蓝的渐变色系
#'
#' @examples
#' zone_colors <- get_zone_colors(10)
#' zone_colors["Zone_0"]  # 深红色
#'
get_zone_colors <- function(n_zones = 10) {
  # 从深红到深蓝的渐变
  zone_colors <- colorRampPalette(c(
    "#67001f",  # 深红（Zone_0，核心，高密度）
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061"   # 深蓝（Zone_N-1，外围，低密度）
  ))(n_zones)
  
  # Zone_0 对应第一个颜色（深红）
  zone_names <- sprintf("Zone_%d", 0:(n_zones - 1))
  names(zone_colors) <- zone_names
  
  return(zone_colors)
}


#' 生成统一的细胞类型颜色
#'
#' @param celltypes 细胞类型向量
#' @return 命名的颜色向量
#'
#' @details
#' 根据细胞类型数量自动选择合适的调色板：
#' - ≤8种：使用 RColorBrewer Set2
#' - ≤12种：使用 RColorBrewer Set3
#' - >12种：组合 Set1 + Set2 + Set3
#'
#' @examples
#' celltype_colors <- get_celltype_colors(c("T cell", "B cell", "Macrophage"))
#'
get_celltype_colors <- function(celltypes) {
  n_celltypes <- length(celltypes)
  
  if (n_celltypes <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n_celltypes), "Set2")
  } else if (n_celltypes <= 12) {
    colors <- RColorBrewer::brewer.pal(n_celltypes, "Set3")
  } else {
    colors <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3")
    )[1:n_celltypes]
  }
  
  names(colors) <- celltypes
  return(colors)
}


#' 获取等高线颜色渐变
#'
#' @param n_breaks 等高线断点数量
#' @return 颜色向量（从深蓝到深红）
#'
#' @examples
#' contour_colors <- get_contour_colors(11)
#'
get_contour_colors <- function(n_breaks) {
  colorRampPalette(c(
    "#053061",  # 深蓝 (低密度)
    "#2166ac",
    "#4393c3",
    "#92c5de",
    "#d1e5f0",
    "#fddbc7",
    "#f4a582",
    "#d6604d",
    "#b2182b",
    "#67001f"   # 深红 (高密度)
  ))(n_breaks)
}
```

---

### 08_plot_celltype_utils/02_density_zones.R

- **大小**: 6.67 KB
- **修改时间**: 2025-11-07 10:54:32

```r
# ===================================================================
# 02_density_zones.R (完整修复版)
# 密度分区计算（真正支持边界扩展）
# Author: Assistant (Fixed Version)
# Date: 2025-11-07
# ===================================================================

#' 计算密度分区（基于等距切分，支持边界扩展）
#'
#' @param df 数据框，必须包含 col, row, ClockGene_High 列
#' @param density_bins 等高线分级数量，默认 10（对应0.1间隔）
#' @param expand_margin 边界扩展比例，默认 0.1 (10%)
#'
#' @return 包含以下元素的列表：
#'   - grid: 密度网格数据框（包含扩展区域）
#'   - spot_zones: 每个spot的zone分配
#'   - kde_result: KDE计算结果
#'   - equal_breaks: 等距断点
#'   - col_range: 原始切片列坐标范围
#'   - row_range: 原始切片行坐标范围
#'   - col_range_expanded: 扩展后的列坐标范围
#'   - row_range_expanded: 扩展后的行坐标范围
#'
calculate_density_zones <- function(df, density_bins = 10, expand_margin = 0.1) {
  
  require(dplyr)
  require(MASS)
  
  # ========================================
  # 1. 提取高表达点
  # ========================================
  high_points <- df %>% dplyr::filter(ClockGene_High)
  
  if (nrow(high_points) < 10) {
    warning("高表达点数量不足（< 10），无法计算密度")
    return(NULL)
  }
  
  # ========================================
  # 2. 计算原始范围和扩展范围
  # ========================================
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  
  # ✅ 关键修复：计算扩展后的范围（用于网格生成）
  col_margin <- diff(col_range) * expand_margin
  row_margin <- diff(row_range) * expand_margin
  
  col_range_expanded <- c(col_range[1] - col_margin, col_range[2] + col_margin)
  row_range_expanded <- c(row_range[1] - row_margin, row_range[2] + row_margin)
  
  cat(sprintf("   ✅ 密度计算范围:\n"))
  cat(sprintf("      原始: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range[1], col_range[2], row_range[1], row_range[2]))
  cat(sprintf("      扩展: col [%.1f, %.1f], row [%.1f, %.1f] (边距=%.0f%%)\n",
              col_range_expanded[1], col_range_expanded[2], 
              row_range_expanded[1], row_range_expanded[2],
              expand_margin * 100))
  
  # ========================================
  # 3. 使用扩展范围进行 KDE 计算
  # ========================================
  kde_result <- tryCatch({
    MASS::kde2d(
      x = high_points$col,
      y = high_points$row,
      n = 200,  # 网格分辨率
      lims = c(col_range_expanded, row_range_expanded)  # ✅ 使用扩展范围
    )
  }, error = function(e) {
    warning(sprintf("密度估计失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(kde_result)) return(NULL)
  
  # ========================================
  # 4. 转换为 data frame（保留扩展区域）
  # ========================================
  density_df <- expand.grid(
    col = kde_result$x,
    row = kde_result$y
  )
  density_df$density <- as.vector(kde_result$z)
  
  # 归一化密度到 [0, 1]
  max_density <- max(density_df$density, na.rm = TRUE)
  if (max_density > 0) {
    density_df$density_norm <- density_df$density / max_density
  } else {
    density_df$density_norm <- 0
  }
  
  # ✅ 关键修复：不要过滤掉扩展区域！注释掉原来的过滤代码
  # density_df <- density_df %>%
  #   dplyr::filter(
  #     col >= col_range[1] & col <= col_range[2],
  #     row >= row_range[1] & row <= row_range[2]
  #   )
  
  cat(sprintf("   ✅ 密度网格包含扩展区域: %d x %d = %d 个点\n",
              length(kde_result$x), length(kde_result$y), nrow(density_df)))
  
  # ========================================
  # 5. 等距切分密度区域
  # ========================================
  equal_breaks <- seq(0, 1, length.out = density_bins + 1)
  
  # 打印边界信息
  cat(sprintf("   ✅ Zone边界（等距切分，%d个区域）:\n", density_bins))
  for (i in 1:(length(equal_breaks) - 1)) {
    cat(sprintf("      Zone_%d: [%.2f, %.2f)\n", 
                density_bins - i, 
                equal_breaks[i], 
                equal_breaks[i + 1]))
  }
  
  # 分级（Zone_0 = 最高密度核心，Zone_9 = 最低密度外围）
  density_df$density_zone <- cut(
    density_df$density_norm,
    breaks = equal_breaks,
    labels = sprintf("Zone_%d", (density_bins - 1):0),
    include.lowest = TRUE,
    right = TRUE
  )
  
  # ========================================
  # 6. 为每个 spot 分配密度区域
  # ========================================
  spot_zones <- df %>%
    dplyr::select(col, row) %>%
    dplyr::mutate(
      col_idx = sapply(col, function(x) which.min(abs(kde_result$x - x))),
      row_idx = sapply(row, function(y) which.min(abs(kde_result$y - y)))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      grid_col = kde_result$x[col_idx],
      grid_row = kde_result$y[row_idx]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      density_df %>% dplyr::select(col, row, density_zone, density_norm),
      by = c("grid_col" = "col", "grid_row" = "row")
    ) %>%
    dplyr::select(col, row, density_zone, density_value = density_norm)
  
  # 处理 NA（使用最近邻填充）
  if (any(is.na(spot_zones$density_zone))) {
    na_count <- sum(is.na(spot_zones$density_zone))
    cat(sprintf("   ⚠️  %d 个 spots 需要最近邻填充\n", na_count))
    
    na_spots <- which(is.na(spot_zones$density_zone))
    
    for (idx in na_spots) {
      spot_col <- spot_zones$col[idx]
      spot_row <- spot_zones$row[idx]
      
      distances <- sqrt((density_df$col - spot_col)^2 + (density_df$row - spot_row)^2)
      valid_idx <- which(!is.na(density_df$density_zone))
      
      if (length(valid_idx) > 0) {
        nearest_valid <- valid_idx[which.min(distances[valid_idx])]
        spot_zones$density_zone[idx] <- density_df$density_zone[nearest_valid]
        spot_zones$density_value[idx] <- density_df$density_norm[nearest_valid]
      }
    }
  }
  
  # ========================================
  # 7. 返回结果
  # ========================================
  return(list(
    grid = density_df,                        # 包含扩展区域的完整网格
    spot_zones = spot_zones,
    kde_result = kde_result,
    equal_breaks = equal_breaks,
    col_range = col_range,                    # 原始范围
    row_range = row_range,
    col_range_expanded = col_range_expanded,  # 扩展范围
    row_range_expanded = row_range_expanded
  ))
}
```

---

### 08_plot_celltype_utils/03_plot_overlay.R

- **大小**: 10.01 KB
- **修改时间**: 2025-11-07 11:24:05

```r
# ===================================================================
# 03_plot_overlay.R (修复版)
# 细胞类型+密度叠加图（使用raster，无网格线）
# Author: Assistant (Fixed Version)
# Date: 2025-11-07
# ===================================================================

#' 绘制细胞类型和密度区域叠加图
#'
#' @param df 数据框，包含细胞类型和坐标信息
#' @param density_data 密度计算结果（来自 calculate_density_zones）
#' @param sample_id 样本ID
#' @param CONFIG 配置列表
#'
#' @return ggplot对象
#'


plot_celltype_density_overlay <- function(df, density_data, sample_id, CONFIG) {
  
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  require(RANN)
  
  # ========================================
  # 1. 准备数据
  # ========================================
  
  n_zones <- length(unique(density_data$grid$density_zone))
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  zone_colors <- CONFIG$colors$zone_colors %||% get_zone_colors(n_zones)
  
  # 清理 celltype
  df$celltype_clean <- as.character(df$celltype_clean)
  df$celltype_clean[is.na(df$celltype_clean) | df$celltype_clean == ""] <- "Unknown"
  all_celltypes <- sort(unique(df$celltype_clean))
  
  # 获取配置的颜色
  celltype_colors <- CONFIG$colors$celltype_colors
  
  # 确保所有类型都有颜色
  missing_types <- setdiff(all_celltypes, names(celltype_colors))
  if (length(missing_types) > 0) {
    extra_colors <- rainbow(length(missing_types))
    names(extra_colors) <- missing_types
    celltype_colors <- c(celltype_colors, extra_colors)
  }
  
  # 只保留实际存在的类型，并确保顺序一致
  celltype_colors <- celltype_colors[all_celltypes]
  
  cat("   📊 Celltype 颜色映射:\n")
  for (ct in all_celltypes) {
    cat(sprintf("      %s → %s\n", ct, celltype_colors[ct]))
  }
  
  # ========================================
  # 2. 坐标范围
  # ========================================
  
  col_range_raw <- density_data$col_range
  row_range_raw <- density_data$row_range
  
  if (!is.null(density_data$col_range_expanded)) {
    col_limits <- density_data$col_range_expanded
    row_limits <- density_data$row_range_expanded
  } else {
    expand_margin <- 0.1
    col_margin <- diff(col_range_raw) * expand_margin
    row_margin <- diff(row_range_raw) * expand_margin
    col_limits <- c(col_range_raw[1] - col_margin, col_range_raw[2] + col_margin)
    row_limits <- c(row_range_raw[1] - row_margin, row_range_raw[2] + row_margin)
  }
  
  cat(sprintf("   ✅ 绘图范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_limits[1], col_limits[2], row_limits[1], row_limits[2]))
  
  # ========================================
  # 3. 准备等高线数据
  # ========================================
  
  zone_density_ranges <- density_data$grid %>%
    group_by(density_zone) %>%
    summarise(
      min_density = min(density_norm, na.rm = TRUE),
      max_density = max(density_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(density_zone = factor(density_zone, levels = zone_levels)) %>%
    arrange(density_zone)
  
  zone_labels <- zone_density_ranges %>%
    mutate(zone_label = sprintf("%s (%.2f-%.2f)", density_zone, min_density, max_density)) %>%
    pull(zone_label)
  names(zone_labels) <- as.character(zone_density_ranges$density_zone)
  
  contour_breaks <- density_data$equal_breaks
  contour_data <- density_data$grid %>%
    mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # ========================================
  # 4. 准备细胞数据
  # ========================================
  
  df_filtered <- df %>% 
    filter(!is.na(density_zone)) %>%
    mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      celltype_clean = factor(celltype_clean, levels = all_celltypes)
    )
  
  # 计算细胞大小
  if (nrow(df_filtered) > 10000) {
    sample_idx <- sample(nrow(df_filtered), 10000)
    coords_sample <- df_filtered[sample_idx, c("col", "row")]
  } else {
    coords_sample <- df_filtered[, c("col", "row")]
  }
  
  nn_dist <- RANN::nn2(coords_sample, k = 2)$nn.dists[, 2]
  median_dist <- median(nn_dist, na.rm = TRUE)
  square_size <- median_dist * 1.0
  
  cat(sprintf("   📏 细胞正方形大小: %.3f\n", square_size))
  
  # 等高线颜色
  contour_colors <- get_contour_colors(length(contour_breaks))
  
  # 统一的图例样式参数
  legend_key_width <- 1.0
  legend_key_height <- 0.7
  legend_text_size <- 10
  legend_title_size <- 11
  
  # ========================================
  # 5. 绘图
  # ========================================
  
  p <- ggplot() +
    # ========================================
    # Layer 1: 细胞类型（底层）- 使用 FILL
    # ========================================
    geom_tile(
      data = df_filtered,
      aes(x = col, y = row, fill = celltype_clean),  # ✅ 使用 fill
      width = square_size,
      height = square_size,
      color = NA,  # ✅ 不要边框
      alpha = 1
    ) +
    scale_fill_manual(
      values = celltype_colors,  # ✅ 必须是命名向量
      name = "Cell Type",
      breaks = all_celltypes,
      drop = TRUE,
      na.value = "gray50",
      guide = guide_legend(
        order = 2,
        override.aes = list(
          alpha = 1,
          color = NA  # ✅ 图例中也不要边框
        ),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keywidth = unit(legend_key_width, "cm"),
        keyheight = unit(legend_key_height, "cm")
      )
    ) +
    
    # ========================================
    # 新的scale用于density zones
    # ========================================
    ggnewscale::new_scale_fill() +
    
    # Layer 2: Zone填充（半透明覆盖层）
    geom_raster(
      data = contour_data,
      aes(x = col, y = row, fill = density_zone),
      alpha = 0.3,
      interpolate = TRUE
    ) +
    scale_fill_manual(
      values = zone_colors,
      labels = zone_labels,
      name = "Density Zones\n(Zone_0 = Core Red → Zone_9 = Outer Blue)",
      breaks = zone_levels,
      na.value = "transparent",
      drop = FALSE,
      guide = guide_legend(
        order = 1,
        override.aes = list(
          alpha = 0.8,
          color = "gray40",
          linewidth = 0.3
        ),
        title.position = "top",
        title.hjust = 0,
        ncol = 1,
        byrow = TRUE,
        keywidth = unit(legend_key_width, "cm"),
        keyheight = unit(legend_key_height, "cm")
      )
    ) +
    
    # ========================================
    # 为等高线准备新的color scale
    # ========================================
    ggnewscale::new_scale_color() +
    
    # Layer 3: 等高线边界
    geom_contour(
      data = contour_data,
      aes(x = col, y = row, z = density_norm, color = after_stat(level)),
      breaks = contour_breaks,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    scale_color_gradientn(
      colors = contour_colors,
      limits = c(min(contour_breaks), max(contour_breaks)),
      guide = "none"
    ) +
    
    # ========================================
    # 坐标系统
    # ========================================
    scale_x_continuous(limits = col_limits, expand = c(0, 0)) +
    scale_y_reverse(limits = rev(row_limits), expand = c(0, 0)) +
    coord_fixed(ratio = 1, xlim = col_limits, ylim = rev(row_limits), clip = "off") +
    
    # ========================================
    # 标题和主题
    # ========================================
    labs(
      title = sprintf("Cell Type Distribution in Density Zones - %s", sample_id),
      subtitle = sprintf("Bottom = Cell types | Middle = Density zones (α=0.3) | Top = %d contour lines", 
                        length(contour_breaks)),
      x = NULL, y = NULL
    ) +
    
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30", margin = margin(b = 10)),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = unit(0.6, "cm"),
      legend.title = element_text(size = legend_title_size, face = "bold", hjust = 0, margin = margin(b = 8)),
      legend.text = element_text(size = legend_text_size, lineheight = 1.2, margin = margin(l = 3, r = 5, t = 2, b = 2)),
      legend.key = element_rect(color = "gray60", fill = NA, linewidth = 0.3),
      legend.key.spacing.y = unit(0.2, "cm"),
      legend.background = element_rect(fill = "white", color = "gray50", linewidth = 0.5),
      legend.margin = margin(10, 10, 10, 10),
      plot.margin = margin(15, 20, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# ========================================
# 辅助函数
# ========================================

#' 生成 zone 颜色（红到蓝渐变）
get_zone_colors <- function(n_zones) {
  colorRampPalette(c(
    "#B2182B",  # 深红（核心高密度）
    "#EF8A62",  # 浅红
    "#FDDBC7",  # 粉色
    "#F7F7F7",  # 白色（中间）
    "#D1E5F0",  # 浅蓝
    "#67A9CF",  # 蓝色
    "#2166AC"   # 深蓝（外围低密度）
  ))(n_zones)
}

#' 生成等高线颜色（紫色渐变）
get_contour_colors <- function(n_contours) {
  colorRampPalette(c(
    "#542788",  # 深紫
    "#8073AC",  # 中紫
    "#B2ABD2",  # 浅紫
    "#D8DAEB"   # 淡紫
  ))(n_contours)
}

#' 为细胞类型生成颜色
get_celltype_colors <- function(celltypes) {
  require(RColorBrewer)
  n <- length(celltypes)
  
  if (n <= 3) {
    colors <- brewer.pal(3, "Set2")[1:n]
  } else if (n <= 12) {
    colors <- brewer.pal(n, "Set3")
  } else {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n)
  }
  
  names(colors) <- celltypes
  return(colors)
}

```

---

### 08_plot_celltype_utils/04_plot_composition.R

- **大小**: 3.01 KB
- **修改时间**: 2025-11-06 22:25:19

```r
# ===================================================================
# 04_plot_composition.R
# 区域组成柱状图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制zone组成柱状图
#'
#' @param zone_composition zone组成数据框
#' @param sample_id 样本ID
#' @param CONFIG 配置列表
#'
#' @return patchwork组合图（细胞类型组成 + spot数量）
#'
#' @examples
#' p <- plot_zone_composition(zone_comp, "Sample_01", CONFIG)
#'
plot_zone_composition <- function(zone_composition, sample_id, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # 使用统一的颜色方案
  n_zones <- length(unique(zone_composition$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  celltype_colors <- get_celltype_colors(unique(zone_composition$celltype_clean))
  
  # 确保zone按 Zone_0, Zone_1, ... 排序
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  zone_composition <- zone_composition %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 图1：细胞类型组成堆叠柱状图
  p1 <- ggplot(zone_composition, aes(x = density_zone, y = percentage, fill = celltype_clean)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = sprintf("Cell Type Composition by Density Zone - %s", sample_id),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  # 图2：Zone的spot数量
  zone_totals <- zone_composition %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(total = sum(count), .groups = "drop")
  
  p2 <- ggplot(zone_totals, aes(x = density_zone, y = total, fill = density_zone)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.5) +
    geom_text(aes(label = total), vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Total Spots per Density Zone",
      x = "Density Zone (Zone_0=Core → Higher=Outer)",
      y = "Count"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # 合并
  p_combined <- p1 / p2 + plot_layout(heights = c(2, 1))
  
  return(p_combined)
}
```

---

### 08_plot_celltype_utils/05_plot_heatmap.R

- **大小**: 3.32 KB
- **修改时间**: 2025-11-06 22:25:31

```r
# ===================================================================
# 05_plot_heatmap.R
# 合并热图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制合并热图（所有样本）
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表
#'
#' @return patchwork组合图（zone颜色条 + 热图）
#'
#' @examples
#' p <- plot_combined_heatmap(combined_data, CONFIG)
#'
plot_combined_heatmap <- function(combined_data, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # 计算平均百分比
  heatmap_data <- combined_data %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(
      mean_pct = mean(percentage),
      sd_pct = sd(percentage),
      n_samples = n(),
      .groups = "drop"
    )
  
  # 确保zone按顺序排列
  n_zones <- length(unique(heatmap_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  heatmap_data <- heatmap_data %>%
    dplyr::mutate(density_zone = factor(density_zone, levels = zone_levels))
  
  # 热图主体
  p <- ggplot(heatmap_data, aes(x = density_zone, y = celltype_clean, fill = mean_pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_pct)), size = 3.5, color = "black", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white", "#fee090", "#fc8d59", "#d73027"),
      name = "Mean %",
      limits = c(0, NA),
      guide = guide_colorbar(
        barwidth = 1.5,
        barheight = 15,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = "Cell Type Composition Across Density Zones (All Samples)",
      subtitle = sprintf("Averaged across %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30", margin = margin(b = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 1)
    )
  
  # 添加zone颜色参考条（顶部）
  zone_bar_data <- data.frame(
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 1
  )
  
  p_zone_bar <- ggplot(zone_bar_data, aes(x = density_zone, y = y_position, fill = density_zone)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # 合并图形
  p_final <- p_zone_bar / p + plot_layout(heights = c(0.05, 1))
  
  return(p_final)
}
```

---

### 08_plot_celltype_utils/06_plot_combined.R

- **大小**: 4.39 KB
- **修改时间**: 2025-11-06 22:25:45

```r
# ===================================================================
# 06_plot_combined.R
# 综合分析图绘制
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 绘制综合分析图（箱线图 + 趋势图）
#'
#' @param combined_data 合并的zone组成数据
#' @param CONFIG 配置列表
#'
#' @return patchwork组合图
#'
#' @examples
#' p <- plot_combined_analysis(combined_data, CONFIG)
#'
plot_combined_analysis <- function(combined_data, CONFIG) {
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  # 获取统一的颜色方案
  n_zones <- length(unique(combined_data$density_zone))
  zone_colors <- get_zone_colors(n_zones)
  zone_levels <- sprintf("Zone_%d", 0:(n_zones - 1))
  celltype_colors <- get_celltype_colors(unique(combined_data$celltype_clean))
  
  # 确保zone按顺序排列
  combined_data <- combined_data %>%
    dplyr::mutate(
      density_zone = factor(density_zone, levels = zone_levels),
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    )
  
  # 1. 箱线图
  p1 <- ggplot(combined_data, aes(x = density_zone, y = percentage, fill = density_zone)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    facet_wrap(~celltype_clean, scales = "free_y", ncol = 4) +
    labs(
      title = "Cell Type Percentage Distribution by Density Zone",
      subtitle = sprintf("Data from %d samples", length(unique(combined_data$sample))),
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
  
  # 2. 趋势图
  trend_data <- combined_data %>%
    dplyr::group_by(celltype_clean, zone_numeric, density_zone) %>%
    dplyr::summarise(
      mean_pct = mean(percentage),
      se_pct = sd(percentage) / sqrt(n()),
      .groups = "drop"
    )
  
  p2 <- ggplot(trend_data, aes(x = zone_numeric, y = mean_pct, color = celltype_clean, group = celltype_clean)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct), 
      width = 0.2, 
      linewidth = 0.8,
      alpha = 0.7
    ) +
    scale_color_manual(values = celltype_colors, name = "Cell Type") +
    scale_x_continuous(
      breaks = 0:(n_zones - 1),
      labels = zone_levels
    ) +
    labs(
      title = "Cell Type Enrichment Trend Across Density Zones",
      subtitle = "Mean ± SE across all samples",
      x = "Density Zone (Zone_0=Core/High → Higher=Outer/Low)",
      y = "Mean Percentage (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  # 3. 添加zone颜色参考条
  zone_ref_data <- data.frame(
    zone_numeric = 0:(n_zones - 1),
    density_zone = factor(zone_levels, levels = zone_levels),
    y_position = 0
  )
  
  p2 <- p2 +
    geom_tile(
      data = zone_ref_data,
      aes(x = zone_numeric, y = y_position, fill = density_zone),
      height = max(trend_data$mean_pct) * 0.05,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = zone_colors, guide = "none")
  
  # 合并
  p_combined <- p1 / p2 + plot_layout(heights = c(2, 1.2))
  
  return(p_combined)
}
```

---

### 08_plot_celltype_utils/07_statistics.R

- **大小**: 1.47 KB
- **修改时间**: 2025-11-06 22:25:57

```r
# ===================================================================
# 07_statistics.R
# 统计摘要生成
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 生成统计摘要
#'
#' @param combined_data 合并的zone组成数据
#'
#' @return 统计摘要数据框
#'
#' @details
#' 计算内容：
#' - 每种细胞类型的平均百分比和标准差
#' - 富集最多/最少的zone
#' - 核心区（Zone_0和Zone_1）vs外围区的富集差异
#' - 样本数量
#'
#' @examples
#' summary <- generate_summary_statistics(combined_data)
#'
generate_summary_statistics <- function(combined_data) {
  
  require(dplyr)
  
  # 计算每种细胞类型在不同区域的富集情况
  summary <- combined_data %>%
    dplyr::mutate(zone_numeric = as.numeric(gsub("Zone_", "", density_zone))) %>%
    dplyr::group_by(celltype_clean) %>%
    dplyr::summarise(
      mean_pct_all = mean(percentage),
      sd_pct_all = sd(percentage),
      max_zone = density_zone[which.max(percentage)],
      max_pct = max(percentage),
      min_zone = density_zone[which.min(percentage)],
      min_pct = min(percentage),
      # Zone_0和Zone_1是核心区，其他是外围
      core_enrichment = mean(percentage[zone_numeric <= 1]) - mean(percentage[zone_numeric > 1]),
      n_samples = length(unique(sample)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(core_enrichment))
  
  return(summary)
}
```

---

### 08_plot_celltype_utils/08_validation.R

- **大小**: 3.71 KB
- **修改时间**: 2025-11-07 14:03:40

```r
#!/usr/bin/env Rscript
# ===================================================================
# 验证模块
# ===================================================================

#' 验证输入参数
#' 
#' @param sample_list 样本列表
#' @param CONFIG 配置对象
validate_inputs <- function(sample_list, CONFIG) {
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  # 验证必需目录
  required_dirs <- c("overlay", "celltype", "composition", "heatmaps", "combined")
  
  for (dir_name in required_dirs) {
    if (is.null(CONFIG$dirs[[dir_name]])) {
      stop(sprintf("❌ CONFIG$dirs$%s 未定义", dir_name))
    }
    if (!dir.exists(CONFIG$dirs[[dir_name]])) {
      dir.create(CONFIG$dirs[[dir_name]], recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  invisible(TRUE)
}


#' 验证必需函数
validate_required_functions <- function() {
  
  required_functions <- c(
    "calculate_density_zones",
    "plot_celltype_density_overlay",
    "plot_zone_composition",
    "plot_combined_heatmap",
    "plot_combined_analysis",
    "generate_summary_statistics"
  )
  
  missing_funcs <- required_functions[!sapply(required_functions, exists)]
  
  if (length(missing_funcs) > 0) {
    stop(sprintf("❌ 缺少必需函数: %s", paste(missing_funcs, collapse = ", ")))
  }
  
  invisible(TRUE)
}


#' 设置颜色方案
#' 
#' @param first_sample 第一个样本（用于提取细胞类型）
#' @param CONFIG 配置对象
#' @param celltype_col 细胞类型列名
#' @param density_bins 密度分区数量
setup_colors <- function(first_sample, CONFIG, celltype_col, density_bins) {
  
  # 从第一个样本获取所有细胞类型
  all_celltypes <- sort(unique(as.character(first_sample[[celltype_col]][,1])))
  
  if (is.null(CONFIG$colors$celltype_colors)) {
    CONFIG$colors$celltype_colors <- get_celltype_colors(all_celltypes)
    cat(sprintf("🎨 已生成 %d 种细胞类型颜色方案\n", length(CONFIG$colors$celltype_colors)))
  }
  
  if (is.null(CONFIG$colors$zone_colors)) {
    CONFIG$colors$zone_colors <- get_zone_colors(density_bins)
  }
  
  invisible(CONFIG)
}


#' 验证样本数据
#' 
#' @param seurat_subset Seurat 对象
#' @param sample_id 样本 ID
#' @param celltype_col 细胞类型列名
#' 
#' @return 验证结果列表
validate_sample_data <- function(seurat_subset, sample_id, celltype_col) {
  
  # 检查数据量
  if (ncol(seurat_subset) == 0) {
    return(list(valid = FALSE, error = "No data"))
  }
  
  # 获取坐标
  coords <- tryCatch({
    Seurat::GetTissueCoordinates(
      seurat_subset,
      cols = c("row", "col"),
      scale = NULL
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(coords)) {
    return(list(valid = FALSE, error = "Cannot get coordinates"))
  }
  
  # 合并元数据
  df <- seurat_subset@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(
      coords %>% tibble::rownames_to_column("barcode"), 
      by = "barcode"
    ) %>%
    dplyr::filter(!is.na(col), !is.na(row))
  
  if (nrow(df) == 0) {
    return(list(valid = FALSE, error = "No valid coordinates"))
  }
  
  # 检查必需列
  if (!celltype_col %in% colnames(df)) {
    return(list(valid = FALSE, error = "Missing celltype column"))
  }
  
  if (!"ClockGene_High" %in% colnames(df)) {
    return(list(valid = FALSE, error = "Missing ClockGene_High column"))
  }
  
  # 清理细胞类型
  df$celltype_clean <- as.character(df[[celltype_col]])
  df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
  
  return(list(valid = TRUE, df = df))
}

cat("✅ 08_validation.R 已加载\n")
```

---

### 08_plot_celltype_utils/09_save_plots.R

- **大小**: 5.5 KB
- **修改时间**: 2025-11-07 14:04:12

```r
#!/usr/bin/env Rscript
# ===================================================================
# 图形保存模块
# ===================================================================

#' 处理单个样本
#' 
#' @param sample_id 样本 ID
#' @param sample_list 样本列表
#' @param CONFIG 配置对象
#' @param celltype_col 细胞类型列名
#' @param density_bins 密度分区数量
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param progressor 进度条对象
#' 
#' @return 处理结果
process_single_sample <- function(sample_id, sample_list, CONFIG, 
                                  celltype_col, density_bins,
                                  plot_overlay, plot_composition,
                                  progressor) {
  
  tryCatch({
    
    # 1. 获取并验证数据
    seurat_subset <- sample_list[[sample_id]]
    validation <- validate_sample_data(seurat_subset, sample_id, celltype_col)
    
    if (!validation$valid) {
      progressor(message = sprintf("⚠️  %s - %s", sample_id, validation$error))
      return(list(sample = sample_id, success = FALSE, error = validation$error))
    }
    
    df <- validation$df
    
    # 统计基本信息
    n_spots <- nrow(df)
    n_high <- sum(df$ClockGene_High, na.rm = TRUE)
    high_pct <- 100 * mean(df$ClockGene_High, na.rm = TRUE)
    
    # 2. 计算密度区域
    density_data <- calculate_density_zones(
      df = df,
      density_bins = density_bins,
      expand_margin = CONFIG$plot$expand_margin %||% 0.1
    )
    
    if (is.null(density_data)) {
      progressor(message = sprintf("⚠️  %s - 密度计算失败", sample_id))
      return(list(sample = sample_id, success = FALSE, error = "Density calculation failed"))
    }
    
    # 合并密度信息
    df <- df %>%
      dplyr::left_join(
        density_data$spot_zones %>% 
          dplyr::select(col, row, density_zone, density_value),
        by = c("col", "row")
      )
    
    # 3. 计算区域组成
    zone_composition <- df %>%
      dplyr::filter(!is.na(density_zone)) %>%
      dplyr::group_by(density_zone, celltype_clean) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(density_zone) %>%
      dplyr::mutate(
        total = sum(count),
        percentage = 100 * count / total
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sample = sample_id)
    
    # 4. 绘制并保存图形
    plot_result <- save_sample_plots(
      df = df,
      density_data = density_data,
      zone_composition = zone_composition,
      sample_id = sample_id,
      CONFIG = CONFIG,
      plot_overlay = plot_overlay,
      plot_composition = plot_composition
    )
    
    # 更新进度
    progressor(message = sprintf("✅ %s (%.2f MB)", sample_id, plot_result$total_size_mb))
    
    # 5. 返回结果
    return(list(
      sample = sample_id,
      success = TRUE,
      zone_composition = zone_composition,
      n_spots = n_spots,
      n_high = n_high,
      high_pct = high_pct,
      n_zones = length(unique(zone_composition$density_zone)),
      n_celltypes = length(unique(zone_composition$celltype_clean)),
      n_na_zones = sum(is.na(df$density_zone)),
      output_files = plot_result$output_files,
      total_size_mb = plot_result$total_size_mb
    ))
    
  }, error = function(e) {
    progressor(message = sprintf("❌ %s - %s", sample_id, e$message))
    return(list(
      sample = sample_id,
      success = FALSE,
      error = as.character(e$message)
    ))
  })
}


#' 保存样本图形
#' 
#' @param df 数据框
#' @param density_data 密度数据
#' @param zone_composition 区域组成
#' @param sample_id 样本 ID
#' @param CONFIG 配置对象
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' 
#' @return 保存结果
save_sample_plots <- function(df, density_data, zone_composition, sample_id, CONFIG,
                              plot_overlay, plot_composition) {
  
  output_files <- list()
  total_size <- 0
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # 叠加图
  if (plot_overlay) {
    p_overlay <- plot_celltype_density_overlay(
      df = df,
      density_data = density_data,
      sample_id = sample_id,
      CONFIG = CONFIG
    )
    
    overlay_file <- file.path(
      CONFIG$dirs$overlay, 
      sprintf("celltype_overlay_%s.pdf", safe_name)
    )
    
    ggplot2::ggsave(
      overlay_file,
      plot = p_overlay,
      width = 12, 
      height = 10,
      dpi = CONFIG$plot$dpi %||% 300,
      bg = "white"
    )
    
    output_files$overlay <- overlay_file
    total_size <- total_size + file.size(overlay_file)
  }
  
  # 组成图
  if (plot_composition) {
    p_comp <- plot_zone_composition(
      zone_composition = zone_composition,
      sample_id = sample_id,
      CONFIG = CONFIG
    )
    
    composition_file <- file.path(
      CONFIG$dirs$composition, 
      sprintf("composition_%s.pdf", safe_name)
    )
    
    ggplot2::ggsave(
      composition_file,
      plot = p_comp,
      width = 12, 
      height = 6,
      dpi = CONFIG$plot$dpi %||% 300,
      bg = "white"
    )
    
    output_files$composition <- composition_file
    total_size <- total_size + file.size(composition_file)
  }
  
  return(list(
    output_files = output_files,
    total_size_mb = total_size / 1024^2
  ))
}

cat("✅ 09_save_plots.R 已加载\n")
```

---

### 08_plot_celltype_utils/10_summary.R

- **大小**: 8.85 KB
- **修改时间**: 2025-11-07 14:04:35

```r
#!/usr/bin/env Rscript
# ===================================================================
# 汇总统计模块
# ===================================================================

#' 打印样本汇总
#' 
#' @param results 结果列表
#' @param sample_list 样本列表
#' @param elapsed 耗时
print_sample_summary <- function(results, sample_list, elapsed) {
  
  n_success <- sum(sapply(results, function(x) x$success))
  n_failed <- length(results) - n_success
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本处理完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d (%.1f%%)\n", 
              n_success, 
              length(sample_list),
              100 * n_success / length(sample_list)))
  
  if (n_failed > 0) {
    cat(sprintf("❌ 失败: %d/%d\n\n", n_failed, length(sample_list)))
    cat("失败样本:\n")
    for (res in results) {
      if (!res$success) {
        cat(sprintf("  • %s: %s\n", res$sample, res$error))
      }
    }
    cat("\n")
  }
  
  if (n_success > 0) {
    cat("成功样本:\n")
    cat(sprintf("%-30s %8s %7s %8s %7s %8s %7s %10s\n",
                "样本", "Spots", "High", "High%", "Zones", "Types", "NA", "大小(MB)"))
    cat(paste(rep("-", 100), collapse = ""), "\n")
    
    total_size <- 0
    total_spots <- 0
    
    for (res in results) {
      if (res$success) {
        cat(sprintf("%-30s %8d %7d %7.2f%% %7d %8d %7d %10.2f\n",
                    res$sample,
                    res$n_spots,
                    res$n_high,
                    res$high_pct,
                    res$n_zones,
                    res$n_celltypes,
                    res$n_na_zones,
                    res$total_size_mb))
        
        total_size <- total_size + res$total_size_mb
        total_spots <- total_spots + res$n_spots
      }
    }
    
    if (n_success > 1) {
      cat(paste(rep("-", 100), collapse = ""), "\n")
      cat(sprintf("%-30s %8d %7s %8s %7s %8s %7s %10.2f\n",
                  "总计",
                  total_spots,
                  "-", "-", "-", "-", "-",
                  total_size))
    }
    
    cat("\n")
  }
  
  cat(sprintf("⏱️  样本处理耗时: %.2f 秒 (平均 %.2f 秒/样本)\n\n", 
              as.numeric(elapsed),
              as.numeric(elapsed) / length(sample_list)))
  
  invisible(NULL)
}


#' 收集合并数据
#' 
#' @param results 结果列表
#' @return 合并的数据框
collect_combined_data <- function(results) {
  
  combined_data <- data.frame()
  
  for (res in results) {
    if (res$success) {
      combined_data <- dplyr::bind_rows(combined_data, res$zone_composition)
    }
  }
  
  return(combined_data)
}


#' 生成综合分析
#' 
#' @param combined_data 合并数据
#' @param CONFIG 配置对象
#' @param seurat_basename 基础名
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
generate_combined_analysis <- function(combined_data, CONFIG, seurat_basename,
                                       plot_heatmap, plot_combined) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   生成综合统计图\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  combined_start <- Sys.time()
  
  main_title <- seurat_basename %||% "Seurat Object"
  
  # 热图
  if (plot_heatmap) {
    cat("📊 生成细胞类型热图...\n")
    
    tryCatch({
      p_heatmap <- plot_combined_heatmap(
        combined_data = combined_data, 
        CONFIG = CONFIG
      ) + ggplot2::ggtitle(main_title)
      
      heatmap_file <- file.path(
        CONFIG$dirs$heatmaps, 
        "celltype_heatmap_all_samples.pdf"
      )
      
      ggplot2::ggsave(
        heatmap_file,
        plot = p_heatmap, 
        width = 14, 
        height = 10, 
        dpi = CONFIG$plot$dpi %||% 300, 
        bg = "white"
      )
      
      cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                  basename(heatmap_file),
                  file.size(heatmap_file) / 1024^2))
    }, error = function(e) {
      cat(sprintf("   ⚠️  热图生成失败: %s\n", e$message))
    })
  }
  
  # 综合分析图
  if (plot_combined) {
    cat("📊 生成综合分析图...\n")
    
    tryCatch({
      p_combined <- plot_combined_analysis(
        combined_data = combined_data, 
        CONFIG = CONFIG
      ) + ggplot2::ggtitle(main_title)
      
      combined_file <- file.path(
        CONFIG$dirs$combined, 
        "combined_analysis.pdf"
      )
      
      ggplot2::ggsave(
        combined_file,
        plot = p_combined, 
        width = 16, 
        height = 12, 
        dpi = CONFIG$plot$dpi %||% 300, 
        bg = "white"
      )
      
      cat(sprintf("   ✅ 保存: %s (%.2f MB)\n", 
                  basename(combined_file),
                  file.size(combined_file) / 1024^2))
    }, error = function(e) {
      cat(sprintf("   ⚠️  综合图生成失败: %s\n", e$message))
    })
  }
  
  # 保存数据
  cat("💾 保存统计数据...\n")
  
  composition_csv <- file.path(
    CONFIG$dirs$composition, 
    "celltype_composition_all_samples.csv"
  )
  write.csv(combined_data, composition_csv, row.names = FALSE)
  cat(sprintf("   ✅ 组成数据: %s\n", basename(composition_csv)))
  
  tryCatch({
    summary_stats <- generate_summary_statistics(combined_data)
    summary_csv <- file.path(
      CONFIG$dirs$composition, 
      "summary_statistics.csv"
    )
    write.csv(summary_stats, summary_csv, row.names = FALSE)
    cat(sprintf("   ✅ 汇总统计: %s\n", basename(summary_csv)))
  }, error = function(e) {
    cat(sprintf("   ⚠️  统计计算失败: %s\n", e$message))
  })
  
  combined_end <- Sys.time()
  combined_elapsed <- difftime(combined_end, combined_start, units = "secs")
  
  cat(sprintf("\n⏱️  综合图生成耗时: %.2f 秒\n", as.numeric(combined_elapsed)))
  
  invisible(NULL)
}


#' 打印最终汇总
#' 
#' @param results 结果列表
#' @param sample_list 样本列表
#' @param start_time 开始时间
#' @param combined_data 合并数据
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param plot_heatmap 是否绘制热图
#' @param plot_combined 是否绘制综合图
#' @param CONFIG 配置对象
print_final_summary <- function(results, sample_list, start_time, combined_data,
                               plot_overlay, plot_composition, plot_heatmap, 
                               plot_combined, CONFIG) {
  
  total_elapsed <- difftime(Sys.time(), start_time, units = "secs")
  n_success <- sum(sapply(results, function(x) x$success))
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("✅ 成功: %d/%d\n", n_success, length(sample_list)))
  cat(sprintf("⏱️  总耗时: %.2f 秒 (%.2f 分钟)\n", 
              as.numeric(total_elapsed),
              as.numeric(total_elapsed) / 60))
  
  if (n_success > 0) {
    cat("\n📊 生成内容:\n")
    if (plot_overlay) 
      cat(sprintf("   • 叠加图: %d 个\n", n_success))
    if (plot_composition) 
      cat(sprintf("   • 组成图: %d 个\n", n_success))
    if (plot_heatmap && nrow(combined_data) > 0) 
      cat("   • 热图: 1 个\n")
    if (plot_combined && nrow(combined_data) > 0) 
      cat("   • 综合图: 1 个\n")
  }
  
  cat("\n📁 输出目录:\n")
  cat(sprintf("   • Overlay:     %s\n", CONFIG$dirs$overlay))
  cat(sprintf("   • Composition: %s\n", CONFIG$dirs$composition))
  cat(sprintf("   • Heatmaps:    %s\n", CONFIG$dirs$heatmaps))
  cat(sprintf("   • Combined:    %s\n", CONFIG$dirs$combined))
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(NULL)
}

cat("✅ 10_summary.R 已加载\n")
```

---

### 09_save_results.R

- **大小**: 1022 B
- **修改时间**: 2025-11-05 17:10:12

```r
#!/usr/bin/env Rscript
# ===================================================================
# 保存结果
# ===================================================================

save_results <- function(seurat_obj, config) {
  cat("💾 保存结果...\n")
  
  # 保存metadata
  write.csv(
    seurat_obj@meta.data, 
    file.path(config$metadata_dir, "Lymph_2-25M_clockgene_metadata.csv"),
    row.names = TRUE
  )
  
  # 可选：保存完整对象
  if (config$save_full_object) {
    saveRDS(
      seurat_obj, 
      file.path(config$metadata_dir, "Lymph_2-25M_with_clockgene_niche.rds")
    )
  }
  
  cat("✅ 结果保存完成\n\n")
}

print_summary <- function(config) {
  cat("📊 文件统计:\n")
  cat(sprintf("   图形文件夹: %s\n", config$figure_dir))
  cat(sprintf("   - Isoheight: %d 个文件\n", length(list.files(config$dirs$isoheight))))
  cat(sprintf("   - Spatial: %d 个文件\n", length(list.files(config$dirs$spatial))))
  cat("\n✅ 全部完成！\n")
}
```

---

### 10_batch_processing.R

- **大小**: 3.66 KB
- **修改时间**: 2025-11-07 14:17:44

```r
#!/usr/bin/env Rscript
# ===================================================================
# 批量处理核心模块
# ===================================================================

#' 处理所有文件
#'
#' @param seurat_files 文件列表
#' @param gene_list 基因列表
#' @param CONFIG 配置对象
#' 
#' @return 处理结果列表
#'
process_all_files <- function(seurat_files, gene_list, CONFIG) {
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   开始批量处理\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  results <- list()
  
  for (i in seq_along(seurat_files)) {
    
    # 打印进度
    print_progress_header(i, length(seurat_files))
    
    # 处理单个文件
    result <- process_seurat_file(
      seurat_path = seurat_files[i],
      gene_list = gene_list,
      base_config = CONFIG
    )
    
    results[[i]] <- result
    
    # 估计剩余时间
    if (i < length(seurat_files)) {
      estimate_remaining_time(results, seurat_files, i)
    }
    
    # 强制内存清理
    gc(verbose = FALSE)
  }
  
  return(results)
}


#' 打印进度头部
#'
#' @param current 当前序号
#' @param total 总数
#'
print_progress_header <- function(current, total) {
  cat(sprintf("\n╔═══════════════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  进度: [%d/%d] (%.1f%%)%*s║\n", 
              current, total, (current/total)*100,
              60 - nchar(sprintf("  进度: [%d/%d] (%.1f%%)", current, total, (current/total)*100)), ""))
  cat(sprintf("╚═══════════════════════════════════════════════════════════╝\n"))
}


#' 估算剩余时间
#'
#' @param results 已完成的结果
#' @param seurat_files 所有文件
#' @param current_idx 当前索引
#'
estimate_remaining_time <- function(results, seurat_files, current_idx) {
  avg_time <- mean(sapply(results, function(x) x$processing_time), na.rm = TRUE)
  remaining_time <- avg_time * (length(seurat_files) - current_idx)
  
  cat(sprintf("\n📊 预计剩余时间: %.2f 分钟 (%.2f 小时)\n", 
              remaining_time, remaining_time/60))
}


#' 确认批量处理
#'
#' @param seurat_files 文件列表
#' @param CONFIG 配置对象
#' 
#' @return 逻辑值
#'
confirm_batch_processing <- function(seurat_files, CONFIG) {
  
  if (!CONFIG$batch_mode) return(TRUE)
  if (length(seurat_files) <= 1) return(TRUE)
  if (!interactive()) return(TRUE)
  
  response <- readline(prompt = sprintf(
    "即将处理 %d 个文件，是否继续? (y/n): ", length(seurat_files)))
  
  return(tolower(response) == "y")
}


#' 创建汇总对象
#'
#' @param results 结果列表
#' @param total_elapsed 总耗时
#' @param log_files 日志文件
#' 
#' @return 汇总对象
#'
create_summary_object <- function(results, total_elapsed, log_files) {
  list(
    total = length(results),
    success = sum(sapply(results, function(x) x$success)),
    failed = sum(sapply(results, function(x) !x$success)),
    total_time = as.numeric(total_elapsed),
    log_file = log_files$log,
    csv_file = log_files$csv
  )
}

cat("✅ 10_batch_processing.R 已加载\n")
```

---

### 11_sample_preprocessing.R

- **大小**: 7.07 KB
- **修改时间**: 2025-11-07 14:18:05

```r
#!/usr/bin/env Rscript
# ===================================================================
# 样本预处理模块
# ===================================================================

#' 预处理样本（一次性切分所有样本）
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 要处理的样本列表
#' @param config 配置对象
#' 
#' @return 切分后的样本列表
#'
preprocess_samples <- function(seurat_obj, samples_to_plot, config) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   样本预处理（统一切分）\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 验证样本
  validation <- validate_samples(seurat_obj, samples_to_plot)
  samples_to_plot <- validation$valid_samples
  
  if (length(samples_to_plot) == 0) {
    stop("❌ 没有有效的样本")
  }
  
  print_preprocessing_info(seurat_obj, samples_to_plot)
  
  # 开始切分
  cat("🔧 切分样本...\n")
  start_time <- Sys.time()
  
  sample_list <- list()
  sample_stats <- initialize_sample_stats()
  
  # 切分每个样本
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    
    result <- split_single_sample(seurat_obj, sample_id, i, length(samples_to_plot))
    
    if (!is.null(result$seurat_subset)) {
      sample_list[[sample_id]] <- result$seurat_subset
      sample_stats <- rbind(sample_stats, result$stats)
    }
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  # 打印汇总
  print_preprocessing_summary(sample_list, sample_stats, elapsed, config)
  
  # 动态调整线程数
  recommended_workers <- calculate_safe_workers(sample_stats, config)
  
  # 保存统计到属性
  attr(sample_list, "stats") <- sample_stats
  attr(sample_list, "recommended_workers") <- recommended_workers
  
  return(sample_list)
}


#' 验证样本
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 样本列表
#' 
#' @return 验证结果
#'
validate_samples <- function(seurat_obj, samples_to_plot) {
  
  available_samples <- unique(seurat_obj$orig.ident)
  invalid_samples <- setdiff(samples_to_plot, available_samples)
  
  if (length(invalid_samples) > 0) {
    warning(sprintf("⚠️  以下样本不存在，将跳过: %s", 
                    paste(invalid_samples, collapse = ", ")))
  }
  
  valid_samples <- intersect(samples_to_plot, available_samples)
  
  return(list(
    valid_samples = valid_samples,
    invalid_samples = invalid_samples
  ))
}


#' 初始化样本统计表
#'
#' @return 空数据框
#'
initialize_sample_stats <- function() {
  data.frame(
    sample_id = character(),
    n_spots = integer(),
    n_high = integer(),
    high_pct = numeric(),
    size_mb = numeric(),
    stringsAsFactors = FALSE
  )
}


#' 切分单个样本
#'
#' @param seurat_obj Seurat 对象
#' @param sample_id 样本 ID
#' @param idx 当前索引
#' @param total 总数
#' 
#' @return 切分结果
#'
split_single_sample <- function(seurat_obj, sample_id, idx, total) {
  
  seurat_subset <- tryCatch({
    subset(seurat_obj, subset = orig.ident == sample_id)
  }, error = function(e) {
    seurat_obj[, seurat_obj$orig.ident == sample_id]
  })
  
  if (ncol(seurat_subset) == 0) {
    warning(sprintf("⚠️  样本 %s 无数据，已跳过", sample_id))
    return(list(seurat_subset = NULL, stats = NULL))
  }
  
  # 统计信息
  n_spots <- ncol(seurat_subset)
  n_high <- sum(seurat_subset$ClockGene_High, na.rm = TRUE)
  high_pct <- 100 * mean(seurat_subset$ClockGene_High, na.rm = TRUE)
  size_mb <- as.numeric(object.size(seurat_subset)) / 1024^2
  
  stats <- data.frame(
    sample_id = sample_id,
    n_spots = n_spots,
    n_high = n_high,
    high_pct = high_pct,
    size_mb = size_mb
  )
  
  cat(sprintf("  [%2d/%2d] ✅ %-30s | %6d spots | %4d high (%.2f%%) | %.2f MB\n",
              idx, total, sample_id, n_spots, n_high, high_pct, size_mb))
  
  return(list(seurat_subset = seurat_subset, stats = stats))
}


#' 打印预处理信息
#'
#' @param seurat_obj Seurat 对象
#' @param samples_to_plot 样本列表
#'
print_preprocessing_info <- function(seurat_obj, samples_to_plot) {
  available_samples <- unique(seurat_obj$orig.ident)
  
  cat(sprintf("📊 原始数据: %d spots, %d 个样本\n", 
              ncol(seurat_obj), length(available_samples)))
  cat(sprintf("📊 将处理: %d 个样本\n\n", length(samples_to_plot)))
}


#' 打印预处理汇总
#'
#' @param sample_list 样本列表
#' @param sample_stats 统计数据
#' @param elapsed 耗时
#' @param config 配置对象
#'
print_preprocessing_summary <- function(sample_list, sample_stats, elapsed, config) {
  
  total_spots <- sum(sample_stats$n_spots)
  total_size_mb <- sum(sample_stats$size_mb)
  avg_size_mb <- mean(sample_stats$size_mb)
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat(sprintf("✅ 成功切分 %d 个样本\n", length(sample_list)))
  cat(sprintf("📊 总计: %d spots (%.2f MB)\n", total_spots, total_size_mb))
  cat(sprintf("📊 平均: %.0f spots/样本 (%.2f MB/样本)\n", 
              total_spots / length(sample_list), avg_size_mb))
  cat(sprintf("⏱️  耗时: %.2f 秒\n", as.numeric(elapsed)))
  
  # 动态调整线程数建议
  max_memory_gb <- config$max_memory_gb %||% 100
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  recommended_workers <- min(safe_workers, length(sample_list), config$n_workers)
  
  cat(sprintf("\n💡 推荐线程数: %d (基于内存 %.0f GB)\n", 
              recommended_workers, max_memory_gb))
  
  if (recommended_workers < config$n_workers) {
    cat(sprintf("⚠️  原配置 %d 线程可能导致内存不足，已自动调整\n", 
                config$n_workers))
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 计算安全线程数
#'
#' @param sample_stats 样本统计
#' @param config 配置对象
#' 
#' @return 推荐线程数
#'
calculate_safe_workers <- function(sample_stats, config) {
  
  max_memory_gb <- config$max_memory_gb %||% 100
  avg_size_mb <- mean(sample_stats$size_mb)
  
  # 每个线程需要约 1.5 倍样本大小
  safe_workers <- floor(max_memory_gb * 1024 / (avg_size_mb * 1.5))
  
  # 不超过样本数和配置的线程数
  recommended_workers <- min(
    safe_workers, 
    nrow(sample_stats), 
    config$n_workers
  )
  
  return(max(1, recommended_workers))
}

cat("✅ 11_sample_preprocessing.R 已加载\n")
```

---

### 12_file_utils.R

- **大小**: 6.54 KB
- **修改时间**: 2025-11-07 14:18:32

```r
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
update_config_for_file <- function(seurat_path, base_config) {
  
  config <- base_config
  config$seurat_path <- seurat_path
  
  # 提取文件名
  seurat_basename <- tools::file_path_sans_ext(basename(seurat_path))
  config$output_dir <- file.path(config$output_base_dir, seurat_basename)
  
  # 更新所有目录路径
  config <- update_config_paths(config)
  
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
```

---

### 13_reporting.R

- **大小**: 10.21 KB
- **修改时间**: 2025-11-07 14:18:51

```r
#!/usr/bin/env Rscript
# ===================================================================
# 报告生成模块
# ===================================================================

#' 打印批量处理头部
#'
print_batch_header <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║        Clock Gene Niche Analysis - Batch Processing       ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n")
}


#' 打印文件处理头部
#'
#' @param seurat_basename 文件基础名
#'
print_file_header <- function(seurat_basename) {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  处理文件: %-46s ║\n", seurat_basename))
  cat("╚═══════════════════════════════════════════════════════════╝\n")
}


#' 打印文件处理成功
#'
#' @param seurat_basename 文件基础名
#' @param n_samples 样本数
#' @param elapsed 耗时
#' @param config 配置对象
#'
print_file_success <- function(seurat_basename, n_samples, elapsed, config) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    处理完成                                ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("✅ 文件: %s\n", seurat_basename))
  cat(sprintf("📊 处理样本: %d\n", n_samples))
  cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(elapsed)))
  cat(sprintf("📁 输出: %s\n", config$output_dir))
  
  print_summary(config)
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 打印文件处理失败
#'
#' @param seurat_basename 文件基础名
#' @param error_msg 错误信息
#' @param elapsed 耗时
#'
print_file_failure <- function(seurat_basename, error_msg, elapsed) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    处理失败                                ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("❌ 文件: %s\n", seurat_basename))
  cat(sprintf("❌ 错误: %s\n", error_msg))
  cat(sprintf("⏱️  耗时: %.2f 分钟\n", as.numeric(elapsed)))
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}


#' 打印批量处理总结
#'
#' @param results 结果列表
#' @param total_elapsed 总耗时
#' @param config 配置对象
#'
print_batch_summary <- function(results, total_elapsed, config) {
  
  success_count <- sum(sapply(results, function(x) x$success))
  fail_count <- length(results) - success_count
  
  cat("\n\n")
  cat("╔═══════════════════════════════════════════════════════════╗\n")
  cat("║                    批量处理总结                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("📊 总文件数: %d\n", length(results)))
  cat(sprintf("✅ 成功: %d (%.1f%%)\n", 
              success_count, (success_count/length(results))*100))
  cat(sprintf("❌ 失败: %d (%.1f%%)\n", 
              fail_count, (fail_count/length(results))*100))
  cat(sprintf("⏱️  总耗时: %.2f 分钟 (%.2f 小时)\n", 
              as.numeric(total_elapsed), as.numeric(total_elapsed)/60))
  
  if (success_count > 0) {
    print_success_statistics(results)
  }
  
  cat(sprintf("📁 输出目录: %s\n\n", config$output_base_dir))
  
  if (success_count > 0) {
    print_success_table(results)
  }
  
  if (fail_count > 0) {
    print_failure_table(results)
  }
}


#' 打印成功统计
#'
#' @param results 结果列表
#'
print_success_statistics <- function(results) {
  
  successful_results <- results[sapply(results, function(x) x$success)]
  avg_time <- mean(sapply(successful_results, function(x) x$processing_time))
  total_samples <- sum(sapply(successful_results, function(x) x$n_samples))
  
  cat(sprintf("📈 平均耗时: %.2f 分钟/文件\n", avg_time))
  cat(sprintf("📊 总样本数: %d\n", total_samples))
}


#' 打印成功文件表格
#'
#' @param results 结果列表
#'
print_success_table <- function(results) {
  
  cat("✅ 成功处理的文件:\n")
  cat(sprintf("%-4s %-40s %10s %10s\n", "No.", "文件名", "耗时(分)", "样本数"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  j <- 1
  for (i in seq_along(results)) {
    if (results[[i]]$success) {
      cat(sprintf("%3d. %-40s %10.2f %10d\n", 
                  j,
                  results[[i]]$file,
                  results[[i]]$processing_time,
                  results[[i]]$n_samples))
      j <- j + 1
    }
  }
  cat("\n")
}


#' 打印失败文件表格
#'
#' @param results 结果列表
#'
print_failure_table <- function(results) {
  
  cat("❌ 失败的文件:\n")
  cat(sprintf("%-4s %-40s %s\n", "No.", "文件名", "错误信息"))
  cat(paste(rep("-", 100), collapse = ""), "\n")
  
  j <- 1
  for (i in seq_along(results)) {
    if (!results[[i]]$success) {
      cat(sprintf("%3d. %-40s %s\n", 
                  j,
                  results[[i]]$file,
                  substr(results[[i]]$error, 1, 50)))
      j <- j + 1
    }
  }
  cat("\n")
}


#' 保存批量处理日志
#'
#' @param results 结果列表
#' @param start_time 开始时间
#' @param end_time 结束时间
#' @param config 配置对象
#' 
#' @return 日志文件路径
#'
save_batch_logs <- function(results, start_time, end_time, config) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  log_file <- file.path(config$output_base_dir, 
                       sprintf("batch_processing_log_%s.txt", timestamp))
  
  csv_file <- file.path(config$output_base_dir, 
                       sprintf("batch_summary_%s.csv", timestamp))
  
  # 保存文本日志
  save_text_log(results, start_time, end_time, log_file)
  
  # 保存 CSV
  save_csv_summary(results, csv_file)
  
  return(list(log = log_file, csv = csv_file))
}


#' 保存文本日志
#'
#' @param results 结果列表
#' @param start_time 开始时间
#' @param end_time 结束时间
#' @param log_file 日志文件路径
#'
save_text_log <- function(results, start_time, end_time, log_file) {
  
  tryCatch({
    sink(log_file)
    
    cat("═══════════════════════════════════════════════════════════\n")
    cat("           Batch Processing Log\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    total_time <- difftime(end_time, start_time, units = "mins")
    
    cat(sprintf("Start time: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("End time:   %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Total time: %.2f minutes (%.2f hours)\n\n", 
                as.numeric(total_time), as.numeric(total_time)/60))
    
    # 详细结果
    for (i in seq_along(results)) {
      result <- results[[i]]
      status <- if(result$success) "SUCCESS" else "FAILED"
      
      cat(sprintf("[%s] File %2d/%d: %s\n", 
                  status, i, length(results), result$file))
      
      if (result$success) {
        cat(sprintf("           Time: %.2f min, Samples: %d\n", 
                    result$processing_time, result$n_samples))
      } else {
        cat(sprintf("           Error: %s\n", result$error))
      }
      cat("\n")
    }
    
    sink()
    
    cat(sprintf("📝 日志已保存:\n   %s\n", log_file))
    
  }, error = function(e) {
    sink()
    warning(sprintf("⚠️  无法保存日志: %s", e$message))
  })
}


#' 保存 CSV 汇总
#'
#' @param results 结果列表
#' @param csv_file CSV 文件路径
#'
save_csv_summary <- function(results, csv_file) {
  
  tryCatch({
    summary_df <- data.frame(
      File_Number = seq_along(results),
      File_Name = sapply(results, function(x) x$file),
      Status = sapply(results, function(x) ifelse(x$success, "Success", "Failed")),
      Processing_Time_Minutes = sapply(results, function(x) round(x$processing_time, 2)),
      Number_of_Samples = sapply(results, function(x) x$n_samples),
      Error_Message = sapply(results, function(x) ifelse(!x$success, x$error, "")),
      stringsAsFactors = FALSE
    )
    
    write.csv(summary_df, csv_file, row.names = FALSE, quote = TRUE)
    cat(sprintf("📊 CSV已保存:\n   %s\n\n", csv_file))
    
  }, error = function(e) {
    warning(sprintf("⚠️  无法保存CSV: %s", e$message))
  })
}

cat("✅ 13_reporting.R 已加载\n")
```

---

### AI_trans.R

- **大小**: 6.52 KB
- **修改时间**: 2025-11-07 14:31:56

```r

#!/usr/bin/env Rscript
# 项目信息导出脚本
# 将当前目录的项目结构和内容导出为Markdown文件

# 加载必要的包
if (!require("tools", quietly = TRUE)) install.packages("tools")

# 配置参数
output_file <- "project_export.md"
max_file_size <- 1024 * 1024  # 1MB, 超过此大小的文件不读取内容
ignore_dirs <- c(".git", ".Rproj.user", "node_modules", "__pycache__", ".venv", "venv")
ignore_files <- c(".DS_Store", "Thumbs.db", ".gitignore")
code_extensions <- c("R", "r", "py", "js", "html", "css", "java", "cpp", "c", "h", 
                     "sh", "sql", "md", "rmd", "yml", "yaml", "json", "xml", "txt")

# 函数：获取文件大小的可读格式
format_file_size <- function(size) {
  if (size < 1024) {
    return(paste0(size, " B"))
  } else if (size < 1024^2) {
    return(paste0(round(size / 1024, 2), " KB"))
  } else if (size < 1024^3) {
    return(paste0(round(size / 1024^2, 2), " MB"))
  } else {
    return(paste0(round(size / 1024^3, 2), " GB"))
  }
}

# 函数：生成目录树
generate_tree <- function(path, prefix = "", is_last = TRUE) {
  tree_lines <- c()
  
  files <- list.files(path, all.files = FALSE, include.dirs = TRUE, no.. = TRUE)
  files <- files[!files %in% ignore_dirs & !files %in% ignore_files]
  
  if (length(files) == 0) return(tree_lines)
  
  for (i in seq_along(files)) {
    file_path <- file.path(path, files[i])
    is_last_item <- (i == length(files))
    
    connector <- if (is_last_item) "└── " else "├── "
    tree_lines <- c(tree_lines, paste0(prefix, connector, files[i]))
    
    if (dir.exists(file_path)) {
      new_prefix <- paste0(prefix, if (is_last_item) "    " else "│   ")
      tree_lines <- c(tree_lines, generate_tree(file_path, new_prefix, is_last_item))
    }
  }
  
  return(tree_lines)
}

# 函数：递归获取所有文件
get_all_files <- function(path) {
  all_files <- c()
  
  files <- list.files(path, all.files = FALSE, include.dirs = TRUE, 
                     full.names = TRUE, no.. = TRUE)
  
  for (file_path in files) {
    file_name <- basename(file_path)
    
    # 跳过忽略的目录和文件
    if (file_name %in% ignore_dirs || file_name %in% ignore_files) next
    
    if (dir.exists(file_path)) {
      all_files <- c(all_files, get_all_files(file_path))
    } else {
      all_files <- c(all_files, file_path)
    }
  }
  
  return(all_files)
}

# 函数：判断是否为文本文件
is_text_file <- function(file_path) {
  ext <- tools::file_ext(file_path)
  return(tolower(ext) %in% tolower(code_extensions))
}

# 主函数
export_project <- function() {
  cat("开始导出项目信息...\n")
  
  # 创建Markdown内容
  md_content <- c()
  
  # 标题和基本信息
  project_name <- basename(getwd())
  md_content <- c(md_content, paste0("# 项目导出: ", project_name))
  md_content <- c(md_content, "")
  md_content <- c(md_content, paste0("**导出时间**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  md_content <- c(md_content, paste0("**项目路径**: ", getwd()))
  md_content <- c(md_content, "")
  md_content <- c(md_content, "---")
  md_content <- c(md_content, "")
  
  # 目录结构
  cat("生成目录树...\n")
  md_content <- c(md_content, "## 目录结构")
  md_content <- c(md_content, "")
  md_content <- c(md_content, "```")
  md_content <- c(md_content, project_name)
  tree_lines <- generate_tree(getwd())
  md_content <- c(md_content, tree_lines)
  md_content <- c(md_content, "```")
  md_content <- c(md_content, "")
  
  # 获取所有文件
  cat("收集文件信息...\n")
  all_files <- get_all_files(getwd())
  
  # 文件统计
  md_content <- c(md_content, "## 文件统计")
  md_content <- c(md_content, "")
  md_content <- c(md_content, paste0("- **总文件数**: ", length(all_files)))
  
  # 按扩展名统计
  extensions <- sapply(all_files, tools::file_ext)
  ext_table <- table(extensions)
  ext_table <- sort(ext_table, decreasing = TRUE)
  
  md_content <- c(md_content, "- **文件类型分布**:")
  for (i in seq_along(ext_table)) {
    ext_name <- if (names(ext_table)[i] == "") "(无扩展名)" else paste0(".", names(ext_table)[i])
    md_content <- c(md_content, paste0("  - ", ext_name, ": ", ext_table[i], " 个"))
  }
  
  # 总大小
  total_size <- sum(sapply(all_files, file.size))
  md_content <- c(md_content, paste0("- **项目总大小**: ", format_file_size(total_size)))
  md_content <- c(md_content, "")
  
  # 文件详细内容
  md_content <- c(md_content, "## 文件内容")
  md_content <- c(md_content, "")
  
  cat("读取文件内容...\n")
  for (file_path in all_files) {
    relative_path <- sub(paste0("^", getwd(), "/"), "", file_path)
    file_size <- file.size(file_path)
    
    md_content <- c(md_content, paste0("### ", relative_path))
    md_content <- c(md_content, "")
    md_content <- c(md_content, paste0("- **大小**: ", format_file_size(file_size)))
    md_content <- c(md_content, paste0("- **修改时间**: ", format(file.mtime(file_path), "%Y-%m-%d %H:%M:%S")))
    md_content <- c(md_content, "")
    
    # 如果是文本文件且大小合适，读取内容
    if (is_text_file(file_path) && file_size <= max_file_size) {
      tryCatch({
        content <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
        ext <- tools::file_ext(file_path)
        lang <- if (ext == "") "" else tolower(ext)
        
        md_content <- c(md_content, paste0("```", lang))
        md_content <- c(md_content, content)
        md_content <- c(md_content, "```")
        md_content <- c(md_content, "")
      }, error = function(e) {
        md_content <<- c(md_content, "*无法读取文件内容*")
        md_content <<- c(md_content, "")
      })
    } else if (file_size > max_file_size) {
      md_content <- c(md_content, "*文件过大，跳过内容显示*")
      md_content <- c(md_content, "")
    } else {
      md_content <- c(md_content, "*二进制文件，跳过内容显示*")
      md_content <- c(md_content, "")
    }
    
    md_content <- c(md_content, "---")
    md_content <- c(md_content, "")
  }
  
  # 写入文件
  cat(paste0("写入文件: ", output_file, "\n"))
  writeLines(md_content, output_file, useBytes = TRUE)
  
  cat(paste0("✓ 导出完成！文件已保存到: ", output_file, "\n"))
  cat(paste0("  共处理 ", length(all_files), " 个文件\n"))
}

# 执行导出
tryCatch({
  export_project()
}, error = function(e) {
  cat("错误:", conditionMessage(e), "\n")
})
```

---

### main.R

- **大小**: 7.94 KB
- **修改时间**: 2025-11-07 14:17:20

```r
#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis - Main Script (Optimized)
# Author: Zhangbin
# Date: 2024-11-04
# Optimized: 2024-11-07
#   - 模块化拆分
#   - 统一环境初始化
#   - 内存管理优化
# ===================================================================

# ===================================================================
# 加载配置和模块
# ===================================================================

source("00_config.R")
source("01_setup.R")
source("02_utils.R")
source("03_load_data.R")
source("04_module_score.R")
source("05_niche_analysis.R")
source("06_plot_isoheight.R")
source("07_plot_spatial.R")
source("08_plot_celltype.R")
source("09_save_results.R")
source("10_batch_processing.R")       # 批量处理
source("11_sample_preprocessing.R")   # 样本预处理
source("12_file_utils.R")             # 文件工具
source("13_reporting.R")              # 报告生成


# ===================================================================
# 单文件处理函数（核心流程）
# ===================================================================

#' 处理单个 Seurat 文件
#'
#' @param seurat_path Seurat 文件路径
#' @param gene_list 基因列表
#' @param base_config 基础配置
#' 
#' @return 处理结果列表
#'
process_seurat_file <- function(seurat_path, gene_list, base_config) {
  
  # 1. 更新配置
  config <- update_config_for_file(seurat_path, base_config)
  seurat_basename <- tools::file_path_sans_ext(basename(seurat_path))
  
  # 2. 打印处理信息
  print_file_header(seurat_basename)
  
  file_start_time <- Sys.time()
  
  tryCatch({
    
    # ----------------------------------------
    # 步骤 1-5: 数据准备和分析
    # ----------------------------------------
    cat("\n【步骤 1/9】环境设置\n")
    setup_environment(config)
    
    cat("\n【步骤 2/9】加载 Seurat 对象\n")
    seurat_obj <- load_seurat_object(config)
    genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
    
    cat("\n【步骤 3/9】计算 Clock Gene Score\n")
    seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, config)
    
    cat("\n【步骤 4/9】识别高表达区域\n")
    result <- define_high_expression(seurat_obj, config)
    seurat_obj <- result$seurat_obj
    threshold <- result$threshold
    
    cat("\n【步骤 5/9】Niche 分析\n")
    seurat_obj <- perform_niche_analysis(seurat_obj, threshold, config)
    
    # ----------------------------------------
    # 步骤 6: 样本预处理（统一切分）
    # ----------------------------------------
    cat("\n【步骤 6/9】样本预处理\n")
    
    samples <- unique(seurat_obj$orig.ident)
    samples_to_plot <- if (config$debug_mode) {
      head(samples, config$debug_sample_limit %||% 3)
    } else {
      samples
    }
    
    # 一次性切分所有样本
    sample_list <- preprocess_samples(seurat_obj, samples_to_plot, config)
    
    # 更新配置中的线程数（基于内存估算）
    recommended_workers <- attr(sample_list, "recommended_workers")
    if (!is.null(recommended_workers)) {
      config$n_workers <- recommended_workers
    }
    
    # ----------------------------------------
    # 步骤 7-9: 可视化分析
    # ----------------------------------------
    cat("\n【步骤 7/9】绘制等高线密度图\n")
    iso_results <- plot_isoheight(
      sample_list = sample_list,
      CONFIG = config
    )
    
    cat("\n【步骤 8/9】绘制空间梯度图\n")
    spatial_results <- plot_spatial_gradient(
      sample_list = sample_list,
      CONFIG = config
    )
    
    cat("\n【步骤 9/9】细胞类型 Niche 分析\n")
    celltype_results <- analyze_celltype_niche(
      sample_list = sample_list,
      CONFIG = config,
      seurat_basename = seurat_basename
    )
    
    # ----------------------------------------
    # 保存结果
    # ----------------------------------------
    save_results(seurat_obj, config)
    
    # ----------------------------------------
    # 完成
    # ----------------------------------------
    file_end_time <- Sys.time()
    file_elapsed <- difftime(file_end_time, file_start_time, units = "mins")
    
    print_file_success(seurat_basename, length(sample_list), file_elapsed, config)
    
    # 清理内存
    rm(seurat_obj, sample_list)
    gc(verbose = FALSE)
    
    return(list(
      success = TRUE,
      file = seurat_basename,
      processing_time = as.numeric(file_elapsed),
      n_samples = length(samples_to_plot),
      error = NULL
    ))
    
  }, error = function(e) {
    
    file_end_time <- Sys.time()
    file_elapsed <- difftime(file_end_time, file_start_time, units = "mins")
    
    print_file_failure(seurat_basename, e$message, file_elapsed)
    
    # 清理内存
    gc(verbose = FALSE)
    
    return(list(
      success = FALSE,
      file = seurat_basename,
      processing_time = as.numeric(file_elapsed),
      n_samples = 0,
      error = e$message
    ))
  })
}


# ===================================================================
# 批量处理主函数（简化版）
# ===================================================================

#' 批量处理主函数
#'
#' @return 批量处理结果
#'
main_batch <- function() {
  
  batch_start_time <- Sys.time()
  
  print_batch_header()
  
  # ----------------------------------------
  # 1. 统一初始化环境
  # ----------------------------------------
  cat("\n【初始化】环境设置\n")
  
  init_result <- initialize_environment(
    config = CONFIG,
    custom_scripts = c("niche_marker.R", "SSS_isoheight_plot.R")
  )
  
  if (length(init_result$packages$failed) > 0) {
    warning("⚠️  部分包加载失败，可能影响分析")
  }
  
  # ----------------------------------------
  # 2. 验证输出目录
  # ----------------------------------------
  validate_output_directory(CONFIG)
  
  # ----------------------------------------
  # 3. 扫描输入文件
  # ----------------------------------------
  seurat_files <- scan_seurat_files(CONFIG)
  
  if (length(seurat_files) == 0) {
    stop("❌ 未找到可处理的文件")
  }
  
  print_file_list(seurat_files)
  
  # 确认处理
  if (!confirm_batch_processing(seurat_files, CONFIG)) {
    cat("❌ 已取消处理\n")
    return(invisible(NULL))
  }
  
  # ----------------------------------------
  # 4. 加载基因列表（只加载一次）
  # ----------------------------------------
  gene_list <- load_gene_list_once(CONFIG)
  
  # ----------------------------------------
  # 5. 批量处理文件
  # ----------------------------------------
  results <- process_all_files(seurat_files, gene_list, CONFIG)
  
  # ----------------------------------------
  # 6. 生成总结报告
  # ----------------------------------------
  batch_end_time <- Sys.time()
  total_elapsed <- difftime(batch_end_time, batch_start_time, units = "mins")
  
  print_batch_summary(results, total_elapsed, CONFIG)
  
  log_files <- save_batch_logs(results, batch_start_time, batch_end_time, CONFIG)
  
  cat("\n🎉 批量处理完成！\n\n")
  
  return(invisible(list(
    results = results,
    summary = create_summary_object(results, total_elapsed, log_files)
  )))
}


# ===================================================================
# 辅助操作符
# ===================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}


# ===================================================================
# 运行主流程
# ===================================================================

if (!interactive()) {
  main_batch()
}

cat("✅ main.R 已加载\n")
cat("📚 使用 main_batch() 开始批量处理\n\n")
```

---

### niche_grade_entropy.R

- **大小**: 5.24 KB
- **修改时间**: 2025-10-30 11:38:58

```r
# niche_grade_entropy.R

GetAllCoordinates <- function(.data) {
    .data@images %>%
        names() %>%
        unique() %>%
        map_dfr(~{
            GetTissueCoordinates(
                    .data,
                    image = gsub('-', '.', .x),
                    cols = c("row", "col"),
                    scale = NULL
                ) %>%
            rownames_to_column(var = "cellid")
        })
}

# ################################################################
#
# Neighbors Count
#
# ################################################################

image_spot_neighbors_count <- function(meta.data, celltype_col, neighbor_range, ...) {
    .sample_coord <- meta.data %>%
        select(row,col,!!celltype_col, ...) %>%
        mutate(.celltype = as.numeric(factor(!!celltype_col)))
    # print(head(.sample_coord))
        
    # build celltype matrix
    .celltype_mat <- Matrix::sparseMatrix(
        .sample_coord$row,
        .sample_coord$col,
        x = .sample_coord$.celltype)

    row_r <- dim(.celltype_mat)[1]
    col_r <- dim(.celltype_mat)[2]
    print(str_c("row range: ", row_r, " col range: ", col_r))

    # count
    .sample_coord %>%
        mutate(
            neighbors = map2(row, col, ~{
                .celltype_code <- .celltype_mat[
                        max(0, .x - neighbor_range):min(.x + neighbor_range, dim(.celltype_mat)[1]),
                        max(0, .y - neighbor_range):min(.y + neighbor_range, dim(.celltype_mat)[2])
                    ] %>%
                    as.vector()

                .celltype_code <- .celltype_code[which(.celltype_code > 0)]
                .celltype_code <- .celltype_code[-match(.celltype_mat[.x,.y], .celltype_code)]

                table(.celltype_code)
            })
        )
}

# ################################################################
#
# Neighbors Count
#
# ################################################################


#' SSS niche gradient entropy definition
#'
#' @param .data tibble obj
#' @param celltype_col column in meta.data, use for classifying spots
#' @param ...  columns in meta.data, and will be reserve in result
#' @param roi_col column in meta.data, use for classifying niche gradient
#' @param neighbor_range neighborhood range
#' @param R resample times
#' @param n_work number of threads
#'
niche_grade_entropy <- function(
    .data, 
    ..., 
    celltype_col = seurat_clusters, 
    neighbor_range = 1, 
    roi_col, 
    R = 100, 
    n_work = 3
) {

    celltype_col = enquo(celltype_col)
    roi_col = enquo(roi_col)
    group_vars <- enquos(..., .named = TRUE)

    library(future)
    library(future.apply)

    plan(multisession, workers=n_work)
    options(future.globals.maxSize= Inf)
    options(future.seed=TRUE)
    message("outside >> how many cores can use now: ", nbrOfWorkers())

    set.seed(2023)

    df <- .data@meta.data %>%
        as_tibble(rownames = "cellid") 
    # Get Coordinates
    if(!("col" %in% colnames(df) && "row" %in% colnames(df))) {
        df <- df %>%
            left_join(
                GetAllCoordinates(.data)
            )
    }
    # spot2niche
    df <- df %>%
        group_by(age, orig.ident) %>%
        group_nest() %>%
        mutate(
            data = future_lapply(data, function(df) {
                image_spot_neighbors_count(
                    df,
                    celltype_col,
                    neighbor_range,
                    cellid, !!roi_col
                )
            }, future.chunk.size = Inf)
        ) %>%
        unnest(data) %>%
        filter(!is.na(!!roi_col)) %>%
        group_by(age, !!roi_col) %>%
        mutate(min_roi_spot_num = n())
    # get roi_min_spot_num
    df$min_roi_spot_num <- min(df$min_roi_spot_num)

    # get entropy
    df <- df %>%
        group_by(age, !!roi_col, min_roi_spot_num) %>%
        group_nest() %>%
        mutate(
            data = map2(data, min_roi_spot_num, ~{
                message("")
                message(str_c("## sample spot num: ", .y))

                data <- .x
                min_roi_spot_num <- .y
                # resample
                res <- seq(R) %>%
                    as.list() %>%
                    future_lapply(function(idx) {
                        data %>%
                            slice_sample(n = min_roi_spot_num) %>%
                            group_by(!!celltype_col, .celltype, neighbors) %>%
                            summarise(f_ij = n()) %>%
                            ungroup() %>%
                            mutate(
                                p_ij = f_ij / sum(f_ij),
                                entropy_ij = - p_ij * log2(p_ij),
                                condition = n()
                            ) %>%
                            group_by(condition) %>%
                            summarise(entropy = sum(entropy_ij)) %>%
                            mutate(
                                rep_idx = idx,
                                Pielou = entropy / log(condition))
                    }, future.chunk.size = Inf, future.seed = TRUE) %>%
                    bind_rows() 
                    
            })
        ) %>%
        unnest(data)
}

```

---

### niche_marker.R

- **大小**: 15.37 KB
- **修改时间**: 2025-11-07 12:33:24

```r
# ========================================================================
# niche_marker.R - 完整修复版
# 移除 GetTissueCoordinates 依赖，直接从 @images 提取坐标
# ========================================================================

library(dplyr)
library(tibble)
library(purrr)
library(proxy)
library(future)
library(future.apply)


# ========================================================================
# 1. 坐标提取函数（主方法）
# ========================================================================

GetAllCoordinates <- function(.data) {
  cat("\n🔍 提取空间坐标（直接从 @images 读取）...\n")
  
  image_names <- names(.data@images)
  
  if (length(image_names) == 0) {
    stop("❌ Seurat 对象中没有空间图像数据")
  }
  
  cat(sprintf(">> 发现 %d 个图像\n", length(image_names)))
  
  # 逐个提取坐标
  coords_list <- list()
  
  for (img_name in image_names) {
    cat(sprintf("  >> 提取: %s ... ", img_name))
    
    tryCatch({
      # 直接从 coordinates 槽提取
      img_obj <- .data@images[[img_name]]
      coords_df <- img_obj@coordinates
      
      # 转换为 data.frame
      if (!is.data.frame(coords_df)) {
        coords_df <- as.data.frame(coords_df)
      }
      
      # 获取细胞ID
      cell_ids <- rownames(coords_df)
      if (is.null(cell_ids) || all(is.na(cell_ids))) {
        stop("坐标数据缺少行名（细胞ID）")
      }
      
      # 查找坐标列
      col_names <- colnames(coords_df)
      
      # 可能的列名
      row_candidates <- c("row", "imagerow", "array_row", "tissue_row", "pxl_row_in_fullres")
      col_candidates <- c("col", "imagecol", "array_col", "tissue_col", "pxl_col_in_fullres")
      
      # 找到实际的列名
      row_col <- intersect(col_names, row_candidates)
      col_col <- intersect(col_names, col_candidates)
      
      if (length(row_col) == 0 || length(col_col) == 0) {
        stop(sprintf("未找到坐标列。可用列: %s", paste(col_names, collapse=", ")))
      }
      
      # 使用第一个匹配的列名
      row_col_name <- row_col[1]
      col_col_name <- col_col[1]
      
      # 提取坐标
      result <- data.frame(
        cellid = cell_ids,
        row = as.numeric(coords_df[[row_col_name]]),
        col = as.numeric(coords_df[[col_col_name]]),
        stringsAsFactors = FALSE
      )
      
      # 检查 NA
      n_na <- sum(is.na(result$row) | is.na(result$col))
      if (n_na > 0) {
        warning(sprintf("%s: %d 个细胞的坐标为 NA", img_name, n_na))
      }
      
      coords_list[[img_name]] <- result
      cat(sprintf("✓ %d 个细胞\n", nrow(result)))
      
    }, error = function(e) {
      cat(sprintf("❌ 失败: %s\n", e$message))
      stop(sprintf("样本 %s 坐标提取失败", img_name))
    })
  }
  
  # 合并所有坐标
  all_coords <- bind_rows(coords_list)
  
  if (nrow(all_coords) == 0) {
    stop("❌ 未能提取任何坐标数据")
  }
  
  cat(sprintf("✅ 总共提取 %d 个细胞的坐标\n\n", nrow(all_coords)))
  
  return(all_coords)
}


# ========================================================================
# 2. 单个样本的距离计算
# ========================================================================

single_marker <- function(df, intra_df, spot_type, dist_method, FUN) {
  
  if (nrow(intra_df) > 0) {
    # 准备所有细胞的坐标
    all_df <- df %>%
      column_to_rownames("cellid") %>%
      select(row, col)

    cat(sprintf("  计算距离矩阵: %d 个查询点 × %d 个标记点\n", 
                nrow(all_df), nrow(intra_df)))

    # 计算距离矩阵
    mat <- proxy::dist(all_df, intra_df, method = dist_method) %>%
      as.matrix()

    # 计算每个细胞到最近标记点的距离
    spot_dist <- tibble(cellid = rownames(mat))
    
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      spot_dist[[spot_type]] <- matrixStats::rowMins(mat, na.rm = TRUE)
    } else {
      spot_dist[[spot_type]] <- apply(mat, 1, min, na.rm = TRUE)
    }

    # 应用转换函数（如果提供）
    if (!is.na(FUN)) {
      spot_dist[[spot_type]] <- FUN(spot_dist[[spot_type]])
    }

    # 合并回原始数据
    res <- df %>%
      left_join(spot_dist, by = "cellid")

  } else {
    # 没有标记点，所有距离设为 Inf
    cat("  ⚠️ 警告：没有找到标记点，Distance 设为 Inf\n")
    res <- df %>%
      mutate(!!spot_type := Inf)
  }

  # 移除坐标列
  res %>% select(-c(row, col))
}


# ========================================================================
# 3. 主函数：Niche Marker 分析
# ========================================================================

niche_marker <- function(
  .data,
  marker,
  spot_type,
  slide = "orig.ident",
  dist_method = "Euclidean",
  FUN = NA,
  n_work = 3
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   开始 Niche Marker 分析\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 转换参数为字符串
  marker <- as.character(substitute(marker))
  spot_type <- as.character(substitute(spot_type))
  slide <- as.character(substitute(slide))
  
  cat(sprintf("参数配置:\n"))
  cat(sprintf("  Marker 列: %s\n", marker))
  cat(sprintf("  输出列: %s\n", spot_type))
  cat(sprintf("  样本列: %s\n", slide))
  cat(sprintf("  距离方法: %s\n", dist_method))
  cat(sprintf("  工作线程: %d\n\n", n_work))

  # 加载必要的包
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("需要安装 future 包: install.packages('future')")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("需要安装 future.apply 包: install.packages('future.apply')")
  }

  library(future)
  library(future.apply)

  # 设置并行计算
  plan(multisession, workers = n_work)
  options(future.globals.maxSize = Inf)
  cat(sprintf(">> 并行核心数: %d\n\n", nbrOfWorkers()))

  # 数据统计
  n_total <- ncol(.data)
  n_marker <- sum(.data@meta.data[[marker]], na.rm = TRUE)
  cat(sprintf("数据概况:\n"))
  cat(sprintf("  总细胞数: %d\n", n_total))
  cat(sprintf("  标记细胞: %d (%.1f%%)\n\n", n_marker, 100 * n_marker / n_total))

  # 保存原始细胞顺序（关键！）
  original_cell_order <- colnames(.data)
  cat(sprintf(">> 保存原始细胞顺序: %d 个细胞\n\n", length(original_cell_order)))

  # ========== 提取空间坐标 ==========
  cat("🔄 提取空间坐标...\n")
  all_coords <- tryCatch({
    GetAllCoordinates(.data)
  }, error = function(e) {
    stop(sprintf("❌ 坐标提取失败: %s", e$message))
  })
  
  # 验证坐标完整性
  if (nrow(all_coords) != n_total) {
    stop(sprintf("❌ 坐标数量 (%d) 与细胞数量 (%d) 不匹配", 
                 nrow(all_coords), n_total))
  }
  
  # 检查是否所有细胞都有坐标
  missing_cells <- setdiff(original_cell_order, all_coords$cellid)
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ %d 个细胞缺少坐标数据", length(missing_cells)))
  }

  # ========== 合并 metadata 和坐标 ==========
  cat("\n🔄 合并 metadata 和坐标...\n")
  meta_with_coords <- .data@meta.data %>%
    rownames_to_column(var = "cellid") %>%
    left_join(all_coords, by = "cellid")
  
  # 验证合并结果
  n_missing_coords <- sum(is.na(meta_with_coords$row) | is.na(meta_with_coords$col))
  if (n_missing_coords > 0) {
    stop(sprintf("❌ %d 个细胞在合并后缺少坐标", n_missing_coords))
  }
  cat("✅ 合并完成\n")

  # ========== 按样本分组并计算距离 ==========
  cat("\n🔄 按样本计算距离...\n")
  
  # 分组
  sample_groups <- meta_with_coords %>%
    group_by(.data[[slide]]) %>%
    group_split()
  
  cat(sprintf(">> 将处理 %d 个样本\n\n", length(sample_groups)))

  # 并行处理每个样本
  results_list <- future_lapply(sample_groups, function(df) {
    
    slide_name <- df[[slide]][1]
    cat(sprintf("处理样本: %s\n", slide_name))
    
    # 提取标记点
    intra_df <- df %>%
      filter(!is.na(.data[[marker]]) & .data[[marker]] == TRUE) %>%
      column_to_rownames("cellid") %>%
      select(row, col)
    
    n_sample <- nrow(df)
    n_marker_sample <- nrow(intra_df)
    
    cat(sprintf("  样本细胞数: %d\n", n_sample))
    cat(sprintf("  标记细胞数: %d (%.1f%%)\n", 
                n_marker_sample, 100 * n_marker_sample / n_sample))
    
    # 计算距离
    result <- single_marker(
      df = df, 
      intra_df = intra_df, 
      spot_type = spot_type,
      dist_method = dist_method, 
      FUN = FUN
    )
    
    return(result)
    
  }, future.seed = TRUE, future.chunk.size = Inf)
  
  cat("\n🔄 合并所有样本结果...\n")

  # 合并结果
  combined_results <- bind_rows(results_list)
  
  # 将结果转换为以 cellid 为行名的 data.frame
  combined_results <- combined_results %>%
    column_to_rownames(var = "cellid")

  # ========== 恢复原始细胞顺序 ==========
  cat("\n🔄 恢复原始细胞顺序...\n")
  
  current_cells <- rownames(combined_results)
  missing_cells <- setdiff(original_cell_order, current_cells)
  extra_cells <- setdiff(current_cells, original_cell_order)
  
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ 结果中缺少 %d 个细胞！", length(missing_cells)))
  }
  
  if (length(extra_cells) > 0) {
    warning(sprintf("⚠️ 结果中有 %d 个多余细胞，将被移除", length(extra_cells)))
    combined_results <- combined_results[original_cell_order, ]
  } else {
    # 按原始顺序重新排列
    combined_results <- combined_results[original_cell_order, ]
  }
  
  # 验证顺序
  if (!identical(rownames(combined_results), original_cell_order)) {
    stop("❌ 严重错误：细胞顺序恢复失败！")
  }
  
  cat("✅ 细胞顺序已恢复并验证\n")

  # ========== 将结果添加到 Seurat 对象 ==========
  .data@meta.data <- combined_results

  # ========== 结果验证 ==========
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   结果验证\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  dist_vals <- .data@meta.data[[spot_type]]
  
  cat(sprintf("Distance 统计:\n"))
  cat(sprintf("  最小值: %.2f\n", min(dist_vals, na.rm = TRUE)))
  cat(sprintf("  最大值: %.2f\n", max(dist_vals, na.rm = TRUE)))
  cat(sprintf("  平均值: %.2f\n", mean(dist_vals, na.rm = TRUE)))
  cat(sprintf("  中位数: %.2f\n", median(dist_vals, na.rm = TRUE)))
  cat(sprintf("  NA 数量: %d\n\n", sum(is.na(dist_vals))))

  # 验证标记点的距离
  marker_cells <- .data@meta.data[[marker]]
  marker_dist <- dist_vals[marker_cells]
  
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- sum(!is.na(marker_dist))
  pct_zero <- 100 * n_marker_zero / n_marker_total
  
  cat(sprintf("标记细胞验证:\n"))
  cat(sprintf("  标记细胞总数: %d\n", n_marker_total))
  cat(sprintf("  Distance=0: %d (%.1f%%)\n", n_marker_zero, pct_zero))
  cat(sprintf("  Distance>0: %d (%.1f%%)\n", 
              n_marker_total - n_marker_zero, 
              100 - pct_zero))
  
  if (pct_zero < 95) {
    warning(sprintf("⚠️ 警告：只有 %.1f%% 的标记细胞 Distance=0！预期应接近 100%%", pct_zero))
    
    # 显示异常的标记细胞
    abnormal <- which(marker_dist > 0.1)
    if (length(abnormal) > 0) {
      cat(sprintf("\n前 5 个异常标记细胞:\n"))
      print(head(marker_dist[abnormal], 5))
    }
  } else {
    cat("\n✅ 验证通过：几乎所有标记细胞的 Distance = 0\n")
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("   Niche Marker 分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 关闭并行
  plan(sequential)

  return(.data)
}


# ========================================================================
# 4. 辅助函数：诊断空间坐标
# ========================================================================

diagnose_spatial_coordinates <- function(.data) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("           空间坐标诊断报告\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  image_names <- names(.data@images)
  
  if (length(image_names) == 0) {
    cat("❌ 未找到空间图像数据\n\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("总图像数: %d\n\n", length(image_names)))
  
  for (i in seq_along(image_names)) {
    img_name <- image_names[i]
    cat(sprintf("[%d/%d] 图像: %s\n", i, length(image_names), img_name))
    cat("─────────────────────────────────────────────────────────\n")
    
    img_obj <- .data@images[[img_name]]
    coords <- img_obj@coordinates
    
    cat(sprintf("   细胞数: %d\n", nrow(coords)))
    cat(sprintf("   坐标列: %s\n", paste(colnames(coords), collapse=", ")))
    
    # 检查标准列
    has_row <- "row" %in% colnames(coords)
    has_col <- "col" %in% colnames(coords)
    
    cat(sprintf("   标准列: row=%s, col=%s\n", 
                ifelse(has_row, "✓", "✗"),
                ifelse(has_col, "✓", "✗")))
    
    # 如果有坐标，显示范围
    row_col <- intersect(colnames(coords), 
                        c("row", "imagerow", "array_row", "tissue_row"))
    col_col <- intersect(colnames(coords), 
                        c("col", "imagecol", "array_col", "tissue_col"))
    
    if (length(row_col) > 0 && length(col_col) > 0) {
      cat(sprintf("   坐标范围:\n"))
      cat(sprintf("      %s: [%.1f, %.1f]\n", 
                  row_col[1],
                  min(coords[[row_col[1]]], na.rm=TRUE), 
                  max(coords[[row_col[1]]], na.rm=TRUE)))
      cat(sprintf("      %s: [%.1f, %.1f]\n", 
                  col_col[1],
                  min(coords[[col_col[1]]], na.rm=TRUE), 
                  max(coords[[col_col[1]]], na.rm=TRUE)))
    }
    
    cat("\n")
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
}
```

---

### SSS_isoheight_plot.R

- **大小**: 3.93 KB
- **修改时间**: 2025-11-04 15:14:57

```r
# SSS_isoheight_plot.R (完整版)

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
            rownames_to_column(var = "cellid")
        })
}

celltype_isoheight_plot <- function(
    .data,
    density_top,
    col_bg = "gray90",
    col_top = "darkred",
    col_isoheight = "white",
    col_white_ratio = 0.2,
    cols_fill_isoheight = c(
        rep("white", round(100 * col_white_ratio)),
        colorRampPalette(brewer.pal(5, "YlOrRd")[2:5])(round(100 * (1 - col_white_ratio)))
    ),
    size_bg = 0.1,
    size_top = size_bg,
    nrow = 2,
    expand_margin = 0.05  # ✅ 新增参数：边缘扩展比例（5%）
) {

    density_top  <- enquo(density_top)

    df <- .data@meta.data %>%
        rownames_to_column("cellid") %>%
        inner_join(GetAllCoordinates(.data)) %>%
        as_tibble()

    # ✅ 计算坐标范围并扩展
    col_range <- range(df$col, na.rm = TRUE)
    row_range <- range(df$row, na.rm = TRUE)
    
    col_expand <- diff(col_range) * expand_margin
    row_expand <- diff(row_range) * expand_margin
    
    col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
    row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
    
    cat(sprintf("✅ 坐标范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
                col_limits[1], col_limits[2], row_limits[1], row_limits[2]))

    # ✅ 绘图
    p <- ggplot(mapping = aes(x = col, y = row)) +
        # 1. 背景点
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        ) +
        ggnewscale::new_scale_fill() +
        
        # 2. 密度热图 (关键：设置 bins 和 expand)
        stat_density_2d_filled(
            data = df %>% filter(!!density_top),
            mapping = aes(fill = after_stat(ndensity)),
            geom = "raster", 
            contour = FALSE,
            alpha = 0.8,
            bins = 100,  # ✅ 增加分辨率
            show.legend = TRUE,
            n = 200      # ✅ 增加网格密度
        ) +
        scale_fill_gradientn(
            colours = cols_fill_isoheight,
            name = "Density"
        ) +
        guides(alpha = "none") +
        ggnewscale::new_scale_fill() +
        
        # 3. 等高线
        geom_density_2d(
            data = df %>% filter(!!density_top),
            color = col_isoheight,
            contour_var = "ndensity",
            bins = 10,    # ✅ 等高线数量
            show.legend = FALSE
        ) +
        
        # 4. 高亮点
        geom_point(
            data = df %>% filter(!!density_top),
            color = col_top, alpha = 0.5, size = size_top
        ) +
        
        # ✅ 关键：手动设置坐标范围
        scale_x_continuous(
            limits = col_limits,
            expand = expansion(mult = 0.02)  # 额外 2% 边距
        ) +
        scale_y_reverse(
            limits = rev(row_limits),  # ✅ 注意反转顺序
            expand = expansion(mult = 0.02)
        ) +
        
        # 坐标系统
        coord_fixed(ratio = 1) +  # ✅ 保持宽高比
        
        # 主题
        NoAxes() +
        theme(
            aspect.ratio = 1,
            panel.background = element_blank(),
            plot.background = element_blank(),  # ✅ 移除绘图背景
            strip.background = element_blank(),
            legend.position = "right",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.margin = margin(5, 5, 5, 5)  # ✅ 统一边距
        ) +
        facet_wrap(vars(sample), nrow = nrow)
    
    return(p)
}
```

---

