
# ===================================================================
# 08_plot_celltype.R (带调试缓存 - 完整版)
# 细胞类型在密度区域中的分布分析（全局统一配色版 + 调试缓存）
# Author: Assistant | Date: 2025-11-07 | Version: 2.3
# ===================================================================

cat("🔧 加载 08_plot_celltype.R (带调试缓存)...\n")

# ===================================================================
# 自动检测script_dir
# ===================================================================

if (!exists("script_dir")) {
  # 从当前脚本路径推断
  current_script <- tryCatch({
    # 方法1: 使用sys.frame
    normalizePath(sys.frame(1)$ofile, winslash = "/")
  }, error = function(e) {
    # 方法2: 使用commandArgs (适用于Rscript)
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      sub("^--file=", "", file_arg)
    } else {
      # 方法3: 使用当前工作目录
      file.path(getwd(), "08_plot_celltype.R")
    }
  })
  
  script_dir <- dirname(current_script)
  cat(sprintf("   📂 脚本目录: %s\n", script_dir))
}

# 检查工具函数目录是否存在
utils_dir <- file.path(script_dir, "08_plot_celltype_utils")

if (!dir.exists(utils_dir)) {
  stop(sprintf("❌ 工具函数目录不存在: %s\n请确保 08_plot_celltype_utils 文件夹在正确位置", utils_dir))
}

# 加载工具函数
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
  source(file.path(utils_dir, "08_validation.R"))
  source(file.path(utils_dir, "10_summary.R"))
  
  cat("   ✅ 所有工具已加载\n\n")
  validate_required_functions()
  
}, error = function(e) {
  stop(sprintf("❌ 工具函数加载失败: %s", e$message))
})


# ===================================================================
# 缓存管理函数
# ===================================================================

#' 生成绘图数据的缓存key
#' 
#' @param sample_id 样本ID
#' @param CONFIG 配置对象
#' @return 缓存key字符串
#'
generate_plot_cache_key <- function(sample_id, CONFIG) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("需要安装 digest 包: install.packages('digest')")
  }
  
  # 提取影响绘图的关键参数
  key_params <- list(
    sample_id = sample_id,
    density_threshold = CONFIG$params$density_threshold_percentile,
    n_zones = CONFIG$params$n_zones,
    grid_resolution = CONFIG$params$grid_resolution,
    celltype_col = CONFIG$params$celltype_col,
    version = "v2.3"  # 版本标识，可在逻辑变化时强制更新缓存
  )
  
  cache_key <- digest::digest(key_params, algo = "md5")
  return(cache_key)
}


#' 保存绘图数据到缓存
#' 
#' @param sample_id 样本ID
#' @param plot_data 绘图数据列表
#' @param CONFIG 配置对象
#'
save_plot_cache <- function(sample_id, plot_data, CONFIG) {
  
  if (is.null(CONFIG$cache_dir) || !CONFIG$debug_mode) {
    return(invisible(NULL))
  }
  
  # 确保缓存目录存在
  if (!dir.exists(CONFIG$cache_dir)) {
    dir.create(CONFIG$cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 生成缓存文件名
  cache_key <- generate_plot_cache_key(sample_id, CONFIG)
  cache_file <- file.path(CONFIG$cache_dir, sprintf("celltype_plot_%s.rds", cache_key))
  
  tryCatch({
    saveRDS(plot_data, cache_file)
    file_size_mb <- file.size(cache_file) / 1024^2
    cat(sprintf("      💾 缓存已保存: %.2f MB\n", file_size_mb))
  }, error = function(e) {
    warning(sprintf("      ⚠️  缓存保存失败: %s", e$message))
  })
  
  invisible(cache_file)
}


#' 从缓存加载绘图数据
#' 
#' @param sample_id 样本ID
#' @param CONFIG 配置对象
#' @return 绘图数据列表或NULL
#'
load_plot_cache <- function(sample_id, CONFIG) {
  
  if (is.null(CONFIG$cache_dir) || !CONFIG$debug_mode) {
    return(NULL)
  }
  
  cache_key <- generate_plot_cache_key(sample_id, CONFIG)
  cache_file <- file.path(CONFIG$cache_dir, sprintf("celltype_plot_%s.rds", cache_key))
  
  if (!file.exists(cache_file)) {
    return(NULL)
  }
  
  tryCatch({
    plot_data <- readRDS(cache_file)
    file_size_mb <- file.size(cache_file) / 1024^2
    cat(sprintf("      📂 从缓存加载: %.2f MB\n", file_size_mb))
    return(plot_data)
  }, error = function(e) {
    warning(sprintf("      ⚠️  缓存加载失败: %s", e$message))
    return(NULL)
  })
}


#' 清理过期缓存
#' 
#' @param CONFIG 配置对象
#'
clean_expired_cache <- function(CONFIG) {
  
  if (is.null(CONFIG$cache_dir) || is.null(CONFIG$cache_max_age_hours)) {
    return(invisible(NULL))
  }
  
  cache_files <- list.files(
    CONFIG$cache_dir, 
    pattern = "^celltype_plot_.*\\.rds$", 
    full.names = TRUE
  )
  
  if (length(cache_files) == 0) {
    return(invisible(NULL))
  }
  
  current_time <- Sys.time()
  max_age_secs <- CONFIG$cache_max_age_hours * 3600
  
  expired_files <- c()
  
  for (cache_file in cache_files) {
    file_time <- file.info(cache_file)$mtime
    age_secs <- as.numeric(difftime(current_time, file_time, units = "secs"))
    
    if (age_secs > max_age_secs) {
      expired_files <- c(expired_files, cache_file)
    }
  }
  
  if (length(expired_files) > 0) {
    cat(sprintf("   🗑️  清理 %d 个过期缓存文件\n", length(expired_files)))
    unlink(expired_files)
  }
  
  invisible(expired_files)
}


#' 列出所有缓存文件信息
#' 
#' @param CONFIG 配置对象
#' @return 缓存文件信息数据框
#'
list_cache_info <- function(CONFIG) {
  
  if (is.null(CONFIG$cache_dir) || !dir.exists(CONFIG$cache_dir)) {
    cat("❌ 缓存目录不存在\n")
    return(invisible(NULL))
  }
  
  cache_files <- list.files(
    CONFIG$cache_dir, 
    pattern = "^celltype_plot_.*\\.rds$", 
    full.names = TRUE
  )
  
  if (length(cache_files) == 0) {
    cat("📂 缓存目录为空\n")
    return(invisible(NULL))
  }
  
  cache_info <- data.frame(
    file = basename(cache_files),
    size_mb = sapply(cache_files, function(f) file.size(f) / 1024^2),
    modified = sapply(cache_files, function(f) as.character(file.info(f)$mtime)),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("📂 缓存文件列表 (%s):\n", CONFIG$cache_dir))
  print(cache_info)
  cat(sprintf("\n总计: %d 个文件, %.2f MB\n", 
              nrow(cache_info), sum(cache_info$size_mb)))
  
  return(invisible(cache_info))
}


#' 清除所有celltype绘图缓存
#' 
#' @param CONFIG 配置对象
#' @param confirm 是否需要确认
#'
clear_all_cache <- function(CONFIG, confirm = TRUE) {
  
  if (is.null(CONFIG$cache_dir) || !dir.exists(CONFIG$cache_dir)) {
    cat("❌ 缓存目录不存在\n")
    return(invisible(NULL))
  }
  
  cache_files <- list.files(
    CONFIG$cache_dir, 
    pattern = "^celltype_plot_.*\\.rds$", 
    full.names = TRUE
  )
  
  if (length(cache_files) == 0) {
    cat("📂 缓存目录为空，无需清理\n")
    return(invisible(NULL))
  }
  
  total_size_mb <- sum(sapply(cache_files, file.size)) / 1024^2
  
  if (confirm) {
    cat(sprintf("⚠️  将删除 %d 个缓存文件 (%.2f MB)\n", length(cache_files), total_size_mb))
    response <- readline(prompt = "确认删除? (yes/no): ")
    
    if (tolower(trimws(response)) != "yes") {
      cat("❌ 已取消\n")
      return(invisible(NULL))
    }
  }
  
  unlink(cache_files)
  cat(sprintf("✅ 已删除 %d 个缓存文件\n", length(cache_files)))
  
  invisible(cache_files)
}


# ===================================================================
# 细胞类型名称标准化
# ===================================================================

#' 标准化细胞类型名称
#' 
#' @param names 细胞类型名称向量
#' @param mode 标准化模式 ("underscore"/"hyphen"/"space")
#' @param title_case 是否首字母大写
#' @return 标准化后的名称向量
#'
standardize_celltype_names <- function(names, mode = "underscore", title_case = TRUE) {
  
  # 基础清理
  names <- as.character(names)
  names <- trimws(names)
  names[is.na(names) | names == ""] <- "Unknown"
  
  # 统一分隔符
  if (mode == "underscore") {
    names <- gsub("-", "_", names)
    names <- gsub("\\s+", "_", names)
  } else if (mode == "hyphen") {
    names <- gsub("_", "-", names)
    names <- gsub("\\s+", "-", names)
  } else if (mode == "space") {
    names <- gsub("_", " ", names)
    names <- gsub("-", " ", names)
    names <- gsub("\\s+", " ", names)
  }
  
  # 首字母大写
  if (title_case) {
    separator <- if(mode == "underscore") "_" else if(mode == "hyphen") "-" else " "
    
    names <- sapply(names, function(name) {
      if (name == "Unknown") return("Unknown")
      parts <- strsplit(name, sprintf("[%s]", separator))[[1]]
      parts <- tolower(parts)
      parts <- paste0(toupper(substring(parts, 1, 1)), substring(parts, 2))
      paste(parts, collapse = separator)
    }, USE.NAMES = FALSE)
  }
  
  # 去除前后多余分隔符
  names <- gsub("^[_\\-\\s]+|[_\\-\\s]+$", "", names)
  
  # 保留常见缩写大写
  names <- gsub("\\b(Smc)\\b", "SMC", names)
  names <- gsub("\\b(Pp)\\b", "PP", names)
  names <- gsub("\\b(Bv)\\b", "BV", names)
  names <- gsub("\\b(Mv)\\b", "MV", names)
  names <- gsub("\\b(Tv)\\b", "TV", names)
  
  return(names)
}


# ===================================================================
# 单样本处理 (带缓存支持)
# ===================================================================

#' 处理单个样本 (支持调试缓存)
#'
#' @param df 数据框（包含坐标和细胞类型）
#' @param sample_id 样本ID
#' @param CONFIG 配置列表
#' @return 处理结果列表
#'
process_single_sample <- function(df, sample_id, CONFIG) {
  
  cat(sprintf("\n[%s]\n", sample_id))
  
  # ===================================================================
  # 1. 尝试加载缓存
  # ===================================================================
  
  cached_data <- load_plot_cache(sample_id, CONFIG)
  
  if (!is.null(cached_data)) {
    cat("      🎨 使用缓存数据直接绘图...\n")
    
    # 使用缓存数据绘图
    p_overlay <- plot_celltype_density_overlay(
      cached_data$df, 
      cached_data$density_data, 
      sample_id, 
      CONFIG
    )
    
    p_composition <- plot_zone_composition(
      cached_data$zone_composition, 
      sample_id, 
      CONFIG
    )
    
    # 保存图形
    overlay_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_overlay.png", sample_id))
    composition_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_composition.png", sample_id))
    
    ggsave(overlay_file, plot = p_overlay, width = 16, height = 12, dpi = 300, bg = "white")
    ggsave(composition_file, plot = p_composition, width = 14, height = 10, dpi = 300, bg = "white")
    
    # 统计信息
    n_spots <- nrow(cached_data$df)
    n_high_density <- sum(!is.na(cached_data$df$density_zone))
    n_celltypes <- length(setdiff(unique(cached_data$df$celltype_clean), "Unknown"))
    
    cat(sprintf("  ✅ %d spots | %d high | %d celltypes (缓存)\n", 
                n_spots, n_high_density, n_celltypes))
    
    return(list(
      density_data = cached_data$density_data,
      zone_composition = cached_data$zone_composition,
      plots = list(overlay = p_overlay, composition = p_composition),
      stats = list(
        n_spots = n_spots,
        n_high_density = n_high_density,
        n_celltypes = n_celltypes
      ),
      from_cache = TRUE
    ))
  }
  
  # ===================================================================
  # 2. 缓存不存在,正常计算
  # ===================================================================
  
  # 验证颜色方案
  if (is.null(CONFIG$colors$celltype)) {
    stop("❌ 全局颜色方案未初始化！请先调用 create_global_color_scheme()")
  }
  
  # 标准化细胞类型名称
  raw_celltypes <- df[[CONFIG$params$celltype_col]]
  df$celltype_clean <- standardize_celltype_names(raw_celltypes, mode = "underscore", title_case = TRUE)
  
  # 打印标准化示例
  unique_raw <- unique(raw_celltypes)
  unique_clean <- unique(df$celltype_clean)
  n_show <- min(5, length(unique_raw))
  
  cat("  🔄 细胞类型标准化:\n")
  for (i in 1:n_show) {
    if (unique_raw[i] != unique_clean[i]) {
      cat(sprintf("     '%s' → '%s'\n", unique_raw[i], unique_clean[i]))
    }
  }
  if (length(unique_raw) > 5) {
    cat(sprintf("     ... 还有 %d 个\n", length(unique_raw) - 5))
  }
  
  # 检查未知细胞类型
  all_celltypes_global <- names(CONFIG$colors$celltype)
  sample_celltypes <- setdiff(unique(df$celltype_clean), "Unknown")
  missing_types <- setdiff(sample_celltypes, all_celltypes_global)
  
  if (length(missing_types) > 0) {
    warning(sprintf("  ⚠️  未知细胞类型: %s", paste(missing_types, collapse = ", ")))
  }
  
  # 计算密度区域
  density_data <- calculate_density_zones(
    df = df,
    col_col = CONFIG$params$col_col,
    row_col = CONFIG$params$row_col,
    density_threshold_percentile = CONFIG$params$density_threshold_percentile,
    n_zones = CONFIG$params$n_zones,
    grid_resolution = CONFIG$params$grid_resolution
  )
  
  df$density_zone <- density_data$cell_zones
  
  # 计算zone组成
  zone_composition <- df %>%
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(density_zone) %>%
    dplyr::mutate(
      total = sum(count),
      percentage = (count / total) * 100
    ) %>%
    dplyr::ungroup()
  
  # ===================================================================
  # 3. 保存到缓存 (仅调试模式)
  # ===================================================================
  
  if (CONFIG$debug_mode) {
    plot_data <- list(
      df = df,
      density_data = density_data,
      zone_composition = zone_composition,
      # 保存关键参数便于验证
      params = list(
        sample_id = sample_id,
        n_spots = nrow(df),
        n_zones = CONFIG$params$n_zones,
        celltype_col = CONFIG$params$celltype_col,
        cache_time = Sys.time()
      )
    )
    
    save_plot_cache(sample_id, plot_data, CONFIG)
  }
  
  # ===================================================================
  # 4. 绘制图形
  # ===================================================================
  
  p_overlay <- plot_celltype_density_overlay(df, density_data, sample_id, CONFIG)
  p_composition <- plot_zone_composition(zone_composition, sample_id, CONFIG)
  
  # 保存图形
  overlay_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_overlay.png", sample_id))
  composition_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_composition.png", sample_id))
  
  ggsave(overlay_file, plot = p_overlay, width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(composition_file, plot = p_composition, width = 14, height = 10, dpi = 300, bg = "white")
  
  # 保存数据
  zone_comp_file <- file.path(CONFIG$output$data_dir, sprintf("%s_zone_composition.csv", sample_id))
  write.csv(zone_composition, zone_comp_file, row.names = FALSE)
  
  # 统计信息
  n_spots <- nrow(df)
  n_high_density <- sum(!is.na(df$density_zone))
  n_celltypes <- length(setdiff(unique(df$celltype_clean), "Unknown"))
  
  cat(sprintf("  ✅ %d spots | %d high | %d celltypes\n", n_spots, n_high_density, n_celltypes))
  
  return(list(
    density_data = density_data,
    zone_composition = zone_composition,
    plots = list(overlay = p_overlay, composition = p_composition),
    stats = list(
      n_spots = n_spots,
      n_high_density = n_high_density,
      n_celltypes = n_celltypes
    ),
    from_cache = FALSE
  ))
}


# ===================================================================
# 创建全局颜色方案
# ===================================================================

#' 创建全局统一颜色方案
#'
#' @param data_list 数据框列表
#' @param celltype_col 细胞类型列名
#' @param n_zones 密度区域数量
#' @return 颜色方案列表
#'
create_global_color_scheme <- function(data_list, celltype_col, n_zones = 10) {
  
  cat("\n🎨 生成全局颜色方案...\n")
  
  # 收集所有细胞类型
  all_celltypes <- lapply(data_list, function(df) {
    ct <- df[[celltype_col]]
    ct <- standardize_celltype_names(ct, mode = "underscore", title_case = TRUE)
    ct <- ct[ct != "Unknown"]
    return(unique(ct))
  })
  
  all_celltypes_unique <- sort(unique(unlist(all_celltypes)))
  n_celltypes <- length(all_celltypes_unique)
  
  cat(sprintf("  📊 发现 %d 个细胞类型（标准化后）\n", n_celltypes))
  
  # 打印列表
  if (n_celltypes <= 10) {
    for (ct in all_celltypes_unique) {
      cat(sprintf("     • %s\n", ct))
    }
  } else {
    for (i in 1:10) {
      cat(sprintf("     • %s\n", all_celltypes_unique[i]))
    }
    cat(sprintf("     ... 还有 %d 个\n", n_celltypes - 10))
  }
  
  # 生成颜色
  celltype_colors <- get_celltype_colors(all_celltypes_unique)
  zone_colors <- get_zone_colors(n_zones)
  names(zone_colors) <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf("  ✅ 颜色方案完成 (%d 细胞类型 + %d 区域)\n", n_celltypes, n_zones))
  
  return(list(
    celltype = celltype_colors,
    density_zone = zone_colors
  ))
}


# ===================================================================
# 综合分析 (带缓存统计)
# ===================================================================

#' 运行细胞类型分布综合分析
#'
#' @param data_list 数据框列表
#' @param sample_ids 样本ID向量
#' @param CONFIG 配置列表
#' @return 综合分析结果
#'
run_celltype_analysis <- function(data_list, sample_ids, CONFIG) {
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║  细胞类型在密度区域中的分布分析                          ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n")
  
  # 清理过期缓存
  if (CONFIG$debug_mode && !is.null(CONFIG$cache_max_age_hours)) {
    clean_expired_cache(CONFIG)
  }
  
  # 显示缓存模式
  if (CONFIG$debug_mode && !is.null(CONFIG$cache_dir)) {
    cat(sprintf("\n🔧 调试模式: 开启\n"))
    cat(sprintf("📂 缓存目录: %s\n", CONFIG$cache_dir))
  }
  
  # 创建全局颜色方案
  CONFIG$colors <- create_global_color_scheme(
    data_list = data_list,
    celltype_col = CONFIG$params$celltype_col,
    n_zones = CONFIG$params$n_zones
  )
  
  # 处理每个样本
  cat("\n🔬 处理样本...\n")
  
  results_list <- list()
  n_from_cache <- 0
  
  for (i in seq_along(data_list)) {
    cat(sprintf("\n[%2d/%2d]", i, length(data_list)))
    
    result <- process_single_sample(
      df = data_list[[i]],
      sample_id = sample_ids[i],
      CONFIG = CONFIG
    )
    
    if (!is.null(result$from_cache) && result$from_cache) {
      n_from_cache <- n_from_cache + 1
    }
    
    results_list[[sample_ids[i]]] <- result
  }
  
  # 打印缓存统计
  if (CONFIG$debug_mode && n_from_cache > 0) {
    cat(sprintf("\n💾 缓存命中: %d/%d (%.1f%%)\n", 
                n_from_cache, length(data_list), 
                100 * n_from_cache / length(data_list)))
  }
  
  # 合并数据
  cat("\n\n📊 合并数据...\n")
  
  combined_data <- do.call(rbind, lapply(names(results_list), function(sid) {
    comp <- results_list[[sid]]$zone_composition
    comp$sample <- sid
    return(comp)
  }))
  
  combined_file <- file.path(CONFIG$output$data_dir, "combined_zone_composition.csv")
  write.csv(combined_data, combined_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(combined_file)))
  
  # 绘制热图
  cat("\n📊 生成热图...\n")
  
  p_heatmap <- plot_combined_heatmap(combined_data, CONFIG)
  
  heatmap_file <- file.path(CONFIG$output$plot_dir, "combined_heatmap.png")
  ggsave(heatmap_file, plot = p_heatmap, width = 18, height = 14, dpi = 300, bg = "white")
  cat(sprintf("  ✅ %s\n", basename(heatmap_file)))
  
  # 绘制综合分析图
  cat("\n📊 生成综合分析图...\n")
  
  p_combined <- plot_combined_analysis(combined_data, CONFIG)
  
  combined_plot_file <- file.path(CONFIG$output$plot_dir, "combined_analysis.png")
  ggsave(combined_plot_file, plot = p_combined, width = 20, height = 16, dpi = 300, bg = "white")
  cat(sprintf("  ✅ %s\n", basename(combined_plot_file)))
  
  # 生成统计摘要
  cat("\n📊 生成统计摘要...\n")
  
  summary_stats <- generate_summary_statistics(combined_data, CONFIG)
  
  summary_file <- file.path(CONFIG$output$data_dir, "summary_statistics.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(summary_file)))
  
  # 完成
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║  ✅ 分析完成！                                            ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  return(list(
    individual_results = results_list,
    combined_data = combined_data,
    summary_statistics = summary_stats,
    combined_plots = list(
      heatmap = p_heatmap,
      analysis = p_combined
    ),
    config = CONFIG,
    cache_stats = list(
      n_from_cache = n_from_cache,
      total = length(data_list),
      cache_hit_rate = 100 * n_from_cache / length(data_list)
    )
  ))
}


# ===================================================================
# 主接口函数（兼容原有调用）
# ===================================================================

#' 细胞类型Niche分析主函数
#'
#' @param sample_list 样本列表
#' @param CONFIG 配置对象
#' @param seurat_basename Seurat对象基础名称
#' @return 分析结果
#'
analyze_celltype_niche <- function(sample_list, CONFIG, seurat_basename = NULL) {
  
  # 提取样本ID
  sample_ids <- names(sample_list)
  
  # 确保输出目录存在
  if (!dir.exists(CONFIG$output$plot_dir)) {
    dir.create(CONFIG$output$plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(CONFIG$output$data_dir)) {
    dir.create(CONFIG$output$data_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 运行分析
  results <- run_celltype_analysis(
    data_list = sample_list,
    sample_ids = sample_ids,
    CONFIG = CONFIG
  )
  
  return(results)
}


# ===================================================================
# 导出可用函数列表
# ===================================================================

cat("✅ 08_plot_celltype.R 已加载 (支持调试缓存)\n")
cat("📚 可用函数:\n")
cat("  主函数:\n")
cat("    - analyze_celltype_niche(sample_list, CONFIG, seurat_basename)\n")
cat("    - run_celltype_analysis(data_list, sample_ids, CONFIG)\n")
cat("  缓存管理:\n")
cat("    - list_cache_info(CONFIG)           # 列出缓存信息\n")
cat("    - clear_all_cache(CONFIG)           # 清除所有缓存\n")
cat("    - clean_expired_cache(CONFIG)       # 清理过期缓存\n")
cat("  辅助函数:\n")
cat("    - create_global_color_scheme(data_list, celltype_col, n_zones)\n")
cat("    - standardize_celltype_names(names, mode, title_case)\n")
cat("\n")
