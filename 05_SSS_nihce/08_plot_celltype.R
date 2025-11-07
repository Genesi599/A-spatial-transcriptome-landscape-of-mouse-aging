# ==================================================================
# 08_plot_celltype.R (不标准化细胞类型名称)
# 细胞类型在密度区域中的分布分析
# Author: Assistant | Date: 2025-11-07 | Version: 2.6
# ==================================================================

cat("🔧 加载 08_plot_celltype.R (无标准化版)...\n")

# ==================================================================
# 加载依赖工具函数
# ==================================================================

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
  source(file.path(utils_dir, "10_summary.R"))
}, error = function(e) {
  stop(sprintf("❌ 工具函数加载失败: %s", e$message))
})


# ==================================================================
# 缓存管理
# ==================================================================

generate_plot_cache_key <- function(sample_id, CONFIG) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("需要安装 digest 包: install.packages('digest')")
  }
  
  key_params <- list(
    sample_id = sample_id,
    n_zones = CONFIG$params$n_zones,
    celltype_col = CONFIG$params$celltype_col,
    version = "v2.6"
  )
  
  digest::digest(key_params, algo = "md5")
}

save_plot_cache <- function(sample_id, plot_data, CONFIG) {
  if (is.null(CONFIG$cache_dir) || !CONFIG$debug_mode) {
    return(invisible(NULL))
  }
  
  if (!dir.exists(CONFIG$cache_dir)) {
    dir.create(CONFIG$cache_dir, 
               recursive = TRUE, 
               showWarnings = FALSE)
  }
  
  cache_key <- generate_plot_cache_key(sample_id, CONFIG)
  cache_file <- file.path(
    CONFIG$cache_dir, 
    sprintf("celltype_plot_%s.rds", cache_key)
  )
  
  tryCatch({
    saveRDS(plot_data, cache_file)
    file_size_mb <- file.size(cache_file) / 1024^2
    cat(sprintf("      💾 缓存已保存: %.2f MB\n", file_size_mb))
  }, error = function(e) {
    warning(sprintf("      ⚠️  缓存保存失败: %s", e$message))
  })
  
  invisible(cache_file)
}

load_plot_cache <- function(sample_id, CONFIG) {
  if (is.null(CONFIG$cache_dir) || !CONFIG$debug_mode) {
    return(NULL)
  }
  
  cache_key <- generate_plot_cache_key(sample_id, CONFIG)
  cache_file <- file.path(
    CONFIG$cache_dir, 
    sprintf("celltype_plot_%s.rds", cache_key)
  )
  
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

clean_expired_cache <- function(CONFIG) {
  if (is.null(CONFIG$cache_dir) || 
      is.null(CONFIG$cache_max_age_hours)) {
    return(invisible(NULL))
  }
  
  cache_files <- list.files(
    CONFIG$cache_dir, 
    pattern = "^celltype_plot_.*\\.rds$", 
    full.names = TRUE
  )
  
  if (length(cache_files) == 0) return(invisible(NULL))
  
  current_time <- Sys.time()
  max_age_secs <- CONFIG$cache_max_age_hours * 3600
  
  expired_files <- cache_files[sapply(cache_files, function(f) {
    age_secs <- as.numeric(
      difftime(current_time, file.info(f)$mtime, units = "secs")
    )
    age_secs > max_age_secs
  })]
  
  if (length(expired_files) > 0) {
    cat(sprintf(
      "   🗑️  清理 %d 个过期缓存\n", 
      length(expired_files)
    ))
    unlink(expired_files)
  }
  
  invisible(expired_files)
}

list_cache_info <- function(CONFIG) {
  if (is.null(CONFIG$cache_dir) || 
      !dir.exists(CONFIG$cache_dir)) {
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
    size_mb = sapply(
      cache_files, 
      function(f) file.size(f) / 1024^2
    ),
    modified = sapply(
      cache_files, 
      function(f) as.character(file.info(f)$mtime)
    ),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("📂 缓存列表 (%s):\n", CONFIG$cache_dir))
  print(cache_info)
  cat(sprintf(
    "\n总计: %d 文件, %.2f MB\n", 
    nrow(cache_info), 
    sum(cache_info$size_mb)
  ))
  
  invisible(cache_info)
}

clear_all_cache <- function(CONFIG, confirm = TRUE) {
  if (is.null(CONFIG$cache_dir) || 
      !dir.exists(CONFIG$cache_dir)) {
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
  
  total_size_mb <- sum(sapply(cache_files, file.size)) / 1024^2
  
  if (confirm) {
    cat(sprintf(
      "⚠️  将删除 %d 文件 (%.2f MB)\n", 
      length(cache_files), 
      total_size_mb
    ))
    response <- readline(prompt = "确认删除? (yes/no): ")
    if (tolower(trimws(response)) != "yes") {
      cat("❌ 已取消\n")
      return(invisible(NULL))
    }
  }
  
  unlink(cache_files)
  cat(sprintf("✅ 已删除 %d 文件\n", length(cache_files)))
  invisible(cache_files)
}


# ==================================================================
# 单样本处理（保持原始细胞类型名称）
# ==================================================================

process_single_sample <- function(df, sample_id, CONFIG) {
  # ------------------------------
  # 小工具: 安全拿到 data.frame + 坐标列
  # ------------------------------
  pick_coord_cols <- function(md, coords = NULL) {
    # 常见候选
    cands <- list(
      list(col = "col",                row = "row"),
      list(col = "x",                  row = "y"),
      list(col = "imagecol",           row = "imagerow"),
      list(col = "pxl_col_in_fullres", row = "pxl_row_in_fullres"),
      list(col = "array_col",          row = "array_row")
    )
    if (!is.null(coords)) {
      if (all(c("x","y") %in% colnames(coords))) return(list(col = "x", row = "y", src = "coords"))
      if (all(c("imagecol","imagerow") %in% colnames(coords))) return(list(col = "imagecol", row = "imagerow", src = "coords"))
    }
    for (c in cands) {
      if (all(c(c$col, c$row) %in% colnames(md))) return(list(col = c$col, row = c$row, src = "meta"))
    }
    return(NULL)
  }

  get_df_std <- function(x) {
    if (inherits(x, "Seurat")) {
      seu <- x
      md <- seu@meta.data
      coords <- NULL
      if (requireNamespace("Seurat", quietly = TRUE)) {
        coords <- tryCatch(Seurat::GetTissueCoordinates(seu), error = function(e) NULL)
        if (!is.null(coords) && is.null(rownames(coords)) && "barcode" %in% colnames(coords)) {
          rownames(coords) <- coords$barcode
        }
        if (!is.null(coords)) {
          common <- intersect(rownames(md), rownames(coords))
          md <- md[common, , drop = FALSE]
          coords <- coords[common, , drop = FALSE]
        }
      }
      info <- pick_coord_cols(md, coords)
      if (is.null(info)) {
        stop("无法在 Seurat 对象中识别坐标列，请检查 meta.data 或 GetTissueCoordinates 结果")
      }
      if (identical(info$src, "coords")) {
        out <- cbind(md, coords[, c(info$col, info$row), drop = FALSE])
        colnames(out)[(ncol(out)-1):ncol(out)] <- c("col", "row")
      } else {
        out <- md
        colnames(out)[match(c(info$col, info$row), colnames(out))] <- c("col", "row")
      }
      out <- as.data.frame(out, stringsAsFactors = FALSE)
      return(out)
    } else if (is.data.frame(x)) {
      out <- x
      # 若无标准列名，尝试自动映射
      if (!all(c("col","row") %in% colnames(out))) {
        info <- pick_coord_cols(out, NULL)
        if (is.null(info)) {
          stop("数据中缺少坐标列（col/row 或常见别名 x/y, imagecol/imagerow, pxl_col_in_fullres/pxl_row_in_fullres）")
        }
        colnames(out)[match(c(info$col, info$row), colnames(out))] <- c("col", "row")
      }
      return(out)
    } else {
      stop("df 必须是 Seurat 或 data.frame")
    }
  }

  qout <- function(...) {
    # 受全局 CONFIG$quiet 控制
    if (!isTRUE(CONFIG$quiet)) cat(sprintf(...))
  }

  # ------------------------------
  # 标题输出（受 quiet 控制）
  # ------------------------------
  qout("\n[%s]\n", sample_id)

  # ------------------------------
  # 加载缓存
  # ------------------------------
  cached_data <- load_plot_cache(sample_id, CONFIG)
  if (!is.null(cached_data)) {
    qout("      🎨 使用缓存数据...\n")

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

    overlay_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_overlay.png", sample_id))
    composition_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_composition.png", sample_id))

    ggsave(overlay_file, p_overlay, width = 16, height = 12, dpi = 300, bg = "white")
    ggsave(composition_file, p_composition, width = 14, height = 10, dpi = 300, bg = "white")

    n_spots <- nrow(cached_data$df)
    n_high  <- sum(!is.na(cached_data$df$density_zone))
    n_types <- length(setdiff(unique(cached_data$df$celltype_clean), "Unknown"))

    qout("  ✅ %d spots | %d high | %d types (缓存)\n", n_spots, n_high, n_types)

    return(list(
      density_data = cached_data$density_data,
      zone_composition = cached_data$zone_composition,
      plots = list(overlay = p_overlay, composition = p_composition),
      stats = list(n_spots = n_spots, n_high_density = n_high, n_celltypes = n_types),
      from_cache = TRUE
    ))
  }

  # ------------------------------
  # 颜色方案验证
  # ------------------------------
  if (is.null(CONFIG$colors$celltype)) {
    stop("❌ 全局颜色方案未初始化！")
  }

  # ------------------------------
  # 标准化 df：确保是 data.frame 且有 col/row
  # ------------------------------
  df <- get_df_std(df)

  # ------------------------------
  # 清洗细胞类型（仅 NA/空串 -> Unknown）
  # ------------------------------
  if (is.null(CONFIG$params$celltype_col) || !(CONFIG$params$celltype_col %in% colnames(df))) {
    stop(sprintf("❌ 细胞类型列 '%s' 不存在", CONFIG$params$celltype_col))
  }
  raw_celltypes <- df[[CONFIG$params$celltype_col]]
  if (is.null(raw_celltypes) || length(raw_celltypes) == 0) {
    stop(sprintf("❌ 细胞类型列 '%s' 为空", CONFIG$params$celltype_col))
  }
  df$celltype_clean <- ifelse(is.na(raw_celltypes) | raw_celltypes == "", "Unknown", as.character(raw_celltypes))

  # 打印简洁的类型信息（静默模式不打印）
  unique_types <- sort(unique(df$celltype_clean[df$celltype_clean != "Unknown"]))
  n_types <- length(unique_types)
  qout("  📊 细胞类型: %d 个\n", n_types)
  if (n_types > 0 && !isTRUE(CONFIG$quiet)) {
    head_n <- min(10, n_types)
    for (i in seq_len(head_n)) qout("     • %s\n", unique_types[i])
    if (n_types > head_n) qout("     ... 还有 %d 个\n", n_types - head_n)
  }

  # 检查未知类型（不终止，仅 warn）
  all_types_global <- names(CONFIG$colors$celltype)
  sample_types <- setdiff(unique(df$celltype_clean), "Unknown")
  missing_types <- setdiff(sample_types, all_types_global)
  if (length(missing_types) > 0) {
    warning(sprintf("  ⚠️  未知类型: %s", paste(missing_types, collapse = ", ")))
  }

  # ------------------------------
  # 计算密度区域（静默）
  # ------------------------------
  expand_margin <- if (!is.null(CONFIG$params$expand_margin)) CONFIG$params$expand_margin else 0.1
  # 确保有 ClockGene_High 列（calculate_density_zones 需要）
  if (!("ClockGene_High" %in% colnames(df))) {
    stop("缺少列 'ClockGene_High'，无法进行密度计算")
  }
  density_data <- calculate_density_zones(
    df = df,
    density_bins  = CONFIG$params$n_zones,
    expand_margin = expand_margin,
    quiet = TRUE
  )
  if (is.null(density_data)) {
    warning(sprintf("[%s] 密度计算失败或高表达点不足，跳过该样本", sample_id))
    return(NULL)
  }

  # 将密度 zone 合并回 df
  # density_data$spot_zones: col,row,density_zone,density_value
  df <- df %>%
    dplyr::left_join(
      density_data$spot_zones,
      by = c("col", "row")
    )

  # ------------------------------
  # 计算每个 zone 的细胞组成
  # ------------------------------
  zone_composition <- df %>%
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(density_zone) %>%
    dplyr::mutate(
      total = sum(count),
      percentage = (count / total) * 100
    ) %>%
    dplyr::ungroup()

  # ------------------------------
  # 保存缓存（仅 debug 模式）
  # ------------------------------
  if (isTRUE(CONFIG$debug_mode)) {
    plot_data <- list(
      df = df,
      density_data = density_data,
      zone_composition = zone_composition,
      params = list(
        sample_id = sample_id,
        n_spots = nrow(df),
        n_zones = CONFIG$params$n_zones,
        cache_time = Sys.time()
      )
    )
    save_plot_cache(sample_id, plot_data, CONFIG)
  }

  # ------------------------------
  # 绘图
  # ------------------------------
  p_overlay <- plot_celltype_density_overlay(df, density_data, sample_id, CONFIG)
  p_composition <- plot_zone_composition(zone_composition, sample_id, CONFIG)

  # ------------------------------
  # 保存图与数据
  # ------------------------------
  overlay_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_overlay.png", sample_id))
  composition_file <- file.path(CONFIG$output$plot_dir, sprintf("%s_composition.png", sample_id))
  ggsave(overlay_file, p_overlay, width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(composition_file, p_composition, width = 14, height = 10, dpi = 300, bg = "white")

  zone_comp_file <- file.path(CONFIG$output$data_dir, sprintf("%s_zone_composition.csv", sample_id))
  utils::write.csv(zone_composition, zone_comp_file, row.names = FALSE)

  # ------------------------------
  # 统计输出（静默）
  # ------------------------------
  n_spots <- nrow(df)
  n_high  <- sum(!is.na(df$density_zone))
  n_types <- length(setdiff(unique(df$celltype_clean), "Unknown"))
  qout("  ✅ %d spots | %d high | %d types\n", n_spots, n_high, n_types)

  return(list(
    density_data = density_data,
    zone_composition = zone_composition,
    plots = list(overlay = p_overlay, composition = p_composition),
    stats = list(n_spots = n_spots, n_high_density = n_high, n_celltypes = n_types),
    from_cache = FALSE
  ))
}


# ==================================================================
# 创建全局颜色方案（保持原始名称）
# ==================================================================

create_global_color_scheme <- function(
    data_list, 
    celltype_col, 
    n_zones = 10) {
  
  cat("\n🎨 生成全局颜色方案...\n")
  
  all_celltypes <- unique(unlist(lapply(data_list, function(df) {
    ct <- df[[celltype_col]]
    ct <- as.character(ct)
    unique(ct[!is.na(ct) & ct != "" & ct != "Unknown"])
  })))
  
  all_celltypes <- sort(all_celltypes)
  n_celltypes <- length(all_celltypes)
  
  cat(sprintf("  📊 发现 %d 个细胞类型\n", n_celltypes))
  
  if (n_celltypes <= 10) {
    for (ct in all_celltypes) {
      cat(sprintf("     • %s\n", ct))
    }
  } else {
    for (i in 1:10) {
      cat(sprintf("     • %s\n", all_celltypes[i]))
    }
    cat(sprintf("     ... 还有 %d 个\n", n_celltypes - 10))
  }
  
  # 生成颜色
  celltype_colors <- get_celltype_colors(all_celltypes)
  zone_colors <- get_zone_colors(n_zones)
  names(zone_colors) <- sprintf("Zone_%d", 0:(n_zones - 1))
  
  cat(sprintf(
    "  ✅ 完成 (%d 类型 + %d 区域)\n", 
    n_celltypes, n_zones
  ))
  
  list(celltype = celltype_colors, density_zone = zone_colors)
}


# ==================================================================
# 综合分析
# ==================================================================

run_celltype_analysis <- function(data_list, sample_ids, CONFIG) {
  
  cat(paste0(
    "\n╔═══════════════════════════════════════════════════",
    "═════════╗\n"
  ))
  cat("║  细胞类型分布分析                                ")
  cat("        ║\n")
  cat(paste0(
    "╚═══════════════════════════════════════════════════",
    "═════════╝\n"
  ))
  
  # 清理过期缓存
  if (CONFIG$debug_mode && 
      !is.null(CONFIG$cache_max_age_hours)) {
    clean_expired_cache(CONFIG)
  }
  
  if (CONFIG$debug_mode && !is.null(CONFIG$cache_dir)) {
    cat(sprintf(
      "\n🔧 调试模式开启\n📂 缓存: %s\n", 
      CONFIG$cache_dir
    ))
  }
  
  # 创建颜色方案
  CONFIG$colors <- create_global_color_scheme(
    data_list, 
    CONFIG$params$celltype_col, 
    CONFIG$params$n_zones
  )
  
  # 处理样本
  cat("\n🔬 处理样本...\n")
  
  results_list <- list()
  n_from_cache <- 0
  
  for (i in seq_along(data_list)) {
    cat(sprintf("\n[%2d/%2d]", i, length(data_list)))
    
    result <- process_single_sample(
      data_list[[i]], 
      sample_ids[i], 
      CONFIG
    )
    
    if (!is.null(result$from_cache) && result$from_cache) {
      n_from_cache <- n_from_cache + 1
    }
    
    results_list[[sample_ids[i]]] <- result
  }
  
  # 缓存统计
  if (CONFIG$debug_mode && n_from_cache > 0) {
    cat(sprintf(
      "\n💾 缓存命中: %d/%d (%.1f%%)\n", 
      n_from_cache, 
      length(data_list), 
      100 * n_from_cache / length(data_list)
    ))
  }
  
  # 合并数据
  cat("\n\n📊 合并数据...\n")
  
  combined_data <- do.call(rbind, lapply(
    names(results_list), 
    function(sid) {
      comp <- results_list[[sid]]$zone_composition
      comp$sample <- sid
      comp
    }
  ))
  
  combined_file <- file.path(
    CONFIG$output$data_dir, 
    "combined_zone_composition.csv"
  )
  write.csv(combined_data, combined_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(combined_file)))
  
  # 热图
  cat("\n📊 生成热图...\n")
  p_heatmap <- plot_combined_heatmap(combined_data, CONFIG)
  heatmap_file <- file.path(
    CONFIG$output$plot_dir, 
    "combined_heatmap.png"
  )
  ggsave(heatmap_file, p_heatmap, 
         width = 18, height = 14, dpi = 300, bg = "white")
  cat(sprintf("  ✅ %s\n", basename(heatmap_file)))
  
  # 综合图
  cat("\n📊 生成综合图...\n")
  p_combined <- plot_combined_analysis(combined_data, CONFIG)
  combined_plot_file <- file.path(
    CONFIG$output$plot_dir, 
    "combined_analysis.png"
  )
  ggsave(combined_plot_file, p_combined, 
         width = 20, height = 16, dpi = 300, bg = "white")
  cat(sprintf("  ✅ %s\n", basename(combined_plot_file)))
  
  # 统计摘要
  cat("\n📊 生成统计...\n")
  summary_stats <- generate_summary_statistics(
    combined_data, CONFIG
  )
  summary_file <- file.path(
    CONFIG$output$data_dir, 
    "summary_statistics.csv"
  )
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat(sprintf("  ✅ %s\n", basename(summary_file)))
  
  cat(paste0(
    "\n╔═══════════════════════════════════════════════════",
    "═════════╗\n"
  ))
  cat("║  ✅ 分析完成！                                    ")
  cat("        ║\n")
  cat(paste0(
    "╚═══════════════════════════════════════════════════",
    "═════════╝\n\n"
  ))
  
  list(
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
  )
}


# ==================================================================
# 主接口函数
# ==================================================================

analyze_celltype_niche <- function(
    sample_list, 
    CONFIG, 
    seurat_basename = NULL) {
  
  sample_ids <- names(sample_list)
  
  # 自适应配置输出目录
  output_needs_init <- is.null(CONFIG$output) || 
    is.null(CONFIG$output$plot_dir) || 
    is.null(CONFIG$output$data_dir)
  
  if (output_needs_init) {
    
    cat("   🔧 初始化输出目录...\n")
    
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
      base_metadata_dir <- file.path(
        CONFIG$output_dir, "metadata"
      )
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
  }
  
  # 创建目录
  for (dir_path in CONFIG$output) {
    if (!is.null(dir_path) && !dir.exists(dir_path)) {
      dir.create(dir_path, 
                 recursive = TRUE, 
                 showWarnings = FALSE)
    }
  }
  
  # 运行分析
  run_celltype_analysis(
    data_list = sample_list, 
    sample_ids = sample_ids, 
    CONFIG = CONFIG
  )
}


# ==================================================================
# 导出
# ==================================================================

cat("✅ 08_plot_celltype.R 加载完成 (v2.6 - 无标准化)\n")
cat("📚 主要函数:\n")
cat(paste0(
  "  - analyze_celltype_niche(",
  "sample_list, CONFIG, seurat_basename)\n"
))
cat(paste0(
  "  - run_celltype_analysis(",
  "data_list, sample_ids, CONFIG)\n"
))
cat(paste0(
  "  - create_global_color_scheme(",
  "data_list, celltype_col, n_zones)\n"
))
cat("  - list_cache_info(CONFIG)\n")
cat("  - clear_all_cache(CONFIG)\n\n")