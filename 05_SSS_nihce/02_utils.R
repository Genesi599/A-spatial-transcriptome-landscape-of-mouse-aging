

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


