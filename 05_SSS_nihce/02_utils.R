

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

is_cache_valid <- function(cache_file, source_file = NULL, max_age_hours = NULL) {
  if (!file.exists(cache_file)) return(FALSE)
  
  if (!is.null(source_file) && file.exists(source_file)) {
    if (file.mtime(cache_file) < file.mtime(source_file)) {
      cat("⚠️ 源文件已更新，缓存失效\n")
      return(FALSE)
    }
  }
  
  if (!is.null(max_age_hours)) {
    age <- difftime(Sys.time(), file.mtime(cache_file), units = "hours")
    if (age > max_age_hours) {
      cat(sprintf("⚠️ 缓存已过期 (%.1f 小时)\n", age))
      return(FALSE)
    }
  }
  
  TRUE
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