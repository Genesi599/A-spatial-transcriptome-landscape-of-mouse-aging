# ==================================================
# 09_cache.R
# 缓存管理工具函数
# ==================================================

generate_plot_cache_key <- function(sample_id, CONFIG) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("需要安装 digest 包: install.packages('digest')")
  }
  
  key_params <- list(
    sample_id = sample_id,
    n_zones = CONFIG$params$n_zones,
    celltype_col = CONFIG$params$celltype_col,
    version = "v2.7"
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