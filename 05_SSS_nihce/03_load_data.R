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
  
  # ✅ 确保目录存在
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