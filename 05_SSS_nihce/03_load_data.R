#!/usr/bin/env Rscript
# ===================================================================
# 数据加载模块
# ===================================================================

load_gene_list <- function(config) {
  cat("📄 读取基因列表...\n")
  
  cache_file <- file.path(config$cache_dir, "gene_list.rds")
  
  if (is_cache_valid(cache_file, config$gene_list_path)) {
    gene_list <- load_cache(cache_file, "基因列表")
  } else {
    gene_list <- read.table(config$gene_list_path, header = TRUE, stringsAsFactors = FALSE)[[1]]
    gene_list <- trimws(gene_list[gene_list != ""])
    save_cache(gene_list, cache_file, "基因列表")
  }
  
  cat(sprintf("✅ 共读取 %d 个基因\n\n", length(gene_list)))
  return(gene_list)
}

load_seurat_object <- function(config) {
  cat("🧠 加载 Seurat 对象...\n")
  
  seurat_obj <- readRDS(config$seurat_path)
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  
  cat(sprintf("✅ Spots: %d, Genes: %d\n\n", ncol(seurat_obj), nrow(seurat_obj)))
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