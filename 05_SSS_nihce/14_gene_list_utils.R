# ==================================================
# 14_gene_list_utils.R
# 基因列表扫描和加载工具
# ==================================================

scan_gene_lists <- function(path, pattern = "\\.txt$") {
  if (file.exists(path) && !dir.exists(path)) {
    # 单个文件
    return(normalizePath(path, winslash = "/"))
  }
  
  if (!dir.exists(path)) {
    stop(sprintf("❌ 基因列表路径不存在: %s", path))
  }
  
  # 目录扫描
  gene_files <- list.files(
    path,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(gene_files) == 0) {
    stop(sprintf(
      "❌ 未在 %s 中找到匹配 %s 的基因列表文件",
      path, pattern
    ))
  }
  
  gene_files <- normalizePath(gene_files, winslash = "/")
  
  return(gene_files)
}

load_gene_list_safe <- function(gene_file) {
  if (!file.exists(gene_file)) {
    stop(sprintf("❌ 基因列表文件不存在: %s", gene_file))
  }
  
  genes <- readLines(gene_file, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[genes != "" & !startsWith(genes, "#")]
  
  if (length(genes) == 0) {
    stop(sprintf("❌ 基因列表为空: %s", gene_file))
  }
  
  return(genes)
}

get_genelist_name <- function(gene_file) {
  basename <- tools::file_path_sans_ext(basename(gene_file))
  # 移除可能的前缀（如 "NET_gene_list_mouse" -> "mouse"）
  basename <- gsub("^(gene_list_|genelist_)", "", 
                   basename, ignore.case = TRUE)
  return(basename)
}

print_genelist_info <- function(gene_files) {
  cat("\n📋 基因列表文件:\n")
  for (i in seq_along(gene_files)) {
    genelist_name <- get_genelist_name(gene_files[i])
    n_genes <- length(load_gene_list_safe(gene_files[i]))
    cat(sprintf(
      "   [%2d] %s (%d genes)\n",
      i, genelist_name, n_genes
    ))
  }
  cat("\n")
}

cat("✅ 14_gene_list_utils.R 已加载\n")