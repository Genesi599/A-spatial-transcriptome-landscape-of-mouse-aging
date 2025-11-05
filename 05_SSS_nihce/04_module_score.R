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