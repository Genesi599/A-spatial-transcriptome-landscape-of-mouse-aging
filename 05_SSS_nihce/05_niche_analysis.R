#!/usr/bin/env Rscript
# ===================================================================
# Niche 距离计算
# ===================================================================

perform_niche_analysis <- function(seurat_obj, threshold, config) {
  cat("📈 开始 Niche 分析...\n")
  
  # 生成缓存key
  cache_key <- generate_cache_key(
    threshold, 
    sum(seurat_obj$ClockGene_High), 
    ncol(seurat_obj), 
    config$niche_dist_method
  )
  cache_file <- file.path(config$cache_dir, sprintf("niche_analysis_%s.rds", cache_key))
  
  if (file.exists(cache_file)) {
    niche_data <- load_cache(cache_file, "Niche 距离")
    seurat_obj$ClockGene_Distance <- niche_data$ClockGene_Distance
  } else {
    cat("🔄 正在进行 Niche 分析（多线程）...\n")
    plan(multisession, workers = config$n_workers)
    
    seurat_obj <- niche_marker(
      .data = seurat_obj,
      marker = ClockGene_High,
      spot_type = ClockGene_Distance,
      slide = orig.ident,
      dist_method = config$niche_dist_method,
      FUN = NA,
      n_work = config$n_workers
    )
    
    niche_data <- data.frame(ClockGene_Distance = seurat_obj$ClockGene_Distance)
    save_cache(niche_data, cache_file, "Niche 距离")
  }
  
  cat(sprintf("✅ 距离范围: %.2f ~ %.2f\n\n",
              min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
              max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))
  
  return(seurat_obj)
}