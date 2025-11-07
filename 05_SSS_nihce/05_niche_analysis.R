#!/usr/bin/env Rscript
# ===================================================================
# Niche 距离计算
# ===================================================================

perform_niche_analysis <- function(seurat_obj, threshold, config) {
  cat("📈 开始 Niche 分析...\n")
  
  # ========== 新增：坐标列检查和标准化 ==========
  cat("🔍 检查并标准化空间坐标...\n")
  seurat_obj <- standardize_spatial_coordinates(seurat_obj)
  cat("✅ 坐标标准化完成\n\n")
  
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
    cat(sprintf(">> 使用核心数: %d\n", config$n_workers))
    cat(sprintf(">> 总点数: %d, 标记点数: %d (%.1f%%)\n", 
                ncol(seurat_obj), 
                sum(seurat_obj$ClockGene_High),
                100 * sum(seurat_obj$ClockGene_High) / ncol(seurat_obj)))
    
    # 设置并行计划
    plan(multisession, workers = config$n_workers)
    
    # ========== 新增：添加错误捕获 ==========
    result <- tryCatch({
      niche_marker(
        .data = seurat_obj,
        marker = ClockGene_High,
        spot_type = ClockGene_Distance,
        slide = orig.ident,
        dist_method = config$niche_dist_method,
        FUN = NA,
        n_work = config$n_workers
      )
    }, error = function(e) {
      cat("\n❌ Niche 分析出错！\n")
      cat(sprintf("错误信息: %s\n", e$message))
      
      # 诊断信息
      cat("\n🔍 诊断信息:\n")
      sample_names <- unique(seurat_obj$orig.ident)
      cat(sprintf("   样本数: %d\n", length(sample_names)))
      cat(sprintf("   样本列表: %s\n", paste(head(sample_names, 3), collapse=", ")))
      
      # 检查第一个样本的坐标
      if (length(sample_names) > 0) {
        first_sample <- sample_names[1]
        if (first_sample %in% names(seurat_obj@images)) {
          coords <- seurat_obj@images[[first_sample]]@coordinates
          cat(sprintf("   第一个样本的坐标列: %s\n", 
                      paste(colnames(coords), collapse=", ")))
        }
      }
      
      stop(sprintf("Niche 分析失败: %s", e$message))
    })
    
    seurat_obj <- result
    
    niche_data <- data.frame(ClockGene_Distance = seurat_obj$ClockGene_Distance)
    save_cache(niche_data, cache_file, "Niche 距离")
  }
  
  cat(sprintf("✅ 距离范围: %.2f ~ %.2f\n\n",
              min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
              max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))
  
  return(seurat_obj)
}