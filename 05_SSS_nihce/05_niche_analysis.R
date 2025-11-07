# 05_niche_analysis.R

#!/usr/bin/env Rscript
# ===================================================================
# Niche 距离计算
# ===================================================================

perform_niche_analysis <- function(seurat_obj, threshold, config) {


  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   Niche 距离分析\n")
  cat("═══════════════════════════════════════════════════════════\n\n")

     # ========================================
  # 🔍 类型验证（防止 filter 错误）
  # ========================================
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf(
      "❌ seurat_obj 必须是 Seurat 对象，实际类型: %s",
      class(seurat_obj)[1]
    ))
  }
  
  cat(sprintf("   ✅ 输入对象类型: %s\n", class(seurat_obj)[1]))
  cat(sprintf("   ✅ 细胞数量: %d\n", ncol(seurat_obj)))
  
  # 生成缓存key
  cache_key <- generate_cache_key(
    threshold, 
    sum(seurat_obj$ClockGene_High, na.rm = TRUE), 
    ncol(seurat_obj), 
    config$niche_dist_method
  )
  cache_file <- file.path(config$cache_dir, sprintf("niche_analysis_%s.rds", cache_key))
  
  # 检查缓存
  if (file.exists(cache_file)) {
    cat("📦 从缓存加载 Niche 距离数据...\n")
    niche_data <- load_cache(cache_file, "Niche 距离")
    
    # 验证缓存数据
    if (length(niche_data$ClockGene_Distance) != ncol(seurat_obj)) {
      warning("⚠️ 缓存数据大小不匹配，将重新计算")
      file.remove(cache_file)
    } else {
      seurat_obj$ClockGene_Distance <- niche_data$ClockGene_Distance
      
      cat(sprintf("✅ 距离范围: %.2f ~ %.2f\n",
                  min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
                  max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))
      cat("✅ Niche 分析完成（从缓存加载）\n\n")
      
      return(seurat_obj)
    }
  }
  
  # 如果没有缓存或缓存无效，进行计算
  cat("🔄 开始计算 Niche 距离...\n\n")
  
  # 数据统计
  n_total <- ncol(seurat_obj)
  n_marker <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
  cat(sprintf("数据概况:\n"))
  cat(sprintf("  总细胞数: %d\n", n_total))
  cat(sprintf("  标记细胞: %d (%.1f%%)\n", n_marker, 100 * n_marker / n_total))
  cat(sprintf("  使用核心数: %d\n", config$n_workers))
  cat(sprintf("  距离方法: %s\n\n", config$niche_dist_method))
  
  # 验证必需的列
  if (!"ClockGene_High" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Seurat 对象中缺少 'ClockGene_High' 列，请先运行 define_high_expression()")
  }
  
  if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Seurat 对象中缺少 'orig.ident' 列")
  }
  
  # 验证空间数据
  if (length(names(seurat_obj@images)) == 0) {
    stop("❌ Seurat 对象中没有空间图像数据")
  }
  
  # 调用 niche_marker 函数
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
    # 详细的错误诊断
    cat("\n")
    cat("═══════════════════════════════════════════════════════════\n")
    cat("   ❌ Niche 分析失败\n")
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    cat(sprintf("错误信息: %s\n\n", e$message))
    
    # 诊断信息
    cat("🔍 诊断信息:\n\n")
    
    # 样本信息
    sample_names <- unique(seurat_obj$orig.ident)
    cat(sprintf("1. 样本数量: %d\n", length(sample_names)))
    cat(sprintf("   样本列表（前5个）: %s\n\n", 
                paste(head(sample_names, 5), collapse=", ")))
    
    # 空间数据信息
    image_names <- names(seurat_obj@images)
    cat(sprintf("2. 空间图像数: %d\n", length(image_names)))
    
    if (length(image_names) > 0) {
      cat("   图像列表（前5个）:\n")
      for (img in head(image_names, 5)) {
        cat(sprintf("     - %s\n", img))
        if (img %in% names(seurat_obj@images)) {
          coords <- seurat_obj@images[[img]]@coordinates
          cat(sprintf("       细胞数: %d\n", nrow(coords)))
          cat(sprintf("       坐标列: %s\n", paste(colnames(coords), collapse=", ")))
        }
      }
    }
    cat("\n")
    
    # 标记细胞信息
    cat(sprintf("3. 标记细胞统计:\n"))
    marker_table <- table(seurat_obj$ClockGene_High)
    print(marker_table)
    cat("\n")
    
    # 按样本统计标记细胞
    cat("4. 各样本标记细胞数:\n")
    marker_by_sample <- table(
      seurat_obj$orig.ident[seurat_obj$ClockGene_High]
    )
    print(head(marker_by_sample, 10))
    cat("\n")
    
    cat("═══════════════════════════════════════════════════════════\n\n")
    
    # 抛出原始错误
    stop(sprintf("Niche 分析失败: %s", e$message))
  })
  
  seurat_obj <- result
  
  # 验证结果
  if (!"ClockGene_Distance" %in% colnames(seurat_obj@meta.data)) {
    stop("❌ Niche 分析未能生成 'ClockGene_Distance' 列")
  }
  
  if (any(is.na(seurat_obj$ClockGene_Distance))) {
    n_na <- sum(is.na(seurat_obj$ClockGene_Distance))
    warning(sprintf("⚠️ 警告：%d 个细胞的距离值为 NA", n_na))
  }
  
  # 保存缓存
  cat("\n💾 保存结果到缓存...\n")
  niche_data <- data.frame(
    ClockGene_Distance = seurat_obj$ClockGene_Distance,
    stringsAsFactors = FALSE
  )
  save_cache(niche_data, cache_file, "Niche 距离")
  
  # 输出结果摘要
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   结果摘要\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  dist_vals <- seurat_obj$ClockGene_Distance
  cat(sprintf("Distance 统计:\n"))
  cat(sprintf("  最小值: %.2f\n", min(dist_vals, na.rm = TRUE)))
  cat(sprintf("  第25百分位: %.2f\n", quantile(dist_vals, 0.25, na.rm = TRUE)))
  cat(sprintf("  中位数: %.2f\n", median(dist_vals, na.rm = TRUE)))
  cat(sprintf("  第75百分位: %.2f\n", quantile(dist_vals, 0.75, na.rm = TRUE)))
  cat(sprintf("  最大值: %.2f\n", max(dist_vals, na.rm = TRUE)))
  cat(sprintf("  平均值: %.2f\n", mean(dist_vals, na.rm = TRUE)))
  cat(sprintf("  标准差: %.2f\n", sd(dist_vals, na.rm = TRUE)))
  
  # 标记细胞的距离验证
  marker_dist <- dist_vals[seurat_obj$ClockGene_High]
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- sum(!is.na(marker_dist))
  
  cat(sprintf("\n标记细胞验证:\n"))
  cat(sprintf("  标记细胞数: %d\n", n_marker_total))
  cat(sprintf("  Distance=0: %d (%.1f%%)\n", 
              n_marker_zero, 
              100 * n_marker_zero / n_marker_total))
  
  if (n_marker_zero / n_marker_total < 0.95) {
    cat("\n⚠️ 警告：少于95%的标记细胞距离为0，可能存在问题\n")
  } else {
    cat("\n✅ 验证通过：标记细胞距离计算正确\n")
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("   Niche 分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  return(seurat_obj)
}


# ===================================================================
# 辅助函数：快速诊断（可选）
# ===================================================================

quick_diagnose_niche <- function(seurat_obj) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   快速 Niche 诊断\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 基本信息
  cat("1. 基本信息:\n")
  cat(sprintf("   总细胞数: %d\n", ncol(seurat_obj)))
  cat(sprintf("   总基因数: %d\n", nrow(seurat_obj)))
  
  # 检查必需的列
  cat("\n2. 必需列检查:\n")
  required_cols <- c("ClockGene_High", "orig.ident")
  for (col in required_cols) {
    exists <- col %in% colnames(seurat_obj@meta.data)
    cat(sprintf("   %s: %s\n", col, ifelse(exists, "✓", "✗")))
  }
  
  # 空间数据
  cat("\n3. 空间数据:\n")
  n_images <- length(names(seurat_obj@images))
  cat(sprintf("   图像数: %d\n", n_images))
  
  if (n_images > 0) {
    cat("   样本列表:\n")
    for (img in names(seurat_obj@images)) {
      coords <- seurat_obj@images[[img]]@coordinates
      cat(sprintf("     - %s: %d 个细胞, 坐标列: [%s]\n", 
                  img, 
                  nrow(coords),
                  paste(colnames(coords), collapse=", ")))
    }
  }
  
  # 标记细胞
  if ("ClockGene_High" %in% colnames(seurat_obj@meta.data)) {
    cat("\n4. 标记细胞:\n")
    n_marked <- sum(seurat_obj$ClockGene_High, na.rm = TRUE)
    cat(sprintf("   总数: %d (%.1f%%)\n", 
                n_marked, 
                100 * n_marked / ncol(seurat_obj)))
    
    # 按样本统计
    if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      cat("\n   各样本标记细胞数:\n")
      marker_by_sample <- table(
        seurat_obj$orig.ident[seurat_obj$ClockGene_High]
      )
      for (i in seq_along(marker_by_sample)) {
        sample_name <- names(marker_by_sample)[i]
        n_marker <- marker_by_sample[i]
        n_total <- sum(seurat_obj$orig.ident == sample_name)
        cat(sprintf("     - %s: %d/%d (%.1f%%)\n", 
                    sample_name, n_marker, n_total,
                    100 * n_marker / n_total))
      }
    }
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
  
  invisible(NULL)
}