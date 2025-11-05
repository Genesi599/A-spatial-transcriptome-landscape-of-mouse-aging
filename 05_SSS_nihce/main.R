#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis - Main Script
# Author: Zhangbin (optimized)
# Date: 2024-11-04
# ===================================================================

# 加载配置
source("00_config.R")

# 加载模块
source("01_setup.R")
source("02_utils.R")
source("03_load_data.R")
source("04_module_score.R")
source("05_niche_analysis.R")
source("06_plot_isoheight.R")
source("07_plot_spatial.R")
source("08_plot_celltype.R")
source("09_save_results.R")

# ===================================================================
# 主流程
# ===================================================================

main <- function() {
  # 1. 环境设置
  setup_environment(CONFIG)
  load_packages()
  load_custom_functions()
  
  # 2. 加载数据
  gene_list <- load_gene_list(CONFIG)
  seurat_obj <- load_seurat_object(CONFIG)
  genes_in_data <- check_gene_overlap(gene_list, seurat_obj)
  
  # 3. 计算评分
  seurat_obj <- calculate_module_score(seurat_obj, genes_in_data, CONFIG)
  result <- define_high_expression(seurat_obj, CONFIG)
  seurat_obj <- result$seurat_obj
  threshold <- result$threshold
  
  # 4. Niche分析
  seurat_obj <- perform_niche_analysis(seurat_obj, threshold, CONFIG)
  
  # 5. 确定要处理的样本
  samples <- unique(seurat_obj$orig.ident)
  samples_to_plot <- if (CONFIG$debug_mode) {
    head(samples, CONFIG$debug_sample_limit)
  } else {
    samples
  }
  cat(sprintf("📋 将处理 %d 个样本\n\n", length(samples_to_plot)))
  
  # 6. 绘图
  plot_isoheight_all(seurat_obj, samples_to_plot, CONFIG)
  plot_spatial_gradient(seurat_obj, samples_to_plot, CONFIG)
  plot_celltype_analysis(seurat_obj, samples_to_plot, CONFIG)
  
  # 7. 保存结果
  save_results(seurat_obj, CONFIG)
  print_summary(CONFIG)
}

# 运行主流程
if (!interactive()) {
  main()
}