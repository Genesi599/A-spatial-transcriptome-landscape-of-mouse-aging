# 00_config.R

#!/usr/bin/env Rscript
# ===================================================================
# 配置参数
# ===================================================================

CONFIG <- list(
  # ===== 路径设置 =====
  work_dir = "/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce",
  output_base_dir = "/dellstorage09/quj_lab/yanghang/spatial",
  
  # 数据路径
  gene_list_path = "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt",
  
  # ===== 批量处理设置 =====
  batch_mode = FALSE,  # TRUE: 批量处理模式, FALSE: 单文件模式
  
  # 单文件模式（当 batch_mode = FALSE 时使用）
  seurat_path = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lung_2-25M.rds",
  
  # 批量模式（当 batch_mode = TRUE 时使用）
  seurat_dir = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds",
  seurat_pattern = "\\.rds$",  # 匹配所有 .rds 文件
  recursive_search = FALSE,    # 是否递归搜索子目录
  
  # 可选：指定要处理的文件（留空则处理所有匹配的文件）
  # 例如: c("Lung_2-25M.rds", "Heart_2-25M.rds")
  specific_files = NULL,
  
  # 可选：排除某些文件
  # 例如: c("test.rds", "backup.rds")
  exclude_files = NULL,

  # ===== 分析参数 =====
  threshold_quantile = 0.95,  # Top 5%
  niche_dist_method = "Euclidean",
  n_workers = 6,
  
  # ===== 绘图参数 =====
  plot = list(
    contour_bins = 8,
    point_size_bg = 0.3,
    point_size_top = 1.2,
    point_size_scatter = 2.5,
    contour_alpha = 0.25,
    interp_resolution = 200,
    expand_margin = 0.05,
    dpi = 300
  ),
  
  # ===== 调试参数 =====
  debug_mode = FALSE,
  debug_sample_limit = 3,
  save_full_object = FALSE,
  
  # ===== 缓存参数 =====
  cache_max_age_hours = NULL
)

# ===================================================================
# 打印配置信息
# ===================================================================

cat(sprintf("\n╔═══════════════════════════════════════════════════════════╗\n"))
cat(sprintf("║                    配置信息                                ║\n"))
cat(sprintf("╚═══════════════════════════════════════════════════════════╝\n\n"))

if (CONFIG$batch_mode) {
  cat(sprintf("运行模式: 批量处理\n"))
  cat(sprintf("输入目录: %s\n", CONFIG$seurat_dir))
  cat(sprintf("文件模式: %s\n", CONFIG$seurat_pattern))
  cat(sprintf("递归搜索: %s\n", ifelse(CONFIG$recursive_search, "是", "否")))
  
  if (!is.null(CONFIG$specific_files)) {
    cat(sprintf("指定文件: %s\n", paste(CONFIG$specific_files, collapse = ", ")))
  }
  
  if (!is.null(CONFIG$exclude_files)) {
    cat(sprintf("排除文件: %s\n", paste(CONFIG$exclude_files, collapse = ", ")))
  }
} else {
  cat(sprintf("运行模式: 单文件处理\n"))
  cat(sprintf("Seurat 文件: %s\n", basename(CONFIG$seurat_path)))
}

cat(sprintf("输出目录: %s\n", CONFIG$output_base_dir))
cat(sprintf("阈值分位数: %.2f (Top %.0f%%)\n", 
            CONFIG$threshold_quantile, 
            (1 - CONFIG$threshold_quantile) * 100))
cat(sprintf("并行工作数: %d\n", CONFIG$n_workers))
cat(sprintf("调试模式: %s\n", ifelse(CONFIG$debug_mode, "开启", "关闭")))
cat(sprintf("\n"))