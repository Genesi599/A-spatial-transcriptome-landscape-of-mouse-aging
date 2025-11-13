# 00_config.R

CONFIG <- list(
  work_dir = "/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce",
  output_base_dir = "/dellstorage09/quj_lab/yanghang/spatial",
  
  gene_list_path = "/dellstorage09/quj_lab/yanghang/spatial/ref/clock_genes.txt",
  
  cache_dir = "/dellstorage09/quj_lab/yanghang/spatial/cache",
  
  batch_mode = TRUE,
  seurat_path = NULL,
  seurat_dir = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds",
  seurat_pattern = "\\.rds$",
  recursive_search = FALSE,
  specific_files = NULL,
  exclude_files = NULL,
  score_column_name = "ClockGene_Score1",
  threshold_quantile = 0.95,
  niche_dist_method = "Euclidean",
  n_workers = 10,
  
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
  
  debug_mode = FALSE,
  debug_sample_limit = 2,
  save_full_object = FALSE,
  
  cache_max_age_hours = NULL
)

if (!is.null(CONFIG$cache_dir)) {
  if (!dir.exists(CONFIG$cache_dir)) {
    dir.create(CONFIG$cache_dir, recursive = TRUE, 
               showWarnings = FALSE)
    cat(sprintf("✅ 创建缓存目录: %s\n", CONFIG$cache_dir))
  }
}

cat("\n", strrep("=", 60), "\n")
cat("配置信息\n")
cat(strrep("=", 60), "\n\n")

if (CONFIG$batch_mode) {
  cat("运行模式: 批量处理\n")
  cat(sprintf("输入目录: %s\n", CONFIG$seurat_dir))
  cat(sprintf("基因列表: %s\n", 
              ifelse(file.exists(CONFIG$gene_list_path),
                     basename(CONFIG$gene_list_path),
                     "⚠️  文件不存在")))
} else {
  cat("运行模式: 单文件处理\n")
  cat(sprintf("Seurat: %s\n", basename(CONFIG$seurat_path)))
  cat(sprintf("基因列表: %s\n", basename(CONFIG$gene_list_path)))
}

cat(sprintf("输出目录: %s\n", CONFIG$output_base_dir))
cat(sprintf("缓存目录: %s\n", 
            ifelse(is.null(CONFIG$cache_dir), 
                   "未设置", CONFIG$cache_dir)))
cat(sprintf("调试模式: %s\n\n", 
            ifelse(CONFIG$debug_mode, "开启", "关闭")))