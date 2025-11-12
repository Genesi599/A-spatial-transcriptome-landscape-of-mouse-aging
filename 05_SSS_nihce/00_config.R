# 00_config.R (支持多基因列表)

CONFIG <- list(
  # ===== 路径设置 =====
  work_dir = "/data/home/quj_lab/yanghang/...",
  output_base_dir = "/dellstorage09/quj_lab/yanghang/spatial",
  
  # ✅ 修改：基因列表可以是单个文件或目录
  gene_list_path = "/dellstorage09/quj_lab/yanghang/spatial/ref",
  gene_list_pattern = "\\.txt$",  # 匹配 .txt 文件
  
  cache_dir = "/dellstorage09/quj_lab/yanghang/spatial/cache",
  
  # ===== 批量处理设置 =====
  batch_mode = TRUE,
  seurat_path = "...",
  seurat_dir = "...",
  seurat_pattern = "\\.rds$",
  recursive_search = FALSE,
  specific_files = NULL,
  exclude_files = NULL,

  # ===== 分析参数 =====
  threshold_quantile = 0.95,
  niche_dist_method = "Euclidean",
  n_workers = 10,
  
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
  debug_mode = TRUE,
  debug_sample_limit = 2,
  save_full_object = FALSE,
  
  # ===== 缓存参数 =====
  cache_max_age_hours = NULL
)

# 缓存目录初始化
if (!is.null(CONFIG$cache_dir)) {
  if (!dir.exists(CONFIG$cache_dir)) {
    dir.create(CONFIG$cache_dir, recursive = TRUE, 
               showWarnings = FALSE)
    cat(sprintf("✅ 创建缓存目录: %s\n", CONFIG$cache_dir))
  }
}

# 配置打印
cat("\n╔════════════════════════════════════════════════╗\n")
cat("║                  配置信息                      ║\n")
cat("╚════════════════════════════════════════════════╝\n\n")

if (CONFIG$batch_mode) {
  cat("运行模式: 批量处理\n")
  cat(sprintf("输入目录: %s\n", CONFIG$seurat_dir))
  cat(sprintf("基因列表: %s\n", CONFIG$gene_list_path))  # ✅ 修改
} else {
  cat("运行模式: 单文件处理\n")
  cat(sprintf("Seurat: %s\n", basename(CONFIG$seurat_path)))
  cat(sprintf("基因列表: %s\n", CONFIG$gene_list_path))  # ✅ 修改
}

cat(sprintf("输出目录: %s\n", CONFIG$output_base_dir))
cat(sprintf("缓存目录: %s\n", 
            ifelse(is.null(CONFIG$cache_dir), 
                   "未设置", CONFIG$cache_dir)))
cat(sprintf("调试模式: %s\n\n", 
            ifelse(CONFIG$debug_mode, "开启", "关闭")))