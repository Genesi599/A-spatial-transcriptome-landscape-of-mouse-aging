#!/usr/bin/env Rscript
# ===================================================================
# 配置参数
# ===================================================================

CONFIG <- list(
  # ===== 路径设置 =====
  work_dir = "/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce",
  output_dir = "/dellstorage09/quj_lab/yanghang/spatial",
  
  # 数据路径
  gene_list_path = "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt",
  seurat_path = "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lung_2-25M.rds",
  
  # ===== 分析参数 =====
  threshold_quantile = 0.95,  # Top 10%
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

# 生成目录路径
CONFIG$cache_dir <- file.path(CONFIG$output_dir, "cache")
CONFIG$figure_dir <- file.path(CONFIG$output_dir, "figure")
CONFIG$metadata_dir <- file.path(CONFIG$output_dir, "metadata")

# 子目录
CONFIG$dirs <- list(
  cache = CONFIG$cache_dir,
  figure = CONFIG$figure_dir,
  metadata = CONFIG$metadata_dir,
  isoheight = file.path(CONFIG$figure_dir, "isoheight"),
  spatial = file.path(CONFIG$figure_dir, "spatial"),
  overlay = file.path(CONFIG$figure_dir, "isoheight", "01_overlay_plots"),
  celltype = file.path(CONFIG$figure_dir, "isoheight", "02_celltype_only"),
  composition = file.path(CONFIG$figure_dir, "isoheight", "03_composition_stats"),
  heatmaps = file.path(CONFIG$figure_dir, "isoheight", "04_heatmaps"),
  combined = file.path(CONFIG$figure_dir, "isoheight", "05_combined_analysis")
)