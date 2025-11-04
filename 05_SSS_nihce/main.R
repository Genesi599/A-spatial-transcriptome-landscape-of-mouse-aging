#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis Workflow (with Caching)
# Author: Zhangbin (adapted by ChatGPT)
# Date: 2025-11-04
# ===================================================================

# -----------------------------
# 0. 环境设置
# -----------------------------
# 设置当前工作目录为脚本所在路径
setwd("/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce")

output_dir <- "/dellstorage09/quj_lab/yanghang/spatial"
cache_dir <- "/dellstorage09/quj_lab/yanghang/spatial/cache"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cat("✅ 工作目录:", getwd(), "\n")
cat("✅ 输出目录:", output_dir, "\n")
cat("✅ 缓存目录:", cache_dir, "\n")

# -----------------------------
# 1. 加载 R 包
# -----------------------------
cat("\n📦 加载所需R包...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(proxy)
  library(future)
  library(future.apply)
  library(ggnewscale)
  library(RColorBrewer)
  library(patchwork)
  library(digest)  # 用于生成缓存 key
})

# -----------------------------
# 🔧 缓存工具函数
# -----------------------------

# 生成缓存 key (基于参数的 MD5 哈希)
generate_cache_key <- function(...) {
  params <- list(...)
  key <- digest::digest(params, algo = "md5")
  return(key)
}

# 保存缓存
save_cache <- function(obj, cache_file, description = "") {
  tryCatch({
    saveRDS(obj, cache_file)
    cat(sprintf("💾 缓存已保存: %s\n", basename(cache_file)))
    if (description != "") cat(sprintf("   描述: %s\n", description))
  }, error = function(e) {
    warning(sprintf("⚠️ 缓存保存失败: %s\n", e$message))
  })
}

# 加载缓存
load_cache <- function(cache_file, description = "") {
  if (file.exists(cache_file)) {
    cat(sprintf("📂 发现缓存文件: %s\n", basename(cache_file)))
    if (description != "") cat(sprintf("   描述: %s\n", description))
    obj <- readRDS(cache_file)
    cat("✅ 缓存加载成功！\n")
    return(obj)
  } else {
    return(NULL)
  }
}

# 检查缓存是否有效
is_cache_valid <- function(cache_file, source_file = NULL, max_age_hours = NULL) {
  if (!file.exists(cache_file)) return(FALSE)
  
  # 检查源文件是否更新
  if (!is.null(source_file) && file.exists(source_file)) {
    if (file.mtime(cache_file) < file.mtime(source_file)) {
      cat("⚠️ 源文件已更新，缓存失效\n")
      return(FALSE)
    }
  }
  
  # 检查缓存年龄
  if (!is.null(max_age_hours)) {
    cache_age <- difftime(Sys.time(), file.mtime(cache_file), units = "hours")
    if (cache_age > max_age_hours) {
      cat(sprintf("⚠️ 缓存已过期 (%.1f 小时)\n", cache_age))
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# -----------------------------
# 2. 加载函数脚本
# -----------------------------
source("niche_marker.R")
source("SSS_isoheight_plot.R")

# -----------------------------
# 3. 读取基因列表 (缓存)
# -----------------------------
cat("\n📄 读取基因列表...\n")
gene_list_path <- "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt"

# 缓存文件名
gene_list_cache <- file.path(cache_dir, "gene_list.rds")

if (is_cache_valid(gene_list_cache, gene_list_path)) {
  gene_list <- load_cache(gene_list_cache, "基因列表")
} else {
  gene_list <- read.table(gene_list_path, header = TRUE, stringsAsFactors = FALSE)[[1]]
  gene_list <- trimws(gene_list)
  gene_list <- gene_list[gene_list != ""]
  save_cache(gene_list, gene_list_cache, "基因列表")
}

cat(sprintf("✅ 共读取 %d 个基因。\n", length(gene_list)))
print(head(gene_list))

# -----------------------------
# 4. 加载 Seurat 对象 (缓存)
# -----------------------------
cat("\n🧠 加载 Seurat 对象...\n")
seurat_path <- "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lymph_2-25M.rds"

# 缓存文件名
seurat_cache <- file.path(cache_dir, "Lymph_2-25M_original.rds")

if (is_cache_valid(seurat_cache, seurat_path)) {
  seurat_obj <- load_cache(seurat_cache, "Seurat 对象")
} else {
  seurat_obj <- readRDS(seurat_path)
  save_cache(seurat_obj, seurat_cache, "Seurat 对象")
}

cat(sprintf("✅ Spots 数量: %d, 基因数量: %d\n", ncol(seurat_obj), nrow(seurat_obj)))

# --------------------------------------------------
# 🧬 验证 Seurat 对象基因名格式
# --------------------------------------------------
cat("\n🔍 检查 Seurat 对象中的基因名格式...\n")
gene_names_preview <- head(rownames(seurat_obj), 100)
cat("📄 前 10 个基因名示例:\n")
print(head(gene_names_preview, 10))
cat("✅ 基因名检查完成。\n")

# -----------------------------
# 5. 检查基因
# -----------------------------
genes_in_data <- intersect(gene_list, rownames(seurat_obj))
genes_missing <- setdiff(gene_list, rownames(seurat_obj))

cat(sprintf("✅ 匹配上的基因数: %d / %d (%.1f%%)\n",
            length(genes_in_data),
            length(gene_list),
            100 * length(genes_in_data) / length(gene_list)))

if (length(genes_in_data) < length(gene_list) * 0.3) {
  upper_match <- sum(toupper(gene_list) %in% toupper(rownames(seurat_obj)))
  if (upper_match > length(genes_in_data)) {
    cat("💡 提示：基因名大小写可能不一致，可尝试统一大写。\n")
  }
}

if (length(genes_missing) > 0) {
  cat("⚠️ 以下部分基因不在数据中（将被忽略）:\n")
  print(utils::head(genes_missing, 15))
  if (length(genes_missing) > 15)
    cat(sprintf("... 其余 %d 个未显示\n", length(genes_missing) - 15))
}

# -----------------------------
# 6. 计算基因集评分 (Module Score) - 缓存
# -----------------------------
cat("\n🧮 计算 Clock Gene Module Score...\n")

# 生成缓存 key
score_cache_key <- generate_cache_key(
  genes = genes_in_data,
  seurat_dims = dim(seurat_obj),
  method = "AddModuleScore"
)
score_cache_file <- file.path(cache_dir, sprintf("module_score_%s.rds", score_cache_key))

if (file.exists(score_cache_file)) {
  cat("📂 发现 Module Score 缓存...\n")
  score_data <- load_cache(score_cache_file, "Module Score 评分")
  seurat_obj$ClockGene_Score1 <- score_data$ClockGene_Score1
} else {
  cat("🔄 正在计算 Module Score (可能需要几分钟)...\n")
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = list(clock_gene_set = genes_in_data),
    name = "ClockGene_Score"
  )
  
  # 保存缓存
  score_data <- data.frame(ClockGene_Score1 = seurat_obj$ClockGene_Score1)
  save_cache(score_data, score_cache_file, "Module Score 评分")
}

# 定义全局阈值
threshold <- quantile(seurat_obj$ClockGene_Score1, 0.7, na.rm = TRUE)
cat(sprintf("✅ 高表达阈值: %.3f (Top 30%%)\n", threshold))

# 创建辅助列
seurat_obj$ClockGene_High <- seurat_obj$ClockGene_Score1 > threshold
cat("✅ 高/低表达分组:\n")
print(table(seurat_obj$ClockGene_High))

# -----------------------------
# 7. Niche 分析 - 缓存
# -----------------------------
cat("\n📈 开始 Niche 分析...\n")

# 生成缓存 key
niche_cache_key <- generate_cache_key(
  threshold = threshold,
  high_spots = sum(seurat_obj$ClockGene_High),
  total_spots = ncol(seurat_obj),
  dist_method = "Euclidean"
)
niche_cache_file <- file.path(cache_dir, sprintf("niche_analysis_%s.rds", niche_cache_key))

if (file.exists(niche_cache_file)) {
  cat("📂 发现 Niche 分析缓存...\n")
  niche_data <- load_cache(niche_cache_file, "Niche 距离数据")
  seurat_obj$ClockGene_Distance <- niche_data$ClockGene_Distance
} else {
  cat("🔄 正在进行 Niche 分析 (可能需要较长时间)...\n")
  
  library(future)
  plan(multisession, workers = 6)
  
  seurat_obj <- niche_marker(
    .data = seurat_obj,
    marker = ClockGene_High,
    spot_type = ClockGene_Distance,
    slide = orig.ident,
    dist_method = "Euclidean",
    FUN = NA,
    n_work = 6
  )
  
  # 保存缓存
  niche_data <- data.frame(ClockGene_Distance = seurat_obj$ClockGene_Distance)
  save_cache(niche_data, niche_cache_file, "Niche 距离数据")
}

cat("✅ Niche 分析完成。\n")

# -----------------------------
# 8. 绘制 Isoheight 等高线图 - 缓存
# -----------------------------
cat("\n🎨 绘制 Isoheight 图...\n")

iso_plot_cache <- file.path(cache_dir, sprintf("isoheight_plot_%s.rds", niche_cache_key))

if (file.exists(iso_plot_cache)) {
  cat("📂 发现等高线图缓存...\n")
  p_iso_col <- load_cache(iso_plot_cache, "等高线图对象")
} else {
  cat("🔄 正在生成等高线图...\n")
  
  p_iso_col <- celltype_isoheight_plot(
    .data = seurat_obj,
    density_top = ClockGene_High,
    col_bg = "gray92",
    col_top = "#d62728",
    col_isoheight = "white",
    col_white_ratio = 0.25,
    cols_fill_isoheight = c(
      rep("white", 25),
      colorRampPalette(brewer.pal(9, "YlOrRd")[3:9])(75)
    ),
    size_bg = 0.3,
    size_top = 0.7,
    nrow = 2
  )
  
  save_cache(p_iso_col, iso_plot_cache, "等高线图对象")
}

# 保存 PDF
ggsave(
  file.path(output_dir, "ClockGene_niche_isoheight.pdf"),
  plot = p_iso_col,
  width = 14, height = 10, dpi = 300
)

cat("✅ 等高线图已保存。\n")

# -----------------------------
# 9. 可视化 Niche 距离梯度
# -----------------------------
cat("\n🔥 绘制空间梯度图...\n")

p_score <- SpatialFeaturePlot(
  seurat_obj,
  features = "ClockGene_Score1",
  pt.size.factor = 1.5,
  alpha = c(0.1, 1)
) + scale_fill_gradientn(
  colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
  name = "Clock Gene\nScore"
)

p_niche <- SpatialFeaturePlot(
  seurat_obj,
  features = "ClockGene_Distance",  # ✅ 修改：改为 ClockGene_Distance
  pt.size.factor = 1.5,
  alpha = c(0.1, 1)
) + scale_fill_gradientn(
  colors = rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582",
                 "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de")),
  name = "Distance to\nHigh Score Region"
)

p_combined <- (p_score | p_niche) +
  plot_annotation(title = "Clock Gene Niche Analysis (Lymph_2-25M)")

ggsave(file.path(output_dir, "ClockGene_combined_spatial.pdf"),
       plot = p_combined, width = 18, height = 9, dpi = 300)
# -----------------------------
# 10. 保存结果
# -----------------------------
cat("\n💾 保存结果...\n")

saveRDS(seurat_obj, file.path(output_dir, "Lymph_2-25M_with_clockgene_niche.rds"))

write.csv(seurat_obj@meta.data,
          file.path(output_dir, "Lymph_2-25M_clockgene_metadata.csv"),
          row.names = TRUE)

# -----------------------------
# 11. 缓存管理信息
# -----------------------------
cat("\n📊 缓存统计:\n")
cache_files <- list.files(cache_dir, full.names = TRUE)
if (length(cache_files) > 0) {
  cache_info <- file.info(cache_files)
  total_size <- sum(cache_info$size) / 1024^2  # MB
  cat(sprintf("   缓存文件数: %d\n", length(cache_files)))
  cat(sprintf("   总大小: %.1f MB\n", total_size))
  cat(sprintf("   位置: %s\n", cache_dir))
  cat("\n💡 提示: 如需清除缓存，运行: unlink(cache_dir, recursive = TRUE)\n")
}

cat("\n✅ 全部完成！结果保存在:\n")
cat("   输出目录:", output_dir, "\n")
cat("   缓存目录:", cache_dir, "\n")