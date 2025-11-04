#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis Workflow (Optimized)
# Author: Zhangbin (optimized by Assistant)
# Date: 2025-11-04
# 优化点：
#   1. 所有图片保存到 figure 子文件夹
#   2. 简化代码结构，保持功能不变
#   3. 优化缓存机制
# ===================================================================

# -----------------------------
# 0. 环境设置
# -----------------------------
setwd("/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce")

# 主输出目录
output_dir <- "/dellstorage09/quj_lab/yanghang/spatial"
cache_dir <- file.path(output_dir, "cache")
figure_dir <- file.path(output_dir, "figure")  # ✅ 统一的图形文件夹

# ✅ 创建目录结构
dirs <- list(
  cache = cache_dir,
  figure = figure_dir,
  isoheight = file.path(figure_dir, "isoheight"),
  spatial = file.path(figure_dir, "spatial"),
  sss_niche = file.path(figure_dir, "sss_niche"),
  metadata = file.path(output_dir, "metadata")
)

lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

cat("✅ 工作目录:", getwd(), "\n")
cat("✅ 输出目录:", output_dir, "\n")
cat("✅ 图形目录:", figure_dir, "\n")
cat("✅ 缓存目录:", cache_dir, "\n")

# -----------------------------
# 1. 加载 R 包
# -----------------------------
cat("\n📦 加载所需R包...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(future)
  library(future.apply)
  library(ggnewscale)
  library(RColorBrewer)
  library(patchwork)
  library(digest)
  library(akima)
})

# -----------------------------
# 2. 缓存工具函数
# -----------------------------
generate_cache_key <- function(...) digest::digest(list(...), algo = "md5")

save_cache <- function(obj, file, desc = "") {
  tryCatch({
    saveRDS(obj, file)
    cat(sprintf("💾 缓存已保存: %s (%.2f MB) %s\n", 
                basename(file), file.size(file)/1024^2, 
                ifelse(desc != "", paste0("- ", desc), "")))
  }, error = function(e) warning(sprintf("⚠️ 缓存保存失败: %s\n", e$message)))
}

load_cache <- function(file, desc = "") {
  if (file.exists(file)) {
    cat(sprintf("📂 加载缓存: %s (%.2f MB) %s\n", 
                basename(file), file.size(file)/1024^2,
                ifelse(desc != "", paste0("- ", desc), "")))
    return(readRDS(file))
  }
  NULL
}

is_cache_valid <- function(cache_file, source_file = NULL, max_age_hours = NULL) {
  if (!file.exists(cache_file)) return(FALSE)
  if (!is.null(source_file) && file.exists(source_file)) {
    if (file.mtime(cache_file) < file.mtime(source_file)) {
      cat("⚠️ 源文件已更新，缓存失效\n")
      return(FALSE)
    }
  }
  if (!is.null(max_age_hours)) {
    age <- difftime(Sys.time(), file.mtime(cache_file), units = "hours")
    if (age > max_age_hours) {
      cat(sprintf("⚠️ 缓存已过期 (%.1f 小时)\n", age))
      return(FALSE)
    }
  }
  TRUE
}

# -----------------------------
# 3. 加载函数脚本
# -----------------------------
cat("\n📚 加载分析函数...\n")
source("niche_marker.R")
source("SSS_isoheight_plot.R")

# -----------------------------
# 4. 读取基因列表
# -----------------------------
cat("\n📄 读取基因列表...\n")
gene_list_path <- "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt"
gene_list_cache <- file.path(cache_dir, "gene_list.rds")

if (is_cache_valid(gene_list_cache, gene_list_path)) {
  gene_list <- load_cache(gene_list_cache, "基因列表")
} else {
  gene_list <- read.table(gene_list_path, header = TRUE, stringsAsFactors = FALSE)[[1]]
  gene_list <- trimws(gene_list[gene_list != ""])
  save_cache(gene_list, gene_list_cache, "基因列表")
}
cat(sprintf("✅ 共读取 %d 个基因\n", length(gene_list)))

# -----------------------------
# 5. 加载 Seurat 对象
# -----------------------------
cat("\n🧠 加载 Seurat 对象...\n")
seurat_path <- "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lymph_2-25M.rds"

seurat_obj <- readRDS(seurat_path)
seurat_obj <- UpdateSeuratObject(seurat_obj)
cat(sprintf("✅ Spots: %d, Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))

# -----------------------------
# 6. 检查基因匹配
# -----------------------------
cat("\n🔍 检查基因匹配情况...\n")
genes_in_data <- intersect(gene_list, rownames(seurat_obj))
genes_missing <- setdiff(gene_list, rownames(seurat_obj))

cat(sprintf("✅ 匹配基因: %d / %d (%.1f%%)\n",
            length(genes_in_data), length(gene_list),
            100 * length(genes_in_data) / length(gene_list)))

if (length(genes_missing) > 0 && length(genes_missing) <= 10) {
  cat("⚠️ 缺失基因:", paste(genes_missing, collapse = ", "), "\n")
} else if (length(genes_missing) > 10) {
  cat(sprintf("⚠️ 缺失 %d 个基因 (前10个): %s ...\n", 
              length(genes_missing), paste(head(genes_missing, 10), collapse = ", ")))
}

# -----------------------------
# 7. 计算基因集评分
# -----------------------------
cat("\n🧮 计算 Clock Gene Module Score...\n")
score_cache_key <- generate_cache_key(genes_in_data, dim(seurat_obj), "AddModuleScore")
score_cache_file <- file.path(cache_dir, sprintf("module_score_%s.rds", score_cache_key))

if (file.exists(score_cache_file)) {
  score_data <- load_cache(score_cache_file, "Module Score")
  seurat_obj$ClockGene_Score1 <- score_data$ClockGene_Score1
} else {
  cat("🔄 正在计算 Module Score...\n")
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = list(clock_gene_set = genes_in_data),
    name = "ClockGene_Score"
  )
  score_data <- data.frame(ClockGene_Score1 = seurat_obj$ClockGene_Score1)
  save_cache(score_data, score_cache_file, "Module Score")
}

cat(sprintf("✅ 评分范围: %.3f ~ %.3f\n", 
            min(seurat_obj$ClockGene_Score1, na.rm = TRUE),
            max(seurat_obj$ClockGene_Score1, na.rm = TRUE)))

# -----------------------------
# 8. 设置阈值
# -----------------------------
THRESHOLD_QUANTILE <- 0.90  # Top 10%
threshold <- quantile(seurat_obj$ClockGene_Score1, THRESHOLD_QUANTILE, na.rm = TRUE)
threshold_pct <- (1 - THRESHOLD_QUANTILE) * 100
threshold_desc <- sprintf("Top %.1f%%", threshold_pct)

cat(sprintf("✅ 高表达阈值: %.3f (%s)\n", threshold, threshold_desc))

seurat_obj$ClockGene_High <- seurat_obj$ClockGene_Score1 > threshold
cat("✅ 高/低表达分组:\n")
print(table(seurat_obj$ClockGene_High))

# -----------------------------
# 9. Niche 分析
# -----------------------------
cat("\n📈 开始 Niche 分析...\n")
niche_cache_key <- generate_cache_key(threshold, sum(seurat_obj$ClockGene_High), 
                                      ncol(seurat_obj), "Euclidean")
niche_cache_file <- file.path(cache_dir, sprintf("niche_analysis_%s.rds", niche_cache_key))

if (file.exists(niche_cache_file)) {
  niche_data <- load_cache(niche_cache_file, "Niche 距离")
  seurat_obj$ClockGene_Distance <- niche_data$ClockGene_Distance
} else {
  cat("🔄 正在进行 Niche 分析（多线程）...\n")
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
  
  niche_data <- data.frame(ClockGene_Distance = seurat_obj$ClockGene_Distance)
  save_cache(niche_data, niche_cache_file, "Niche 距离")
}

cat(sprintf("✅ 距离范围: %.2f ~ %.2f\n",
            min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
            max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))

# -----------------------------
# 10. 绘图配置
# -----------------------------
DEBUG_MODE <- TRUE  # 改为 FALSE 绘制所有样本
DEBUG_SAMPLE_LIMIT <- 3

samples <- unique(seurat_obj$orig.ident)
cat(sprintf("\n✅ 检测到 %d 个样本\n", length(samples)))

samples_to_plot <- if (DEBUG_MODE) {
  cat(sprintf("🔧 调试模式：只处理前 %d 个样本\n", DEBUG_SAMPLE_LIMIT))
  head(samples, DEBUG_SAMPLE_LIMIT)
} else {
  cat("🚀 生产模式：处理所有样本\n")
  samples
}

# -----------------------------
# 11. 绘制 Isoheight 图
# -----------------------------
cat("\n🎨 绘制 Isoheight 图...\n")

for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
  
  seurat_subset <- tryCatch(
    subset(seurat_obj, subset = orig.ident == sample_id),
    error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
  )
  
  p_iso <- celltype_isoheight_plot(
    .data = seurat_subset,
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
    size_top = 1.2,
    nrow = 1
  )
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  ggsave(file.path(dirs$isoheight, sprintf("ClockGene_isoheight_%s.pdf", safe_name)),
         plot = p_iso, width = 8, height = 8, dpi = 300)
}

# -----------------------------
# 12. 绘制空间梯度图（修复版）
# -----------------------------
cat("\n🔥 绘制空间梯度图...\n")

for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
  
  seurat_subset <- tryCatch(
    subset(seurat_obj, subset = orig.ident == sample_id),
    error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
  )
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # Score 图
  p_score <- SpatialFeaturePlot(
    seurat_subset, features = "ClockGene_Score1",
    pt.size.factor = 1.5, alpha = c(0.1, 1)
  ) + scale_fill_gradientn(
    colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
    name = "Clock Gene\nScore"
  ) + ggtitle(sample_id) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Distance 图（✅ 修复版：使用单向渐变 + 标记高表达区）
  p_niche <- SpatialFeaturePlot(
    seurat_subset, features = "ClockGene_Distance",
    pt.size.factor = 1.5, alpha = c(0.1, 1)
  ) + scale_fill_gradient(
    low = "#d73027",   # 红色 = 近（Distance=0，高表达核心）
    high = "#313695",  # 深蓝 = 远（Distance大）
    name = "Distance\nto High\nScore Region"
  ) + ggtitle(sample_id) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # 合并图
  p_combined <- (p_score | p_niche) +
    plot_annotation(
      title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
    )
  
  ggsave(file.path(dirs$spatial, sprintf("ClockGene_spatial_%s.pdf", safe_name)),
         plot = p_combined, width = 18, height = 9, dpi = 300)
}

# -----------------------------
# 13. 绘制 SSS Niche 热图
# -----------------------------
cat("\n🎨 绘制 SSS Niche 热图（平滑插值）...\n")

for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  tryCatch({
    # 提取数据
    sample_meta <- seurat_obj@meta.data %>%
      filter(orig.ident == sample_id) %>%
      rownames_to_column("cellid")
    
    # 检查必需列
    if (!all(c("col", "row", "ClockGene_High", "ClockGene_Distance") %in% colnames(sample_meta))) {
      cat("   ⚠️ 缺少必需列，跳过\n")
      next
    }
    
    # 空间插值
    col_range <- range(sample_meta$col, na.rm = TRUE)
    row_range <- range(sample_meta$row, na.rm = TRUE)
    
    interp_result <- akima::interp(
      x = sample_meta$col, y = sample_meta$row, z = sample_meta$ClockGene_Distance,
      xo = seq(col_range[1], col_range[2], length.out = 200),
      yo = seq(row_range[1], row_range[2], length.out = 200),
      linear = FALSE, extrap = FALSE
    )
    
    interp_df <- expand.grid(col = interp_result$x, row = interp_result$y) %>%
      mutate(distance = as.vector(interp_result$z)) %>%
      filter(!is.na(distance))
    
    # 绘图
    n_high <- sum(sample_meta$ClockGene_High, na.rm = TRUE)
    
    p_sss <- ggplot() +
      geom_raster(data = interp_df, aes(x = col, y = row, fill = distance), 
                  interpolate = TRUE) +
      scale_fill_gradientn(
        colours = c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", 
                    "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac"),
        name = "Distance\n(bins)", na.value = "gray95",
        guide = guide_colorbar(barwidth = 1.5, barheight = 10)
      ) +
      geom_point(data = sample_meta, aes(x = col, y = row), 
                 color = "white", size = 0.8, alpha = 0.6) +
      geom_point(data = filter(sample_meta, ClockGene_High), aes(x = col, y = row),
                 color = "black", size = 2.5, alpha = 0.9) +
      scale_x_continuous(expand = expansion(mult = 0.02)) +
      scale_y_reverse(expand = expansion(mult = 0.02)) +
      coord_fixed(ratio = 1) +
      labs(title = sample_id,
           subtitle = sprintf("⚫ SSS: %d spots (%.1f%%) | ⚪ All spots: %d",
                              n_high, 100 * n_high / nrow(sample_meta), nrow(sample_meta))) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(dirs$sss_niche, sprintf("ClockGene_SSS_niche_%s.pdf", safe_name)),
           plot = p_sss, width = 10, height = 10, dpi = 300)
    
  }, error = function(e) {
    cat(sprintf("   ❌ 绘制失败: %s\n", conditionMessage(e)))
  })
}

# -----------------------------
# 14. 保存结果
# -----------------------------
cat("\n💾 保存结果...\n")
write.csv(seurat_obj@meta.data, 
          file.path(dirs$metadata, "Lymph_2-25M_clockgene_metadata.csv"),
          row.names = TRUE)

# 可选：保存完整对象
save_full_object <- FALSE
if (save_full_object) {
  saveRDS(seurat_obj, file.path(dirs$metadata, "Lymph_2-25M_with_clockgene_niche.rds"))
}

# -----------------------------
# 15. 统计信息
# -----------------------------
cat("\n📊 文件统计:\n")
cat(sprintf("   图形文件夹: %s\n", figure_dir))
cat(sprintf("   - Isoheight: %d 个文件\n", length(list.files(dirs$isoheight))))
cat(sprintf("   - Spatial: %d 个文件\n", length(list.files(dirs$spatial))))
cat(sprintf("   - SSS Niche: %d 个文件\n", length(list.files(dirs$sss_niche))))

cat("\n✅ 全部完成！\n")
cat(sprintf("📁 所有图形已保存到: %s\n", figure_dir))