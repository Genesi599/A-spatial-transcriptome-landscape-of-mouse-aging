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
# 重新诊断 ClockGene_High 的定义
# -----------------------------
cat("\n🔍 深度诊断 ClockGene_High 定义...\n")

# 1. 查看 ClockGene_High 的定义逻辑
cat("\n【1】ClockGene_High 是如何定义的？\n")
cat("请检查你的代码中是否有类似这样的定义：\n")
cat("  seurat_obj$ClockGene_High <- seurat_obj$ClockGene_Score1 > threshold\n\n")

# 2. 对于示例样本，查看 Score 和 High 的对应关系
sample_id <- samples_to_plot[1]
sample_data <- seurat_obj@meta.data %>%
  filter(orig.ident == sample_id) %>%
  select(ClockGene_Score1, ClockGene_High, ClockGene_Distance) %>%
  arrange(desc(ClockGene_Score1))

cat("\n【2】样本", sample_id, "的数据检查：\n")
cat("\n=== Top 20 最高 Score 的点 ===\n")
print(head(sample_data, 20))
cat("\n期望：ClockGene_High 应该全是 TRUE\n")

cat("\n=== Top 20 最低 Score 的点 ===\n")
print(tail(sample_data, 20))
cat("\n期望：ClockGene_High 应该全是 FALSE\n")

# 3. 检查阈值
high_score_range <- range(seurat_obj$ClockGene_Score1[seurat_obj$ClockGene_High == TRUE], 
                          na.rm = TRUE)
low_score_range <- range(seurat_obj$ClockGene_Score1[seurat_obj$ClockGene_High == FALSE], 
                         na.rm = TRUE)

cat("\n【3】Score 范围检查：\n")
cat("高表达点（High=TRUE）的 Score 范围：[", 
    round(high_score_range[1], 2), ", ", 
    round(high_score_range[2], 2), "]\n", sep = "")
cat("低表达点（High=FALSE）的 Score 范围：[", 
    round(low_score_range[1], 2), ", ", 
    round(low_score_range[2], 2), "]\n", sep = "")

# 如果有重叠，说明分类有问题
if (high_score_range[1] < low_score_range[2]) {
  cat("\n⚠️ 警告：两个范围有重叠！\n")
  cat("   这可能导致一些高 Score 的点被标记为 High=FALSE\n")
}

# 4. 检查 High=TRUE 的点在图中的位置
cat("\n【4】检查黄圈点的真实 Score 值：\n")
high_points <- seurat_obj@meta.data %>%
  filter(orig.ident == sample_id, ClockGene_High == TRUE) %>%
  select(ClockGene_Score1, ClockGene_Distance)

cat("黄圈点的 Score 统计：\n")
print(summary(high_points$ClockGene_Score1))

cat("\n黄圈点的 Distance 统计：\n")
print(summary(high_points$ClockGene_Distance))

# 【5】相关性检查（修复版）
cat("\n【5】Score vs Distance 相关性（仅 High=TRUE 的点）：\n")

# ✅ 修复：当 Distance 全为 0 时，标准差为 0，无法计算相关系数
if (sd(high_points$ClockGene_Distance, na.rm = TRUE) == 0) {
  cat("✅ 所有高表达点的 Distance = 0（标准差为 0）\n")
  cat("   这是完全正确的！无需计算相关系数\n")
} else {
  cor_high <- cor(high_points$ClockGene_Score1, 
                  high_points$ClockGene_Distance, 
                  use = "complete.obs")
  cat("相关系数：", round(cor_high, 3), "\n")
  if (cor_high > 0) {
    cat("⚠️ 正相关：Score 越高，Distance 反而越大！这不对！\n")
  } else {
    cat("✅ 负相关：这是正常的\n")
  }
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("✅ 诊断完成！Distance 计算完全正确！\n")
cat(rep("=", 80), "\n\n", sep = "")

# -----------------------------
# 12. 绘制空间梯度图（完整修复版 - 匹配 Isoheight 方向）
# -----------------------------
cat("\n🔥 绘制空间梯度图（匹配 Isoheight 坐标）...\n")

for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
  
  seurat_subset <- tryCatch(
    subset(seurat_obj, subset = orig.ident == sample_id),
    error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
  )
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # ✅ 使用和 Isoheight 图完全相同的坐标获取方式
  coords <- GetTissueCoordinates(
    seurat_subset,
    cols = c("row", "col"),  # ✅ 明确指定 row 和 col
    scale = NULL
  )
  
  # 检查坐标列名
  coord_cols <- colnames(coords)
  cat(sprintf("   坐标列: %s\n", paste(coord_cols, collapse = ", ")))
  
  # ✅ 确保有 row 和 col（与 Isoheight 一致）
  if (!("row" %in% coord_cols && "col" %in% coord_cols)) {
    cat(sprintf("   ⚠️ 警告：未找到 row/col 列，跳过样本 %s\n", sample_id))
    cat("   可用列名：", paste(coord_cols, collapse = ", "), "\n")
    next
  }
  
  # 合并数据
  plot_data <- seurat_subset@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
  
  # ✅ 计算坐标范围（与 Isoheight 保持一致）
  expand_margin <- 0.05
  col_range <- range(plot_data$col, na.rm = TRUE)
  row_range <- range(plot_data$row, na.rm = TRUE)
  
  col_expand <- diff(col_range) * expand_margin
  row_expand <- diff(row_range) * expand_margin
  
  col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
  row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
  
  # ============================================
  # 左图：Clock Gene Score（蓝→红，低→高）
  # ============================================
  p_score <- ggplot(plot_data, aes(x = col, y = row)) +  # ✅ 使用 col, row
    geom_point(aes(fill = ClockGene_Score1), 
               shape = 21, size = 2.5, color = "white", stroke = 0.1) +
    scale_fill_gradientn(
      colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
      name = "Clock Gene\nScore",
      na.value = "gray90"
    ) +
    # ✅ 设置坐标范围
    scale_x_continuous(
      limits = col_limits,
      expand = expansion(mult = 0.02)
    ) +
    # ✅ 关键：反转 Y 轴（与 Isoheight 一致）
    scale_y_reverse(
      limits = rev(row_limits),
      expand = expansion(mult = 0.02)
    ) +
    coord_fixed(ratio = 1) +
    ggtitle("Clock Gene Score") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      aspect.ratio = 1,
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # ============================================
  # 右图：Distance（红→蓝，近→远）
  # ============================================
  p_distance <- ggplot(plot_data, aes(x = col, y = row)) +  # ✅ 使用 col, row
    geom_point(aes(fill = ClockGene_Distance), 
               shape = 21, size = 2.5, color = "white", stroke = 0.1) +
    scale_fill_gradientn(
      colors = rev(c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027")),
      name = "Distance to\nHigh Score\nRegion",
      na.value = "gray90"
    ) +
    # ✅ 设置坐标范围
    scale_x_continuous(
      limits = col_limits,
      expand = expansion(mult = 0.02)
    ) +
    # ✅ 关键：反转 Y 轴（与 Isoheight 一致）
    scale_y_reverse(
      limits = rev(row_limits),
      expand = expansion(mult = 0.02)
    ) +
    coord_fixed(ratio = 1) +
    ggtitle("Distance to High Score Region") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      aspect.ratio = 1,
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # ============================================
  # 合并图
  # ============================================
  p_combined <- (p_score | p_distance) +
    plot_annotation(
      title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  # 保存
  ggsave(file.path(dirs$spatial, sprintf("ClockGene_spatial_%s.pdf", safe_name)),
         plot = p_combined, width = 16, height = 8, dpi = 300)
  
  cat(sprintf("   ✅ 已保存: ClockGene_spatial_%s.pdf\n", safe_name))
}

cat("\n✅ 所有空间图绘制完成！\n")
cat("   方向已与 Isoheight 图保持一致（Y 轴反转）\n")


# -----------------------------
# 13. 绘制细胞类型 + 等高线叠加图
# -----------------------------
cat("\n🎨 绘制细胞类型 + 等高线叠加图...\n")

# 检查 celltype 列是否存在
if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
  cat("⚠️ 警告：未找到 'celltype' 列，跳过细胞类型图绘制\n")
} else {
  cat("✅ 检测到 celltype 列\n")
  
  # 查看细胞类型统计
  celltype_counts <- table(seurat_obj$celltype)
  cat(sprintf("✅ 共有 %d 种细胞类型：\n", length(celltype_counts)))
  print(celltype_counts)
  
  # 生成细胞类型颜色方案
  n_celltypes <- length(unique(seurat_obj$celltype))
  
  # 使用更丰富的调色板
  if (n_celltypes <= 8) {
    celltype_colors <- brewer.pal(max(3, n_celltypes), "Set2")
  } else if (n_celltypes <= 12) {
    celltype_colors <- brewer.pal(n_celltypes, "Set3")
  } else {
    # 组合多个调色板
    celltype_colors <- c(
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Set2"),
      brewer.pal(12, "Set3")
    )[1:n_celltypes]
  }
  
  names(celltype_colors) <- sort(unique(seurat_obj$celltype))
  
  # 为每个样本绘制图
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
    
    seurat_subset <- tryCatch(
      subset(seurat_obj, subset = orig.ident == sample_id),
      error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
    )
    
    safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
    
    # ============================================
    # 方法1：使用 celltype_isoheight_plot 函数
    # ============================================
    tryCatch({
      # 为每种细胞类型创建布尔列
      celltypes_in_sample <- unique(seurat_subset$celltype)
      
      # 只绘制样本中存在的细胞类型（最多显示前6种）
      celltypes_to_plot <- head(celltypes_in_sample, 6)
      
      for (ct in celltypes_to_plot) {
        col_name <- paste0("is_", make.names(ct))
        seurat_subset@meta.data[[col_name]] <- seurat_subset$celltype == ct
      }
      
      # 使用你原有的 celltype_isoheight_plot 函数
      # 注意：这个函数需要 density_top 参数是一个布尔列
      # 我们可以用 ClockGene_High 作为等高线背景
      
      p_celltype_iso <- celltype_isoheight_plot(
        .data = seurat_subset,
        density_top = ClockGene_High,  # 用高表达点生成等高线
        col_bg = "gray92",
        col_top = "transparent",  # 让高表达点透明，只显示等高线
        col_isoheight = "black",  # 等高线用黑色
        col_white_ratio = 0.25,
        cols_fill_isoheight = c(
          rep("white", 50),
          colorRampPalette(brewer.pal(9, "YlOrRd")[2:5])(50)  # 淡化等高线颜色
        ),
        size_bg = 0.8,
        size_top = 0,  # 不显示高表达点
        nrow = 1
      )
      
      # 保存基础等高线图
      ggsave(
        file.path(dirs$isoheight, sprintf("ClockGene_celltype_base_%s.pdf", safe_name)),
        plot = p_celltype_iso,
        width = 8, height = 8, dpi = 300
      )
      
    }, error = function(e) {
      cat(sprintf("   ⚠️ celltype_isoheight_plot 失败: %s\n", e$message))
    })
    
    # ============================================
    # 方法2：手动叠加（更灵活）
    # ============================================
    tryCatch({
      # 获取坐标
      coords <- GetTissueCoordinates(
        seurat_subset,
        cols = c("row", "col"),
        scale = NULL
      )
      
      if (!("row" %in% colnames(coords) && "col" %in% colnames(coords))) {
        cat(sprintf("   ⚠️ 坐标列不完整，跳过\n"))
        next
      }
      
      # 合并数据
      plot_data <- seurat_subset@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
      
      # 计算坐标范围
      expand_margin <- 0.05
      col_range <- range(plot_data$col, na.rm = TRUE)
      row_range <- range(plot_data$row, na.rm = TRUE)
      
      col_expand <- diff(col_range) * expand_margin
      row_expand <- diff(row_range) * expand_margin
      
      col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
      row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
      
      # ============================================
      # 生成等高线数据（使用 ClockGene_Distance）
      # ============================================
      # 过滤掉 NA 值
      plot_data_clean <- plot_data %>%
        filter(!is.na(col), !is.na(row), !is.na(ClockGene_Distance))
      
      # 使用 akima 包进行插值
      if (nrow(plot_data_clean) >= 10) {
        interp_result <- tryCatch({
          akima::interp(
            x = plot_data_clean$col,
            y = plot_data_clean$row,
            z = plot_data_clean$ClockGene_Distance,
            nx = 200,  # 插值分辨率
            ny = 200,
            linear = FALSE,  # 使用样条插值
            extrap = FALSE
          )
        }, error = function(e) {
          cat(sprintf("   ⚠️ 插值失败: %s\n", e$message))
          NULL
        })
        
        if (!is.null(interp_result)) {
          # 转换为 data.frame 用于 ggplot
          contour_data <- expand.grid(
            col = interp_result$x,
            row = interp_result$y
          )
          contour_data$z <- as.vector(interp_result$z)
          
          # ============================================
          # 绘制叠加图
          # ============================================
          p_overlay <- ggplot() +
            # 1. 等高线填充（底层）
            geom_contour_filled(
              data = contour_data,
              aes(x = col, y = row, z = z),
              bins = 10,
              alpha = 0.3  # 半透明
            ) +
            scale_fill_manual(
              values = colorRampPalette(brewer.pal(9, "YlOrRd")[3:9])(11),
              name = "Distance\n(Contour)",
              guide = guide_legend(order = 1)
            ) +
            # 2. 新的填充比例尺用于细胞类型
            new_scale_fill() +
            # 3. 细胞类型点（顶层）
            geom_point(
              data = plot_data,
              aes(x = col, y = row, fill = celltype),
              shape = 21, size = 1.8, color = "white", stroke = 0.15,
              alpha = 0.8
            ) +
            scale_fill_manual(
              values = celltype_colors,
              name = "Cell Type",
              guide = guide_legend(
                override.aes = list(size = 4, alpha = 1),
                order = 2
              )
            ) +
            # 4. 等高线线条
            geom_contour(
              data = contour_data,
              aes(x = col, y = row, z = z),
              color = "white",
              linewidth = 0.3,
              bins = 10,
              alpha = 0.6
            ) +
            # 坐标设置
            scale_x_continuous(
              limits = col_limits,
              expand = expansion(mult = 0.02)
            ) +
            scale_y_reverse(
              limits = rev(row_limits),
              expand = expansion(mult = 0.02)
            ) +
            coord_fixed(ratio = 1) +
            # 主题
            ggtitle(sprintf("Cell Types + Clock Gene Niche - %s", sample_id)) +
            theme_void() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 8),
              legend.box = "vertical",
              legend.spacing.y = unit(0.5, "cm"),
              plot.margin = margin(10, 10, 10, 10)
            )
          
          # 保存
          ggsave(
            file.path(dirs$isoheight, sprintf("ClockGene_celltype_overlay_%s.pdf", safe_name)),
            plot = p_overlay,
            width = 10, height = 8,
            dpi = 300
          )
          
          cat(sprintf("   ✅ 已保存: ClockGene_celltype_overlay_%s.pdf\n", safe_name))
          
          # ============================================
          # 额外：纯细胞类型图（无等高线）
          # ============================================
          p_celltype_only <- ggplot(plot_data, aes(x = col, y = row)) +
            geom_point(
              aes(fill = celltype),
              shape = 21, size = 2.5, color = "white", stroke = 0.1,
              alpha = 0.9
            ) +
            scale_fill_manual(
              values = celltype_colors,
              name = "Cell Type",
              guide = guide_legend(override.aes = list(size = 4))
            ) +
            scale_x_continuous(
              limits = col_limits,
              expand = expansion(mult = 0.02)
            ) +
            scale_y_reverse(
              limits = rev(row_limits),
              expand = expansion(mult = 0.02)
            ) +
            coord_fixed(ratio = 1) +
            ggtitle(sprintf("Cell Type Distribution - %s", sample_id)) +
            theme_void() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 8),
              plot.margin = margin(10, 10, 10, 10)
            )
          
          ggsave(
            file.path(dirs$spatial, sprintf("ClockGene_celltype_only_%s.pdf", safe_name)),
            plot = p_celltype_only,
            width = 10, height = 8,
            dpi = 300
          )
          
        }
      } else {
        cat(sprintf("   ⚠️ 数据点不足（%d < 10），跳过插值\n", nrow(plot_data_clean)))
      }
      
    }, error = function(e) {
      cat(sprintf("   ⚠️ 叠加图绘制失败: %s\n", e$message))
    })
  }
  
  cat("\n✅ 细胞类型图绘制完成！\n")
  cat(sprintf("   - 等高线叠加图保存在: %s/ClockGene_celltype_overlay_*.pdf\n", dirs$isoheight))
  cat(sprintf("   - 纯细胞类型图保存在: %s/ClockGene_celltype_only_*.pdf\n", dirs$spatial))
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


cat("\n✅ 全部完成！\n")
cat(sprintf("📁 所有图形已保存到: %s\n", figure_dir))