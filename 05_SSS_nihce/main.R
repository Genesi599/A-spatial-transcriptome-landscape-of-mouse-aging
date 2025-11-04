#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis Workflow (Optimized Caching)
# Author: Zhangbin (adapted by ChatGPT)
# Date: 2025-11-04
# 优化点：只缓存计算结果，不重复保存大对象
# ===================================================================

# -----------------------------
# 0. 环境设置
# -----------------------------
setwd("/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce")

# ✅ 主输出目录
output_dir <- "/dellstorage09/quj_lab/yanghang/spatial"
cache_dir <- file.path(output_dir, "cache")

# ✅ 创建子文件夹结构
output_subdirs <- list(
  metadata = file.path(output_dir, "metadata"),
  isoheight = file.path(output_dir, "isoheight_plots"),
  spatial = file.path(output_dir, "spatial_plots"),
  score = file.path(output_dir, "score_plots"),
  distance = file.path(output_dir, "distance_plots"),
  sss_niche = file.path(output_dir, "sss_niche_plots")
)

# 创建所有子目录
lapply(c(output_dir, cache_dir, output_subdirs), function(d) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
})

cat("✅ 工作目录:", getwd(), "\n")
cat("✅ 输出目录:", output_dir, "\n")
cat("✅ 缓存目录:", cache_dir, "\n")
cat("✅ 子文件夹:\n")
for (name in names(output_subdirs)) {
  cat(sprintf("   - %s: %s\n", name, basename(output_subdirs[[name]])))
}

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
    size_mb <- file.size(cache_file) / 1024^2
    cat(sprintf("💾 缓存已保存: %s (%.2f MB)\n", basename(cache_file), size_mb))
    if (description != "") cat(sprintf("   📝 %s\n", description))
  }, error = function(e) {
    warning(sprintf("⚠️ 缓存保存失败: %s\n", e$message))
  })
}

# 加载缓存
load_cache <- function(cache_file, description = "") {
  if (file.exists(cache_file)) {
    size_mb <- file.size(cache_file) / 1024^2
    cat(sprintf("📂 加载缓存: %s (%.2f MB)\n", basename(cache_file), size_mb))
    if (description != "") cat(sprintf("   📝 %s\n", description))
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
cat("\n📚 加载分析函数...\n")
source("niche_marker.R")
source("SSS_isoheight_plot.R")
cat("✅ 函数加载完成\n")

# -----------------------------
# 3. 读取基因列表 (缓存)
# -----------------------------
cat("\n📄 读取基因列表...\n")
gene_list_path <- "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt"
gene_list_cache <- file.path(cache_dir, "gene_list.rds")

if (is_cache_valid(gene_list_cache, gene_list_path)) {
  gene_list <- load_cache(gene_list_cache, "基因列表")
} else {
  gene_list <- read.table(gene_list_path, header = TRUE, stringsAsFactors = FALSE)[[1]]
  gene_list <- trimws(gene_list)
  gene_list <- gene_list[gene_list != ""]
  save_cache(gene_list, gene_list_cache, "基因列表")
}

cat(sprintf("✅ 共读取 %d 个基因\n", length(gene_list)))
cat("📋 前 10 个基因:\n")
print(head(gene_list, 10))

# -----------------------------
# 4. 加载 Seurat 对象 (只读取一次)
# -----------------------------
cat("\n🧠 加载 Seurat 对象...\n")
seurat_path <- "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lymph_2-25M.rds"

cat("🔄 从原始路径读取（整个流程只读一次）...\n")
seurat_obj <- readRDS(seurat_path)

# ✅ 添加：更新对象到当前 Seurat 版本
cat("🔧 更新 Seurat 对象版本...\n")
seurat_obj <- UpdateSeuratObject(seurat_obj)
cat("✅ 对象更新完成\n")

cat(sprintf("✅ Spots: %d, Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))

# 基本信息
cat(sprintf("   样本数: %d\n", length(unique(seurat_obj$orig.ident))))
cat("   样本列表:\n")
print(table(seurat_obj$orig.ident))

# -----------------------------
# 5. 检查基因
# -----------------------------
cat("\n🔍 检查基因匹配情况...\n")
genes_in_data <- intersect(gene_list, rownames(seurat_obj))
genes_missing <- setdiff(gene_list, rownames(seurat_obj))

cat(sprintf("✅ 匹配基因: %d / %d (%.1f%%)\n",
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
  cat("⚠️ 以下基因不在数据中（将被忽略）:\n")
  print(head(genes_missing, 10))
  if (length(genes_missing) > 10) {
    cat(sprintf("   ... 其余 %d 个未显示\n", length(genes_missing) - 10))
  }
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
  
# 保存缓存（只保存计算结果列）
  score_data <- data.frame(ClockGene_Score1 = seurat_obj$ClockGene_Score1)
  save_cache(score_data, score_cache_file, "Module Score 评分")
}

cat("✅ Module Score 计算完成\n")
cat(sprintf("   评分范围: %.3f ~ %.3f\n", 
            min(seurat_obj$ClockGene_Score1, na.rm = TRUE),
            max(seurat_obj$ClockGene_Score1, na.rm = TRUE)))

# ✅ 灵活的阈值设置
THRESHOLD_QUANTILE <- 0.90  # ← 修改这里：0.90 = Top 10%, 0.95 = Top 5%, 0.99 = Top 1%

# 计算阈值
threshold <- quantile(seurat_obj$ClockGene_Score1, THRESHOLD_QUANTILE, na.rm = TRUE)

# 自动生成描述
threshold_pct <- (1 - THRESHOLD_QUANTILE) * 100
if (threshold_pct < 1) {
  threshold_desc <- sprintf("Top %.1f%%", threshold_pct)
} else if (threshold_pct == round(threshold_pct)) {
  threshold_desc <- sprintf("Top %d%%", as.integer(threshold_pct))
} else {
  threshold_desc <- sprintf("Top %.1f%%", threshold_pct)
}

cat(sprintf("✅ 高表达阈值: %.3f (%s, quantile=%.2f)\n", 
            threshold, threshold_desc, THRESHOLD_QUANTILE))

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
  cat("   使用多线程加速...\n")
  
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
  
  # 保存缓存（只保存计算结果列）
  niche_data <- data.frame(ClockGene_Distance = seurat_obj$ClockGene_Distance)
  save_cache(niche_data, niche_cache_file, "Niche 距离数据")
}

cat("✅ Niche 分析完成\n")
cat(sprintf("   距离范围: %.2f ~ %.2f\n",
            min(seurat_obj$ClockGene_Distance, na.rm = TRUE),
            max(seurat_obj$ClockGene_Distance, na.rm = TRUE)))


# -----------------------------
# 7.5. 绘图配置
# -----------------------------
cat("\n🔧 配置绘图参数...\n")

# ✅ 1. 调试模式开关
DEBUG_MODE <- TRUE  # ← 改为 FALSE 绘制所有样本
DEBUG_SAMPLE_LIMIT <- 3  # 调试模式下只画前 N 个样本

# ✅ 2. 获取所有样本名称
samples <- unique(seurat_obj$orig.ident)
cat(sprintf("✅ 检测到 %d 个样本\n", length(samples)))

# ✅ 3. 打印样本列表（便于检查）
if (length(samples) <= 10) {
  cat("📋 样本列表:\n")
  print(samples)
} else {
  cat("📋 前 10 个样本:\n")
  print(head(samples, 10))
  cat(sprintf("   ... 其余 %d 个未显示\n", length(samples) - 10))
}

# ✅ 4. 根据调试模式决定处理哪些样本
if (DEBUG_MODE) {
  samples_to_plot <- head(samples, min(DEBUG_SAMPLE_LIMIT, length(samples)))
  cat(sprintf("\n🔧 调试模式已启用：只处理前 %d 个样本\n", length(samples_to_plot)))
  cat("📋 待处理样本:", paste(samples_to_plot, collapse = ", "), "\n")
  cat("💡 关闭调试模式: 设置 DEBUG_MODE <- FALSE\n")
} else {
  samples_to_plot <- samples
  cat(sprintf("\n🚀 生产模式：将处理全部 %d 个样本\n", length(samples_to_plot)))
}


# -----------------------------
# 8. 绘制 Isoheight 图 - 分样本保存
# -----------------------------
cat("\n🎨 绘制 Isoheight 图（分样本）...\n")


# 为每个样本单独绘图
for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("\n📊 [%d/%d] 正在处理: %s\n", i, length(samples_to_plot), sample_id))
  
  # 提取单个样本
  tryCatch({
    seurat_subset <- subset(seurat_obj, subset = orig.ident == sample_id)
  }, error = function(e) {
    cat("   ⚠️ subset 失败，使用索引方法\n")
    sample_cells <- colnames(seurat_obj)[seurat_obj$orig.ident == sample_id]
    seurat_subset <<- seurat_obj[, sample_cells]
  })
  
  cat(sprintf("   Spots 数: %d\n", ncol(seurat_subset)))
  
  # 绘制等高线图
  cat("   🔄 绘制等高线图...\n")
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
  
  # ✅ 保存到 isoheight_plots 子文件夹
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  output_file <- file.path(output_subdirs$isoheight, 
                           sprintf("ClockGene_isoheight_%s.pdf", safe_name))
  ggsave(output_file, plot = p_iso, width = 8, height = 8, dpi = 300)
  cat(sprintf("✅ 已保存: %s\n", basename(output_file)))
}

if (DEBUG_MODE) {
  cat(sprintf("\n⚠️ 调试模式：已完成 %d/%d 个样本\n", length(samples_to_plot), length(samples)))
  cat("💡 关闭调试模式: 设置 DEBUG_MODE <- FALSE\n")
} else {
  cat("\n✅ 所有样本的等高线图已保存\n")
}


# -----------------------------
# 9. 可视化 Niche 距离梯度 - 分样本保存
# -----------------------------
cat("\n🔥 绘制空间梯度图（分样本）...\n")

# 使用相同的调试设置
if (DEBUG_MODE) {
  samples_to_plot <- head(samples, DEBUG_SAMPLE_LIMIT)
  cat(sprintf("🔧 调试模式：只处理前 %d 个样本\n", length(samples_to_plot)))
} else {
  samples_to_plot <- samples
  cat("🚀 生产模式：处理所有样本\n")
}

# 为每个样本单独绘图
for (i in seq_along(samples_to_plot)) {
  sample_id <- samples_to_plot[i]
  cat(sprintf("\n📊 [%d/%d] 正在处理: %s\n", i, length(samples_to_plot), sample_id))
  
  # 提取单个样本
  tryCatch({
    seurat_subset <- subset(seurat_obj, subset = orig.ident == sample_id)
  }, error = function(e) {
    cat("   ⚠️ subset 失败，使用索引方法\n")
    sample_cells <- colnames(seurat_obj)[seurat_obj$orig.ident == sample_id]
    seurat_subset <<- seurat_obj[, sample_cells]
  })
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # 绘制 Score 图
  cat("   🔄 绘制 Score 图...\n")
  p_score <- SpatialFeaturePlot(
    seurat_subset,
    features = "ClockGene_Score1",
    pt.size.factor = 1.5,
    alpha = c(0.1, 1)
  ) + scale_fill_gradientn(
    colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
    name = "Clock Gene\nScore"
  ) + ggtitle(sample_id) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # ✅ 修复 Distance 图配色
  cat("   🔄 绘制 Distance 图...\n")
  p_niche <- SpatialFeaturePlot(
    seurat_subset,
    features = "ClockGene_Distance",
    pt.size.factor = 1.5,
    alpha = c(0.1, 1)
  ) + scale_fill_gradientn(
    # ✅ 修复：移除 rev()，让小值（近）= 红色，大值（远）= 蓝色
    colors = c("#67001f", "#b2182b", "#d6604d", "#f4a582",
               "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac"),
    name = "Distance to\nHigh Score Region",
    # ✅ 添加清晰的图例标签
    labels = function(x) sprintf("%.0f", x)
  ) + ggtitle(sample_id) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # 合并图（左右对比）
  p_combined <- (p_score | p_niche) +
    plot_annotation(
      title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
    )
  
  # 保存到各自的子文件夹
  ggsave(
    file.path(output_subdirs$spatial, sprintf("ClockGene_spatial_%s.pdf", safe_name)),
    plot = p_combined,
    width = 18, height = 9, dpi = 300
  )
  
  ggsave(
    file.path(output_subdirs$score, sprintf("ClockGene_score_%s.pdf", safe_name)),
    plot = p_score,
    width = 10, height = 9, dpi = 300
  )
  
  ggsave(
    file.path(output_subdirs$distance, sprintf("ClockGene_distance_%s.pdf", safe_name)),
    plot = p_niche,
    width = 10, height = 9, dpi = 300
  )
  
  cat(sprintf("✅ 已保存 3 个图到: spatial/score/distance 文件夹\n"))
}

if (DEBUG_MODE) {
  cat(sprintf("\n⚠️ 调试模式：已完成 %d/%d 个样本\n", length(samples_to_plot), length(samples)))
  cat("💡 关闭调试模式: 设置 DEBUG_MODE <- FALSE\n")
} else {
  cat("\n✅ 所有空间梯度图已保存\n")
}

# -----------------------------
# 9.5. SSS Niche 热图可视化
# -----------------------------
cat("\n🎨 绘制 SSS Niche 热图...\n")

# 根据调试模式决定绘制的样本
if (DEBUG_MODE) {
  samples_to_plot_sss <- head(samples, DEBUG_SAMPLE_LIMIT)
  cat(sprintf("🔧 调试模式：只绘制前 %d 个样本的 SSS 热图\n", length(samples_to_plot_sss)))
} else {
  samples_to_plot_sss <- samples
  cat(sprintf("🚀 生产模式：绘制所有 %d 个样本的 SSS 热图\n", length(samples_to_plot_sss)))
}

# 为每个样本单独绘图
for (i in seq_along(samples_to_plot_sss)) {
  sample_id <- samples_to_plot_sss[i]
  cat(sprintf("\n📊 [%d/%d] 正在绘制: %s\n", i, length(samples_to_plot_sss), sample_id))
  
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # 输出到 sss_niche_plots 子文件夹
  output_file <- file.path(output_subdirs$sss_niche, 
                           sprintf("ClockGene_SSS_niche_%s.pdf", safe_name))
  
  # 提取单个样本数据
  tryCatch({
    sample_meta <- seurat_obj@meta.data %>%
      filter(orig.ident == sample_id) %>%
      rownames_to_column("cellid")
    
    # 检查坐标信息
    if (!all(c("col", "row") %in% colnames(sample_meta))) {
      cat("   🔄 获取空间坐标...\n")
      sample_coords <- GetAllCoordinates(seurat_obj[, seurat_obj$orig.ident == sample_id])
      sample_meta <- sample_meta %>%
        left_join(sample_coords, by = "cellid")
    }
    
    # 检查必需列
    required_cols <- c("col", "row", "ClockGene_High", "ClockGene_Distance")
    missing_cols <- setdiff(required_cols, colnames(sample_meta))
    
    if (length(missing_cols) > 0) {
      cat(sprintf("   ⚠️ 警告：缺少列 %s，跳过该样本\n", paste(missing_cols, collapse = ", ")))
      next
    }
    
    # 数据统计
    n_high <- sum(sample_meta$ClockGene_High, na.rm = TRUE)
    n_low <- sum(!sample_meta$ClockGene_High, na.rm = TRUE)
    cat(sprintf("   📊 SSS: %d spots (%.1f%%) | Others: %d spots (%.1f%%)\n", 
                n_high, 100 * n_high / nrow(sample_meta),
                n_low, 100 * n_low / nrow(sample_meta)))
    
    # ✅ 绘制 SSS 热图（修复配色）
    cat("   🔄 绘制 SSS 热图...\n")
    p_sss_niche <- ggplot(sample_meta, aes(x = col, y = row)) +
      # 1. 背景热图（显示 niche 距离）
      geom_tile(
        aes(fill = ClockGene_Distance), 
        width = 1, 
        height = 1
      ) +
      # ✅ 修复配色：距离近（小值）= 红色，距离远（大值）= 蓝色
      scale_fill_gradientn(
        colours = c(
          "#67001f", "#b2182b", "#d6604d", "#f4a582",  # 红色系（近）
          "#fddbc7", "#f7f7f7",                        # 白色过渡
          "#d1e5f0", "#92c5de", "#4393c3", "#2166ac"   # 蓝色系（远）
        ),
        name = "Distance\n(to Niche)",
        na.value = "white",
        # ✅ 添加更清晰的图例
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = 1.5,
          barheight = 10
        )
      ) +
      
      # 2. 叠加背景点 (Others - 低表达)
      geom_point(
        data = sample_meta %>% filter(ClockGene_High == FALSE),
        aes(x = col, y = row),
        color = "gray70",
        size = 0.3,
        alpha = 0.5
      ) +
      
      # 3. 高亮点 (SSS - 高表达，距离应该为 0，显示为深红色)
      geom_point(
        data = sample_meta %>% filter(ClockGene_High == TRUE),
        aes(x = col, y = row),
        color = "black",
        size = 0.8,
        alpha = 0.8
      ) +
      
      # 4. 坐标和主题
      scale_y_reverse() +
      coord_fixed(ratio = 1) +
      labs(
        title = sample_id,
        subtitle = sprintf(
          "🔴 SSS (High): %d spots (%.1f%%) | ⚪ Others: %d spots (%.1f%%)",
          n_high, 100 * n_high / nrow(sample_meta),
          n_low, 100 * n_low / nrow(sample_meta)
        )
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    # 保存 PDF
    ggsave(
      output_file, 
      plot = p_sss_niche, 
      width = 10, 
      height = 10, 
      dpi = 300
    )
    
    cat(sprintf("   ✅ 已保存: %s\n", basename(output_file)))
    
    # 调试模式下保存 PNG 预览
    if (DEBUG_MODE) {
      output_png <- file.path(output_subdirs$sss_niche, 
                              sprintf("ClockGene_SSS_niche_%s.png", safe_name))
      ggsave(output_png, plot = p_sss_niche, width = 10, height = 10, dpi = 150)
      cat(sprintf("   ✅ 已保存预览: %s\n", basename(output_png)))
    }
    
  }, error = function(e) {
    cat(sprintf("   ❌ 绘制失败: %s\n", conditionMessage(e)))
    cat("   跳过该样本...\n")
  })
}

cat("\n✅ SSS Niche 热图绘制完成\n")
if (DEBUG_MODE) {
  cat(sprintf("⚠️ 调试模式：只生成了 %d 张图\n", length(samples_to_plot_sss)))
  cat("💡 关闭调试模式以生成所有样本的图\n")
}

# -----------------------------
# 10. 保存结果（优化版 - 不保存大对象）
# -----------------------------
cat("\n💾 保存结果...\n")

# 1️⃣ 保存元数据到 metadata 子文件夹
cat("📝 保存元数据表格...\n")
write.csv(
  seurat_obj@meta.data,
  file.path(output_subdirs$metadata, "Lymph_2-25M_clockgene_metadata.csv"),
  row.names = TRUE
)
cat("✅ 元数据已保存 (CSV格式)\n")

# 2️⃣ 可选：保存完整对象
save_full_object <- FALSE  # ← 改为 TRUE 时才保存完整对象

if (save_full_object) {
  cat("\n⚠️ 正在保存完整 Seurat 对象（较慢，文件较大）...\n")
  saveRDS(
    seurat_obj, 
    file.path(output_subdirs$metadata, "Lymph_2-25M_with_clockgene_niche.rds")
  )
  cat("✅ 完整对象已保存\n")
} else {
  cat("\n💡 提示：完整 Seurat 对象未保存（节省时间和空间）\n")
  cat("   所有计算结果已缓存，重新运行脚本可快速恢复\n")
  cat("   如需保存完整对象用于分享，设置 save_full_object <- TRUE\n")
}

# -----------------------------
# 11. 缓存管理信息
# -----------------------------
cat("\n📊 缓存统计:\n")
cache_files <- list.files(cache_dir, full.names = TRUE, pattern = "\\.rds$")
if (length(cache_files) > 0) {
  cache_info <- file.info(cache_files)
  cache_sizes <- cache_info$size / 1024^2  # MB
  total_size <- sum(cache_sizes)
  
  cat(sprintf("   📦 缓存文件数: %d\n", length(cache_files)))
  cat(sprintf("   💾 总大小: %.1f MB\n", total_size))
  cat(sprintf("   📁 位置: %s\n", cache_dir))
  
  # 显示最大的缓存文件
  if (length(cache_files) >= 3) {
    top_idx <- order(cache_sizes, decreasing = TRUE)[1:min(3, length(cache_files))]
    cat("\n   🔝 最大的缓存文件:\n")
    for (i in seq_along(top_idx)) {
      idx <- top_idx[i]
      cat(sprintf("      %d. %s (%.1f MB)\n", 
                  i, basename(cache_files[idx]), cache_sizes[idx]))
    }
  }
  
  cat("\n💡 缓存管理命令:\n")
  cat("   查看: list.files(cache_dir)\n")
  cat("   清除: unlink(file.path(cache_dir, '*.rds'))\n")
  cat("   清空: unlink(cache_dir, recursive = TRUE)\n")
}

# 输出文件统计
cat("\n📂 输出文件统计:\n")
output_files <- list.files(output_dir, pattern = "\\.(pdf|csv|rds)$", full.names = TRUE)
if (length(output_files) > 0) {
  output_info <- file.info(output_files)
  output_sizes <- output_info$size / 1024^2
  cat(sprintf("   📄 输出文件数: %d\n", length(output_files)))
  cat(sprintf("   💾 总大小: %.1f MB\n", sum(output_sizes)))
  
  # 按类型统计
  pdf_files <- grep("\\.pdf$", output_files, value = TRUE)
  csv_files <- grep("\\.csv$", output_files, value = TRUE)
  rds_files <- grep("\\.rds$", output_files, value = TRUE)
  
  cat(sprintf("   - PDF 图形: %d 个\n", length(pdf_files)))
  cat(sprintf("   - CSV 表格: %d 个\n", length(csv_files)))
  cat(sprintf("   - RDS 对象: %d 个\n", length(rds_files)))
}

# -----------------------------
# 12. 完成总结
# -----------------------------
cat("\n" , rep("=", 60), "\n", sep = "")
cat("✅ 全部完成！\n")
cat(rep("=", 60), "\n", sep = "")

cat("\n📍 结果位置:\n")
cat(sprintf("   输出目录: %s\n", output_dir))
cat(sprintf("   缓存目录: %s\n", cache_dir))

cat("\n📊 生成的文件:\n")
cat("   - 元数据: Lymph_2-25M_clockgene_metadata.csv\n")
cat(sprintf("   - 等高线图: ClockGene_isoheight_*.pdf (%d 个)\n", length(samples)))
cat(sprintf("   - 空间梯度图: ClockGene_spatial_*.pdf (%d 个)\n", length(samples)))
cat(sprintf("   - Score 图: ClockGene_score_*.pdf (%d 个)\n", length(samples)))
cat(sprintf("   - Distance 图: ClockGene_distance_*.pdf (%d 个)\n", length(samples)))

cat("\n🚀 性能提示:\n")
cat("   首次运行: 完整计算（约 30-40 分钟）\n")
cat("   后续运行: 缓存加速（约 1-2 分钟）\n")

cat("\n💡 下一步:\n")
cat("   1. 查看元数据: read.csv('Lymph_2-25M_clockgene_metadata.csv')\n")
cat("   2. 查看图形: 打开 output_dir 中的 PDF 文件\n")
cat("   3. 修改参数: 调整 threshold (Top 5%) 或其他参数后重新运行\n")

cat("\n")