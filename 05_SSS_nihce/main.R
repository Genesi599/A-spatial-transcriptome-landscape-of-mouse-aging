#!/usr/bin/env Rscript
# ===================================================================
# Clock Gene Niche Analysis Workflow
# Author: Zhangbin (adapted by ChatGPT)
# Date: 2025-11-04
# ===================================================================

# -----------------------------
# 0. 环境设置
# -----------------------------
# 设置当前工作目录为脚本所在路径
setwd("/data/home/quj_lab/yanghang/A-spatial-transcriptome-landscape-of-mouse-aging/05_SSS_nihce")

output_dir <- "/dellstorage09/quj_lab/yanghang/out_file"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("✅ 工作目录:", getwd(), "\n")
cat("✅ 输出目录:", output_dir, "\n")

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
})

# -----------------------------
# 2. 加载函数脚本
# -----------------------------
# 请将这些R脚本放在当前目录 (或指定完整路径)
source("niche_marker.R")
source("SSS_isoheight_plot.R")
# source("niche_grade_entropy.R")  # 如需熵分析可启用

# -----------------------------
# 3. 读取基因列表
# -----------------------------
cat("\n📄 读取基因列表...\n")
gene_list_path <- "/dellstorage09/quj_lab/yanghang/spatial/ref/NET_gene_list_mouse.txt"

gene_list <- read.table(gene_list_path, header = TRUE, stringsAsFactors = FALSE)[[1]]
gene_list <- trimws(gene_list)
gene_list <- gene_list[gene_list != ""]
cat(sprintf("✅ 共读取 %d 个基因。\n", length(gene_list)))
print(head(gene_list))

# -----------------------------
# 4. 加载 Seurat 对象
# -----------------------------
cat("\n🧠 加载 Seurat 对象...\n")
seurat_path <- "/dellstorage01/quj_lab/zhangbin/published_project/mouse_spatial_transcriptome_2024/stereo_seq_data/seurat_rds/Lymph_2-25M.rds"
seurat_obj <- readRDS(seurat_path)

cat(sprintf("✅ Spots 数量: %d, 基因数量: %d\n", ncol(seurat_obj), nrow(seurat_obj)))


# --------------------------------------------------
# 🧬 验证 Seurat 对象基因名格式
# --------------------------------------------------

cat("\n🔍 检查 Seurat 对象中的基因名格式...\n")

# 1. 查看前 10 个基因名
gene_names_preview <- head(rownames(seurat_obj), 100)
cat("📄 前 10 个基因名示例:\n")
print(gene_names_preview)

# 2. 统计特征
has_ensembl <- any(grepl("^ENSG", rownames(seurat_obj)))
has_version <- any(grepl("\\.", rownames(seurat_obj)))
has_symbol_suffix <- any(grepl(".*[_|][A-Z]+$", rownames(seurat_obj)))

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
    cat("💡 提示：基因名大小写可能不一致，可尝试统一大写：\n")
    cat("   gene_list <- toupper(gene_list)\n")
    cat("   rownames(seurat_obj) <- toupper(rownames(seurat_obj))\n")
  }
}

if (length(genes_missing) > 0) {
  cat("⚠️ 以下部分基因不在数据中（将被忽略）:\n")
  print(utils::head(genes_missing, 15))
  if (length(genes_missing) > 15)
    cat(sprintf("... 其余 %d 个未显示\n", length(genes_missing) - 15))
} else {
  cat("🎉 所有基因均成功匹配！\n")
}

# -----------------------------
# 6. 计算基因集评分 (Module Score)
# -----------------------------
cat("\n🧮 计算 Clock Gene Module Score...\n")

# Step 1️⃣ 计算 Score
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(clock_gene_set = genes_in_data),
  name = "ClockGene_Score"
)

# Step 2️⃣ 阈值计算（Top 30%）
threshold_value <- quantile(seurat_obj$ClockGene_Score1, 0.7, na.rm = TRUE)
cat(sprintf("✅ 高表达阈值设定为: %.3f (Top 30%%)\n", threshold_value))

# Step 3️⃣ 创建高低群组列
seurat_obj$ClockGene_niche <- ifelse(
  seurat_obj$ClockGene_Score1 > threshold_value,
  "ClockGene_High", "ClockGene_Low"
)
cat("✅ 已自动生成字段: ClockGene_niche\n")
print(table(seurat_obj$ClockGene_niche))

# Step 4️⃣ 在 meta.data 中添加布尔列（供 niche_marker 使用）
seurat_obj$Marker_Boolean <- seurat_obj$ClockGene_Score1 > threshold_value
cat("✅ 已在 meta.data 中创建布尔列: Marker_Boolean\n")

# Step 5️⃣ 并行设置
library(future)
plan(multisession, workers = 6)
cat("\n📈 开始 Niche 分析...\n")
cat(">> 可使用核心数: ", nbrOfWorkers(), "\n")

# Step 6️⃣ 调用修改后的 niche_marker 函数
seurat_obj <- niche_marker(
  .data = seurat_obj,
  marker = "Marker_Boolean",       # ⚠️ 注意：传入字符串列名
  spot_type = "ClockGene_niche",   # ⚠️ 也是列名字符串
  slide = "orig.ident",            # ⚠️ 同样列名字符串
  dist_method = "Euclidean",
  FUN = ceiling,
  n_work = 6
)

cat("✅ Niche 分析完成。\n")

# -----------------------------
# 8. 绘制 Isoheight 等高线图
# -----------------------------
cat("\n🎨 绘制 Isoheight 图...\n")

p_iso <- celltype_isoheight_plot(
  .data = seurat_obj,
  density_top = ClockGene_Score1 > threshold,
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

ggsave(
  file.path(output_dir, "ClockGene_niche_isoheight.pdf"),
  plot = p_iso,
  width = 14, height = 10, dpi = 300
)

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
  features = "ClockGene_niche",
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

cat("\n✅ 全部完成！结果保存在:\n")
cat(file.path(getwd(), output_dir), "\n")