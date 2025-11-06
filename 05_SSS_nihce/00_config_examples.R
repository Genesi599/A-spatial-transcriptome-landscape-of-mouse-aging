#!/usr/bin/env Rscript
# ===================================================================
# 配置示例
# ===================================================================

# ===================================================================
# 示例 1: 批量处理目录中的所有 .rds 文件
# ===================================================================
CONFIG_EXAMPLE_1 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  recursive_search = FALSE,
  specific_files = NULL,
  exclude_files = NULL
)

# ===================================================================
# 示例 2: 批量处理指定的几个文件
# ===================================================================
CONFIG_EXAMPLE_2 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  specific_files = c("Lung_2-25M.rds", "Heart_2-25M.rds", "Brain_2-25M.rds"),
  exclude_files = NULL
)

# ===================================================================
# 示例 3: 批量处理，但排除某些文件
# ===================================================================
CONFIG_EXAMPLE_3 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  exclude_files = c("test.rds", "backup.rds", "old_version.rds")
)

# ===================================================================
# 示例 4: 递归搜索子目录
# ===================================================================
CONFIG_EXAMPLE_4 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  recursive_search = TRUE  # 会搜索所有子目录
)

# ===================================================================
# 示例 5: 单文件模式
# ===================================================================
CONFIG_EXAMPLE_5 <- list(
  batch_mode = FALSE,
  seurat_path = "/path/to/single_file.rds"
)

# ===================================================================
# 示例 6: 匹配特定命名模式的文件
# ===================================================================
CONFIG_EXAMPLE_6 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "^Lung.*\\.rds$",  # 只处理以 Lung 开头的文件
  recursive_search = FALSE
)

# ===================================================================
# 示例 7: 调试模式（只处理每个文件的前几个样本）
# ===================================================================
CONFIG_EXAMPLE_7 <- list(
  batch_mode = TRUE,
  seurat_dir = "/path/to/seurat_files",
  seurat_pattern = "\\.rds$",
  debug_mode = TRUE,
  debug_sample_limit = 2  # 每个文件只处理前2个样本
)