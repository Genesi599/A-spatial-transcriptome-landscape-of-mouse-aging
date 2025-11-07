# ===================================================================
# 01_color_schemes.R
# 统一的颜色方案管理（全局配色版）
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

library(RColorBrewer)

# ===================================================================
# 原有函数（保持不变）
# ===================================================================

#' 生成统一的zone颜色方案（支持任意数量的区域）
#'
#' @param n_zones zone数量，默认10
#' @return 命名的颜色向量
#'
#' @details
#' Zone_0 = 核心（深红色，高密度）
#' Zone_N = 外围（深蓝色，低密度）
#' 使用深红到深蓝的渐变色系
#'
#' @examples
#' zone_colors <- get_zone_colors(10)
#' zone_colors["Zone_0"]  # 深红色
#'
get_zone_colors <- function(n_zones = 10) {
  # 从深红到深蓝的渐变
  zone_colors <- colorRampPalette(c(
    "#67001f",  # 深红（Zone_0，核心，高密度）
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061"   # 深蓝（Zone_N-1，外围，低密度）
  ))(n_zones)
  
  # Zone_0 对应第一个颜色（深红）
  zone_names <- sprintf("Zone_%d", 0:(n_zones - 1))
  names(zone_colors) <- zone_names
  
  return(zone_colors)
}


#' 生成统一的细胞类型颜色
#'
#' @param celltypes 细胞类型向量
#' @return 命名的颜色向量
#'
#' @details
#' 根据细胞类型数量自动选择合适的调色板：
#' - ≤8种：使用 RColorBrewer Set2
#' - ≤12种：使用 RColorBrewer Set3
#' - >12种：组合 Set1 + Set2 + Set3
#'
#' @examples
#' celltype_colors <- get_celltype_colors(c("T cell", "B cell", "Macrophage"))
#'
get_celltype_colors <- function(celltypes) {
  n_celltypes <- length(celltypes)
  
  if (n_celltypes <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n_celltypes), "Set2")
  } else if (n_celltypes <= 12) {
    colors <- RColorBrewer::brewer.pal(n_celltypes, "Set3")
  } else {
    colors <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3")
    )[1:n_celltypes]
  }
  
  names(colors) <- celltypes
  return(colors)
}


#' 获取等高线颜色渐变
#'
#' @param n_breaks 等高线断点数量
#' @return 颜色向量（从深蓝到深红）
#'
#' @examples
#' contour_colors <- get_contour_colors(11)
#'
get_contour_colors <- function(n_breaks) {
  colorRampPalette(c(
    "#053061",  # 深蓝 (低密度)
    "#2166ac",
    "#4393c3",
    "#92c5de",
    "#d1e5f0",
    "#fddbc7",
    "#f4a582",
    "#d6604d",
    "#b2182b",
    "#67001f"   # 深红 (高密度)
  ))(n_breaks)
}


# ===================================================================
# 新增：全局统一配色功能
# ===================================================================

#' 生成全局统一的颜色方案（核心函数）
#'
#' @param sample_list 所有样本列表
#' @param celltype_col 细胞类型列名
#' @param density_bins 密度分区数量
#' 
#' @return 包含所有颜色映射的列表
#'
#' @details
#' 这个函数会：
#' 1. 遍历所有样本，收集所有独特的细胞类型
#' 2. 为每个细胞类型分配固定颜色
#' 3. 为密度区域生成渐变颜色
#' 4. 返回统一的颜色方案供所有图表使用
#'
#' @examples
#' color_scheme <- create_global_color_scheme(sample_list, "celltype", 10)
#' CONFIG$colors <- color_scheme
#'
create_global_color_scheme <- function(sample_list, celltype_col, density_bins) {
  
  cat("\n🎨 生成全局统一颜色方案...\n")
  
  # ========================================
  # 1. 收集所有细胞类型
  # ========================================
  
  all_celltypes <- character(0)
  
  for (i in seq_along(sample_list)) {
    
    sample_id <- names(sample_list)[i]
    seurat_obj <- sample_list[[sample_id]]
    
    # 自动检测细胞类型列
    detected_col <- detect_celltype_column(seurat_obj, celltype_col)
    
    if (is.null(detected_col)) {
      warning(sprintf("样本 %s 无法找到细胞类型列，跳过", sample_id))
      next
    }
    
    # 提取并清理细胞类型名称
    celltypes <- seurat_obj@meta.data[[detected_col]]
    celltypes_clean <- clean_celltype_names(celltypes)
    
    all_celltypes <- c(all_celltypes, celltypes_clean)
  }
  
  # 去重并排序
  unique_celltypes <- sort(unique(all_celltypes))
  n_celltypes <- length(unique_celltypes)
  
  if (n_celltypes == 0) {
    stop("❌ 未找到任何细胞类型")
  }
  
  cat(sprintf("   📊 发现 %d 个独特细胞类型\n", n_celltypes))
  
  # ========================================
  # 2. 生成细胞类型颜色映射
  # ========================================
  
  celltype_colors <- get_celltype_colors(unique_celltypes)
  
  cat(sprintf("   ✅ 为 %d 个细胞类型分配了固定颜色\n", n_celltypes))
  
  # 打印颜色映射（前10个）
  if (n_celltypes <= 10) {
    cat("\n   细胞类型颜色映射:\n")
    for (ct in unique_celltypes) {
      cat(sprintf("      • %-25s → %s\n", ct, celltype_colors[ct]))
    }
  } else {
    cat("\n   细胞类型颜色映射（前10个）:\n")
    for (i in 1:10) {
      ct <- unique_celltypes[i]
      cat(sprintf("      • %-25s → %s\n", ct, celltype_colors[ct]))
    }
    cat(sprintf("      ... 还有 %d 个细胞类型\n", n_celltypes - 10))
  }
  
  # ========================================
  # 3. 生成密度区域颜色映射
  # ========================================
  
  zone_colors <- get_zone_colors(density_bins)
  
  cat(sprintf("\n   ✅ 为 %d 个密度区域分配了渐变颜色\n", density_bins))
  cat(sprintf("      Zone_0 (核心) → %s (深红)\n", zone_colors["Zone_0"]))
  cat(sprintf("      Zone_%d (外围) → %s (深蓝)\n", 
              density_bins - 1, zone_colors[sprintf("Zone_%d", density_bins - 1)]))
  
  # ========================================
  # 4. 返回颜色配置
  # ========================================
  
  color_scheme <- list(
    celltype = celltype_colors,
    density_zone = zone_colors,
    n_celltypes = n_celltypes,
    n_zones = density_bins,
    celltype_names = unique_celltypes,
    zone_names = names(zone_colors),
    celltype_col = celltype_col
  )
  
  cat("   ✅ 全局颜色方案已创建\n\n")
  
  return(color_scheme)
}


#' 自动检测细胞类型列
#'
#' @param seurat_obj Seurat 对象
#' @param preferred_col 首选列名
#' 
#' @return 检测到的列名或 NULL
#'
#' @details
#' 按以下顺序尝试：
#' 1. 首选列名（如 "celltype"）
#' 2. 常见的细胞类型列名
#' 3. 包含关键词的列（模糊匹配）
#'
detect_celltype_column <- function(seurat_obj, preferred_col = "celltype") {
  
  meta_cols <- colnames(seurat_obj@meta.data)
  
  # 1. 如果首选列存在，直接使用
  if (preferred_col %in% meta_cols) {
    return(preferred_col)
  }
  
  # 2. 尝试常见的细胞类型列名（按优先级排序）
  candidate_cols <- c(
    "celltype", "cell_type", "CellType",
    "predicted.id", "predicted.celltype",
    "annotation", "Annotation",
    "celltype.l1", "celltype.l2",
    "cluster", "seurat_clusters"
  )
  
  for (col in candidate_cols) {
    if (col %in% meta_cols) {
      return(col)
    }
  }
  
  # 3. 模糊匹配（包含关键词）
  pattern_matches <- grep("type|cluster|annotation|label|class", 
                          meta_cols, ignore.case = TRUE, value = TRUE)
  
  if (length(pattern_matches) > 0) {
    return(pattern_matches[1])
  }
  
  # 4. 未找到
  return(NULL)
}


#' 清理细胞类型名称
#'
#' @param celltypes 原始细胞类型向量
#' 
#' @return 清理后的细胞类型向量
#'
#' @details
#' - 移除特殊字符，替换为下划线
#' - 移除多余的下划线
#' - 移除首尾下划线
#'
clean_celltype_names <- function(celltypes) {
  
  celltypes_clean <- gsub("[^[:alnum:]_]", "_", celltypes)
  celltypes_clean <- gsub("_{2,}", "_", celltypes_clean)
  celltypes_clean <- gsub("^_|_$", "", celltypes_clean)
  
  return(celltypes_clean)
}


#' 从颜色方案中获取细胞类型颜色
#'
#' @param celltype 细胞类型名称（可以是向量）
#' @param color_scheme 全局颜色方案（来自 CONFIG$colors）
#' 
#' @return 颜色值（向量）
#'
#' @examples
#' color <- get_color_for_celltype("T_cells", CONFIG$colors)
#' colors <- get_color_for_celltype(c("T_cells", "B_cells"), CONFIG$colors)
#'
get_color_for_celltype <- function(celltype, color_scheme) {
  
  sapply(celltype, function(ct) {
    if (ct %in% names(color_scheme$celltype)) {
      return(color_scheme$celltype[ct])
    } else {
      warning(sprintf("细胞类型 '%s' 未找到，使用灰色", ct))
      return("#CCCCCC")
    }
  }, USE.NAMES = FALSE)
}


#' 从颜色方案中获取区域颜色
#'
#' @param zone 区域名称（可以是向量）
#' @param color_scheme 全局颜色方案（来自 CONFIG$colors）
#' 
#' @return 颜色值（向量）
#'
#' @examples
#' color <- get_color_for_zone("Zone_0", CONFIG$colors)
#' colors <- get_color_for_zone(c("Zone_0", "Zone_1"), CONFIG$colors)
#'
get_color_for_zone <- function(zone, color_scheme) {
  
  sapply(zone, function(z) {
    if (z %in% names(color_scheme$density_zone)) {
      return(color_scheme$density_zone[z])
    } else {
      warning(sprintf("区域 '%s' 未找到，使用灰色", z))
      return("#CCCCCC")
    }
  }, USE.NAMES = FALSE)
}


#' 向后兼容函数：setup_colors
#'
#' @param seurat_obj 单个 Seurat 对象（已弃用）
#' @param CONFIG 配置对象
#' @param celltype_col 细胞类型列名（已弃用）
#' @param density_bins 密度分区数量（已弃用）
#' 
#' @return 无返回值
#'
#' @details
#' 这个函数保持向后兼容，但实际上不再做任何事情。
#' 请在分析前使用 create_global_color_scheme() 生成全局颜色方案。
#'
setup_colors <- function(seurat_obj, CONFIG, celltype_col, density_bins) {
  
  # 如果 CONFIG 中已经有全局颜色方案，直接返回
  if (!is.null(CONFIG$colors$celltype) && !is.null(CONFIG$colors$density_zone)) {
    # 静默返回（颜色已设置）
    return(invisible(NULL))
  }
  
  # 否则警告
  warning("⚠️  setup_colors() 已弃用。", 
          "请在 analyze_celltype_niche() 开始时自动调用 create_global_color_scheme()")
  
  return(invisible(NULL))
}


# ===================================================================
# 工具函数：验证颜色方案
# ===================================================================

#' 验证颜色方案完整性
#'
#' @param color_scheme 颜色方案对象
#' 
#' @return 逻辑值，TRUE 表示有效
#'
validate_color_scheme <- function(color_scheme) {
  
  required_fields <- c("celltype", "density_zone", "n_celltypes", "n_zones")
  
  missing_fields <- setdiff(required_fields, names(color_scheme))
  
  if (length(missing_fields) > 0) {
    warning(sprintf("颜色方案缺少字段: %s", paste(missing_fields, collapse = ", ")))
    return(FALSE)
  }
  
  if (length(color_scheme$celltype) != color_scheme$n_celltypes) {
    warning("细胞类型颜色数量与 n_celltypes 不匹配")
    return(FALSE)
  }
  
  if (length(color_scheme$density_zone) != color_scheme$n_zones) {
    warning("区域颜色数量与 n_zones 不匹配")
    return(FALSE)
  }
  
  return(TRUE)
}


#' 打印颜色方案摘要
#'
#' @param color_scheme 颜色方案对象
#' 
#' @return 无返回值
#'
print_color_scheme <- function(color_scheme) {
  
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("   颜色方案摘要\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  cat(sprintf("细胞类型数量: %d\n", color_scheme$n_celltypes))
  cat(sprintf("密度区域数量: %d\n", color_scheme$n_zones))
  
  if (color_scheme$n_celltypes <= 15) {
    cat("\n细胞类型颜色:\n")
    for (ct in names(color_scheme$celltype)) {
      cat(sprintf("  %-30s %s\n", ct, color_scheme$celltype[ct]))
    }
  } else {
    cat("\n细胞类型颜色（前15个）:\n")
    for (i in 1:15) {
      ct <- names(color_scheme$celltype)[i]
      cat(sprintf("  %-30s %s\n", ct, color_scheme$celltype[ct]))
    }
    cat(sprintf("  ... 还有 %d 个\n", color_scheme$n_celltypes - 15))
  }
  
  cat("\n密度区域颜色:\n")
  for (zn in names(color_scheme$density_zone)) {
    cat(sprintf("  %-10s %s\n", zn, color_scheme$density_zone[zn]))
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
}


cat("✅ 01_color_schemes.R 已加载（全局配色版）\n")