# ===================================================================
# 08_plot_celltype_utils/01_color_schemes.R
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

