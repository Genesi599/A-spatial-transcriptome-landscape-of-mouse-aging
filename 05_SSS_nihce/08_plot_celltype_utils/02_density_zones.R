# ===================================================================
# 02_density_zones.R
# 密度分区计算（0.1等距切分）
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 计算密度分区（基于等距切分）
#'
#' @param df 数据框，必须包含 col, row, ClockGene_High 列
#' @param density_bins 等高线分级数量，默认 10（对应0.1间隔）
#' @param expand_margin KDE计算时的扩展边距，默认 0.05
#'
#' @return 包含以下元素的列表：
#'   - grid: 密度网格数据框
#'   - spot_zones: 每个spot的zone分配
#'   - kde_result: KDE计算结果
#'   - equal_breaks: 等距断点
#'   - col_range: 切片列坐标范围
#'   - row_range: 切片行坐标范围
#'
#' @details
#' 使用MASS::kde2d计算密度，然后按0.1等距切分为N个区域：
#' - Zone_0: [0.9, 1.0] 核心高密度区
#' - Zone_1: [0.8, 0.9)
#' - ...
#' - Zone_9: [0.0, 0.1) 外围低密度区
#'
#' @examples
#' density_data <- calculate_density_zones(df, density_bins = 10)
#'
calculate_density_zones <- function(df, density_bins = 10, expand_margin = 0.05) {
  
  # 只使用高表达点计算密度
  high_points <- df %>% dplyr::filter(ClockGene_High)
  
  if (nrow(high_points) < 10) {
    warning("高表达点数量不足（< 10），无法计算密度")
    return(NULL)
  }
  
  # 计算坐标范围（不扩展，严格限制在切片范围内）
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  
  # 用于KDE计算的范围（稍微扩展以避免边界效应）
  col_expand <- diff(col_range) * expand_margin
  row_expand <- diff(row_range) * expand_margin
  
  col_limits_kde <- c(col_range[1] - col_expand, col_range[2] + col_expand)
  row_limits_kde <- c(row_range[1] - row_expand, row_range[2] + row_expand)
  
  # 使用 MASS::kde2d 计算密度
  kde_result <- tryCatch({
    MASS::kde2d(
      x = high_points$col,
      y = high_points$row,
      n = 200,
      lims = c(col_limits_kde, row_limits_kde)
    )
  }, error = function(e) {
    warning(sprintf("密度估计失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(kde_result)) return(NULL)
  
  # 转换为 data frame
  density_df <- expand.grid(
    col = kde_result$x,
    row = kde_result$y
  )
  density_df$density <- as.vector(kde_result$z)
  
  # 归一化密度
  density_df$density_norm <- density_df$density / max(density_df$density, na.rm = TRUE)
  
  # 只保留切片范围内的密度网格
  density_df <- density_df %>%
    dplyr::filter(
      col >= col_range[1] & col <= col_range[2],
      row >= row_range[1] & row <= row_range[2]
    )
  
  cat(sprintf("   ✅ 密度网格限制在切片范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range[1], col_range[2], row_range[1], row_range[2]))
  
  # 使用0.1等距切分
  equal_breaks <- seq(0, 1, length.out = density_bins + 1)
  
  # 打印边界信息
  cat(sprintf("   ✅ Zone边界（等距切分，%d个区域）:\n", density_bins))
  for (i in 1:(length(equal_breaks) - 1)) {
    cat(sprintf("      Zone_%d: [%.2f, %.2f)\n", 
                length(equal_breaks) - 1 - i, 
                equal_breaks[i], 
                equal_breaks[i + 1]))
  }
  
  # 分级（反转：Zone_0 = 最高密度核心）
  density_df$density_zone <- cut(
    density_df$density_norm,
    breaks = equal_breaks,
    labels = sprintf("Zone_%d", (length(equal_breaks) - 2):0),
    include.lowest = TRUE,
    right = TRUE
  )
  
  # 为每个spot分配最近的密度区域
  spot_zones <- df %>%
    dplyr::select(col, row) %>%
    dplyr::mutate(
      col_idx = sapply(col, function(x) which.min(abs(kde_result$x - x))),
      row_idx = sapply(row, function(y) which.min(abs(kde_result$y - y)))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      grid_col = kde_result$x[col_idx],
      grid_row = kde_result$y[row_idx]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      density_df %>% dplyr::select(col, row, density_zone, density_norm),
      by = c("grid_col" = "col", "grid_row" = "row")
    ) %>%
    dplyr::select(col, row, density_zone, density_value = density_norm)
  
  # 检查NA并填充（使用最近邻）
  if (any(is.na(spot_zones$density_zone))) {
    na_spots <- which(is.na(spot_zones$density_zone))
    
    for (idx in na_spots) {
      spot_col <- spot_zones$col[idx]
      spot_row <- spot_zones$row[idx]
      
      distances <- sqrt((density_df$col - spot_col)^2 + (density_df$row - spot_row)^2)
      valid_idx <- which(!is.na(density_df$density_zone))
      
      if (length(valid_idx) > 0) {
        nearest_valid <- valid_idx[which.min(distances[valid_idx])]
        spot_zones$density_zone[idx] <- density_df$density_zone[nearest_valid]
        spot_zones$density_value[idx] <- density_df$density_norm[nearest_valid]
      }
    }
  }
  
  return(list(
    grid = density_df,
    spot_zones = spot_zones,
    kde_result = kde_result,
    equal_breaks = equal_breaks,
    col_range = col_range,
    row_range = row_range
  ))
}