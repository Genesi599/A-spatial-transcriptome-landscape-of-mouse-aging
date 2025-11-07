# ===================================================================
# 02_density_zones.R (完整修复版)
# 密度分区计算（真正支持边界扩展）
# Author: Assistant (Fixed Version)
# Date: 2025-11-07
# ===================================================================

#' 计算密度分区（基于等距切分，支持边界扩展）
#'
#' @param df 数据框，必须包含 col, row, ClockGene_High 列
#' @param density_bins 等高线分级数量，默认 10（对应0.1间隔）
#' @param expand_margin 边界扩展比例，默认 0.1 (10%)
#'
#' @return 包含以下元素的列表：
#'   - grid: 密度网格数据框（包含扩展区域）
#'   - spot_zones: 每个spot的zone分配
#'   - kde_result: KDE计算结果
#'   - equal_breaks: 等距断点
#'   - col_range: 原始切片列坐标范围
#'   - row_range: 原始切片行坐标范围
#'   - col_range_expanded: 扩展后的列坐标范围
#'   - row_range_expanded: 扩展后的行坐标范围
#'
calculate_density_zones <- function(df, density_bins = 10, expand_margin = 0.1) {
  
  require(dplyr)
  require(MASS)
  
  # ========================================
  # 1. 提取高表达点
  # ========================================
  high_points <- df %>% dplyr::filter(ClockGene_High)
  
  if (nrow(high_points) < 10) {
    warning("高表达点数量不足（< 10），无法计算密度")
    return(NULL)
  }
  
  # ========================================
  # 2. 计算原始范围和扩展范围
  # ========================================
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  
  # ✅ 关键修复：计算扩展后的范围（用于网格生成）
  col_margin <- diff(col_range) * expand_margin
  row_margin <- diff(row_range) * expand_margin
  
  col_range_expanded <- c(col_range[1] - col_margin, col_range[2] + col_margin)
  row_range_expanded <- c(row_range[1] - row_margin, row_range[2] + row_margin)
  
  cat(sprintf("   ✅ 密度计算范围:\n"))
  cat(sprintf("      原始: col [%.1f, %.1f], row [%.1f, %.1f]\n",
              col_range[1], col_range[2], row_range[1], row_range[2]))
  cat(sprintf("      扩展: col [%.1f, %.1f], row [%.1f, %.1f] (边距=%.0f%%)\n",
              col_range_expanded[1], col_range_expanded[2], 
              row_range_expanded[1], row_range_expanded[2],
              expand_margin * 100))
  
  # ========================================
  # 3. 使用扩展范围进行 KDE 计算
  # ========================================
  kde_result <- tryCatch({
    MASS::kde2d(
      x = high_points$col,
      y = high_points$row,
      n = 200,  # 网格分辨率
      lims = c(col_range_expanded, row_range_expanded)  # ✅ 使用扩展范围
    )
  }, error = function(e) {
    warning(sprintf("密度估计失败: %s", e$message))
    return(NULL)
  })
  
  if (is.null(kde_result)) return(NULL)
  
  # ========================================
  # 4. 转换为 data frame（保留扩展区域）
  # ========================================
  density_df <- expand.grid(
    col = kde_result$x,
    row = kde_result$y
  )
  density_df$density <- as.vector(kde_result$z)
  
  # 归一化密度到 [0, 1]
  max_density <- max(density_df$density, na.rm = TRUE)
  if (max_density > 0) {
    density_df$density_norm <- density_df$density / max_density
  } else {
    density_df$density_norm <- 0
  }
  
  # ✅ 关键修复：不要过滤掉扩展区域！注释掉原来的过滤代码
  # density_df <- density_df %>%
  #   dplyr::filter(
  #     col >= col_range[1] & col <= col_range[2],
  #     row >= row_range[1] & row <= row_range[2]
  #   )
  
  cat(sprintf("   ✅ 密度网格包含扩展区域: %d x %d = %d 个点\n",
              length(kde_result$x), length(kde_result$y), nrow(density_df)))
  
  # ========================================
  # 5. 等距切分密度区域
  # ========================================
  equal_breaks <- seq(0, 1, length.out = density_bins + 1)
  
  # 打印边界信息
  cat(sprintf("   ✅ Zone边界（等距切分，%d个区域）:\n", density_bins))
  for (i in 1:(length(equal_breaks) - 1)) {
    cat(sprintf("      Zone_%d: [%.2f, %.2f)\n", 
                density_bins - i, 
                equal_breaks[i], 
                equal_breaks[i + 1]))
  }
  
  # 分级（Zone_0 = 最高密度核心，Zone_9 = 最低密度外围）
  density_df$density_zone <- cut(
    density_df$density_norm,
    breaks = equal_breaks,
    labels = sprintf("Zone_%d", (density_bins - 1):0),
    include.lowest = TRUE,
    right = TRUE
  )
  
  # ========================================
  # 6. 为每个 spot 分配密度区域
  # ========================================
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
  
  # 处理 NA（使用最近邻填充）
  if (any(is.na(spot_zones$density_zone))) {
    na_count <- sum(is.na(spot_zones$density_zone))
    cat(sprintf("   ⚠️  %d 个 spots 需要最近邻填充\n", na_count))
    
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
  
  # ========================================
  # 7. 返回结果
  # ========================================
  return(list(
    grid = density_df,                        # 包含扩展区域的完整网格
    spot_zones = spot_zones,
    kde_result = kde_result,
    equal_breaks = equal_breaks,
    col_range = col_range,                    # 原始范围
    row_range = row_range,
    col_range_expanded = col_range_expanded,  # 扩展范围
    row_range_expanded = row_range_expanded
  ))
}