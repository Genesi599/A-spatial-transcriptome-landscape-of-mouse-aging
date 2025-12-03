# ==================================================================
# 02_density_zones.R (优化版)
# 密度分区计算（真正支持边界扩展）
# Author: Assistant
# Date: 2025-11-07
# ==================================================================

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


# ==================================================================
# 02_density_zones.R (兼容 Seurat + 可静默)
# 密度分区计算（真正支持边界扩展）
# ==================================================================

calculate_density_zones <- function(
  df,
  density_bins = 10,
  expand_margin = 0.1,
  # 手动列名映射（若与你的数据不同可以在调用处修改）
  col_col = "col",
  row_col = "row",
  flag_col = "ClockGene_High",
  quiet = FALSE
) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("需要安装 MASS 包")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("需要安装 dplyr 包")

  # 小工具：选择坐标列
  pick_coord_cols <- function(md, coords = NULL) {
    # 候选命名（从常见到少见）
    candidates <- list(
      list(col = "col",                  row = "row"),
      list(col = "x",                    row = "y"),
      list(col = "imagecol",             row = "imagerow"),
      list(col = "pxl_col_in_fullres",   row = "pxl_row_in_fullres"),
      list(col = "array_col",            row = "array_row")
    )
    # 若有 Seurat::GetTissueCoordinates 返回的 coords
    if (!is.null(coords)) {
      if (all(c("x","y") %in% colnames(coords))) {
        return(list(col = "x", row = "y", from = "coords"))
      }
      if (all(c("imagecol","imagerow") %in% colnames(coords))) {
        return(list(col = "imagecol", row = "imagerow", from = "coords"))
      }
    }
    # 在 meta.data 里找
    for (cand in candidates) {
      if (all(c(cand$col, cand$row) %in% colnames(md))) {
        return(list(col = cand$col, row = cand$row, from = "meta.data"))
      }
    }
    return(NULL)
  }

  # 若传入 Seurat，则从 meta.data/coords 取列；否则按 data.frame 处理
  if (inherits(df, "Seurat")) {
    seu <- df
    md <- seu@meta.data
    coords <- NULL
    if (requireNamespace("Seurat", quietly = TRUE)) {
      coords <- tryCatch(Seurat::GetTissueCoordinates(seu),
                         error = function(e) NULL)
      # 让 coords 行名与细胞一致，以便后续 join
      if (!is.null(coords) && is.null(rownames(coords)) &&
          "barcode" %in% colnames(coords)) {
        rownames(coords) <- coords$barcode
      }
      # 截齐到共同细胞
      if (!is.null(coords)) {
        common <- intersect(rownames(md), rownames(coords))
        md <- md[common, , drop = FALSE]
        coords <- coords[common, , drop = FALSE]
      }
    }
    # 优先使用用户显式指定的列名；否则自动探测
    coord_info <- NULL
    if (all(c(col_col, row_col) %in% colnames(md))) {
      coord_info <- list(col = col_col, row = row_col, from = "meta.data")
    } else if (!is.null(coords)) {
      coord_info <- pick_coord_cols(md = md, coords = coords)
      if (is.null(coord_info)) coord_info <- pick_coord_cols(md = md, coords = NULL)
    } else {
      coord_info <- pick_coord_cols(md = md, coords = NULL)
    }
    if (is.null(coord_info)) {
      stop("Seurat@meta.data 或坐标中找不到 col/row（常见别名也未找到），",
           "请在调用时显式指定 col_col/row_col 对应的列名")
    }

    # 组装 df：坐标来自 coords 或 meta.data
    if (coord_info$from == "coords") {
      tmp <- coords[, c(coord_info$col, coord_info$row), drop = FALSE]
    } else {
      tmp <- md[, c(coord_info$col, coord_info$row), drop = FALSE]
    }
    colnames(tmp) <- c("col", "row")

    if (!(flag_col %in% colnames(md))) {
      stop(sprintf("Seurat@meta.data 缺少标记列: %s", flag_col))
    }
    tmp[[flag_col]] <- md[[flag_col]]
    df <- tmp
  } else if (is.data.frame(df)) {
    # data.frame：重命名到标准列名
    req <- c(col_col, row_col, flag_col)
    if (!all(req %in% colnames(df))) {
      stop(sprintf("df 缺少列: %s",
                   paste(setdiff(req, colnames(df)), collapse = ", ")))
    }
    df <- df[, req, drop = TRUE]
    colnames(df) <- c("col", "row", "ClockGene_High")
  } else {
    stop("df 必须是 Seurat 或 data.frame")
  }

  # 标准化类型并去 NA
  df$col <- as.numeric(df$col)
  df$row <- as.numeric(df$row)
  if (!("ClockGene_High" %in% colnames(df))) {
    # 如果还没改名，补一手（Seurat分支已经保证存在）
    df$ClockGene_High <- df[[flag_col]]
  }
  if (is.numeric(df$ClockGene_High))
    df$ClockGene_High <- df$ClockGene_High > 0
  if (!is.logical(df$ClockGene_High))
    stop("ClockGene_High 必须是逻辑或数值(0/1)")
  df <- df[!is.na(df$col) & !is.na(df$row) & !is.na(df$ClockGene_High), ,
           drop = FALSE]

  # 1) 取高表达点
  high_points <- df %>% dplyr::filter(ClockGene_High)
  if (nrow(high_points) < 10) {
    warning("高表达点数量不足（< 10），无法计算密度")
    return(NULL)
  }

  # 2) 范围与扩展
  col_range <- range(df$col, na.rm = TRUE)
  row_range <- range(df$row, na.rm = TRUE)
  if (diff(col_range) == 0) col_range <- col_range + c(-0.5, 0.5)
  if (diff(row_range) == 0) row_range <- row_range + c(-0.5, 0.5)
  col_margin <- diff(col_range) * expand_margin
  row_margin <- diff(row_range) * expand_margin
  col_range_expanded <- c(col_range[1] - col_margin, col_range[2] + col_margin)
  row_range_expanded <- c(row_range[1] - row_margin, row_range[2] + row_margin)

  if (!quiet) {
    cat("   ✅ 密度计算范围:\n")
    cat(sprintf("      原始: col [%.1f, %.1f], row [%.1f, %.1f]\n",
                col_range[1], col_range[2], row_range[1], row_range[2]))
    cat(sprintf(paste0("      扩展: col [%.1f, %.1f], row [%.1f, %.1f] ",
                       "(边距=%.0f%%)\n"),
                col_range_expanded[1], col_range_expanded[2],
                row_range_expanded[1], row_range_expanded[2],
                expand_margin * 100))
  }

  # 3) KDE
  kde_result <- tryCatch({
    MASS::kde2d(
      x = high_points$col,
      y = high_points$row,
      n = 200,
      lims = c(col_range_expanded, row_range_expanded)
    )
  }, error = function(e) {
    warning(sprintf("密度估计失败: %s", e$message))
    return(NULL)
  })
  if (is.null(kde_result)) return(NULL)

  # 4) 网格
  density_df <- expand.grid(
    col = kde_result$x,
    row = kde_result$y
  )
  density_df$density <- as.vector(kde_result$z)
  max_density <- max(density_df$density, na.rm = TRUE)
  density_df$density_norm <- if (max_density > 0) density_df$density / max_density else 0
  if (!quiet) {
    cat(sprintf("   ✅ 密度网格包含扩展区域: %d x %d = %d 个点\n",
                length(kde_result$x), length(kde_result$y), nrow(density_df)))
  }

  # 5) 等距切分
  equal_breaks <- seq(0, 1, length.out = density_bins + 1)
  if (!quiet) {
    cat(sprintf("   ✅ Zone边界（等距切分，%d个区域）:\n", density_bins))
    for (i in seq_len(length(equal_breaks) - 1)) {
      cat(sprintf("      Zone_%d: [%.2f, %.2f)\n",
                  density_bins - i, equal_breaks[i], equal_breaks[i + 1]))
    }
  }
  density_df$density_zone <- cut(
    density_df$density_norm,
    breaks = equal_breaks,
    labels = sprintf("Zone_%d", (density_bins - 1):0),
    include.lowest = TRUE,
    right = TRUE
  )

  # 6) 分配到每个 spot
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

  # 最近邻填充 NA
  if (any(is.na(spot_zones$density_zone))) {
    na_spots <- which(is.na(spot_zones$density_zone))
    if (!quiet) cat(sprintf("   ⚠️  %d 个 spots 需要最近邻填充\n", length(na_spots)))
    for (idx in na_spots) {
      sc <- spot_zones$col[idx]; sr <- spot_zones$row[idx]
      d2 <- (density_df$col - sc)^2 + (density_df$row - sr)^2
      vi <- which(!is.na(density_df$density_zone))
      if (length(vi) > 0) {
        nv <- vi[which.min(d2[vi])]
        spot_zones$density_zone[idx] <- density_df$density_zone[nv]
        spot_zones$density_value[idx] <- density_df$density_norm[nv]
      }
    }
  }

  list(
    grid = density_df,
    spot_zones = spot_zones,
    kde_result = kde_result,
    equal_breaks = equal_breaks,
    col_range = col_range,
    row_range = row_range,
    col_range_expanded = col_range_expanded,
    row_range_expanded = row_range_expanded
  )
}
