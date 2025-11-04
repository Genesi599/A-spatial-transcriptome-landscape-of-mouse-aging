GetAllCoordinates <- function(.data) {
  .data@images %>%
    names() %>%
    unique() %>%
    map_dfr(~ {
      GetTissueCoordinates(
        .data,
        image = .x,
        cols = c("row", "col"),
        scale = NULL
      ) %>%
        rownames_to_column(var = "cellid")
    })
}

single_marker <- function(df, intra_df, spot_type, dist_method, FUN, zero_check = FALSE) {
  if (nrow(intra_df) > 0) {  # ✅ 修复：使用 nrow() 而不是 length()
    all_df <- df %>%
      column_to_rownames("cellid") %>%
      select(row, col)

    # ✅ 添加调试信息
    cat(sprintf("  计算距离矩阵: %d 个查询点 × %d 个目标点\n", 
                nrow(all_df), nrow(intra_df)))

    mat <- proxy::dist(all_df, intra_df, method = dist_method) %>%
      as.matrix()

    spot_dist <- tibble(cellid = rownames(mat))
    
    # ✅ 修复：确保 else 分支可用
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      spot_dist[[spot_type]] <- matrixStats::rowMins(mat, na.rm = TRUE)
      cat("  使用 matrixStats::rowMins 计算最小距离\n")
    } else {
      # ✅ 取消注释 else 分支
      spot_dist[[spot_type]] <- apply(mat, 1, min, na.rm = TRUE)
      cat("  使用 apply(min) 计算最小距离\n")
    }

    # ✅ 添加验证
    cat(sprintf("  Distance 范围: %.2f ~ %.2f\n",
                min(spot_dist[[spot_type]], na.rm = TRUE),
                max(spot_dist[[spot_type]], na.rm = TRUE)))

    if (!is.na(FUN)) {
      spot_dist[[spot_type]] <- FUN(spot_dist[[spot_type]])
    }

    res <- df %>%
      left_join(spot_dist, by = "cellid")

  } else {
    # ✅ 修复：当没有标记点时，所有距离应该是 Inf
    cat("  ⚠️ 警告：没有找到标记点（intra_df 为空），Distance 设为 Inf\n")
    res <- df %>%
      mutate(!!spot_type := Inf)
  }

  res %>% select(-c(row, col))
}

# ------------------ 主函数（完全修复版）------------------ #

niche_marker <- function(
  .data,
  marker,
  spot_type,
  slide = "orig.ident",
  dist_method = "Euclidean",
  FUN = NA,
  n_work = 3
) {
  # 列名字符串
  marker <- as.character(substitute(marker))
  spot_type <- as.character(substitute(spot_type))
  slide <- as.character(substitute(slide))

  library(future)
  library(future.apply)

  plan(multisession, workers = n_work)
  options(future.globals.maxSize = Inf)
  message(">> 使用核心数: ", nbrOfWorkers())

  # ✅ 添加全局统计
  n_total <- ncol(.data)
  n_marker <- sum(.data@meta.data[[marker]], na.rm = TRUE)
  message(sprintf(">> 总点数: %d, 标记点数: %d (%.1f%%)",
                  n_total, n_marker, 100 * n_marker / n_total))

  # ✅✅✅ 关键修复：保存原始的细胞顺序
  original_cell_order <- colnames(.data)
  message(sprintf(">> 保存原始细胞顺序: %d 个细胞", length(original_cell_order)))

  .data@meta.data <-
    .data@meta.data %>%
    rownames_to_column(var = "cellid") %>%
    left_join(GetAllCoordinates(.data), by = "cellid") %>%
    group_by(.data[[slide]]) %>%
    group_split() %>%
    future_lapply(function(df) {
      slide_name <- df[[slide]][1]
      cat(sprintf("\n处理样本: %s\n", slide_name))

      # ✅ 修复：过滤时需要处理 NA 值
      intra_df <- df %>%
        filter(!is.na(.data[[marker]]) & .data[[marker]] == TRUE) %>%
        column_to_rownames("cellid") %>%
        select(row, col)

      cat(sprintf("  样本总点数: %d, 标记点数: %d\n", 
                  nrow(df), nrow(intra_df)))

      single_marker(df, intra_df, spot_type = spot_type,
                    dist_method = dist_method, FUN = FUN)
    }, future.chunk.size = Inf) %>%
    bind_rows() %>%
    column_to_rownames(var = "cellid")

  # ✅✅✅ 关键修复：严格按原始顺序重新排列
  message("\n>> 重新排序 meta.data 以匹配 Seurat object...")
  
  # 检查是否所有细胞都存在
  current_cells <- rownames(.data@meta.data)
  missing_cells <- setdiff(original_cell_order, current_cells)
  extra_cells <- setdiff(current_cells, original_cell_order)
  
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ 错误：meta.data 中缺少 %d 个细胞！", length(missing_cells)))
  }
  
  if (length(extra_cells) > 0) {
    warning(sprintf("⚠️ meta.data 中有 %d 个多余细胞，将被移除", length(extra_cells)))
    .data@meta.data <- .data@meta.data[original_cell_order, ]
  } else {
    # 强制按原始顺序重新排列
    .data@meta.data <- .data@meta.data[original_cell_order, ]
  }
  
  # 验证排序结果
  if (!identical(rownames(.data@meta.data), original_cell_order)) {
    stop("❌ 严重错误：重新排序后仍不匹配！")
  }
  
  message("✅ meta.data 行顺序已修正并验证")

  # ✅ 最终验证
  message("\n>> Distance 计算完成！")
  dist_vals <- .data@meta.data[[spot_type]]
  message(sprintf(">> Distance 统计: 最小=%.2f, 最大=%.2f, 均值=%.2f",
                  min(dist_vals, na.rm = TRUE),
                  max(dist_vals, na.rm = TRUE),
                  mean(dist_vals, na.rm = TRUE)))

  # ✅ 关键验证：标记点的 Distance 应该是 0
  marker_dist <- dist_vals[.data@meta.data[[marker]]]
  message(sprintf(">> 标记点的 Distance: 最小=%.2f, 最大=%.2f, 均值=%.2f",
                  min(marker_dist, na.rm = TRUE),
                  max(marker_dist, na.rm = TRUE),
                  mean(marker_dist, na.rm = TRUE)))
  
  # ✅✅✅ 增强验证：检查标记点 Distance = 0 的比例
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- length(marker_dist[!is.na(marker_dist)])
  pct_marker_zero <- 100 * n_marker_zero / n_marker_total
  
  message(sprintf(">> 标记点中 Distance=0 的比例: %d/%d (%.1f%%)",
                  n_marker_zero, n_marker_total, pct_marker_zero))
  
  if (pct_marker_zero < 99) {
    warning(sprintf("⚠️ 警告：只有 %.1f%% 的标记点 Distance=0！应该接近 100%%", pct_marker_zero))
    
    # 找出异常的标记点
    abnormal_marker_cells <- names(marker_dist)[marker_dist > 0]
    if (length(abnormal_marker_cells) > 0) {
      message(sprintf(">> 前 10 个异常标记点（Distance > 0）:"))
      abnormal_info <- .data@meta.data[head(abnormal_marker_cells, 10), 
                                       c(slide, marker, spot_type)]
      print(abnormal_info)
    }
  } else {
    message("✅ 验证通过：几乎所有标记点的 Distance = 0")
  }

  return(.data)
}