GetAllCoordinates <- function(.data) {
  cat("🔍 提取所有样本的空间坐标...\n")
  
  image_names <- names(.data@images)
  
  if (length(image_names) == 0) {
    stop("❌ Seurat 对象中没有空间图像数据")
  }
  
  # 使用 map_dfr 合并所有样本的坐标
  all_coords <- purrr::map_dfr(image_names, function(img_name) {
    cat(sprintf("  >> 提取样本: %s\n", img_name))
    
    tryCatch({
      # 方法1：直接从 coordinates 槽提取
      coords <- .data@images[[img_name]]@coordinates
      
      # 确保有 row 和 col 列
      if (!all(c("row", "col") %in% colnames(coords))) {
        # 尝试使用其他列名
        possible_row <- c("row", "imagerow", "array_row", "tissue_row")
        possible_col <- c("col", "imagecol", "array_col", "tissue_col")
        
        actual_row <- intersect(colnames(coords), possible_row)[1]
        actual_col <- intersect(colnames(coords), possible_col)[1]
        
        if (is.na(actual_row) || is.na(actual_col)) {
          stop(sprintf("未找到有效坐标列。可用列: %s", 
                      paste(colnames(coords), collapse=", ")))
        }
        
        coords$row <- coords[[actual_row]]
        coords$col <- coords[[actual_col]]
      }
      
      # 提取需要的列
      result <- coords %>%
        as.data.frame() %>%
        select(row, col) %>%
        rownames_to_column(var = "cellid")
      
      cat(sprintf("     ✓ 提取 %d 个细胞\n", nrow(result)))
      
      return(result)
      
    }, error = function(e) {
      cat(sprintf("     ❌ 提取失败: %s\n", e$message))
      
      # 尝试方法2：使用 GetTissueCoordinates（如果可用）
      tryCatch({
        cat("     🔄 尝试使用 GetTissueCoordinates...\n")
        
        # 不同的参数组合
        result <- tryCatch({
          # Seurat v4 风格
          Seurat::GetTissueCoordinates(
            object = .data,
            image = img_name,
            cols = c("row", "col")
          ) %>%
            rownames_to_column(var = "cellid")
        }, error = function(e2) {
          # Seurat v5 风格
          Seurat::GetTissueCoordinates(
            object = .data,
            image = img_name
          ) %>%
            rownames_to_column(var = "cellid")
        })
        
        cat(sprintf("     ✓ 提取 %d 个细胞\n", nrow(result)))
        return(result)
        
      }, error = function(e2) {
        cat(sprintf("     ❌ GetTissueCoordinates 也失败: %s\n", e2$message))
        return(NULL)
      })
    })
  })
  
  if (is.null(all_coords) || nrow(all_coords) == 0) {
    stop("❌ 无法提取任何坐标数据")
  }
  
  cat(sprintf("✅ 总共提取 %d 个细胞的坐标\n", nrow(all_coords)))
  
  # 验证数据
  if (any(is.na(all_coords$row)) || any(is.na(all_coords$col))) {
    warning("⚠️ 坐标中包含 NA 值")
  }
  
  return(all_coords)
}

single_marker <- function(df, intra_df, spot_type, dist_method, FUN, zero_check = FALSE) {
  if (nrow(intra_df) > 0) {
    all_df <- df %>%
      column_to_rownames("cellid") %>%
      select(row, col)

    cat(sprintf("  计算距离矩阵: %d 个查询点 × %d 个目标点\n", 
                nrow(all_df), nrow(intra_df)))

    mat <- proxy::dist(all_df, intra_df, method = dist_method) %>%
      as.matrix()

    spot_dist <- tibble(cellid = rownames(mat))
    
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      spot_dist[[spot_type]] <- matrixStats::rowMins(mat, na.rm = TRUE)
    } else {
      spot_dist[[spot_type]] <- apply(mat, 1, min, na.rm = TRUE)
    }

    if (!is.na(FUN)) {
      spot_dist[[spot_type]] <- FUN(spot_dist[[spot_type]])
    }

    res <- df %>%
      left_join(spot_dist, by = "cellid")

  } else {
    cat("  ⚠️ 警告：没有找到标记点，Distance 设为 Inf\n")
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
  library(dplyr)
  library(tibble)

  plan(multisession, workers = n_work)
  options(future.globals.maxSize = Inf)
  message(">> 使用核心数: ", nbrOfWorkers())

  # 全局统计
  n_total <- ncol(.data)
  n_marker <- sum(.data@meta.data[[marker]], na.rm = TRUE)
  message(sprintf(">> 总点数: %d, 标记点数: %d (%.1f%%)",
                  n_total, n_marker, 100 * n_marker / n_total))

  # 保存原始细胞顺序
  original_cell_order <- colnames(.data)
  message(sprintf(">> 保存原始细胞顺序: %d 个细胞", length(original_cell_order)))

  # ========== 关键修复：使用修复版的 GetAllCoordinates ==========
  all_coords <- tryCatch({
    GetAllCoordinates(.data)
  }, error = function(e) {
    message("⚠️ GetAllCoordinates 失败，尝试简单版本...")
    GetAllCoordinates_Simple(.data)
  })
  
  # 验证坐标提取结果
  if (nrow(all_coords) != n_total) {
    warning(sprintf("⚠️ 坐标数量 (%d) 与细胞数量 (%d) 不匹配", 
                   nrow(all_coords), n_total))
  }

  .data@meta.data <-
    .data@meta.data %>%
    rownames_to_column(var = "cellid") %>%
    left_join(all_coords, by = "cellid") %>%
    group_by(.data[[slide]]) %>%
    group_split() %>%
    future_lapply(function(df) {
      slide_name <- df[[slide]][1]
      cat(sprintf("\n处理样本: %s\n", slide_name))

      # 过滤标记点
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

  # 严格按原始顺序重新排列
  message("\n>> 重新排序 meta.data 以匹配 Seurat object...")
  
  current_cells <- rownames(.data@meta.data)
  missing_cells <- setdiff(original_cell_order, current_cells)
  extra_cells <- setdiff(current_cells, original_cell_order)
  
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ 错误：meta.data 中缺少 %d 个细胞！", length(missing_cells)))
  }
  
  if (length(extra_cells) > 0) {
    warning(sprintf("⚠️ meta.data 中有 %d 个多余细胞，将被移除", length(extra_cells)))
  }
  
  .data@meta.data <- .data@meta.data[original_cell_order, ]
  
  if (!identical(rownames(.data@meta.data), original_cell_order)) {
    stop("❌ 严重错误：重新排序后仍不匹配！")
  }
  
  message("✅ meta.data 行顺序已修正并验证")

  # 最终验证
  message("\n>> Distance 计算完成！")
  dist_vals <- .data@meta.data[[spot_type]]
  message(sprintf(">> Distance 统计: 最小=%.2f, 最大=%.2f, 均值=%.2f",
                  min(dist_vals, na.rm = TRUE),
                  max(dist_vals, na.rm = TRUE),
                  mean(dist_vals, na.rm = TRUE)))

  # 验证标记点
  marker_dist <- dist_vals[.data@meta.data[[marker]]]
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- length(marker_dist[!is.na(marker_dist)])
  pct_marker_zero <- 100 * n_marker_zero / n_marker_total
  
  message(sprintf(">> 标记点中 Distance=0 的比例: %d/%d (%.1f%%)",
                  n_marker_zero, n_marker_total, pct_marker_zero))
  
  if (pct_marker_zero < 99) {
    warning(sprintf("⚠️ 警告：只有 %.1f%% 的标记点 Distance=0！", pct_marker_zero))
  }

  return(.data)
}
