# ========================================================================
# niche_marker.R - 完整修复版
# 移除 GetTissueCoordinates 依赖，直接从 @images 提取坐标
# ========================================================================

library(dplyr)
library(tibble)
library(purrr)
library(proxy)
library(future)
library(future.apply)


# ========================================================================
# 1. 坐标提取函数（主方法）
# ========================================================================

GetAllCoordinates <- function(.data) {
  cat("\n🔍 提取空间坐标（直接从 @images 读取）...\n")
  
  image_names <- names(.data@images)
  
  if (length(image_names) == 0) {
    stop("❌ Seurat 对象中没有空间图像数据")
  }
  
  cat(sprintf(">> 发现 %d 个图像\n", length(image_names)))
  
  # 逐个提取坐标
  coords_list <- list()
  
  for (img_name in image_names) {
    cat(sprintf("  >> 提取: %s ... ", img_name))
    
    tryCatch({
      # 直接从 coordinates 槽提取
      img_obj <- .data@images[[img_name]]
      coords_df <- img_obj@coordinates
      
      # 转换为 data.frame
      if (!is.data.frame(coords_df)) {
        coords_df <- as.data.frame(coords_df)
      }
      
      # 获取细胞ID
      cell_ids <- rownames(coords_df)
      if (is.null(cell_ids) || all(is.na(cell_ids))) {
        stop("坐标数据缺少行名（细胞ID）")
      }
      
      # 查找坐标列
      col_names <- colnames(coords_df)
      
      # 可能的列名
      row_candidates <- c("row", "imagerow", "array_row", "tissue_row", "pxl_row_in_fullres")
      col_candidates <- c("col", "imagecol", "array_col", "tissue_col", "pxl_col_in_fullres")
      
      # 找到实际的列名
      row_col <- intersect(col_names, row_candidates)
      col_col <- intersect(col_names, col_candidates)
      
      if (length(row_col) == 0 || length(col_col) == 0) {
        stop(sprintf("未找到坐标列。可用列: %s", paste(col_names, collapse=", ")))
      }
      
      # 使用第一个匹配的列名
      row_col_name <- row_col[1]
      col_col_name <- col_col[1]
      
      # 提取坐标
      result <- data.frame(
        cellid = cell_ids,
        row = as.numeric(coords_df[[row_col_name]]),
        col = as.numeric(coords_df[[col_col_name]]),
        stringsAsFactors = FALSE
      )
      
      # 检查 NA
      n_na <- sum(is.na(result$row) | is.na(result$col))
      if (n_na > 0) {
        warning(sprintf("%s: %d 个细胞的坐标为 NA", img_name, n_na))
      }
      
      coords_list[[img_name]] <- result
      cat(sprintf("✓ %d 个细胞\n", nrow(result)))
      
    }, error = function(e) {
      cat(sprintf("❌ 失败: %s\n", e$message))
      stop(sprintf("样本 %s 坐标提取失败", img_name))
    })
  }
  
  # 合并所有坐标
  all_coords <- bind_rows(coords_list)
  
  if (nrow(all_coords) == 0) {
    stop("❌ 未能提取任何坐标数据")
  }
  
  cat(sprintf("✅ 总共提取 %d 个细胞的坐标\n\n", nrow(all_coords)))
  
  return(all_coords)
}


# ========================================================================
# 2. 单个样本的距离计算
# ========================================================================

single_marker <- function(df, intra_df, spot_type, dist_method, FUN) {
  
  if (nrow(intra_df) > 0) {
    # 准备所有细胞的坐标
    all_df <- df %>%
      tibble::column_to_rownames("cellid") %>%
      dplyr::select(row, col)

    cat(sprintf("  计算距离矩阵: %d 个查询点 × %d 个标记点\n", 
                nrow(all_df), nrow(intra_df)))

    # 计算距离矩阵
    mat <- proxy::dist(all_df, intra_df, method = dist_method) %>%
      as.matrix()

    # 计算每个细胞到最近标记点的距离
    spot_dist <- tibble(cellid = rownames(mat))
    
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      spot_dist[[spot_type]] <- matrixStats::rowMins(mat, na.rm = TRUE)
    } else {
      spot_dist[[spot_type]] <- apply(mat, 1, min, na.rm = TRUE)
    }

    # 应用转换函数（如果提供）
    if (!is.na(FUN)) {
      spot_dist[[spot_type]] <- FUN(spot_dist[[spot_type]])
    }

    # 合并回原始数据
    res <- df %>%
      dplyr::left_join(spot_dist, by = "cellid")

  } else {
    # 没有标记点，所有距离设为 Inf
    cat("  ⚠️ 警告：没有找到标记点，Distance 设为 Inf\n")
    res <- df %>%
      dplyr::mutate(!!spot_type := Inf)
  }

  # 移除坐标列
  res %>% dplyr::select(-c(row, col))
}


# ========================================================================
# 3. 主函数：Niche Marker 分析
# ========================================================================

niche_marker <- function(
  .data,
  marker,
  spot_type,
  slide = "orig.ident",
  dist_method = "Euclidean",
  FUN = NA,
  n_work = 3
) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   开始 Niche Marker 分析\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 转换参数为字符串
  marker <- as.character(substitute(marker))
  spot_type <- as.character(substitute(spot_type))
  slide <- as.character(substitute(slide))
  
  cat(sprintf("参数配置:\n"))
  cat(sprintf("  Marker 列: %s\n", marker))
  cat(sprintf("  输出列: %s\n", spot_type))
  cat(sprintf("  样本列: %s\n", slide))
  cat(sprintf("  距离方法: %s\n", dist_method))
  cat(sprintf("  工作线程: %d\n\n", n_work))

  # 加载必要的包
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("需要安装 future 包: install.packages('future')")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("需要安装 future.apply 包: install.packages('future.apply')")
  }

  library(future)
  library(future.apply)

  # 设置并行计算
  plan(multisession, workers = n_work)
  options(future.globals.maxSize = Inf)
  cat(sprintf(">> 并行核心数: %d\n\n", nbrOfWorkers()))

  # 数据统计
  n_total <- ncol(.data)
  n_marker <- sum(.data@meta.data[[marker]], na.rm = TRUE)
  cat(sprintf("数据概况:\n"))
  cat(sprintf("  总细胞数: %d\n", n_total))
  cat(sprintf("  标记细胞: %d (%.1f%%)\n\n", n_marker, 100 * n_marker / n_total))

  # 保存原始细胞顺序（关键！）
  original_cell_order <- colnames(.data)
  cat(sprintf(">> 保存原始细胞顺序: %d 个细胞\n\n", length(original_cell_order)))

  # ========== 提取空间坐标 ==========
  cat("🔄 提取空间坐标...\n")
  all_coords <- tryCatch({
    GetAllCoordinates(.data)
  }, error = function(e) {
    stop(sprintf("❌ 坐标提取失败: %s", e$message))
  })
  
  # 验证坐标完整性
  if (nrow(all_coords) != n_total) {
    stop(sprintf("❌ 坐标数量 (%d) 与细胞数量 (%d) 不匹配", 
                 nrow(all_coords), n_total))
  }
  
  # 检查是否所有细胞都有坐标
  missing_cells <- setdiff(original_cell_order, all_coords$cellid)
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ %d 个细胞缺少坐标数据", length(missing_cells)))
  }

  # ========== 合并 metadata 和坐标 ==========
  cat("\n🔄 合并 metadata 和坐标...\n")
  meta_with_coords <- .data@meta.data %>%
    tibble::rownames_to_column(var = "cellid") %>%
    dplyr::left_join(all_coords, by = "cellid")
  
  # 验证合并结果
  n_missing_coords <- sum(is.na(meta_with_coords$row) | is.na(meta_with_coords$col))
  if (n_missing_coords > 0) {
    stop(sprintf("❌ %d 个细胞在合并后缺少坐标", n_missing_coords))
  }
  cat("✅ 合并完成\n")

  # ========== 按样本分组并计算距离 ==========
  cat("\n🔄 按样本计算距离...\n")
  
  # 分组
  sample_groups <- meta_with_coords %>%
    dplyr::group_by(.data[[slide]]) %>%
    group_split()
  
  cat(sprintf(">> 将处理 %d 个样本\n\n", length(sample_groups)))

  # 并行处理每个样本
  results_list <- future_lapply(sample_groups, function(df) {
    
    slide_name <- df[[slide]][1]
    cat(sprintf("处理样本: %s\n", slide_name))
    
    # 提取标记点
    intra_df <- df %>%
      dplyr::filter(!is.na(.data[[marker]]) & .data[[marker]] == TRUE) %>%
      tibble::column_to_rownames("cellid") %>%
      dplyr::select(row, col)
    
    n_sample <- nrow(df)
    n_marker_sample <- nrow(intra_df)
    
    cat(sprintf("  样本细胞数: %d\n", n_sample))
    cat(sprintf("  标记细胞数: %d (%.1f%%)\n", 
                n_marker_sample, 100 * n_marker_sample / n_sample))
    
    # 计算距离
    result <- single_marker(
      df = df, 
      intra_df = intra_df, 
      spot_type = spot_type,
      dist_method = dist_method, 
      FUN = FUN
    )
    
    return(result)
    
  }, future.seed = TRUE, future.chunk.size = Inf)
  
  cat("\n🔄 合并所有样本结果...\n")

  # 合并结果
  combined_results <- bind_rows(results_list)
  
  # 将结果转换为以 cellid 为行名的 data.frame
  combined_results <- combined_results %>%
    tibble::column_to_rownames(var = "cellid")

  # ========== 恢复原始细胞顺序 ==========
  cat("\n🔄 恢复原始细胞顺序...\n")
  
  current_cells <- rownames(combined_results)
  missing_cells <- setdiff(original_cell_order, current_cells)
  extra_cells <- setdiff(current_cells, original_cell_order)
  
  if (length(missing_cells) > 0) {
    stop(sprintf("❌ 结果中缺少 %d 个细胞！", length(missing_cells)))
  }
  
  if (length(extra_cells) > 0) {
    warning(sprintf("⚠️ 结果中有 %d 个多余细胞，将被移除", length(extra_cells)))
    combined_results <- combined_results[original_cell_order, ]
  } else {
    # 按原始顺序重新排列
    combined_results <- combined_results[original_cell_order, ]
  }
  
  # 验证顺序
  if (!identical(rownames(combined_results), original_cell_order)) {
    stop("❌ 严重错误：细胞顺序恢复失败！")
  }
  
  cat("✅ 细胞顺序已恢复并验证\n")

  # ========== 将结果添加到 Seurat 对象 ==========
  .data@meta.data <- combined_results

  # ========== 结果验证 ==========
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("   结果验证\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  dist_vals <- .data@meta.data[[spot_type]]
  
  cat(sprintf("Distance 统计:\n"))
  cat(sprintf("  最小值: %.2f\n", min(dist_vals, na.rm = TRUE)))
  cat(sprintf("  最大值: %.2f\n", max(dist_vals, na.rm = TRUE)))
  cat(sprintf("  平均值: %.2f\n", mean(dist_vals, na.rm = TRUE)))
  cat(sprintf("  中位数: %.2f\n", median(dist_vals, na.rm = TRUE)))
  cat(sprintf("  NA 数量: %d\n\n", sum(is.na(dist_vals))))

  # 验证标记点的距离
  marker_cells <- .data@meta.data[[marker]]
  marker_dist <- dist_vals[marker_cells]
  
  n_marker_zero <- sum(marker_dist == 0, na.rm = TRUE)
  n_marker_total <- sum(!is.na(marker_dist))
  pct_zero <- 100 * n_marker_zero / n_marker_total
  
  cat(sprintf("标记细胞验证:\n"))
  cat(sprintf("  标记细胞总数: %d\n", n_marker_total))
  cat(sprintf("  Distance=0: %d (%.1f%%)\n", n_marker_zero, pct_zero))
  cat(sprintf("  Distance>0: %d (%.1f%%)\n", 
              n_marker_total - n_marker_zero, 
              100 - pct_zero))
  
  if (pct_zero < 95) {
    warning(sprintf("⚠️ 警告：只有 %.1f%% 的标记细胞 Distance=0！预期应接近 100%%", pct_zero))
    
    # 显示异常的标记细胞
    abnormal <- which(marker_dist > 0.1)
    if (length(abnormal) > 0) {
      cat(sprintf("\n前 5 个异常标记细胞:\n"))
      print(head(marker_dist[abnormal], 5))
    }
  } else {
    cat("\n✅ 验证通过：几乎所有标记细胞的 Distance = 0\n")
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("   Niche Marker 分析完成\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # 关闭并行
  plan(sequential)

  return(.data)
}


# ========================================================================
# 4. 辅助函数：诊断空间坐标
# ========================================================================

diagnose_spatial_coordinates <- function(.data) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("           空间坐标诊断报告\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  image_names <- names(.data@images)
  
  if (length(image_names) == 0) {
    cat("❌ 未找到空间图像数据\n\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("总图像数: %d\n\n", length(image_names)))
  
  for (i in seq_along(image_names)) {
    img_name <- image_names[i]
    cat(sprintf("[%d/%d] 图像: %s\n", i, length(image_names), img_name))
    cat("─────────────────────────────────────────────────────────\n")
    
    img_obj <- .data@images[[img_name]]
    coords <- img_obj@coordinates
    
    cat(sprintf("   细胞数: %d\n", nrow(coords)))
    cat(sprintf("   坐标列: %s\n", paste(colnames(coords), collapse=", ")))
    
    # 检查标准列
    has_row <- "row" %in% colnames(coords)
    has_col <- "col" %in% colnames(coords)
    
    cat(sprintf("   标准列: row=%s, col=%s\n", 
                ifelse(has_row, "✓", "✗"),
                ifelse(has_col, "✓", "✗")))
    
    # 如果有坐标，显示范围
    row_col <- intersect(colnames(coords), 
                        c("row", "imagerow", "array_row", "tissue_row"))
    col_col <- intersect(colnames(coords), 
                        c("col", "imagecol", "array_col", "tissue_col"))
    
    if (length(row_col) > 0 && length(col_col) > 0) {
      cat(sprintf("   坐标范围:\n"))
      cat(sprintf("      %s: [%.1f, %.1f]\n", 
                  row_col[1],
                  min(coords[[row_col[1]]], na.rm=TRUE), 
                  max(coords[[row_col[1]]], na.rm=TRUE)))
      cat(sprintf("      %s: [%.1f, %.1f]\n", 
                  col_col[1],
                  min(coords[[col_col[1]]], na.rm=TRUE), 
                  max(coords[[col_col[1]]], na.rm=TRUE)))
    }
    
    cat("\n")
  }
  
  cat("═══════════════════════════════════════════════════════════\n\n")
}
