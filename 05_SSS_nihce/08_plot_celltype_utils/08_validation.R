#!/usr/bin/env Rscript
# ===================================================================
# 验证模块
# ===================================================================

#' 验证输入参数
#' 
#' @param sample_list 样本列表
#' @param CONFIG 配置对象
validate_inputs <- function(sample_list, CONFIG) {
  
  if (!is.list(sample_list) || length(sample_list) == 0) {
    stop("❌ sample_list 必须是非空列表")
  }
  
  # 验证必需目录
  required_dirs <- c("overlay", "celltype", "composition", "heatmaps", "combined")
  
  for (dir_name in required_dirs) {
    if (is.null(CONFIG$dirs[[dir_name]])) {
      stop(sprintf("❌ CONFIG$dirs$%s 未定义", dir_name))
    }
    if (!dir.exists(CONFIG$dirs[[dir_name]])) {
      dir.create(CONFIG$dirs[[dir_name]], recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  invisible(TRUE)
}


#' 验证必需函数
validate_required_functions <- function() {
  
  required_functions <- c(
    "calculate_density_zones",
    "plot_celltype_density_overlay",
    "plot_zone_composition",
    "plot_combined_heatmap",
    "plot_combined_analysis",
    "generate_summary_statistics"
  )
  
  missing_funcs <- required_functions[!sapply(required_functions, exists)]
  
  if (length(missing_funcs) > 0) {
    stop(sprintf("❌ 缺少必需函数: %s", paste(missing_funcs, collapse = ", ")))
  }
  
  invisible(TRUE)
}


#' 设置颜色方案
#' 
#' @param first_sample 第一个样本（用于提取细胞类型）
#' @param CONFIG 配置对象
#' @param celltype_col 细胞类型列名
#' @param density_bins 密度分区数量
setup_colors <- function(first_sample, CONFIG, celltype_col, density_bins) {
  
  # 从第一个样本获取所有细胞类型
  all_celltypes <- sort(unique(as.character(first_sample[[celltype_col]][,1])))
  
  if (is.null(CONFIG$colors$celltype_colors)) {
    CONFIG$colors$celltype_colors <- get_celltype_colors(all_celltypes)
    cat(sprintf("🎨 已生成 %d 种细胞类型颜色方案\n", length(CONFIG$colors$celltype_colors)))
  }
  
  if (is.null(CONFIG$colors$zone_colors)) {
    CONFIG$colors$zone_colors <- get_zone_colors(density_bins)
  }
  
  invisible(CONFIG)
}


#' 验证样本数据
#' 
#' @param seurat_subset Seurat 对象
#' @param sample_id 样本 ID
#' @param celltype_col 细胞类型列名
#' 
#' @return 验证结果列表
validate_sample_data <- function(seurat_subset, sample_id, celltype_col) {
  
  # 检查数据量
  if (ncol(seurat_subset) == 0) {
    return(list(valid = FALSE, error = "No data"))
  }
  
  # 获取坐标
  coords <- tryCatch({
    Seurat::GetTissueCoordinates(
      seurat_subset,
      cols = c("row", "col"),
      scale = NULL
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(coords)) {
    return(list(valid = FALSE, error = "Cannot get coordinates"))
  }
  
  # 合并元数据
  df <- seurat_subset@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(
      coords %>% tibble::rownames_to_column("barcode"), 
      by = "barcode"
    ) %>%
    dplyr::filter(!is.na(col), !is.na(row))
  
  if (nrow(df) == 0) {
    return(list(valid = FALSE, error = "No valid coordinates"))
  }
  
  # 检查必需列
  if (!celltype_col %in% colnames(df)) {
    return(list(valid = FALSE, error = "Missing celltype column"))
  }
  
  if (!"ClockGene_High" %in% colnames(df)) {
    return(list(valid = FALSE, error = "Missing ClockGene_High column"))
  }
  
  # 清理细胞类型
  df$celltype_clean <- as.character(df[[celltype_col]])
  df$celltype_clean[is.na(df$celltype_clean)] <- "Unknown"
  
  return(list(valid = TRUE, df = df))
}

cat("✅ 08_validation.R 已加载\n")