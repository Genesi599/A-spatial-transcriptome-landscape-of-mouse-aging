# ===================================================================
# 08_validation.R (简化版)
# 数据验证模块
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 验证配置对象
#' 
#' @param CONFIG 配置对象
#' @return invisible(TRUE)
#'
#' @details
#' 验证内容：
#' - 必需参数是否存在
#' - 输出目录是否定义
#' - 自动创建不存在的目录
#'
validate_config <- function(CONFIG) {
  
  cat("\n🔍 验证配置...\n")
  
  # ========================================
  # 1. 验证必需参数
  # ========================================
  
  required_params <- c(
    "col_col",                        # 列坐标列名
    "row_col",                        # 行坐标列名
    "celltype_col",                   # 细胞类型列名
    "density_threshold_percentile",   # 密度阈值百分位
    "n_zones"                         # 密度区域数量
  )
  
  if (is.null(CONFIG$params)) {
    stop("❌ CONFIG$params 未定义")
  }
  
  missing_params <- required_params[!required_params %in% names(CONFIG$params)]
  
  if (length(missing_params) > 0) {
    stop(sprintf("❌ 缺少必需参数: %s", paste(missing_params, collapse = ", ")))
  }
  
  cat("   ✅ 必需参数完整\n")
  
  # ========================================
  # 2. 验证参数值的合理性
  # ========================================
  
  # 验证密度阈值百分位
  if (CONFIG$params$density_threshold_percentile < 0 || 
      CONFIG$params$density_threshold_percentile > 1) {
    stop("❌ density_threshold_percentile 必须在 [0, 1] 范围内")
  }
  
  # 验证zone数量
  if (CONFIG$params$n_zones < 2 || CONFIG$params$n_zones > 20) {
    warning("⚠️  n_zones 通常在 2-20 之间，当前值可能不合理")
  }
  
  cat("   ✅ 参数值合理\n")
  
  # ========================================
  # 3. 验证输出目录
  # ========================================
  
  if (is.null(CONFIG$output)) {
    stop("❌ CONFIG$output 未定义")
  }
  
  if (is.null(CONFIG$output$plot_dir) || is.null(CONFIG$output$data_dir)) {
    stop("❌ 必须指定 CONFIG$output$plot_dir 和 CONFIG$output$data_dir")
  }
  
  # 创建输出目录
  for (dir_path in c(CONFIG$output$plot_dir, CONFIG$output$data_dir)) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("   📁 创建目录: %s\n", dir_path))
    } else {
      cat(sprintf("   📁 目录已存在: %s\n", dir_path))
    }
  }
  
  cat("   ✅ 输出目录就绪\n")
  
  invisible(TRUE)
}


#' 验证数据列表
#' 
#' @param data_list 数据框列表
#' @param sample_ids 样本ID向量
#' @param CONFIG 配置对象
#' @return invisible(TRUE)
#'
#' @details
#' 验证内容：
#' - 数据列表结构
#' - 样本ID数量匹配
#' - 每个数据框的必需列
#' - 数据量检查
#'
validate_data_list <- function(data_list, sample_ids, CONFIG) {
  
  cat("\n🔍 验证数据列表...\n")
  
  # ========================================
  # 1. 验证基本结构
  # ========================================
  
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("❌ data_list 必须是非空列表")
  }
  
  if (length(data_list) != length(sample_ids)) {
    stop(sprintf("❌ data_list 长度 (%d) 与 sample_ids 长度 (%d) 不匹配",
                length(data_list), length(sample_ids)))
  }
  
  cat(sprintf("   📦 数据列表包含 %d 个样本\n", length(data_list)))
  
  # ========================================
  # 2. 验证每个数据框
  # ========================================
  
  n_valid <- 0
  n_warnings <- 0
  
  required_cols <- c(
    CONFIG$params$col_col, 
    CONFIG$params$row_col, 
    CONFIG$params$celltype_col
  )
  
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    sid <- sample_ids[i]
    
    # 检查是否为数据框
    if (!is.data.frame(df)) {
      warning(sprintf("   ⚠️  样本 %s: 不是数据框，跳过", sid))
      n_warnings <- n_warnings + 1
      next
    }
    
    # 检查必需列
    missing_cols <- required_cols[!required_cols %in% colnames(df)]
    
    if (length(missing_cols) > 0) {
      warning(sprintf("   ⚠️  样本 %s: 缺少列 %s", 
                     sid, paste(missing_cols, collapse = ", ")))
      n_warnings <- n_warnings + 1
      next
    }
    
    # 检查数据量
    if (nrow(df) == 0) {
      warning(sprintf("   ⚠️  样本 %s: 数据为空", sid))
      n_warnings <- n_warnings + 1
      next
    }
    
    # 检查坐标是否有效
    n_na_coords <- sum(is.na(df[[CONFIG$params$col_col]]) | 
                      is.na(df[[CONFIG$params$row_col]]))
    
    if (n_na_coords > 0) {
      warning(sprintf("   ⚠️  样本 %s: %d 个spots坐标缺失", sid, n_na_coords))
      n_warnings <- n_warnings + 1
    }
    
    n_valid <- n_valid + 1
  }
  
  # ========================================
  # 3. 汇总验证结果
  # ========================================
  
  cat(sprintf("   ✅ 有效样本: %d/%d\n", n_valid, length(data_list)))
  
  if (n_warnings > 0) {
    cat(sprintf("   ⚠️  警告: %d 个\n", n_warnings))
  }
  
  if (n_valid == 0) {
    stop("❌ 没有有效的样本数据")
  }
  
  cat("   ✅ 数据验证通过\n")
  
  invisible(TRUE)
}


#' 验证必需函数是否加载
#' 
#' @return invisible(TRUE)
#'
#' @details
#' 检查所有必需的函数是否已加载到环境中
#'
validate_required_functions <- function() {
  
  cat("\n🔍 验证必需函数...\n")
  
  required_functions <- c(
    # 核心功能
    "standardize_celltype_names",      # 名称标准化
    "create_global_color_scheme",       # 全局颜色方案
    "calculate_density_zones",          # 密度区域计算
    
    # 绘图函数
    "plot_celltype_density_overlay",    # 叠加图
    "plot_zone_composition",            # 组成图
    "plot_combined_heatmap",            # 热图
    "plot_combined_analysis",           # 综合分析图
    
    # 统计分析
    "generate_summary_statistics",      # 统计摘要
    
    # 操作符
    "%||%"                              # 空值默认值
  )
  
  missing_funcs <- character()
  
  for (func_name in required_functions) {
    if (!exists(func_name, mode = "function")) {
      missing_funcs <- c(missing_funcs, func_name)
    }
  }
  
  if (length(missing_funcs) > 0) {
    stop(sprintf(
      "❌ 缺少必需函数: %s\n\n请检查是否加载了所有工具函数文件:\n%s",
      paste(missing_funcs, collapse = ", "),
      paste(
        "  • 00_operators.R",
        "  • 01_color_schemes.R",
        "  • 02_density_zones.R",
        "  • 03_plot_overlay.R",
        "  • 04_plot_composition.R",
        "  • 05_plot_heatmap.R",
        "  • 06_plot_combined.R",
        "  • 07_statistics.R",
        sep = "\n"
      )
    ))
  }
  
  cat(sprintf("   ✅ 所有 %d 个必需函数已加载\n", length(required_functions)))
  
  invisible(TRUE)
}


#' 快速数据质量检查
#' 
#' @param df 数据框
#' @param sample_id 样本ID
#' @param CONFIG 配置对象
#' @return 质量检查结果列表
#'
quick_quality_check <- function(df, sample_id, CONFIG) {
  
  result <- list(
    sample_id = sample_id,
    passed = TRUE,
    warnings = character(),
    info = list()
  )
  
  # 检查数据量
  result$info$n_rows <- nrow(df)
  
  if (nrow(df) < 100) {
    result$warnings <- c(result$warnings, "数据量较少（< 100 spots）")
  }
  
  # 检查坐标范围
  col_range <- range(df[[CONFIG$params$col_col]], na.rm = TRUE)
  row_range <- range(df[[CONFIG$params$row_col]], na.rm = TRUE)
  
  result$info$col_range <- col_range
  result$info$row_range <- row_range
  
  # 检查细胞类型数量
  n_celltypes <- length(unique(df[[CONFIG$params$celltype_col]]))
  result$info$n_celltypes <- n_celltypes
  
  if (n_celltypes < 2) {
    result$warnings <- c(result$warnings, "细胞类型少于2种")
  }
  
  # 检查缺失值
  n_missing_celltype <- sum(is.na(df[[CONFIG$params$celltype_col]]))
  
  if (n_missing_celltype > 0) {
    result$warnings <- c(result$warnings, 
                        sprintf("%d spots 缺少细胞类型信息", n_missing_celltype))
  }
  
  return(result)
}

cat("✅ 08_validation.R 已加载（简化版）\n")