# ===================================================================
# 07_statistics.R (全局统一配色版)
# 统计摘要生成
# Author: Assistant
# Date: 2025-11-07
# ===================================================================

#' 生成统计摘要
#'
#' @param combined_data 合并的zone组成数据（必须包含 celltype_clean, density_zone, percentage, sample）
#' @param CONFIG 配置列表（可选，用于获取全局细胞类型颜色）
#'
#' @return 统计摘要数据框，包含以下列：
#'   - celltype_clean: 细胞类型名称
#'   - mean_pct_all: 跨所有zone的平均百分比
#'   - sd_pct_all: 标准差
#'   - median_pct_all: 中位数
#'   - max_zone: 富集最高的zone
#'   - max_pct: 最高百分比
#'   - min_zone: 富集最低的zone
#'   - min_pct: 最低百分比
#'   - core_enrichment: 核心区（Zone_0, Zone_1）与外围区的差异
#'   - enrichment_type: 富集类型（Core-enriched/Peripheral-enriched/Evenly-distributed）
#'   - n_samples: 包含该细胞类型的样本数
#'   - n_observations: 观测次数
#'
#' @details
#' 富集类型判定标准：
#' - Core-enriched: core_enrichment > 5
#' - Peripheral-enriched: core_enrichment < -5
#' - Evenly-distributed: -5 <= core_enrichment <= 5
#'
#' @examples
#' summary_stats <- generate_summary_statistics(combined_data)
#' write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)
#'
generate_summary_statistics <- function(combined_data, CONFIG = NULL) {
  
  require(dplyr)
  
  cat("\n📊 生成统计摘要...\n")
  
  # ========================================
  # 1. 数据验证
  # ========================================
  
  required_cols <- c("celltype_clean", "density_zone", "percentage", "sample")
  missing_cols <- required_cols[!required_cols %in% colnames(combined_data)]
  
  if (length(missing_cols) > 0) {
    stop(sprintf("❌ 数据缺少必需列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  if (nrow(combined_data) == 0) {
    stop("❌ combined_data 为空")
  }
  
  # ========================================
  # 2. 计算统计摘要
  # ========================================
  
  summary <- combined_data %>%
    dplyr::mutate(
      # 提取zone编号
      zone_numeric = as.numeric(gsub("Zone_", "", density_zone))
    ) %>%
    dplyr::group_by(celltype_clean) %>%
    dplyr::summarise(
      # 基础统计量
      mean_pct_all = mean(percentage, na.rm = TRUE),
      sd_pct_all = sd(percentage, na.rm = TRUE),
      median_pct_all = median(percentage, na.rm = TRUE),
      
      # 富集最高的区域
      max_zone = density_zone[which.max(percentage)],
      max_pct = max(percentage, na.rm = TRUE),
      
      # 富集最低的区域
      min_zone = density_zone[which.min(percentage)],
      min_pct = min(percentage, na.rm = TRUE),
      
      # 核心区（Zone_0 和 Zone_1）vs 外围区富集差异
      # 正值表示核心富集，负值表示外围富集
      core_enrichment = mean(percentage[zone_numeric <= 1], na.rm = TRUE) - 
                       mean(percentage[zone_numeric > 1], na.rm = TRUE),
      
      # 样本覆盖度
      n_samples = length(unique(sample)),
      n_observations = n(),
      
      .groups = "drop"
    ) %>%
    # ========================================
    # 3. 分类富集类型
    # ========================================
    dplyr::mutate(
      enrichment_type = dplyr::case_when(
        core_enrichment > 5 ~ "Core-enriched",
        core_enrichment < -5 ~ "Peripheral-enriched",
        TRUE ~ "Evenly-distributed"
      )
    ) %>%
    # 按核心富集度排序（核心富集的细胞类型排在前面）
    dplyr::arrange(desc(core_enrichment))
  
  # ========================================
  # 4. 添加颜色信息（如果提供了CONFIG）
  # ========================================
  
  if (!is.null(CONFIG) && !is.null(CONFIG$colors) && !is.null(CONFIG$colors$celltype)) {
    celltype_colors <- CONFIG$colors$celltype
    
    summary <- summary %>%
      dplyr::mutate(
        color = sapply(celltype_clean, function(ct) {
          celltype_colors[ct] %||% "#CCCCCC"
        })
      )
  }
  
  # ========================================
  # 5. 打印摘要
  # ========================================
  
  cat(sprintf("   📊 分析了 %d 种细胞类型\n", nrow(summary)))
  cat(sprintf("   📈 总观测次数: %d\n", sum(summary$n_observations)))
  cat(sprintf("   📦 来自 %d 个样本\n", length(unique(combined_data$sample))))
  cat("\n")
  
  # 富集类型统计
  n_core <- sum(summary$enrichment_type == "Core-enriched")
  n_periph <- sum(summary$enrichment_type == "Peripheral-enriched")
  n_even <- sum(summary$enrichment_type == "Evenly-distributed")
  
  cat("   📍 富集类型分布:\n")
  cat(sprintf("      🔴 核心富集 (Core):        %2d (%.1f%%)\n", 
              n_core, 100 * n_core / nrow(summary)))
  cat(sprintf("      🔵 外围富集 (Peripheral):  %2d (%.1f%%)\n", 
              n_periph, 100 * n_periph / nrow(summary)))
  cat(sprintf("      ⚪ 均匀分布 (Even):         %2d (%.1f%%)\n", 
              n_even, 100 * n_even / nrow(summary)))
  cat("\n")
  
  # ========================================
  # 6. 打印核心富集 TOP 5
  # ========================================
  
  if (n_core > 0) {
    cat("   🔴 核心富集 TOP 5:\n")
    top_core <- head(summary %>% dplyr::filter(enrichment_type == "Core-enriched"), 5)
    
    for (i in 1:nrow(top_core)) {
      cat(sprintf("      %d. %-30s: +%6.2f%% (最高: %s, %.2f%%)\n", 
                  i, 
                  top_core$celltype_clean[i], 
                  top_core$core_enrichment[i],
                  top_core$max_zone[i],
                  top_core$max_pct[i]))
    }
    cat("\n")
  }
  
  # ========================================
  # 7. 打印外围富集 TOP 5
  # ========================================
  
  if (n_periph > 0) {
    cat("   🔵 外围富集 TOP 5:\n")
    top_periph <- head(
      summary %>% 
        dplyr::filter(enrichment_type == "Peripheral-enriched") %>%
        dplyr::arrange(core_enrichment), 
      5
    )
    
    for (i in 1:nrow(top_periph)) {
      cat(sprintf("      %d. %-30s: %7.2f%% (最高: %s, %.2f%%)\n", 
                  i, 
                  top_periph$celltype_clean[i], 
                  top_periph$core_enrichment[i],
                  top_periph$max_zone[i],
                  top_periph$max_pct[i]))
    }
    cat("\n")
  }
  
  # ========================================
  # 8. 打印均匀分布细胞类型
  # ========================================
  
  if (n_even > 0 && n_even <= 5) {
    cat("   ⚪ 均匀分布细胞类型:\n")
    even_types <- summary %>% dplyr::filter(enrichment_type == "Evenly-distributed")
    
    for (i in 1:nrow(even_types)) {
      cat(sprintf("      • %-30s: %7.2f%% (范围: %.2f%% - %.2f%%)\n", 
                  even_types$celltype_clean[i], 
                  even_types$core_enrichment[i],
                  even_types$min_pct[i],
                  even_types$max_pct[i]))
    }
    cat("\n")
  }
  
  cat("   ✅ 统计摘要生成完成\n")
  
  return(summary)
}


#' 生成zone级别的统计摘要
#'
#' @param combined_data 合并的zone组成数据
#'
#' @return zone级别的统计摘要
#'
generate_zone_summary <- function(combined_data) {
  
  require(dplyr)
  
  cat("\n📊 生成zone级别统计...\n")
  
  zone_summary <- combined_data %>%
    dplyr::group_by(density_zone) %>%
    dplyr::summarise(
      n_celltypes = length(unique(celltype_clean)),
      n_samples = length(unique(sample)),
      n_observations = n(),
      mean_diversity = mean(percentage),  # 平均富集度（多样性指标）
      sd_diversity = sd(percentage),
      .groups = "drop"
    ) %>%
    dplyr::arrange(density_zone)
  
  cat(sprintf("   📊 分析了 %d 个密度区域\n", nrow(zone_summary)))
  
  for (i in 1:nrow(zone_summary)) {
    cat(sprintf("   • %-10s: %2d 种细胞类型, %2d 个样本\n",
                zone_summary$density_zone[i],
                zone_summary$n_celltypes[i],
                zone_summary$n_samples[i]))
  }
  
  cat("   ✅ Zone统计完成\n")
  
  return(zone_summary)
}

cat("✅ 07_statistics.R 已加载（全局统一配色版）\n")