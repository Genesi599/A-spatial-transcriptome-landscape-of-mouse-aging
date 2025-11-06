# ===================================================================
# 07_statistics.R
# 统计摘要生成
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' 生成统计摘要
#'
#' @param combined_data 合并的zone组成数据
#'
#' @return 统计摘要数据框
#'
#' @details
#' 计算内容：
#' - 每种细胞类型的平均百分比和标准差
#' - 富集最多/最少的zone
#' - 核心区（Zone_0和Zone_1）vs外围区的富集差异
#' - 样本数量
#'
#' @examples
#' summary <- generate_summary_statistics(combined_data)
#'
generate_summary_statistics <- function(combined_data) {
  
  require(dplyr)
  
  # 计算每种细胞类型在不同区域的富集情况
  summary <- combined_data %>%
    dplyr::mutate(zone_numeric = as.numeric(gsub("Zone_", "", density_zone))) %>%
    dplyr::group_by(celltype_clean) %>%
    dplyr::summarise(
      mean_pct_all = mean(percentage),
      sd_pct_all = sd(percentage),
      max_zone = density_zone[which.max(percentage)],
      max_pct = max(percentage),
      min_zone = density_zone[which.min(percentage)],
      min_pct = min(percentage),
      # Zone_0和Zone_1是核心区，其他是外围
      core_enrichment = mean(percentage[zone_numeric <= 1]) - mean(percentage[zone_numeric > 1]),
      n_samples = length(unique(sample)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(core_enrichment))
  
  return(summary)
}