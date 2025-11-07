#!/usr/bin/env Rscript
# ===================================================================
# 图形保存模块（支持串联和并行模式）
# ===================================================================

#' 处理单个样本
#' 
#' @param sample_id 样本 ID
#' @param sample_list 样本列表
#' @param CONFIG 配置对象
#' @param celltype_col 细胞类型列名
#' @param density_bins 密度分区数量
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' @param progressor 进度条对象（可选，串联模式下为 NULL）
#' 
#' @return 处理结果
process_single_sample <- function(sample_id, sample_list, CONFIG, 
                                  celltype_col, density_bins,
                                  plot_overlay, plot_composition,
                                  progressor = NULL) {
  
  tryCatch({
    
    # 1. 获取并验证数据
    seurat_subset <- sample_list[[sample_id]]
    validation <- validate_sample_data(seurat_subset, sample_id, celltype_col)
    
    if (!validation$valid) {
      # ✅ 只在 progressor 存在时调用
      if (!is.null(progressor)) {
        progressor(message = sprintf("⚠️  %s - %s", sample_id, validation$error))
      }
      return(list(sample = sample_id, success = FALSE, error = validation$error))
    }
    
    df <- validation$df
    
    # 统计基本信息
    n_spots <- nrow(df)
    n_high <- sum(df$ClockGene_High, na.rm = TRUE)
    high_pct <- 100 * mean(df$ClockGene_High, na.rm = TRUE)
    
    # 2. 计算密度区域
    density_data <- calculate_density_zones(
      df = df,
      density_bins = density_bins,
      expand_margin = CONFIG$plot$expand_margin %||% 0.1
    )
    
    if (is.null(density_data)) {
      # ✅ 只在 progressor 存在时调用
      if (!is.null(progressor)) {
        progressor(message = sprintf("⚠️  %s - 密度计算失败", sample_id))
      }
      return(list(sample = sample_id, success = FALSE, error = "Density calculation failed"))
    }
    
    # 合并密度信息
    df <- df %>%
      dplyr::left_join(
        density_data$spot_zones %>% 
          dplyr::select(col, row, density_zone, density_value),
        by = c("col", "row")
      )
    
    # 3. 计算区域组成
    zone_composition <- df %>%
      dplyr::filter(!is.na(density_zone)) %>%
      dplyr::group_by(density_zone, celltype_clean) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(density_zone) %>%
      dplyr::mutate(
        total = sum(count),
        percentage = 100 * count / total
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sample = sample_id)
    
    # 4. 绘制并保存图形
    plot_result <- save_sample_plots(
      df = df,
      density_data = density_data,
      zone_composition = zone_composition,
      sample_id = sample_id,
      CONFIG = CONFIG,
      plot_overlay = plot_overlay,
      plot_composition = plot_composition
    )
    
    # ✅ 只在 progressor 存在时更新进度
    if (!is.null(progressor)) {
      progressor(message = sprintf("✅ %s (%.2f MB)", sample_id, plot_result$total_size_mb))
    }
    
    # 5. 返回结果
    return(list(
      sample = sample_id,
      success = TRUE,
      zone_composition = zone_composition,
      n_spots = n_spots,
      n_high = n_high,
      high_pct = high_pct,
      n_zones = length(unique(zone_composition$density_zone)),
      n_celltypes = length(unique(zone_composition$celltype_clean)),
      n_na_zones = sum(is.na(df$density_zone)),
      output_files = plot_result$output_files,
      total_size_mb = plot_result$total_size_mb
    ))
    
  }, error = function(e) {
    # ✅ 只在 progressor 存在时更新进度
    if (!is.null(progressor)) {
      progressor(message = sprintf("❌ %s - %s", sample_id, e$message))
    }
    return(list(
      sample = sample_id,
      success = FALSE,
      error = as.character(e$message)
    ))
  })
}


#' 保存样本图形
#' 
#' @param df 数据框
#' @param density_data 密度数据
#' @param zone_composition 区域组成
#' @param sample_id 样本 ID
#' @param CONFIG 配置对象
#' @param plot_overlay 是否绘制叠加图
#' @param plot_composition 是否绘制组成图
#' 
#' @return 保存结果
save_sample_plots <- function(df, density_data, zone_composition, sample_id, CONFIG,
                              plot_overlay, plot_composition) {
  
  output_files <- list()
  total_size <- 0
  safe_name <- gsub("[^[:alnum:]]", "_", sample_id)
  
  # 叠加图
  if (plot_overlay) {
    p_overlay <- plot_celltype_density_overlay(
      df = df,
      density_data = density_data,
      sample_id = sample_id,
      CONFIG = CONFIG
    )
    
    overlay_file <- file.path(
      CONFIG$dirs$overlay, 
      sprintf("celltype_overlay_%s.pdf", safe_name)
    )
    
    ggplot2::ggsave(
      overlay_file,
      plot = p_overlay,
      width = 12, 
      height = 10,
      dpi = CONFIG$plot$dpi %||% 300,
      bg = "white"
    )
    
    output_files$overlay <- overlay_file
    total_size <- total_size + file.size(overlay_file)
  }
  
  # 组成图
  if (plot_composition) {
    p_comp <- plot_zone_composition(
      zone_composition = zone_composition,
      sample_id = sample_id,
      CONFIG = CONFIG
    )
    
    composition_file <- file.path(
      CONFIG$dirs$composition, 
      sprintf("composition_%s.pdf", safe_name)
    )
    
    ggplot2::ggsave(
      composition_file,
      plot = p_comp,
      width = 12, 
      height = 6,
      dpi = CONFIG$plot$dpi %||% 300,
      bg = "white"
    )
    
    output_files$composition <- composition_file
    total_size <- total_size + file.size(composition_file)
  }
  
  return(list(
    output_files = output_files,
    total_size_mb = total_size / 1024^2
  ))
}

cat("✅ 09_save_plots.R 已加载\n")