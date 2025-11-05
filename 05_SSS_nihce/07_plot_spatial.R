#!/usr/bin/env Rscript
# ===================================================================
# 空间梯度图绘制
# ===================================================================

plot_spatial_gradient <- function(seurat_obj, samples_to_plot, config) {
  cat("🔥 绘制空间梯度图...\n")
  
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
    
    seurat_subset <- tryCatch(
      subset(seurat_obj, subset = orig.ident == sample_id),
      error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
    )
    
    tryCatch({
      # 获取坐标
      coords <- get_coordinates_safely(seurat_subset)
      
      # 合并数据
      plot_data <- seurat_subset@meta.data %>%
        rownames_to_column("barcode") %>%
        left_join(coords %>% rownames_to_column("barcode"), by = "barcode")
      
      # 计算坐标范围
      limits <- calculate_coord_limits(plot_data, config$plot$expand_margin)
      
      # 左图：Score
      p_score <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_point(aes(fill = ClockGene_Score1), 
                   shape = 21, size = config$plot$point_size_scatter, 
                   color = "white", stroke = 0.1) +
        scale_fill_gradientn(
          colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
          name = "Clock Gene\nScore"
        ) +
        scale_x_continuous(limits = limits$col, expand = expansion(mult = 0.02)) +
        scale_y_reverse(limits = rev(limits$row), expand = expansion(mult = 0.02)) +
        coord_fixed(ratio = 1) +
        ggtitle("Clock Gene Score") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right"
        )
      
      # 右图：Distance
      p_distance <- ggplot(plot_data, aes(x = col, y = row)) +
        geom_point(aes(fill = ClockGene_Distance), 
                   shape = 21, size = config$plot$point_size_scatter, 
                   color = "white", stroke = 0.1) +
        scale_fill_gradientn(
          colors = rev(c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027")),
          name = "Distance"
        ) +
        scale_x_continuous(limits = limits$col, expand = expansion(mult = 0.02)) +
        scale_y_reverse(limits = rev(limits$row), expand = expansion(mult = 0.02)) +
        coord_fixed(ratio = 1) +
        ggtitle("Distance to High Score Region") +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right"
        )
      
      # 合并
      p_combined <- (p_score | p_distance) +
        plot_annotation(
          title = sprintf("Clock Gene Niche Analysis - %s", sample_id),
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )
      
      # 保存
      safe_name <- safe_filename(sample_id)
      ggsave(
        file.path(config$dirs$spatial, sprintf("ClockGene_spatial_%s.pdf", safe_name)),
        plot = p_combined, 
        width = 16, height = 8, 
        dpi = config$plot$dpi
      )
      
    }, error = function(e) {
      cat(sprintf("   ⚠️ 绘图失败: %s\n", e$message))
    })
  }
  
  cat("✅ 空间梯度图绘制完成\n\n")
}