#!/usr/bin/env Rscript
# ===================================================================
# Isoheight 图绘制
# ===================================================================

plot_isoheight_all <- function(seurat_obj, samples_to_plot, config) {
  cat("🎨 绘制 Isoheight 图...\n")
  
  for (i in seq_along(samples_to_plot)) {
    sample_id <- samples_to_plot[i]
    cat(sprintf("[%d/%d] %s\n", i, length(samples_to_plot), sample_id))
    
    seurat_subset <- tryCatch(
      subset(seurat_obj, subset = orig.ident == sample_id),
      error = function(e) seurat_obj[, seurat_obj$orig.ident == sample_id]
    )
    
    p_iso <- celltype_isoheight_plot(
      .data = seurat_subset,
      density_top = ClockGene_High,
      col_bg = "gray92",
      col_top = "#d62728",
      col_isoheight = "white",
      col_white_ratio = 0.25,
      cols_fill_isoheight = c(
        rep("white", 25),
        colorRampPalette(brewer.pal(9, "YlOrRd")[3:9])(75)
      ),
      size_bg = config$plot$point_size_bg,
      size_top = config$plot$point_size_top,
      nrow = 1
    )
    
    safe_name <- safe_filename(sample_id)
    ggsave(
      file.path(config$dirs$isoheight, sprintf("ClockGene_isoheight_%s.pdf", safe_name)),
      plot = p_iso, 
      width = 8, height = 8, 
      dpi = config$plot$dpi
    )
  }
  
  cat("✅ Isoheight 图绘制完成\n\n")
}