#!/usr/bin/env Rscript
# ===================================================================
# 保存结果
# ===================================================================

save_results <- function(seurat_obj, config) {
  cat("💾 保存结果...\n")
  
  # 保存metadata
  write.csv(
    seurat_obj@meta.data, 
    file.path(config$metadata_dir, "Lymph_2-25M_clockgene_metadata.csv"),
    row.names = TRUE
  )
  
  # 可选：保存完整对象
  if (config$save_full_object) {
    saveRDS(
      seurat_obj, 
      file.path(config$metadata_dir, "Lymph_2-25M_with_clockgene_niche.rds")
    )
  }
  
  cat("✅ 结果保存完成\n\n")
}

print_summary <- function(config) {
  cat("📊 文件统计:\n")
  cat(sprintf("   图形文件夹: %s\n", config$figure_dir))
  cat(sprintf("   - Isoheight: %d 个文件\n", length(list.files(config$dirs$isoheight))))
  cat(sprintf("   - Spatial: %d 个文件\n", length(list.files(config$dirs$spatial))))
  cat("\n✅ 全部完成！\n")
}