# ===================================================================
# 批量修复所有 filter() 调用
# ===================================================================

fix_filter_in_files <- function() {
  
  cat("🔍 开始批量修复 filter() 调用...\n\n")
  
  # 需要修复的文件列表
  files_to_fix <- c(
    "07_plot_spatial.R",
    "niche_grade_entropy.R",
    "niche_marker.R",
    "SSS_isoheight_plot.R",
    "08_plot_celltype.R",
    "08_plot_celltype_utils/02_density_zones.R",
    "08_plot_celltype_utils/03_plot_overlay.R",
    "08_plot_celltype_utils/07_statistics.R"
  )
  
  fixed_count <- 0
  skipped_count <- 0
  
  for (file in files_to_fix) {
    
    if (!file.exists(file)) {
      cat(sprintf("⚠️  跳过: %s (文件不存在)\n", file))
      skipped_count <- skipped_count + 1
      next
    }
    
    # 备份原文件
    backup_file <- paste0(file, ".backup.", format(Sys.time(), "%Y%m%d_%H%M%S"))
    file.copy(file, backup_file, overwrite = FALSE)
    
    # 读取文件内容
    content <- readLines(file, warn = FALSE)
    original_content <- content
    
    # 修复模式
    patterns <- list(
      # filter( -> dplyr::filter(
      list(
        pattern = "([^:_a-zA-Z0-9])filter\\s*\\(",
        replacement = "\\1dplyr::filter("
      ),
      # left_join( -> dplyr::left_join(
      list(
        pattern = "([^:_a-zA-Z0-9])left_join\\s*\\(",
        replacement = "\\1dplyr::left_join("
      ),
      # right_join( -> dplyr::right_join(
      list(
        pattern = "([^:_a-zA-Z0-9])right_join\\s*\\(",
        replacement = "\\1dplyr::right_join("
      ),
      # inner_join( -> dplyr::inner_join(
      list(
        pattern = "([^:_a-zA-Z0-9])inner_join\\s*\\(",
        replacement = "\\1dplyr::inner_join("
      ),
      # rownames_to_column( -> tibble::rownames_to_column(
      list(
        pattern = "([^:_a-zA-Z0-9])rownames_to_column\\s*\\(",
        replacement = "\\1tibble::rownames_to_column("
      ),
      # column_to_rownames( -> tibble::column_to_rownames(
      list(
        pattern = "([^:_a-zA-Z0-9])column_to_rownames\\s*\\(",
        replacement = "\\1tibble::column_to_rownames("
      ),
      # select( -> dplyr::select(
      list(
        pattern = "([^:_a-zA-Z0-9])select\\s*\\(",
        replacement = "\\1dplyr::select("
      ),
      # mutate( -> dplyr::mutate(
      list(
        pattern = "([^:_a-zA-Z0-9])mutate\\s*\\(",
        replacement = "\\1dplyr::mutate("
      ),
      # group_by( -> dplyr::group_by(
      list(
        pattern = "([^:_a-zA-Z0-9])group_by\\s*\\(",
        replacement = "\\1dplyr::group_by("
      ),
      # summarize( -> dplyr::summarize(
      list(
        pattern = "([^:_a-zA-Z0-9])summarize\\s*\\(",
        replacement = "\\1dplyr::summarize("
      ),
      # summarise( -> dplyr::summarise(
      list(
        pattern = "([^:_a-zA-Z0-9])summarise\\s*\\(",
        replacement = "\\1dplyr::summarise("
      ),
      # arrange( -> dplyr::arrange(
      list(
        pattern = "([^:_a-zA-Z0-9])arrange\\s*\\(",
        replacement = "\\1dplyr::arrange("
      )
    )
    
    # 应用所有修复
    for (pat in patterns) {
      content <- gsub(pat$pattern, pat$replacement, content, perl = TRUE)
    }
    
    # 检查是否有修改
    if (identical(content, original_content)) {
      cat(sprintf("✓  无需修改: %s\n", file))
    } else {
      # 写回文件
      writeLines(content, file)
      
      # 统计修改次数
      n_changes <- sum(content != original_content)
      
      cat(sprintf("✅ 已修复: %s (%d 行修改, 备份: %s)\n", 
                  file, n_changes, basename(backup_file)))
      fixed_count <- fixed_count + 1
    }
  }
  
  cat("\n")
  cat(paste(rep("═", 70), collapse = ""), "\n")
  cat(sprintf("✅ 修复完成: %d 个文件已修改\n", fixed_count))
  if (skipped_count > 0) {
    cat(sprintf("⚠️  跳过: %d 个文件\n", skipped_count))
  }
  cat(paste(rep("═", 70), collapse = ""), "\n\n")
  
  cat("📝 下一步:\n")
  cat("   1. 重启 R 会话: .rs.restartR()\n")
  cat("   2. 重新运行: source('main.R'); main_batch()\n\n")
  
  invisible(list(
    fixed = fixed_count,
    skipped = skipped_count,
    total = length(files_to_fix)
  ))
}

# 运行修复
fix_filter_in_files()