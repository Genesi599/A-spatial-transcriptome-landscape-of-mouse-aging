# ===================================================================
# 精确修复 calculate_density_zones 中的 filter
# ===================================================================

fix_density_zones <- function() {
  
  file <- "08_plot_celltype_utils/02_density_zones.R"
  
  if (!file.exists(file)) {
    stop("❌ 文件不存在: ", file)
  }
  
  cat(sprintf("🔍 检查文件: %s\n\n", file))
  
  # 备份
  backup <- paste0(file, ".backup.", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(file, backup)
  cat(sprintf("💾 备份: %s\n\n", basename(backup)))
  
  # 读取
  lines <- readLines(file, warn = FALSE)
  original_lines <- lines
  
  # 显示所有包含 filter 的行（修复前）
  cat("🔍 修复前的 filter 调用:\n")
  cat(paste(rep("─", 70), collapse = ""), "\n")
  
  for (i in seq_along(lines)) {
    if (grepl("filter\\s*\\(", lines[i], perl = TRUE) && 
        !grepl("^\\s*#", lines[i])) {
      cat(sprintf("%4d: %s\n", i, trimws(lines[i])))
    }
  }
  cat(paste(rep("─", 70), collapse = ""), "\n\n")
  
  # 修复模式
  fixes <- list(
    # filter( -> dplyr::filter(
    list(
      name = "filter",
      pattern = "([^:_a-zA-Z0-9])filter\\s*\\(",
      replacement = "\\1dplyr::filter("
    ),
    # select( -> dplyr::select(
    list(
      name = "select",
      pattern = "([^:_a-zA-Z0-9])select\\s*\\(",
      replacement = "\\1dplyr::select("
    ),
    # mutate( -> dplyr::mutate(
    list(
      name = "mutate",
      pattern = "([^:_a-zA-Z0-9])mutate\\s*\\(",
      replacement = "\\1dplyr::mutate("
    ),
    # left_join( -> dplyr::left_join(
    list(
      name = "left_join",
      pattern = "([^:_a-zA-Z0-9])left_join\\s*\\(",
      replacement = "\\1dplyr::left_join("
    ),
    # group_by( -> dplyr::group_by(
    list(
      name = "group_by",
      pattern = "([^:_a-zA-Z0-9])group_by\\s*\\(",
      replacement = "\\1dplyr::group_by("
    ),
    # summarize( -> dplyr::summarize(
    list(
      name = "summarize",
      pattern = "([^:_a-zA-Z0-9])summarize\\s*\\(",
      replacement = "\\1dplyr::summarize("
    ),
    # arrange( -> dplyr::arrange(
    list(
      name = "arrange",
      pattern = "([^:_a-zA-Z0-9])arrange\\s*\\(",
      replacement = "\\1dplyr::arrange("
    )
  )
  
  # 应用所有修复
  for (fix in fixes) {
    lines <- gsub(fix$pattern, fix$replacement, lines, perl = TRUE)
  }
  
  # 显示修复后的结果
  cat("✅ 修复后的 filter 调用:\n")
  cat(paste(rep("─", 70), collapse = ""), "\n")
  
  for (i in seq_along(lines)) {
    if (grepl("filter\\s*\\(", lines[i], perl = TRUE) && 
        !grepl("^\\s*#", lines[i])) {
      cat(sprintf("%4d: %s\n", i, trimws(lines[i])))
    }
  }
  cat(paste(rep("─", 70), collapse = ""), "\n\n")
  
  # 统计变化
  n_changes <- sum(lines != original_lines)
  
  if (n_changes == 0) {
    cat("ℹ️  文件无需修改\n\n")
    return(invisible(NULL))
  }
  
  # 写回文件
  writeLines(lines, file)
  
  cat(sprintf("✅ 已修复 %d 行\n", n_changes))
  cat(sprintf("📁 文件: %s\n", file))
  cat(sprintf("💾 备份: %s\n\n", backup))
  
  return(invisible(list(
    file = file,
    changes = n_changes,
    backup = backup
  )))
}

# 执行修复
fix_density_zones()

cat("🔧 修复完成！现在重新运行:\n")
cat("   .rs.restartR()\n")
cat("   source('main.R')\n")
cat("   main_batch()\n\n")