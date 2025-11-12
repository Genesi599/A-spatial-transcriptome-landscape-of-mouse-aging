#!/usr/bin/env Rscript
# ===================================================================
# 10_batch_processing.R
# ===================================================================

process_all_files <- function(seurat_files, gene_list, CONFIG) {
  is_multi <- is.character(gene_list) && length(gene_list) > 1 && 
              all(file.exists(gene_list))
  
  if (is_multi) {
    return(process_all_files_multi_genelist(
      seurat_files, gene_list, CONFIG))
  } else {
    return(process_all_files_single_genelist(
      seurat_files, gene_list, CONFIG))
  }
}

process_all_files_single_genelist <- function(
    seurat_files, gene_list, CONFIG) {
  
  cat(strrep("=", 60), "\n")
  cat("开始批量处理（单基因列表）\n")
  cat(strrep("=", 60), "\n\n")
  
  results <- list()
  
  for (i in seq_along(seurat_files)) {
    print_progress_header(i, length(seurat_files))
    
    result <- process_seurat_file(
      seurat_path = seurat_files[i],
      gene_list = gene_list,
      base_config = CONFIG
    )
    
    seurat_basename <- tools::file_path_sans_ext(
      basename(seurat_files[i]))
    result$task_id <- seurat_basename
    result$genelist <- NULL
    results[[seurat_basename]] <- result
    
    if (i < length(seurat_files)) {
      estimate_remaining_time(results, seurat_files, i)
    }
    gc(verbose = FALSE)
  }
  return(results)
}

process_all_files_multi_genelist <- function(
    seurat_files, gene_files, CONFIG) {
  
  total_tasks <- length(seurat_files) * length(gene_files)
  completed <- 0
  results <- list()
  
  cat(strrep("=", 60), "\n")
  cat("开始批量处理（多基因列表）\n")
  cat(sprintf("%d Seurat × %d 基因列表 = %d 任务\n",
              length(seurat_files), length(gene_files), total_tasks))
  cat(strrep("=", 60), "\n\n")
  
  for (gene_file in gene_files) {
    genelist_name <- get_genelist_name(gene_file)
    
    cat("\n", strrep("=", 60), "\n")
    cat("基因列表:", genelist_name, "\n")
    cat(strrep("=", 60), "\n\n")
    
    gene_list <- load_gene_list_safe(gene_file)
    cat(sprintf("📋 基因数: %d\n", length(gene_list)))
    
    config_current <- CONFIG
    config_current$output_base_dir <- file.path(
      CONFIG$output_base_dir, genelist_name)
    
    cat(sprintf("📂 输出目录: %s\n\n", 
                config_current$output_base_dir))
    
    for (i in seq_along(seurat_files)) {
      seurat_file <- seurat_files[i]
      seurat_basename <- tools::file_path_sans_ext(
        basename(seurat_file))
      
      task_id <- sprintf("%s_%s", genelist_name, seurat_basename)
      completed <- completed + 1
      
      print_progress_header_multi(
        completed, total_tasks, genelist_name, seurat_basename)
      
      result <- process_seurat_file(
        seurat_path = seurat_file,
        gene_list = gene_list,
        base_config = config_current
      )
      
      result$genelist <- genelist_name
      result$task_id <- task_id
      results[[task_id]] <- result
      
      if (completed < total_tasks) {
        estimate_remaining_time_multi(results, total_tasks, 
                                      completed)
      }
      gc(verbose = FALSE)
    }
    
    merge_genelist_statistics(
      genelist_name, config_current$output_base_dir, seurat_files)
  }
  
  return(results)
}

merge_genelist_statistics <- function(
    genelist_name, base_output_dir, seurat_files) {
  
  cat(sprintf("\n📊 合并 %s 统计数据...\n", genelist_name))
  
  all_stats <- lapply(seurat_files, function(sf) {
    seurat_name <- tools::file_path_sans_ext(basename(sf))
    stat_file <- file.path(
      base_output_dir, seurat_name, 
      "celltype_analysis", "zone_celltype_counts.csv")
    
    if (file.exists(stat_file)) {
      read.csv(stat_file, stringsAsFactors = FALSE)
    } else NULL
  })
  
  all_stats <- do.call(rbind, Filter(Negate(is.null), all_stats))
  
  if (is.null(all_stats) || nrow(all_stats) == 0) {
    warning("无统计数据可合并")
    return(invisible(NULL))
  }
  
  output_file <- file.path(
    base_output_dir, 
    sprintf("%s_zone_celltype_counts_merged.csv", genelist_name))
  
  write.csv(all_stats, output_file, row.names = FALSE)
  cat(sprintf("  ✅ %s (%d 行)\n", basename(output_file), 
              nrow(all_stats)))
  
  invisible(output_file)
}

print_progress_header <- function(current, total) {
  pct <- sprintf("%.1f%%", (current/total)*100)
  cat(sprintf("\n%s\n进度: [%d/%d] %s\n%s\n",
              strrep("=", 60), current, total, pct, 
              strrep("=", 60)))
}

print_progress_header_multi <- function(
    current, total, genelist_name, seurat_basename) {
  
  pct <- sprintf("%.1f%%", (current/total)*100)
  cat(sprintf("\n%s\n任务: [%d/%d] %s\n基因列表: %s\nSeurat: %s\n%s\n",
              strrep("=", 60), current, total, pct, 
              genelist_name, seurat_basename, strrep("=", 60)))
}

estimate_remaining_time <- function(
    results, seurat_files, current_idx) {
  
  avg_time <- mean(sapply(results, function(x) {
    x$processing_time %||% 0
  }), na.rm = TRUE)
  
  remaining <- avg_time * (length(seurat_files) - current_idx)
  cat(sprintf("\n📊 预计剩余: %.2f 分钟 (%.2f 小时)\n", 
              remaining, remaining/60))
}

estimate_remaining_time_multi <- function(
    results, total_tasks, completed) {
  
  avg_time <- mean(sapply(results, function(x) {
    x$processing_time %||% 0
  }), na.rm = TRUE)
  
  remaining <- avg_time * (total_tasks - completed)
  cat(sprintf("\n📊 预计剩余: %.2f 分钟 (%.2f 小时)\n", 
              remaining, remaining/60))
}

confirm_batch_processing <- function(seurat_files, CONFIG) {
  if (!CONFIG$batch_mode) return(TRUE)
  if (length(seurat_files) <= 1) return(TRUE)
  if (!interactive()) return(TRUE)
  
  response <- readline(sprintf(
    "即将处理 %d 个文件，是否继续? (y/n): ", 
    length(seurat_files)))
  
  return(tolower(response) == "y")
}

create_summary_object <- function(
    results, total_elapsed, log_files) {
  
  list(
    total = length(results),
    success = sum(sapply(results, function(x) x$success)),
    failed = sum(sapply(results, function(x) !x$success)),
    total_time = as.numeric(total_elapsed),
    log_file = log_files$log,
    csv_file = log_files$csv
  )
}

cat("✅ 10_batch_processing.R 已加载\n")