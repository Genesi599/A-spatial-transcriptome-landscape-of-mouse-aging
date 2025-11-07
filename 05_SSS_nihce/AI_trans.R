
#!/usr/bin/env Rscript
# 项目信息导出脚本
# 将当前目录的项目结构和内容导出为Markdown文件

# 加载必要的包
if (!require("tools", quietly = TRUE)) install.packages("tools")

# 配置参数
output_file <- "project_export.md"
max_file_size <- 1024 * 1024  # 1MB, 超过此大小的文件不读取内容
ignore_dirs <- c(".git", ".Rproj.user", "node_modules", "__pycache__", ".venv", "venv")
ignore_files <- c(".DS_Store", "Thumbs.db", ".gitignore")
code_extensions <- c("R", "r", "py", "js", "html", "css", "java", "cpp", "c", "h", 
                     "sh", "sql", "md", "rmd", "yml", "yaml", "json", "xml", "txt")

# 函数：获取文件大小的可读格式
format_file_size <- function(size) {
  if (size < 1024) {
    return(paste0(size, " B"))
  } else if (size < 1024^2) {
    return(paste0(round(size / 1024, 2), " KB"))
  } else if (size < 1024^3) {
    return(paste0(round(size / 1024^2, 2), " MB"))
  } else {
    return(paste0(round(size / 1024^3, 2), " GB"))
  }
}

# 函数：生成目录树
generate_tree <- function(path, prefix = "", is_last = TRUE) {
  tree_lines <- c()
  
  files <- list.files(path, all.files = FALSE, include.dirs = TRUE, no.. = TRUE)
  files <- files[!files %in% ignore_dirs & !files %in% ignore_files]
  
  if (length(files) == 0) return(tree_lines)
  
  for (i in seq_along(files)) {
    file_path <- file.path(path, files[i])
    is_last_item <- (i == length(files))
    
    connector <- if (is_last_item) "└── " else "├── "
    tree_lines <- c(tree_lines, paste0(prefix, connector, files[i]))
    
    if (dir.exists(file_path)) {
      new_prefix <- paste0(prefix, if (is_last_item) "    " else "│   ")
      tree_lines <- c(tree_lines, generate_tree(file_path, new_prefix, is_last_item))
    }
  }
  
  return(tree_lines)
}

# 函数：递归获取所有文件
get_all_files <- function(path) {
  all_files <- c()
  
  files <- list.files(path, all.files = FALSE, include.dirs = TRUE, 
                     full.names = TRUE, no.. = TRUE)
  
  for (file_path in files) {
    file_name <- basename(file_path)
    
    # 跳过忽略的目录和文件
    if (file_name %in% ignore_dirs || file_name %in% ignore_files) next
    
    if (dir.exists(file_path)) {
      all_files <- c(all_files, get_all_files(file_path))
    } else {
      all_files <- c(all_files, file_path)
    }
  }
  
  return(all_files)
}

# 函数：判断是否为文本文件
is_text_file <- function(file_path) {
  ext <- tools::file_ext(file_path)
  return(tolower(ext) %in% tolower(code_extensions))
}

# 主函数
export_project <- function() {
  cat("开始导出项目信息...\n")
  
  # 创建Markdown内容
  md_content <- c()
  
  # 标题和基本信息
  project_name <- basename(getwd())
  md_content <- c(md_content, paste0("# 项目导出: ", project_name))
  md_content <- c(md_content, "")
  md_content <- c(md_content, paste0("**导出时间**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  md_content <- c(md_content, paste0("**项目路径**: ", getwd()))
  md_content <- c(md_content, "")
  md_content <- c(md_content, "---")
  md_content <- c(md_content, "")
  
  # 目录结构
  cat("生成目录树...\n")
  md_content <- c(md_content, "## 目录结构")
  md_content <- c(md_content, "")
  md_content <- c(md_content, "```")
  md_content <- c(md_content, project_name)
  tree_lines <- generate_tree(getwd())
  md_content <- c(md_content, tree_lines)
  md_content <- c(md_content, "```")
  md_content <- c(md_content, "")
  
  # 获取所有文件
  cat("收集文件信息...\n")
  all_files <- get_all_files(getwd())
  
  # 文件统计
  md_content <- c(md_content, "## 文件统计")
  md_content <- c(md_content, "")
  md_content <- c(md_content, paste0("- **总文件数**: ", length(all_files)))
  
  # 按扩展名统计
  extensions <- sapply(all_files, tools::file_ext)
  ext_table <- table(extensions)
  ext_table <- sort(ext_table, decreasing = TRUE)
  
  md_content <- c(md_content, "- **文件类型分布**:")
  for (i in seq_along(ext_table)) {
    ext_name <- if (names(ext_table)[i] == "") "(无扩展名)" else paste0(".", names(ext_table)[i])
    md_content <- c(md_content, paste0("  - ", ext_name, ": ", ext_table[i], " 个"))
  }
  
  # 总大小
  total_size <- sum(sapply(all_files, file.size))
  md_content <- c(md_content, paste0("- **项目总大小**: ", format_file_size(total_size)))
  md_content <- c(md_content, "")
  
  # 文件详细内容
  md_content <- c(md_content, "## 文件内容")
  md_content <- c(md_content, "")
  
  cat("读取文件内容...\n")
  for (file_path in all_files) {
    relative_path <- sub(paste0("^", getwd(), "/"), "", file_path)
    file_size <- file.size(file_path)
    
    md_content <- c(md_content, paste0("### ", relative_path))
    md_content <- c(md_content, "")
    md_content <- c(md_content, paste0("- **大小**: ", format_file_size(file_size)))
    md_content <- c(md_content, paste0("- **修改时间**: ", format(file.mtime(file_path), "%Y-%m-%d %H:%M:%S")))
    md_content <- c(md_content, "")
    
    # 如果是文本文件且大小合适，读取内容
    if (is_text_file(file_path) && file_size <= max_file_size) {
      tryCatch({
        content <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
        ext <- tools::file_ext(file_path)
        lang <- if (ext == "") "" else tolower(ext)
        
        md_content <- c(md_content, paste0("```", lang))
        md_content <- c(md_content, content)
        md_content <- c(md_content, "```")
        md_content <- c(md_content, "")
      }, error = function(e) {
        md_content <<- c(md_content, "*无法读取文件内容*")
        md_content <<- c(md_content, "")
      })
    } else if (file_size > max_file_size) {
      md_content <- c(md_content, "*文件过大，跳过内容显示*")
      md_content <- c(md_content, "")
    } else {
      md_content <- c(md_content, "*二进制文件，跳过内容显示*")
      md_content <- c(md_content, "")
    }
    
    md_content <- c(md_content, "---")
    md_content <- c(md_content, "")
  }
  
  # 写入文件
  cat(paste0("写入文件: ", output_file, "\n"))
  writeLines(md_content, output_file, useBytes = TRUE)
  
  cat(paste0("✓ 导出完成！文件已保存到: ", output_file, "\n"))
  cat(paste0("  共处理 ", length(all_files), " 个文件\n"))
}

# 执行导出
tryCatch({
  export_project()
}, error = function(e) {
  cat("错误:", conditionMessage(e), "\n")
})
