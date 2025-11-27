# ==================================================
# 11_sample_process.R
# 单样本处理函数
# ==================================================

pick_coord_cols <- function(md, coords = NULL) {
  cands <- list(
    list(col = "col",                row = "row"),
    list(col = "x",                  row = "y"),
    list(col = "imagecol",           row = "imagerow"),
    list(col = "pxl_col_in_fullres", row = "pxl_row_in_fullres"),
    list(col = "array_col",          row = "array_row")
  )
  
  if (!is.null(coords)) {
    if (all(c("x","y") %in% colnames(coords))) 
      return(list(col = "x", row = "y", src = "coords"))
    if (all(c("imagecol","imagerow") %in% colnames(coords))) 
      return(list(col = "imagecol", row = "imagerow", 
                  src = "coords"))
  }
  
  for (c in cands) {
    if (all(c(c$col, c$row) %in% colnames(md))) 
      return(list(col = c$col, row = c$row, src = "meta"))
  }
  
  return(NULL)
}

get_df_std <- function(x) {
  if (inherits(x, "Seurat")) {
    seu <- x
    md <- seu@meta.data
    coords <- NULL
    
    if (requireNamespace("Seurat", quietly = TRUE)) {
      coords <- tryCatch(
        Seurat::GetTissueCoordinates(seu), 
        error = function(e) NULL
      )
      
      if (!is.null(coords) && is.null(rownames(coords)) && 
          "barcode" %in% colnames(coords)) {
        rownames(coords) <- coords$barcode
      }
      
      if (!is.null(coords)) {
        common <- intersect(rownames(md), rownames(coords))
        md <- md[common, , drop = FALSE]
        coords <- coords[common, , drop = FALSE]
      }
    }
    
    info <- pick_coord_cols(md, coords)
    if (is.null(info)) {
      stop(paste0(
        "无法在 Seurat 对象中识别坐标列，",
        "请检查 meta.data 或 GetTissueCoordinates 结果"
      ))
    }
    
    if (identical(info$src, "coords")) {
      out <- cbind(
        md, 
        coords[, c(info$col, info$row), drop = FALSE]
      )
      colnames(out)[(ncol(out)-1):ncol(out)] <- c("col", "row")
    } else {
      out <- md
      colnames(out)[match(
        c(info$col, info$row), 
        colnames(out)
      )] <- c("col", "row")
    }
    
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    return(out)
    
  } else if (is.data.frame(x)) {
    out <- x
    
    if (!all(c("col","row") %in% colnames(out))) {
      info <- pick_coord_cols(out, NULL)
      if (is.null(info)) {
        stop(paste0(
          "数据中缺少坐标列（col/row 或常见别名 ",
          "x/y, imagecol/imagerow, ",
          "pxl_col_in_fullres/pxl_row_in_fullres）"
        ))
      }
      colnames(out)[match(
        c(info$col, info$row), 
        colnames(out)
      )] <- c("col", "row")
    }
    return(out)
    
  } else {
    stop("df 必须是 Seurat 或 data.frame")
  }
}

qout <- function(CONFIG, ...) {
  if (!isTRUE(CONFIG$quiet)) cat(sprintf(...))
}

process_single_sample <- function(df, sample_id, CONFIG) {
  qout(CONFIG, "\n[%s]\n", sample_id)
  
  cached_data <- load_plot_cache(sample_id, CONFIG)
  if (!is.null(cached_data)) {
    qout(CONFIG, "      🎨 使用缓存数据...\n")
    
    p_overlay <- plot_celltype_density_overlay(
      cached_data$df, 
      cached_data$density_data, 
      sample_id, 
      CONFIG
    )
    p_composition <- plot_zone_composition(
      cached_data$zone_composition, 
      sample_id, 
      CONFIG
    )
    
    overlay_png <- file.path(
      CONFIG$output$plot_dir, 
      sprintf("%s_overlay.png", sample_id)
    )
    overlay_pdf <- file.path(
      CONFIG$output$plot_dir, 
      sprintf("%s_overlay.pdf", sample_id)
    )
    composition_png <- file.path(
      CONFIG$output$plot_dir, 
      sprintf("%s_composition.png", sample_id)
    )
    composition_pdf <- file.path(
      CONFIG$output$plot_dir, 
      sprintf("%s_composition.pdf", sample_id)
    )
    
    ggsave(overlay_png, p_overlay, 
           width = 16, height = 12, dpi = 300, bg = "white")
    ggsave(overlay_pdf, p_overlay, 
           width = 16, height = 12, bg = "white")
    ggsave(composition_png, p_composition, 
           width = 14, height = 10, dpi = 300, bg = "white")
    ggsave(composition_pdf, p_composition, 
           width = 14, height = 10, bg = "white")
    
    n_spots <- nrow(cached_data$df)
    n_high  <- sum(!is.na(cached_data$df$density_zone))
    n_types <- length(setdiff(
      unique(cached_data$df$celltype_clean), 
      "Unknown"
    ))
    
    qout(CONFIG, "  ✅ %d spots | %d high | %d types (缓存)\n", 
         n_spots, n_high, n_types)
    
    return(list(
      density_data = cached_data$density_data,
      zone_composition = cached_data$zone_composition,
      plots = list(
        overlay = p_overlay, 
        composition = p_composition
      ),
      stats = list(
        n_spots = n_spots, 
        n_high_density = n_high, 
        n_celltypes = n_types
      ),
      from_cache = TRUE
    ))
  }
  
  if (is.null(CONFIG$colors$celltype)) {
    stop("❌ 全局颜色方案未初始化！")
  }
  
  df <- get_df_std(df)
  
  if (is.null(CONFIG$params$celltype_col) || 
      !(CONFIG$params$celltype_col %in% colnames(df))) {
    stop(sprintf(
      "❌ 细胞类型列 '%s' 不存在", 
      CONFIG$params$celltype_col
    ))
  }
  
  raw_celltypes <- df[[CONFIG$params$celltype_col]]
  if (is.null(raw_celltypes) || length(raw_celltypes) == 0) {
    stop(sprintf(
      "❌ 细胞类型列 '%s' 为空", 
      CONFIG$params$celltype_col
    ))
  }
  
  df$celltype_clean <- ifelse(
    is.na(raw_celltypes) | raw_celltypes == "", 
    "Unknown", 
    as.character(raw_celltypes)
  )
  
  unique_types <- sort(unique(
    df$celltype_clean[df$celltype_clean != "Unknown"]
  ))
  n_types <- length(unique_types)
  qout(CONFIG, "  📊 细胞类型: %d 个\n", n_types)
  
  if (n_types > 0 && !isTRUE(CONFIG$quiet)) {
    head_n <- min(10, n_types)
    for (i in seq_len(head_n)) 
      qout(CONFIG, "     • %s\n", unique_types[i])
    if (n_types > head_n) 
      qout(CONFIG, "     ... 还有 %d 个\n", n_types - head_n)
  }
  
  all_types_global <- names(CONFIG$colors$celltype)
  sample_types <- setdiff(unique(df$celltype_clean), "Unknown")
  missing_types <- setdiff(sample_types, all_types_global)
  if (length(missing_types) > 0) {
    warning(sprintf(
      "  ⚠️  未知类型: %s", 
      paste(missing_types, collapse = ", ")
    ))
  }
  
  expand_margin <- if (!is.null(CONFIG$params$expand_margin)) 
    CONFIG$params$expand_margin else 0.1
  
  if (!("ClockGene_High" %in% colnames(df))) {
    stop("缺少列 'ClockGene_High'，无法进行密度计算")
  }
  
  density_data <- calculate_density_zones(
    df = df,
    density_bins  = CONFIG$params$n_zones,
    expand_margin = expand_margin,
    quiet = TRUE
  )
  
  if (is.null(density_data)) {
    warning(sprintf(
      "[%s] 密度计算失败或高表达点不足，跳过该样本", 
      sample_id
    ))
    return(NULL)
  }
  
  df <- df %>%
    dplyr::left_join(
      density_data$spot_zones,
      by = c("col", "row")
    )
  
  zone_composition <- df %>%
    dplyr::filter(!is.na(density_zone)) %>%
    dplyr::group_by(density_zone, celltype_clean) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(density_zone) %>%
    dplyr::mutate(
      total = sum(count),
      percentage = (count / total) * 100
    ) %>%
    dplyr::ungroup()
  
  if (isTRUE(CONFIG$debug_mode)) {
    plot_data <- list(
      df = df,
      density_data = density_data,
      zone_composition = zone_composition,
      params = list(
        sample_id = sample_id,
        n_spots = nrow(df),
        n_zones = CONFIG$params$n_zones,
        cache_time = Sys.time()
      )
    )
    save_plot_cache(sample_id, plot_data, CONFIG)
  }
  
  p_overlay <- plot_celltype_density_overlay(
    df, density_data, sample_id, CONFIG
  )
  p_composition <- plot_zone_composition(
    zone_composition, sample_id, CONFIG
  )
  
  # 定义文件路径
  overlay_file_png <- file.path(
    CONFIG$output$plot_dir, 
    sprintf("%s_overlay.png", sample_id)
  )
  overlay_file_pdf <- file.path(
    CONFIG$output$plot_dir, 
    sprintf("%s_overlay.pdf", sample_id)
  )

  composition_file_png <- file.path(
    CONFIG$output$plot_dir, 
    sprintf("%s_composition.png", sample_id)
  )
  composition_file_pdf <- file.path(
    CONFIG$output$plot_dir, 
    sprintf("%s_composition.pdf", sample_id)
  )

  # 保存 PNG 版本
  ggsave(overlay_file_png, p_overlay, 
        width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(composition_file_png, p_composition, 
        width = 14, height = 10, dpi = 300, bg = "white")

  # 保存 PDF 版本
  ggsave(overlay_file_pdf, p_overlay, 
        width = 16, height = 12, bg = "white")
  ggsave(composition_file_pdf, p_composition, 
        width = 14, height = 10, bg = "white")
  
  zone_comp_file <- file.path(
    CONFIG$output$data_dir, 
    sprintf("%s_zone_composition.csv", sample_id)
  )
  utils::write.csv(zone_composition, zone_comp_file, 
                   row.names = FALSE)
  
  n_spots <- nrow(df)
  n_high  <- sum(!is.na(df$density_zone))
  n_types <- length(setdiff(
    unique(df$celltype_clean), 
    "Unknown"
  ))
  qout(CONFIG, "  ✅ %d spots | %d high | %d types\n", 
       n_spots, n_high, n_types)
  
  return(list(
    density_data = density_data,
    zone_composition = zone_composition,
    plots = list(
      overlay = p_overlay, 
      composition = p_composition
    ),
    stats = list(
      n_spots = n_spots, 
      n_high_density = n_high, 
      n_celltypes = n_types
    ),
    from_cache = FALSE
  ))
}