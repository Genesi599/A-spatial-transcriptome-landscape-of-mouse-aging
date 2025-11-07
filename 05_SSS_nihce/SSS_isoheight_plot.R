# SSS_isoheight_plot.R (完整版)

GetAllCoordinates <- function(.data) {
    .data@images %>%
        names() %>%
        unique() %>%
        map_dfr(~{
            GetTissueCoordinates(
                    .data,
                    image = .x,
                    cols = c("row", "col"),
                    scale = NULL
                ) %>%
            tibble::rownames_to_column(var = "cellid")
        })
}

celltype_isoheight_plot <- function(
    .data,
    density_top,
    col_bg = "gray90",
    col_top = "darkred",
    col_isoheight = "white",
    col_white_ratio = 0.2,
    cols_fill_isoheight = c(
        rep("white", round(100 * col_white_ratio)),
        colorRampPalette(brewer.pal(5, "YlOrRd")[2:5])(round(100 * (1 - col_white_ratio)))
    ),
    size_bg = 0.1,
    size_top = size_bg,
    nrow = 2,
    expand_margin = 0.05  # ✅ 新增参数：边缘扩展比例（5%）
) {

    density_top  <- enquo(density_top)

    df <- .data@meta.data %>%
        tibble::rownames_to_column("cellid") %>%
        dplyr::inner_join(GetAllCoordinates(.data)) %>%
        as_tibble()

    # ✅ 计算坐标范围并扩展
    col_range <- range(df$col, na.rm = TRUE)
    row_range <- range(df$row, na.rm = TRUE)
    
    col_expand <- diff(col_range) * expand_margin
    row_expand <- diff(row_range) * expand_margin
    
    col_limits <- c(col_range[1] - col_expand, col_range[2] + col_expand)
    row_limits <- c(row_range[1] - row_expand, row_range[2] + row_expand)
    
    cat(sprintf("✅ 坐标范围: col [%.1f, %.1f], row [%.1f, %.1f]\n",
                col_limits[1], col_limits[2], row_limits[1], row_limits[2]))

    # ✅ 绘图
    p <- ggplot(mapping = aes(x = col, y = row)) +
        # 1. 背景点
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        ) +
        ggnewscale::new_scale_fill() +
        
        # 2. 密度热图 (关键：设置 bins 和 expand)
        stat_density_2d_filled(
            data = df %>% dplyr::filter(!!density_top),
            mapping = aes(fill = after_stat(ndensity)),
            geom = "raster", 
            contour = FALSE,
            alpha = 0.8,
            bins = 100,  # ✅ 增加分辨率
            show.legend = TRUE,
            n = 200      # ✅ 增加网格密度
        ) +
        scale_fill_gradientn(
            colours = cols_fill_isoheight,
            name = "Density"
        ) +
        guides(alpha = "none") +
        ggnewscale::new_scale_fill() +
        
        # 3. 等高线
        geom_density_2d(
            data = df %>% dplyr::filter(!!density_top),
            color = col_isoheight,
            contour_var = "ndensity",
            bins = 10,    # ✅ 等高线数量
            show.legend = FALSE
        ) +
        
        # 4. 高亮点
        geom_point(
            data = df %>% dplyr::filter(!!density_top),
            color = col_top, alpha = 0.5, size = size_top
        ) +
        
        # ✅ 关键：手动设置坐标范围
        scale_x_continuous(
            limits = col_limits,
            expand = expansion(mult = 0.02)  # 额外 2% 边距
        ) +
        scale_y_reverse(
            limits = rev(row_limits),  # ✅ 注意反转顺序
            expand = expansion(mult = 0.02)
        ) +
        
        # 坐标系统
        coord_fixed(ratio = 1) +  # ✅ 保持宽高比
        
        # 主题
        NoAxes() +
        theme(
            aspect.ratio = 1,
            panel.background = element_blank(),
            plot.background = element_blank(),  # ✅ 移除绘图背景
            strip.background = element_blank(),
            legend.position = "right",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.margin = margin(5, 5, 5, 5)  # ✅ 统一边距
        ) +
        facet_wrap(vars(sample), nrow = nrow)
    
    return(p)
}
