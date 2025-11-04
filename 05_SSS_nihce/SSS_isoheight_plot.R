# SSS_isoheight_plot.R (修复版)

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
            rownames_to_column(var = "cellid")
        })
}

#' high light spot isoheight plot
#' 
#' @import ggnewscale
#' 
#' @param .data seurat obj
#' @param density_top expressions use for filter to select the high light spot
#' @param col_bg background spot color
#' @param col_top high light spot color
#' @param col_isoheight isoheight line color
#' @param col_white_ratio white ratio in isoheight fill colors
#' @param cols_fill_isoheight isoheight fill colors
#' @param size_bg background spot size
#' @param size_top high light spot size
#' @param nrow number of rows use for facet_wrap
#' 
#' @return ggplot obj
#' 
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
    nrow = 2
) {

    density_top  <- enquo(density_top)

    df <- .data@meta.data %>%
        rownames_to_column("cellid") %>%
        inner_join(
            GetAllCoordinates(.data)
        ) %>%
        as_tibble()

    # ✅ 修复：使用更宽松的坐标范围（增加 5% 边距）
    xy_r_max <- max(df$col, df$row)
    xy_r_min <- min(df$col, df$row)
    margin <- (xy_r_max - xy_r_min) * 0.05
    
    xlim_range <- c(max(0, min(df$col) - margin), max(df$col) + margin)
    ylim_range <- c(max(0, min(df$row) - margin), max(df$row) + margin)

    p <- ggplot(mapping = aes(x = col, y = row)) +
        # 背景点
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        )

    p <- p +
        ggnewscale::new_scale_fill() +
        # 密度填充（等高线填充）
        stat_density_2d_filled(
            data = df %>% filter(!!density_top),
            mapping = aes(fill = after_stat(ndensity), alpha = after_stat(ndensity)),
            geom = "raster", contour = FALSE,
            show.legend = TRUE  # ✅ 显示填充图例
        ) +
        scale_fill_gradientn(
            colours = cols_fill_isoheight,
            name = "Density"  # ✅ 统一图例标题
        ) +
        # ✅ 修复：移除 alpha 的图例
        guides(alpha = "none") +
        
        ggnewscale::new_scale_fill() +
        # 等高线线条
        geom_density_2d(
            data = df %>% filter(!!density_top),
            color = col_isoheight,
            contour_var = "ndensity",
            show.legend = FALSE  # ✅ 关闭线条图例（避免重复）
        )

    # 高亮点
    p <- p +
        geom_point(
            data = df %>% filter(!!density_top),
            color = col_top, alpha = 0.5, size = size_top
        )

    # 坐标和主题
    p <- p +
        scale_y_reverse() +
        NoAxes() +
        coord_fixed(
            ratio = 1, 
            xlim = xlim_range,  # ✅ 使用动态范围
            ylim = ylim_range,
            expand = FALSE
        ) +
        theme(
            aspect.ratio = 1,
            panel.background = element_blank(),
            strip.background = element_blank(),
            legend.position = "right",  # ✅ 图例位置
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10)
        ) +
        facet_wrap(vars(sample), nrow = nrow)
    
    return(p)
}