# SSS_isoheight_plot.R (最彻底的修复)

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
        inner_join(GetAllCoordinates(.data)) %>%
        as_tibble()

    # ✅ 完全不设置坐标限制
    p <- ggplot(mapping = aes(x = col, y = row)) +
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        ) +
        ggnewscale::new_scale_fill() +
        stat_density_2d_filled(
            data = df %>% filter(!!density_top),
            mapping = aes(fill = after_stat(ndensity)),
            geom = "raster", 
            contour = FALSE,
            alpha = 0.8,
            show.legend = TRUE
        ) +
        scale_fill_gradientn(
            colours = cols_fill_isoheight,
            name = "Density"
        ) +
        guides(alpha = "none") +
        ggnewscale::new_scale_fill() +
        geom_density_2d(
            data = df %>% filter(!!density_top),
            color = col_isoheight,
            contour_var = "ndensity",
            show.legend = FALSE
        ) +
        geom_point(
            data = df %>% filter(!!density_top),
            color = col_top, alpha = 0.5, size = size_top
        ) +
        scale_y_reverse() +  # ✅ 只保留 y 轴翻转
        coord_fixed(ratio = 1) +  # ✅ 只保持宽高比
        NoAxes() +
        theme(
            aspect.ratio = 1,
            panel.background = element_blank(),
            strip.background = element_blank(),
            legend.position = "right",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10)
        ) +
        facet_wrap(vars(sample), nrow = nrow)
    
    return(p)
}