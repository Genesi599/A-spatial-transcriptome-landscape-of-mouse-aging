
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

    xy_r <- max(df$col, df$row)

    p <- ggplot(mapping = aes(x = col, y = row)) +
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        )

    p <- p +
        ggnewscale::new_scale_fill() +
        stat_density_2d_filled(
            data = df %>% filter(!!density_top),
            mapping = aes(fill = ..ndensity.., alpha = ..ndensity.. ),
            geom = "raster", contour = F
        ) +
        scale_fill_gradientn(colours = cols_fill_isoheight) +
        ggnewscale::new_scale_fill() +
        geom_density_2d(
            data = df %>% filter(!!density_top),
            color = col_isoheight,
            contour_var	= "ndensity",
            show.legend = T
        )

    p <- p +
        geom_point(
            data = df %>% filter(!!density_top),
            color = col_top, alpha = 0.5, size = size_top
        )

    p <- p +
        scale_y_reverse() +
        NoAxes() +
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r)) +
        theme(
            aspect.ratio = 1,
            panel.background = element_blank(),
            strip.background = element_blank()
        ) +
        facet_wrap(vars(sample), nrow = nrow)
}
