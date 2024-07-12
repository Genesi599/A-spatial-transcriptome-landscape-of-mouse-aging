

celltype_isoheight_plot <- function(data,
    density_head,
    density_tail,

    col_bg = "gray90",
    col_cell_head = "darkred",
    col_cell_tail = "black",

    head_isoheight_col = "white",
    h_p_w = 0.2,
    head_cols = c(
        rep("white", round(100 * h_p_w)),
        colorRampPalette(brewer.pal(5, "YlOrRd")[2:5])(round(100 * (1 - h_p_w)))
    ),
    t_p_w = h_p_w,
    tail_cols = c(
        rep("white", round(100 * t_p_w)),
        colorRampPalette(brewer.pal(5, "GnBu")[2:5])(round(100 * (1 - t_p_w)))
    ),
    name_list = NULL,

    size_bg = 0.1,
    size_cell = size_bg,

    nrow = 2) {

    density_head  <- enquo(density_head)

    if (!missing(density_tail)) {
        density_tail  <- enquo(density_tail)
    }

    df <- data@meta.data %>%
        rownames_to_column("cellid") %>%
        inner_join(
            GetAllCoordinates(data)
        ) %>%
        as_tibble()

    if (!is.null(name_list)) {
        obj_rename_list <- name_list$id %>%
            setNames(name_list$sample)

        df <- df %>%
            mutate(
                sample = obj_rename_list[orig.ident],
                sample = factor(sample, levels = obj_rename_list)
            )
    } else {
        df <- df %>%
            arrange(age, sample) %>%
            mutate(sample = factor(sample, levels = unique(sample)))
    }

    xy_r <- max(df$col, df$row)

    p <- ggplot(mapping = aes(x = col, y = row)) +
        geom_point(
            data = df,
            color = col_bg, alpha = 0.5, size = size_bg
        )

    if (!missing(density_tail)) {
        p <- p +
            stat_density_2d_filled(
                data = df %>% filter(!!density_tail),
                mapping = aes(fill = ..ndensity.., alpha = ..ndensity.. ),
                geom = "raster", contour = F
            ) +
            scale_fill_gradientn(colours = tail_cols)
    }

    p <- p +
        ggnewscale::new_scale_fill() +
        stat_density_2d_filled(
            data = df %>% filter(!!density_head),
            mapping = aes(fill = ..ndensity.., alpha = ..ndensity.. ),
            geom = "raster", contour = F
        ) +
        scale_fill_gradientn(colours = head_cols) +
        ggnewscale::new_scale_fill() +
        geom_density_2d(
            data = df %>% filter(!!density_head),
            color = head_isoheight_col,
            contour_var	= "ndensity",
            show.legend = T
        )

    if (!missing(density_tail)) { 
        p <- p +
            geom_point(
                data = df %>% filter(!!density_tail),
                color = col_cell_tail, alpha = 0.5, size = size_cell
            )
    }

    p <- p +
        geom_point(
            data = df %>% filter(!!density_head),
            color = col_cell_head, alpha = 0.5, size = size_cell
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
