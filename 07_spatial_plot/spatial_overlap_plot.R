
GetAllCoordinates <- function(.data) {
    .data@images %>%
        names() %>%
        unique() %>%
        map_dfr(~{
            GetTissueCoordinates(
                    .data,
                    image = gsub('-', '.', .x),
                    cols = c("row", "col"),
                    scale = NULL
                ) %>%
            rownames_to_column(var = "cellid")
        })
}

library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))
PairedColors <- colorRampPalette(colors = (brewer.pal(n = 11, name = "Paired")))

###############################################################

spatial_overlap_plot <- function(data,
    cell_cols,

    col_bg = "#050137",
    col_bg2 = "black",

    size_bg = 0.1,
    size_cell = size_bg,

    nrow = 2) {

    df <- data@meta.data %>%
        rownames_to_column("cellid") %>%
        inner_join(
            GetAllCoordinates(data)
        ) %>%
        as_tibble()

    df <- df %>%
        arrange(age, sample) %>%
        mutate(sample = factor(sample, levels = unique(sample)))

    xy_r <- max(df$col, df$row)

    df$count <- 0
    for(i in names(cell_cols)) {

        df <- df %>%
            mutate(
                count = !!sym(i) + count
            )
    }

    p <- ggplot(mapping = aes(x = col, y = row)) +
        geom_point(
            data = df %>%
                filter(count == 0),
            color = col_bg, 
            alpha = 0.5, 
            size = size_bg
        )

    p <- p +
        geom_point(
            data = df %>%
                filter(count != 0),
            color = col_bg2,
            alpha = 0.5,
            size = size_bg
        )

    point_alpha <- 1 / length(cell_cols)

    print(cell_cols)
    for(i in names(cell_cols)) {
        p <- p +
            geom_point(
                data = df %>% filter(!!sym(i)),
                color = cell_cols[i],
                size = size_cell,
                alpha = point_alpha
            )
    }
    p <- p +
        scale_y_reverse() +
        NoLegend() +
        Seurat::DarkTheme() +
        theme(
            strip.background = element_blank()
        ) +
        NoAxes() +
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r)) +
        facet_wrap(vars(sample), nrow = nrow)
}
