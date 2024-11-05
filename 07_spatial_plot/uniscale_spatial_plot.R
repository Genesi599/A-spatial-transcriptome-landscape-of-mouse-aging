

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

uniscale_SpatialDimPlot <- function(
    .data,
    group.by = NULL,
    cols = NULL,
    images = names(.data@images),
    ncol,
    nrow = ceiling(length(images) / ncol),
    plot_theme = theme(),
    ...
    ) {
    data_co <- GetAllCoordinates(.data)

    xy_r <- max(data_co$col, data_co$row)

    r_r <- .data@images %>%
        map_dbl(Radius) %>%
        min()
    .data@images <- .data@images %>%
        map( ~{
            .x@spot.radius <- r_r
            .x
        })

    plot_list <- SpatialDimPlot(.data,
            group.by = group.by,
            images = images,
            pt.size.factor = 1.2,
            combine = F,
            image.alpha = 0,
            cols = cols,
            ...
        ) %>%
        map(~ .x + coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r)))

    plot_list <- plot_list %>%
        map(~ .x + plot_theme)

    print(length(plot_list))
    ggpubr::ggarrange(
        plotlist = plot_list,
        ncol = ncol,
        nrow = nrow,
        common.legend = T,
        legend = "right"
    )
}


uniscale_SpatialFeaturePlot <- function(
    .data,
    feature,
    cols = SpatialColors(100),
    images = names(.data@images),
    ncol,
    nrow = ceiling(length(images) / ncol),
    limit = TRUE,
    plot_theme = theme(),
    min_vlaue = NA
    ) {
    data_co <- GetAllCoordinates(.data)

    xy_r <- max(data_co$col, data_co$row)

    r_r <- .data@images %>%
        map_dbl(Radius) %>%
        min()
    .data@images <- .data@images %>%
        map( ~{
            .x@spot.radius <- r_r
            .x
        })

    plot_list <- SpatialFeaturePlot(.data,
            features = feature,
            images = images,
            pt.size.factor = 1.2,
            combine = F,
            image.alpha = 0
        ) %>%
        map(~ .x + coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r)))

    if (limit) {
        data <- FetchData(
            object = .data,
            vars = feature
        )
        fill_min <- min(data[[feature]])
        fill_max <- max(data[[feature]])

        if(!is.na(min_vlaue)) {
            fill_min <- min_vlaue
        } 
        plot_list <- plot_list %>%
            map(~ .x + scale_fill_gradientn(colours = cols, limits = c(fill_min, fill_max), na.value = cols[1]))

    } else {
        plot_list <- plot_list %>%
            map(~ .x + scale_fill_gradientn(colours = cols))
    }

    plot_list <- plot_list %>%
        map(~ .x + plot_theme)

    print(length(plot_list))
    ggpubr::ggarrange(
        plotlist = plot_list,
        ncol = ncol,
        nrow = nrow,
        common.legend = limit,
        legend = "right"
    )
}

