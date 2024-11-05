
# ################################################################
#
# plot
#
# ################################################################

#' 
#' @param centers tibble obj
#' @param n_work number of threads
#' 
segment_neighbor_check <- function(centers, n_work) {
    library(future)
    library(future.apply)

    plan(multisession, workers=n_work)
    options(future.globals.maxSize= Inf)
    message(">> how many cores can use now: ", nbrOfWorkers())

    centers %>%
        filter(!is.na(label)) %>%
        select(orig.ident, sample, col, row, label) %>%
        group_by(orig.ident, sample) %>%
        group_nest() %>%
        mutate(
            data = future_lapply(data, function(df) {

                mat <- Matrix::sparseMatrix(
                    df$col,
                    df$row,
                    x = as.numeric(factor(df$label)))

                df %>%
                    pmap_dfr(~{
                        x <- ..1
                        y <- ..2

                        res <- tibble()
                        if (x + 1 <= dim(mat)[1] && mat[x, y] == mat[x + 1,y]) {
                            res <- res %>%
                                bind_rows(tibble(x1 = x, y1 = y, x2 = x + 1, y2 = y))
                        }
                        if (y + 1 <= dim(mat)[2] && mat[x, y] == mat[x, y + 1]) {
                            res <- res %>%
                                bind_rows(tibble(x1 = x, y1 = y, x2 = x, y2 = y + 1))
                        }
                        res
                    }) %>%
                    distinct()

            }, future.chunk.size = Inf)
        ) %>%
        unnest(data)
}

#' plot for tissue_SLIC_entropy reult
#' 
#' @import ggforce
#' @import ggnewscale
#' 
#' @param centers 
#' @param entropy
#' 
#' @param n_work number of threads
#' @param celltype celltype_col
#' @param segment_neighbor spot neighbor information, result of segment_neighbor_check
#' 
#' @param spot_size background spot size
#' @param spot_color
#' @param spot_fill
#' @param spot_guides_hide
#' @param spot_stroke
#' 
#' @param segment_alpha segment connect spots in a superpixel
#' @param segment_color
#' 
#' @param show_entropy entropy point in the center of superpixel
#' @param point_shape
#' @param point_stroke
#' @param point_color
#' @param point_factor
#' @param point_fill_cols
#' 
#' @param show_pie pie plot showing the ratio of spot in a superpixel
#' @param pie_alpha
#' @param pie_color
#' @param pie_stroke
#' 
#' @param samples
#' @param legend_position
#' @param nrow
#' 
#' @return ggplot obj
#' 
SLIC_segment_plot <- function(
    centers, 
    entropy,
    
    n_work = 3,
    segment_neighbor = segment_neighbor_check(centers, n_work),
    celltype = seurat_clusters,

    spot_size = 1,
    spot_color = "black",
    spot_fill = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(centers[[as_label(celltype)]]))),
    spot_stroke = 0.1,
    spot_guides_hide = TRUE,

    segment_alpha = 0.5,
    segment_color = "black",

    show_entropy = TRUE,
    point_shape = 21,
    point_stroke = 1,
    point_color = "white", 
    point_factor = c(2, 10),
    point_fill_cols = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(300),

    show_pie = TRUE,
    pie_alpha = 1,
    pie_color = point_color,
    pie_stroke = 0.1,

    legend_position = "right",
    nrow = 2
) {

    celltype <- enquo(celltype)

    centers <- centers %>%
        mutate(age = factor(as.character(age), levels = c("Young", "Old"))) %>%
        arrange(age, sample) %>%
        mutate(sample = factor(sample, levels = unique(sample)))

    segment_neighbor <- segment_neighbor %>%
        mutate(sample = factor(sample, levels = levels(centers$sample)))

    print(unique(centers$sample))

    p <- ggplot()

    #
    p <- p +
        geom_segment(
            data = segment_neighbor,
            aes(x = x1, y = y1, xend = x2, yend = y2),
            alpha = segment_alpha,
            color = segment_color,
            size = spot_size)

    #
    df_spot <- centers %>%
        distinct(orig.ident, sample, col, row, !!celltype, center_col, center_row)

    p <- p +
        geom_point(
            aes(x = col, y = row, fill = !!celltype),
            data = df_spot,
            size = spot_size,
            shape = 21,
            color = spot_color,
            stroke = spot_stroke
            ) +
        scale_fill_manual(values = spot_fill)
    
    if (spot_guides_hide) {
        p <- p + guides(fill = "none")
    }

    #
    if(show_pie) {
        df_centers_pie <- centers %>%
            filter(!is.na(label)) %>%
            group_by(orig.ident, sample, label, spot_num, center_col, center_row, !!celltype) %>%
            summarise(celltype_num = n()) %>%
            group_by(orig.ident, sample, label, spot_num, center_col, center_row) %>%
            arrange(!!celltype)

        p <- p +
            geom_arc_bar(
                data = df_centers_pie,
                aes(x0 = center_col, y0 = center_row, 
                    r0 = 0, 
                    r = (spot_num - min(spot_num)) / (max(spot_num) - min(spot_num)) * (point_factor[2] - point_factor[1]) / 2 + point_factor[1] * 1.5,
                    amount = celltype_num, 
                    fill = !!celltype),
                color = pie_color,
                alpha = pie_alpha,
                size = pie_stroke,
                stat = "pie"
            )
    }

    #
    if (show_entropy) {
        df_centers <- centers %>%
            filter(!is.na(label)) %>%
            select(orig.ident, sample, label, spot_num, center_col, center_row) %>%
            distinct() %>%
            inner_join(
                entropy %>%
                select(orig.ident, label, entropy) %>%
                rename(entropy = entropy)
            ) %>%
            mutate(
                r = (spot_num - min(spot_num)) / (max(spot_num) - min(spot_num)) * (point_factor[2] - point_factor[1]) / 2 + point_factor[1] / 2
            )

        p <- p +
            ggnewscale::new_scale_fill() +
            geom_circle(
                data = df_centers,
                aes(x0 = center_col, y0 = center_row, r = r, fill = entropy),
                size = point_stroke,
                color = point_color
                ) +
            scale_fill_gradientn(colours = point_fill_cols)

    }

    # theme
    p <- p +
        Seurat::NoAxes() +
        coord_fixed() +
        theme(
            panel.background = element_blank(),
            strip.background = element_blank(),
            legend.position = legend_position
        ) +
        facet_wrap(vars(sample), nrow = nrow) +
        scale_y_reverse()

    p
}
