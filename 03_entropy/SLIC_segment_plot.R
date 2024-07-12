
# ################################################################
#
# plot
#
# ################################################################

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
                        if(x+1 <= dim(mat)[1] && mat[x, y] == mat[x + 1,y]) 
                            res <- res %>%
                                bind_rows(tibble(x1 = x, y1 = y, x2 = x + 1, y2 = y))
                        if(y+1 <= dim(mat)[2] && mat[x, y] == mat[x, y + 1]) 
                            res <- res %>%
                                bind_rows(tibble(x1 = x, y1 = y, x2 = x, y2 = y + 1))
                        res
                    }) %>%
                    distinct()

            }, future.chunk.size = Inf)
        ) %>%
        unnest(data)
}

SLIC_segment_plot <- function(
    centers, 
    entropy,
    
    n_work = 3,
    segment_neighbor = segment_neighbor_check(centers),
    celltype = "celltype",

    spot_size = 1,
    spot_color = "black",
    spot_fill = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(centers[[celltype]]))),
    spot_guides_hide = TRUE,
    spot_stroke = 0.1,
    segment_alpha = 0.5,
    segment_color = "black",

    show_entropy = TRUE,
    point_shape = 21,
    point_stroke = 1,
    point_color = "white", 
    point_factor = c(2, 10),
    showing_mt_mean = TRUE,
    point_fill_cols = ifelse(showing_mt_mean, 
        colorRampPalette(brewer.pal(9, "YlOrRd")),
        colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral"))))(300),
    point_low_fill_color = "gray",

    show_pie = TRUE,
    pie_alpha = 1,
    pie_color = point_color,
    pie_stroke = 0.1,

    samples = NA,
    legend_position = "right",
    nrow = 2) {

    celltype <- sym(celltype)

    # hpts <- chull(X)
    # hpts <- alphahull::ashape(X, alpha = 0.5)

    if(!is.na(samples)) {
        centers <- centers %>%
            filter(sample %in% samples)
        segment_neighbor <- segment_neighbor %>%
            filter(sample %in% samples)
    }

    centers <- centers %>%
        mutate(age = factor(as.character(age), levels = c("Young", "Old"))) %>%
        arrange(age, sample) %>%
        mutate(sample = factor(sample, levels = unique(sample)))

    segment_neighbor <- segment_neighbor %>%
        mutate(sample = factor(sample, levels = levels(centers$sample)))

    print(unique(centers$sample))

    p <- ggplot()

    # 连线
    p <- p +
        geom_segment(
            data = segment_neighbor,
            aes(x = x1, y = y1, xend = x2, yend = y2),
            alpha = segment_alpha,
            color = segment_color,
            size = spot_size)

    # 背景点
    df_spot <- centers %>%
        select(orig.ident, sample, col, row, !!celltype, center_col, center_row) %>%
        distinct()

    p <- p +
        geom_point(
            data = df_spot,
            aes_string(x = "col", y = "row", fill = as_label(celltype)),
            size = spot_size,
            shape = 21,
            color = spot_color,
            stroke = spot_stroke
            ) +
        scale_fill_manual(values = spot_fill)

    if (spot_guides_hide) {
        p <- p + guides(fill = "none")
    }

    # 饼图
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

    # 熵值
    if(show_entropy) {
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

        if (showing_mt_mean) {
            te_mean <- quantile(df_centers$entropy, 0.25) 

            df_centers_h <- df_centers %>%
                filter(entropy > te_mean)

            df_centers_l <- df_centers %>%
                filter(entropy <= te_mean)

            p <- p +
                ggnewscale::new_scale_fill() +
                geom_circle(
                    data = df_centers_h,
                    aes(x0 = center_col, y0 = center_row, r = r, fill = entropy),
                    size = point_stroke,
                    color = point_color
                    ) +
                scale_fill_gradientn(colours = point_fill_cols) +
                geom_circle(
                    data = df_centers_l,
                    aes(x0 = center_col, y0 = center_row, r = r),
                    size = point_stroke,
                    color = point_color,
                    fill = point_low_fill_color
                    )

        } else {
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
