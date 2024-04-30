
# ################################################################
#
# pulic
#
# ################################################################

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

# ################################################################
#
#  Segment
#
# ################################################################

getdist <- function(spots, centers, m, window) {

    dist_s <- proxy::dist(
        spots %>% select(col, row),
        centers %>% select(col, row),
        method = "Euclidean") %>%
        unclass() %>%
        as_tibble(rownames = "cellid") %>%
        pivot_longer(-cellid, names_to = "label", values_to = "d_s") %>%
        filter(d_s <= window)

    # 
    dist_exp <- proxy::dist(
        spots %>% select(starts_with("PC_")),
        centers %>% select(starts_with("PC_")),
        method = "Euclidean") %>%
        unclass() %>%
        as_tibble(rownames = "cellid") %>%
        pivot_longer(-cellid, names_to = "label", values_to = "d_e")

    #
    dist <- dist_s %>%
        inner_join(dist_exp) %>%
        mutate(
            dist = sqrt(m * d_s^2 + d_e^2)
        ) %>%
        group_by(cellid) %>%
        slice_min(dist) %>%
        select(cellid, label, dist)

}


single_SLIC <- function(df, ..., window = 20, min_spot = max(window^2 * 0.1, 5), m = 0.8) {
    print(df[1, "orig.ident"])

    # 1. innitessd
    center_df <- df %>%
        mutate(
            col = ceiling(col /  window),
            row = ceiling(row /  window)
        )  %>%
        group_by(col, row) %>%
        summarise(
            spot_num = n(),
            across(starts_with("PC_"), mean)) %>%
        ungroup() %>%
        mutate(
            col = floor((col - 0.5) * window),
            row = floor((row - 0.5) * window)
        ) %>%
        filter(spot_num > min_spot) %>%
        mutate(
            label = str_c("V", seq(n()))
        ) %>%
        arrange(label)

    df <- df %>%
        mutate(
            label = NaN,
            dist = Inf
        )

    # print(nrow(df) / (window * window))

    res <- Inf
    step <- 1
    while(res > 1) {
        # 2. update pixel
        dists <- getdist(
            df %>% column_to_rownames("cellid"), 
            center_df  %>% column_to_rownames("label"), 
            m, window)

        # 3. update dist
        df <- df %>%
            select(-c(label, dist)) %>%
            left_join(dists)

        # 4. update centers
        center_df_new <- df %>%
            group_by(label) %>%
            summarise(
                spot_num = n(),
                across(c(col, row, starts_with("PC_")), mean)) %>%
            ungroup() %>%
            arrange(label) %>%
            filter(!is.na(label))

        center_df <- center_df %>%
            filter(label %in% center_df_new$label) %>%
            arrange(label)

        res <- as.matrix(center_df %>% select(-label)) - as.matrix(center_df_new %>% select(-label))

        res <- colSums(res)
        # print(res)
        # res <- sum(abs(res[c("col", "row")]))
        res <- sum(abs(res))
        print(str_c("Step: ", step, " Segment num: ", nrow(center_df_new), " NA spot: ", sum(is.na(df$label)), " Diff: ", res))

        center_df <- center_df_new
        step <- step + 1
    }

    df %>%
        select(cellid, orig.ident, col, row, label) %>%
        left_join(
            center_df %>%
                rename(
                    center_col = col,
                    center_row = row
                )
        ) %>%
        select(cellid, orig.ident, 
            col, row, 
            label, center_col, center_row, spot_num,
            ...)
}

STSLIC <- function(.data, nPC = 1:15, window = 20, m = 3, min_spot = max(window^2 * 0.1, 5)) {
    message("run SLIC")

    .data@meta.data %>%
        rownames_to_column(var = "cellid") %>%
        # Get Coordinates
        left_join(
            GetAllCoordinates(.data)
        ) %>%
        # Get PC
        left_join(
            FetchData(
                object = .data,
                vars = str_c("PC_", nPC)
            ) %>%
            # scale() %>%
            as_tibble(rownames = "cellid")
        ) %>%
        # SLIC
        select(cellid, orig.ident, col, row, starts_with("PC_")) %>%
        group_by(orig.ident) %>%
        group_split() %>%
        future_lapply(single_SLIC,
            window = window, m = m, min_spot = min_spot,
            future.chunk.size = Inf) %>%
        bind_rows()
}

# ################################################################
#
# Neighbors Count
#
# ################################################################

ImageSpotNeighborsCount <- function(meta.data, celltype_col, neighbor_range, ...) {
    .sample_coord <- meta.data %>%
        select(row,col,!!celltype_col, ...) %>%
        mutate(.celltype = as.numeric(factor(!!celltype_col)))
    # print(head(.sample_coord))
        
    # build celltype matrix
    .celltype_mat <- Matrix::sparseMatrix(
        .sample_coord$row,
        .sample_coord$col,
        x = .sample_coord$.celltype)

    row_r <- dim(.celltype_mat)[1]
    col_r <- dim(.celltype_mat)[2]
    print(str_c("row range: ", row_r, " col range: ", col_r))

    # count
    .sample_coord %>%
        mutate(
            neighbors = map2(row, col, ~{
                .celltype_code <- .celltype_mat[
                        max(0, .x - neighbor_range):min(.x + neighbor_range, dim(.celltype_mat)[1]),
                        max(0, .y - neighbor_range):min(.y + neighbor_range, dim(.celltype_mat)[2])
                    ] %>%
                    as.vector()

                # remove missing
                .celltype_code <- .celltype_code[which(.celltype_code > 0)]
                # remove self
                .celltype_code <- .celltype_code[-match(.celltype_mat[.x,.y], .celltype_code)]
                # encode
                # index:种类, value:数量比例
                table(.celltype_code)
            })
        )
}

# ################################################################
#
# Neighbors Count
#
# ################################################################

SegmentEntropy <- function(.data, 
    ..., 
    celltype_col = celltype, 
    neighbor_range = 1, 
    roi_col,
    n_work = 3) {

    celltype_col = enquo(celltype_col)
    roi_col = enquo(roi_col)
    group_vars <- enquos(..., .named = TRUE)

    library(future)
    library(future.apply)

    plan(multisession, workers = n_work)
    options(future.globals.maxSize = Inf)
    message("outside >> how many cores can use now: ", nbrOfWorkers())

    set.seed(2023)

    df <- .data@meta.data %>%
        as_tibble(rownames = "cellid") %>%
    # spot2niche
        group_by(orig.ident) %>%
        group_nest() %>%
        mutate(
            data = future_lapply(data, function(df) {
                ImageSpotNeighborsCount(
                    df,
                    celltype_col,
                    neighbor_range,
                    cellid, !!roi_col
                )
            }, future.chunk.size = Inf)
        ) %>%
        unnest(data) %>%
        filter(!is.na(!!roi_col)) %>%
    # get entropy
        group_by(orig.ident, !!roi_col, !!celltype_col, .celltype, neighbors) %>%
        summarise(f_ij = n()) %>%
        group_by(orig.ident, !!roi_col) %>%
        mutate(
            p_ij = f_ij / sum(f_ij),
            entropy_ij = - p_ij * log2(p_ij),
            condition = n()
        ) %>%
        group_by(orig.ident, !!roi_col, condition) %>%
        summarise(entropy = sum(entropy_ij)) %>%
        mutate(Pielou = entropy / log(condition)) %>%
        left_join(
            .data@meta.data %>%
                select(orig.ident, !!!group_vars) %>%
                distinct()
        )
}


coverage_check <- function(df, select_spot_type = c(1,2,3)) {

    df_roi <- df[chull(df$col, df$row), ]

    x_range <- range(df$col)
    y_range <- range(df$row)

    chull_range <- inner_join(
            tibble(col = seq(x_range[1], x_range[2])),
            tibble(row = seq(y_range[1], y_range[2])),
            by = character()
        ) %>%
        mutate(
            roi_type = sp::point.in.polygon(col, row, df_roi$col, df_roi$row)
        ) %>%
        filter(roi_type %in% select_spot_type) %>%
        nrow()

    tibble(spot_num = nrow(df), chull_range = chull_range, coverage = nrow(df)/chull_range)

}

#' tissue SLIC entropy
#'
#' @description 使用SLIC分割样本，根据spot分布，度量切片组织结构的混乱程度
#' 
#' @param .data seurat 对象
#' @param celltype_col meta.data中细胞类别列列名 默认celltype
#' @param ...  meta.data中细胞类别列列名
#' @param neighbor_range spot邻域范围 默认 1
#' 
#' @param nPC 使用的PC 默认 PC1~15
#' @param window 超像素尺寸 
#' @param m 空间距离调节因子
#' @param min_spot 初始超像素最少spot数
#' 
#' @param n_work 线程数
#'
TissueSLICEntropy <- function(.data,
    ..., 
    celltype_col = celltype,
    neighbor_range = 1,

    nPC = 1:15,
    window = 20,
    m = 3,
    min_spot = window^2 * 0.1,

    n_work = 3) {

    library(future)
    library(future.apply)

    plan(multisession, workers = n_work)
    options(future.globals.maxSize = Inf)
    message(">> how many cores can use now: ", nbrOfWorkers())

    # SLIC and Get Coordinates
    .data@meta.data <- .data@meta.data %>%
        rownames_to_column(var = "cellid") %>%
        left_join(
            STSLIC(.data, nPC = nPC, window = window, m = m, min_spot = min_spot)
        ) %>%
        column_to_rownames("cellid")

    res <- SegmentEntropy(.data,
        celltype_col = !!enquo(celltype_col),
        neighbor_range = neighbor_range,
        roi_col = label,
        n_work = n_work)
    
    res <- res %>%
        inner_join(
            .data@meta.data %>%
                select(orig.ident, col, row, label) %>%
                distinct() %>%
                group_nest(orig.ident, label) %>%
                mutate(
                    coverage = map(data, coverage_check)
                ) %>%
                unnest(coverage) %>%
                select(-data)
        ) %>%
        mutate(
            entropy.original = entropy,
            entropy = entropy / log(chull_range)
        ) %>%
        ungroup()

    list(
        "centers" = .data@meta.data %>%
            rownames_to_column(var = "cellid") %>%
            select(orig.ident, ..., !!enquo(celltype_col), cellid, col, row, label, center_col, center_row, spot_num) %>%
            distinct(),
        "entropy" = res,
        "nPC" = nPC,
        "window" = window,
        "m" = m,
        "min_spot" = min_spot
    )
}

