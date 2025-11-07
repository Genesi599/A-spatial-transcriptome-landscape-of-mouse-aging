# niche_grade_entropy.R

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
            tibble::rownames_to_column(var = "cellid")
        })
}

# ################################################################
#
# Neighbors Count
#
# ################################################################

image_spot_neighbors_count <- function(meta.data, celltype_col, neighbor_range, ...) {
    .sample_coord <- meta.data %>%
        dplyr::select(row,col,!!celltype_col, ...) %>%
        dplyr::mutate(.celltype = as.numeric(factor(!!celltype_col)))
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
        dplyr::mutate(
            neighbors = map2(row, col, ~{
                .celltype_code <- .celltype_mat[
                        max(0, .x - neighbor_range):min(.x + neighbor_range, dim(.celltype_mat)[1]),
                        max(0, .y - neighbor_range):min(.y + neighbor_range, dim(.celltype_mat)[2])
                    ] %>%
                    as.vector()

                .celltype_code <- .celltype_code[which(.celltype_code > 0)]
                .celltype_code <- .celltype_code[-match(.celltype_mat[.x,.y], .celltype_code)]

                table(.celltype_code)
            })
        )
}

# ################################################################
#
# Neighbors Count
#
# ################################################################


#' SSS niche gradient entropy definition
#'
#' @param .data tibble obj
#' @param celltype_col column in meta.data, use for classifying spots
#' @param ...  columns in meta.data, and will be reserve in result
#' @param roi_col column in meta.data, use for classifying niche gradient
#' @param neighbor_range neighborhood range
#' @param R resample times
#' @param n_work number of threads
#'
niche_grade_entropy <- function(
    .data, 
    ..., 
    celltype_col = seurat_clusters, 
    neighbor_range = 1, 
    roi_col, 
    R = 100, 
    n_work = 3
) {

    celltype_col = enquo(celltype_col)
    roi_col = enquo(roi_col)
    group_vars <- enquos(..., .named = TRUE)

    library(future)
    library(future.apply)

    plan(multisession, workers=n_work)
    options(future.globals.maxSize= Inf)
    options(future.seed=TRUE)
    message("outside >> how many cores can use now: ", nbrOfWorkers())

    set.seed(2023)

    df <- .data@meta.data %>%
        as_tibble(rownames = "cellid") 
    # Get Coordinates
    if(!("col" %in% colnames(df) && "row" %in% colnames(df))) {
        df <- df %>%
            dplyr::left_join(
                GetAllCoordinates(.data)
            )
    }
    # spot2niche
    df <- df %>%
        dplyr::group_by(age, orig.ident) %>%
        group_nest() %>%
        dplyr::mutate(
            data = future_lapply(data, function(df) {
                image_spot_neighbors_count(
                    df,
                    celltype_col,
                    neighbor_range,
                    cellid, !!roi_col
                )
            }, future.chunk.size = Inf)
        ) %>%
        unnest(data) %>%
        dplyr::filter(!is.na(!!roi_col)) %>%
        dplyr::group_by(age, !!roi_col) %>%
        dplyr::mutate(min_roi_spot_num = n())
    # get roi_min_spot_num
    df$min_roi_spot_num <- min(df$min_roi_spot_num)

    # get entropy
    df <- df %>%
        dplyr::group_by(age, !!roi_col, min_roi_spot_num) %>%
        group_nest() %>%
        dplyr::mutate(
            data = map2(data, min_roi_spot_num, ~{
                message("")
                message(str_c("## sample spot num: ", .y))

                data <- .x
                min_roi_spot_num <- .y
                # resample
                res <- seq(R) %>%
                    as.list() %>%
                    future_lapply(function(idx) {
                        data %>%
                            slice_sample(n = min_roi_spot_num) %>%
                            dplyr::group_by(!!celltype_col, .celltype, neighbors) %>%
                            dplyr::summarise(f_ij = n()) %>%
                            ungroup() %>%
                            dplyr::mutate(
                                p_ij = f_ij / sum(f_ij),
                                entropy_ij = - p_ij * log2(p_ij),
                                condition = n()
                            ) %>%
                            dplyr::group_by(condition) %>%
                            dplyr::summarise(entropy = sum(entropy_ij)) %>%
                            dplyr::mutate(
                                rep_idx = idx,
                                Pielou = entropy / log(condition))
                    }, future.chunk.size = Inf, future.seed = TRUE) %>%
                    bind_rows() 
                    
            })
        ) %>%
        unnest(data)
}

