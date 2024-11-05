

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

#' Adjacency score analysis
#' 
#' @description  https://doi.org/10.1016/j.cmet.2021.07.018
#' 
#' @param .data seurat obj
#' @param ...  columns in meta.data, and will be reserve in result
#' @param spottype use for classifying spots
#' @param d farest distance between two spot which will be calculated
#' @param R resample times
#' @param spot_number_cutoff min spot number in one image
#' @param n_work number of threads
#' 
#' @return tibble obj <orig.ident, spottype, ..., Neighbor,N_spot,Neighbor_spot,K_z,K_o,K_e,K_sd>
#' spottype: spot type A
#' Neighbor: spot type A's neighbor spot type B
#' N_spot: number of A
#' Neighbor_spot: number of B
#' K_z: Adjacency score
#' K_o: observed score
#' K_e: mean of simulations
#' K_sd: sd of simulations
#' 
neighbor_score <- function(
    .data, 
    ..., 
    spottype = celltype, 
    d = 1, R = 50, 
    spot_number_cutoff = 10, 
    n_work = 3
) {
    
    spottype <- enquo(spottype)
    group_vars <- enquos(..., .named = TRUE)

    set.seed(2022)

    plan(multisession, workers = n_work)
    options(future.globals.maxSize = Inf)
    options(future.seed=TRUE)
    message(">> how many cores can use now: ", nbrOfWorkers())

    .data@meta.data %>%
        rownames_to_column(var = "cellid") %>%
        select(cellid, orig.ident, !!spottype, !!!group_vars) %>%
        # get Coordinates
        left_join(
            GetAllCoordinates(.data)
        ) %>%
        # NeighborScore
        group_by(orig.ident) %>%
        group_split() %>%
        future_lapply(function(df) {
            print(df$orig.ident[1])

            mat <- df %>%
                column_to_rownames("cellid") %>%
                select(row, col) %>%
                as.matrix() %>%
                dist(method = "euclidean") %>%
                as.matrix()

            .sample_K <- df %>%
                select(!!spottype) %>%
                distinct() %>%
                inner_join(., ., by = character())
            colnames(.sample_K) <- c("A","B")

            .sample_K <- .sample_K %>%
                split(., row(.)) %>%
                map_dfr(~ {
                    ct <- .x

                    A_spot <- df %>% filter(!!spottype == ct$A) %>% select(cellid)
                    B_spot <- df %>% filter(!!spottype == ct$B) %>% select(cellid)

                    .res <- tibble(
                        !!spottype := ct$A,
                        Neighbor = ct$B,
                        N_spot = nrow(A_spot),
                        Neighbor_spot = nrow(B_spot),
                        Image_spot = nrow(df),
                        d = d)

                    if (nrow(A_spot) > spot_number_cutoff & nrow(B_spot) > spot_number_cutoff) {
                        print(str_c(ct$A, " => ", ct$B))

                        .sample <- seq(R) %>%
                            map_dfr(~{
                                # print(str_c("random sampling: ", .))
                                .sample_A <- df %>%
                                    select(cellid) %>%
                                    slice_sample(n = nrow(A_spot))

                                .sample_B <- df %>%
                                    select(cellid) %>%
                                    slice_sample(n = nrow(B_spot))

                                .sample <- mat[.sample_A$cellid, .sample_B$cellid] %>%
                                    apply(1, function(l) sum(0 < l & l < d + 1, na.rm = TRUE))
                                c("mean" = mean(.sample, na.rm = TRUE), "sd" = sd(.sample, na.rm = TRUE))
                            })

                        .res <- .res %>%
                            bind_cols(
                                tibble(
                                    K_e = mean(.sample$mean, na.rm = TRUE),
                                    K_sd = mean(.sample$sd, na.rm = TRUE),
                                    K_o = mat[A_spot$cellid, B_spot$cellid] %>%
                                        apply(1, function(l) sum(0 < l & l < d + 1, na.rm = TRUE)) %>%
                                        mean(na.rm = TRUE),
                                    K_z = (K_o - K_e) / K_sd
                                )
                            )
                    }

                    .res
                })

            df %>%
                select(orig.ident, !!spottype, !!!group_vars) %>%
                distinct() %>%
                left_join(.sample_K)
        }) %>%
        bind_rows()
}
