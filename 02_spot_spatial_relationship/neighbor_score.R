

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

#' Neighbor Score
#' 
#' @title 邻域指数
#' @description 参考文献 https://doi.org/10.1016/j.cmet.2021.07.018
#' 
#' @import Seurat
#' @import tidyverse
#' @import future
#' @import future.apply
#' 
#' @param .data seurat对象
#' @param ... meta.data中的列, 会merge到结果中
#' @param spottype 分类列名 默认 celltype
#' @param d 统计距离d之内的spot 默认 2
#' @param R 抽样次数 默认 50
#' @param spot_number_cutoff 某类细胞单芯片最少细胞数 默认 10
#' @param n_work 线程数 默认 3
#' 
#' @return tibble <orig.ident, spottype, ..., Neighbor,N_spot,Neighbor_spot,K_z,K_o,K_e,K_sd>
#' 
NeighborScore <- function(.data, ..., spottype = celltype, d = 1, R = 50, spot_number_cutoff = 10, n_work = 3) {
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
