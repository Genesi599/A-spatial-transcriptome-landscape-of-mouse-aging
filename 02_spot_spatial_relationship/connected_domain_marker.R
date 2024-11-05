
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

#' Spatially connected domain identification analysis
#' 
#' @param .data seurat obj
#' @param spottype filter condition, eg: celltype == "SS macro"
#' @param n_work number of threads
#' 
#' @return seurat obj, add result in metadata. <roi_id, roi_size>
#' roi_id: connected range id,
#' roi_size: size of connected range
#' 
connected_domain_marker <- function(
    .data, 
    spottype,
    n_work = 12
) {

    library(future)
    library(future.apply)

    plan(multisession, workers = n_work)
    options(future.globals.maxSize = Inf)
    message(">> how many cores can use now: ", nbrOfWorkers())

    spottype = enquo(spottype)

    df <- .data@meta.data %>%
        as_tibble(rownames = "cellid") %>%
        # get Coordinates
        left_join(
            GetAllCoordinates(.data)
        )  %>%
        filter(!!spottype) %>%
        group_split(orig.ident) %>%
        lapply(function(sample) { 
            print(sample$sample[[1]])

            sample_coord <- sample %>%
                select(cellid, row, col) %>%
                mutate(
                    id = seq(n()),
                    roi_id = 0,
                    roi_size = 1
                )

            cellid_mat <- Matrix::sparseMatrix(
                    sample_coord$row,
                    sample_coord$col,
                    x = sample_coord$id)

            # BFS
            roi_number <- 1
            for(i in 1:nrow(sample_coord)) {

                if(sample_coord[i,]$roi_size > 1) next;

                df_t <- sample_coord[i,]

                roi_size <- 0

                while(roi_size != nrow(df_t)) {
                    roi_size <- nrow(df_t)
                    df_t <- df_t %>%
                        pmap_dfr(~{
                            x <- ..2
                            y <- ..3
                            cellid_list <- cellid_mat[
                                    max(0, x - 1):min(x + 1, dim(cellid_mat)[1]),
                                    max(0, y - 1):min(y + 1, dim(cellid_mat)[2])
                                ] %>%
                                as.vector()
                            cellid_list <- cellid_list[which(cellid_list > 0)]

                            sample_coord[cellid_list,]
                        }) %>%
                        distinct()
                }

                sample_coord[df_t$id,]$roi_id <- roi_number
                sample_coord[df_t$id,]$roi_size <- nrow(df_t)
                roi_number <- roi_number + 1
            }


            sample_coord
        }) %>%
        bind_rows()

    .data@meta.data <- .data@meta.data %>%
        as_tibble(rownames = "cellid") %>%
        left_join(
            df %>%
                select(cellid, roi_id, roi_size)
        ) %>%
        column_to_rownames("cellid")

    .data
}

