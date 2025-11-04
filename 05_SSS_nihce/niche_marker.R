GetAllCoordinates <- function(.data) {
  .data@images %>%
    names() %>%
    unique() %>%
    map_dfr(~ {
      GetTissueCoordinates(
        .data,
        image = .x,
        cols = c("row", "col"),
        scale = NULL
      ) %>%
        rownames_to_column(var = "cellid")
    })
}

single_marker <- function(df, intra_df, spot_type, dist_method, FUN, zero_check = FALSE) {
  if (length(intra_df) != 0) {
    all_df <- df %>%
      column_to_rownames("cellid") %>%
      select(row, col)

    mat <- proxy::dist(all_df, intra_df, method = dist_method) %>%
      as.matrix()

    spot_dist <- tibble(cellid = rownames(mat))
    spot_dist[[spot_type]] <- rowMeans(mat, na.rm = TRUE)

    if (!is.na(FUN)) {
      spot_dist[[spot_type]] <- FUN(spot_dist[[spot_type]])
    }

    res <- df %>%
      left_join(spot_dist, by = "cellid")

  } else {
    res <- df %>%
      mutate(!!spot_type := Inf)
  }

  res %>% select(-c(row, col))
}

# ------------------ 主函数 ------------------ #

niche_marker <- function(
  .data,
  marker,
  spot_type,
  slide = "orig.ident",
  dist_method = "Euclidean",
  FUN = NA,
  n_work = 3
) {
  # 👇 关键修改：列名字符串，不再是 quosure
  marker <- as.character(substitute(marker))
  spot_type <- as.character(substitute(spot_type))
  slide <- as.character(substitute(slide))

  library(future)
  library(future.apply)

  plan(multisession, workers = n_work)
  options(future.globals.maxSize = Inf)
  message(">> how many cores can use now: ", nbrOfWorkers())

  .data@meta.data <-
    .data@meta.data %>%
    rownames_to_column(var = "cellid") %>%
    left_join(GetAllCoordinates(.data), by = "cellid") %>%
    group_by(.data[[slide]]) %>%
    group_split() %>%
    future_lapply(function(df) {
      slide_name <- df[[slide]][1]
      print(slide_name)

      intra_df <- df %>%
        filter(.data[[marker]]) %>%
        column_to_rownames("cellid") %>%
        select(row, col)

      single_marker(df, intra_df, spot_type = spot_type,
                    dist_method = dist_method, FUN = FUN)
    }, future.chunk.size = Inf) %>%
    bind_rows() %>%
    column_to_rownames(var = "cellid")

  .data@meta.data <- .data@meta.data[colnames(.data), ]

  return(.data)
}