
library(Seurat)
library(tidyverse)


wd <- 
tissue <- 
rds_path <- 
n_work <- 20
out <- 

########################################################

O_Y_DEG <- function(data, n_work = 14) {

    plan(multisession, workers=n_work)
    options(future.globals.maxSize= Inf)
    message("outside >> how many cores can use now: ", nbrOfWorkers())

    OvsY_markers_full <- levels(data) %>%
        future_lapply(function(ct)  {
            print(ct)

            t_df <- tryCatch({
                data %>%
                    subset(celltype == ct) %>%
                    FindMarkers(ident.1 = "25M", group.by = "age", test.use = "wilcox") %>%
                    as_tibble(rownames = "gene") %>%
                    mutate(celltype = ct) %>%
                    return()
            }, error = function(err) {
                return(tibble())
            })

        }, future.chunk.size = Inf) %>%
        bind_rows()

    return(OvsY_markers_full)
}

deg_filter <- function(df, logfc_cutoff = 0.25) {
    df_filtered <- df %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > logfc_cutoff) %>%
        filter(pct.1 > 0.1, pct.2 > 0.1) %>%
        filter(!grepl('^Gm|^mt-|^Rps|^Rpl|^Hb[ab]', gene))
    return(df_filtered)
}

########################################################

setwd(wd)

rds <- readRDS(rds_path)

DefaultAssay(rds) <- "Spatial"
deg_df <- O_Y_DEG(rds, n_work) %>%
    deg_filter() %>%
    mutate(tissue = tissue) %>%
    mutate(type = if_else(avg_log2FC > 0, "Up", "Down"))

deg_df %>%
    write_csv(str_c(out, "_stDEGs.csv"))
