
rm(list = ls())
gc()

library(Seurat)
library(tidyverse)

library(future)
library(future.apply)


wd <-
out <-
rds_path <- 
gene_number <- 20

setwd(wd)

#####################################################

cell_identity_score <- function(data, feature) {

    data <- AddModuleScore(data,
        features = feature,
        name = "cell_identity_score__")

    df <- data@meta.data %>%
        as_tibble(rownames = "cellid") %>%
        select(cellid, tissue, age, sample, celltype, starts_with("cell_identity_score__")) %>%
        rename_at(vars(starts_with("cell_identity_score__")), ~{
            t <- names(feature)[as.numeric(gsub("cell_identity_score__", "", .x))]
            str_c("cell_identity_score_", t)
        })
}

get_features <- function(df, n) {
    feature <- df %>%
        mutate(
            feature_name = cluster
        ) %>%
        group_by(feature_name) %>%
        top_n(n, -avg_log2FC) %>%
        select(feature_name, gene) %>%
        ungroup() %>%
        group_nest(feature_name) %>%
        pmap(~ {
            t <- ..2[["gene"]]
            setNames(list(t), ..1)
        }) %>%
        unlist(recursive = F)
}

add_sore_to_rds <- function(rds, score_df) {
    rds@meta.data <- rds@meta.data %>%
        as_tibble(rownames = "cellid") %>%
        left_join(
            score_df %>%
                pivot_longer(
                    starts_with("cell_identity_score"),
                    names_prefix = "cell_identity_score_",
                    names_to = c("cell_identity_score_celltype"),
                    values_to = "cell_identity_score"
                ) %>%
                filter(celltype == cell_identity_score_celltype) %>%
                select(cellid, cell_identity_score),
            by = c("cellid")
        ) %>%
        column_to_rownames("cellid")
    rds
}

#########################################################

plan(multisession, workers=20)
options(future.globals.maxSize= Inf)
message("outside >> how many cores can use now: ", nbrOfWorkers())

rds <- readRDS(rds_path)

rds_y <- rds %>%
    subset(age == "2M" | age == "4M")

df <- FindAllMarkers(rds_y) %>%
    as_tibble() %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0.25) 

df %>%
    as_tibble() %>%
    write_csv(str_c(out, "_cell_identity_gene_all.csv"))

score_df <- cell_identity_score(rds, get_features(df, gene_number))

score_df %>% 
    saveRDS(str_c(out, "_cell_identity_score_top", gene_number,".rds"))
