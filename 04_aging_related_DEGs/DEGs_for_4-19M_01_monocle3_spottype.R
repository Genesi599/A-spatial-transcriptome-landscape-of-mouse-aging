
library(Seurat)
library(monocle3)
library(tidyverse)

source("loess_smooth_heatmap.R")

wd <- 
rds_path <- 

########################################################

deg_filter <- function(gene_df) {
    gene_df_filtered <- gene_df %>%
        mutate(gene = gene_short_name) %>%
        filter(!grepl('^Gm|^mt-|^Rps|^Rpl|^Hb[ab]', gene))
    return(gene_df_filtered)
}

get_age_DEG_mat <- function(meta, cell, exp_mat, out_dir, prefix, n_cores=1){
    pd <- subset(meta, celltype %in% cell)

    # downsample by age group
    cell_min <- min(table(pd$age))
    cell_min <- ifelse(cell_min < 500, cell_min, 500)
    set.seed(2024)
    pd <- pd %>% 
        rownames_to_column('cellid') %>% 
        group_by(age) %>% 
        sample_n(cell_min) %>% 
        as.data.frame() %>% 
        arrange(age) %>% 
        column_to_rownames('cellid')

    exp_mat_tmp <- exp_mat[, rownames(pd)]
    exp_mat_tmp <- exp_mat_tmp[rowSums(exp_mat_tmp) != 0, ]
    fd <- data.frame(gene_short_name = rownames(exp_mat_tmp), row.names = rownames(exp_mat_tmp))

    cds <- new_cell_data_set(exp_mat_tmp, cell_metadata = pd, gene_metadata = fd)
    cds <- cds[, colSums(exprs(cds)) != 0]
    cds <- estimate_size_factors(cds)

    gene_fits <- fit_models(cds, model_formula_str = "~age", verbose = T, cores=n_cores)
    fit_coefs <- coefficient_table(gene_fits)
    gene_res <- subset(fit_coefs, term == 'age' & q_value < 0.05) %>%
            select(gene_short_name, term, q_value, estimate)

    gene_res <- deg_filter(gene_res)
    gene_res$celltype <- cell
    gene_res %>%
        write_csv(paste0(out_dir, '/', prefix, '_', str_replace(cell, "[/&]", "_"), '.csv'))

    mat <- int_obj %>%
        get_log2cpm_mat(
        age, celltype,
        gene_list = gene_res$gene,
        sample_number = 50,
        assay = "Spatial",
        replace = T)
    mat <- smooth_data(mat)

    gene_list <- mat$gene_anno %>%
        column_to_rownames("id") %>%
        select(-slope)

    mat_draw <- mat$mat.smooth[rownames(gene_list), ]
    saveRDS(mat_draw, paste0(out_dir, '/', prefix, '_mat_draw_', str_replace(cell, "[/&]", "_"), '.rds'))
}


########################################################

setwd(wd)

#
int_obj <- readRDS(rds_path)
meta <- int_obj@meta.data %>%
    mutate(
        age = str_extract(age, "\\d+"),
        age = as.numeric(age)
    )
exp_mat <- GetAssayData(int_obj, slot = 'data')

#
for (cell in as.character(unique(meta$celltype)))
try({
    get_age_DEG_mat(meta, cell, exp_mat, out_dir = wd, prefix = 'age_DEG_celltype')
    print(paste('==================', cell, '=================='))
})

