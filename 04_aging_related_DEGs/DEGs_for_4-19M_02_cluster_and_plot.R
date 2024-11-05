#!/data/home/quj_lab/zhangbin/anaconda3/envs/R_sc/bin/R

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

wd <- 
tissue <- 

##########################################

setwd(wd)

all_mat <- list.files(wd, full.names = TRUE) %>%
    .[grepl('age_DEG_celltype_mat_draw_', .)] %>%
    map_dfr(function(x){
        tmp_mat <- readRDS(x)
        colnames(tmp_mat) <- 1:150
        cell <- str_extract(x, '(?<=draw_).*?(?=\\.rds)')
        rownames(tmp_mat) <- paste(rownames(tmp_mat), cell, sep='__')
        as.data.frame(tmp_mat)
    }) %>%
    data.matrix()

# heatmap ####
hm_draw <- function(mat_draw, prefix, num_clusters, plot_height, plot_width){
    grad_cols <- c('4' = 'grey85', '13' = '#F9C670', '19' = '#f39800')
    cell_list <- data.frame(age = rep(c('4', '13', '19'), each = 50))
    p <- pheatmap::pheatmap(
        mat_draw,
        # cluster_rows = F,
        cluster_cols = F,
        show_rownames = F,
        show_colnames = F,
        annotation_names_row = F,
        #annotation_row = gene_list,
        annotation_col = cell_list,
        # annotation_legend = F,
        # legend = F,
        scale = "row",
        color = colorRampPalette(colors = brewer.pal(n = 11, name = "RdYlBu") %>% rev())(100),
        # gaps_row = gaps_row,
        cutree_rows = num_clusters,
        annotation_colors = list(
            age = grad_cols
        )
    )

    gene_list <- data.frame(cluster = factor(cutree(p$tree_row, num_clusters)))

    p <- pheatmap::pheatmap(
        mat_draw,
        # cluster_rows = F,
        cluster_cols = F,
        show_rownames = F,
        show_colnames = F,
        annotation_names_row = F,
        annotation_row = gene_list,
        annotation_col = cell_list,
        # annotation_legend = F,
        # legend = F,
        scale = "row",
        color = colorRampPalette(colors = brewer.pal(n = 11, name = "RdYlBu") %>% rev())(100),
        # gaps_row = gaps_row,
        cutree_rows = num_clusters,
        annotation_colors = list(
            age = grad_cols
        )
    )
    ggsave(paste0(prefix, '_cluster', num_clusters, '.pdf'), plot = p, width = plot_width, height = plot_height)
    ggsave(paste0(prefix, '_cluster', num_clusters, '.png'), plot = p, width = plot_width, height = plot_height)

    gene_list$gene <- str_extract(rownames(gene_list), '.*(?=__)')
    gene_list$celltype <- str_extract(rownames(gene_list), '(?<=__).*')
    gene_list$tissue <- basename(tissue)

    gene_list %>%
        write_csv(str_c(prefix, '_gene_list_cluster', num_clusters, '.csv'))
    
    gene_list
}
age_DEG <- hm_draw(all_mat, str_c(tissue, '_age_DEG_celltype_hm'), 8, plot_height = 8, plot_width = 6) 


#
line_plot_df2 <- all_mat %>%
    as.data.frame() %>% 
    as_tibble(rownames = "gene")  %>%
    left_join(
        age_DEG %>%
            select(cluster) %>%
            as_tibble(rownames = "gene")
    ) %>%
    pivot_longer(`1`:`150`, names_to = 'cellid', values_to = 'value') %>%
    arrange(cellid, cluster, value) %>%
    mutate(
        smooth = 'smooth',
        cellid = as.numeric(cellid))

#
p <- line_plot_df2 %>%
    ggplot(aes(x = cellid, y = value)) +
    geom_smooth(aes(group = gene), se = F, color = 'grey85') +
    geom_smooth(aes(group = smooth), colour = "#D7281B") +
    theme(
        panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'white'),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
    facet_wrap(~cluster, ncol = 2, scale = 'free')

ggsave(str_c(tissue, '_age_DEG_celltype_smooth.pdf'), height = 8, width = 4)
ggsave(str_c(tissue, '_age_DEG_celltype_smooth.png'), height = 8, width = 4)
