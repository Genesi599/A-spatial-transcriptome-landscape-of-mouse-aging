
rm(list = ls())
gc()

######
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
# library(ggsci)
library(ggpubr)

##################################################################################################################

# bin1 => bin_size
bin_data_progress <- function(infile, pro, bs, print_bin_information = FALSE) {
    dat <- fread(file = infile)
    print(paste0("load file: ", infile))

    if (length(grep("MIDCounts", colnames(dat)) > 0)) {
        dat <- dat %>%
            rename(UMICount = MIDCounts)
    }
    if (length(grep("MIDCount", colnames(dat)) > 0)) {
        dat <- dat %>%
            rename(UMICount = MIDCount)
    }

    # (x, y) to bin(x, y)
    dat[, binx := trunc((x - min(x)) / bs + 1)]
    dat[, biny := trunc((y - min(y)) / bs + 1)]

    # write (x, y) to bin(x, y) information
    if (print_bin_information) {
        dat %>%
            rename(
                !!paste0("bin", bs, ".x") := binx,
                !!paste0("bin", bs, ".y") := biny
            ) %>%
            fwrite(
                paste0(pro, "_bin", bs, "_information.txt"),
                col.names = T,
                row.names = F,
                sep = "\t",
                quote = F
            )
    }

    dat[, x := binx][, binx := NULL]
    dat[, y := biny][, biny := NULL]

    #
    dat <- dat[, sum(UMICount), by = .(geneID, x, y)][, bin_ID := max(x) * (y - 1) + x]
    return(dat);
}

generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
    if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

    spot.radius <- unnormalized.radius / max(dim(image))

    return(new(
        Class = "VisiumV1",
        image = image,
        scale.factors = scalefactors(
            spot = scale.factors$tissue_hires_scalef,
            fiducial = scale.factors$fiducial_diameter_fullres,
            hires = scale.factors$tissue_hires_scalef,
            lowres = scale.factors$tissue_lowres_scalef
        ),
        coordinates = tissue.positions,
        spot.radius = spot.radius
    ))
}

generate_seurat_spatialObj <- function(infile, pro, bs, obj_name, print_bin_information = FALSE) {


    dat <- bin_data_progress(infile, pro, bs, print_bin_information = FALSE)
    print(paste0("bin", bs))


    bin.coor <- dat[, sum(V1), by = .(x, y)]

    ##
    hash.G <- unique(dat$geneID) %>%
        data.frame(
            row.names = .,
            values = seq_along(.)
        )
    gen <- hash.G[dat$geneID, "values"]

    ##
    hash.B <- unique(dat$bin_ID) %>%
        data.frame(
            row.names = sprintf("%d", .),
            values = .
        )
    bin <- hash.B[sprintf("%d", dat$bin_ID), "values"]

    ##
    cnt <- dat$V1

    rm(dat)
    gc()

    ##
    tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))

    tissue_positions_list <- data.frame(
        row.names = paste("BIN", rownames(hash.B), sep = "."),
        tissue = 1,
        row = bin.coor$y,
        col = bin.coor$x,
        imagerow = bin.coor$y,
        imagecol = bin.coor$x
    )

    scalefactors_json <- toJSON(list(
        fiducial_diameter_fullres = 1,
        tissue_hires_scalef = 1,
        tissue_lowres_scalef = 1
    ))

    ##
    print(paste0("build mat start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    mat <- sparseMatrix(i = gen, j = bin, x = cnt)
    rownames(mat) <- rownames(hash.G)
    colnames(mat) <- paste("BIN", sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")
    print(paste0("build mat end: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    ##
    print(paste0("Create Seurat Object start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    seurat_spatialObj <- CreateSeuratObject(
        mat,
        project = obj_name,
        assay = "Spatial",
        min.cells = 5,
        min.features = 5
    )
    print(paste0("Create Seurat Object end: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    ##
    print(paste0("generate spatial Object start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    spatialObj <- generate_spatialObj(
        image = tissue_lowres_image,
        scale.factors = fromJSON(scalefactors_json),
        tissue.positions = tissue_positions_list
    )
    print(paste0("generate spatial Object end: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    ##
    spatialObj <- spatialObj[Cells(seurat_spatialObj)]
    DefaultAssay(spatialObj) <- "Spatial"
    seurat_spatialObj[[obj_name]] <- spatialObj

    return(seurat_spatialObj)
}

qc_information <- function(seurat_spatialObj, dir_name, pro, bs, nFeature_threshold = 0, scale_factor = 1) {

    pro_name <- paste0(dir_name, "/", pro, "_bin_", bs, "_Feature_threshold_", nFeature_threshold)
    csv_name <- paste0(pro_name, "_meta_data.csv")
    pdf_name <- paste0(pro_name, "_figures.pdf")
    pdf_all_in_one <- paste0(pro_name, "_figures_all.pdf")
    note_message <- paste0(pro, "_bin", bs, "_Feature_threshold_", nFeature_threshold)

    seurat_spatialObj@images[[1]]@spot.radius <- seurat_spatialObj@images[[1]]@spot.radius * scale_factor

    Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
    Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
    upper <- as.numeric(Q3 + 1.5 * (Q3 - Q1))
    lower <- as.numeric(Q1 - 1.5 * (Q3 - Q1))
    mean_gene <- as.integer(mean(seurat_spatialObj$nFeature_Spatial))
    median_gene <- as.integer(median(seurat_spatialObj$nFeature_Spatial))
    mean_UMI <- as.integer(mean(seurat_spatialObj$nCount_Spatial))
    median_UMI <- as.integer(median(seurat_spatialObj$nCount_Spatial))
    mean_mito <- as.integer(mean(seurat_spatialObj$percent.mt))
    median_mito <- as.integer(median(seurat_spatialObj$percent.mt))
    bin_number <- dim(seurat_spatialObj)[2]

    preqc_meta_data <- tibble(
        sample = pro,
        bin_number = bin_number,
        mean_gene = mean_gene,
        median_gene = median_gene,
        mean_UMI = mean_UMI,
        median_UMI = median_UMI,
        mean_mito = mean_mito,
        median_mito = median_mito,
        Feature_threshold = nFeature_threshold
    )

    write.csv(preqc_meta_data, csv_name, quote = F, row.names = F)

    p1 <- VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3, pt.size = 0)

    p2 <- ggplot(seurat_spatialObj@meta.data, aes(x = nFeature_Spatial)) +
        geom_density(colour = "black") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold.italic"), 
            legend.position = "none", 
            axis.title = element_text(size = 15, face = "bold.italic"), 
            axis.text.x = element_text(size = 12), 
            axis.ticks.x = element_blank()) +
        geom_vline(aes(xintercept = 500, colour = "#999999", linetype = "twodash")) +
        geom_vline(aes(xintercept = lower, colour = "#377EB8", linetype = "twodash")) +
        geom_vline(aes(xintercept = upper, colour = "#E41A1C", linetype = "twodash")) +
        xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial), max(seurat_spatialObj@meta.data$nFeature_Spatial)) +
        ggtitle(paste0(
            pro, ".BIN_", bs, ":", "\n", 
            "nGene:", dim(seurat_spatialObj)[1], "; ",
            "nBIN:", dim(seurat_spatialObj)[2], "\n",
            "Gene:", mean_gene, " ", median_gene, "\n",
            "UMI:", mean_UMI, " ", median_UMI, "\n",
            "Mt ratio:", mean_mito, " ", median_mito
        ))
    
    xy_r <- GetTissueCoordinates(
            seurat_spatialObj,
            cols = c("row", "col"),
            scale = NULL
        )

    xy_r <- max(xy_r$col, xy_r$row)

    p3 <- SpatialFeaturePlot(seurat_spatialObj, features = "nFeature_Spatial") + theme(legend.position = "right") +
        scale_y_reverse() &
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r))
    
    p4 <- SpatialFeaturePlot(seurat_spatialObj, features = "nCount_Spatial") + theme(legend.position = "right") +
        scale_y_reverse() &
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r))
    
    pdf(pdf_name, width = 8, height = 8)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    dev.off()

    pdf(file = pdf_all_in_one, width = 20, height = 20)
    p_per <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    annotate_figure(p_per, top = text_grob(note_message, face = "bold", size = 20)) %>%
    print()
    dev.off()
}

qc_plot_png <- function(seurat_spatialObj, dir_name, pro, bs, nFeature_threshold = 0, scale_factor = 0.8) {

    pro_name <- paste0(dir_name, "/", pro, "_bin_", bs, "_Feature_threshold_", nFeature_threshold)
    png_all_in_one <- paste0(pro_name, "_figures_all.png")

    seurat_spatialObj@images[[1]]@spot.radius <- seurat_spatialObj@images[[1]]@spot.radius * scale_factor

    Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
    Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
    upper <- as.numeric(Q3 + 1.5 * (Q3 - Q1))
    lower <- as.numeric(Q1 - 1.5 * (Q3 - Q1))
    mean_gene <- as.integer(mean(seurat_spatialObj$nFeature_Spatial))
    median_gene <- as.integer(median(seurat_spatialObj$nFeature_Spatial))
    mean_UMI <- as.integer(mean(seurat_spatialObj$nCount_Spatial))
    median_UMI <- as.integer(median(seurat_spatialObj$nCount_Spatial))
    mean_mito <- as.integer(mean(seurat_spatialObj$percent.mt))
    median_mito <- as.integer(median(seurat_spatialObj$percent.mt))


    p1 <- VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3, pt.size = 0)

    p2 <- ggplot(seurat_spatialObj@meta.data, aes(x = nFeature_Spatial)) +
        geom_density(colour = "black") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold.italic"), 
            legend.position = "none", 
            axis.title = element_text(size = 15, face = "bold.italic"), 
            axis.text.x = element_text(size = 12), 
            axis.ticks.x = element_blank()) +
        geom_vline(aes(xintercept = 500, colour = "#999999", linetype = "twodash")) +
        geom_vline(aes(xintercept = lower, colour = "#377EB8", linetype = "twodash")) +
        geom_vline(aes(xintercept = upper, colour = "#E41A1C", linetype = "twodash")) +
        xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial), max(seurat_spatialObj@meta.data$nFeature_Spatial)) +
        ggtitle(paste0(
            pro, ".BIN_", bs, ":", "\n", 
            "nGene:", dim(seurat_spatialObj)[1], "; ",
            "nBIN:", dim(seurat_spatialObj)[2], "\n",
            "Gene:", mean_gene, " ", median_gene, "\n",
            "UMI:", mean_UMI, " ", median_UMI, "\n",
            "Mt ratio:", mean_mito, " ", median_mito
        ))

    
    xy_r <- GetTissueCoordinates(
        seurat_spatialObj,
        cols = c("row", "col"),
        scale = NULL
    )

    xy_r <- max(xy_r$col, xy_r$row)


    
    p3 <- SpatialFeaturePlot(seurat_spatialObj, features = "nFeature_Spatial") + theme(legend.position = "right") +
        scale_y_reverse() &
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r))
    
    p4 <- SpatialFeaturePlot(seurat_spatialObj, features = "nCount_Spatial") + theme(legend.position = "right") +
        scale_y_reverse() &
        coord_fixed(ratio = 1, xlim = c(1, xy_r), ylim = c(1, xy_r))

    ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2) +
    ggtitle(pro)

    ggsave(png_all_in_one, width = 20, height = 20)
}

##################################################################################################################

args <- commandArgs(T)
infile <- args[1]
bs <- as.numeric(args[2])
outdir <- args[3]
obj_name <- args[4]
nFeature_threshold <- as.numeric(args[5])

if (!dir.exists(outdir)) {
    dir.create(outdir)
}
setwd(outdir)

pro <- infile %>%
    strsplit("/") %>%
    unlist() %>%
    tail(1) %>%
    gsub(".bin50|.txt|.tsv|.gz|_filterd|.gem|_bin1", "", .)


# check perQC rds
preqc_rds <- paste0(obj_name, "_bin", bs, "_preQC.rds")
dir_name <- "qc_Feature_threshold_0"
if (file.exists(preqc_rds)) {
    seurat_spatialObj <- readRDS(preqc_rds)
    qc_plot_png(seurat_spatialObj, dir_name, obj_name, bs)
} else {
    dir.create(dir_name)
    
    print("creat Spatial Object")
    # 1. bin data
    # 2. creat Spatial Object
    seurat_spatialObj <- generate_seurat_spatialObj(infile, obj_name, bs, obj_name)
    gc()

    print("Spatial Analyse")
    #  3. Spatial Analyse
    seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^mt-")
    qc_information(seurat_spatialObj, dir_name, obj_name, bs)
    saveRDS(seurat_spatialObj, file = preqc_rds)
}

# check nFeature_threshold QC rds
dir_name <- paste0("qc_Feature_threshold_", nFeature_threshold)
after_qc_rds <- paste0(obj_name, "_bin", bs, "_Feature_threshold_", nFeature_threshold, "_qc.rds")
if (file.exists(after_qc_rds) || nFeature_threshold == 0 || is.na(nFeature_threshold)) {
    print(paste0("Dir exists::", dir_name, "!"))
    seurat_spatialObj <- subset(seurat_spatialObj, nFeature_Spatial > nFeature_threshold)
    qc_plot_png(seurat_spatialObj, dir_name, obj_name, bs, nFeature_threshold = nFeature_threshold)
} else {
    dir.create(dir_name)
    #  4. quality check
    seurat_spatialObj <- subset(seurat_spatialObj, nFeature_Spatial > nFeature_threshold)
    qc_information(seurat_spatialObj, dir_name, obj_name, bs, nFeature_threshold = nFeature_threshold)
    saveRDS(seurat_spatialObj, file = after_qc_rds)
}
