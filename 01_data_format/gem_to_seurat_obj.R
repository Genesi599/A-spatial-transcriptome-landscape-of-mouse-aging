
# bin1 gem to seurat rds 
# filted bin50 gem + annotation to seurat rds

rm(list = ls())
gc()

#
library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(Matrix)
library(ggplot2)
library(ggpubr)

#
args <- commandArgs(T)
infile <- args[1] # input gem file, XXX.gem.gz / XXX.gem
bs <- as.numeric(args[2]) # bin size, if use for rebuilding filtered bin50 gem, bs should be 1
outdir <- args[3] # output dir
obj_name <- args[4] # output file name
nFeature_threshold <- as.numeric(args[5]) # qc threshold, if use for rebuilding filtered bin50 gem, nFeature_threshold should be 0 
anno_csv <- args[6] # annotation file path, only use for rebuilding filtered bin50 gem

##################################################################################################################

bin_data_progress <- function(infile, pro, bs) {
    dat <- fread(file = infile)
    print(paste0("load file: ", infile)) 

    if (length(grep("MIDCounts", colnames(dat)) > 0)) {
        dat <- dat %>%
            rename(MIDCount = MIDCounts)
    }

    if(bs != 1) {
        # (x, y) to bin(x, y)
        dat[, binx := trunc((x - min(x)) / bs + 1)]
        dat[, biny := trunc((y - min(y)) / bs + 1)]
        dat[, x := binx][, binx := NULL]
        dat[, y := biny][, biny := NULL]
    }

    #
    dat <- dat[, sum(MIDCount), by = .(geneID, x, y)][, bin_ID := max(x) * (y - 1) + x]
    return(dat)
}

generate_spatialObj <- function(image, scale.factors, tissue.positions) {

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

generate_seurat_spatialObj <- function(infile, pro, bs, obj_name) {

    dat <- bin_data_progress(infile, pro, bs)
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
        row = bin.coor$y,
        col = bin.coor$x,
        imagerow = bin.coor$y,
        imagecol = bin.coor$x
    )

    ##
    mat <- sparseMatrix(i = gen, j = bin, x = cnt)
    rownames(mat) <- rownames(hash.G)
    colnames(mat) <- paste("BIN", sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")

    ##
    seurat_spatialObj <- CreateSeuratObject(
        mat,
        project = obj_name,
        assay = "Spatial",
        min.cells = 5,
        min.features = 5
    )

    ##
    spatialObj <- generate_spatialObj(
        image = tissue_lowres_image,
        scale.factors = list(
            fiducial_diameter_fullres = 1,
            tissue_hires_scalef = 1,
            tissue_lowres_scalef = 1
        ),
        tissue.positions = tissue_positions_list
    )

    ##
    spatialObj <- spatialObj[Cells(seurat_spatialObj)]
    DefaultAssay(spatialObj) <- "Spatial"
    seurat_spatialObj[[obj_name]] <- spatialObj

    return(seurat_spatialObj)
}

add_ct_anno <- function(seurat_spatialObj, anno_csv) {

    seurat_spatialObj@meta.data <- seurat_spatialObj@meta.data %>%
        as_tibble(rownames = "cellid") %>%
        # get Coordinates
        left_join(
            GetTissueCoordinates(
                seurat_spatialObj,
                cols = c("row", "col"),
                scale = NULL
            ) %>%
            as_tibble(rownames = "cellid")
        )  %>%
        left_join(
            fread(anno_csv) %>%
                rename(row = y, col = x) %>%
                select(row, col, celltype)
        ) %>%
        select(-c(row, col)) %>%
        column_to_rownames("cellid")

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

##################################################################################################################

#
dir.create(outdir)
setwd(outdir)

# check perQC rds
preqc_rds <- paste0(obj_name, "_bin", bs, "_preQC.rds")
dir_name <- "qc_Feature_threshold_0"
if (grepl("_filtered_bin50.gem", basename(infile))) {
    bs <- 1
    preqc_rds <- paste0(obj_name, "_filtered_bin50.rds")
    dir_name <- "filtered_bin50_qc"
    nFeature_threshold <- 0
}

if (file.exists(preqc_rds)) {
    seurat_spatialObj <- readRDS(preqc_rds)
} else {
    dir.create(dir_name)
    
    print("creat Spatial Object")
    # 1. bin data
    # 2. creat Spatial Object
    seurat_spatialObj <- generate_seurat_spatialObj(infile, obj_name, bs, obj_name)
    gc()

    # add anno
    if (grepl("_filtered_bin50.gem", basename(infile)) & file.exists(anno_csv)) {
        print("add celltype annotation to Spatial Object")
        print(paste0("load file: ", anno_csv))

        seurat_spatialObj <- add_ct_anno(seurat_spatialObj, anno_csv)
    }

    print("Spatial Analyse")
    #  3. Spatial Analyse
    seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^mt-")
    qc_information(seurat_spatialObj, dir_name, obj_name, bs)
    saveRDS(seurat_spatialObj, file = preqc_rds)
}

#############################################
# check nFeature_threshold QC rds
dir_name <- paste0("qc_Feature_threshold_", nFeature_threshold)
after_qc_rds <- paste0(obj_name, "_bin", bs, "_Feature_threshold_", nFeature_threshold, "_qc.rds")

if (file.exists(after_qc_rds) || nFeature_threshold == 0 || is.na(nFeature_threshold)) {
    print(paste0("Dir exists::", dir_name, "!"))
} else {
    dir.create(dir_name)
    #  4. quality check
    seurat_spatialObj <- subset(seurat_spatialObj, nFeature_Spatial > nFeature_threshold)
    qc_information(seurat_spatialObj, dir_name, obj_name, bs, nFeature_threshold = nFeature_threshold)
    saveRDS(seurat_spatialObj, file = after_qc_rds)
}
