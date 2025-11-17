# Psilocybin snRNAseq
# Filter low quality cells. Remove lowly expressed or noisy genes.
# Kennedi Todd
# November 7, 2025

# load libraries
library(dotenv) # load_dot_env
setwd(".")

# source thresholds and output paths
source("../../refs/thresholds_and_outs.R")
out <- paste0(out, "pass1/")

# load the environment variables
load_dot_env(file = "../../refs/.env")
ann_dir <- Sys.getenv("ANNOTATION_REFS")

# read in annotation file
genes <- readRDS("../../rObjects/annotation.rds")

# load data
mouse <- readRDS("../../rObjects/cellranger_filtered_seurat_obj.rds")

# filter
mouse_filtered <- subset(mouse, subset = (nCount_RNA > nCount.min) &
                           (nCount_RNA < nCount.max) &
                           (nFeature_RNA > nFeature.min) &
                           (cell_complexity > complexity.cutoff) &
                           (percent_mt < mt.cutoff) &
                           (percent_hb < hb.cutoff) &
                           (percent_nuclear > nuclear.cutoff))

# print cells removed
print(paste0(dim(mouse)[2] - dim(mouse_filtered)[2]," cells removed"))

# filter lowly expressed genes
counts <- GetAssayData(object = mouse_filtered, layer = "counts")
nonzero <- counts > 0  # produces logical
keep <- Matrix::rowSums(nonzero) >= 10  # sum the true/false
counts_filtered <- counts[keep,]  # keep certain genes

# print features removed
print(paste0(dim(counts)[1] - dim(counts_filtered)[1], " lowly expressed features removed"))

# get mt.genes
gene.names <- genes[genes$type == "gene",]
gene.names <- gene.names[gene.names$gene_type == "protein_coding",]
mt.genes <- gene.names[gene.names$seqnames == "chrM", "gene_name"]

# remove mt.genes
keep <- !rownames(counts_filtered) %in% mt.genes # false when mt.gene
counts_filtered2 <- counts_filtered[keep,]

# create new obj
options(Seurat.object.assay.calcn = TRUE)
mouse_filtered <- CreateSeuratObject(counts_filtered2)
mouse_filtered <- AddMetaData(mouse_filtered, metadata = mouse@meta.data[colnames(mouse_filtered), c(1,4:25)])

# print features removed
print(paste0(dim(counts)[1] - dim(counts_filtered2)[1], " features removed"))

# cleanup data
remove(mouse,counts,counts_filtered,counts_filtered2,nonzero,genes,gene.names)

# save
saveRDS(object = mouse_filtered, 
        file = paste0( "../../rObjects/pass1_filtered_seurat_obj.rds"), 
        compress = FALSE)
