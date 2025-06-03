# Psilocybin snRNAseq
# Create merged seurat object from cellranger output
# Kennedi Todd
# 05/21/2025

# load libraries
library(dotenv)      # load_dot_env()
library(dplyr)       # left_join()
library(rtracklayer) # import()
library(scCustomize) # Merge_Seurat_List()
library(Seurat)      # CreateSeuratObject()
library(stringr)     # str_match()

# load the environment variables
load_dot_env(file = "../../refs/.env")
ann_dir <- Sys.getenv("ANNOTATION_REFS")

# get sample names
samples <- readLines("../../refs/sample_list.tsv")
samples <- gtools::mixedsort(samples)

# read in annotation file
if (file.exists("../../rObjects/annotation.rds")) {
  genes <- readRDS("../../rObjects/annotation.rds")
} else {
  file <- paste0(ann_dir, "/refdata-gex-GRCm39-2024-A/genes/genes.gtf")
  genes <- rtracklayer::import(file)
  genes <- as.data.frame(genes)
  saveRDS(genes, "../../rObjects/annotation.rds")
}

# read in meta
meta <- readRDS("../../rObjects/meta.rds")

# read counts and create seurat obj
prefix <- "../../counts/"
suffix <- "/outs/filtered_feature_bc_matrix.h5"

# Initialize an empty list
seurat_obj_list <- list()
  
# create list of individual seurat objects
for (i in 1:length(samples)) {
  print(i)
  sample <- samples[i]
    
  # Create Seurat object with PIPseeker output
  obj <- CreateSeuratObject(
    Read10X_h5(paste0(prefix, sample, suffix))
    )
  
  # Get sample ID
  id <- meta[meta$filename == sample,]
  id <- id$sample_id
  
  # Add sample ID as prefix to cell names
  obj <- RenameCells(obj, add.cell.id = id)
    
  # Add Seurat object to the list with the sample name as the key
  seurat_obj_list[[id]] <- obj
    
  # cleanup - helps with memory
  remove(obj)
  gc()
}
  
# Merge all Seurat objects
mouse <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])
  
# Set project name
mouse@project.name <- "Psilocybin Project 1 snRNAseq"
  
# Join layers
mouse$orig.ident <- colnames(mouse)
mouse <- JoinLayers(mouse)
  
# Extract sample_id, H.24h.M.31_AAACCAAAGCATCCGG-1
mouse$sample_id <- str_match(colnames(mouse), "(.+)_[ACTG]")[,2]
  
# Check
table(mouse$sample_id)
  
# Add meta
mouse@meta.data <- left_join(x = mouse@meta.data, y = meta, by = "sample_id")
rownames(mouse@meta.data) <- mouse$orig.ident

# Check
all.equal(colnames(mouse), rownames(mouse@meta.data))
  
# cleanup
remove(meta, seurat_obj_list)
gc()

# preview
mouse

# set idents
Idents(mouse) <- "sample"

# add meta columns
# cell.complexity
mouse$cell_complexity <- log10(mouse$nFeature_RNA) / log10(mouse$nCount_RNA)
# percent.mt
mt.genes <- genes[genes$type == "gene",]
mt.genes <- mt.genes[mt.genes$gene_type == "protein_coding",]
mt.genes <- mt.genes[mt.genes$seqnames == "chrM",]
mt.genes <- mt.genes$gene_name
mouse$percent_mt <- PercentageFeatureSet(mouse, features = mt.genes)
mt.genes
# percent.ribo
# ribosomal proteins begin with 'Rps' or 'Rpl' in this annotation file
# mitochondrial ribosomes start with 'Mrps' or 'Mrpl'
ribo <- genes[genes$type == "gene",]
ribo <- ribo[ribo$gene_type == "protein_coding",]
mt.ribo <- ribo[grep("^Mrp[sl]", ribo$gene_name),"gene_name"]
ribo <- ribo[grep("^Rp[sl]", ribo$gene_name), "gene_name"]
ribo.combined <- c(mt.ribo,ribo)
ribo.combined <- ribo.combined[ribo.combined %in% rownames(mouse)]
mouse$percent_ribo <- PercentageFeatureSet(mouse, features = ribo.combined)
ribo.combined
# percent.hb
hb.genes <- c("Hba-x","Hba-a1","Hba-a2","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
mouse$percent_hb <- PercentageFeatureSet(mouse, features = hb.genes)
hb.genes

# set levels
mouse$group <- factor(mouse$group, levels = c("S.24h", "L.24h", "H.24h"))
mouse$group2 <- factor(mouse$group2, levels = c("S.24h.F", "S.24h.M", 
                                                "L.24h.F", "L.24h.M", 
                                                "H.24h.F", "H.24h.M"))
mouse$treatment <- factor(mouse$treatment, levels = c("saline","psilocybin"))

# set sample levels
new_order <- mouse@meta.data %>%
  arrange(group2) %>%
  select(sample_id)
new_order <- unique(new_order$sample_id)
mouse$sample_id <- factor(mouse$sample_id, levels = new_order)

# save with compression
saveRDS(mouse, "../../rObjects/cellranger_filtered_seurat_obj.rds", compress = FALSE)
