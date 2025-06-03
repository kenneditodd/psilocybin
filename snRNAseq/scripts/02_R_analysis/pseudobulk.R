# load libraries
library(Seurat)
library(dplyr)
library(DESeq2)

# Source thresholds and output paths
source("../../refs/thresholds_and_outs.R")
out <- paste0(out, "pass2/")

# Source custom functions
files <- list.files("../../functions", full.names = TRUE)
invisible(lapply(files, source))

# Load environment variables
load_dot_env(file = "../../refs/.env")
project_dir <- Sys.getenv("PROJECT_DIR")
de_method <- Sys.getenv("DE_METHOD", unset = "DESeq2")

# Read Seurat object
mouse.annotated <- readRDS(paste0(project_dir, "/rObjects/", filtering_method,
                                  "_pass2_annotated_seurat_obj.rds"))

# Step 1: Aggregate pseudobulk counts by cell type and sample

# Step 2: Get list of unique cell types
cell_types <- unique(mouse.annotated$annotated_clusters)

# Step 3: Split the pseudobulk matrix by cell type
split_pb_counts <- lapply(celltypes, function(ct) {
  cols <- grep(paste0("^", ct, "_"), colnames(pb_counts), value = TRUE)
  pb_counts[, cols, drop = FALSE]
})
names(split_pb_counts) <- celltypes

# Step 4: Extract metadata mapping sample_id to group
sample_meta <- unique(mouse.annotated@meta.data[, c("sample_id", "group", "group2", "treatment")])

# Step 6: Run pseudobulk DESeq2 for all cell types
de_results <- lapply(celltypes, function(ct) {
  run_pseudobulk_deseq2(ct, split_pb_counts, sample_meta, group1 = "L.24h", group2 = "S.24h")
})
names(de_results) <- celltypes

# Step 7: Combine all DE results into one data frame
all_de <- do.call(rbind, de_results)
