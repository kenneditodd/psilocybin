# load libraries
library(dotenv)      # load_dot_env()
library(dplyr)       # left_join()
library(rtracklayer) # import()
library(scCustomize) # Read_CellBender_h5_Mat
library(Seurat)      # CreateSeuratObject()
library(stringr)     # str_match()

# get sample names
samples <- readLines("../../refs/sample_list.tsv")
samples <- gtools::mixedsort(samples)

# read in meta
meta <- readRDS("../../rObjects/meta.rds")

# read counts and create seurat obj
seurat_obj_list <- list()
if (file.exists(paste0("../../rObjects/cellbender_obj_merged.rds"))) {
  mouse <- readRDS(paste0("../../rObjects/cellbender_obj_merged.rds"))
} else {
  
  # path info
  prefix <- "../../counts/"
  suffix <- "_cellbender_filtered.h5"
  
  # create list of individual seurat objects
  for (i in 1:length(samples)) {
    print(i)
    sample <- samples[i]
    
    # Create Seurat object with PIPseeker output
    obj <- CreateSeuratObject(
      Read_CellBender_h5_Mat(paste0(prefix, sample, "/outs/", sample, suffix))
    )
    
    # Add sample ID as prefix to cell names
    obj <- RenameCells(obj, add.cell.id = sample)
    
    # Add Seurat object to the list with the sample name as the key
    seurat_obj_list[[sample]] <- obj
    
    # cleanup - helps with memory
    remove(obj)
    gc()
  }
  
  # Merge all Seurat objects
  mouse <- merge(seurat_obj_list[[1]], 
                 y = seurat_obj_list[-1])
  
  # Set project name
  mouse@project.name <- "E.coli Mice scRNAseq"
  
  # Join layers
  mouse$orig.ident <- colnames(mouse)
  mouse <- JoinLayers(mouse)
  
  # Extract animal_id
  mouse$animal_id <- str_match(colnames(mouse), "[oldyung]+_[fesm]+_([0-9]+)_.+")[,2]
  
  # Check
  table(mouse$animal_id)
  
  # Add meta
  mouse@meta.data <- left_join(x = mouse@meta.data, y = meta, by = "animal_id")
  rownames(mouse@meta.data) <- mouse$orig.ident
  
  # save
  saveRDS(mouse, "../../rObjects/cellbender_obj_merged.rds")
  
  # cleanup
  remove(meta, seurat_obj_list)
  gc()
  
} # end of else statement

# preview
mouse

# extract barcodes that survived cell bender filtering
barcodes <- colnames(mouse)
saveRDS(barcodes, file = "../../rObjects/cellbender_passed_barcodes.rds", compress = FALSE)
