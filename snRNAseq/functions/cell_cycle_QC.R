#' cell_cycle_QC
#'
#' @param obj A v5 Seurat object.
#' @param markersPath Path to table with cell markers.
#' @param species A character with value "mouse" or "human". Necessary for subsetting cell cycle phase marker table.
#' @param outDir Path to output directory where QC plots will be saved. The default directory is the current one.
#' @param verbose Display messages
#' @return A character vector containing cell phases.
#' @export
#' @examples
#' #' mySeuratObject[["phase"]] <- cell_cycle_QC(mySeuratObject, "mouse")
cell_cycle_QC <- function(obj, markersPath, species, sampleCol = "sample", outDir = NULL, verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object")}
  if(!species %in% c("mouse","human")){stop("species argument is not a valid option")}
  if(!is(markersPath,"character")){stop("markers is not a character")}
  if(!is(sampleCol,"character")){stop("sampleCol is not a character")}
  if(!is(verbose,"logical")) { stop("verbose argument is not a logical")}
  
  # set output dir
  output <- "./"
  if(!is.null(outDir)) {output <- outDir}
  if(!endsWith(output,"/")) { output <- paste0(output,"/")}
  
  # Raw counts are not comparable between cells. Each cell has a different nCount_RNA. 
  # The log normalization, ((nCount_RNA / nFeature_RNA) * log1p, is taken in order to explore variation
  if(verbose){print("Normalizing data")}
  obj <- NormalizeData(obj)
  
  # load cell cycle phase markers
  phase.markers <- read.delim(markersPath, header = TRUE, sep = "\t")
  phase.markers <- phase.markers[phase.markers$species == species,]
  
  # save the file
  write.table(phase.markers, 
              paste0(output, species, "_cell_cycle_phase_markers.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # subset cell cycle markers
  g2m <- phase.markers[phase.markers$phase == "G2/M", "gene_name"]
  s <- phase.markers[phase.markers$phase == "S", "gene_name"]
  
  # score cells
  if(verbose){print("Cell cycle scoring")}
  obj <- CellCycleScoring(obj,
                          g2m.features = g2m,
                          s.features = s,
                          set.ident = TRUE)
  
  # find variable features
  if(verbose) {print("Finding variable features")}
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  
  # scale
  if(verbose){print("Scaling data")}
  obj <- ScaleData(obj)
  
  # run PCA
  if(verbose){print("Running PCA")}
  obj <- RunPCA(obj)
  
  # plot and save PCA plot
  if(verbose){print("Outputing plots")}
  pdf(paste0(output, "cell_cycle_pca.pdf"), height = 4, width = 6)
  pca <- DimPlot(obj,
                 reduction = "pca",
                 group.by = "Phase",
                 split.by = "Phase")
  print(pca)
  dev.off()
  
  # update metadata with sample column name given
  meta <- dplyr::rename(obj@meta.data, cell_cycle_sample_split=sampleCol)
  obj@meta.data <- meta
  
  # percent cells per phase plot
  pdf(paste0(output, "cell_cycle_percent_phase_per_sample.pdf"), height = 4, width = 6)
  percent.phase <- obj@meta.data %>%
    group_by(cell_cycle_sample_split, Phase) %>%
    dplyr::count() %>%
    group_by(cell_cycle_sample_split) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = cell_cycle_sample_split, y = percent, fill = Phase)) +
    geom_col() +
    ggtitle("Percentage of phase per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Sample")
  print(percent.phase)
  dev.off()
  
  # return phase vector
  return(obj$Phase)
}
