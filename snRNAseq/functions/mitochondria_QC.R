#' mitochondria_QC
#'
#' @param obj 
#' @param outDir 
#' @param verbose 
#' @param mitoCol
#' @return
#' @export
#'
#' @examples
mitochondria_QC <- function(obj, outDir = NULL, mitoCol = NULL, sampleCol = "sample", verbose = TRUE) {
  
  # required packages
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # check class
  if(!is(obj,"Seurat")) {stop("obj argument is not a Seurat object.")}
  if(!is(outDir,"character")) {stop("outDir argument is not a character.")}
  if(!is(mitoCol,"character") && !is(mitoCol, "NULL")) { stop("mitoCol argument is not a character")}
  if(!is(sampleCol,"character")){stop("sampleCol is not a character")}
  if(!is(verbose,"logical")) { stop("verbose argument is not a logical")}
  
  # set output directory
  output <- "./"
  if(!is.null(outDir)) {output <- outDir}
  if(!endsWith(output,"/")) { output <- paste0(output,"/")}
  
  # Raw counts are not comparable between cells. Each cell has a different nCount_RNA. 
  # The log normalization, ((nCount_RNA / nFeature_RNA) * log1p, is taken in order to explore variation
  if(verbose){print("Normalizing data")}
  obj <- NormalizeData(obj)
  
  # find variable features
  if(verbose) {print("Finding variable features")}
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  
  # scale
  if(verbose){print("Scaling data")}
  obj <- ScaleData(obj)
  
  # run PCA
  if(verbose){print("Running PCA")}
  obj <- RunPCA(obj)
  
  # set percent.mt column name if named alternatively
  mito_column <- "percent_mt"
  if(!is.null(mitoCol)) {mito_column <- mitoCol}
  
  # set quartile values
  meta <- obj@meta.data
  first <- as.numeric(summary(meta[,mito_column])[2])
  mean <- as.numeric(summary(meta[,mito_column])[4])
  third <- as.numeric(summary(meta[,mito_column])[5])
  
  # turn percent.mt into factor based on quartile value
  if(verbose){print("Factoring percent.mt")}
  obj[["mito_factor"]] <- cut(meta[,mito_column],
                              breaks = c(-Inf, first, mean, third, Inf),
                              labels = c("Low", "Medium", "Medium high", "High"))
  
  # plot and save PCA
  if(verbose){print("Outputing plots")}
  pdf(paste0(output, "mitocondria_mito_factor_pca.pdf"), height = 4, width = 6)
  pca <- DimPlot(obj,
                 reduction = "pca",
                 group.by = "mito_factor",
                 split.by = "mito_factor")
  print(pca)
  dev.off()
  
  # update metadata with sample column name given
  meta <- dplyr::rename(obj@meta.data, mito_sample_split=sampleCol)
  obj@meta.data <- meta
  
  # percent cells per mito.factor per sample plot
  pdf(paste0(output, "mitochondria_percent_mito_factor_per_sample.pdf"), height = 4, width = 6)
  percent.mito.factor <- obj@meta.data %>%
    group_by(mito_sample_split, mito_factor) %>%
    dplyr::count() %>%
    group_by(mito_sample_split) %>%
    dplyr::mutate(percent = 100*n/sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = mito_sample_split, y = percent, fill = mito_factor)) +
    geom_col() +
    ggtitle("Percentage of mito_factor per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Sample")
  print(percent.mito.factor)
  dev.off()
  
  # return phase vector
  return(obj$mito_factor)
}
