#' DE_within_each_cluster
#'
#' @param obj 
#' @param outDir 
#' @param clusterCol 
#' @param groupCol
#' @return
#' @export
#'
#' @examples
DE_within_each_cluster <- function(obj, outDir, clusterCol = "annotated_clusters", groupCol = "group", group1, group2, method = "DESeq2") {
  
  # add column
  obj$clusterCol <- obj[[clusterCol]]
  obj$groupCol <- obj[[groupCol]]
  
  # initialize variables
  cell_types <- levels(obj$clusterCol)
  if(is.null(cell_types)) {cell_types <- unique(obj$clusterCol)}
  master.df <- data.frame()
  
  # loop through clusters
  for (i in cell_types) {
    
    # print cell type you're on
    print(i)
    
    # subset object by cell type
    cluster <- subset(obj, clusterCol == i)
    
    # set idents based on group
    Idents(cluster) <- cluster$groupCol
    
    # differential expression
    markers <- FindMarkers(object = cluster,
                           ident.1 = group1,
                           ident.2 = group2,
                           only.pos = FALSE, # default
                           min.pct = 0.10,
                           test.use = method,
                           verbose = TRUE,
                           assay = "RNA")
    
    # check if no DEGs
    if(nrow(markers) == 0) {
      next
    }
    
    # add markers to master table
    markers$cluster <- i
    markers$gene <- rownames(markers)
    master.df <- rbind(master.df, markers)
  }
  
  # rename dynamically using rename_with()
  master.df <- master.df %>%
    rename_with(~paste0("percent_", c(group1, group2)), c(pct.1, pct.2))
  
  # set row names
  rownames(master.df) <- 1:nrow(master.df)
  
  # calculate difference
  master.df$percent_difference <- 
    abs(master.df[[paste0("percent_", group1)]] - master.df[[paste0("percent_", group2)]])
  
  # reorder columns
  master.df <- master.df[,c(6,7,1,5,2,3,4,8)]
  
  # check if the outDir ends with a slash
  if(!stringr::str_ends(outDir, "/")) {
    outDir <- paste0(outDir, "/")
  }
  
  # write table
  write.table(x = master.df, 
              file = paste0(outDir, group1, "_vs_", group2, "_DEGs.tsv"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}
