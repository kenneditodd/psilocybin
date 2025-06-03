# variables
filtering_method <- "thresh1_no_integration"
out <- paste0("../../results/", filtering_method, "/")

# thresholds
if (filtering_method == "thresh1_no_integration") {
  # Cell Ranger output
  nCount.min <- 5000
  nCount.max <- 50000
  nFeature.min <- 2000
  complexity.cutoff <- 0.8
  mt.cutoff <- 1
  hb.cutoff <- 1
} else if (filtering_method == "thresh1_cellbender") {
  # Cell Ranger output
  nCount.min <- 5000
  nCount.max <- 50000
  nFeature.min <- 2000
  complexity.cutoff <- 0.8
  mt.cutoff <- 1
  hb.cutoff <- 1
} 
