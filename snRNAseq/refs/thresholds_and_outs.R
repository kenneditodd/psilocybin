# variables
filtering_method <- "thresh1_soupx"
out <- paste0("../../results/", filtering_method, "/")

# thresholds
if (filtering_method == "thresh1_soupx") {
  # Cell Ranger output
  nCount.min <- 7000
  nCount.max <- 80000
  nFeature.min <- 2000
  complexity.cutoff <- 0.78
  mt.cutoff <- 1
  hb.cutoff <- 1
  nuclear.cutoff <- 60
} else if (filtering_method == "thresh2_cellbender") {
  # cellbender output
  nCount.min <- 7000
  nCount.max <- 80000
  nFeature.min <- 2000
  complexity.cutoff <- 0.8
  mt.cutoff <- 1
  hb.cutoff <- 1
}

