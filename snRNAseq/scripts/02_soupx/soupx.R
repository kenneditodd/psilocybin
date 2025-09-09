library(Seurat)
library(SoupX)

sc = SoupChannel(table_of_droplets,table_of_counts)
sc = setClusters(sc,cluster_labels)
sc = autoEstCont(sc)
out = adjustCounts(sc)
