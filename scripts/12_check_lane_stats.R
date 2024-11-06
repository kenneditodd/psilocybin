setwd("/research/labs/neurology/fryer/m214960/psi1")

# read and reformat lane 3
lane3 <- read.delim2("rawQC/lane3_multiqc_data/multiqc_general_stats.txt")
colnames(lane3) <- gsub("FastQC_mqc.generalstats.fastqc.","", colnames(lane3))
rownames(lane3) <- lane3$Sample
lane3$Sample <- NULL
lane3$avg_sequence_length <- NULL
lane3$median_sequence_length <- NULL
lane3[] <- lapply(lane3, as.numeric)
summary(lane3)

# read and reformat lane4
lane4 <- read.delim2("rawQC/lane4_multiqc_data/multiqc_general_stats.txt")
colnames(lane4) <- gsub("FastQC_mqc.generalstats.fastqc.","", colnames(lane4))
rownames(lane4) <- lane4$Sample
lane4$Sample <- NULL
lane4$avg_sequence_length <- NULL
lane4$median_sequence_length <- NULL
lane4[] <- lapply(lane4, as.numeric)
summary(lane4)

# T-test
# Assuming lane3 and lane4 have the same column names
results <- lapply(1:4, function(i) {
  t.test(lane3[[i]], lane4[[i]], var.equal = TRUE) # Set var.equal=FALSE if variances are unequal
})

# Name the results by column names for easy reference
names(results) <- colnames(lane3)

# Extract p-values for each T-test
p_values <- sapply(results, function(x) x$p.value)
p_values

adjusted_p_values <- p.adjust(p_values, method = "BH")
adjusted_p_values

