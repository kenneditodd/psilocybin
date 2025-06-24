DESeq2_per_cell_type <- function(cell_type, pb_counts, sample_meta, group_column, group1, group2, min_counts = 20, min_samples = 3) {
  message("==========")
  message(paste0("Running DE for ", group1, " vs ", group2, " in ", cell_type, " (", group_column, ")"))
  
  # Get sample IDs from count matrix
  sample_ids <- gsub(paste0(cell_type, "_"), "", colnames(pb_counts))
  
  # Subset and align metadata
  meta <- sample_meta %>%
    filter(sample_id %in% sample_ids) %>%
    arrange(match(sample_id, sample_ids))
  
  # Subset count matrix to those samples
  count_cols <- paste0(cell_type, "_", meta$sample_id)
  counts <- pb_counts[, count_cols, drop = FALSE]
  
  # Subset genes: keep if expressed in ≥ min_samples with ≥ min_counts
  keep_genes <- rowSums(counts >= min_counts) >= min_samples
  counts <- counts[keep_genes, ]
  
  # Skip if too few genes remain
  if (nrow(counts) < 10) {
    message("⚠️ Not enough genes after filtering. Skipping ", cell_type)
    return(NULL)
  }
  
  # Create condition factor
  meta$condition <- factor(meta[[group_column]], levels = c(group2, group1))
  if (length(unique(meta$condition)) < 2) {
    message("❌ Only one group present in filtered samples. Skipping DE.")
    return(NULL)
  }
  
  rownames(meta) <- colnames(counts)
  
  # Run DESeq2
  message("Creating DESeq2 object...")
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~condition)
  message("Running DESeq...")
  dds <- DESeq(dds)
  
  message("Extracting results...")
  res <- results(dds, contrast = c("condition", group1, group2))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res$cell_type <- cell_type
  res$comparison <- paste0(group1, "_vs_", group2)
  res$group_by <- group_column
  
  message("✅ Done with ", cell_type)
  return(res)
}
