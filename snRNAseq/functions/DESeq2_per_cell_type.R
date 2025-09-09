#' Differential expression for pseudobulk data within a cell type
#'
#' Run DESeq2 on pseudobulk counts for a given cell type and comparison.
#'
#' @param cell_type Character. Name of the cell type being tested.
#' @param pb_counts Matrix. Pseudobulk count matrix (genes x samples), already subset to the current cell type and comparison samples (columns like "<cell_type>_<sample_id>").
#' @param sample_meta Data frame. Metadata for samples, must include `sample_id` and the grouping column specified in `group_column`. Should already be filtered to this comparison.
#' @param group_column Character. Column in `sample_meta` defining groups (e.g. "group2").
#' @param group1 Character. Name of the first group (numerator in the contrast).
#' @param group2 Character. Name of the second group (reference in the contrast).
#' @param min_counts Integer. Minimum per-sample count threshold to consider a gene expressed (default 20).
#' @param samples_per_group Integer or NULL. If NULL (default), automatically set to the size of the smaller group in this comparison. The filter then keeps genes with counts ≥ `min_counts` in at least `samples_per_group` samples total (across both groups combined).
#'
#' @return Data frame of DESeq2 results with added columns: gene, cell_type, comparison, group_by.
#'   Returns NULL if too few genes remain after filtering or only one group is present.
#'
#' @examples
#' \dontrun{
#' res <- DESeq2_per_cell_type(
#'   cell_type    = "Astrocytes",
#'   pb_counts    = pb_mat_astros,   # genes x samples for this cell type & comparison
#'   sample_meta  = meta_cmp,
#'   group_column = "group2",
#'   group1       = "H.24h.F",
#'   group2       = "S.24h.F",
#'   min_counts   = 20,
#'   samples_per_group = NULL
#' )
#' }
DESeq2_per_cell_type <- function(cell_type, pb_counts, sample_meta, group_column, group1, group2, min_counts = 20, samples_per_group = NULL) {
  
  message("==========")
  message(paste0("Running DE for ", group1, " vs ", group2, " in ", cell_type, " (", group_column, ")"))
  
  # Align meta rows with pb_counts columns
  sample_ids <- gsub(paste0("^", cell_type, "_"), "", colnames(pb_counts))
  meta <- sample_meta
  rownames(meta) <- paste0(ct, "_", rownames(meta))
  meta <- meta[colnames(pb_counts),]
  
  # Check
  stopifnot( all.equal(colnames(pb_counts), rownames(meta)))
  
  # Create condition factor for DESeq2
  meta$condition <- factor(meta[[group_column]], levels = c(group2, group1))
  if (length(unique(meta$condition)) < 2) {
    message("Only one group present in filtered samples. Skipping DE.")
    return(NULL)
  }
  
  # Auto-calculate samples_per_group if not provided: size of the smaller group
  if (is.null(samples_per_group)) {
    samples_per_group <- min(
      sum(meta$condition == group1),
      sum(meta$condition == group2)
    )
  }
  
  # Filter genes: keep if expressed (>= min_counts) in >= samples_per_group samples (total)
  keep_genes <- rowSums(pb_counts >= min_counts) >= samples_per_group
  pb_counts  <- pb_counts[keep_genes, , drop = FALSE]
  
  # Run DESeq2
  message("Creating DESeq2 object...")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb_counts, colData = meta, design = ~ condition)
  
  message("Running DESeq...")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  message("Extracting results...")
  res <- DESeq2::results(
    object = dds, 
    contrast = c("condition", group1, group2))
  res <- as.data.frame(res)
  res$gene       <- rownames(res)
  res$cell_type  <- cell_type
  res$comparison <- paste0(group1, "_vs_", group2)
  res$group_by   <- group_column
  
  message("Done with ", cell_type)
  return(res)
}
