# Kennedi Todd
# December 29, 2025
# Sex-specific model

# Load libraries
library(dotenv)
library(dplyr)
library(rtracklayer)
library(sleuth)

# Load environment variables
print("Sex-specific analysis")
out <- "../results/"
load_dot_env("../refs/.env")

# Read annotation file
annotation_rds <- "../rObjects/annotation.rds"
if (file.exists(annotation_rds)) {
  transcripts <- readRDS(annotation_rds)
} else {
  transcripts <- rtracklayer::import(Sys.getenv("MMUSCULUS_GTF")) %>%
    as.data.frame() %>%
    filter(type == "transcript") %>%
    dplyr::rename("target_id" = "transcript_id")
  saveRDS(transcripts, annotation_rds)
}

# Read meta
meta <- read.delim2(file = "../refs/metadata.tsv")
meta$notes <- NULL

# Add kallisto quant output dir path
meta$path <- paste0(
  Sys.getenv("PROJECT_DIR"), "/kallisto/", meta$filename, 
  c(rep("_L3_", 45), rep("_L4_", 45)),"kallisto"
)

# Rename sample column
meta$sample <- paste0(meta$filename, c(rep("_L3_", 45), rep("_L4_", 45)),"kallisto")

# Set factors
meta$group <- factor(meta$group)
meta$sex <- factor(meta$sex)
meta$group2 <- factor(meta$group2)

# Create sleuth object
# Loads your Kallisto data, merges metadata, filters low-expressed transcripts, 
# normalizes expression, summarizes uncertainty from bootstraps, attaches annotation, 
# transforms the data, and builds a single object that Sleuth can model.
so <- sleuth_prep(
  sample_to_covariates = meta, 
  full_model = ~ group2,
  num_cores = 12,
  read_bootstrap_tpm = TRUE,
  target_mapping = transcripts)

# transcripts that passed filtering
sleuth_keep <- so$filter_df$target_id

# Extract Hbb-bs TPM and add to meta
# est_counts = Kallisto’s estimated number of fragments assigned to this transcript
# obs_norm = est_counts but with but normalized across samples using DESeq-style size factors
# eff_len = effective transcript length used by Kallisto to compute abundances adjusted for fragment length distribution
tpm <- so$obs_norm_filt %>% 
  left_join(so$target_mapping, by = "target_id")
hbb_log2_tpm <- tpm %>%
  filter(target_id == "ENSMUST00000023934") %>%
  mutate(hbb_log2_tpm = log2(tpm + 0.5)) %>%
  select(sample, hbb_log2_tpm)
meta <- meta %>% left_join(hbb_log2_tpm, by = "sample")

# Filter by transcript_type and remove MT genes
custom_keep <- transcripts %>%
  dplyr::filter(
    transcript_type %in% c("protein_coding",
                           "retained_intron",
                           "nonsense_mediated_decay"),
    seqnames != "chrM"
  ) %>%
  dplyr::pull(target_id)
final_keep <- intersect(sleuth_keep, custom_keep)

# Re-create sleuth object with new design
so <- sleuth_prep(
  sample_to_covariates = meta, 
  full_model = ~ group2 + hbb_log2_tpm,
  num_cores = 12,
  read_bootstrap_tpm = TRUE,
  target_mapping = transcripts,
  filter_target_id = final_keep
)

# Set comparisons
myContrasts <- list(
  c("L.8h.F", "S.8h.F"), 
  c("L.8h.M", "S.8h.M"),
  c("L.24h.F", "S.24h.F"),
  c("L.24h.M", "S.24h.M"),
  c("L.7d.F", "S.7d.F"),
  c("L.7d.M", "S.7d.M"),
  c("H.8h.F", "S.8h.F"),
  c("H.8h.M", "S.8h.M"),
  c("H.24h.F", "S.24h.F"),
  c("H.24h.M", "S.24h.M"),
  c("H.7d.F", "S.7d.F"),
  c("H.7d.M", "S.7d.M"),
  c("S.8h.F", "S.8h.M"),
  c("S.24h.F", "S.24h.M"),
  c("S.7d.F", "S.7d.M")
)

# Differential expression loop
orig_groups <- so$sample_to_covariates$group2

for (ctr in myContrasts) {
  num <- ctr[1]   # numerator
  ref <- ctr[2]   # reference
  
  message("Running contrast: ", num, " vs ", ref)
  
  # Reset group2 to original levels each time
  so$sample_to_covariates$group2 <- orig_groups
  so$sample_to_covariates$group2 <- relevel(so$sample_to_covariates$group2, ref = ref)
  
  # Give model a unique name
  fit_name  <- paste0("full_", num, "_vs_", ref)
  so <- sleuth_fit(so, ~ group2 + hbb_log2_tpm, fit_name = fit_name)
  
  # Construct expected beta name
  design_cols <- colnames(models(so)$full$design_matrix)
  beta_name <- paste0("group2", num)
  
  # Wald test for this contrast
  so <- sleuth_wt(so, which_beta = beta_name, which_model = fit_name)
  
  # Extract results
  res <- sleuth_results(so, test = beta_name, test_type = "wt", which_model = fit_name)
  
  # Save
  out_file <- paste0(out, "DEG_tables/", num, "_vs_", ref, "_DEGs.tsv")
  write.table(res, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  message("Saved: ", out_file)
}

saveRDS(object = so, file = "../rObjects/so_group2_hbb.rds")
