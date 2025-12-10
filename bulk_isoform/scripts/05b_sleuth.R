# Kennedi Todd
# November 18, 2025
# Sex-specific model

# Load libraries
library(dotenv)
library(dplyr)
library(rtracklayer)
library(sleuth)

# Load environment variables
print("Both sexes analysis")
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
  full_model = ~ group + sex,
  num_cores = 12,
  read_bootstrap_tpm = TRUE,
  target_mapping = transcripts)

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
# cleanup
remove(hbb_log2_tpm, tpm)

# Re-create sleuth object with new design
so <- sleuth_prep(
  sample_to_covariates = meta, 
  full_model = ~ group + sex + hbb_log2_tpm,
  num_cores = 12,
  read_bootstrap_tpm = TRUE,
  target_mapping = transcripts)

# Calculate percentage isoform used (PIU)
piu <- so$obs_norm_filt %>% 
  left_join(so$target_mapping, by = "target_id") %>%
  group_by(gene_id, sample) %>%
  filter(sum(tpm) > 0) %>% 
  mutate(piu = tpm / sum(tpm)) %>%
  ungroup() %>% 
  select(target_id, sample, tpm, piu)
so$piu <- piu

# Set comparisons
myContrasts <- list(
  c("L.8h", "S.8h"), 
  c("L.24h", "S.24h"),
  c("L.7d", "S.7d"),
  c("H.8h", "S.8h"),
  c("H.24h", "S.24h"),
  c("H.7d", "S.7d")
)

# Differential expression loop
orig_groups <- so$sample_to_covariates$group

for (ctr in myContrasts) {
  num <- ctr[1]   # numerator
  ref <- ctr[2]   # reference
  
  message("Running contrast: ", num, " vs ", ref)
  
  # Reset group to original levels each time
  so$sample_to_covariates$group <- orig_groups
  so$sample_to_covariates$group <- relevel(so$sample_to_covariates$group, ref = ref)
  
  # Give model a unique name
  fit_name  <- paste0("full_", num, "_vs_", ref)
  so <- sleuth_fit(so, ~ group + sex + hbb_log2_tpm, fit_name = fit_name)
  
  # Construct expected beta name
  design_cols <- colnames(models(so)$full$design_matrix)
  beta_name <- paste0("group", num)
  
  if (!beta_name %in% design_cols) {
    stop("Could not find beta '", beta_name,
         "' in design matrix.")
  }
  
  # Wald test for this contrast
  so <- sleuth_wt(so, which_beta = beta_name, which_model = fit_name)
  
  # Extract results
  res <- sleuth_results(so, test = beta_name, test_type = "wt", which_model = fit_name)
  
  # Filter by transcript_type
  res <- subset(res, transcript_type %in% c("protein_coding",
                                            "retained_intron",
                                            "nonsense_mediated_decay"))
  
  # Save
  out_file <- paste0(out, "DEG_tables/", num, "_vs_", ref, "_DEGs.tsv")
  write.table(res, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  message("Saved: ", out_file)
}

saveRDS(object = so, file = "../rObjects/so_group_sex_hbb.rds")
