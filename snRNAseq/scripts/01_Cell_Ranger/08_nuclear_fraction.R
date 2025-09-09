# working directory and libraries
setwd(".")
library(dotenv)
library(DropletQC) # package installed from https://github.com/powellgenomicslab/DropletQC

# load the environment variables
load_dot_env(file = "../../refs/.env")

# access the environment variables
project_dir <- Sys.getenv("PROJECT_DIR")
annotation_dir <- Sys.getenv("ANNOTATION_REFS")
sample <- Sys.getenv("SAMPLE")

# annotation file location
gtf_file <- paste0(annotation_dir, "/refdata-gex-GRCm39-2024-A/genes/genes.gtf")

# get nuclear fraction
nf <- nuclear_fraction_annotation(
  annotation_path = gtf_file,
  bam = paste0(project_dir, "/counts/", sample, "/outs/per_sample_outs/", sample, "/count/sample_alignments.bam"),
  barcodes = paste0(project_dir, "/counts/", sample, "/outs/per_sample_outs/", sample, "/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"),
  verbose = FALSE
  )

# reformat
nf$barcode <- rownames(nf)

# write output
write.table(x = nf,
            file = paste0(project_dir, "/counts/", sample, "/", sample, "_nuclear_fraction.tsv"),
            quote = FALSE,
            row.names = FALSE)
