# load libraries
library(stringr)
library(dotenv)
library(dplyr)

# load the environment variables
load_dot_env(file = "../../refs/.env")

# access the environment variables
project_dir <- Sys.getenv("PROJECT_DIR")

# file locations
samples <- readLines(paste0(project_dir, "/refs/sample_list.tsv"))
locations <- paste0(project_dir, "/counts/", samples, "/outs/metrics_summary.csv")

# initialize df and loop through files
df <- data.frame()
for (i in 1:length(locations)) {
  if (i == 1) {
    df <- read.csv(locations[i])
  } else {
    row <- read.csv(locations[i])[1,]
    df <- rbind(df,row)
  }
}

# rename rows
rownames(df) <- samples

# rename columns
df <- df %>%
  rename(
    estimated_cells = Estimated.Number.of.Cells,
    mean_reads_per_cell = Mean.Reads.per.Cell,
    median_genes_per_cell = Median.Genes.per.Cell,
    number_reads = Number.of.Reads,
    valid_barcodes = Valid.Barcodes,
    sequencing_saturation = Sequencing.Saturation,
    Q30_bases_barcode = Q30.Bases.in.Barcode,
    Q30_bases_read = Q30.Bases.in.RNA.Read,
    Q30_bases_UMI = Q30.Bases.in.UMI,
    reads_mapped_genome = Reads.Mapped.to.Genome,
    confident_reads_mapped_genome = Reads.Mapped.Confidently.to.Genome,
    confident_intergenic_reads_mapped = Reads.Mapped.Confidently.to.Intergenic.Regions,
    confident_intronic_reads_mapped = Reads.Mapped.Confidently.to.Intronic.Regions,
    confident_exonic_reads_mapped = Reads.Mapped.Confidently.to.Exonic.Regions,
    confident_reads_mapped_transcriptome = Reads.Mapped.Confidently.to.Transcriptome,
    reads_mapped_antisense = Reads.Mapped.Antisense.to.Gene,
    fraction_reads = Fraction.Reads.in.Cells,
    total_genes = Total.Genes.Detected,
    median_UMI_per_cell = Median.UMI.Counts.per.Cell
  )

write.table(df, 
            paste0(project_dir, "/counts/web_summaries/overall_metrics.tsv"),
            sep = "\t",
            quote = FALSE)
