# read files
df <- read.delim2("../../refs/metadata.tsv")

# SRA biosample
biosample <- data.frame(sample_name = df$sample_id,
                        organism = "Mus musculus",
                        strain = "C57BL/6",
                        age = "8-12 weeks at injection",
                        collection_date = "01-Oct-2023",
                        geo_loc_name = "USA",
                        tissue = "prefrontal cortex")
# Combine and reorder
biosample <- cbind(biosample, df[, c(
  "group",
  "group2",
  "treatment",
  "dose",
  "timepoint",
  "sex",
  "animal_id",
  "RNA_extraction_batch",
  "sacrifice_batch",
  "age_at_treatment",
  "weight_at_treatment",
  "project_id",
  "filename"
)])
write.table(biosample, "../../refs/SRA_BioSample.tsv", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# SRA metadata
sra.meta <- data.frame(sample_name = df$sample_id,
                       library_ID = df$animal_id,
                       title = paste0("Bulk RNA-Seq of mouse prefrontal cortex"),
                       library_strategy = "RNA-Seq",
                       library_source = "TRANSCRIPTOMIC",
                       library_selection = "PolyA",
                       library_layout = "paired",
                       platform = "ILLUMINA",
                       instrument_model = "Illumina NovaSeq 6000",
                       design_description = "Illumina Stranded mRNA Prep",
                       filetype = "fastq",
                       filename = paste0(df$filename, "_L", df$lane, "_R1.fastq.gz"),
                       filename2 = paste0(df$filename, "_L", df$lane, "_R2.fastq.gz"))
# Combine and reorder
sra.meta <- cbind(sra.meta, df[, c(
  "group",
  "group2",
  "treatment",
  "dose",
  "timepoint",
  "sex",
  "animal_id",
  "RNA_extraction_batch",
  "sacrifice_batch",
  "age_at_treatment",
  "instrument",
  "run",
  "flow_cell",
  "lane"
)])
write.table(sra.meta, "../../refs/SRA_Metadata.tsv", 
            quote = FALSE, row.names = FALSE, sep = "\t")
