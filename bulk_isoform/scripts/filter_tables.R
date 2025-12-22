
# List files
out <- "../results/sex_specific/"
files <- list.files(paste0(out, "DEGs/DEG_tables/"))


# Loop through files, filter, resave
for (i in seq_along(files)) {
  my_table <- read.delim2(paste0(out, "DEGs/DEG_tables/", files[i]))
  my_table <- subset(my_table, transcript_type %in% c("protein_coding",
                                                      "retained_intron",
                                                      "nonsense_mediated_decay"))
  write.table(my_table,
              paste0(out, "DEGs_filtered_transcript_type/DEG_tables/", files[i]), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
}
