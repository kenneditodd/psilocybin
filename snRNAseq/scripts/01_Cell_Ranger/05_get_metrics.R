# load libraries
library(stringr)
library(dotenv)
library(dplyr)
library(readr)
library(tidyr)

# load the environment variables
load_dot_env(file = "../../refs/.env")

# access the environment variables
project_dir <- Sys.getenv("PROJECT_DIR")

# file locations
samples <- readLines(paste0(project_dir, "/refs/sample_list.tsv"))
locations <- paste0(project_dir, "/counts/", samples, "/outs/per_sample_outs/", samples, "/metrics_summary.csv")

# initialize df and loop through files
df <- data.frame()
for (i in 1:length(locations)) {
  if (i == 1) {
    df <- readr::read_csv(locations[i], show_col_types = FALSE)
    df <- df %>%
      filter(`Grouped By` == "Physical library ID") %>%
      select(`Group Name`, `Metric Name`, `Metric Value`) %>%
      pivot_wider(names_from = `Metric Name`, values_from = `Metric Value`)
    df$tgen_subject_name <- samples[i]
  } else {
    row <- readr::read_csv(locations[i], show_col_types = FALSE)
    row <- row %>%
      filter(`Grouped By` == "Physical library ID") %>%
      select(`Group Name`, `Metric Name`, `Metric Value`) %>%
      pivot_wider(names_from = `Metric Name`, values_from = `Metric Value`)
    row$tgen_subject_name <- samples[i]
    row <- row[1,]
    df <- rbind(df,row)
  }
}
remove(row)

# read metadata
meta <- readRDS("../../rObjects/meta.rds")

# add sample_id by joining with tgen_subject_name
df <- left_join(x = df, y = meta[,c(1,16)], by = "tgen_subject_name")
remove(meta)

# rename rows
df <- as.data.frame(df)
rownames(df) <- df$sample_id

# rename/reorder columns
colnames(df) <- gsub(" ", "_", tolower(colnames(df)))
df <- df[,c(18,17,2,11:16,3:10)]

# sort by ascending # of cells
df <- df[order(df$cells),]

write.table(df, 
            paste0(project_dir, "/counts/web_summaries/overall_metrics.tsv"),
            sep = "\t",
            quote = FALSE)
