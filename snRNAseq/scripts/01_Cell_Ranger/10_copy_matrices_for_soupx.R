# load libraries
library(dotenv)
library(fs)

# load the environment variables
load_dot_env(file = "../../refs/.env")

# access the environment variables
project_dir <- Sys.getenv("PROJECT_DIR")

# get raw matrix directory locations
samples <- readLines(paste0(project_dir, "/refs/sample_list.tsv"))
raw_locations <- paste0(project_dir, "/counts/", samples, "/outs/multi/count/raw_feature_bc_matrix")

# get filtered matrix directory locations
samples <- readLines(paste0(project_dir, "/refs/sample_list.tsv"))
filtered_locations <- paste0(project_dir, "/counts/", samples, "/outs/per_sample_outs/", 
                             samples, "/count/sample_filtered_feature_bc_matrix")

for (i in seq_along(samples)) {
  sample <- samples[i]
  raw <- raw_locations[i]
  filt <- filtered_locations[i]
  
  # Make new soupx subdir
  soupx_dir <- path(project_dir, "soupx", sample)
  dir_create(soupx_dir)
  
  # Copy raw matrix info into new subdir
  dir_copy(raw, path(soupx_dir, "raw_feature_bc_matrix"))
  
  # Copy filtered matrix info into new subdir
  dir_copy(filt, path(soupx_dir, "filtered_feature_bc_matrix"))
  
  cat("✅ Copied for sample:", sample, "\n")
}
