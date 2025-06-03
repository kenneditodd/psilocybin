#!/usr/bin/env Rscript

# Set working directory
setwd(".")

# Load packages
library(dotenv)
library(dplyr)
library(gtools)
library(DESeq2)
library(parallel)
library(Seurat)
library(stringr)
library(tidyr)

# Source thresholds and output paths
source("../../refs/thresholds_and_outs.R")
out <- paste0(out, "pass2/")

# Source custom functions
files <- list.files("../../functions", full.names = TRUE)
invisible(lapply(files, source))

# Load environment variables
load_dot_env(file = "../../refs/.env")
project_dir <- Sys.getenv("PROJECT_DIR")
de_method <- Sys.getenv("DE_METHOD", unset = "DESeq2")

# Read Seurat object
mouse.annotated <- readRDS(paste0(project_dir, "/rObjects/", filtering_method,
                                  "_pass2_annotated_seurat_obj.rds"))

# Define comparisons
comparisons_no_sex <- list(
  c("L.24h", "S.24h"), 
  c("H.24h", "S.24h")
)

comparisons_with_sex <- list(
  c("L.24h.F", "S.24h.F"), 
  c("L.24h.M", "S.24h.M"),
  c("H.24h.F", "S.24h.F"), 
  c("H.24h.M", "S.24h.M")
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 09a_DE.R [mode: no_sex|sex] [index]")

mode <- args[1]
index <- as.integer(args[2])

# Select comparison and metadata column
if (mode == "no_sex") {
  comparison <- comparisons_no_sex[[index]]
  group_col <- "group"
  out_dir <- paste0(out, "DEGs_", de_method, "/both_sexes/DEG_tables")
} else if (mode == "sex") {
  comparison <- comparisons_with_sex[[index]]
  group_col <- "group2"
  out_dir <- paste0(out, "DEGs_", de_method, "/sex_specific/DEG_tables")
} else {
  stop("Invalid mode. Use 'no_sex' or 'sex'")
}

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Run DE
group1 <- comparison[1]
group2 <- comparison[2]
message(sprintf("Comparing %s vs. %s using %s (%s)", group1, group2, group_col, de_method))

DE_within_each_cluster(
  obj = mouse.annotated,
  outDir = out_dir,
  clusterCol = "annotated_clusters",
  groupCol = group_col,
  method = de_method,
  group1 = group1,
  group2 = group2
)
