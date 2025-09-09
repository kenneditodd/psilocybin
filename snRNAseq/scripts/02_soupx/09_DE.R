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

# source thresholds and output paths
source("../../refs/thresholds_and_outs.R")
out <- paste0(out, "pass3_downsampled/")

# Source custom functions
files <- list.files("../../functions", full.names = TRUE)
invisible(lapply(files, source))

# Load environment variables
load_dot_env(file = "../../refs/.env")
project_dir <- Sys.getenv("PROJECT_DIR")
de_method <- Sys.getenv("DE_METHOD")

# Read Seurat object
mouse.annotated <- readRDS(
  paste0(project_dir, "/rObjects/", filtering_method, "_pass3_downsampled_annotated_seurat_obj.rds"))

# Define comparisons
comparisons_both_sexes <- list(
  c("L.24h", "S.24h"), 
  c("H.24h", "S.24h")
)
comparisons_sex_specific <- list(
  c("L.24h.F", "S.24h.F"), 
  c("L.24h.M", "S.24h.M"),
  c("H.24h.F", "S.24h.F"), 
  c("H.24h.M", "S.24h.M")
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript DE.R [mode: both_sexes|sex_specific] [index]")

# Get arguments for DE function
mode <- args[1]
index <- as.integer(args[2])
pct <- 0.6

# Set and create output directory
out_dir <- paste0(out, "DEGs_", de_method, "_pct_", pct, "/DEG_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Select comparison and metadata column
if (mode == "both_sexes") {
  comparison <- comparisons_both_sexes[[index]]
  group_col <- "group"
  latent_vars <- "sex"
} else if (mode == "sex_specific") {
  comparison <- comparisons_sex_specific[[index]]
  group_col <- "group2"
  latent_vars <- NULL
} else {
  stop("Invalid mode. Use 'both_sexes' or 'sex_specific'")
}

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
  group2 = group2,
  pct = pct,
  latentVars = latent_vars
)
