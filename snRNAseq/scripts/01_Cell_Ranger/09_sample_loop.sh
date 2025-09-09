#!/bin/bash

# Load Singularity
module load singularity

# Load environment variables
source ../../refs/.env

# Parameters
SIMG_IMAGE="/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif"
R_SCRIPT="${PROJECT_DIR}/scripts/01_Cell_Ranger/08_nuclear_fraction.R"

# Submit no_sex comparisons (2 total)
while read -r SAMPLE; do
  sbatch --job-name=${SAMPLE}_nuclear_fraction \
         --output=logs/%x.%j.stdout\
         --error=logs/%x.%j.stderr \
         --time=24:00:00 \
         --mem=100G \
         --cpus-per-task=8 \
         --export=ALL,SAMPLE=$SAMPLE,R_SCRIPT=$R_SCRIPT,SIMG_IMAGE=$SIMG_IMAGE \
         --wrap="singularity exec --bind /tgen_labs \$SIMG_IMAGE Rscript \$R_SCRIPT"
done < ../../refs/sample_list.tsv
