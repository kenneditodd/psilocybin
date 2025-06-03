#!/bin/bash

# Load Singularity
module load singularity

# Load environment variables
source ../../refs/.env

# Parameters
SIMG_IMAGE="/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif"
R_SCRIPT="${PROJECT_DIR}/scripts/02_R_analysis/09a_DE.R"
DE_METHOD="DESeq2"  # <-- Set method here: "DESeq2" or "MAST"

# Submit no_sex comparisons (2 total)
for i in {1..2}; do
  sbatch --job-name=DE_${DE_METHOD}_no_sex_${i} \
         --output=logs/DE_${DE_METHOD}_no_sex_${i}.out \
         --error=logs/DE_${DE_METHOD}_no_sex_${i}.err \
         --time=24:00:00 \
         --mem=500G \
         --export=ALL,MODE=no_sex,INDEX=$i,DE_METHOD=$DE_METHOD,R_SCRIPT=$R_SCRIPT,SIMG_IMAGE=$SIMG_IMAGE \
         --wrap="singularity exec --bind /tgen_labs \$SIMG_IMAGE Rscript \$R_SCRIPT \$MODE \$INDEX"
done

# Submit sex comparisons (4 total)
for i in {1..4}; do
  sbatch --job-name=DE_${DE_METHOD}_sex_${i} \
         --output=logs/DE_${DE_METHOD}_sex_${i}.out \
         --error=logs/DE_${DE_METHOD}_sex_${i}.err \
         --time=24:00:00 \
         --mem=500G \
         --export=ALL,MODE=sex,INDEX=$i,DE_METHOD=$DE_METHOD,R_SCRIPT=$R_SCRIPT,SIMG_IMAGE=$SIMG_IMAGE \
         --wrap="singularity exec --bind /tgen_labs \$SIMG_IMAGE Rscript \$R_SCRIPT \$MODE \$INDEX"
done
