#!/bin/bash

# Load Singularity
module load singularity

# Load environment variables
source ../../refs/.env

# Parameters
IMAGE="/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif"
R_SCRIPT="${PROJECT_DIR}/scripts/02_soupx/14_recluster_neurons_DE.R"
DE_METHOD="MAST"

# Submit no_sex comparisons (2 total)
for i in {1..2}; do
  sbatch --job-name=DE_${DE_METHOD}_both_sexes_${i} \
         --output=logs/DE_${DE_METHOD}_both_sexes_${i}.out \
         --error=logs/DE_${DE_METHOD}_both_sexes_${i}.err \
         --time=24:00:00 \
         --mem=300G \
         --export=ALL,MODE=both_sexes,INDEX=$i,DE_METHOD=$DE_METHOD,R_SCRIPT=$R_SCRIPT,IMAGE=$IMAGE \
         --wrap="singularity exec --bind /tgen_labs \$IMAGE Rscript \$R_SCRIPT \$MODE \$INDEX"
done

# Submit sex comparisons (4 total)
for i in {1..4}; do
  sbatch --job-name=DE_${DE_METHOD}_sex_specific_${i} \
         --output=logs/DE_${DE_METHOD}_sex_specific_${i}.out \
         --error=logs/DE_${DE_METHOD}_sex_specific_${i}.err \
         --time=24:00:00 \
         --mem=300G \
         --export=ALL,MODE=sex_specific,INDEX=$i,DE_METHOD=$DE_METHOD,R_SCRIPT=$R_SCRIPT,IMAGE=$IMAGE \
         --wrap="singularity exec --bind /tgen_labs \$IMAGE Rscript \$R_SCRIPT \$MODE \$INDEX"
done
