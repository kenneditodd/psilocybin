#!/bin/bash
#SBATCH --job-name=sleuth
#SBATCH --cpus-per-task=12
#SBATCH --mem=60GB
#SBATCH --time=2:00:00
#SBATCH --output=logs/%x.%A_%a.stdout
#SBATCH --error=logs/%x.%A_%a.stderr

# Source shared variables
source ../refs/.env

export SIMG_PATH="/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif"
export SCRIPT_PATH="${PROJECT_DIR}/scripts/05a_sleuth.R"

# Load Singularity
module load singularity/3.11.5

# Run script inside container
singularity exec -B /tgen_labs/jfryer:/tgen_labs/jfryer \
  "$SIMG_PATH" /usr/local/bin/Rscript "$SCRIPT_PATH"
