#!/bin/bash
#SBATCH --job-name=tune_umap_pass1_harmony
#SBATCH --ntasks=1
#SBATCH --mem=60GB
#SBATCH --time=1:00:00
#SBATCH --output=logs/%x.%A_%a.stdout
#SBATCH --error=logs/%x.%A_%a.stderr
#SBATCH --array=2-211

# Source shared variables
source ../../refs/.env

export SIMG_PATH="/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif"
export SCRIPT_PATH="/tgen_labs/jfryer/ktodd/psi1/snRNAseq/scripts/02_soupx/tune_umap.R"

# Load Singularity
module load singularity/3.11.5

# Read parameter line corresponding to this array task ID
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../refs/umap_params.csv)
DIMS=$(echo $PARAMS | cut -d',' -f1 | tr -d '\r')
MIN_DIST=$(echo $PARAMS | cut -d',' -f2 | tr -d '\r')
N_NEIGHBORS=$(echo $PARAMS | cut -d',' -f3 | tr -d '\r')

# Run script inside container
singularity exec -B /tgen_labs/jfryer:/tgen_labs/jfryer \
  "$SIMG_PATH" /usr/local/bin/Rscript "$SCRIPT_PATH" "$DIMS" "$MIN_DIST" "$N_NEIGHBORS"
