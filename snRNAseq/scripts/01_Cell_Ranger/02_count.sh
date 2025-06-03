#!/bin/sh
#SBATCH --job-name=n10x_multi       # Name of the job
#SBATCH --mem=50G                   # Amount of memory allocated for the job
#SBATCH --tasks=32                  # Number of tasks (or CPU cores)
#SBATCH --output=logs/%x.%j.stdout  # File for standard output
#SBATCH --error=logs/%x.%j.stderr   # File for standard error output
#SBATCH --time=24:00:00             # Maximum time the job is allowed to run, HH:MM:SS

# source settings and environment variables
source $HOME/.bash_profile
source ../../refs/.env

# print cellranger version
cellranger -V

# change directory to your desired output folder
cd $PROJECT_DIR/counts
echo "Project directory: $PROJECT_DIR"

# print fastq directory
echo "FASTQ directory: $FASTQ_DIR"

# Config CSV path
CONFIG_FILE="${PROJECT_DIR}/refs/config.csv"
echo "Using config: $CONFIG_FILE"

# print sample variable passed from 03_sample_loop.sh script
SAMPLE=$1
echo "Sample: $SAMPLE"

# run cellranger
cellranger multi \
  --id="$SAMPLE" \
  --csv="$CONFIG_FILE" \
  --localcores=$SLURM_NTASKS \
  --localmem=$(($SLURM_MEM_PER_NODE / 1024))

# key:
# --id, sets the name of the output directory for this run (e.g., counts/JDFSeq_0169)
# --csv, points to the config CSV that includes fastq_ids and sample_id for this sample
# --localcores, number of CPU cores to use (matches --cpus-per-task from SLURM)
# --localmem, memory in GB (matches --mem from SLURM)
# [gene-expression] section in config file defines the transcriptome path and options
# [libraries] section lists the two flowcell-specific FASTQ prefixes for the sample
# [samples] section defines the sample_id (should match $SAMPLE)
# Output goes to: $PROJECT_DIR/counts/$SAMPLE/outs/
