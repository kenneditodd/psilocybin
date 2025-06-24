#!/bin/sh
#SBATCH --job-name=cellranger_multi_part1 # Name of the job
#SBATCH --mem=8G                          # Amount of memory allocated for the job
#SBATCH --cpus-per-task=1                 # Number of tasks (or CPU cores)
#SBATCH --output=logs/%x.%j.stdout        # File for standard output
#SBATCH --error=logs/%x.%j.stderr         # File for standard error output
#SBATCH --time=5-00:00:00                 # Maximum time the job is allowed to run, DD-HH:MM:SS

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

# print sample variable passed from 02b_sample_loop.sh script
SAMPLE=$1
echo "Sample: $SAMPLE"

# Config CSV path
CONFIG_FILE="${PROJECT_DIR}/refs/configs/${SAMPLE}.csv"
echo "Using config: $CONFIG_FILE"

# run cellranger
cellranger multi \
  --id="$SAMPLE" \
  --csv="$CONFIG_FILE" \
  --jobmode=../refs/slurm.template \
  --disable-ui \
  --maxjobs=20 \
  --mempercore=35000

# key:
# --id, sets the name of the output directory for this run (e.g., counts/JDFSeq_0169)
# --csv, points to the config CSV that includes fastq_ids and sample_id for this sample
# --jobmode, points to the SLURM template used to manage sub-job execution
# --disable-ui, turns off the graphical web-based UI (CLI only mode)
# --maxjobs, limits the number of SLURM sub-jobs Cell Ranger can run in parallel (20 here)
# --mempercore, maximum memory (in MB) allowed per SLURM task/core (35 GB in this case)
# [gene-expression] section in config file defines the transcriptome path and options
# [libraries] section lists the two flowcell-specific FASTQ prefixes for the sample
# [samples] section defines the sample_id (should match $SAMPLE)
# Output goes to: $PROJECT_DIR/counts/$SAMPLE/outs/
