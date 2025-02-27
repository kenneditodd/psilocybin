#!/bin/sh
#SBATCH --job-name=run_snakemake              # Job name
#SBATCH --mem=5G                              # Memory allocation
#SBATCH --mail-type=END,FAIL                  # Notify on job end or failure
#SBATCH --output=logs/%x.%j.stdout            # Standard output log
#SBATCH --error=logs/%x.%j.stderr             # Standard error log
#SBATCH --time=24:00:00                       # Time limit (HH:MM:SS)

# Activate conda environment
source $HOME/.bash_profile
conda activate psilo

# Change directory to where the Snakefile is located
cd ../

# Before running Snakemake, you should do a dry run with the command below
#snakemake --dry-run --printshellcmds

# Run Snakemake
snakemake --snakefile Snakefile --jobs 20 --rerun-incomplete --latency-wait 60 --cluster "sbatch --mem=60G --output=scripts/logs/snakemake_job_logs/%x.%N.%j.stdout --error=scripts/logs/snakemake_job_logs/%x.%N.%j.stderr --tasks=20 --time=05:00:00"

# Key for the Snakemake command:
# --snakefile Snakefile: Specifies the Snakefile to use.
# --jobs 20            : Runs up to 20 parallel jobs.
# --rerun-incomplete   : Reruns jobs that were incomplete.
# --latency-wait 60    : Waits 60 seconds for files to appear before failing.
# --cluster            : Submits jobs using the sbatch command with the following arguments:

# Key for the Snakemake SLURM job:
# --mem=60G             : Allocates 60 GB of memory for each job.
# --output              : Specifies the path for standard output logs.
# --error               : Specifies the path for standard error logs.
# --tasks=20            : Runs 20 tasks for each job.
# --time=05:00:00       : Sets a time limit of 5 hours for each job.
