#!/bin/sh
#SBATCH --job-name fastqc
#SBATCH --mem 50G
#SBATCH --tasks 30
#SBATCH --mail-user todd.kennedi@mayo.edu
#SBATCH --mail-type END,FAIL
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-short
#SBATCH --time 4:00:00 ## HH:MM:SS
#SBATCH --propagate=NONE

# activate conda environment
source $HOME/.bash_profile
conda activate psilo

# change directory to raw reads
#cd ../../rawReads

# run raw fastqc
#fastqc --threads 30 --outdir ../rawQC *.fastq.gz

# change directory to trimmed reads
cd ../../trimmedReads

# run trimmed fastqc
fastqc --threads 30 --outdir ../trimmedQC *.fastq.gz

