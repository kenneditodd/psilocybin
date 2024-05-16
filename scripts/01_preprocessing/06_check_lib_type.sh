#!/bin/sh
#SBATCH --job-name check_lib_type
#SBATCH --mem 10G
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-med
#SBATCH --tasks 15
#SBATCH --time 02:00:00 ## HH:MM:SS

# salmon was installed in it's own environment to work properly
# activate environment
source $HOME/.bash_profile
conda activate salmon

# salmon version
salmon -v

# go to raw reads dir
cd ../../rawReads

# validate mappings
# note this same sample has other lanes but this should be sufficient to check
salmon quant --libType A \
             --index /research/labs/neurology/fryer/projects/references/mouse/salmonIndexGRCm39 \
             --mates1 Psi1_A10_High_Female_L3_R1.fastq.gz \
             --mates2 Psi1_A10_High_Female_L3_R2.fastq.gz \
             --output ../refs/transcript_quant \
             --threads 15 \
             --validateMappings

# Key
# --libType A is for autodetect library type
# --index Salmon index
# --mates1
# --mates2

# Results
# Automatically detected most likely library type as ISR
# ISR = inward stranded reverse
# This is the -s2 argument for featureCounts (reversely stranded)

# Notes
# You can also check your library type in IGV


