#!/bin/sh
#SBATCH --partition=gpu-a100        # Explicitly set GPU partition
#SBATCH --nodes=1                   # Requests 1 compute node
#SBATCH --ntasks=1                  # Requests 1 task/process
#SBATCH --gres=gpu:1                # Requests 1 GPU
#SBATCH --time=8:00:00              # Maximum time the job is allowed to run, HH:MM:SS
#SBATCH --job-name=cellbender       # Name of the job
#SBATCH --mem=50G                   # Memory allocation
#SBATCH --output=logs/%x.%j.stdout  # Stdout log file
#SBATCH --error=logs/%x.%j.stderr   # Stderr log file

# source settings
source $HOME/.bash_profile
source ../../refs/.env

# activate conda environment
conda activate cellbender

# check gpu
nvidia-smi

# cellbender version
cellbender --version

# print sample passed from sample_loop.sh script
SAMPLE=$1
echo "sample: $SAMPLE"

# in order for ckpt.tar.gz in cellbender to save correctly, go to where data is located
cd ${PROJECT_DIR}/counts/${SAMPLE}/outs

# run cellbender remove-background
cellbender remove-background \
      --input raw_feature_bc_matrix.h5 \
      --output ${SAMPLE}_cellbender.h5 \
      --cuda

# Key
# --input       :un-filtered data file, .h5 or matrix directory from 10x
# --output      :output file location containing .h5 extension
# --cuda        :for use with GPU
