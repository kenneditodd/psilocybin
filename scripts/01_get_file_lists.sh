#!/bin/bash

# Load environment variables
source ../refs/.env

# Get fastq file list
# There are 2 fastq files per sample
cd $RAW_READS_DIR
out="${PROJECT_DIR}/refs/fastq_file_list.txt"
ls -1 | grep .fastq.gz > $out

# Get sample file list
out="${PROJECT_DIR}/refs/sample_file_list.txt"
ls -1 | grep _R1 > $out
