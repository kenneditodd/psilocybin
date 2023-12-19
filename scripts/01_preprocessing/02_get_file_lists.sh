#!/bin/bash

# Get fastq file list
# There are 2 fastq files per sample
cd ../../rawReads
out=../refs/fastq_file_list.txt
ls -1 | grep .fastq.gz > $out

# Get sample file list
out=../refs/sample_file_list.txt
ls -1 | grep _R1_ > $out
