#!/bin/bash

# Get fastq file list
# There are ? fastq files per sample
# Each sample has R1 and R2 on all ? lanes
cd /research/labs/neurology/fryer/projects/psilocybin/psil1
out=/research/labs/neurology/fryer/m214960/psilocybin/refs/fastq_file_list.txt
ls -1 | grep .fastq.gz > $out

# Get sample file list
cd /research/labs/neurology/fryer/projects/psilocybin/psil1
out=/research/labs/neurology/fryer/m214960/psilocybin/refs/sample_file_list.txt
ls -1 | grep _R1_ > $out
