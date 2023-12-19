#!/bin/bash
# This script will print to standard out.
# Redirect output to save.

# set variables
files=/research/labs/neurology/fryer/m214960/psilocybin/refs/fastq_file_list.txt

# go to fastq dir
cd /research/labs/neurology/fryer/projects/psilocybin/psil1

# print fastq file name + header
cat $files | while read file
do
  header=$(zcat $file | head -1)
  echo -n $file && echo -ne '\t' && echo $header
done
