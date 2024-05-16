#!/bin/bash
# This script will print to standard out.
# Redirect output to save.

# go to fastq dir
cd ../../rawReads

# set variables
files=../refs/fastq_file_list.txt

# print fastq file name + header
cat $files | while read file
do
  header=$(zcat $file | head -1)
  echo -n $file && echo -ne '\t' && echo $header
done
