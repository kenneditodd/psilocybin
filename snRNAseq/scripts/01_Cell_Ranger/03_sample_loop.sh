#!/bin/bash
while read -r sample; do
    sbatch 02_count.sh "$sample"
done < ../../refs/sample_list.tsv
