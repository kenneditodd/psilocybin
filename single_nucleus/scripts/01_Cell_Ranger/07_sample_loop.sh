#!/bin/bash
while read -r sample; do
    sbatch 06_cellbender.sh "$sample"
done < ../../refs/sample_list.tsv
