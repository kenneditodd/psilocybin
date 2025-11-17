#!/bin/bash
while read -r sample; do
    sbatch 02_cellranger_multi.sh "$sample"
done < ../../refs/sample_list.tsv
