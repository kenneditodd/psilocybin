#!/bin/bash

# source environment variables
source ../../refs/.env

# extract webs summaries
while read -r sample; do
	cd $PROJECT_DIR/counts/$sample/outs/per_sample_outs/$sample
	cp web_summary.html ../../../../web_summaries/"$sample"_web_summary.html
done < $PROJECT_DIR/refs/sample_list.tsv

