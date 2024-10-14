#!/usr/bin/env python3

import os
import json

# Set file paths
sample_list_path = '../../refs/sample_file_list.txt'
output_file_path = '../../refs/config.json'

# Initialize directories and sample data
directories = {
    "rawReads": "rawReads/",
    "rawQC": "rawQC/",
    "trimmedReads": "trimmedReads/",
    "trimmedQC": "trimmedQC/",
    "starAligned": "starAligned/",
    "featureCounts": "featureCounts/",
    "genomeDir": "refs/starGenomeDir/"
}

cluster_info = {
    "threads": "20"
}

# Load sample names and strip _R1 suffix
with open(sample_list_path, 'r') as infile:
    all_samples = [line.strip().replace("_R1.fastq.gz", "") for line in infile]

# Create the config structure
config = {
    "DIRECTORIES": directories,
    "SAMPLE_INFORMATION": {
        "allSamples": all_samples
    },
    "CLUSTER_INFORMATION": cluster_info,
    "SAMPLES": {}
}

# Populate sample information for read1 and read2
with open(sample_list_path, 'r') as infile:
    for line in infile:
        sample = line.strip()
        base_name = sample.replace("_R1.fastq.gz", "")
        read1 = sample.replace(".fastq.gz", "")
        read2 = sample.replace("_R1.fastq.gz", "_R2")
        config["SAMPLES"][base_name] = {
            "read1": read1,
            "read2": read2
        }

# Write config to JSON file
with open(output_file_path, 'w') as outfile:
    json.dump(config, outfile, indent=4)

print(f"Config file saved to {output_file_path}")
