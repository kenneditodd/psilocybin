#!/usr/bin/env python3

import os
import json

# Set file paths
sample_list_path = '../refs/sample_file_list.txt'
output_file_path = '../refs/config.json'

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

# Open the sample list file and process each line
with open(sample_list_path, 'r') as infile:
  
    # Read each line, strip whitespace, and remove the _R1 suffix
    sample_base_names = []  # Create an empty list

    for line in infile:
        line = line.strip()
        line = line.split(".")[0]
        line = line.replace("_R1", "")
        sample_base_names.append(line)
    
    # Loop through each sample name in sample_base_names
    male_samples = []
    female_samples = []
    for sample in sample_base_names:
        if "_Male" in sample:
            male_samples.append(sample)
        if "_Female" in sample:
            female_samples.append(sample)

# Create the config structure
config = {
    "DIRECTORIES": directories,
    "SAMPLE_INFORMATION": {
        "allSamples": sample_base_names,
        "femaleSamples": female_samples,
        "maleSamples": male_samples
    },
    "CLUSTER_INFORMATION": cluster_info,
    "SAMPLES": {}
}

# Populate sample information for read1 and read2
with open(sample_list_path, 'r') as infile:
    for line in infile:
        sample = line.strip()
        base_name = sample.split(".")[0]
        read1 = base_name
        read2 = base_name.replace("_R1", "_R2")
        base_name = base_name.replace("_R1", "")

        config["SAMPLES"][base_name] = {
            "read1": read1,
            "read2": read2
        }

# Write config to JSON file
with open(output_file_path, 'w') as outfile:
    json.dump(config, outfile, indent=4)

print(f"Config file saved to {output_file_path}")
