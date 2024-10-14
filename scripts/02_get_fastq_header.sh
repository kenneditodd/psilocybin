#!/bin/bash
# This script prints each fastq file name along with its header and saves to fastq_headers.tsv.

# Go to fastq directory
cd ../rawReads || exit 1  # Exit if directory change fails

# Set variables
fastq_files="../refs/fastq_file_list.txt"
output_file="../refs/fastq_headers.tsv"

# Check if the file list exists
if [[ ! -f "$fastq_files" ]]; then
  echo "File list not found: $fastq_files"
  exit 1
fi

# Clear the output file if it exists, or create a new one
> "$output_file"

# Print fastq file name + header and save to output file
while IFS= read -r file; do
  if [[ -f "$file" ]]; then
    header=$(zcat "$file" | head -1)
    echo -e "${file}\t${header}" >> "$output_file"  # Append to output file
  else
    echo "File not found: $file" >> "$output_file"  # Log missing files
  fi
done < "$fastq_files"

# Message after completion
echo "Headers saved to $output_file"
