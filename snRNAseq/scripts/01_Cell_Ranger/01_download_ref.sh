#!/bin/bash

# From cellranger
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"

# Extract files
tar -xvzf refdata-gex-GRCm39-2024-A.tar.gz
