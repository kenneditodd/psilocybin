# setup
setwd(".")
library('stringr')

# read fastq names
file <- read.delim2("refs/original_fastq_names.txt",
                    header = FALSE)
colnames(file) <- "original_names"

# create new file names
new <- 
  str_match(file$original,
            "(Psi1_A[0-9]+_[HighLowNS]+_[FmMale]+)\\.[A-Z0-9]+_(L[0-9]_R[12])_[A-Z]+-[A-Z]+(.fastq.gz)")
new <- as.data.frame(new)
file$new_name <- paste0(new$V2, "_", new$V3, new$V4)

# check all names are unique
table(duplicated(file$new_name))

# create bash script
header <- c("#!/bin/bash","cd /research/labs/neurology/fryer/projects/psilocybin/psi1")
path <- "/research/labs/neurology/fryer/m214960/psilocybin/rawReads/"
bsh <- paste0("cp ", file$original_names, " ", path, file$new_name)
script <- as.data.frame(c(header,bsh))

# save script
write.table(x = script,
            file = "scripts/01_preprocessing/03_rename_fastq_files.sh",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

