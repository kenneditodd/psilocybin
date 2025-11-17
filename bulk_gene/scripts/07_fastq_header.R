# read files
meta <- read.delim2("refs/metadata.tsv")
seq.info <- read.delim2("refs/fastq_headers.tsv", header = FALSE)

# reformat seq.info
colnames(seq.info) <- c("filename","fastq_header")
seq.info$fastq_header <- gsub("@","",seq.info$fastq_header)
seq.info <- tidyr::separate_wider_delim(seq.info, 
                                        cols = fastq_header, 
                                        delim = " ", 
                                        names = c("header1","header2"))
seq.info <- tidyr::separate_wider_delim(seq.info, 
                                        cols = header1, 
                                        delim = ":", 
                                        names = c("instrument","run","flow_cell","lane",
                                                  "xpos","ypos","tile"))

# remove unnecessary info
seq.info$header2 <- NULL
seq.info$xpos <- NULL
seq.info$ypos <- NULL
seq.info$tile <- NULL

# remove R2 files
keep <- is.na(stringr::str_match(seq.info$filename, ".+_R2"))
seq.info <- seq.info[keep,]

# check - all on the same instrument, run, and flow cell
table(seq.info$instrument)
table(seq.info$run)
table(seq.info$flow_cell)

# get the base of the filename
seq.info$filename <- 
  stringr::str_match(seq.info$filename, 
                     "(Psi1_A[0-9]+_[HighLowNS]+_[FemMale]+)_L[34]_R1.fastq.gz")[,2]

# merge
df <- dplyr::left_join(meta, seq.info, by = "filename")

# check lane distribution with group
table(paste0("Lane", df$lane, "-Group", df$group))
table(paste0("Lane", df$lane, "-Group", df$group2))

# save joined meta
write.table(x = df, file = "refs/metadata.tsv", sep = "\t",
            quote = FALSE)
