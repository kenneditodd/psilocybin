# read table
mouse.meta <- read.delim2("../../refs/joes_meta.tsv", sep = "\t", header = TRUE)

# rename columns
names <- c("animal_id", "filename", "RLIMS_number", "RLIMS_name", "sacrifice_batch", 
           "group", "treatment", "dose", "timepoint", "sex", "age_at_treatment",
           "weight_at_treatment","RIN","genotype","notes")
colnames(mouse.meta) <- names

# remove unnecessary columns
mouse.meta$RLIMS_name <- NULL
mouse.meta$RLIMS_number <- NULL

# add project_id
mouse.meta$project_id <- "Psi1"

# fix column units
mouse.meta$age_at_treatment <- paste0(mouse.meta$age_at_treatment, " weeks")
mouse.meta$weight_at_treatment <- paste0(mouse.meta$weight_at_treatment, " grams")
mouse.meta$treatment <- tolower(mouse.meta$treatment)

# make the group col more informative
short.treatment <- gsub("psilocybin","psilo",mouse.meta$treatment)
short.treatment <- gsub("saline","sal",short.treatment)
dose.group <- gsub("1 mg / kg", "high", mouse.meta$dose)
dose.group <- gsub("0.25 mg / kg", "low", dose.group)
time <- gsub(" hour", "h", mouse.meta$timepoint)
time <- gsub(" day", "d", time)
group <- paste0(short.treatment, ".", dose.group, ".", time)
group <- gsub(".0.9%","",group)
mouse.meta$group <- group

# sample_id
sample.id <- paste0(short.treatment, ".", dose.group, ".", time,".",
                    mouse.meta$sex,".",mouse.meta$animal_id)
sample.id <- gsub(".0.9%","",sample.id)
mouse.meta$sample_id <- sample.id

# RNA extraction batch (not library prep)
# From Joe's email
# Group 1: A81, A86
# Group 2: A77-A80, A82-A85
# Group 3: A66-A76, A87-A90
# Group 4: A51-A65
# Group 5: A37-A50
# Group 6: A23-A36
# Group 7: A9-A22
# Group 8: A1-A8
mouse.meta$RNA_extraction_batch <- ""
mouse.meta[c(81,86),"RNA_extraction_batch"] <- 1
mouse.meta[c(77:80,82:85),"RNA_extraction_batch"] <- 2
mouse.meta[c(66:76,87:90),"RNA_extraction_batch"] <- 3
mouse.meta[c(51:65),"RNA_extraction_batch"] <- 4
mouse.meta[c(37:50),"RNA_extraction_batch"] <- 5
mouse.meta[c(23:36),"RNA_extraction_batch"] <- 6
mouse.meta[c(9:22),"RNA_extraction_batch"] <- 7
mouse.meta[c(1:8),"RNA_extraction_batch"] <- 8

# rearrange columns
colnames(mouse.meta)[c(15,4:8,1,9:14,3,16,2)]
meta <- mouse.meta[,c(15,4:8,1,9:14,3,16,2)]

# remove time point 28 (wasn't sent for sequencing)
meta <- meta[!meta$timepoint == "28 day",]

# save meta
write.table(x = meta, file = "../../refs/metadata.tsv", sep = "\t",
            quote = FALSE)



