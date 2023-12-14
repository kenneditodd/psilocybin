# read table
mouse.meta <- read.delim2("../../refs/mouse_meta.tsv", sep = "\t", header = TRUE)

# rename columns
names <- c("animal_id", "filename", "RLIMS_number", "RLIMS_name", "sacrifice_batch", 
           "group", "treatment", "dose", "timepoint", "sex", "age_at_treatment",
           "weight_at_treatment","RIN","genotype","notes")
colnames(mouse.meta) <- names

# add project_id
mouse.meta$project_id <- "Psil1"

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

# rearrange columns
colnames(mouse.meta)[c(17,6,7:10,1,11:16,2:5)]
meta <- mouse.meta[,c(17,6,7:10,1,11:16,2:5)]

# save meta
write.table(x = mouse.meta, file = "../../refs/metadata.tsv", sep = "\t",
            quote = FALSE)



