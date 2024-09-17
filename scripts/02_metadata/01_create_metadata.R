# Load libraries
library(dplyr)

# Read table
mouse.meta <- read.delim2("../../refs/joes_meta.tsv", sep = "\t", header = TRUE)

# Rename columns explicitly
mouse.meta <- mouse.meta %>%
  rename(
    animal_id = Animal_Num,
    filename = Unique_ID,
    RLIMS_number = RLIMS_Num,
    RLIMS_name = RLIMS_Name,
    sacrifice_batch = Sacrifice_Group_Num,
    group = Group_Num,
    treatment = Treatment,
    dose = Dose,
    timepoint = Timepoint,
    sex = Sex,
    age_at_treatment = Age_Weeks,
    weight_at_treatment = Weight_grams_day_of_treatment,
    RIN = RIN,
    genotype = Genotype,
    notes = Notes
  )

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
group <- paste0(gsub(" ", "", mouse.meta$dose),
                ".",
                gsub(" hour", "h", mouse.meta$timepoint))
group <- gsub(" day", 'd', group)
group <- gsub("0.9%", 'S', group)
group <- gsub("0.25mg/kg", 'L', group)
group <- gsub("1mg/kg", 'H', group)
mouse.meta$group <- group
group2 <- paste0(group, ".", mouse.meta$sex)
mouse.meta$group2 <- group2

# sample_id
sample.id <- paste0(group2, ".", mouse.meta$animal_id)
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
meta <- mouse.meta[, c(
  "sample_id",
  "group",
  "group2",
  "treatment",
  "dose",
  "timepoint",
  "sex",
  "animal_id",
  "RNA_extraction_batch",
  "sacrifice_batch",
  "age_at_treatment",
  "weight_at_treatment",
  "RIN",
  "genotype",
  "notes",
  "project_id",
  "filename"
)]

# remove 28 day time point (wasn't sent for sequencing)
meta <- meta[!meta$timepoint == "28 day",]

# save meta
write.table(x = meta, 
            file = "../../refs/metadata.tsv", 
            sep = "\t",
            quote = FALSE)
