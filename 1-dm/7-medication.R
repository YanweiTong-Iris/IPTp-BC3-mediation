################################################################
# IPTp and child growth
# Script for maternal medication
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

all_visits <- readRDS((paste0(data_path, "outcome2_zscores_data.RDS")))

# Categorizing medicine codes in the dataset into 8 categories of medicine 
antibacterial <-c(8, 12, 17, 99, 253, 19, 32, 276, 9, 98, 13, 11, 16, 43, 52, 53, 184, 199)
betalactam <-c(8, 12, 17, 253, 19, 32, 276, 9, 13, 11, 16)
fluoroquinolone <- c(43, 184)
sulfonamide <- c(52, 53)
antibacterial_other <- c(99, 98, 199)
antimalarial <- c(50, 183, 124, 125, 137, 255)
antiparasitic <- c(96)
antiviral <- c(1)

# Find the minimum and maximum values for codes
nmed_combined <- unlist(all_visits[, paste0("nmed", 1:7)])
min_value <- min(nmed_combined, na.rm = TRUE) # 1
max_value <- max(nmed_combined, na.rm = TRUE) # 280

anymed <- c(min_value:max_value)

# Total counts of medicine taken for each category, within each row of the dataset 
get_category_counts <- function(row, category) {
    sum(row %in% category)
  }
  
  # Apply the function to each row in the dataset
  antibacterial_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = antibacterial)
  betalactam_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = betalactam)
  fluoroquinolone_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = fluoroquinolone)
  sulfonamide_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = sulfonamide)
  antibacterial_other_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = antibacterial_other)
  antimalarial_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = antimalarial)
  antiparasitic_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = antiparasitic)
  antiviral_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = antiviral)
  anymed_counts <- apply(all_visits[, paste0("nmed", 1:7)], 1, get_category_counts, category = anymed)


# Creating 8 new binary columns that take the value 1 for rows with at least one count for each category
  all_visits <- all_visits %>%
    mutate(
      antibacterial = ifelse(antibacterial_counts > 0, 1, 0),
      betalactam = ifelse(betalactam_counts > 0, 1, 0),
      fluoroquinolone = ifelse(fluoroquinolone_counts > 0, 1, 0),
      sulfonamide = ifelse(sulfonamide_counts > 0, 1, 0),
      antibacterial_other = ifelse(antibacterial_other_counts > 0, 1, 0),
      antimalarial = ifelse(antimalarial_counts > 0, 1, 0),
      antiparasitic = ifelse(antiparasitic_counts > 0, 1, 0),
      antiviral = ifelse(antiviral_counts > 0, 1, 0),
      anymed = ifelse(anymed_counts > 0, 1, 0)
    )
  
#temp = all_visits %>%  select(id, nmed1:nmed7, antibacterial,betalactam, fluoroquinolone, sulfonamide, antibacterial_other, antimalarial, antiparasitic, antiviral), checking if the code worked 

saveRDS(all_visits, paste0(data_path, "maternal_medication.RDS")) 

