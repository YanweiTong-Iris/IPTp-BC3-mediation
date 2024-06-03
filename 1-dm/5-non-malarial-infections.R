################################################################
# IPTp and child growth
# Code maternal non-malarial outcomes 
# Adapted from code by Jordan Lee
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# read in data ------------------------------------------------------------
BC3.all.visits = read_dta(paste0(box_path, "Maternal/", "BC-3 mothers all visits database FINAL.dta"))
BC3.expanded = read_dta(paste0(box_path, "Maternal/", "BC-3 mothers expanded database FINAL.dta"))

diagnostic.codes = read_excel(paste0(box_path, "Codes/", "Diagnostic_codes.xlsx"))
medication.codes = read_excel(paste0(box_path, "Codes/", "Medication_codes.xlsx"))


# Frequency of Diagnoses =======================================================

apply((BC3.all.visits[c('diag1', 'diag2', 'diag3', 'diag4',
                        'diag5', 'diag6', 'diag7')]), 2, table) # No diagnoses available for 'diag5', 'diag6', 'diag7'

diagnoses <- dplyr::select(BC3.all.visits, diag1:diag7)
diagnoses.lvls <- unique(unlist(diagnoses))

diagnoses.freq <- sapply(diagnoses,
                         function(x) table(factor(x, levels = diagnoses.lvls,
                                                  ordered=TRUE)))
diagnoses.freq <- data.frame(diagnoses.freq)
diagnoses.freq$cum_diagnoses = rowSums(diagnoses.freq, na.rm = TRUE)
diagnoses.freq.sorted <- rownames_to_column(diagnoses.freq[order(diagnoses.freq$cum_diagnoses,
                                                                 decreasing=TRUE),], "Code")
diagnoses.freq.sorted <- transform(diagnoses.freq.sorted, Code = as.numeric(Code))

combined.diagnoses <- diagnoses.freq.sorted
combined.diagnoses <- merge(diagnoses.freq.sorted, diagnostic.codes, by = "Code", all = FALSE)
combined.diagnoses <- combined.diagnoses[order(combined.diagnoses$cum_diagnoses, decreasing=TRUE),] 

table(diagnoses.freq.sorted$Code)
table(combined.diagnoses$Code) # Diagnostic Code '12' does not have a description.

write_xlsx(combined.diagnoses, 
           path = paste0(box_path, "Codes/", "Combined Diagnoses.xlsx")) # Exporting combined_diagnoses data frame to Excel spreadsheet

# Creating the Different Classifications of Diagnoses ==========================

BC3.all.visits$febrile.illness <- ifelse(BC3.all.visits$temp>=38.0 | BC3.all.visits$fever==1, 1, 0)
sum(is.na(BC3.all.visits$febrile.illness)==TRUE)

summary(freqlist(~febrile.illness+temp+fever+febrile, data = BC3.all.visits, addNA = TRUE))

BC3.all.visits$malaria.diagnosed <- 0 # Assuming participant was not diagnosed with malaria at each clinical visit

BC3.all.visits$malaria.diagnosed[BC3.all.visits$diag1==202 | BC3.all.visits$diag2==202] <- 1 
# Diagnostic Code for Malaria = '202'

table(BC3.all.visits$malaria.diagnosed)

summary(freqlist(~mstatus+malaria.diagnosed, data = BC3.all.visits, addNA = TRUE))
filter(BC3.all.visits, malaria.diagnosed==1, mstatus==0)
filter(BC3.all.visits, malaria.diagnosed==0, mstatus==1)
BC3.all.visits$malaria.diagnosed[BC3.all.visits$mstatus == 0 & BC3.all.visits$malaria.diagnosed == 1] <- 0
BC3.all.visits$malaria.diagnosed[BC3.all.visits$mstatus == 1 & BC3.all.visits$malaria.diagnosed == 0] <- 1

# Observations:

# ID 30524 and 30561 were diagnosed with malaria but had 'mstatus' = 0
# They were febrile but had a negative blood smear --> They were treated for malaria, however.

# ID 30709 was not diagnosed with malaria but had 'mstatus' = 1

# Diagnostic code does not match 'mstatus' --> 'mstatus' takes precedence because data have been embedded already 
# 4 times where participants were diagnosed with malaria but had missing data on 'mstatus' --> All post-partum visits, so will not be included in final analysis

summary(freqlist(~mstatus+febrile+malaria.diagnosed, data = BC3.all.visits, addNA = TRUE))
summary(freqlist(~BSdich+malaria.diagnosed+mstatus+febrile.illness, data = BC3.all.visits, addNA = TRUE))
val_lab(BC3.all.visits$BSdich)

BC3.all.visits$non.malarial.febrile = ifelse(BC3.all.visits$mstatus!=0, 0,
                                             ifelse(BC3.all.visits$febrile==1, 1, 0))
table(BC3.all.visits$fever)
summary(freqlist(~mstatus+temp+fever+febrile+non.malarial.febrile, data = BC3.all.visits, addNA = TRUE))
summary(freqlist(~mstatus+febrile+non.malarial.febrile, data = BC3.all.visits, addNA = TRUE))
summary(freqlist(~febrile.illness+malaria.diagnosed+non.malarial.febrile, data = BC3.all.visits, addNA = TRUE))

# Observations: 

# 288 observations where participant had no presence of malaria but WAS febrile
# 674 observations where participant had no presence of malaria but had missing data on whether they were febrile

missing.febrile <- BC3.all.visits %>%
  dplyr::select(id, mstatus, temp, fever, febrile) %>%
  filter(is.na(febrile)==TRUE)

summary(freqlist(~mstatus+temp+fever+febrile, data = missing.febrile, addNA = TRUE))

### For above observations, temperature number must be above threshold (>= 38.0) OR presence of fever ('fever'==1).
### Override missing data for 'febrile' to '0' or '1' depending on whether participant meets either category.
### However, none of the observations contain such data.

respiratory.codes <- combined.diagnoses %>%
  filter(Code %in% c(183, 25, 174, 142, 139, 121, 122))
sum(respiratory.codes$cum_diagnoses)
# Observation: Total number of diagnoses = 1556

genitourinary.codes <- combined.diagnoses %>%
  filter(Code %in% c(184, 192, 168, 259, 148, 258, 260, 263))
sum(genitourinary.codes$cum_diagnoses)
# Observation: Total number of diagnoses = 820

gastrointestinal.codes <- combined.diagnoses %>%
  filter(Code %in% c(47, 43, 62, 46, 23, 44, 72, 196))
sum(gastrointestinal.codes$cum_diagnoses)
# Observation: Total number of diagnoses = 368

integumentary.codes <- combined.diagnoses %>%
  filter(Code %in% c(172, 1, 48, 206, 22, 207, 333, 20, 101, 111, 120, 305, 339, 361))
sum(integumentary.codes$cum_diagnoses)
# Observation: Total number of diagnoses = 128

abdpain.code <- combined.diagnoses %>%
  filter(Code %in% c(194))
sum(abdpain.code$cum_diagnoses)
# Observation: Total number of diagnoses = 611

headache.code <- combined.diagnoses %>%
  filter(Code %in% c(69))
sum(headache.code$cum_diagnoses)
# Observation: Total number of diagnoses = 654

other.pain.code <- combined.diagnoses %>%
  filter(Code %in% c(123))
sum(other.pain.code$cum_diagnoses)
# Observation: Total number of diagnoses = 593

BC3.all.visits <- BC3.all.visits %>%
  rowwise %>%
  mutate(respiratory.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% respiratory.codes$Code), na.rm=T),
         genitourinary.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% genitourinary.codes$Code), na.rm=T),
         gastrointestinal.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% gastrointestinal.codes$Code), na.rm=T),
         integumentary.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% integumentary.codes$Code), na.rm=T),
         abdpain.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% abdpain.code$Code), na.rm=T),
         headache.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% headache.code$Code), na.rm=T),
         other.pain.diagnosed = rowSums(across(diag1:diag4, ~ .x %in% other.pain.code$Code), na.rm=T)) %>%
  relocate(febrile.illness:other.pain.diagnosed, .before = nmed1)

# ==============================================================================================================
# New code beyond what Jordan originally wrote 

# Calculate gestational age at each measure  =======================================================
indiv <- read_dta(paste0(box_path, "Maternal/", "Maternal BC3 individual level database FINAL.dta")) %>% 
  dplyr::select(id, dod, finaledd, GAdelivery) %>% 
  mutate(date_preg_start = dod - GAdelivery*7)

merged <- left_join(BC3.all.visits, indiv, by = "id") %>% 
  mutate(GAweeks = as.numeric((date - date_preg_start)/7)) 

# create gestational age categories 
merged <- merged %>% mutate(GA_visit_cat = case_when(
  GAweeks <21 ~ "0-20 weeks",
  GAweeks >=21 & GAweeks <=29 ~ "21-28 weeks",
  GAweeks >=29 & GAweeks <=37 ~ "29-36 weeks",
  GAweeks >=37 ~ "37+ weeks"
))

merged = merged %>%
  filter(id %in% dfz$motherid)


saveRDS(merged, file = paste0(data_path, "non_malarial_infections.RDS"))
View(merged)
