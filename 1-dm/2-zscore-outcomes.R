################################################################
# IPTp and child growth
# Script that generates data that will be used to calculate 
# incidence and prevalence
################################################################


rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

#Loading the dataset 
zscores = readRDS(paste0(data_path, "cleaned_zscores_data.RDS"))

#Recategorizing age into the following categories 
zscores = zscores %>% mutate(agecat_birth = case_when(
  age == 0 ~ "Birth",
  age >0 & age <= 3 ~ "1 day-3 months",
  age >3 & age <= 6 ~ ">3-6 months",
  age >6 & age <= 9 ~ ">6-9 months",
  age >9 & age <= 12 ~ ">9-12 months")) %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "1 day-3 months", ">3-6 months",
    ">6-9 months", ">9-12 months"
  )))

saveRDS(zscores, paste0(data_path, "cleaned_zscores_agecat_birth.RDS"))

#Creating variables for moderate to severe and severe stunting and wasting; underweight 
zscores_data = zscores %>% mutate(haz_ms_stunt = ifelse(haz < -2 , 1, 0), 
                                     haz_s_stunt = ifelse(haz < - 3, 1, 0), 
                                     whz_ms_waste = ifelse(whz < -2, 1, 0), 
                                     whz_s_waste = ifelse(whz < -3, 1, 0),
                                     waz_underwt = ifelse(waz < -2, 1, 0),
                                     maternal_agecat = case_when(enrollage <20 ~ "less than 20",
                                                                 enrollage<25 ~ "20-24",
                                                                 enrollage<30 ~ "25-29",
                                                                 enrollage>29 ~ "30+")#recategorizing maternal age categories for sparsity 
)  

# Read the SGA csv file from Jordan
SGA_data <- read.csv((paste0(data_path, "Output for INTERGROWTH21_BW for gest age.csv"))) 

# Merge analysis data and SGA data
merged_data <- merge(zscores_data, SGA_data, by.x = "id", by.y = "Id", all.x = TRUE)


saveRDS(merged_data, paste0(data_path, "outcome_zscores_covariates.RDS"))


# Scramble the treatment arm
set.seed(10)
shuffled_Txarm <- sample(zscores_data$Txarm)
zscores_data$rand_Txarm <- shuffled_Txarm[match(zscores_data$id, unique(zscores_data$id))]
zscores_data$rand_Txarm <- relevel(factor(zscores_data$rand_Txarm), "SP")

zscores_data = zscores_data %>% 
  mutate(APdichenroll = factor(APdichenroll), educdich = factor(educdich),
         placentalmal = factor(placentalmal), placentalBSdich = factor(placentalBSdich), 
         placentalLAMPdich = factor(placentalLAMPdich), BSdichenroll = factor(BSdichenroll)) %>%
  mutate(enrollage_binary = factor(ifelse(enrollage< median(enrollage), 0 , 1))) %>%
  mutate(enrollage_binary_30 = factor(ifelse(enrollage< 30, 0 , 1))) %>%
  mutate(wealth_binary = factor(ifelse(wealthcat=="Least poor", 0, 1)))

#This dataset will be used to calculate prevalence and incidence in subsequent scripts 
saveRDS(zscores_data, paste0(data_path, "outcome_zscores_data.RDS"))
#View(readRDS(paste0(data_path, "outcome_zscores_data.RDS")))
