################################################################
# IPTp and child growth
# Table of total effects on incidence of binary outcomes
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Read in aim1 incidence ratio data 
incidence_ratio = readRDS(paste0(results_path, "aim1-stratified/aim1_incidence_3mo_results_stratified.RDS")) %>% 
  filter(is.na(modifier_level), gravidity_strata == "all") %>% 
  rename(agecat_birth = age_group) %>% 
  dplyr::select(outcome, agecat_birth, point_estimate, lower_95CI, upper_95CI) %>% 
  distinct() %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "1 day-3 months", ">3-6 months", ">6-9 months",">9-12 months"
  )))

# Read in quarterly incidence data 
incidence = readRDS(paste0(data_path, "analysis_data_incidence_quarterly.RDS"))[1:12]
# Variable names
current_names <- names(incidence)

# Renaming variables for ease of using dynamic variable names 
new_names <- gsub("^atrisk_", "atrisk_incident_", current_names)

# Replacing the incidence dataset with new names using setNames()
incidence <- setNames(incidence, new_names)


# Read in outcome zscores data to add Txarm column to incidence data
zscore_data = readRDS(paste0(data_path, "outcome_zscores_data.RDS")) %>% 
  dplyr::select(id, Txarm) %>% 
  group_by(id) %>% #since the outcome zscores data is long
  summarise(Txarm = first(Txarm)) 

# Merge Txarm to incidence dataframe 
incidence <- merge(incidence, zscore_data, by = "id", all.x = TRUE)

# Function for calculating incidence of each outcome  
incident_function <- function(dataset, incident_vars) {
  results <- list() # creating an empty list called results 
  for (var in incident_vars) {
    
    atrisk_var <- paste0("atrisk_", var) 
    
    dataset_clean = dataset[!is.na(dataset[[atrisk_var]]), ]
    
    #grouping, filtering for atrisk == TRUE, and then calculating incidence and ci 
    temp_result <- dataset_clean %>%
      group_by(Txarm, agecat_birth) %>%
      filter(get(atrisk_var) == 1, !is.na(get(var))) %>%
      summarise(incident = sum(get(var), na.rm = TRUE),
                rows = n(),
                .groups = 'drop') %>%
      mutate(var_name = var) %>%
      rowwise() %>%
      do(data.frame(.,
                    binconf(.$incident, .$rows, method = "wilson"))) %>%
      ungroup()
    #appending the results in temp_results to the list
    results[[var]] <- temp_result
  }
  return(results)
}
#looping through the outcome variables and applying the function to the incidence dataset
incident_vars <- c("incident_haz_ms_stunt_agecat_birth",
                   "incident_haz_s_stunt_agecat_birth",
                   "incident_whz_ms_waste_agecat_birth",
                   "incident_whz_s_waste_agecat_birth",
                   "incident_waz_underwt_agecat_birth")
results <- incident_function(incidence, incident_vars)

incidence_results <- dplyr::bind_rows(results) 

# Pivot the dataset to create a table that matches the table shell 
wide_incidence_results <- incidence_results %>%
  pivot_wider(
    names_from = Txarm,
    values_from = c(incident, rows, PointEst, Lower, Upper),
    names_sep = "_"
  ) %>% 
  rename(outcome = var_name)


# Merge the datasets for incidence ratio and incidence 
merged_data <- merge(incidence_ratio, wide_incidence_results, by = c("outcome", "agecat_birth")) %>%
  mutate(across(c(PointEst_DP, PointEst_SP, Lower_DP, Upper_DP, Lower_SP, Upper_SP), ~ . * 100))

# Format the table 
merged_data <- merged_data %>% 
  rename("Outcome" = outcome, "Age category" = agecat_birth, "N DP" = rows_DP, "N SP" = rows_SP) %>% 
  mutate("Incidence ratio (95% CI)" = paste0(sprintf("%.2f", point_estimate), " (", 
                             sprintf("%.2f", lower_95CI), ", ", sprintf("%.2f", upper_95CI), ")"), 
         "Incidence (95% CI) SP" = paste0(sprintf("%.1f", PointEst_SP), " (", sprintf("%.1f", Lower_SP), ", ", sprintf("%.1f", Upper_SP), ")"), 
         "Incidence (95% CI) DP" = paste0(sprintf("%.1f", PointEst_DP), " (", sprintf("%.1f", Lower_DP), ", ", sprintf("%.1f", Upper_DP), ")")) %>% 
  dplyr::select("Outcome", "Age category", "N DP", "Incidence (95% CI) DP", "N SP", "Incidence (95% CI) SP", "Incidence ratio (95% CI)") 

# Clean the Outcome column 
merged_data <- merged_data %>%
  mutate("Outcome" = case_when(
    Outcome == "incident_haz_ms_stunt_agecat_birth" ~ "Moderate to severe stunting",
    Outcome == "incident_haz_s_stunt_agecat_birth" ~ "Severe stunting",
    Outcome == "incident_whz_ms_waste_agecat_birth" ~ "Moderate to severe wasting",
    Outcome == "incident_whz_s_waste_agecat_birth" ~ "Severe wasting",
    Outcome == "incident_waz_underwt_agecat_birth" ~ "Underweight"
  ))  %>% 
  mutate(Outcome = factor(Outcome, levels = c(
    "Moderate to severe stunting","Severe stunting",
    "Moderate to severe wasting","Severe wasting","Underweight"
  ))) %>% 
  arrange(`Outcome`, `Age category`)


write.csv(merged_data, here::here(tables_path, "table_aim1_incidence.csv")) 



#-----------------------------------------------------
# Primigravidae
#-----------------------------------------------------

# Read in aim1 incidence ratio data (primigravidae)
incidence_ratio_single = readRDS(paste0(results_path, "aim1-stratified/aim1_incidence_3mo_results_stratified.RDS")) %>% 
  filter(is.na(modifier_level), gravidity_strata == "single") %>% 
  rename(agecat_birth = age_group) %>% 
  dplyr::select(outcome, agecat_birth, point_estimate, lower_95CI, upper_95CI) %>% 
  distinct() %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "1 day-3 months", ">3-6 months", ">6-9 months",">9-12 months"
  )))


#subset primigravidae
incidence_single <- incidence %>% filter(Gravidity == 1)

results_single <- incident_function(incidence_single, incident_vars)

incidence_results_single <- dplyr::bind_rows(results_single) 

# Pivot the dataset to create a table that matches the table shell 
wide_incidence_results_single <- incidence_results_single %>%
  pivot_wider(
    names_from = Txarm,
    values_from = c(incident, rows, PointEst, Lower, Upper),
    names_sep = "_"
  ) %>% 
  rename(outcome = var_name)


# Merge the datasets for incidence ratio and incidence 
merged_data_single <- merge(incidence_ratio_single, wide_incidence_results_single, by = c("outcome", "agecat_birth")) %>%
  mutate(across(c(PointEst_DP, PointEst_SP, Lower_DP, Upper_DP, Lower_SP, Upper_SP), ~ . * 100))

# Format the table 
merged_data_single = data_format(merged_data_single)

write.csv(merged_data_single, here::here(tables_path, "table_aim1_incidence_primigravidae.csv")) 


#-----------------------------------------------------
# Multigravidae
#-----------------------------------------------------

# Read in aim1 incidence ratio data (multigravidae)
incidence_ratio_multi = readRDS(paste0(results_path, "aim1-stratified/aim1_incidence_3mo_results_stratified.RDS")) %>% 
  filter(is.na(modifier_level), gravidity_strata == "multi") %>% 
  rename(agecat_birth = age_group) %>% 
  dplyr::select(outcome, agecat_birth, point_estimate, lower_95CI, upper_95CI) %>% 
  distinct() %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "1 day-3 months", ">3-6 months", ">6-9 months",">9-12 months"
  )))


#subset multigravidae
incidence_multi <- incidence %>% filter(Gravidity > 1)

results_multi <- incident_function(incidence_multi, incident_vars)

incidence_results_multi <- dplyr::bind_rows(results_multi) 

# Pivot the dataset to create a table that matches the table shell 
wide_incidence_results_multi <- incidence_results_multi %>%
  pivot_wider(
    names_from = Txarm,
    values_from = c(incident, rows, PointEst, Lower, Upper),
    names_sep = "_"
  ) %>% 
  rename(outcome = var_name)


# Merge the datasets for incidence ratio and incidence 
merged_data_multi <- merge(incidence_ratio_multi, wide_incidence_results_multi, by = c("outcome", "agecat_birth")) %>%
  mutate(across(c(PointEst_DP, PointEst_SP, Lower_DP, Upper_DP, Lower_SP, Upper_SP), ~ . * 100))

# Format the table 
merged_data_multi = data_format(merged_data_multi)

write.csv(merged_data_multi, here::here(tables_path, "table_aim1_incidence_multigravidae.csv")) 


