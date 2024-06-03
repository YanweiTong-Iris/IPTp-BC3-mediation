################################################################
# IPTp and child growth
# Create binary indicators for incidence variables
# Within different age categories
################################################################
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

zscore_data = readRDS(paste0(data_path, "outcome_zscores_data.RDS"))

# Subset to relevant variables; create age categories  --------------------
incidence_data = zscore_data %>%  dplyr::select(
  id, date, agedays, agecat, agecat_birth, agemonthcat,
  haz_ms_stunt, haz_s_stunt, whz_ms_waste, whz_s_waste, waz_underwt
)  %>% 
  # create age categories
  mutate(agecat = factor(agecat, levels = c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>%  
  
  mutate(age_1_12 = case_when(
    agedays >= 1 & agedays <=365.25 ~ "1 day- 12 months",
    agedays <1 ~ "birth")) %>% #any incidence after birth 
  
  mutate(age_6_12_month = case_when(
    agedays <1 ~ "birth",
    agedays >= 1 & agedays <182.625 ~ "1 day- <6 months",
    agedays >= 182.625 & agedays <=365.25 ~ "6-12 months")) %>% #incidence below 6 months and 6 months and above 
  
  mutate(age_1_12 = factor(age_1_12, levels = c("birth", "1 day- 12 months")),
         age_6_12_month = factor(age_6_12_month, levels = c("birth", "1 day- <6 months", "6-12 months")),
         agemonthcat = as.factor(agemonthcat))

# List containing all outcomes and ages    --------------------
outcome_list = c("haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", "waz_underwt")
agecat_list = c("agemonthcat", "agecat_birth", "age_1_12", "age_6_12_month") 

# Function for calculating the first occurrence of each outcome in outcome_list     --------------------
# within different age categories 
find_index_of_first_outcome = function(outcome, agecat_col){ 
  incidence_data %>% 
    arrange(id, agedays) %>%
    group_by(id) %>% 
    mutate(Index = which(!!sym(outcome) == 1)[1], 
           ROW = row_number(), 
           FIRST = (Index == ROW), 
           FIRST = ifelse(is.na(FIRST), FALSE, FIRST))  %>% 
    ungroup() %>%
    dplyr::select(!!glue::glue("first_{outcome}_{agecat_col}") := FIRST) ## squiggly bracket are for dynamic strings
}

# Looping through the function using outcome_list and age category list    --------------------
incidence_list = list()
for(i in 1:length(agecat_list)){
  incidence_list[[i]] = bind_cols(lapply(outcome_list, function(x) 
    find_index_of_first_outcome(outcome = x, 
                                agecat_col = agecat_list[i]))) 
}

incidence_all_ages = incidence_data %>% bind_cols(incidence_list)


# Function for defining at risk population   --------------------
incidence_outcomes = function(outcome, agecat_col) { 
  incidence_all_ages %>%
    arrange(id, !!sym(agecat_col)) %>% 
    group_by(id, !!sym(agecat_col)) %>% 
    dplyr::summarize(incident_case = sum(!!sym(glue::glue("first_{outcome}_agecat_birth")), na.rm =TRUE)) %>%
    mutate(cum_outcome = stats::lag(cumsum(incident_case), k=1),
           cum_outcome = replace_na(cum_outcome, 0),
           at_risk = cum_outcome == 0) %>% 
    mutate(at_risk= ifelse(incident_case == 1, T, at_risk)) %>% 
    rename(!!glue::glue("incident_{outcome}_{agecat_col}") := incident_case,
           !!glue::glue("atrisk_{outcome}_{agecat_col}") := at_risk) %>% 
    dplyr::select(id, !!sym(agecat_col), 
           !!glue::glue("incident_{outcome}_{agecat_col}"), 
           !!glue::glue("atrisk_{outcome}_{agecat_col}"))# := string on the left will take the value of the right side
}

# Looping through the function using outcome_list    --------------------
incidence_data_age_1_12 = lapply(outcome_list, incidence_outcomes, agecat_col = "age_1_12") %>% 
  Reduce(left_join, .) #used to join datasets with unequal lengths 

incidence_data_age_6_12_month = lapply(outcome_list, incidence_outcomes, agecat_col = "age_6_12_month") %>% 
  Reduce(left_join, .) #used to join datasets with unequal lengths 

incidence_data_agecat_birth = lapply(outcome_list, incidence_outcomes, agecat_col = "agecat_birth") %>% 
  Reduce(left_join, .) #used to join datasets with unequal lengths 

incidence_data_agecat_month = lapply(outcome_list, incidence_outcomes, agecat_col = "agemonthcat") %>% 
  Reduce(left_join, .) #used to join datasets with unequal lengths 

# Save data    --------------------
saveRDS(incidence_data_age_1_12, paste0(data_path, "incidence_data_annual.RDS"))
saveRDS(incidence_data_age_6_12_month, paste0(data_path, "incidence_data_biannual.RDS"))
saveRDS(incidence_data_agecat_birth, paste0(data_path, "incidence_data_quarterly.RDS"))
saveRDS(incidence_data_agecat_month, paste0(data_path, "incidence_data_monthly.RDS"))

