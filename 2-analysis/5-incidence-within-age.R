################################################################
# IPTp and child growth
# Calculate incidence within age and covariate levels
################################################################
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Read in data   --------------------
incidence_data_age_1_12 = readRDS(paste0(data_path, "incidence_data_annual.RDS"))
incidence_data_age_6_12_month = readRDS(paste0(data_path, "incidence_data_biannual.RDS"))
incidence_data_agecat_birth = readRDS(paste0(data_path, "incidence_data_quarterly.RDS"))
incidence_data_agecatmonth = readRDS(paste0(data_path, "incidence_data_monthly.RDS"))

zscore_data = readRDS(paste0(data_path, "outcome_zscores_data.RDS"))

# List containing all outcomes and ages   --------------------
outcome_list = c("haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", "waz_underwt")

# List of covariates of interest    --------------------
covar_list = c("LBW", "GAcat", "mparasitestatus", "sex", "maternal_agecat", "Txarm", 
               "Gravidcat", "anemiaenroll","LBWpreterm", "anyHP", "BSdichenroll", "APdichenroll")
covar_df = zscore_data %>% dplyr::select(all_of(c("id", covar_list))) %>% distinct()

# Merge in covariates of interest    --------------------
incidence_data_age_1_12 = left_join(incidence_data_age_1_12, covar_df, by = "id")
incidence_data_age_6_12_month = left_join(incidence_data_age_6_12_month, covar_df, by = "id")
incidence_data_agecat_birth = left_join(incidence_data_agecat_birth, covar_df, by = "id")
incidence_data_agecatmonth = left_join(incidence_data_agecatmonth, covar_df, by = "id")


# Calculates incidence and CIs for each age, covariate, and outcome   --------------------
## Age in annual categories  --------------------

#Creating a list that will contain CIs 
inc_results_covar_list1 = list()
inc_results_outcome_list1 = list()

covar_list1 = c("age_1_12", covar_list)

for (o in 1:length(outcome_list)){
  for(c in 1:length(covar_list1)){
    inc_results_covar_list1[[c]] = lapply(levels(incidence_data_age_1_12$age_1_12), function(x) 
      get_Inc_CI(data = incidence_data_age_1_12, 
                 outcome_col = outcome_list[o], 
                 strat_varname = covar_list1[c],
                 agecat_col = "age_1_12",
                 age_level = x) 
    ) %>% bind_rows()
  }
  inc_results_outcome_list1[[o]] = inc_results_covar_list1 %>% bind_rows()
}

inc_results_annual = inc_results_outcome_list1 %>% bind_rows()

saveRDS(inc_results_annual, paste0(results_path, "incidence-age/incidence_estimates_annual.RDS")) 

## Age in biannual categories  --------------------

#Creating a list that will contain CIs 
inc_results_covar_list2 = list()
inc_results_outcome_list2 = list()

covar_list2 = c("age_6_12_month", covar_list)

for (o in 1:length(outcome_list)){
  for(c in 1:length(covar_list2)){
    inc_results_covar_list2[[c]] = lapply(levels(incidence_data_age_6_12_month$age_6_12_month), function(x) 
      get_Inc_CI(data = incidence_data_age_6_12_month, 
                 outcome_col = outcome_list[o], 
                 strat_varname = covar_list2[c],
                 agecat_col = "age_6_12_month",
                 age_level = x) 
    ) %>% bind_rows()
  }
  inc_results_outcome_list2[[o]] = inc_results_covar_list2 %>% bind_rows()
}

inc_results_biannual = inc_results_outcome_list2 %>% bind_rows()

saveRDS(inc_results_biannual, paste0(results_path, "incidence-age/incidence_estimates_biannual.RDS")) 


## Age in quarterly categories  --------------------

#Creating a list that will contain CIs 
inc_results_covar_list3 = list()
inc_results_outcome_list3 = list()

covar_list3 = c("agecat_birth", covar_list)

for (o in 1:length(outcome_list)){
  for(c in 1:length(covar_list3)){
    inc_results_covar_list3[[c]] = lapply(levels(incidence_data_agecat_birth$agecat_birth), function(x) 
      get_Inc_CI(data = incidence_data_agecat_birth, 
                 outcome_col = outcome_list[o], 
                 strat_varname = covar_list3[c],
                 agecat_col = "agecat_birth",
                 age_level = x) 
    ) %>% bind_rows()
  }
  inc_results_outcome_list3[[o]] = inc_results_covar_list3 %>% bind_rows()
}

inc_results_quarterly = inc_results_outcome_list3 %>% bind_rows()

saveRDS(inc_results_quarterly, paste0(results_path, "incidence-age/incidence_estimates_quarterly.RDS")) 


## Age in monthly categories  --------------------

#Creating a list that will contain CIs 
inc_results_covar_list4 = list()
inc_results_outcome_list4 = list()

covar_list4 = c("agemonthcat", covar_list)

for (o in 1:length(outcome_list)){
  for(c in 1:length(covar_list4)){
    inc_results_covar_list4[[c]] = lapply(levels(incidence_data_agecatmonth$agemonthcat), function(x) 
      get_Inc_CI(data = incidence_data_agecatmonth, 
                 outcome_col = outcome_list[o], 
                 strat_varname = covar_list4[c],
                 agecat_col = "agemonthcat",
                 age_level = x) 
    ) %>% bind_rows()
  }
  inc_results_outcome_list4[[o]] = inc_results_covar_list4 %>% bind_rows()
}

inc_results_monthly = inc_results_outcome_list4 %>% bind_rows()

saveRDS(inc_results_monthly, paste0(results_path, "incidence-age/incidence_estimates_monthly.RDS")) 

