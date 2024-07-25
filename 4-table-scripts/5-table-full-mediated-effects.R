################################################################
# IPTp and child growth
# Full ACME and ADE table 
# (quarterly zscore, velocity, and incidence; with and without interaction)
# Last updated: July 19, 2024
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library('data.table')
library(purrr)

#Loading the data
zscore_results = readRDS(paste0(results_path, "aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS"))
velocity_results = readRDS(paste0(results_path, "aim2-stratified/aim2_single_mediator_velocity_results_3mo.RDS"))
incidence_results = readRDS(paste0(results_path, "aim2-stratified/aim2_single_mediator_incidence_results_3mo.RDS"))


full_results = rbind(zscore_results[c(1:9, 11:25)], velocity_results[c(1:9, 11:25)], incidence_results[c(1:9, 11:25)]) %>%
  dplyr::select(-time_unit, -ACME_p_val) %>% 
  mutate_if(is.numeric, ~ ifelse(abs(.) > 1000, NA, .)) %>%
  mutate(across(c(total_effect, total_effect_CI, total_effect_lower_CI, total_effect_upper_CI),
      ~ ifelse(if_any(c(total_effect, total_effect_lower_CI, total_effect_upper_CI), is.na), NA, .))) %>%
  mutate(across(c(ACME_average, ACME_average_CI, ACME_average_lower_CI, ACME_average_upper_CI),
                ~ ifelse(if_any(c(ACME_average, ACME_average_lower_CI, ACME_average_upper_CI), is.na), NA, .))) %>%
  mutate(across(c(ADE_average, ADE_average_CI, ADE_average_lower_CI, ADE_average_upper_CI),
                ~ ifelse(if_any(c(ADE_average, ADE_average_lower_CI, ADE_average_upper_CI), is.na), NA, .))) %>%
  dplyr::select(-total_effect_lower_CI, -total_effect_upper_CI, -ACME_average_lower_CI, -ACME_average_upper_CI,
                -ADE_average_lower_CI, -ADE_average_upper_CI)

View(full_results)


write.csv(full_results, here::here(tables_path, "table_full_mediation_results.csv")) 

