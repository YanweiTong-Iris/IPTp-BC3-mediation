################################################################
# IPTp and child growth
# Table of mediator-outcome effects
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

MO_incidence_3mo = readRDS(paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_3mo_stratified.RDS")) %>%
  mutate_if(is.numeric, round, digits=2) %>%
  dplyr::select(independent_variable, dependent_variable, N_from_analysis, age_group, gravidae, 
                point_estimate, lower_95CI,upper_95CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>%
  filter(gravidae == "all") %>%
  filter(dependent_variable %in% c("incident_haz_ms_stunt_agecat_birth", "incident_whz_ms_waste_agecat_birth")) %>%
  filter(independent_variable %in% c(
    "anemia_28binary", "gestational_weightchange",
    "placentalmal", "preterm",
    "birthlength", "birthweight",
    "birthweight_kg"
  )) %>%
  mutate(mediator_label = case_when(
    independent_variable == "anemia_28binary" ~ "Anemia",
    independent_variable == "gestational_weightchange" ~ "Gestational weight change (kg)",
    independent_variable == "placentalmal" ~ "Placental malaria",
    independent_variable == "preterm" ~ "Pre-term birth",
    independent_variable == "birthlength" ~ "Birth length (cm)",
    independent_variable== "birthweight" ~ "Birth weight (g)",
    independent_variable== "birthweight_kg" ~ "Birth weight (kg)"
  )) %>% 
  mutate(
    mediator_label = factor(mediator_label, levels = c(
      "Anemia", 
      "Gestational weight change (kg)",
      "Placental malaria",
      "Pre-term birth",
      "Birth length (cm)",
      "Birth weight (g)",
      "Birth weight (kg)"
    ))) %>%
  mutate(outcome = dependent_variable) %>%
  mutate(effect_CI = paste0(point_estimate, " (", lower_95CI, ", ", upper_95CI, ")")) %>%
  dplyr::select(gravidae, outcome, mediator_label, age_group, N_from_analysis, effect_CI) %>%
  mutate(N_from_analysis = as.integer(N_from_analysis)) 


MO_incidence_3mo = with(MO_incidence_3mo, MO_incidence_3mo[order(outcome, mediator_label, age_group),])

write.csv(MO_incidence_3mo,  paste0(tables_path,"table_MO_incidence_3mo.csv"))
