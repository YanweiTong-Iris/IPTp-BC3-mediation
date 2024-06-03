################################################################
# IPTp and child growth
# Script for sensitivity analysis to check unmeasured 
# mediator-outcome confounding
# Last updated: May 27, 2024
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#prevent using scientific notations
options(scipen = 999)

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_zscore_quarterly = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS"))
data_incidence_3month = readRDS(paste0(data_path,"analysis_data_incidence_quarterly.RDS"))

result_single_mediator_zscore_3mo = 
  readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS")) %>%
  filter(interaction == 0) %>%
  dplyr::select(mediator, outcome, age_group, gravidae,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>% 
  filter(gravidae == "all") %>%
  filter(mediator %in% c("birthweight_kg", "birthlength", "SCF", "IL18", "CDCP1", "CD6", "OPG", "DNER"))

result_single_mediator_incidence_3mo = readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_incidence_results_3mo.RDS")) %>%
  filter(interaction == 0) %>%
  dplyr::select(mediator, outcome, age_group, gravidae,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>% 
  filter(gravidae == "all") %>%
  filter(mediator %in% c("LBW","preterm","birthweight_kg", "birthlength", "SCF"))


#--------------------------------------------------------
# Mediational e-value for ACME and ADE (incidence outcome)
#--------------------------------------------------------

results_single_mediator_ms_stunt = subset(result_single_mediator_incidence_3mo, 
                                          outcome == "incident_haz_ms_stunt_agecat_birth")
results_single_mediator_ms_waste = subset(result_single_mediator_incidence_3mo, 
                                          outcome == "incident_whz_ms_waste_agecat_birth")


# Function to compute the E-value for RR > 1
compute_e_val_rr_gt_1 <- function(RR, LL, UL) {
  e_value <- RR + sqrt(RR * (RR - 1))
  e_value_bound <- ifelse(LL <= 1, 1, LL + sqrt(LL * (LL - 1)))
  return(list(E_value = e_value, E_value_bound = e_value_bound))
}

# Function to compute the E-value for RR < 1
compute_e_val_rr_lt_1 <- function(RR, LL, UL) {
  RR_star <- 1 / RR
  e_value <- RR_star + sqrt(RR_star * (RR_star - 1))
  e_value_bound <- ifelse(UL >= 1, 1, {
    UL_star <- 1 / UL
    UL_star + sqrt(UL_star * (UL_star - 1))
  })
  return(list(E_value = e_value, E_value_bound = e_value_bound))
}

# Wrapper function to add E-value and bounds to the dataset
add_mediational_e_value <- function(data) {
  data <- data %>%
    rowwise() %>%
    mutate(
      ACME_e_val_result = ifelse(ACME_average > 1,
                                 list(compute_e_val_rr_gt_1(ACME_average, ACME_average_lower_CI, ACME_average_upper_CI)),
                                 list(compute_e_val_rr_lt_1(ACME_average, ACME_average_lower_CI, ACME_average_upper_CI))),
      ADE_e_val_result = ifelse(ADE_average > 1,
                                list(compute_e_val_rr_gt_1(ADE_average, ADE_average_lower_CI, ADE_average_upper_CI)),
                                list(compute_e_val_rr_lt_1(ADE_average, ADE_average_lower_CI, ADE_average_upper_CI))),
      ACME_mediational_eval = ACME_e_val_result$E_value,
      ACME_mediational_eval_bound = ACME_e_val_result$E_value_bound,
      ADE_mediational_eval = ADE_e_val_result$E_value,
      ADE_mediational_eval_bound = ADE_e_val_result$E_value_bound
    ) %>%
    ungroup() %>%
    dplyr::select(-ACME_e_val_result, -ADE_e_val_result)
  
  return(data)
}



results_single_mediator_ms_stunt_eval = add_mediational_e_value(results_single_mediator_ms_stunt) 
results_single_mediator_ms_waste_eval = add_mediational_e_value(results_single_mediator_ms_waste)


sensitivity_eval_incidence = rbind(results_single_mediator_ms_stunt_eval, results_single_mediator_ms_waste_eval)
saveRDS(sensitivity_eval_incidence, paste0(results_path,"sensitivity/sensitivity_eval_incidence.RDS"))


#--------------------------------------------------------
# Mediational e-value for ACME and ADE (zscore outcome)
#--------------------------------------------------------

results_single_mediator_laz = subset(result_single_mediator_zscore_3mo, outcome == "haz_quarter")
results_single_mediator_wlz = subset(result_single_mediator_zscore_3mo, outcome == "whz_quarter")


z_data_long <- data_zscore_quarterly %>%
  pivot_longer(cols = c(haz_quarter, whz_quarter), names_to = "outcome", values_to = "zscore")


sample_size_sd <- z_data_long %>%
  group_by(agecat_birth, outcome) %>%
  summarise(
    sample_size = sum(!is.na(zscore)), 
    outcome_sd = sd(zscore, na.rm = TRUE)
  )

sample_size_sd_laz <- sample_size_sd %>%
  filter(outcome == "haz_quarter") %>%
  dplyr::select(-outcome)

sample_size_sd_wlz <- sample_size_sd %>%
  filter(outcome == "whz_quarter") %>%
  dplyr::select(-outcome)

results_single_mediator_laz <- results_single_mediator_laz %>%
  left_join(sample_size_sd_laz, by = c("age_group" = "agecat_birth")) %>%
  #rename(sample_size = sample_size, outcome_sd = outcome_sd) %>%
  mutate(ACME_d = ACME_average / outcome_sd,
         standardized_ACME_SE = (ACME_average_upper_CI - ACME_average_lower_CI) / (3.92*outcome_sd),
         ADE_d = ADE_average / outcome_sd,
         standardized_ADE_SE = (ADE_average_upper_CI - ADE_average_lower_CI) / (3.92*outcome_sd)) %>%
  mutate(approx_ACME_RR = exp(0.91*ACME_d),
         approx_ACME_RR_lower_CI = exp(0.91*ACME_d - 1.78*standardized_ACME_SE),
         approx_ACME_RR_upper_CI = exp(0.91*ACME_d + 1.78*standardized_ACME_SE),
         approx_ADE_RR = exp(0.91*ADE_d),
         approx_ADE_RR_lower_CI = exp(0.91*ADE_d - 1.78*standardized_ADE_SE),
         approx_ADE_RR_upper_CI = exp(0.91*ADE_d + 1.78*standardized_ADE_SE))

results_single_mediator_wlz <- results_single_mediator_wlz %>%
  left_join(sample_size_sd_wlz, by = c("age_group" = "agecat_birth")) %>%
  mutate(ACME_d = ACME_average / outcome_sd,
         standardized_ACME_SE = (ACME_average_upper_CI - ACME_average_lower_CI) / (3.92*outcome_sd),
         ADE_d = ADE_average / outcome_sd,
         standardized_ADE_SE = (ADE_average_upper_CI - ADE_average_lower_CI) / (3.92*outcome_sd)) %>%
  mutate(approx_ACME_RR = exp(0.91*ACME_d),
         approx_ACME_RR_lower_CI = exp(0.91*ACME_d - 1.78*standardized_ACME_SE),
         approx_ACME_RR_upper_CI = exp(0.91*ACME_d + 1.78*standardized_ACME_SE),
         approx_ADE_RR = exp(0.91*ADE_d),
         approx_ADE_RR_lower_CI = exp(0.91*ADE_d - 1.78*standardized_ADE_SE),
         approx_ADE_RR_upper_CI = exp(0.91*ADE_d + 1.78*standardized_ADE_SE))
  
# Wrapper function to add approx E-value and bounds to the dataset
add_approx_mediational_e_value <- function(data) {
    data <- data %>%
      rowwise() %>%
      mutate(
        ACME_e_val_result = ifelse(approx_ACME_RR > 1,
                                   list(compute_e_val_rr_gt_1(approx_ACME_RR, approx_ACME_RR_lower_CI, approx_ACME_RR_upper_CI)),
                                   list(compute_e_val_rr_lt_1(approx_ACME_RR, approx_ACME_RR_lower_CI, approx_ACME_RR_upper_CI))),
        ADE_e_val_result = ifelse(approx_ADE_RR > 1,
                                  list(compute_e_val_rr_gt_1(approx_ADE_RR, approx_ADE_RR_lower_CI, approx_ADE_RR_upper_CI)),
                                  list(compute_e_val_rr_lt_1(approx_ADE_RR, approx_ADE_RR_lower_CI, approx_ADE_RR_upper_CI))),
        approx_ACME_mediational_eval = ACME_e_val_result$E_value,
        approx_ACME_mediational_eval_bound = ACME_e_val_result$E_value_bound,
        approx_ADE_mediational_eval = ADE_e_val_result$E_value,
        approx_ADE_mediational_eval_bound = ADE_e_val_result$E_value_bound
      ) %>%
      ungroup() %>%
      dplyr::select(-ACME_e_val_result, -ADE_e_val_result)

    return(data)
}

results_single_mediator_laz_eval = add_approx_mediational_e_value(results_single_mediator_laz)
results_single_mediator_wlz_eval = add_approx_mediational_e_value(results_single_mediator_wlz)
sensitivity_eval_zscore = rbind(results_single_mediator_laz_eval, results_single_mediator_wlz_eval)
saveRDS(sensitivity_eval_zscore, paste0(results_path,"sensitivity/sensitivity_eval_zscore.RDS"))




