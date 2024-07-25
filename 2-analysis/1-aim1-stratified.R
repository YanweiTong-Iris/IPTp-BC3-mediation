################################################################
# IPTp and child growth
# Script for Aim1: estimate the impact of IPTp on child growth
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


#--------------------------------------------------------
# create outcome, age, and modifier list 
#--------------------------------------------------------

outcome_zscore_quarter = c("haz_quarter", "whz_quarter", "waz_quarter")
outcome_incidence_3mo = c(
  "incident_haz_ms_stunt_agecat_birth",
  "incident_haz_s_stunt_agecat_birth",
  "incident_whz_ms_waste_agecat_birth",
  "incident_whz_s_waste_agecat_birth",
  "incident_waz_underwt_agecat_birth"
)

age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))

#--------------------------------------------------------
# aim1 analysis wrapper function
#--------------------------------------------------------

# Documentation: aim1_analysis
# Usage: aim1_analysis(data, time_unit, age_group, outcome, outcome_type, treatment, modifier)
# Description: Estimate the impact of IPTp on child growth by each age category (stratified by gravidity)
# For continuous outcomes, use GLM with a Gaussian family and identity link.
# For binary outcomes, use GLM with a Poisson family and log link.
# If the sample size is sufficient, investigate whether effects are modified by
# the modifiers below in both a relative and additive scale.
# Return an outcome data frame with outcome name, age group, N from analysis,
# modifier name (or “none” if no modifier), point estimate, SE, lower bound of 95% CI of the point estimate,
# upper bound of 95% CI of the point estimate, additive modifier level (or NA if no modifier),
# lower bound of 95% CI of the additive modifier level, upper bound of 95% CI of the additive modifier level,
# relative modifier level (or NA if no modifier), lower bound of 95% CI of the relative modifier level,
# upper bound of 95% CI of the relative modifier level

# Args/Options:
# time_unit: "quarter"
# age_group: 
#     - data_zscore_quarterly: agecat_birth
#     - data_incidence_3month: agecat_birth
# outcome: outcome name, as a string (w/ different time units at the end)
# continuous: "haz_quarter", "whz_quarter", "waz_quarter"
# binary: "incident_haz_ms_stunt", "incident_haz_s_stunt", "incident_whz_ms_waste", 
#         "incident_whz_s_waste", "incident_waz_underwt","SGA"
# outcome_type: type of outcome, as a string (either "continuous" or "binary")
# treatment: either "Txarm" or "rand_Txarm")
# modifier_name: name of effect modifier, as a string (default = NULL)
#   "enrollage" (maternal age), "gravidity_cat" (first pregnancy or not), 
#   "gestational_weightchange_cat" (weight difference at gestaional week 36 and 20), 
#   "mother_heightcat" (the mother's height measured during pregnancy),
#   "preterm" (pre-term delivery), "sex" (child sex), "anymalaria" (child malaria infection),
#   and "wet_season" (wet season or not, not available for biannual and annual datasets),
# Note: when gravidity = all, the model is adjusted for gravidity; 
# when the data is stratified by multi or primigravida, then it's not adjusted for gravidity

# Returns: the data frame as described above
# Output: prints the data frame described above

aim1_analysis <-
  function(data,
           time_unit,
           age_group,
           age_category,
           outcome,
           outcome_type,
           treatment = "Txarm",
           modifier = "NA") {
    
    print(paste0("outcome: ", outcome, ", age: ", age_group, ", modifier: ", modifier))
    
    data_age = data[data[[age_category]] == age_group, ]
    
    output_list = list()
    strata_list = list(
      "all" = c("single", "multi"),
      "single" = c("single"),
      "multi" = c("multi")
    )
    
    for (i in 1:3) {
      strata_name = names(strata_list)[[i]]
      strata = strata_list[[strata_name]]
      data_stratified = data_age %>% filter(gravidity_cat %in% strata)
      
      modifier_counts <- table(data_stratified[[modifier]])
      
      if (modifier == "NA") {
        formula <- as.formula(paste(outcome, "~", treatment))
      } else if (length(modifier_counts) > 1 && all(modifier_counts > 3)){
        data_stratified <- data_stratified[!is.na(data_stratified[[modifier]]), ]
        formula <-
          as.formula(paste(outcome, "~", treatment, "+", modifier, "+", treatment, ":", modifier))
      } else {
        message("One of the levels doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
      print(formula)
      
      if (outcome_type == "continuous") {
        model <- glm(formula, data = data_stratified, family = gaussian(link = "identity"), na.action = na.omit)
        
      } else if (outcome_type == "binary") {
        model <- glm(formula, data = data_stratified, family = poisson(link = "log"), na.action = na.omit)
      }
      
      eligible_N_data <- data_stratified[!is.na(data_stratified[[outcome]]), ]
      N <- length(unique(eligible_N_data$id))
      
      # model without interaction tx by gravidity
      estimates <- summary(model)$coefficients
      coefs = names(coef(model))
      treatment_row = as.numeric(grep("armDP$", coefs))
      modifier_row = as.numeric(grep(paste0("^", modifier), coefs))
      modifier_levels = grep(paste0("^", modifier), coefs, value = TRUE)
      interaction_row = as.numeric(grep(paste0("^", treatment, "DP:", modifier), coefs))
      
      point_estimate <- estimates[treatment_row, 1]
      point_estimate_modified <- ifelse(outcome_type == "continuous", point_estimate, exp(point_estimate))
      SE <- estimates[treatment_row, 2]
      lower_CI <- ifelse(outcome_type == "continuous", point_estimate - qnorm(0.975) * SE, exp(point_estimate - qnorm(0.975) * SE))
      upper_CI <- ifelse(outcome_type == "continuous", point_estimate + qnorm(0.975) * SE, exp(point_estimate + qnorm(0.975) * SE))
      
      # Assessment of effect modification
      if (modifier != "NA" && length(modifier_counts) > 1 && all(modifier_counts > 3)) {
        beta2 = coef(model)[modifier_row]
        beta3 = coef(model)[interaction_row]
        # Additive interaction level
        RERI = 1 + exp(beta3 + beta2 + point_estimate) - exp(point_estimate) - exp(beta2)
        addi_g <-
          as.formula(~ exp(x3 + x2 + x1) + 1 - exp(x1) - exp(x2))
        RERI_SE = deltamethod(g = addi_g,
                              mean = coef(model)[c(treatment_row, modifier_row, interaction_row)],
                              cov = vcov(model)[c(treatment_row, modifier_row, interaction_row), c(treatment_row, modifier_row, interaction_row)])
        RERI_CI <-
          cbind(RERI - qnorm(0.975) * RERI_SE, RERI + qnorm(0.975) * RERI_SE)
        # Relative (multiplicative) interaction level
        rel_int <- exp(beta3)
        rel_g = as.formula(~ exp(x1))
        rel_SE = deltamethod(g = rel_g,
                             mean = coef(model)[interaction_row],
                             cov = vcov(model)[interaction_row, interaction_row])
        rel_CI <-
          cbind(rel_int - qnorm(0.975) * rel_SE,
                rel_int + qnorm(0.975) * rel_SE)
        modifier_levels = levels(data_stratified[[modifier]])[c(2:length(levels(data_stratified[[modifier]])))]
      } else {
        RERI <- NA
        rel_int <- NA
        RERI_CI <- matrix(NA, nrow = 1, ncol = 2)
        rel_CI <- matrix(NA, nrow = 1, ncol = 2)
        modifier_levels = NA
      }
      
      outcome_remark = case_when(
        outcome == "haz_quarter" ~ "height-for-age z score",
        outcome == "whz_quarter" ~ "weight-for-height z score",
        outcome == "waz_quarter" ~ "weight-for-age z score",
        grepl("incident_haz_ms_stunt", outcome) ~ "incidence: moderate to severe stunting",
        grepl("incident_haz_s_stunt", outcome) ~ "incidence: severe stunting",
        grepl("incident_whz_ms_waste", outcome) ~ "incidence: moderate to severe wasting",
        grepl("incident_whz_s_waste", outcome) ~ "incidence: severe wasting",
        grepl("incident_waz_underwt", outcome) ~ "incidence: underweight",
        TRUE ~ NA
      )
      
      # Build output data frame
      output <- data.frame(
        outcome = outcome,
        outcome_remark = outcome_remark,
        time_unit = time_unit,
        age_group = age_group,
        gravidity_strata = strata_name,
        N_from_analysis = N,
        point_estimate = point_estimate_modified,
        SE = SE,
        lower_95CI = lower_CI,
        upper_95CI = upper_CI,
        modifier_name = ifelse(modifier=="NA", NA, modifier),
        modifier_level = modifier_levels,
        additive_effect = RERI,
        additive_CI = RERI_CI,
        relative_effect = rel_int,
        relative_CI = rel_CI
      )
      
      rownames(output) <- NULL
      output_list[[i]] = output
      
    }
    
    output_df = bind_rows(output_list)
    return(output_df)
  }


#--------------------------------------------------------
# apply the function and generate results
#--------------------------------------------------------

full_modifier_list = c("NA", "anymalaria", "maternal_agecat", "mother_heightcat", "wet_season")


crossing_zscore_3mo = crossing(
  outcome = outcome_zscore_quarter,
  time_unit = "3 month",
  age_group = age_list_3mo_birth,
  age_category = "agecat_birth"
) %>%
  mutate(age_group = as.character(age_group))

crossing_incidence_3mo = rbind(
  crossing(
    outcome = outcome_incidence_3mo,
    time_unit = "3 months",
    age_group = age_list_3mo_birth
  ),
  c("SGA_quarterly", "3 months", "Birth")
) %>% mutate(age_category = "agecat_birth", age_group = as.character(age_group))


# --------------------------------------------------------

aim1_application <-
  function(data_set,
           crossing_set,
           modifier_list,
           outcome_type) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      results_list = lapply(modifier_list, function(x)
        aim1_analysis(
          outcome = as.character(crossing_set[i, "outcome"]),
          data = data_set,
          time_unit = as.character(crossing_set[i, "time_unit"]),
          age_category = as.character(crossing_set[i, "age_category"]),
          age_group = as.character(crossing_set[i, "age_group"]),
          outcome_type = outcome_type,
          treatment = "Txarm",
          modifier = x
        ))
      full_results_list = append(full_results_list, results_list)
    }
    
    combined_results_df <- bind_rows(full_results_list)
    combined_results_df = data.frame(lapply(combined_results_df, function(x) if(is.numeric(x)) round(x, 4) else x))
    
    View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# apply the wrapper function and save the results
#--------------------------------------------------------

z_score_quarterly_results_stratified = aim1_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_3mo, modifier_list = full_modifier_list, outcome_type = "continuous")
#View(z_score_quarterly_results_stratified)
saveRDS(z_score_quarterly_results_stratified, paste0(results_path,"aim1-stratified/aim1_zscore_quarterly_results_stratified.RDS"))

incidence_3mo_results_stratified = aim1_application(data_set = data_incidence_3month, crossing_set = crossing_incidence_3mo, modifier_list = full_modifier_list, outcome_type = "binary")
#View(incidence_3mo_results_stratified)
saveRDS(incidence_3mo_results_stratified, paste0(results_path,"aim1-stratified/aim1_incidence_3mo_results_stratified.RDS"))

