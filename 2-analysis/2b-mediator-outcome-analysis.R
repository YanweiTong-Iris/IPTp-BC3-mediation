################################################################
# IPTp and child growth
# Script for intervention-mediator and mediator-outcome analysis
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------------------------
# load data
#--------------------------------------------------------

data_zscore_quarterly = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS"))
data_incidence_3month = readRDS(paste0(data_path,"analysis_data_incidence_quarterly.RDS"))
data_velocity_3month = readRDS(paste0(data_path,"analysis_data_velocity_3month.RDS")) 


#--------------------------------------------------------
# create outcome, age, and mediator list 
#--------------------------------------------------------

outcome_z_score_quarter = c("haz_quarter", "whz_quarter", "waz_quarter")
outcome_incidence_3mo = c(
  "incident_haz_ms_stunt_agecat_birth",
  "incident_haz_s_stunt_agecat_birth",
  "incident_whz_ms_waste_agecat_birth",
  "incident_whz_s_waste_agecat_birth",
  "incident_waz_underwt_agecat_birth"
)


age_list_3mo_birth = factor(c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"), levels= c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))


maternal_mediators_main = c("anemia_28binary", "gestational_weightchange",
                            "placentalmal", "preterm", "LBW",
                            "birthweight", "birthlength", "birthweight_kg",
                            # !NOTE: 4E_BP1 was renamed as F4E_BP1!
                            "IL8","VEGFA", "CD8A", "CDCP1", "CD244", 
                            "IL7", "OPG", "LAP_TGF_beta_1", "uPA",	"IL6",	
                            "IL_17C",	"MCP_1",	"CXCL11",	"AXIN1", "TRAIL",	
                            "CXCL9",	"CST5", "OSM", "CXCL1",	"CCL4",	"CD6", 
                            "SCF",	"IL18",	"TGF_alpha",	"MCP_4",	"CCL11",
                            "TNFSF14", "MMP_1", "LIF_R",	"FGF_21",	"CCL19", 
                            "IL_10RB", "IL_18R1","PD_L1",	"CXCL5",	"TRANCE", 
                            "HGF",	"IL_12B",	"MMP_10",	"IL10",	"TNF",
                            "CCL23",	"CD5",	"CCL3", "Flt3L", "CXCL6", 
                            "CXCL10", "F4E_BP1", "SIRT2",	"CCL28",	"DNER",	
                            "EN_RAGE",	"CD40",	"IFN_gamma", "FGF_19",	"MCP_2",	
                            "CASP_8",	"CCL25", "CX3CL1", "TNFRSF9",	"TWEAK",	"CCL20",	
                            "ST1A1",	"STAMBP",	"ADA",	"TNFB",	"CSF_1")

Olink_mediators = c("IL8","VEGFA", "CD8A", "CDCP1", "CD244", "IL7", "OPG", "LAP_TGF_beta_1",
                    "uPA",	"IL6",	"IL_17C",	"MCP_1",	"CXCL11",	"AXIN1", "TRAIL",	"CXCL9",	"CST5",
                    "OSM", "CXCL1",	"CCL4",	"CD6", "SCF",	"IL18",	"TGF_alpha",	"MCP_4",	"CCL11",
                    "TNFSF14", "MMP_1", "LIF_R",	"FGF_21",	"CCL19", "IL_10RB", "IL_18R1",
                    "PD_L1",	"CXCL5",	"TRANCE", "HGF",	"IL_12B",	"MMP_10",	"IL10",	"TNF",
                    "CCL23",	"CD5",	"CCL3", "Flt3L", "CXCL6", "CXCL10", "F4E_BP1", "SIRT2",	
                    "CCL28",	"DNER",	"EN_RAGE",	"CD40",	"IFN_gamma", "FGF_19",	"MCP_2",	
                    "CASP_8",	"CCL25", "CX3CL1", "TNFRSF9",	"TWEAK",	"CCL20",	"ST1A1",	
                    "STAMBP",	"ADA",	"TNFB",	"CSF_1")



#--------------------------------------------------------
# intervention-mediator & mediator-outcome analysis wrapper function
#--------------------------------------------------------

# Documentation: mediator_analysis
# Usage: mediator_analysis(data, time_unit, age_group, model, independent_var, dependent_var, dependent_type)
# Description: intervention-mediator & mediator-outcome analysis
# Intervention: Txarm OR rand_Txarm (for blind testing)
# Mediators of interest: low birth weight (LBW), preterm birth (preterm), 
#                        placental malaria (placentalmal), maternal anemia (anemia_28binary),
#                        gestational weight gain (gestational_weightchange), 
#                        gestational weight change Z score (GWC_Z),
#                        and maternal inflammation (Olink biomarkers). 
# For mediator-outcome models, use GEE with a Gaussian family, identity link, and exchangeable correlation matrix 
# for continuous outcomes since there were repeated measures within each age category. 
# For dichotomous outcomes, use a Poisson family and log link since there are no repeated measures. 
# Mediator-outcome models will be adjusted for confounders of mediators and outcomes by choosing nodes 
# sufficient to block backdoor pathways
# Note: when gravidity = all, the model is adjusted for gravidity; 
# when the data is stratified by multi or primigravida, then it's not adjusted for gravidity


# Args/Options:
# time_unit: "quarter"
# age_category: to indicate which age column to use (as a string)
#     - data_zscore_quarterly: agecat_birth
#     - data_incidence_3month: agecat_birth
# age_group: to subset the full dataset to a single age range
# model: "intervention-mediator" OR "mediator-outcome"
# independent_var: x of the formula (as a string)
# dependent_var: y of the formula (as a string)
# variable to choose: 
# intervention: "Txarm", "rand_Txarm"
# mediator: "anemia_28binary", "placentalmal", "preterm", "LBW",
#           "gestational_weightchange", "GWC_Z","birthweight", "birthlength", "birthweight_kg", 
#           "IL8","VEGFA", "CD8A", "CDCP1", "CD244", "IL7", "OPG", "LAP_TGF_beta_1",
#           "uPA",	"IL6",	"IL_17C",	"MCP_1",	"CXCL11",	"AXIN1", "TRAIL",	"CXCL9",	"CST5",
#           "OSM", "CXCL1",	"CCL4",	"CD6", "SCF",	"IL18",	"TGF_alpha",	"MCP_4",	"CCL11",
#           "TNFSF14", "MMP_1", "LIF_R",	"FGF_21",	"CCL19", "IL_10RB", "IL_18R1",
#           "PD_L1",	"CXCL5",	"TRANCE", "HGF",	"IL_12B",	"MMP_10",	"IL10",	"TNF",
#           "CCL23",	"CD5",	"CCL3", "Flt3L", "CXCL6", "CXCL10", "F4E_BP1", "SIRT2",	
#           "CCL28",	"DNER",	"EN_RAGE",	"CD40",	"IFN_gamma", "FGF_19",	"MCP_2",	
#           "CASP_8",	"CCL25", "CX3CL1", "TNFRSF9",	"TWEAK",	"CCL20",	"ST1A1",	
#           "STAMBP",	"ADA",	"TNFB",	"CSF_1"
# confounders: "sex", "gravidity_cat", "APdichenroll", "enrollage_binary", "educdich", "wealth_binary"
# outcome:
# - continuous: "haz", "whz", "waz", "wgv1", "wlz_gv1", "wgv3", "wlz_gv3",
#             "lgv2", "laz_gv2"， "lgv3", "laz_gv3"
# - binary: "SGA", haz_ms_stunt", "haz_s_stunt", "whz_ms_waste", "whz_s_waste", 
#         "waz_underweight", "incident_haz_ms_stunt", "incident_haz_s_stunt", 
#         "incident_whz_ms_waste", "incident_whz_s_waste", "incident_waz_underwt"
# dependent_type: "continuous" OR "binary" (as a string)

# Returns: the dataframe as described above
# Output: prints the data frame described above

mediator_analysis <- function(data, time_unit, age_category, age_group, model, independent_var, dependent_var, dependent_type, gravidity_strata) {
  print(paste0("model: ", model, ", x: ", independent_var, ", y: ", dependent_var, ", age: ", age_group))
  
  comment({#when mediator = anymalaria, for infant ever infected, only look at data post first infection
  if (independent_var == "new_anymalaria"){
    data = data %>% 
      group_by(id) %>%
      mutate(post_malaria = cumsum(as.numeric(as.character(new_anymalaria))) > 0 | max(as.numeric(as.character(new_anymalaria))) == 0) %>%
      ungroup()
    data = data[data$post_malaria,]
  }})
  
  data_age <- data[data[[age_category]] == age_group, ]
  data_age <- data_age[!is.na(data_age[[dependent_var]]), ]
  print(nrow(data_age))
  
  confounder <- NA
  adjusted <- 0
  
  if(gravidity_strata == "all"){
    data_stratified = data_age
  }else {
    data_stratified = data_age %>% filter(gravidity_cat == gravidity_strata)
    }
    

    # Compose the formula
       adjusted <- 1
       confounder <- "full"
      if (gravidity_strata == "all") {
          formula <-
            as.formula(paste(dependent_var, "~",independent_var,
                             "+ Txarm + sex + gravidity_cat+ GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      } else{
          formula <-as.formula(paste(dependent_var,"~",independent_var,
                                     "+ Txarm + sex + GAenroll + enrollage + APdichenroll + educdich + wealthcat"))
      }
    
    print(formula)
  
    # Check if the outcome has enough variation
    if (dependent_type == "binary") {
      counts = table(data_stratified[[dependent_var]])
      print(counts)
      if (!(length(unique(data_stratified[[dependent_var]])) >= 2 && all(counts >= 10))) {
        message("Variable '", dependent_var, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    }

    if(independent_var %in% c("anemia_28binary", "placentalmal",
                             "preterm", "LBW", "anymalaria")){
      counts = table(data_stratified[[independent_var]])
      print(counts)
      if (!(length(unique(data_stratified[[independent_var]])) >= 2 && all(counts >= 10))) {
        message("Variable '", independent_var, "' doesn't have enough variability. Skipping analysis.")
        return(NULL)
      }
    }
  
  
  # Run GLM model
  if (dependent_type == "continuous") {
    model.fit <- glm(formula, data = data_stratified, family = gaussian(link = "identity"), na.action = na.omit)
  } else if (dependent_type == "binary") {
    data_stratified <- data_stratified[!is.na(data_stratified[[independent_var]]), ]
    model.fit <- glm(formula, data = data_stratified, family = poisson(link = "log"), na.action = na.omit)
  }
  #print(summary(model.fit))
  
  N <- data_stratified %>% filter(!is.na(dependent_var)) %>% distinct(id) %>% count()
  N_tx <- data_stratified %>% 
    group_by(Txarm) %>% 
    summarise(N = n())
  
  N_DP <- N_tx$N[N_tx$Txarm == "DP"]
  N_SP <- N_tx$N[N_tx$Txarm == "SP"]
  
  estimates <- summary(model.fit)$coefficients
  coefs <- names(coef(model.fit))
  independent_var_row <- as.numeric(grep(paste0("^", independent_var), coefs))
  
  point_estimate <- estimates[independent_var_row, 1]
  point_estimate_modified <- ifelse(dependent_type == "continuous", point_estimate, exp(point_estimate))
  SE <- estimates[independent_var_row, 2]
  lower_CI <- ifelse(dependent_type == "continuous", point_estimate - qnorm(0.975) * SE, exp(point_estimate - qnorm(0.975) * SE))
  upper_CI <- ifelse(dependent_type == "continuous", point_estimate + qnorm(0.975) * SE, exp(point_estimate + qnorm(0.975) * SE))
  p_value <- estimates[independent_var_row, 4]
  
  outcome_remark <- case_when(
    dependent_var == "haz_quarter" | dependent_var == "haz" ~ "height-for-age z score",
    dependent_var == "whz_quarter" | dependent_var == "whz"~ "weight-for-height z score",
    dependent_var == "waz_quarter" | dependent_var == "waz" ~ "weight-for-age z score",
    dependent_var == "wgv1" ~ "absolute weight growth velocity (1 month increment)",
    dependent_var == "wgv3" ~ "absolute weight growth velocity (3 month increment)",
    dependent_var == "wlz_gv1" ~ "weight velocity based on weight-for-length Z-score (1 month increment)",
    dependent_var == "wlz_gv3" ~ "weight velocity based on weight-for-length Z-score (3 month increment)",
    dependent_var == "lgv2" ~ "absolute length growth velocity (2 month increment)",
    dependent_var == "lgv3" ~ "absolute length growth velocity (3 month increment)",
    dependent_var == "laz_gv2" ~ "length velocity based on length-for-age Z-score (2 month increment)",
    dependent_var == "laz_gv3" ~ "length velocity based on length-for-age Z-score (3 month increment)",
    grepl("SGA", dependent_var) ~ "small for gestational weight",
    dependent_var == "haz_ms_stunt_quarter" | dependent_var == "haz_ms_stunt" ~ "prevalence: moderate to severe stunting",
    dependent_var == "haz_s_stunt_quarter"| dependent_var == "haz_s_stunt" ~ "prevalence: severe stunting",
    dependent_var == "whz_ms_waste_quarter" | dependent_var == "whz_ms_waste" ~ "prevalence: moderate to severe wasting",
    dependent_var == "whz_s_waste_quarter" | dependent_var == "whz_s_waste" ~ "prevalence: severe wasting",
    dependent_var == "waz_underwt_quarter" | dependent_var == "waz_underwt" ~ "prevalence: underweight",
    grepl("incident_haz_ms_stunt", dependent_var) ~ "incidence: moderate to severe stunting",
    grepl("incident_haz_s_stunt", dependent_var) ~ "incidence: severe stunting",
    grepl("incident_whz_ms_waste", dependent_var) ~ "incidence: moderate to severe wasting",
    grepl("incident_whz_s_waste", dependent_var) ~ "incidence: severe wasting",
    grepl("incident_waz_underwt", dependent_var) ~ "incidence: underweight",
    TRUE ~ NA
  )
  
  mediator_remark <- case_when(
    independent_var == "anemia_28binary" | dependent_var == "anemia_28binary" ~ "Anemia at gestational week 28",
    independent_var == "placentalmal" | dependent_var == "placentalmal" ~ "Placental malaria",
    independent_var == "LBW" | dependent_var == "LBW" ~ "Low birth weight",
    independent_var == "preterm" | dependent_var == "preterm" ~ "Pre-term birth",
    independent_var == "gestational_weightchange" | dependent_var == "gestational_weightchange" ~ "Gestational weight change (kg)",
    independent_var == "SCF" | dependent_var == "betalactam_binary" ~ "Stem cell factor",
    independent_var == "birthweight" | dependent_var == "birthweight" ~ "Birth weight",
    independent_var == "birthweight_kg" | dependent_var == "birthweight_kg" ~ "Birth weight (kg)",
    independent_var == "birthlength" | dependent_var == "birthlength" ~ "Birth length",
    TRUE ~ NA
  )
  
  # Build output data frame
  if (model == "intervention-mediator") {
    output <- data.frame(
      model = model,
      gravidae = gravidity_strata,
      adjusted = adjusted,
      confounder = confounder,
      independent_variable = independent_var,
      dependent_variable = dependent_var,
      dependent_type = dependent_type,
      mediator_remark = mediator_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      N_DP = N_DP,
      N_SP = N_SP,
      point_estimate = point_estimate_modified,
      SE = SE,
      lower_95CI = lower_CI,
      upper_95CI = upper_CI,
      p_value = p_value
    )
  } else {
    output <- data.frame(
      model = model,
      gravidae = gravidity_strata,
      adjusted = adjusted,
      confounder = confounder,
      independent_variable = independent_var,
      dependent_variable = dependent_var,
      dependent_type = dependent_type,
      mediator_remark = mediator_remark,
      outcome_remark = outcome_remark,
      time_unit = time_unit,
      age_group = age_group,
      N_from_analysis = N,
      N_DP = N_DP,
      N_SP = N_SP,
      point_estimate = point_estimate_modified,
      SE = SE,
      lower_95CI = lower_CI,
      upper_95CI = upper_CI,
      p_value = p_value
    )
  }
  rownames(output) <- NULL
  return(output)
}


#--------------------------------------------------------
# function application
#--------------------------------------------------------

mediator_analysis_application <-
  function(data_set,
           crossing_set,
           model) {
    full_results_list = list()
    
    for (i in 1:nrow(crossing_set)) {
      try(full_results_list[[i]] <- mediator_analysis(
        data = data_set,
        time_unit = as.character(crossing_set[i, "time_unit"]),
        age_category = as.character(crossing_set[i, "age_category"]),
        age_group =  as.character(crossing_set[i, "age_group"]),
        gravidity_strata = as.character(crossing_set[i, "gravidity_strata"]),
        model = model,
        independent_var = as.character(crossing_set[i, "independent_var"]),
        dependent_var = as.character(crossing_set[i, "dependent_var"]),
        dependent_type = as.character(crossing_set[i, "dependent_type"])))
    }
    
    combined_results_df <- do.call(rbind, full_results_list)
    combined_results_df = data.frame(lapply(combined_results_df, function(x) if(is.numeric(x)) round(x, 4) else x))
    #View(combined_results_df)
    return(combined_results_df)
  }


#--------------------------------------------------------
# create crossing sets 
#--------------------------------------------------------

crossing_zscore_MO_3mo =crossing(
    independent_var = maternal_mediators_main,
    dependent_var = outcome_z_score_quarter,
    time_unit = "3 month",
    age_group = age_list_3mo_birth,
    age_category = "agecat_birth",
    gravidity_strata = c("all", "single", "multi"),
    dependent_type = "continuous"
  ) %>% mutate(age_group = as.character(age_group))


crossing_incidence_MO_3mo = rbind(
  crossing(independent_var = maternal_mediators_main,
           dependent_var = outcome_incidence_3mo,
           time_unit = "3 month",
           age_group = age_list_3mo_birth,
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary"),
  crossing(independent_var = maternal_mediators_main,
           dependent_var = "SGA_quarterly",
           time_unit = "3 month",
           age_group = "Birth",
           age_category = "agecat_birth",
           gravidity_strata = c("all", "single", "multi"),
           dependent_type = "binary")) %>%
  mutate(age_group = as.character(age_group))


#--------------------------------------------------------
# save the results with BH-adjusted p-value for Olink markers
#--------------------------------------------------------

#prevent using scientific notations
options(scipen = 999)

MO_zscore_3mo_stratified = mediator_analysis_application(data_set = data_zscore_quarterly, crossing_set = crossing_zscore_MO_3mo, model = "mediator-outcome")
MO_zscore_3mo_stratified = add_column(MO_zscore_3mo_stratified, adj_p_value = NA, .after="p_value" )
MO_zscore_3mo_stratified_non_Olink = MO_zscore_3mo_stratified %>%
  filter(!independent_variable %in% Olink_mediators)
MO_zscore_3mo_stratified_Olink = MO_zscore_3mo_stratified %>%
  filter(independent_variable %in% Olink_mediators) %>%
  group_by(gravidae, dependent_variable, age_group) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>% 
  ungroup()
MO_zscore_3mo_stratified_save = rbind(MO_zscore_3mo_stratified_non_Olink, MO_zscore_3mo_stratified_Olink)
sig_mediator_results = MO_zscore_3mo_stratified %>% 
  filter(!(independent_variable %in% Olink_mediators) | (independent_variable %in% Olink_mediators & p_value <= 0.05)) 
sig_mediators = unique(sig_mediator_results$independent_variable)
MO_zscore_3mo_stratified_save = MO_zscore_3mo_stratified_save %>% filter(independent_variable %in% sig_mediators)
#View(MO_zscore_3mo_stratified_save)
saveRDS(MO_zscore_3mo_stratified_save, paste0(results_path,"IM-MO-stratified/mediator_outcome_zscore_results_3mo_stratified.RDS"))


MO_incidence_3mo_stratified = mediator_analysis_application(data_set = data_incidence_3month, crossing_set = crossing_incidence_MO_3mo, model = "mediator-outcome")
#View(MO_incidence_3mo_stratified)
MO_incidence_3mo_stratified = add_column(MO_incidence_3mo_stratified, adj_p_value = NA, .after="p_value" )
MO_incidence_3mo_stratified_non_Olink = MO_incidence_3mo_stratified %>%
  filter(!independent_variable %in% Olink_mediators)
MO_incidence_3mo_stratified_Olink = MO_incidence_3mo_stratified %>%
  filter(independent_variable %in% Olink_mediators) %>%
  group_by(gravidae, dependent_variable, age_group) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>% 
  ungroup()
MO_incidence_3mo_stratified_save = rbind(MO_incidence_3mo_stratified_non_Olink, MO_incidence_3mo_stratified_Olink)
sig_mediator_results = MO_incidence_3mo_stratified %>% 
  filter(!(independent_variable %in% Olink_mediators) | (independent_variable %in% Olink_mediators & p_value <= 0.05)) 
sig_mediators = unique(sig_mediator_results$independent_variable)
MO_incidence_3mo_stratified_save = MO_incidence_3mo_stratified_save %>% filter(independent_variable %in% sig_mediators)
#View(MO_incidence_3mo_stratified_save)
saveRDS(MO_incidence_3mo_stratified_save, paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_3mo_stratified.RDS"))

