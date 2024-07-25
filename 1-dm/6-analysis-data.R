################################################################
# IPTp and child growth
# Script for preparing the data sets for aim1 & 2 analysis
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

# -----------------------------------------
# load data
# -----------------------------------------

dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))
df_season = readRDS(paste0(box_path, "/processed_data/outcome_zscores_covariates.RDS"))[, c("uniqueid", "wet_season", "incidentmalaria_season")]
dfz = merge(dfz, df_season, by = "uniqueid")
df_wgv1 = readRDS(paste0(box_path,"/processed_data/weight_linear_growth_velocity_1month.RDS"))[, c("id", "month1_interval","age_diff1","wgv1","wlz_gv1")]
df_lgv2 = readRDS(paste0(box_path,"/processed_data/length_linear_growth_velocity_2month.RDS"))[, c("id", "month2_interval", "age_diff2","lgv2", "laz_gv2")]
df_gv3 = readRDS(paste0(box_path,"/processed_data/weight&length_linear_growth_velocity_3month.RDS"))[, c("id", "month3_interval", "age_diff3", "wgv3", "wlz_gv3", "lgv3", "laz_gv3")]
df_incidence_monthly_original = readRDS(paste0(box_path, "/processed_data/incidence_data_monthly.RDS"))
df_incidence_quarterly_original = readRDS(paste0(box_path, "/processed_data/incidence_data_quarterly.RDS"))
df_incidence_biannual_original = readRDS(paste0(box_path, "/processed_data/incidence_data_biannual.RDS"))
df_incidence_annual_original = readRDS(paste0(box_path, "/processed_data/incidence_data_annual.RDS"))
df_maternal = read_dta(paste0(box_path, "Maternal/", "BC-3 mothers expanded database FINAL.dta"))
df_maternal_indi = read_dta(paste0(box_path, "Maternal/", "Maternal BC3 individual level database FINAL.dta"))


Olink_broad = read_excel(paste0(box_path, "Samples/","05-18-23 Olink_Inflammation_BC3 Maternal Speciemen_3 plates_Report_broad.xlsx"))
Olink_broad = Olink_broad %>% 
  mutate(motherid = substr(ID, start = 1, stop = 5)) %>%
  dplyr::select(motherid, "IL8","VEGFA", "CD8A", "CDCP1", "CD244", "IL7", "OPG", "LAP_TGF_beta_1",
                "uPA",	"IL6",	"IL_17C",	"MCP_1",	"CXCL11",	"AXIN1", "TRAIL",	"CXCL9",	"CST5",
                "OSM", "CXCL1",	"CCL4",	"CD6", "SCF",	"IL18",	"TGF_alpha",	"MCP_4",	"CCL11",
                "TNFSF14", "MMP_1", "LIF_R",	"FGF_21",	"CCL19", "IL_10RB", "IL_18R1",
                "PD_L1",	"CXCL5",	"TRANCE", "HGF",	"IL_12B",	"MMP_10",	"IL10",	"TNF",
                "CCL23",	"CD5",	"CCL3", "Flt3L", "CXCL6", "CXCL10", "4E_BP1", "SIRT2",	
                "CCL28",	"DNER",	"EN_RAGE",	"CD40",	"IFN_gamma", "FGF_19",	"MCP_2",	
                "CASP_8",	"CCL25", "CX3CL1", "TNFRSF9",	"TWEAK",	"CCL20",	"ST1A1",	
                "STAMBP",	"ADA",	"TNFB",	"CSF_1") %>%
  mutate(motherid = as.numeric(motherid)) 
#modify Olink marker name (cannot take numerical values as initials)
colnames(Olink_broad)[which(names(Olink_broad) == "4E_BP1")] <- "F4E_BP1"


main_covariate_list = c("sex", "enrollage", "maternal_agecat", "GAenroll", "Gravidity", "gravidity_cat", "LBW", "anyHP", "preterm", "birthweight", "birthweight_kg", 
                        "birthlength", "birthlength_m","mother_height", "mother_heightcat", "gestational_weightchange", "GWC_Z","gestational_weightchange_cat", "wet_season", 
                        "incidentmalaria", "anymalaria","incidentmalaria_season", "anyhb_28", "anemia_28", "anemia_28binary", "anyhb_36", "anemia_36" , 
                        "anemia_36binary", "wbc_t1", "wbc_t2", "wbc_t3", "wbc_t3.t2","neutro_t1", "neutro_t2", "neutro_t3", "neutro_t3.t2",
                        "placentalmal", "placentalBSdich", "placentalLAMPdich", "BSdichenroll","Graviddich", "APdichenroll", "educdich", 
                        "enrollage_binary", "enrollage_binary_30", "wealth_binary", "wealthcat", "IL8","VEGFA", "CD8A", "CDCP1", "CD244", "IL7", 
                        "OPG", "LAP_TGF_beta_1", "uPA",	"IL6",	"IL_17C",	"MCP_1",	"CXCL11",	"AXIN1", "TRAIL",	"CXCL9",	"CST5",
                        "OSM", "CXCL1",	"CCL4",	"CD6", "SCF",	"IL18",	"TGF_alpha",	"MCP_4",	"CCL11", "TNFSF14", "MMP_1", "LIF_R",	"FGF_21",	
                        "CCL19", "IL_10RB", "IL_18R1", "PD_L1",	"CXCL5",	"TRANCE", "HGF",	"IL_12B",	"MMP_10",	"IL10",	"TNF", "CCL23",
                        "CD5",	"CCL3", "Flt3L", "CXCL6", "CXCL10", "F4E_BP1", "SIRT2", "CCL28",	"DNER",	"EN_RAGE",	"CD40",	"IFN_gamma", 
                        "FGF_19",	"MCP_2", "CASP_8",	"CCL25", "CX3CL1", "TNFRSF9",	"TWEAK",	"CCL20",	"ST1A1", "STAMBP",	"ADA",	"TNFB",	"CSF_1") 



# -----------------------------------------
# extract useful columns and create new modified covariates
# -----------------------------------------

ages <- dfz %>% dplyr::select(id, uniqueid, age, agecat, agedays, agecat_birth, agemonthcat) %>%
  mutate(uniqueid = as.factor(uniqueid)) %>%
  
  mutate(agecat = factor(agecat, levels = c("0-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>% 
  
  mutate(age_2_month = case_when(
    agedays <60.875 ~ "<2 months",
    agedays >= 60.875 & agedays <=121.75 ~ "2- 4 months",
    agedays >= 121.75 & agedays <=182.625 ~ "4- 6 months",
    agedays >= 182.625 & agedays <=243.5 ~ "6- 8 months",
    agedays >= 243.5 & agedays <=304.375 ~ "8- 10 months",
    agedays >= 304.375 & agedays <=365.25 ~ "10- 12 months")) %>%
  
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

# weight change = weight of week 36- weight of week 20
# GWC-z at wk 36 based on https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(22)00055-X/fulltext
# excel formula =(1.382972-56.14743*POWER(gestational week,-2)+0.2787683*POWER(gestational week,0.5))
GWC_mean = (1.382972-56.14743*(36^(-2)))+0.2787683*(36^0.5)
# excel formula = 0.2501993731+142.4297879*POWER(gestational week,-2)-61.45345*POWER(gestational week,-2)*LN(gestational week)
GWC_SD = 0.2501993731+142.4297879*(36^(-2))-61.45345*(36^(-2))*log(36)
# height & delta weight cutoffs at 1st, 2nd, and 3rd quantile

maternal <- df_maternal %>%
  dplyr::select(id, gestwks, weight, height) %>%
  fill(height, .direction = "down") %>%
  filter(gestwks == 20 | gestwks == 36) %>%
  mutate(mother_height = height) %>%
  mutate(mother_heightcat = paste0("Q", ntile(mother_height, 4))) %>%
  mutate(mother_heightcat = factor(mother_heightcat, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  group_by(id) %>%
  mutate(gestational_weightchange = as.numeric(ifelse(length(gestwks) == 2, weight[2] - weight[1], NA))) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(gestational_weightchange_cat = paste0("Q", cut(gestational_weightchange, breaks = quantile(gestational_weightchange, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), labels = FALSE, include.lowest = TRUE))) %>% 
  mutate(gestational_weightchange_cat = factor(gestational_weightchange_cat, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  mutate(GWC_tranformation = log(gestational_weightchange+3.3+8.75)) %>%
  mutate(GWC_Z = (GWC_tranformation-GWC_mean)/GWC_SD) %>%
  mutate(motherid = id) %>%
  dplyr::select(motherid, mother_height, mother_heightcat, gestational_weightchange, GWC_Z, gestational_weightchange_cat)


# generate wbc and neutrophill for each trimester
maternal_blood_count <- df_maternal %>% 
  dplyr::select(id, gestwks, wbc, neutro) %>%
  na.omit() %>%
  group_by(id) %>%
  summarise(
    wbc_t1 = if(length(wbc[gestwks >= 0 & gestwks <= 12]) == 0) NA else mean(wbc[gestwks >= 0 & gestwks <= 12], na.rm = TRUE),
    wbc_t2 = if(length(wbc[gestwks >= 13 & gestwks <= 27]) == 0) NA else mean(wbc[gestwks >= 13 & gestwks <= 27], na.rm = TRUE),
    wbc_t3 = if(length(wbc[gestwks >= 28]) == 0) NA else mean(wbc[gestwks >= 28], na.rm = TRUE),
    neutro_t1 = if(length(neutro[gestwks >= 0 & gestwks <= 12]) == 0) NA else mean(neutro[gestwks >= 0 & gestwks <= 12], na.rm = TRUE),
    neutro_t2 = if(length(neutro[gestwks >= 13 & gestwks <= 27]) == 0) NA else mean(neutro[gestwks >= 13 & gestwks <= 27], na.rm = TRUE),
    neutro_t3 = if(length(neutro[gestwks >= 28]) == 0) NA else mean(neutro[gestwks >= 28], na.rm = TRUE)
  ) %>%
  mutate(wbc_t3.t2 = wbc_t3/wbc_t2, neutro_t3.t2 = neutro_t3/neutro_t2) %>%
  ungroup()


# generate anyhb for week 28 and week 36
maternal_hb_28 <- df_maternal %>%
  dplyr::select(id, gestwks, anyhb) %>%
  filter(gestwks == 28) %>%
  mutate(anyhb_28 = anyhb) %>%
  mutate(anemia_28 = case_when(
    anyhb_28 >= 12 ~ NA,
    anyhb_28 >= 11 & anyhb_28 < 12 ~ "mild",
    anyhb_28 >= 8  & anyhb_28 < 11 ~ "moderate",
    anyhb_28 < 8 ~ "severe")) %>%
  mutate(anemia_28 = factor(anemia_28, levels = c("mild", "moderate", "severe"))) %>%
  mutate(anemia_28binary = factor(ifelse(anyhb_28 <11, 1, 0), levels = c(0,1))) %>%
  dplyr::select(id, anyhb_28, anemia_28, anemia_28binary)

  
maternal_hb_36 <- df_maternal %>%
  dplyr::select(id, gestwks, anyhb) %>%
  filter(gestwks == 36) %>%
  mutate(anyhb_36 = anyhb) %>%
  mutate(anemia_36 = case_when(
    anyhb_36 >= 12 ~ NA,
    anyhb_36 >= 11 & anyhb_36 < 12 ~ "mild",
    anyhb_36 >= 8  & anyhb_36 < 11 ~ "moderate",
    anyhb_36 < 8 ~ "severe")) %>%
  mutate(anemia_36 = factor(anemia_36, levels = c("mild", "moderate", "severe"))) %>%
  mutate(anemia_36binary = factor(ifelse(anyhb_36 <11, 1, 0), levels = c(0,1))) %>%
  dplyr::select(id, anyhb_36, anemia_36, anemia_36binary)

maternal_main = merge(maternal, maternal_blood_count, by.x= "motherid", by.y = "id", all.x = TRUE)
maternal_main = merge(maternal_main, maternal_hb_28, by.x= "motherid", by.y = "id", all.x = TRUE)
maternal_main = merge(maternal_main, maternal_hb_36, by.x= "motherid", by.y = "id", all.x = TRUE)



# inflammation
maternal_main = maternal_main %>%
  left_join(Olink_broad, by = c("motherid" = "motherid"))


# check maternal merged
#View(maternal_main)
colnames(maternal_main)


# extract all main covariates
covariates <-
  dfz %>% dplyr::select(uniqueid,id, motherid, Txarm, rand_Txarm, sex, enrollage, 
                        maternal_agecat, Gravidity, preterm, birthweight, height,
                        wet_season, incidentmalaria_season, LBW, anyHP, placentalmal, 
                        placentalBSdich, placentalLAMPdich, BSdichenroll, Graviddich, 
                        APdichenroll, GAenroll, educdich, enrollage_binary, enrollage_binary_30, 
                        wealth_binary, wealthcat, incidentmalaria, anymalaria) %>%
  mutate(uniqueid = factor(uniqueid)) %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP"))) %>%
  mutate(sex = factor(sex, levels = c("F", "M"))) %>%
  mutate(maternal_agecat = factor(
    maternal_agecat,
    levels = c("less than 20", "20-24", "25-29", "30+"))) %>%
  mutate(gravidity_cat = ifelse(Gravidity == 1, "single", "multi")) %>%
  mutate(gravidity_cat = factor(gravidity_cat, levels = c("single", "multi"))) %>%
  mutate(preterm = factor(preterm, levels = c("0", "1"))) %>%
  mutate(
    incidentmalaria = as.factor(incidentmalaria),
    anymalaria = as.factor(anymalaria),
    LBW = as.factor(LBW),
    wet_season = as.factor(wet_season),
    incidentmalaria_season = as.factor(incidentmalaria_season),
    birthweight_kg = birthweight/1000) %>%
  mutate(anyHP = factor(ifelse(anyHP == "Yes", 1, 0))) %>%
  group_by(id) %>%
  mutate(birthlength = first(height)) %>%
  ungroup() %>%
  dplyr::select(-height)

main_covariates = merge(covariates, maternal_main, by = "motherid", all.x = TRUE)

id_list = unique(ages$id)
data = merge(ages, main_covariates, by = c("uniqueid","id")) %>% mutate(uniqueid = as.factor(uniqueid))
colnames(data)[colnames(data)=="agemonthcat"] <- "agemonth_ceiling"

# test a post infant malaria infection flag for mediator-outcome analysis
data = data %>% 
  group_by(id) %>%
  mutate(post_malaria = cumsum(as.numeric(as.character(incidentmalaria))) > 0 | max(as.numeric(as.character(incidentmalaria))) == 0,
         birthlength_m = birthlength/100) %>%
  ungroup()

#View(data)

# -----------------------------------------
# filter incidence data and add SGA
# -----------------------------------------

indi_SGA = dplyr::select(dfz, c("id", "SGA")) %>% group_by(id) %>% slice(1)

df_incidence_monthly = merge(df_incidence_monthly_original, indi_SGA, by = "id") %>% group_by(id) %>%
  mutate(incident_haz_ms_stunt_agemonthcat = ifelse(!atrisk_haz_ms_stunt_agemonthcat, NA, incident_haz_ms_stunt_agemonthcat)) %>%
  mutate(incident_haz_s_stunt_agemonthcat = ifelse(!atrisk_haz_s_stunt_agemonthcat, NA, incident_haz_s_stunt_agemonthcat)) %>%
  mutate(incident_whz_ms_waste_agemonthcat = ifelse(!atrisk_whz_ms_waste_agemonthcat, NA, incident_haz_ms_stunt_agemonthcat)) %>%
  mutate(incident_whz_s_waste_agemonthcat = ifelse(!atrisk_whz_s_waste_agemonthcat, NA, incident_whz_s_waste_agemonthcat)) %>%
  mutate(incident_waz_underwt_agemonthcat = ifelse(!atrisk_waz_underwt_agemonthcat, NA, incident_waz_underwt_agemonthcat)) %>%
  mutate(SGA_monthly = ifelse(row_number() > 1, NA, SGA)) %>% 
  ungroup() %>%
  dplyr::select(- SGA)
#View(df_incidence_monthly)

df_incidence_quarterly = merge(df_incidence_quarterly_original, indi_SGA, by = "id")%>% group_by(id) %>%
  mutate(incident_haz_ms_stunt_agecat_birth = ifelse(!atrisk_haz_ms_stunt_agecat_birth, NA, incident_haz_ms_stunt_agecat_birth)) %>%
  mutate(incident_haz_s_stunt_agecat_birth = ifelse(!atrisk_haz_s_stunt_agecat_birth, NA, incident_haz_s_stunt_agecat_birth)) %>%
  mutate(incident_whz_ms_waste_agecat_birth = ifelse(!atrisk_whz_ms_waste_agecat_birth, NA, incident_whz_ms_waste_agecat_birth)) %>%
  mutate(incident_whz_s_waste_agecat_birth = ifelse(!atrisk_whz_s_waste_agecat_birth, NA, incident_whz_s_waste_agecat_birth)) %>%
  mutate(incident_waz_underwt_agecat_birth = ifelse(!atrisk_waz_underwt_agecat_birth, NA, incident_waz_underwt_agecat_birth)) %>%
  mutate(SGA_quarterly = ifelse(row_number() > 1, NA, SGA)) %>% 
  ungroup() %>%
  dplyr::select(- SGA)
#View(df_incidence_quarterly)

df_incidence_biannual = merge(df_incidence_biannual_original, indi_SGA, by = "id") %>% group_by(id) %>%
  mutate(incident_haz_ms_stunt_age_6_12_month = ifelse(!atrisk_haz_ms_stunt_age_6_12_month, NA, incident_haz_ms_stunt_age_6_12_month)) %>%
  mutate(incident_haz_s_stunt_age_6_12_month = ifelse(!atrisk_haz_s_stunt_age_6_12_month, NA, incident_haz_s_stunt_age_6_12_month)) %>%
  mutate(incident_whz_ms_waste_age_6_12_month = ifelse(!atrisk_whz_ms_waste_age_6_12_month, NA, incident_whz_ms_waste_age_6_12_month)) %>%
  mutate(incident_whz_s_waste_age_6_12_month = ifelse(!atrisk_whz_s_waste_age_6_12_month, NA, incident_whz_s_waste_age_6_12_month)) %>%
  mutate(incident_waz_underwt_age_6_12_month = ifelse(!atrisk_waz_underwt_age_6_12_month, NA, incident_waz_underwt_age_6_12_month)) %>%
  mutate(SGA_biannual = ifelse(row_number() > 1, NA, SGA)) %>% 
  ungroup() %>%
  dplyr::select(- SGA)
#View(df_incidence_biannual)

df_incidence_annual = merge(df_incidence_annual_original, indi_SGA, by = "id") %>% group_by(id) %>%
  mutate(incident_haz_ms_stunt_age_1_12 = ifelse(!atrisk_haz_ms_stunt_age_1_12, NA, incident_haz_ms_stunt_age_1_12)) %>%
  mutate(incident_haz_s_stunt_age_1_12 = ifelse(!atrisk_haz_s_stunt_age_1_12, NA, incident_haz_s_stunt_age_1_12)) %>%
  mutate(incident_whz_ms_waste_age_1_12 = ifelse(!atrisk_whz_ms_waste_age_1_12, NA, incident_whz_ms_waste_age_1_12)) %>%
  mutate(incident_whz_s_waste_age_1_12 = ifelse(!atrisk_whz_s_waste_age_1_12, NA, incident_whz_s_waste_age_1_12)) %>%
  mutate(incident_waz_underwt_age_1_12 = ifelse(!atrisk_waz_underwt_age_1_12, NA, incident_waz_underwt_age_1_12)) %>%
  mutate(SGA_annual = ifelse(row_number() > 1, NA, SGA)) %>% 
  ungroup() %>%
  dplyr::select(- SGA)
#View(df_incidence_annual)


# -----------------------------------------
# continuous
# -----------------------------------------

continuous_prevalence_outcomes <-
  dfz %>% dplyr::select(uniqueid, id, age, agedays, agecat_birth, haz, whz, waz, haz_ms_stunt, haz_s_stunt, whz_ms_waste, whz_s_waste, waz_underwt)

continuous_prevalence_outcomes_mutated <-
  continuous_prevalence_outcomes %>% mutate(agemonth_round = round(age)) %>%
  mutate(agemonth_ceiling = ceiling(age)) %>%
  mutate(uniqueid= as.factor(uniqueid))

#View(continuous_prevalence_outcomes_mutated)

data_continuous = merge(
  dplyr::select(continuous_prevalence_outcomes_mutated, c("uniqueid", "id", "age", "agecat_birth","haz", "whz", "waz", "agemonth_round", "agemonth_ceiling")),
  dplyr::select(data, c("uniqueid", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
  by = c("uniqueid")
) %>% mutate(agemonth_round = as.factor(agemonth_round), agemonth_ceiling= as.factor(agemonth_ceiling)) 

#View(data_continuous)


data_full_longitudinal = merge(
  dplyr::select(continuous_prevalence_outcomes_mutated, 
                c("uniqueid", "id", "age", "agemonth_round", "agecat_birth",
                  "haz", "whz", "waz", "haz_ms_stunt", "haz_s_stunt", 
                  "whz_ms_waste", "whz_s_waste", "waz_underwt")),
  dplyr::select(data, c("uniqueid", "motherid", "Txarm", "rand_Txarm", "agedays","agemonth_ceiling",main_covariate_list)),
  by = c("uniqueid")
) %>% mutate(agemonth_round = as.factor(agemonth_round))

colnames(df_incidence_monthly)[2] = "agemonth_ceiling"

data_full_longitudinal = merge(data_full_longitudinal,
                               df_incidence_monthly,
                               by = c("id", "agemonth_ceiling")) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

data_full_longitudinal <- data_full_longitudinal %>%
  arrange(id, agedays) %>%  
  group_by(id) %>%
  mutate(
    lagged_anymalaria = ifelse(
      any(anymalaria == 1) & (agedays <= agedays[which.max(anymalaria)] + 30),
      1, 
      0
    )
  ) %>%
  ungroup() %>% 
  group_by(id, agemonth_ceiling) %>%
  mutate(malaria_numeric = as.numeric(as.character(anymalaria))) %>%
  mutate(sum_anymalaria = sum(malaria_numeric)) %>%
  ungroup()


# -----------------------------------------
# zscore and prevalence (monthly cutoff)
# -----------------------------------------

data_monthly_2 = data %>% mutate(agemonth_round = round(age)) %>% 
  group_by(id, agemonth_round) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  ungroup() 

prevalence_outcomes_mutated <-
  continuous_prevalence_outcomes_mutated %>% 
  group_by(id, agemonth_round) %>%
  mutate(haz = mean(haz, na.rm = TRUE), 
         waz = mean(waz, na.rm = TRUE), 
         whz = mean(whz, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(haz_ms_stunt = ifelse(haz< -2, 1, 0),
         haz_s_stunt = ifelse(haz< -3, 1, 0),
         waz_underwt = ifelse(waz< -2, 1, 0),
         whz_ms_waste = ifelse(whz< -2, 1, 0),
         whz_s_waste = ifelse(whz< -3, 1, 0))

data_zscore_monthly = merge(
  dplyr::select(prevalence_outcomes_mutated, c("uniqueid", "id", "age", "agemonth_round", "haz" , "whz", "waz")),
  dplyr::select(data_monthly_2, c("id","agemonth_round", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
  by = c("id", "agemonth_round")
) %>% mutate(agemonth_round = as.factor(agemonth_round)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()
#View(data_zscore_monthly)

data_prevalence_monthly = merge(
  dplyr::select(prevalence_outcomes_mutated, c("uniqueid", "id", "age", "agemonth_round", "haz_ms_stunt" ,"haz_s_stunt", "whz_ms_waste", "whz_s_waste", "waz_underwt")),
  dplyr::select(data_monthly_2, c("id","agemonth_round", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
  by = c("id", "agemonth_round")
) %>% mutate(agemonth_round = as.factor(agemonth_round)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

#View(data_prevalence_monthly)


# -----------------------------------------
# velocity, 1 month
# -----------------------------------------

data_monthly = data %>%
  group_by(id, agemonth_ceiling) %>%
  mutate(malaria_numeric = as.numeric(as.character(anymalaria))) %>%
  mutate(sum_anymalaria = sum(malaria_numeric)) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP")))

df_wgv1_new <- df_wgv1 %>%
  separate(
    month1_interval,
    into = c("age_from", "age_to"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(age_from = as.numeric(sub("mo", "", age_from)),
         age_to = as.numeric(sub("mo", "", age_to)))

data_velocity_1month = merge(
  dplyr::select(df_wgv1_new, c("id", "month1_interval", "age_to", "age_diff1", "wgv1", "wlz_gv1")),
  dplyr::select(data_monthly, c("uniqueid", "id", "age", "agemonth_ceiling", "motherid", "Txarm", "rand_Txarm", main_covariate_list, "sum_anymalaria")),
  by.y = c("id", "agemonth_ceiling"),
  by.x = c("id", "age_to"),
) %>% group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

colnames(data_velocity_1month)[2] <- "agemonth_ceiling"
#View(data_velocity_1month)


# -----------------------------------------
# incidence, monthly
# -----------------------------------------

colnames(df_incidence_monthly)[2] = "agemonth_ceiling"

data_incidence_1month = merge(data_incidence_1month,
                              dplyr::select(data_zscore_monthly, c("id", "agemonth_round", "haz", "whz", "waz")),
                              by.x = c("id", "agemonth_ceiling"),
                              by.y = c("id", "agemonth_round")) %>%
  mutate(incident_haz_ms_stunt_agemonthcat = ifelse(is.na(haz), NA, incident_haz_ms_stunt_agemonthcat),
         atrisk_haz_ms_stunt_agemonthcat = ifelse(is.na(haz), NA, atrisk_haz_ms_stunt_agemonthcat),
         incident_haz_s_stunt_agemonthcat = ifelse(is.na(haz), NA, incident_haz_s_stunt_agemonthcat),
         atrisk_haz_s_stunt_agemonthcat = ifelse(is.na(haz), NA, atrisk_haz_s_stunt_agemonthcat),
         incident_whz_ms_waste_agemonthcat = ifelse(is.na(whz), NA, incident_whz_ms_waste_agemonthcat),
         atrisk_whz_ms_waste_agemonthcat = ifelse(is.na(whz), NA, atrisk_whz_ms_waste_agemonthcat),
         incident_whz_s_waste_agemonthcat = ifelse(is.na(whz), NA, incident_whz_s_waste_agemonthcat),
         atrisk_whz_s_waste_agemonthcath = ifelse(is.na(whz), NA, atrisk_whz_s_waste_agemonthcat),
         incident_waz_underwt_agemonthcat = ifelse(is.na(waz), NA, incident_waz_underwt_agemonthcat),
         atrisk_waz_underwt_agemonthcat = ifelse(is.na(waz), NA, atrisk_waz_underwt_agemonthcat)) %>%
  mutate(agemonth_ceiling = factor(agemonth_ceiling))

#View(data_incidence_1month)


# -----------------------------------------
# monthly full outcomes
# (prepare for infant malaria analysis)
# -----------------------------------------

data_full_monthly_round = merge(
  dplyr::select(
    prevalence_outcomes_mutated,
    c("uniqueid", "id", "age",
      "agemonth_round", "haz", "whz", "waz",
      "haz_ms_stunt", "haz_s_stunt", "whz_ms_waste",
      "whz_s_waste", "waz_underwt"
    )),
  dplyr::select(
    data_monthly_2,
    c("id","agemonth_round", "motherid", "Txarm", "rand_Txarm",
      main_covariate_list)),
  by = c("id", "agemonth_round")
) %>% mutate(agemonth_round = as.factor(agemonth_round)) %>%
  group_by(id) %>%
  arrange(id, agemonth_round) %>%
  mutate(malaria_numeric = as.numeric(as.character(anymalaria))) %>%
  mutate(new_anymalaria = ifelse(
      lag(malaria_numeric, default = 0) == 1,
      1, malaria_numeric)) %>%
  ungroup() %>% 
  mutate_if(is.numeric, ~ifelse(is.nan(.), NA, .))

#View(data_full_monthly_round)


zscore_monthly_ceiling <-
  continuous_prevalence_outcomes_mutated %>% 
  group_by(id, agemonth_ceiling) %>%
  mutate(haz = mean(haz, na.rm = TRUE), 
         waz = mean(waz, na.rm = TRUE), 
         whz = mean(whz, na.rm = TRUE)) %>%
  slice(1) %>%
  mutate(haz_ms_stunt = ifelse(haz< -2, 1, 0),
         haz_s_stunt = ifelse(haz< -3, 1, 0),
         waz_underwt = ifelse(waz< -2, 1, 0),
         whz_ms_waste = ifelse(whz< -2, 1, 0),
         whz_s_waste = ifelse(whz< -3, 1, 0)) %>%
  ungroup() 


data_full_monthly_ceiling = merge(
  zscore_monthly_ceiling,
  data_velocity_1month,
  by = c("id", "agemonth_ceiling"),
) %>% group_by(id) %>%
  arrange(id, age.x) %>%
  mutate(malaria_numeric = as.numeric(as.character(anymalaria))) %>%
  mutate(new_anymalaria = ifelse(
    lag(malaria_numeric, default = 0) == 1,
    1, malaria_numeric)) %>%
  ungroup() 

#View(data_full_monthly_ceiling)


# -----------------------------------------
# velocity, 2 month
# -----------------------------------------

data_2month = data %>% group_by(id, age_2_month) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  arrange(id, age)%>%
  ungroup()
#View(data_2month)

df_lgv2_new <- df_lgv2 %>%
  separate(
    month2_interval,
    into = c("age_from", "age_to"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(age_from = as.numeric(sub("mo", "", age_from)),
         age_to = as.numeric(sub("mo", "", age_to))) %>%
  mutate(
    age_2_month = case_when(
      age_to == 2 ~ "<2 months",
      age_to == 4 ~ "2- 4 months",
      age_to == 6 ~ "4- 6 months",
      age_to == 8 ~ "6- 8 months",
      age_to == 10 ~ "8- 10 months",
      age_to == 12 ~ "10- 12 months"
    )
  )
head(df_lgv2_new)

data_velocity_2month = merge(dplyr::select(df_lgv2_new, c("id", "month2_interval", "age_diff2", "lgv2", "laz_gv2", "age_2_month")),
                             dplyr::select(data_2month, c("uniqueid", "id", "age", "age_2_month", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
                             by = c("id", "age_2_month"),) %>%
  mutate(age_2_month = as.factor(age_2_month)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

#View(data_velocity_2month)


# -----------------------------------------
# velocity, 3 month
# -----------------------------------------

data_quarter = data %>% group_by(id, agecat) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  arrange(id, age)%>%
  ungroup()

df_gv3_new <- df_gv3 %>%
  separate(
    month3_interval,
    into = c("age_from", "age_to"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(age_from = as.numeric(sub("mo", "", age_from)),
         age_to = as.numeric(sub("mo", "", age_to))) %>%
  mutate(
    agecat = case_when(
      age_to == 3 ~ "0-3 months",
      age_to == 6 ~ ">3-6 months",
      age_to == 9 ~ ">6-9 months",
      age_to == 12 ~ ">9-12 months",
    )
  )

data_velocity_3month = merge(dplyr::select(df_gv3_new, c("id", "agecat", "month3_interval", "age_diff3", "wgv3", "wlz_gv3", "lgv3", "laz_gv3")),
                             dplyr::select(data_quarter, c("uniqueid", "id", "age", "agecat", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
                             by = c("id", "agecat"),) %>%
  mutate(agecat = as.factor(agecat)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

#View(data_velocity_3month)


# -----------------------------------------
# # zscore and prevalence, quarterly
# -----------------------------------------

preceding_malaria <- dfz %>% 
  dplyr::select(id, age, agecat_birth, haz_ms_stunt, whz_ms_waste, anymalaria) %>%
  group_by(id, agecat_birth) %>%
  mutate(
    previous_malaria = cummax(lag(anymalaria, default = 0)),
    prec_malaria_stunt = ifelse(previous_malaria == 1 & haz_ms_stunt == 1, 1, 0),
    prec_malaria_waste = ifelse(previous_malaria == 1 & whz_ms_waste == 1, 1, 0)
  ) %>%
  summarise(prec_malaria_stunt = max(prec_malaria_stunt),
            prec_malaria_waste = max(prec_malaria_waste)) %>%
  ungroup() %>%
  mutate(prec_malaria_stunt = factor(prec_malaria_stunt),
         prec_malaria_waste = factor(prec_malaria_waste))

data_quarter_birth = data %>% group_by(id, agecat_birth) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  arrange(id, age)%>%
  ungroup()
#View(data_quarter_birth)

prevalence_outcomes_quarterly <-
  continuous_prevalence_outcomes_mutated %>% 
  group_by(id, agecat_birth) %>%
  mutate(haz_quarter = mean(haz, na.rm = TRUE), 
         waz_quarter = mean(waz, na.rm = TRUE), 
         whz_quarter = mean(whz, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(haz_ms_stunt_quarter = ifelse(haz_quarter< -2, 1, 0),
         haz_s_stunt_quarter = ifelse(haz_quarter< -3, 1, 0),
         waz_underwt_quarter = ifelse(waz_quarter< -2, 1, 0),
         whz_ms_waste_quarter = ifelse(whz_quarter< -2, 1, 0),
         whz_s_waste_quarter = ifelse(whz_quarter< -3, 1, 0))
#View(prevalence_outcomes_quarterly)

data_zscore_quarterly = merge(dplyr::select(prevalence_outcomes_quarterly, c("id", "agecat_birth", "haz_quarter", "whz_quarter", "waz_quarter")),
                              dplyr::select(data_quarter_birth, c("uniqueid", "id", "age", "agecat_birth", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
                              by = c("id", "agecat_birth")) %>%
  mutate(agecat_birth = as.factor(agecat_birth)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

data_zscore_quarterly = merge(data_zscore_quarterly,
                                  preceding_malaria,
                                  by = c("id", "agecat_birth"))

#View(data_zscore_quarterly)


data_prevalence_quarterly = merge(dplyr::select(prevalence_outcomes_quarterly, c("id", "agecat_birth", "haz_ms_stunt_quarter" ,"haz_s_stunt_quarter", "whz_ms_waste_quarter", "whz_s_waste_quarter", "waz_underwt_quarter")),
                                  dplyr::select(data_quarter_birth, c("uniqueid", "id", "age", "agecat_birth", "motherid", "Txarm", "rand_Txarm", main_covariate_list)),
                                  by = c("id", "agecat_birth")) %>%
  mutate(agecat_birth = as.factor(agecat_birth)) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

data_prevalence_quarterly = merge(data_prevalence_quarterly,
                              preceding_malaria,
                              by = c("id", "agecat_birth"))

#View(data_prevalence_quarterly)


# -----------------------------------------
# incidence, quarterly
# -----------------------------------------

data_incidence_3month = merge(data_incidence_3month,
                              dplyr::select(data_zscore_quarterly, c("id", "agecat_birth", "haz_quarter", "whz_quarter", "waz_quarter")),
                              by = c("id", "agecat_birth")) %>%
  mutate(agecat_birth = factor(agecat_birth,
                               levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"
                               ))) %>%
  mutate(incident_haz_ms_stunt_agecat_birth = ifelse(is.na(haz_quarter), NA, incident_haz_ms_stunt_agecat_birth),
         atrisk_haz_ms_stunt_agecat_birth = ifelse(is.na(haz_quarter), NA, atrisk_haz_ms_stunt_agecat_birth),
         incident_haz_s_stunt_agecat_birth = ifelse(is.na(haz_quarter), NA, incident_haz_s_stunt_agecat_birth),
         atrisk_haz_s_stunt_agecat_birth = ifelse(is.na(haz_quarter), NA, atrisk_haz_s_stunt_agecat_birth),
         incident_whz_ms_waste_agecat_birth = ifelse(is.na(whz_quarter), NA, incident_whz_ms_waste_agecat_birth),
         atrisk_whz_ms_waste_agecat_birth = ifelse(is.na(whz_quarter), NA, atrisk_whz_ms_waste_agecat_birth),
         incident_whz_s_waste_agecat_birth = ifelse(is.na(whz_quarter), NA, incident_whz_s_waste_agecat_birth),
         atrisk_whz_s_waste_agecat_birth = ifelse(is.na(whz_quarter), NA, atrisk_whz_s_waste_agecat_birth),
         incident_waz_underwt_agecat_birth = ifelse(is.na(waz_quarter), NA, incident_waz_underwt_agecat_birth),
         atrisk_waz_underwt_agecat_birth = ifelse(is.na(waz_quarter), NA, atrisk_waz_underwt_agecat_birth))

data_incidence_3month = merge(data_incidence_3month,
                              preceding_malaria,
                              by = c("id", "agecat_birth"))

#View(data_incidence_3month)


# -----------------------------------------
# incidence, biannual
# -----------------------------------------

data_biannual = dplyr::select(data, c("uniqueid", "id", "age", "age_6_12_month", "motherid", "Txarm", "rand_Txarm", main_covariate_list)) %>% group_by(id, age_6_12_month) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  arrange(id, age)%>%
  ungroup()

data_incidence_6month = merge(df_incidence_biannual,
                                data_biannual,
                                by = c("id", "age_6_12_month")) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

#View(data_incidence_6month)


# -----------------------------------------
# incidence, annual
# -----------------------------------------

data_annual = dplyr::select(data, c("uniqueid", "id", "age", "age_1_12", "motherid", "Txarm", "rand_Txarm", main_covariate_list)) %>% group_by(id, age_1_12) %>%
  arrange(desc(anymalaria), desc(incidentmalaria), desc(wet_season), incidentmalaria_season) %>%
  slice(1) %>%
  arrange(id, age)%>%
  ungroup()

data_incidence_12month = merge(df_incidence_annual,
                              data_annual,
                              by = c("id", "age_1_12"),) %>%
  group_by(id) %>%
  arrange(id, age) %>%
  ungroup()

#View(data_incidence_12month)


# -----------------------------------------
# save data
# -----------------------------------------

saveRDS(data_continuous, paste0(data_path,"analysis_data_continuous.RDS"))
saveRDS(data_zscore_monthly, paste0(data_path,"analysis_data_zscore_monthly.RDS"))
saveRDS(data_zscore_quarterly, paste0(data_path,"analysis_data_zscore_quarterly.RDS"))
saveRDS(data_prevalence_monthly, paste0(data_path,"analysis_data_prevalence_monthly.RDS"))
saveRDS(data_prevalence_quarterly, paste0(data_path,"analysis_data_prevalence_quarterly.RDS"))
saveRDS(data_incidence_1month, paste0(data_path,"analysis_data_incidence_monthly.RDS"))
saveRDS(data_incidence_3month, paste0(data_path,"analysis_data_incidence_quarterly.RDS"))
saveRDS(data_incidence_6month, paste0(data_path,"analysis_data_incidence_biannual.RDS"))
saveRDS(data_incidence_12month, paste0(data_path,"analysis_data_incidence_annual.RDS"))
saveRDS(data_velocity_1month, paste0(data_path,"analysis_data_velocity_1month.RDS"))
saveRDS(data_velocity_2month, paste0(data_path,"analysis_data_velocity_2month.RDS"))
saveRDS(data_velocity_3month, paste0(data_path,"analysis_data_velocity_3month.RDS"))
saveRDS(data_full_monthly_ceiling, paste0(data_path,"analysis_data_monthly_ceiling.RDS"))
saveRDS(data_full_monthly_round, paste0(data_path,"analysis_data_monthly_round.RDS"))


