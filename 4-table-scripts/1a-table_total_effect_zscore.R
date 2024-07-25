################################################################
# IPTp and child growth
# Table of total effects on incidence of binary outcomes
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
options(digits=2)


#-----------------------------------------------------
# All gravidae
#-----------------------------------------------------

# Read in aim1 zscore mean difference data
zscore_mean_difference = readRDS(paste0(results_path, "aim1-stratified/aim1_zscore_quarterly_results_stratified.RDS")) %>% 
  filter(is.na(modifier_level), gravidity_strata == "all") %>% 
  rename(agecat_birth = age_group) %>% 
  dplyr::select(outcome, agecat_birth, point_estimate, lower_95CI, upper_95CI) %>% 
  distinct() %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "1 day-3 months", ">3-6 months", ">6-9 months",">9-12 months"
  ))) %>%
  filter(outcome %in% c("haz_quarter", "whz_quarter")) %>%
  mutate(mean_CI = paste0(point_estimate, " (", lower_95CI, ", ", upper_95CI, ")")) %>%
  dplyr::select(outcome, agecat_birth, mean_CI) %>%
  rename(mean_difference_CI = mean_CI)

zscore_mean_difference_haz = filter(zscore_mean_difference, outcome == "haz_quarter") 
zscore_mean_difference_whz = filter(zscore_mean_difference, outcome == "whz_quarter") 

# Read in quarterly zscore data 
zscore = readRDS(paste0(data_path, "analysis_data_zscore_quarterly.RDS")) %>%
  group_by(Txarm, agecat_birth) %>%
  summarise(N_haz = sum(!is.na(haz_quarter)),
            mean_haz_quarter = mean(haz_quarter, na.rm = TRUE),
            haz_SE = sd(haz_quarter, na.rm = TRUE)/sqrt(N_haz),
            N_whz = sum(!is.na(whz_quarter)),
            mean_whz_quarter = mean(whz_quarter, na.rm = TRUE),
            whz_SE = sd(whz_quarter, na.rm = TRUE)/sqrt(N_whz)) %>%
  mutate(haz_CI = paste0(round(mean_haz_quarter , 2), " (", round(mean_haz_quarter-qnorm(0.975)*haz_SE, 2), ", ", round(mean_haz_quarter+qnorm(0.975)*haz_SE, 2), ")"),
         whz_CI = paste0(round(mean_whz_quarter, 2), " (", round(mean_whz_quarter-qnorm(0.975)*whz_SE, 2), ", ", round(mean_whz_quarter+qnorm(0.975)*whz_SE, 2), ")")) %>%
  dplyr::select(Txarm, agecat_birth, N_haz, haz_CI, N_whz, whz_CI)

zscore_SP = zscore %>% filter(Txarm == "SP") %>%
  rename(N_SP_haz = N_haz,
         SP_haz_CI = haz_CI,
         N_SP_whz = N_whz,
         SP_whz_CI = whz_CI)

zscore_DP = zscore %>% filter(Txarm == "DP") %>%
  rename(N_DP_haz = N_haz,
         DP_haz_CI = haz_CI,
         N_DP_whz = N_whz,
         DP_whz_CI = whz_CI)


#combine the two tables
zscore_mean_difference_haz = merge(zscore_mean_difference_haz, zscore_SP[2:4], by = "agecat_birth")
zscore_mean_difference_haz = merge(zscore_mean_difference_haz, zscore_DP[2:4], by= "agecat_birth") %>%
  rename(N_SP = N_SP_haz,
         SP_CI = SP_haz_CI,
         N_DP = N_DP_haz,
         DP_CI = DP_haz_CI)


zscore_mean_difference_whz = merge(zscore_mean_difference_whz, zscore_SP[c(2, 5:6)], by = "agecat_birth")
zscore_mean_difference_whz = merge(zscore_mean_difference_whz, zscore_DP[c(2, 5:6)], by= "agecat_birth")%>%
  rename(N_SP = N_SP_whz,
         SP_CI = SP_whz_CI,
         N_DP = N_DP_whz,
         DP_CI = DP_whz_CI)


aim1_zscore = rbind(zscore_mean_difference_haz, zscore_mean_difference_whz) %>% 
  group_by(outcome) %>%
  arrange(agecat_birth, .by_group = TRUE) %>%
  ungroup()

write.csv(aim1_zscore, here::here(tables_path, "table_aim1_zscore.csv")) 

View(aim1_zscore)
