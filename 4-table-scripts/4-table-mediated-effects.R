################################################################
# IPTp and child growth
# Table of mediated effects
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

#prevent using scientific notations
options(scipen = 999)

result_single_mediator_zscore = readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS"))%>%
  mutate_if(is.numeric, round, digits=4) %>%
  filter(interaction == 0) %>%
  filter(outcome %in% c("haz_quarter", "whz_quarter")) %>%
  dplyr::select(mediator, outcome, n, age_group, gravidae, 
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>%
  # filter(mediator %in% c(
  #   "anemia_28binary", "gestational_weightchange", 
  #   "placentalmal", "preterm",
  #   "birthlength","birthweight_kg"
  # )) %>% 
  mutate(mediator_label = case_when(
    mediator == "anemia_28binary" ~ "Anemia",
    mediator == "gestational_weightchange" ~ "Gestational weight change (kg)",
    mediator == "placentalmal" ~ "Placental malaria",
    mediator == "preterm" ~ "Pre-term birth",
    mediator == "birthlength" ~ "Birth length (cm)",
    mediator == "birthweight_kg" ~ "Birth weight (kg)"
  )) %>% 
  mutate(
    mediator_label = factor(mediator_label, levels = c(
      "Anemia", 
      "Gestational weight change (kg)",
      "Placental malaria",
      "Pre-term birth",
      "Birth length (cm)",
      "Birth weight (kg)"
    ))
  ) %>% mutate(ACME_CI = paste0(ACME_average, " (", ACME_average_lower_CI, ", ", ACME_average_upper_CI, ")")) %>%
  dplyr::select(gravidae, outcome, mediator, mediator_label, age_group, n, ACME_CI) %>%
  mutate(n = as.integer(n))


result_single_mediator_zscore = with(result_single_mediator_zscore, result_single_mediator_zscore[order(outcome, mediator, age_group),])


table_mediated_effects_primi = result_single_mediator_zscore %>% filter(gravidae == "single")
write.csv(table_mediated_effects_primi,  paste0(tables_path,"table_mediated_effects_primi.csv"))

table_mediated_effects_multi = result_single_mediator_zscore %>% filter(gravidae == "multi")
write.csv(table_mediated_effects_multi,  paste0(tables_path,"table_mediated_effects_multi.csv"))

table_mediated_effects_all = result_single_mediator_zscore %>% filter(gravidae == "all")
write.csv(table_mediated_effects_all,  paste0(tables_path,"table_mediated_effects_all.csv"))
