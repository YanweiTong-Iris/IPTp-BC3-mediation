################################################################
# IPTp and child growth
# Figures showing single mediator results 
# for z-scores and incidence 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# load results -------------------------------------------------------------------
result_single_mediator_zscore_3mo = 
  readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS")) %>%
  filter(interaction == 0) %>%
  dplyr::select(mediator, outcome, age_group, gravidae,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))


result_single_mediator_incidence_3mo = readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_incidence_results_3mo.RDS")) %>%
  filter(interaction == 0) %>%
  dplyr::select(mediator, outcome, age_group, gravidae,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

results_single_mediator_laz = subset(result_single_mediator_zscore_3mo, outcome == "haz_quarter")
results_single_mediator_waz = subset(result_single_mediator_zscore_3mo, outcome == "waz_quarter")
results_single_mediator_wlz = subset(result_single_mediator_zscore_3mo, outcome == "whz_quarter")

results_single_mediator_ms_stunt = subset(result_single_mediator_incidence_3mo, 
                                          outcome == "incident_haz_ms_stunt_agecat_birth")
results_single_mediator_ms_waste = subset(result_single_mediator_incidence_3mo, 
                                          outcome == "incident_whz_ms_waste_agecat_birth")
results_single_mediator_underwt = subset(result_single_mediator_incidence_3mo, 
                                         outcome == "incident_waz_underwt_agecat_birth")

total_effect_res_Z =  readRDS(paste0(results_path,"aim1-stratified/aim1_zscore_quarterly_results_stratified.RDS"))
total_effect_res_inc3m =  readRDS(paste0(results_path,"aim1-stratified/aim1_incidence_3mo_results_stratified.RDS"))

# data processing functions -------------------------------------------------------------------

process_data <- function(input_data, outcome, outcome_type){
  output_acme <- input_data %>% 
    filter(gravidae=="all") %>% 
    filter(mediator!="anemia_36binary") %>% 
    filter(mediator!="birthweight") %>%
    pivot_longer(
      cols = c("ACME_average", "ACME_average_lower_CI", "ACME_average_upper_CI",
               "ADE_average", "ADE_average_lower_CI", "ADE_average_upper_CI"),
      names_to = c("measure", ".value"),
      names_pattern = "^(ACME|ADE)_(\\w+)"
    ) %>% 
    mutate(mediator_label = case_when(
      mediator == "anemia_28binary" ~ "Anemia",
      mediator == "gestational_weightchange" ~ "Gestational\nweight change",
      mediator == "LBW" ~ "Low birth weight",
      mediator == "placentalmal" ~ "Placental malaria",
      mediator == "preterm" ~ "Preterm birth",
      mediator == "SCF" ~ "SCF",
      mediator == "CCL11" ~ "CCL11",
      mediator == "CCL19" ~ "CCL19",
      mediator == "CCL28" ~ "CCL28",
      mediator == "CD244" ~ "CD244",
      mediator == "CD5" ~ "CD5",
      mediator == "CD6" ~ "CD6",
      mediator == "CDCP1" ~ "CDCP1",
      mediator == "CXCL5" ~ "CXCL5",
      mediator == "DNER" ~ "DNER",
      mediator == "IFN_gamma" ~ "IFN_gamma",
      mediator == "IL_12B" ~ "IL_12B",
      mediator == "IL18" ~ "IL18",
      mediator == "LIF_R" ~ "LIF_R",
      mediator == "OPG" ~ "OPG",
      mediator == "PD_L1" ~ "PD_L1",
      mediator == "TNF" ~ "TNF",
      mediator == "TNFRSF9" ~ "TNFRSF9",
      mediator == "birthlength" ~ "Birth length",
      mediator == "birthweight_kg" ~ "Birth weight",
      mediator == "Total effect" ~ "Total effect"
    )) %>% 
    mutate(age_group = factor(age_group, levels = c(
      "Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"
    ))) %>% 
    mutate(mediator_label = factor(mediator_label, levels = c("Anemia", 
                                                              "Gestational\nweight change",
                                                              "Placental malaria",
                                                              "CCL11",
                                                              "CCL19",
                                                              "CCL28",
                                                              "CD244",
                                                              "CD5",
                                                              "CD6",
                                                              "CDCL5",
                                                              "CDCP1",
                                                              "CXCL5",
                                                              "DNER",
                                                              "IFN_gamma",
                                                              "IL_12B",
                                                              "IL18",
                                                              "LIF_R",
                                                              "OPG",
                                                              "PD_L1",
                                                              "TNF",
                                                              "TNFRSF9",
                                                              "SCF",
                                                              "Preterm birth",
                                                              "Birth length",
                                                              "Birth weight",
                                                              "Low birth weight",
                                                              "Total effect"))) %>% 
    filter(!mediator %in% c("antibacterial_binary", "betalactam_binary")) %>% 
    mutate(measure = ifelse(measure=="ACME", "Mediated effect", "Direct effect")) 
  
  # process total effects 
  if(outcome=="haz") output_sub = total_effect_res_Z %>% filter(outcome=="haz_quarter")
  if(outcome=="whz") output_sub = total_effect_res_Z %>% filter(outcome=="whz_quarter")
  if(outcome=="stunting") output_sub = total_effect_res_inc3m %>% filter(outcome=="incident_haz_ms_stunt_agecat_birth")
  if(outcome=="wasting") output_sub = total_effect_res_inc3m %>% filter(outcome=="incident_whz_ms_waste_agecat_birth")
  
  output_total = output_sub %>% filter(is.na(modifier_name)) %>% 
    filter(gravidity_strata=="all") %>% 
    dplyr::select(age_group, point_estimate, lower_95CI, upper_95CI) %>% 
    rename(average = point_estimate, 
           average_lower_CI = lower_95CI,
           average_upper_CI = upper_95CI) %>% 
    mutate(mediator="Total effect",
           measure = "Total effect") 
  
  
  output_data <- bind_rows(output_acme, output_total)  %>% 
    mutate(mediator_label = case_when(
      mediator == "anemia_28binary" ~ "Anemia",
      mediator == "gestational_weightchange" ~ "Gestational\nweight change",
      mediator == "LBW" ~ "Low birth weight",
      mediator == "placentalmal" ~ "Placental malaria",
      mediator == "preterm" ~ "Preterm birth",
      mediator == "SCF" ~ "SCF",
      mediator == "ADA" ~ "ADA",
      mediator == "CCL11" ~ "CCL11",
      mediator == "CCL19" ~ "CCL19",
      mediator == "CCL28" ~ "CCL28",
      mediator == "CD244" ~ "CD244",
      mediator == "CD5" ~ "CD5",
      mediator == "CD6" ~ "CD6",
      mediator == "CDCP1" ~ "CDCP1",
      mediator == "CXCL5" ~ "CXCL5",
      mediator == "DNER" ~ "DNER",
      mediator == "IFN_gamma" ~ "IFN-gamma",
      mediator == "IL10" ~ "IL10",
      mediator == "IL_12B" ~ "IL-12B",
      mediator == "IL18" ~ "IL18",
      mediator == "LIF_R" ~ "LIF-R",
      mediator == "OPG" ~ "OPG",
      mediator == "PD_L1" ~ "PD-L1",
      mediator == "TNF" ~ "TNF",
      mediator == "TNFRSF9" ~ "TNFRSF9",
      mediator == "TWEAK" ~ "TWEAK",
      mediator == "birthlength" ~ "Birth length",
      mediator == "birthweight_kg" ~ "Birth weight",
      mediator == "Total effect" ~ "Total effect"
    )) %>% 
    mutate(age_group = factor(age_group, levels = c(
      "Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"
    ))) %>%
    mutate(mediator_label = factor(mediator_label, levels = c("Anemia", 
                                                              "Gestational\nweight change",
                                                              "Placental malaria",
                                                              "ADA",
                                                              "CCL11",
                                                              "CCL19",
                                                              "CCL28",
                                                              "CD244",
                                                              "CD5",
                                                              "CD6",
                                                              "CDCL5",
                                                              "CDCP1",
                                                              "CXCL5",
                                                              "DNER",
                                                              "IFN_gamma",
                                                              "IL10",
                                                              "IL-12B",
                                                              "IL18",
                                                              "LIF-R",
                                                              "OPG",
                                                              "PD-L1",
                                                              "TNF",
                                                              "TNFRSF9",
                                                              "TWEAK",
                                                              "SCF",
                                                              "Preterm birth",
                                                              "Birth length",
                                                              "Birth weight",
                                                              "Low birth weight",
                                                              "Total effect"))) %>% 
    mutate(measure = ifelse(measure=="ACME", "Mediated effect", measure)) 
  
  if(outcome_type=="binary"){
    
    output_data <- output_data %>% 
      filter(measure!="Direct effect") %>%
      mutate(harmful = ifelse( mediator %in% c("LBW", "anemia_28binary", 
                                               "placentalmal", "preterm"), 1, 0)) %>% 
      mutate(favors = case_when(
        average < 1 & average_upper_CI <1 & harmful==1 ~ "DP promotes\ngrowth",
        average < 1 & average_upper_CI <1 & harmful==0 ~ "SP promotes\ngrowth",
        average_lower_CI <=1 & average_upper_CI >= 1 ~ "Null",
        average > 1 & average_upper_CI >1 & harmful==1 ~ "SP promotes\ngrowth" ,
        average > 1 & average_upper_CI >1 & harmful==0 ~ "DP promotes\ngrowth" 
        
      ))  %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) %>% 
      mutate(signif = case_when(
        average_lower_CI<1 & average_upper_CI<1 ~ 1,
        average_lower_CI>1 & average_upper_CI>1 ~ 1,
        TRUE ~ 0
      )) %>% 
      group_by(mediator_label) %>% 
      mutate(nsignif = sum(signif)) %>% 
      ungroup() %>% 
      filter(mediator_label == "Total effect" |nsignif>0)
  }
  
  
  if(outcome_type=="continuous"){
    output_data <- output_data %>% 
      filter(measure!="Direct effect") %>%
      mutate(favors = case_when(
        average < 0 & average_upper_CI <0 ~ "SP promotes\ngrowth",
        average_lower_CI <=0 & average_upper_CI >= 0 ~ "Null",
        average > 0 & average_upper_CI >0 ~ "DP promotes\ngrowth" 
      )) %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) %>% 
      mutate(signif = case_when(
        average_lower_CI<0 & average_upper_CI<0 ~ 1,
        average_lower_CI>0 & average_upper_CI>0 ~ 1,
        TRUE ~ 0
      )) %>% 
      group_by(mediator_label) %>% 
      mutate(nsignif = sum(signif)) %>% 
      ungroup() %>% 
      filter(mediator_label == "Total effect" | nsignif>0)
  }
  
  
  return(output_data)
  
}


################################
# incidence plot
################################

# stunting -------------------------------------------------------------------

stunt_l <- process_data(input_data = results_single_mediator_ms_stunt,
                        outcome = "stunting",
                        outcome_type = "binary") 
drop_stunt <- which(stunt_l$age_group=="Birth" & stunt_l$mediator_label == "Birth length")
stunt_l <- stunt_l[-drop_stunt,]

# truncate CIs
stunt_l <- stunt_l %>% mutate(
  average_upper_CI = ifelse(average_upper_CI>3, 2.99, average_upper_CI)
)


stunt_plot <- ggplot(stunt_l %>% filter(measure!="Direct effect"), 
                     aes(x = mediator_label, y = average, group = measure)) + 
  geom_hline(yintercept = 1, linewidth=0.3) +
  geom_point(aes(col = favors, shape=measure), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_y_continuous(trans="log", 
                     # limits = c(0.6,5),
                     limits=c(-.6, 3),
                     breaks = c(0.75, 1, 1.5, 2, 3, 5),
                     labels = c(0.75, 1, 1.5, 2, 3, 5)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2),
        plot.title = element_text(hjust = -0.18, size = 12)) +
  ggtitle("A) Stunting")

stunt_plot



# wasting -------------------------------------------------------------------
waste_l <- process_data(input_data = results_single_mediator_ms_waste,
                        outcome = "wasting",
                        outcome_type = "binary")
drop_waste <- which(waste_l$age_group=="Birth" & waste_l$mediator_label == "Birth weight")
waste_l <- waste_l[-drop_waste,]

# truncate CIs
waste_l <- waste_l %>% mutate(
  average_upper_CI = ifelse(average_upper_CI>3, 2.99, average_upper_CI)
)

waste_plot <- ggplot(waste_l  %>% filter(measure!="Direct effect"), 
                     aes(x = mediator_label, y = average, group = measure)) + 
  geom_hline(yintercept = 1, linewidth=0.3) +
  geom_point(aes(col = favors, shape = measure), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_y_continuous(trans="log", 
                     # limits = c(0.6,5),
                     # limits=c(-.6, 3),
                     breaks = c(0.75, 1, 1.5, 2, 3, 5),
                     labels = c(0.75, 1, 1.5, 2, 3, 5)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2),
        plot.title = element_text(hjust = -0.18, size = 12)) +
  
  ggtitle("B) Wasting")

waste_plot


combined_inc_plot <- grid.arrange(stunt_plot, waste_plot, 
                                  nrow=2,ncol=1, heights = c(4,3.4))

ggsave(combined_inc_plot, filename = paste0(figure_path, "plot-mediation-inc.png"),
       width=7, height= 3.5)


################################
# Z-score plot
################################

# LAZ -------------------------------------------------------------------
laz_l <- process_data(input_data = results_single_mediator_laz,
                      outcome = "haz",
                      outcome_type = "continuous")
drop_laz <- which(laz_l$age_group=="Birth" & laz_l$mediator_label == "Birth length")
laz_l <- laz_l[-drop_laz,]

laz_plot <- ggplot(laz_l %>% filter(measure!="Direct effect" &
                                      mediator!="LBW") , 
                   aes(x = mediator_label, y = average, group = measure)) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_point(aes(col = favors, shape=measure), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_y_continuous(
    limits = c(-0.45, 0.45),
    breaks = c( -0.4,-0.2, 0, 0.2, 0.4),
    labels = c( -0.4,-0.2, 0, 0.2, 0.4)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        legend.title = element_blank(),
        legend.position= "none",
        legend.spacing.x = unit(0, "cm"),
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2)) +
  ggtitle("A) Length-for-age Z-score")
laz_plot


# WLZ -------------------------------------------------------------------
wlz_l <- process_data(input_data = results_single_mediator_wlz,
                      outcome = "whz",
                      outcome_type = "continuous")
drop_wlz <- which(wlz_l$age_group=="Birth" & wlz_l$mediator_label == "Birth weight")
wlz_l <- wlz_l[-drop_wlz,]

wlz_plot <- ggplot(wlz_l  %>% filter(measure!="Direct effect" &
                                       mediator!="LBW"),
                   aes(x = mediator_label, y = average, group = measure)) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_point(aes(col = favors, shape=measure), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_y_continuous(
    limits = c(-0.45, 0.45),
    breaks = c( -0.4,-0.2, 0, 0.2, 0.4),
    labels = c( -0.4,-0.2, 0, 0.2, 0.4)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        legend.title = element_blank(),
        legend.position= "none",
        legend.spacing.x = unit(0, "cm"),
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2)) +
  ggtitle("B) Weight-for-length Z-score")

wlz_plot

combined_Z_plot <- grid.arrange(laz_plot,wlz_plot, 
                              nrow=2,ncol=1)

ggsave(combined_Z_plot, filename = paste0(figure_path, "plot-mediation-Z.png"),
       width=7, height=4.5)


### save legend

legend_plot <- ggplot(wlz_l , aes(x = mediator_label, y = average, group = age_group)) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_point(aes(col = favors, shape=measure), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_y_continuous(
    breaks = c(-0.3, -0.2, -0.1, 1, 0.1, 0.2, 0.3, 0.4),
    labels = c(-0.3, -0.2, -0.1, 1, 0.1, 0.2, 0.3, 0.4)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        legend.title = element_blank(),
        legend.position= "bottom",
        legend.spacing.x = unit(0, "cm")) +
  ggtitle("B) Weight-for-length Z-score")

ggsave(legend_plot, filename = paste0(figure_path, "plot-mediation-legend.pdf"),
       width=9, height=5.75)


#grid.text("Title for Plot 1", x = unit(0.25, "npc"), y = unit(0.98, "npc"))



