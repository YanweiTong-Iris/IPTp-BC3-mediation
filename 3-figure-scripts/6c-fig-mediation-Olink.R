################################################################
# IPTp and child growth
# Figures showing single mediator results 
# for z-scores and incidence 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# load results -------------------------------------------------------------------
result_Olink_incidence_3mo <- readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_incidence_results_3mo.RDS")) %>% 
  filter(interaction == 0 & gravidae=="all") %>%
  dplyr::select(mediator, outcome, age_group,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

result_Olink_Zscore_3mo <- readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS")) %>%
  filter(interaction == 0 & gravidae=="all") %>%
  dplyr::select(mediator, outcome, age_group,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

results_Olink_stunt = subset(result_Olink_incidence_3mo, outcome == "incident_haz_ms_stunt_agecat_birth")
results_Olink_waste = subset(result_Olink_incidence_3mo, outcome == "incident_whz_ms_waste_agecat_birth")

results_Olink_laz = subset(result_Olink_Zscore_3mo, outcome == "haz_quarter")
results_Olink_wlz = subset(result_Olink_Zscore_3mo, outcome == "whz_quarter")

# data processing functions -------------------------------------------------------------------

process_data <- function(input_data, outcome_type){
  output_data <- input_data %>% 
    filter(! mediator %in% c("anemia_36binary", "birthweight","anemia_28binary",
                             "birthlength", "birthweight_kg", "gestational_weightchange",
                             "placentalmal","preterm","LBW")) %>% 
    # filter(age_group %in% c("Birth","1 day-3 months")) %>% 
    pivot_longer(
      cols = c("ACME_average", "ACME_average_lower_CI", "ACME_average_upper_CI",
               "ADE_average", "ADE_average_lower_CI", "ADE_average_upper_CI"),
      names_to = c("measure", ".value"),
      names_pattern = "^(ACME|ADE)_(\\w+)"
    ) %>% 
    mutate(mediator = gsub("_", "-", mediator)) %>%
    mutate(mediator = fct_rev(mediator)) %>% 
    mutate(measure = ifelse(measure=="ACME", "Mediated effect", "Direct effect"))  

  if(outcome_type=="binary"){
    output_data <- output_data %>% 
      mutate(favors = case_when(
        average < 1 & average_upper_CI <1 ~ "DP promotes\ngrowth",
        average_lower_CI <=1 & average_upper_CI >= 1 ~ "Null",
        average > 1 & average_upper_CI >1 ~ "SP promotes\ngrowth" 
      )) %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) 
  }
  
  
  if(outcome_type=="continuous"){
    output_data <- output_data %>% 
      mutate(favors = case_when(
        average < 0 & average_upper_CI <0 ~ "SP promotes\ngrowth",
        average_lower_CI <=0 & average_upper_CI >= 0 ~ "Null",
        average > 0 & average_upper_CI >0 ~ "DP promotes\ngrowth" 
      )) %>% 
      mutate(favors = factor(favors, levels = c("SP promotes\ngrowth", "Null", "DP promotes\ngrowth"))) 
  }
  
  
  return(output_data)
  
}


################################
# 3 month
################################

# stunting -------------------------------------------------------------------

stunt_l <- process_data(input_data = results_Olink_stunt,
                        outcome_type = "binary") 

stunt_plot <- ggplot(stunt_l %>% filter(measure=="Mediated effect"), aes(x = mediator, y = average, group = age_group)) + 
  geom_hline(yintercept = 1, linewidth=0.3) +
  geom_point(aes(col = favors), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_y_continuous(trans="log",
                     limits = c(0.6,1.75),
                     breaks = c(0.75, 1, 1.25, 1.5),
                     labels = c(0.75, 1, 1.25, 1.5)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2))

stunt_plot

ggsave(stunt_plot, filename = paste0(figure_path, "plot-mediation-Olink-stunt.pdf"),
       width=8, height=4)

# LAZ -------------------------------------------------------------------
laz_l <- process_data(input_data = results_Olink_laz,
                      outcome_type = "continuous")

laz_plot <- ggplot(laz_l %>% filter(measure=="Mediated effect") , aes(x = mediator, y = average, group = age_group)) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_point(aes(col = favors), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_y_continuous(
    limits = c(-0.23, 0.2),
    breaks = c( -0.2,-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2),
    labels = c( -0.2,-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        # axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position= "bottom",
        legend.spacing.x = unit(0, "cm"),
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2))

laz_plot

ggsave(laz_plot, filename = paste0(figure_path, "plot-mediation-Olink-laz.pdf"),
       width=11, height=4)

# wasting -------------------------------------------------------------------
waste_l <- process_data(input_data = results_Olink_waste,
                        outcome_type = "binary")

waste_plot <- ggplot(waste_l  %>% filter(measure=="Mediated effect"), aes(x = mediator, y = average, group = age_group)) + 
  geom_hline(yintercept = 1, linewidth=0.3) +
  geom_point(aes(col = favors), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_y_continuous(trans="log",
                     limits = c(0.5,3),
                     breaks = c(0.5, 1, 1.5, 2, 3),
                     labels = c(0.5, 1, 1.5, 2, 3)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2)) 

waste_plot

ggsave(waste_plot, filename = paste0(figure_path, "plot-mediation-Olink-waste.pdf"),
       width=8, height=4)


# WLZ -------------------------------------------------------------------
wlz_l <- process_data(input_data = results_Olink_wlz,
                      outcome_type = "continuous")

wlz_plot <- ggplot(wlz_l  %>% filter(measure=="Mediated effect"), aes(x = mediator, y = average, group = age_group)) +
  geom_hline(yintercept = 0, linewidth=0.3) +
  geom_point(aes(col = favors), position = position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(ymin = average_lower_CI, ymax = average_upper_CI,
                     col = favors),
                 position = position_dodge(width=0.5)) + 
  facet_grid(~age_group) +
  scale_color_manual(values = c("#164c9e", "#909190","#9e161d"), drop = FALSE) +
  scale_shape_manual(values = c(1, 16)) + 
  scale_y_continuous(
    limits = c(-0.23, 0.2),
    breaks = c( -0.2,-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2),
    labels = c( -0.2,-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2)) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        # axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position= "bottom",
        legend.spacing.x = unit(0, "cm"),
        panel.grid.minor = element_line(color = "white", size = 0),
        panel.grid.major = element_line(size=0.2)) 

wlz_plot

ggsave(wlz_plot, filename = paste0(figure_path, "plot-mediation-Olink-wlz.pdf"),
       width=11, height=4)



