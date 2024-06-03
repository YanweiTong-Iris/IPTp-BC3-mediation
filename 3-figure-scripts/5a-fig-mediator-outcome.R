################################################################
# IPTp and child growth
# Figures showing mediator-outcome results
# for z-scores  

# Non-olink mediators
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

MO_zscore_3mo = readRDS(paste0(results_path,"IM-MO-stratified/mediator_outcome_zscore_results_3mo_stratified.RDS")) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

# rescale gestational weight change to improve clarity of plot
MO_zscore_3mo = MO_zscore_3mo %>% 
  mutate(point_estimate = ifelse(independent_variable %in% c("gestational_weightchange", "birthlength"), point_estimate*5, point_estimate),
         lower_95CI = ifelse(independent_variable %in% c("gestational_weightchange", "birthlength"), lower_95CI*5, lower_95CI),
         upper_95CI = ifelse(independent_variable %in% c("gestational_weightchange", "birthlength"), upper_95CI*5, upper_95CI))

# define data processing function -------------------------------------------------------------------
process_data <- function(data, outcome_type){
  output_data <- data %>%  filter(independent_variable %in% 
                                    c("anemia_28binary", 
                                      "birthlength",
                                      "birthweight_kg", 
                                      "gestational_weightchange",
                                      "placentalmal",
                                      "preterm"))%>% 
    mutate(independent_variable_label = case_when(
      independent_variable == "anemia_28binary" ~ "Anemia",
      independent_variable == "gestational_weightchange" ~ "Gestational\nweight\nchange (kg)",
      independent_variable == "LBW" ~ "Low birth\nweight",
      independent_variable == "placentalmal" ~ "Placental\nmalaria",
      independent_variable == "preterm" ~ "Preterm\nbirth",
      independent_variable == "birthweight_kg" ~ "Birth\nweight\n(kg)",
      independent_variable == "birthlength" ~ "Birth\nlength\n(cm)",
      TRUE ~ independent_variable
    )) %>%
    mutate(independent_variable_label = factor(
      independent_variable_label,
      levels = c(
        "Anemia",
        "Gestational\nweight\nchange (kg)",
        "Placental\nmalaria",
        "Preterm\nbirth",
        "Birth\nweight\n(kg)",
        "Birth\nlength\n(cm)"
      )
    )) %>%
    mutate(age_group = fct_rev(age_group))  %>% 
    filter(gravidae!="all") %>% 
    mutate(gravidae = ifelse(gravidae=="single", "Primigravidae", "Multigravidae")) %>% 
    mutate(gravidae = factor(gravidae, levels = c("Primigravidae", "Multigravidae")))
  
  return(output_data)
}


# process Z-score data -------------------------------------------------------------------
zscore_plotdata = process_data(data = MO_zscore_3mo, outcome_type="continuous") %>% 
  filter(dependent_variable !="waz_quarter") %>% 
  mutate(dependent_variable_label = case_when(
    dependent_variable == "haz_quarter" ~ "length-for-age Z",
    dependent_variable == "whz_quarter" ~ "weight-for-length Z"
  ))





# make plots -------------------------------------------------------------------
# drop birth length at birth
drop <- which(zscore_plotdata$independent_variable=="birthlength" & zscore_plotdata$age_group=="Birth" &
                zscore_plotdata$dependent_variable=="haz_quarter")
zscore_plotdata <- zscore_plotdata[-drop,]

Zscore_plot <- ggplot(zscore_plotdata %>% 
                        mutate(age_group = fct_rev(age_group)) ,
                      aes(x = age_group, y = point_estimate)) +
  geom_point(aes(col = gravidae), position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = lower_95CI, ymax = upper_95CI,
                     col = gravidae),
                 position = position_dodge(width = 0.6)) +
  facet_grid(dependent_variable_label~independent_variable_label, 
             scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#8a8787","black"), drop = FALSE) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing = unit(0.5, "cm", data = NULL),
        panel.grid.minor.x = element_blank() ,
        panel.grid.major.x = element_blank() ,
        axis.text.x = element_text(size=7,angle = 90, vjust=0.5, hjust=1))+
  xlab("Age Group") +
  ylab("Mean difference in Z-score (95% CI)") 
Zscore_plot

ggsave(Zscore_plot, filename = paste0(figure_path, "plot-med-outcome-zscore.pdf"),
       width=9, height=6)
ggsave(Zscore_plot, filename = paste0(figure_path, "plot-med-outcome-zscore.png"),
       width=9, height=5.5)


