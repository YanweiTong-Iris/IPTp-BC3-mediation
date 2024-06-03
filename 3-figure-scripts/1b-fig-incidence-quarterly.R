################################################################
# IPTp and child growth
# Script for making incidence plots 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))


inc.res.quarterly = readRDS(paste0(results_path, "incidence_estimates_quarterly.RDS")) %>% 
  rename(agevar = agecat_birth,
         agecat = age_level) %>% 
  mutate(strat_var = ifelse(strat_var == "agecat_birth", "agecat", strat_var)) %>% 
  mutate(agecat = factor(agecat, levels = c("Birth","1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>% 
  filter(strat_var=="agecat") %>% 
  filter(outcome!="waz_underwt") %>% 
  mutate(outcome_lab = case_when(
    outcome=="haz_ms_stunt" ~ "Any stunting",
    outcome=="haz_s_stunt" ~ "Severe stunting",
    outcome=="whz_ms_waste" ~ "Any wasting",
    outcome=="whz_s_waste" ~ "Severe wasting"
  )) %>% 
  mutate(outcome_class = case_when(
    outcome=="haz_ms_stunt" ~ "C) Stunting",
    outcome=="haz_s_stunt" ~ "C) Stunting",
    outcome=="whz_ms_waste" ~ "D) Wasting",
    outcome=="whz_s_waste" ~ "D) Wasting"
  )) %>% 
  mutate(outcome_severity = case_when(
    outcome=="haz_ms_stunt" ~ "Any stunting/wasting",
    outcome=="haz_s_stunt" ~ "Severe stunting/wasting",
    outcome=="whz_ms_waste" ~ "Any stunting/wasting",
    outcome=="whz_s_waste" ~ "Severe stunting/wasting"
  )) %>% 
  mutate(agecat_plot = case_when(
    agecat == "Birth"~ "Birth",
    agecat == "1 day-3 months"~ "0-3",
    agecat == ">3-6 months"~ "4-6",
    agecat == ">6-9 months"~ "7-9",
    agecat == ">9-12 months"~ "10-12",
  )) %>% 
  mutate(agecat_plot = factor(agecat_plot, levels = c(
    "Birth", "0-3", "4-6", "7-9", "10-12"
  )))

incidence_plot =  ggplot(inc.res.quarterly , 
                         aes(x = agecat_plot, y = est)) + 
  facet_wrap(~outcome_class) +
  geom_linerange(aes(ymin = lb, ymax = ub, color = outcome_severity), size=1.5, alpha=0.75) +
  geom_point(aes(color = outcome_severity), stat = "identity", size=3) +
  ylab("Incidence (95% CI)") +
  xlab("Age, months") +
  scale_color_manual(values = c("#a18a8a", "black"))+
  scale_y_continuous(labels = c(0, 5, 10, 15, 20, 25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(size=13),
        panel.spacing = unit(3, "lines"),
        strip.text.x = element_text(hjust = 0, margin=margin(l=0, b=3)))


incidence_plot

ggsave(incidence_plot, filename = paste0(here::here(), "/5-figures/plot-incidence-3m.pdf"),
       width=6, height=4)




