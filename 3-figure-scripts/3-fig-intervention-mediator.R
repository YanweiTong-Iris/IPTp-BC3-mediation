################################################################
# IPTp and child growth
# Figures showing mediator-outcome results
# for z-scores and incidence 
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))


IM_results =  readRDS(paste0(results_path,"IM-MO-stratified/intervention_mediator_main_results_stratified.RDS"))  %>% 
  filter(!mediator %in% c(
    "birthweight", "GWC_Z", "antibacterial_binary", "betalactam_binary"
  )) %>% 
  mutate(mediator = ifelse(mediator=="IL_10RB", "IL-10RB", mediator)) %>%
  mutate(mediator = ifelse(mediator=="MMP_1", "MMP-1", mediator)) %>%
  mutate(mediator_remark = mediator) %>%
  mutate(mediator_remark = case_when(
    mediator == "anemia_28binary" ~ "Anemia",
    mediator == "gestational_weightchange" ~ "Gestational weight\nchange (kg)",
    mediator == "placentalmal" ~ "Placental malaria",
    mediator == "birthweight_kg" ~ "Birth weight (kg)",
    mediator == "birthlength" ~ "Birth length (cm)",
    mediator == "preterm" ~ "Pre-term birth",
    mediator == "LBW" ~ "Low birthweight",
    TRUE ~ mediator
  )) %>%
  mutate(
    mediator_remark = factor(mediator_remark, levels = c(
      "Anemia", 
      "Gestational weight\nchange (kg)",
      "Placental malaria",
      "ADA",
      "CCL11",
      "CCL19",
      "CCL28",
      "CD244",
      "CD5",
      "CD6",
      "CDCP1",
      "CX3CL1",
      #"CXCL5",
      "DNER",
      "IL-10RB",
      "IL10",
      "IL18",
      "MMP-1",
      "OPG",
      "TNFRSF9",
      "TWEAK",
      "SCF",
      "Pre-term birth",
      "Birth length (cm)",
      "Birth weight (kg)",
      "Low birthweight"
    ))
  ) %>% 
  mutate(gravidity_strata = case_when(
    gravidity_strata=="single"~ "Primigravidae",
    gravidity_strata=="multi"~ "Multigravidae",
    gravidity_strata=="all"~ "All"
  )) %>% 
  mutate(gravidity_strata = fct_rev(gravidity_strata))


IM_binary = IM_results %>% filter(mediator_remark %in% c("Anemia", 
                                                         "Placental malaria",
                                                         "Pre-term birth",
                                                         "Low birthweight"
))

IM_continuous = IM_results %>% filter(mediator_remark %in% c("Gestational weight\nchange (kg)", 
                                                             "Birth length (cm)",
                                                             "Birth weight (kg)"
))

IM_Olink = IM_results %>% filter(!mediator_remark %in% c(unique(IM_binary$mediator_remark),unique(IM_continuous$mediator_remark)))

# truncate CI bounds
IM_binary = IM_binary %>% mutate(
  unadjusted_lower_95CI = ifelse(unadjusted_lower_95CI<0, 0.101, unadjusted_lower_95CI)
)

binary_plot <- ggplot(IM_binary %>% mutate(mediator_remark = fct_rev(mediator_remark)), 
                      aes(y = unadjusted_point_estimate, x = mediator_remark)) +
  geom_hline(yintercept=1, linewidth = 0.3) +
  geom_point(aes(col = gravidity_strata, shape=gravidity_strata), position = position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin = unadjusted_lower_95CI, 
                     ymax = unadjusted_upper_95CI,
                     col = gravidity_strata),
                 position = position_dodge(width=0.5)) +
  scale_color_manual(values = c("#c29bcc", "#6b4575","black"))+
  scale_shape_manual(values = c(15, 17, 19)) + 
  scale_y_continuous(trans="log",
                     limits = c(0.13, 5),
                     breaks = c(0.25,0.5, 1, 2, 3, 4),
                     labels = c(0.25,0.5, 1, 2, 3, 4)) + 
  coord_flip() +
  theme_minimal() +
  #theme(legend.position = "bottom") +
  theme(legend.position = "none") +
  ylab("Incidence ratio (95% CI)")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = -0.28, size = 12))+
  ggtitle("a) Binary mediators")
binary_plot


continuous_plot <- ggplot(IM_continuous %>% mutate(mediator_remark = fct_rev(mediator_remark)), 
                          aes(y = unadjusted_point_estimate, x = mediator_remark)) +
  geom_hline(yintercept=0, linewidth = 0.3) +
  geom_point(aes(col = gravidity_strata, shape=gravidity_strata), position = position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin = unadjusted_lower_95CI, 
                     ymax = unadjusted_upper_95CI,
                     col = gravidity_strata),
                 position = position_dodge(width=0.5)) +
  scale_color_manual(values = c("#c29bcc", "#6b4575","black"))+
  scale_shape_manual(values = c(15, 17, 19)) + 
  scale_y_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
                     labels = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Mean difference (95% CI)")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = -0.35, size = 12)) +
  ggtitle("b) Continuous mediators")
continuous_plot

olink_plot <- ggplot(IM_Olink %>% mutate(mediator_remark = fct_rev(mediator_remark)), 
                     aes(y = unadjusted_point_estimate, x = mediator_remark)) +
  geom_hline(yintercept=0, linewidth = 0.3) +
  geom_point(aes(col = gravidity_strata, shape=gravidity_strata), position = position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin = unadjusted_lower_95CI, 
                     ymax = unadjusted_upper_95CI,
                     col = gravidity_strata),
                 position = position_dodge(width=0.5)) +
  scale_color_manual(values = c("#c29bcc", "#6b4575","black"))+
  scale_shape_manual(values = c(15, 17, 19)) + 
  scale_y_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
                     labels = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 10, b = 10),
        legend.spacing.x = unit(1, 'cm')) +
  ylab("Mean difference (95% CI)")+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = -0.22, size = 12)) +
  ggtitle("c) Inflammation-related mediators")
olink_plot


spacer <- rectGrob(gp = gpar(col = NA))

# Combine the plots with the spacer
combined_plot <- grid.arrange(binary_plot, spacer, continuous_plot, spacer, olink_plot, nrow = 5, heights = c(1, 0.15, 0.9, 0.15, 3.5))
combined_plot
ggsave(combined_plot, filename = paste0(figure_path, "plot-intervention-mediator-gravidity.pdf"),
       width=6, height= 10)
