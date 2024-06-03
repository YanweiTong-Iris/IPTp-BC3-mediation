################################################################
# IPTp and child growth
# Heatmaps showing mediation analysis results
# for Olink mediators
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(circlize)

mediation_zscore_3mo = readRDS(paste0(results_path,"/aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS")) %>%
  filter(interaction == 0, gravidae=="all") %>%
  filter(outcome %in% c("haz_quarter", "whz_quarter")) %>%
  dplyr::select(mediator, outcome, age_group,
                ACME_average, ACME_average_lower_CI, ACME_average_upper_CI, 
                ACME_p_val, ACME_adj_p,
                ADE_average, ADE_average_lower_CI, ADE_average_upper_CI) %>%
  mutate(signif = ifelse (ACME_adj_p < 0.05, "Significant", "Not significant")) %>% 
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months"))) %>%
  filter(!mediator %in% c("anemia_28binary", "antibacterial_binary", "betalactam_binary", "birthweight", 
                          "birthlength", "birthweight_kg", "gestational_weightchange", "LBW", "placentalmal","preterm")) %>%
  mutate(outcome_remark = case_when(outcome == "haz_quarter" ~ "Length-for-age Z",
                                    outcome == "whz_quarter" ~ "Weight-for-length Z")) %>%
  mutate(mediator = gsub("_", "-", mediator)) %>%
  mutate(mediator = ifelse(mediator == "F4E-BP1", "4E-BP1", mediator)) %>%
  mutate(mediator = ifelse(mediator == "LAP-TGF-beta-1", "LAP TGF-beta-1", mediator)) 


library(ComplexHeatmap)
library(circlize)

make_heatmap <- function(plot_data, outcome_name, outcome_type){
  heatmap_matrix = plot_data %>%
    filter(outcome_remark == outcome_name) %>%
    pivot_wider(values_from = ACME_average, 
                names_from = age_group, 
                values_fill = NA, 
                id_cols = mediator) %>%
    column_to_rownames(var = "mediator") %>%
    as.matrix(ncol = 6) 
  
  print(heatmap_matrix)
  
  if (outcome_type == "continuous") {
    col_fun = colorRamp2(c(-0.1, -0.05, 0, 0.05, 0.1),
                         c("#04388c", "#3b7eeb", "white", "#e38f22", "#6b3d01"))
  } else {
    col_fun = colorRamp2(c(0, 0.25, 1, 1.5, 2),
                         c("#04388c", "#3b7eeb", "white", "#e38f22", "#6b3d01"))
  }

  
  
  # Draw the heatmap with the significance annotations
  Heatmap(heatmap_matrix, 
          name = ifelse(outcome_type == "continuous", "Mean difference\nin z-score", "Incidence ratio"),
          cluster_rows = TRUE, 
          cluster_columns = FALSE,
          clustering_distance_rows = "euclidean",
          clustering_method_rows = "ward.D",
          row_names_gp = gpar(fontsize = 9),
          na_col = "grey90",
          col = col_fun,
          
          ) 

}

pdf(paste0(figure_path, "plot-Olink-mediation-heatmap-LAZ.pdf"), width=5, height=5)
make_heatmap(plot_data = mediation_zscore_3mo, outcome_name = "Length-for-age Z", outcome_type = "continuous")
dev.off()

pdf(paste0(figure_path, "plot-Olink-mediation-heatmap-WLZ.pdf"), width=5, height=5)
make_heatmap(plot_data = mediation_zscore_3mo, outcome_name = "Weight-for-length Z", outcome_type = "continuous")
dev.off()
