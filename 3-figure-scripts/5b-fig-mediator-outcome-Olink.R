################################################################
# IPTp and child growth
# Figures showing mediator-outcome results
# for Olink mediators
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(circlize)

MO_zscore_3mo = readRDS(paste0(results_path,"IM-MO-stratified/mediator_outcome_zscore_results_3mo_stratified.RDS")) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

zscore_plot_data <- MO_zscore_3mo %>% filter(!independent_variable %in% 
                                  c("anemia_28binary", 
                                    "birthweight",
                                    "birthlength",
                                    "birthweight_kg", 
                                    "gestational_weightchange",
                                    "LBW",
                                    "placentalmal",
                                    "preterm")) %>% 
  filter(gravidae == "all") %>% 
  mutate(signif = ifelse (adj_p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(outcome = case_when(
    outcome_remark == "height-for-age z score" ~ "Length-for-age Z",
    outcome_remark == "weight-for-age z score" ~ "Weight-for-age Z",
    outcome_remark == "weight-for-height z score" ~ "Weight-for-length Z"
  )) %>% 
  mutate(independent_variable = gsub("_", "-", independent_variable)) %>%
  mutate(independent_variable = ifelse(independent_variable == "F4E-BP1", "4E-BP1", independent_variable)) %>%
  mutate(independent_variable = ifelse(independent_variable == "LAP-TGF-beta-1", "LAP TGF-beta-1", independent_variable)) %>%
  mutate(independent_variable = fct_rev(independent_variable))
  


MO_binary_3mo = readRDS(paste0(results_path,"IM-MO-stratified/mediator_outcome_incidence_results_3mo_stratified.RDS")) %>%
  mutate(age_group = factor(age_group, levels = c("Birth", "1 day-3 months", ">3-6 months", ">6-9 months", ">9-12 months")))

binary_plot_data <- MO_binary_3mo %>%filter(!independent_variable %in% 
                                              c("anemia_28binary", 
                                                "birthweight",
                                                "birthlength",
                                                "birthweight_kg", 
                                                "gestational_weightchange",
                                                "LBW",
                                                "placentalmal",
                                                "preterm")) %>% 
  filter(gravidae == "all", dependent_variable %in% c("incident_haz_ms_stunt_agecat_birth", "incident_whz_ms_waste_agecat_birth")) %>% 
  mutate(signif = ifelse (adj_p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(outcome = case_when(
    outcome_remark == "incidence: moderate to severe stunting" ~ "Moderate-to-severe stunting",
    outcome_remark == "incidence: moderate to severe wasting" ~ "Moderate-to-severe wasting"
  )) %>% 
  mutate(independent_variable = gsub("_", "-", independent_variable)) %>%
  mutate(independent_variable = ifelse(independent_variable == "F4E-BP1", "4E-BP1", independent_variable)) %>%
  mutate(independent_variable = ifelse(independent_variable == "LAP-TGF-beta-1", "LAP TGF-beta-1", independent_variable)) %>%
  mutate(independent_variable = fct_rev(independent_variable)) %>%
  filter(upper_95CI < 10)



library(ComplexHeatmap)

make_heatmap <- function(plot_data, outcome_name, outcome_type){
  heatmap_matrix = plot_data %>%
    filter(outcome == outcome_name) %>%
    pivot_wider(values_from = point_estimate, 
                names_from = age_group, 
                values_fill = NA, 
                id_cols = independent_variable) %>%
    column_to_rownames(var = "independent_variable") %>%
    as.matrix(ncol = 6) 
  
  if (outcome_type == "continuous") {
    col_fun = colorRamp2(c(-0.4, -0.2, 0, 0.2, 0.4),
                         c("#04388c", "#3b7eeb", "white", "#e38f22", "#6b3d01"))
  } else {
    col_fun = colorRamp2(c(0, 0.25, 1, 1.5, 2),
                         c("#04388c", "#3b7eeb", "white", "#e38f22", "#6b3d01"))
  }
  
  
  signif_matrix = plot_data %>% filter(outcome == outcome_name) %>%
    dplyr::select(independent_variable, age_group, signif) %>%
    pivot_wider(names_from = age_group, 
                values_from = signif, 
                id_cols = independent_variable) %>%
    column_to_rownames(var = "independent_variable") %>%
    as.matrix()
  
  signif_matrix[is.na(signif_matrix)] <- "Not Significant"

  
  # Draw the heatmap with the significance annotations
  Heatmap(heatmap_matrix, 
          name = ifelse(outcome_type == "continuous", "Z-score difference\nper NPX", "Incidence ratio\nper NPX"),
          cluster_rows = TRUE, 
          cluster_columns = FALSE,
          clustering_distance_rows = "euclidean",
          clustering_method_rows = "ward.D",
          row_names_gp = gpar(fontsize = 9),
          na_col = "grey90",
          col = col_fun,
          
          
          # add asterisk to the heat map tiles that are significant
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(signif_matrix[i, j] == "Significant") {
              grid.points(x, y, pch = 8, gp = gpar(fontsize = 7))
            } else {
              grid.text("", x, y)
            }
          }
          
  ) 
  
}

pdf(paste0(figure_path, "plot-med-outcome-heatmap-LAZ.pdf"), width=5, height=10)
make_heatmap(plot_data = zscore_plot_data, outcome_name = "Length-for-age Z", outcome_type = "continuous")
dev.off()

pdf(paste0(figure_path, "plot-med-outcome-heatmap-WLZ.pdf"), width=5, height=10)
make_heatmap(plot_data = zscore_plot_data, outcome_name = "Weight-for-length Z", outcome_type = "continuous")
dev.off()

pdf(paste0(figure_path, "plot-med-outcome-heatmap-stunt.pdf"), width=5, height=7.5)
make_heatmap(plot_data = binary_plot_data, outcome_name = "Moderate-to-severe stunting", outcome_type = "binary")
dev.off()

pdf(paste0(figure_path, "plot-med-outcome-heatmap-waste.pdf"), width=5, height=7.5)
make_heatmap(plot_data = binary_plot_data, outcome_name = "Moderate-to-severe wasting", outcome_type = "binary")
dev.off()
