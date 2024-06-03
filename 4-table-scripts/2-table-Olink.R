################################################################
# IPTp and child growth
# Table of Olink inclusion
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

#-------------------------------------------
# All biomarkers
#-------------------------------------------
Olink_long_df = read_NPX(paste0(box_path, "Samples/","05-18-23 Olink_Inflammation_BC3 Maternal Speciemen_3 plates_Report.xlsx")) %>% 
  mutate(motherID = substr(SampleID, start = 1, stop = 5)) %>%
  filter(!Assay %in% c("Det Ctrl", "Ext Ctrl", "Inc Ctrl 1", "Inc Ctrl 2"))

Olink_all = unique(Olink_long_df$Assay)
length(Olink_all)
# n1 = 92

#-------------------------------------------
# Biomarkers with (NPX>LOD) > 50%
#-------------------------------------------
LOD_check = Olink_long_df %>%
  mutate(below_LOD = NPX < LOD) %>%
  group_by(Assay) %>%
  summarise(percentage_below_LOD = mean(mean(below_LOD, na.rm = TRUE)) * 100) %>%
  filter(percentage_below_LOD < 50)

Olink_LOD = unique(LOD_check$Assay)
length(Olink_LOD)
# n2 = 67


#-------------------------------------------
# Biomarkers with significant IM results
#-------------------------------------------
other_mediators = c("anemia_28binary", "antibacterial_binary", "betalactam_binary",
                          "placentalmal", "preterm", "LBW", "gestational_weightchange", 
                    "GWC_Z","birthweight", "birthlength", "birthweight_kg")

IM_results = readRDS(paste0(results_path,"IM-MO-stratified/intervention_mediator_main_results_stratified.RDS"))
Olink_IM = setdiff(IM_results$mediator, other_mediators)
Olink_IM <- gsub("_", "-", Olink_IM)
length(Olink_IM)
# n3 = 17


#-------------------------------------------
# Biomarkers with significant MO results
#-------------------------------------------
MO_results = readRDS(paste0(results_path,"IM-MO-stratified/mediator_outcome_zscore_results_3mo_stratified.RDS"))
Olink_MO = setdiff(MO_results$independent_variable, other_mediators)
Olink_MO <- gsub("_", "-", Olink_MO)
Olink_MO <- gsub("LAP-TGF-beta-1", "LAP TGF-beta-1", Olink_MO)
Olink_MO <- gsub("F4E-BP1", "4E-BP1", Olink_MO)
length(Olink_MO)
# n4 = 58


#-------------------------------------------
# Biomarkers with both significant MO and IM results
#-------------------------------------------
Olink_IM_MO = intersect(Olink_IM, Olink_MO)
length(Olink_IM_MO)
# n5= 15


#-------------------------------------------
# Biomarkers included in the mediation analyses
#-------------------------------------------
mediation_result = readRDS(paste0(results_path,"aim2-stratified/aim2_single_mediator_zscore_results_3mo.RDS"))
Olink_mediation = setdiff(mediation_result$mediator, other_mediators)
Olink_mediation <- gsub("_", "-", Olink_mediation)
length(Olink_mediation)
# n6 = 15

#-------------------------------------------
# make the table
#-------------------------------------------
table_Olink <- data.frame(
  "Passed LOD check" = rep("", length(Olink_all)),
  "Significant intervention-mediator results" = rep("", length(Olink_all)),
  "Significant mediator-outcome results" = rep("", length(Olink_all)),
  "Included in the mediation analysis" = rep("", length(Olink_all)), 
  row.names = Olink_all,
  check.names = FALSE  
)

# Fill the table with check marks (√) where appropriate
table_Olink["Passed LOD check"][row.names(table_Olink) %in% Olink_LOD,] <- "√"
table_Olink["Significant intervention-mediator results"][row.names(table_Olink) %in% Olink_IM,] <- "√"
table_Olink["Significant mediator-outcome results"][row.names(table_Olink) %in% Olink_MO,] <- "√"
table_Olink["Included in the mediation analysis"][row.names(table_Olink) %in% Olink_mediation,] <- "√"

# View the table
# View(table_Olink)
write.csv(table_Olink,  paste0(tables_path,"table_Olink_inclusion.csv"))

