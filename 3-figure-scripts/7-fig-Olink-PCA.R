################################################################
# IPTp and child growth
# Olink PCA
# Last updated: July 20, 2024
################################################################

dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))
zscore_quarterly = readRDS(paste0(data_path,"analysis_data_zscore_quarterly.RDS"))
df_maternal_indi = read_dta(paste0(box_path, "Maternal/", "Maternal BC3 individual level database FINAL.dta"))
df_maternal_selected = dfz %>% 
  filter(motherid %in% zscore_quarterly$motherid) %>%
  dplyr::select(motherid, Txarm, placentalmal) %>%
  group_by(motherid) %>%
  slice(1) %>%
  ungroup()

Olink_long = read_NPX(paste0(box_path, "Samples/","05-18-23 Olink_Inflammation_BC3 Maternal Speciemen_3 plates_Report.xlsx")) %>% 
  mutate(motherID = substr(SampleID, start = 1, stop = 5)) %>%
  filter(!Assay %in% c("Det Ctrl", "Ext Ctrl", "Inc Ctrl 1", "Inc Ctrl 2"))

Olink_broad = read_excel(
  paste0(box_path, "Samples/","05-18-23 Olink_Inflammation_BC3 Maternal Speciemen_3 plates_Report_broad.xlsx")
)

LOD_check = Olink_long %>%
  mutate(below_LOD = NPX < LOD) %>%
  group_by(Assay) %>%
  summarise(percentage_below_LOD = mean(mean(below_LOD, na.rm = TRUE)) * 100) %>%
  filter(percentage_below_LOD < 50)

Olink_LOD = unique(LOD_check$Assay)

Olink_long = Olink_long %>% filter(Assay %in% Olink_LOD)
olink_pca_plot(Olink_long)

Olink_long_combined <-
  merge(
    Olink_long,
    df_maternal_selected,
    by.x = "motherID",
    by.y = "motherid",
    all.x = TRUE
  )



Olink_LOD = gsub("-", "_", Olink_LOD)
Olink_LOD = gsub("LAP TGF_beta_1", "LAP_TGF_beta_1", Olink_LOD)


Olink_broad = Olink_broad %>% mutate(motherID = as.numeric(substr(ID, start = 1, stop = 5))) %>%
  dplyr::select(c(motherID, all_of(Olink_LOD)))

Olink_broad_combined <-
  merge(
    Olink_broad,
    df_maternal_selected,
    by.x = "motherID",
    by.y = "motherid",
    all.x = TRUE
  )

Olink_broad_combined_clean <- Olink_broad_combined %>%
  filter_all(all_vars(!is.na(.))) %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP")))


PCA_result = prcomp(Olink_broad_combined_clean[2:68], center = TRUE, scale = TRUE)
summary(PCA_result)
biplot(PCA_result) 

screeplot(PCA_result, type = "l")


pca_plot = fviz_pca_ind(
    PCA_result,
    geom.ind = "point",
    pointshape = 21,
    pointsize = 2,
    fill.ind = Olink_broad_combined_clean[["Txarm"]],
    col.ind = "black",
    palette = "jco",
    addEllipses = TRUE,
    label = "var",
    col.var = "black",
    repel = TRUE,
    legend.title = "Intervention"
  ) +
    scale_fill_manual(values = c("SP" = "#164c9e", "DP" = "#9e161d")) +  # Manually set fill colors and legend title
    scale_color_manual(values = c("SP" = "#164c9e", "DP" = "#9e161d")) +  # Manually set border colors for ellipses and legend title
    ggtitle("2D PCA Plot by Intervention") +
    theme(plot.title = element_text(hjust = 0.5))
  
print(pca_plot)
ggsave(pca_plot, filename = paste0(figure_path, "plot-Olink_PCA.pdf"),
       width=6, height=4)





