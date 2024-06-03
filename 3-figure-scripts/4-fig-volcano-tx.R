################################################################
# IPTp and child growth
# Volcano plot figure
################################################################
set.seed(2023)
rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))

Olink_long = read_NPX(paste0(box_path, "Samples/","05-18-23 Olink_Inflammation_BC3 Maternal Speciemen_3 plates_Report.xlsx"))
Olink_long = Olink_long %>% mutate(motherID = substr(SampleID, start = 1, stop = 5))

LOD_check = Olink_long %>%
  mutate(below_LOD = NPX < LOD) %>%
  group_by(Assay) %>%
  summarise(percentage_below_LOD = mean(mean(below_LOD, na.rm = TRUE)) * 100) %>%
  filter(percentage_below_LOD < 50)

Assay_keep = setdiff(unique(LOD_check$Assay), c("Det Ctrl", "Ext Ctrl", "Inc Ctrl 1", "Inc Ctrl 2"))
Olink_long = Olink_long %>% filter(Assay %in% Assay_keep)


df_maternal_selected = dfz %>% dplyr::select(motherid, Txarm, 
                                             sex, enrollage, maternal_agecat, 
                                             Gravidity, preterm, birthweight,
                                             LBW, anyHP, SGA, anyHP, placentalmal, 
                                             placentalBSdich, placentalLAMPdich, 
                                             BSdichenroll, Graviddich, 
                                             APdichenroll, GAenroll, educdich, 
                                             enrollage_binary, enrollage_binary_30, 
                                             wealth_binary, wealthcat) %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP"))) %>%
  mutate(sex = factor(sex, levels = c("F", "M"))) %>%
  mutate(maternal_agecat = factor(
    maternal_agecat,
    levels = c("less than 20", "20-24", "25-29", "30+"))) %>%
  mutate(gravidity_cat = ifelse(Gravidity == 1, "1", ">1")) %>%
  mutate(gravidity_cat = factor(gravidity_cat, levels = c("1", ">1"))) %>%
  mutate(preterm = factor(preterm, levels = c("0", "1"))) %>%
  mutate(LBW = as.factor(LBW)) %>%
  mutate(anyHP = factor(ifelse(anyHP == "Yes", 1, 0))) %>%
  group_by(motherid) %>%
  slice(1) %>%
  ungroup()

Olink_combined <-
  merge(
    Olink_long,
    df_maternal_selected,
    by.x = "motherID",
    by.y = "motherid",
    all.x = TRUE
  ) %>%
  arrange(OlinkID, Index) %>%
  mutate(
    motherID = as.numeric(motherID),
    Txarm = factor(Txarm),
    preterm = factor(preterm),
    SGA = factor(SGA),
    anyHP = factor(anyHP),
    LBW = factor(LBW),
  ) %>% 
  mutate(Txarm = fct_rev(Txarm))

all_OlinkID = unique(Olink_combined$OlinkID)

Olink_mother = na.omit(subset(Olink_combined,!(
  Olink_combined$motherID %in% c("CONS-", "IPC-1", "IPC-2", "IPC-3", "NEG-1", "NEG-2", "NEG-3")
)))


# general risk factors
Olink_mother = subset(Olink_mother,!is.na(NPX)) %>%
  mutate(APdichenroll = factor(APdichenroll) , educdich = factor(educdich)) %>%
  mutate(enrollage_binary = factor(ifelse(enrollage< median(enrollage), 0 , 1))) %>%
  mutate(enrollage_binary_30 = factor(ifelse(enrollage< 30, 0 , 1))) %>%
  mutate(wealth_binary = factor(ifelse(wealthcat==3, 0, 1)))

# customize plot
olink_volcano_plot <- function (p.val_tbl, x_lab = "Estimate", olinkid_list = NULL, ...) {
  if (length(list(...)) > 0) {
    ellipsis_variables <- names(list(...))
    if (length(ellipsis_variables) == 1) {
      if (!(ellipsis_variables == "coloroption")) {
        stop(paste0("The ... option only takes the coloroption argument. ... currently contains the variable ", 
                    ellipsis_variables, "."))
      }
    }
    else {
      stop(paste0("The ... option only takes one argument. ... currently contains the variables ", 
                  paste(ellipsis_variables, collapse = ", "), "."))
    }
  }
  if (is.null(olinkid_list)) {
    olinkid_list <- p.val_tbl %>% dplyr::filter(Threshold == 
                                                  "Significant") %>% dplyr::pull(OlinkID)
  }
  volcano_plot <- p.val_tbl %>% ggplot2::ggplot(ggplot2::aes(x = estimate, 
                                                             y = -log10(p.value), 
                                                             color = Threshold)) + 
    ggplot2::geom_point() + 
    ggplot2::labs(x = x_lab, y = "-log10(p-value)") + 
    ggrepel::geom_label_repel(data = subset(p.val_tbl, OlinkID %in% olinkid_list), 
                              ggplot2::aes(label = Assay), box.padding = 1, show.legend = FALSE) + 
    ggplot2::geom_hline(yintercept = -log10(0.05),  linetype = "dotted") + OlinkAnalyze::set_plot_theme() + 
    OlinkAnalyze::olink_color_discrete(...) + 
    scale_color_manual(values = c("#909190", "#9e161d","#164c9e"))  +
    theme(
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(color = "black"),
      legend.position = "bottom"
    ) +
    xlab("log2(fold change)\nDP vs. SP")  
    
  return(volcano_plot)
}

# make volcano plot ---------------------------------------------------------
# label the top 10 most significant proteins
factor = "Txarm"

ttest_results <- olink_ttest(df = Olink_mother,
                               variable = "Txarm",
                               alternative = 'two.sided')
ttest_results = cbind(Maternal_factor = factor, ttest_results) %>% 
  mutate(Threshold = case_when(
    estimate > 0 & p.value < 0.05  ~"Up-regulated in DP; p<0.05", 
    estimate < 0 & p.value < 0.05  ~"Up-regulated in SP; p<0.05", 
    p.value >= 0.05 ~ "p>=0.05"
  )) %>% 
  mutate(adj_threshold = ifelse(Adjusted_pval < 0.05, "Significant", "Not significant"))

top_10_name <- ttest_results %>%
    slice_head(n = 15) %>%
    pull(OlinkID)

TNF_name <- ttest_results %>% filter(Assay == "DNER") %>% pull(OlinkID)
  
plot = olink_volcano_plot(p.val_tbl = ttest_results,
                            x_lab = factor,
                            olinkid_list = c(top_10_name, TNF_name)
                          #max.overlaps = 20
                          )
plot

ggsave(paste0(figure_path, "plot-volcano-tx.png"), plot, width = 7, height = 7, dpi = 300)

# number of Olink markers with p-value <0.05
ttest_results %>% filter(p.value < 0.05) %>% nrow()

# list of Olink proteins upregulated in DP
ttest_results %>% filter(estimate > 0 & p.value <= 0.05) %>% pull(Assay)

# list of Olink proteins upregulated in SP
ttest_results %>% filter(estimate < 0 & p.value <= 0.05) %>% pull(Assay)


