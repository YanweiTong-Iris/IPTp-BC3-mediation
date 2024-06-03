################################################################
# IPTp and child growth
# Script for making Z-score density plots
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library('data.table')

#Loading the data
z.plot <- readRDS(paste0(data_path, "cleaned_zscores_data.RDS")) %>% 
  mutate(agem_birth = ifelse(agedays==0, 0, agemonthcat)) %>% 
#Recategorizing age into the following categories 
 mutate(agecat_birth = case_when(
  age == 0 ~ "Birth",
  age >0 & age <= 3 ~ "0-3",
  age >3 & age <= 6 ~ "4-6",
  age >6 & age <= 9 ~ "7-9",
  age >9 & age <= 12 ~ "10-12")) %>% 
  mutate(agecat_birth = factor(agecat_birth, levels = c(
    "Birth", "0-3", "4-6",
    "7-9", "10-12"
  )))

median(z.plot$haz[z.plot$agecat_birth=="Birth"], na.rm = T)
median(z.plot$whz[z.plot$agecat_birth=="Birth"], na.rm = T)


library(ggridges)

haz_plot <- ggplot(z.plot, aes(y=agecat_birth, x = haz, fill=stat(x)))+
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2,
                               bandwidth = 0.3) +
  scale_fill_gradient2(low = "#754306", mid = "white", high = "#096a6e") +
  scale_x_continuous(
    limits = c(-5.5, 5.5),
    labels = c(-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5),
    breaks = c(-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5)
  ) + 
  theme_minimal() +
  xlab("Z-score") + 
  coord_flip() + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("A) Length-for-age Z")

wlz_plot <- ggplot(z.plot, aes(y=agecat_birth, x = whz, fill=stat(x)))+
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2,
                               bandwidth = 0.3) +
  scale_fill_gradient2("",low = "#754306", mid = "white", high = "#096a6e") +
  scale_x_continuous(
    limits = c(-5.5, 5.5),
    labels = c(-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5),
    breaks = c(-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5)
  ) + 
  theme_minimal() +
  xlab("Z-score") + 
  coord_flip() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ggtitle("B) Weight-for-length Z")

plot <- grid.arrange(haz_plot, wlz_plot, ncol=2,
                     widths=c(4.3,5))

ggsave(plot, filename = paste0(figure_path, "plot-zscore_density_plot.png"), 
       width = 7, height = 2.5)