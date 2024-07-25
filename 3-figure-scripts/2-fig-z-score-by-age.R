################################################################
# IPTp and child growth
# Z score plots
# Last updated: April 19, 2024
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library('data.table')

#Loading the data
z.plot <- readRDS(paste0(data_path, "cleaned_zscores_data.RDS")) %>% 
  mutate(agem_birth = ifelse(agedays==0, 0, agemonthcat)) %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP")))


z.plot= z.plot %>% mutate(gravid_multi = ifelse(Gravidity ==1, "Primigravidae", "Multigravidae")) %>% 
  mutate(gravid_multi = factor(gravid_multi, levels = c("Primigravidae", "Multigravidae")))


meanz_age_gravid = z.plot %>% group_by(Txarm, gravid_multi, agem_birth) %>% 
  summarise(mean_haz = mean(haz, na.rm=T),
            mean_whz = mean(whz, na.rm=T),
            sd_haz = sd(haz, na.rm=T),
            sd_whz = sd(whz, na.rm=T),
            n=n()) %>% 
  mutate(se_haz = sd_haz/sqrt(n),
         se_whz = sd_whz/sqrt(n)) %>% 
  mutate(lb_haz = mean_haz - qnorm(0.975)*se_haz,
         ub_haz = mean_haz + qnorm(0.975)*se_haz,
         lb_whz = mean_whz - qnorm(0.975)*se_whz,
         ub_whz = mean_whz + qnorm(0.975)*se_whz)

haz_plot <- ggplot(meanz_age_gravid, aes(x = agem_birth, y = mean_haz)) + 
  geom_line(aes(col=Txarm)) +
  geom_point(aes(col=Txarm))+
  geom_ribbon(aes(ymin = lb_haz, ymax=ub_haz, fill=Txarm), alpha=0.3) +
  facet_grid(~gravid_multi, scales="free") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,12,1)) + 
  ylab("Mean length-for-age Z-score") +
  xlab("Age, months") +
  scale_color_manual(values = c("#164c9e", "#9e161d")) + 
  scale_fill_manual(values = c("#164c9e", "#9e161d")) +
  theme_minimal() +
  theme(strip.text = element_text(size=12),
        legend.position = "none",
        legend.title = element_blank()) +
  ggtitle("A)")

whz_plot <- ggplot(meanz_age_gravid, aes(x = agem_birth, y = mean_whz)) + 
  geom_line(aes(col=Txarm)) +
  geom_point(aes(col=Txarm))+
  geom_ribbon(aes(ymin = lb_whz, ymax=ub_whz, fill=Txarm), alpha=0.3) +
  facet_grid(~gravid_multi, scales="free") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,12,1)) + 
  ylab("Mean weight-for-length Z-score") +
  xlab("Age, months") +
  scale_color_manual(values = c("#164c9e", "#9e161d")) + 
  scale_fill_manual(values = c("#164c9e", "#9e161d")) +
  theme_minimal() +
  theme(strip.text = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle("B)")


plot_final <- grid.arrange(haz_plot, whz_plot, ncol=1,
                           heights=c(2,2.65))

ggsave(plot_final, filename = paste0(figure_path, "plot-zscore-tx-age.pdf"),
       width = 8, height = 6)


# *******************************************
# zscores with gestational age-corrected age (CA)
# *******************************************

CA_zscores <- readRDS(paste0(data_path, "cleaned_CA_zscores_data.RDS")) %>% 
  mutate(age_CA_monthcat = ifelse(age_CA_months >= 0, ceiling(age_CA_months), floor(age_CA_months)))%>%
  mutate(age_CA_monthcat = ifelse(age_CA_days ==0, 0, age_CA_monthcat)) %>%
  mutate(Txarm = factor(Txarm, levels = c("SP", "DP")))

CA_zscores = CA_zscores %>% mutate(gravid_multi = ifelse(Gravidity ==1, "Primigravidae", "Multigravidae")) %>% 
  mutate(gravid_multi = factor(gravid_multi, levels = c("Primigravidae", "Multigravidae"))) %>%
  filter(age_CA_days >= 0)


meanz_CA_gravid = CA_zscores %>% group_by(Txarm, gravid_multi, age_CA_monthcat) %>% 
  summarise(mean_haz = mean(CA_laz, na.rm = TRUE),
            mean_whz = mean(CA_wlz, na.rm=TRUE),
            sd_haz = sd(CA_laz, na.rm=T),
            sd_whz = sd(CA_wlz, na.rm=T),
            n=n()) %>% 
  mutate(se_haz = sd_haz/sqrt(n),
         se_whz = sd_whz/sqrt(n)) %>% 
  mutate(lb_haz = mean_haz - qnorm(0.975)*se_haz,
         ub_haz = mean_haz + qnorm(0.975)*se_haz,
         lb_whz = mean_whz - qnorm(0.975)*se_whz,
         ub_whz = mean_whz + qnorm(0.975)*se_whz)


CA_haz_plot <- ggplot(meanz_CA_gravid %>% filter(age_CA_monthcat >= 0), aes(x = age_CA_monthcat, y = mean_haz)) + 
  geom_line(aes(col=Txarm)) +
  geom_point(aes(col=Txarm))+
  geom_ribbon(aes(ymin = lb_haz, ymax=ub_haz, fill=Txarm), alpha=0.3) +
  facet_grid(~gravid_multi, scales="free") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,12,1)) + 
  ylab("Mean length-for-age Z-score") +
  xlab("Gestational age-corrected age (CA), months") +
  scale_color_manual(values = c("#164c9e", "#9e161d")) + 
  scale_fill_manual(values = c("#164c9e", "#9e161d")) + 
  theme_minimal() +
  theme(strip.text = element_text(size=12),
        legend.position = "none",
        legend.title = element_blank()) +
  ggtitle("A)")

CA_whz_plot <- ggplot(meanz_CA_gravid %>% filter(age_CA_monthcat >= 0), aes(x = age_CA_monthcat, y = mean_whz)) + 
  geom_line(aes(col=Txarm)) +
  geom_point(aes(col=Txarm))+
  geom_ribbon(aes(ymin = lb_whz, ymax=ub_whz, fill=Txarm), alpha=0.3) +
  facet_grid(~gravid_multi, scales="free") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,12,1)) + 
  ylab("Mean weight-for-length Z-score") +
  xlab("Gestational age-corrected age (CA), months") +
  scale_color_manual(values = c("#164c9e", "#9e161d")) + 
  scale_fill_manual(values = c("#164c9e", "#9e161d")) + 
  theme_minimal() +
  theme(strip.text = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle("B)")

CA_plot_final <- grid.arrange(CA_haz_plot, CA_whz_plot, ncol=1,
                              heights=c(2,2.65))

ggsave(CA_plot_final, filename = paste0(figure_path, "plot-zscore-tx-CA.pdf"),
       width = 8, height = 6)


