################################################################
# IPTp and child growth
# Extract baseline demographics
# Last updated: July 19, 2024
################################################################

dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))
df_maternal = read_dta(paste0(box_path, "Maternal/", "BC-3 mothers expanded database FINAL.dta"))
df_maternal_indi = read_dta(paste0(box_path, "Maternal/", "Maternal BC3 individual level database FINAL.dta"))
zscore_monthly = readRDS(paste0(data_path,"analysis_data_zscore_monthly.RDS")) %>% mutate(age = as.numeric(age))
zscore_quarterly = readRDS(paste0(data_path,"analysis_data_zscore_quarterly.RDS"))

df_maternal_633 = df_maternal_indi %>% filter(id %in% zscore_quarterly$motherid) %>%
  mutate(anyHP = as.numeric(anyHP), placentalLAMPdich = as.numeric(placentalLAMPdich)) %>%
  mutate(GAenroll_cat = ifelse(GAenroll < 16, 1, 2))
df_maternal_633$placentalmal <- df_maternal_633$anyHP
df_maternal_633$placentalmal[df_maternal_633$placentalLAMPdich == 1 |
                               df_maternal_633$placentalBSdich == 1 |
                               df_maternal_633$anyHP == 1] <- 1
df_maternal_633$placentalmal[is.na(df_maternal_633$anyHP) &
                  is.na(df_maternal_633$placentalBSdich) &
                  is.na(df_maternal_633$placentalLAMPdich)] <- NA

zscore_6month = zscore_monthly %>% filter(age > 6) %>% mutate(age = round(as.numeric(age), 2)) %>%
  filter(if_all(5:7, ~ !is.na(.)))
df_maternal_6 = df_maternal_633 %>% filter(id %in% zscore_6month$motherid)


zscore_12month = zscore_monthly %>% filter(age > 11) %>% mutate(age = round(as.numeric(age), 2)) %>%
  filter(if_all(5:7, ~ !is.na(.)))
df_maternal_12 = df_maternal_633 %>% filter(id %in% zscore_12month$motherid)


################################################################
# placental malaria: histopath vs LAMP
################################################################

placentalmal_corr <- cor(df_maternal_633$anyHP, df_maternal_633$placentalLAMPdich, use = "complete.obs")
print(placentalmal_corr)

df_maternal_633 %>%
  filter(placentalmal == 1) %>%
  nrow()

df_maternal_633 %>%
  filter(anyHP == 1) %>%
  nrow()

df_maternal_633 %>%
  filter(placentalLAMPdich == 1) %>%
  nrow()

df_maternal_633 %>%
  filter(anyHP == 1, placentalLAMPdich == 1) %>%
  nrow()



################################################################
# baseline characteristics
################################################################

table(df_maternal_633$Txarm, df_maternal_633$anyHP)
table(df_maternal_633$Txarm, df_maternal_633$placentalmal)
table(df_maternal_633$Txarm, df_maternal_633$Gravidcat)
table(df_maternal_633$Txarm, df_maternal_633$GAenroll_cat)
table(df_maternal_633$Txarm, df_maternal_633$preterm)

summary_median_IQR <- function(x) {
  med <- round(median(x, na.rm = TRUE), 2)
  q1 <- round(quantile(x, 0.25, na.rm = TRUE), 2)
  q3 <- round(quantile(x, 0.75, na.rm = TRUE), 2)
  return(paste0(med, " (", q1, " – ", q3, ")"))
} 


summary_stats <- df_maternal_633 %>%
  group_by(Txarm) %>%
  summarise(across(c(HBenroll, enrollage, GAenroll), summary_median_IQR))


print(summary_stats)


################################################################
# 12-month completer demographics
################################################################

table(df_maternal_12$Txarm)
table(df_maternal_12$Txarm, df_maternal_12$anyHP)
table(df_maternal_12$Txarm, df_maternal_12$placentalmal)
table(df_maternal_12$Txarm, df_maternal_12$Gravidcat)
table(df_maternal_12$Txarm, df_maternal_12$GAenroll_cat)
table(df_maternal_12$Txarm, df_maternal_12$preterm)
table(df_maternal_12$Txarm, df_maternal_12$gender)
table(df_maternal_12$Txarm, df_maternal_12$educdich)
table(df_maternal_12$Txarm, df_maternal_12$APdichenroll)
table(df_maternal_12$Txarm, df_maternal_12$APdichenroll)


######Categorical#########

characteristics <- c("anyHP", "placentalmal", "Gravidcat", "GAenroll_cat", "preterm", "LBW","gender")

# Function to create a summary table for each characteristic
create_summary_table <- function(df, column) {
  df %>%
    group_by(Txarm, !!sym(column)) %>%
    count() %>%
    ungroup() %>%
    mutate(Characteristic = column) %>%
    rename(Category = !!sym(column))
}


summary_tables <- lapply(characteristics, create_summary_table, df = df_maternal_12)

combined_summary <- bind_rows(summary_tables) %>% na.omit()


######Continuous#########

summary_stats_12month <- df_maternal_12 %>%
  group_by(Txarm) %>%
  summarise(across(c(HBenroll, enrollage, GAenroll, GAdelivery), summary_median_IQR))

print(summary_stats_12month)



################################################################
# 6 months completer vs non-completer demographics
################################################################

table(df_maternal_6$placentalmal)
table(df_maternal_6$Gravidcat)
table(df_maternal_6$preterm)
table(df_maternal_6$gender)
table(df_maternal_6$APdichenroll)
table(df_maternal_6$LBW)

summary_stats_completer_6 <- df_maternal_6 %>%
  summarise(across(c(HBenroll, enrollage, GAenroll, GAdelivery), summary_median_IQR))

print(summary_stats_completer_6)


df_maternal_non_completer_6 = df_maternal_633 %>% filter(!id %in% zscore_6month$motherid)

table(df_maternal_non_completer_6$placentalmal)
table(df_maternal_non_completer_6$Gravidcat)
table(df_maternal_non_completer_6$preterm)
table(df_maternal_non_completer_6$gender)
table(df_maternal_non_completer_6$APdichenroll)
table(df_maternal_non_completer_6$LBW)

summary_stats_non_completer_6 <- df_maternal_non_completer_6 %>%
  summarise(across(c(HBenroll, enrollage, GAenroll, GAdelivery), summary_median_IQR))

print(summary_stats_non_completer_6)





################################################################
# 12 months completer vs non-completer demographics
################################################################

table(df_maternal_12$placentalmal)
table(df_maternal_12$Gravidcat)
table(df_maternal_12$GAenroll_cat)
table(df_maternal_12$preterm)
table(df_maternal_12$gender)
table(df_maternal_12$APdichenroll)
table(df_maternal_12$LBW)

summary_stats_completer <- df_maternal_12 %>%
  summarise(across(c(HBenroll, enrollage, GAenroll, GAdelivery), summary_median_IQR))

print(summary_stats_completer)


df_maternal_non_completer = df_maternal_633 %>% filter(!id %in% zscore_12month$motherid)

table(df_maternal_non_completer$placentalmal)
table(df_maternal_non_completer$Gravidcat)
table(df_maternal_non_completer$preterm)
table(df_maternal_non_completer$gender)
table(df_maternal_non_completer$APdichenroll)
table(df_maternal_non_completer$LBW)

summary_stats_non_completer <- df_maternal_non_completer %>%
  summarise(across(c(HBenroll, enrollage, GAenroll, GAdelivery), summary_median_IQR))

print(summary_stats_non_completer)

