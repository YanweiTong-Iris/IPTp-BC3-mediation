################################################################
# IPTp and child growth
# Script for calculating linear growth velocity
# Last updated: Feburary 27, 2023
################################################################

rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

# load data -----------------------------------------

#df = read.csv(paste0(box_path,"/Pediatric/BC-3 childs all visit database FINAL_withIncidentNMF.csv"))

dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))

data = dfz %>% 
  dplyr::select(id, uniqueid, sex, age, gender, agecat, agemonthcat, weight, height, haz, whz) %>%
  mutate(gram_weight = weight * 1000)

id_list = unique(data$id)


# ---------------------------------------------------
# collapse the data to monthly by keeping the closest age measurement

indices <- 0:12

collapsed_df = list()

for (i in 1:length(id_list)) {
  indi_df = data[which(data$id == id_list[[i]]),]
  ages = indi_df$age
  closest_ages <- sapply(indices, function(x) {
    diff <- abs(ages - x)
    ifelse(min(diff) > 0.5, NA, ages[which.min(diff)])
  })
  collapsed_df[[i]] <-
    indi_df[which(indi_df$age %in% closest_ages),]
}

collapsed_df <- bind_rows(collapsed_df)

age_recat = round(collapsed_df$age)

collapsed_df$agemonth_recat <- age_recat


#----------------------------------------
# growth velocity function: targeted time interval * linear growth velocity per month (diff in growth divided by lapsed age in month)
# data: the data frame for one individual newborn
# interval: measurement interval in month
# yname: the outcome variable - length (height) or weight


growth_velocity = function(data, interval = 1, yname) {
  diff_age = diff(data$age, lag = interval)
  diff_y = diff(data[[yname]], lag = interval)
  LGV = interval * diff_y / diff_age
  return(LGV)
}


#----------------------------------------
# calculate and save linear growth velocities
# wgv1 = weight velocity based on absolute weight in g (1 month increment)
# wgv3 = weight velocity based on absolute weight in g (3 month increment)
# wlz_gv1 = weight velocity based on weight-for-length Z-score (1 month increment)
# wlz_gv3 = weight velocity based on weight-for-length Z-score (3 month increment)
# lgv2 = length velocity based on absolute height in cm (2 month increment)
# lgv3 = length velocity based on absolute height in cm (3 month increment)
# laz_gv2 = length velocity based on length-for-age Z-score (2 month increment)
# laz_gv3 = length velocity based on length-for-age Z-score (3 month increment)


#----------------------------------------
#1-month increment
collapsed_df <-
  collapsed_df %>% mutate(month1_interval = paste0(
    as.character(agemonth_recat),
    "-",
    as.character(agemonth_recat + 1) ,
    "mo"
  ))
#View(collapsed_df)

age_diff_list = list()
weight_results = list()
wlz_results = list()
age_diff_cat_list = list()


for (i in 1:length(id_list)) {
  indi_df = collapsed_df[which(collapsed_df$id == id_list[[i]]), ]
  length_NA = length(indi_df$age) - length(diff(indi_df$age, lag = 1))
  age_diff_list[[i]] <- c(diff(indi_df$age), rep(NA, length_NA))
  weight_results[[i]] <-
    c(
      growth_velocity(
        data = indi_df,
        yname = "gram_weight",
        interval = 1
      ),
      rep(NA, length_NA)
    )
  wlz_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "whz",
      interval = 1
    ),
    rep(NA, length_NA))
  if (length(weight_results[[i]]) >= 2) {
    for (j in 2:length(weight_results[[i]]) - 1) {
      if (!is.na(age_diff_list[[i]][j]) &&
          (age_diff_list[[i]][j] < 0.5 ||
           age_diff_list[[i]][j] > 1.5)) {
        weight_results[[i]][j] <- NA
        wlz_results[[i]][j] <- NA
      }
    }
  }
}


age_diff1 = unlist(age_diff_list)
wgv1 = unlist(weight_results)
wlz_gv1 = unlist(wlz_results)

final_df1 <- cbind(collapsed_df, age_diff1, wgv1, wlz_gv1)
final_df1$month1_interval <-
  factor(
    final_df1$month1_interval,
    levels = c(
      "0-1mo",
      "1-2mo",
      "2-3mo",
      "3-4mo",
      "4-5mo",
      "5-6mo",
      "6-7mo",
      "7-8mo",
      "8-9mo",
      "9-10mo",
      "10-11mo",
      "11-12mo",
      "12-13mo"
    )
  )
final_df1_save <- final_df1 %>% filter(month1_interval != "12-13mo")
#View(final_df1_save)

d1_save <- right_join(dfz, dplyr::select(final_df1_save, c(2, 12:17)), by = "uniqueid")
#View(d1_save)



#----------------------------------------
#2-month increment
final_df1 <-
  final_df1 %>% mutate(month2_interval = paste0(
    as.character(agemonth_recat),
    "-",
    as.character(agemonth_recat + 2) ,
    "mo"
  ))

age_diff_list = list()
length_results = list()
laz_results = list()


for (i in 1:length(id_list)) {
  indi_df = collapsed_df[which(collapsed_df$id == id_list[[i]]),]
  length_NA = length(indi_df$age) - length(diff(indi_df$age, lag = 2))
  age_diff_list[[i]] <-
    c(diff(indi_df$age, lag = 2), rep(NA, length_NA))
  length_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "height",
      interval = 2
    ),
    rep(NA, length_NA))
  laz_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "haz",
      interval = 2
    ),
    rep(NA, length_NA))
  if (length(length_results[[i]]) >= 3) {
    for (j in 2:length(length_results[[i]]) - 1) {
      if (!is.na(age_diff_list[[i]][j]) &&
          (age_diff_list[[i]][j] < 1.5 ||
           age_diff_list[[i]][j] > 2.5)) {
        length_results[[i]][j] <- NA
        laz_results[[i]][j] <- NA
      }
    }
  }
}


age_diff2 = unlist(age_diff_list)
lgv2 = unlist(length_results)
laz_gv2 = unlist(laz_results)

final_df2 <- cbind(final_df1, age_diff2, lgv2, laz_gv2)
final_df2_save <-
  final_df2[final_df2$month2_interval %in% c("0-2mo", "2-4mo", "4-6mo", "6-8mo", "8-10mo", "10-12mo"), ]
final_df2_save$month2_interval <-
  factor(
    final_df2_save$month2_interval,
    levels = c("0-2mo", "2-4mo", "4-6mo", "6-8mo", "8-10mo", "10-12mo")
  )
#View(final_df2_save)

d2_save <-right_join(dfz, dplyr::select(final_df2_save, 2, 12, 13, 18:21), by = "uniqueid")
#View(d2_save)



#----------------------------------------
#3-month increment

final_df2 <-
  final_df2 %>% mutate(month3_interval = paste0(
    as.character(agemonth_recat),
    "-",
    as.character(agemonth_recat + 3) ,
    "mo"
  ))

age_diff_list = list()
weight_results = list()
length_results = list()
wlz_results = list()
laz_results = list()


for (i in 1:length(id_list)) {
  indi_df = collapsed_df[which(collapsed_df$id == id_list[[i]]), ]
  length_NA = length(indi_df$age) - length(diff(indi_df$age, lag = 3))
  age_diff_list[[i]] <-
    c(diff(indi_df$age, lag = 3), rep(NA, length_NA))
  weight_results[[i]] <-
    c(
      growth_velocity(
        data = indi_df,
        yname = "gram_weight",
        interval = 3
      ),
      rep(NA, length_NA)
    )
  wlz_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "whz",
      interval = 3
    ),
    rep(NA, length_NA))
  length_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "height",
      interval = 3
    ),
    rep(NA, length_NA))
  laz_results[[i]] <-
    c(growth_velocity(
      data = indi_df,
      yname = "haz",
      interval = 3
    ),
    rep(NA, length_NA))
  if (length(length_results[[i]]) >= 3) {
    for (j in 2:length(length_results[[i]]) - 1) {
      if (!is.na(age_diff_list[[i]][j]) &&
          (age_diff_list[[i]][j] < 2.5 ||
           age_diff_list[[i]][j] > 3.5)) {
        weight_results[[i]][j] <- NA
        wlz_results[[i]][j] <- NA
        length_results[[i]][j] <- NA
        laz_results[[i]][j] <- NA
      }
    }
  }
}


age_diff3 = unlist(age_diff_list)
wgv3 = unlist(weight_results)
wlz_gv3 = unlist(wlz_results)
lgv3 = unlist(length_results)
laz_gv3 = unlist(laz_results)

final_df3 <-
  cbind(final_df2, age_diff3, wgv3, wlz_gv3, lgv3, laz_gv3)
final_df3_save <-
  final_df3[final_df3$month3_interval %in% c("0-3mo", "3-6mo", "6-9mo", "9-12mo"), ]
final_df3_save$month3_interval <-
  factor(final_df3_save$month3_interval,
         levels = c("0-3mo", "3-6mo", "6-9mo", "9-12mo"))
#View(final_df3_save)

d3_save <-right_join(dfz, dplyr::select(final_df3_save, 2, 12, 13, 22:27), by = "uniqueid")
#View(d3_save)



#----------------------------------------
#save df
#----------------------------------------

saveRDS(d1_save, paste0(data_path, "weight_linear_growth_velocity_1month.RDS"))
saveRDS(d2_save, paste0(data_path, "length_linear_growth_velocity_2month.RDS"))
saveRDS(d3_save, paste0(data_path,"weight&length_linear_growth_velocity_3month.RDS"))

