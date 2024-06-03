################################################################
# IPTp and child growth
# Script for calculating Zscores using anthro_zscores package
################################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# load data -----------------------------------------

d1 = read_dta(paste0(box_path, "Maternal/", "Maternal BC3 individual level database FINAL.dta"))  
d = read.csv(paste0(box_path, "Pediatric/BC-3 childs all visit database FINAL_withIncidentNMF.csv"))
colnames(d)

# exclude part of indi maternal data ----------------
d1$placentalmal <- d1$anyHP
d1$placentalmal[d1$placentalLAMPdich == 1 |
                  d1$placentalBSdich == 1 |
                  d1$anyHP == 1] <- 1
d1$placentalmal[is.na(d1$anyHP) &
                  is.na(d1$placentalBSdich) &
                  is.na(d1$placentalLAMPdich)] <- NA

# exclude stillbirth, twin, spont abortion, no histopath
d1 = d1 %>% mutate(include_flag = ifelse(is.na(placentalmal), 0, 1)) %>% 
  mutate(include_flag = ifelse(withdrawalanalysis %in% c(2,3), 0, include_flag)) %>%
  mutate(stillbirthdich = as.factor(stillbirthdich)) %>%
  mutate(include_flag = ifelse(!is.na(stillbirthdich) & stillbirthdich == 1, 0, include_flag))

d1_include = d1[d1$include_flag == 1,]
d1_include = d1_include %>% 
  dplyr::select(id, placentalmal, placentalBSdich, placentalLAMPdich, BSdichenroll)

#Joining the child and maternal datasets (633 mother-infant pairs left)
d = d %>% inner_join(d1_include, by = c("motherid" = "id"))


#sum(is.na(d$agecat))   # 0, I checked that there are no NA values for agecat 


# add GA corrected ages and zscores
d = d %>% mutate(agedays = age* 30.4375,
                 DOB = lubridate::dmy(d$dob),
                 Date = lubridate::dmy(d$date)) %>% 
  mutate(age_exact_days = as.numeric(difftime(Date, DOB, units="days"))) %>%
  mutate(age_exact_wk = age_exact_days/7) %>%
  #mutate(age_exact_month = as.numeric(difftime(Date, DOB, units="months"))) %>%
  mutate(age_GA_wk = age_exact_wk - 40 + GAdelivery) %>%
  mutate(age_GA_days = age_GA_wk * 7) %>%
  mutate(age_GA_months = age_GA_days/30.4375) %>%
  mutate(age_CA_days = ifelse(GAdelivery< 37, age_GA_days, age_exact_days),
         age_CA_wk = ifelse(GAdelivery< 37, age_GA_wk, age_exact_wk),
         age_CA_months = ifelse(GAdelivery< 37, age_GA_months, age))
  #mutate(age_GA_month = lubridate::interval(DOB, Date-lubridate::weeks(40-GAdelivery)) %/% months(1.00))
  

# creating WHO Z scores
d2 = d %>% mutate(sex = ifelse(gender=="Female", "F", "M")) 

#anthro_zscores is a package that calculates Z scores for outcomes related to stuning, wasting and underweight 
zscores <- anthro_zscores(
  sex = d2$sex,
  age = d2$agedays,
  is_age_in_month = FALSE, 
  weight = d2$weight,
  lenhei = d2$height,
  measure = "L",
  armc = NA_real_,
  triskin = NA_real_,
  subskin = NA_real_,
  oedema = "n"
) %>%  dplyr::select(zlen, flen, zwei, fwei, zwfl, fwfl)

table(zscores$flen)
table(zscores$fwei)
table(zscores$fwfl)


d3 = d2 %>% filter(age_CA_days >= 0)

CA_zscores <- anthro_zscores(
  sex = d3$sex,
  age = d3$age_CA_days,
  is_age_in_month = FALSE, 
  weight = d3$weight,
  lenhei = d3$height,
  measure = "L",
  armc = NA_real_,
  triskin = NA_real_,
  subskin = NA_real_,
  oedema = "n"
) %>%  dplyr::select(zlen, flen, zwei, fwei, zwfl, fwfl)

table(CA_zscores$flen)
table(CA_zscores$fwei)
table(CA_zscores$fwfl)


#Dropping zscores outside of standard range 
zscores <- zscores %>% mutate(
  zlen = ifelse(flen==1, NA, zlen),
  zwei = ifelse(fwei==1, NA, zwei),
  zwfl = ifelse(fwfl==1, NA, zwfl)
)

CA_zscores <- CA_zscores %>% mutate(
  zlen = ifelse(flen==1, NA, zlen),
  zwei = ifelse(fwei==1, NA, zwei),
  zwfl = ifelse(fwfl==1, NA, zwfl)
)

#abs(zlen) >6

summary(zscores$zlen) # max 230.63 before dropping z scores outside of standard range in above step
summary(zscores$zwei) # 75th percentile still stays the same 
summary(zscores$zwfl)

#Renaming to height for age, weight for height and weight for age z scores 
zscores <- zscores %>% rename(
  haz = zlen,
  whz = zwfl,
  waz = zwei
) %>% dplyr::select(haz, whz, waz) #use these for the z score plots 

CA_zscores <- CA_zscores %>% rename(
  CA_laz = zlen,
  CA_wlz = zwfl,
  CA_waz = zwei
) %>% dplyr::select(CA_laz, CA_wlz, CA_waz) 

assert_that(nrow(d2)==nrow(zscores)) 
assert_that(nrow(d3) == nrow(CA_zscores))

zscores <- bind_cols(d2, zscores) 
CA_zscores <- bind_cols(d3, CA_zscores)
CA_zscores_full <- left_join(d2, CA_zscores[c(1,5,146:148)], by = c("id", "age"))

saveRDS(zscores, paste0(data_path, "cleaned_zscores_data.RDS"))
saveRDS(CA_zscores_full, paste0(data_path, "cleaned_CA_zscores_data.RDS"))
saveRDS(d1, paste0(data_path, "Maternal_BC3_individual_level_database_Roh.RDS"))

