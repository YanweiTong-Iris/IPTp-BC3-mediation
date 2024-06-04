########################################
#Coding covariates needed for analysis
#Creates variables for wet season and 
# month of measurement 
########################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

#data <- readRDS(paste0(data_path, "cleaned_zscores_agecat_birth.RDS"))
dfz = readRDS(paste0(box_path, "/processed_data/outcome_zscores_data.RDS"))
data = read.csv(paste0(box_path, "Pediatric/BC-3 childs all visit database FINAL_withIncidentNMF.csv"))

data = data %>%
  filter(id %in% dfz$id)

# taking only the months from the date string
data$month_string <- gsub("[^A-Za-z]+", "", data$date) 

# Create new variable for wet_season using month_string
data$wet_season <- ifelse(data$month_string %in% c("mar", "apr", "may", "oct", "nov", "dec"), 1, 0)


# Create new variable for season of malaria infection 
data$incidentmalaria_season <- 0

# Replace incidentmalaria_season values based on incidentmalaria and wet_season
data$incidentmalaria_season[data$incidentmalaria == 1 & data$wet_season == 1] <- 1
data$incidentmalaria_season[data$incidentmalaria == 1 & data$wet_season == 0] <- 2

data = data %>% mutate(incidentmalaria_season = case_when(
  incidentmalaria_season == 0 ~ "no incidence",
  incidentmalaria_season == 1 ~ "incidence in wet season",
  incidentmalaria_season == 2 ~ "incidence in dry season",
))
# check props of incident malaria in each season 
#prop.table(table(data$incidentmalaria_season))

data_diff = read_dta(paste0(box_path, "Maternal/", "BC-3 mothers expanded database FINAL.dta"))
data_diff <- data_diff %>%
  group_by(id) %>%
  mutate(weight_20wks = mean(weight[gestwks == 20], na.rm = TRUE), 
         #calculating mean weight at 20wks since there are multiple observations at 20 wks
         weight_diff = weight_20wks - last(weight)) %>% 
  #difference between mean weight at 20 weeks and last weight measured for each id  
  filter(!is.na(weight_diff)) %>% 
  ungroup()

#Collapsing the dataset before merging it with the pediatric dataset
data_diff_collapsed <- distinct(data_diff, id, weight_diff)

#Merging the paediatric and maternal dataset (all.x = TRUE keeps all rows in the pediatric dataset)
comb_data <- merge(data, data_diff_collapsed, by.x = "motherid", by.y = "id", all.x=TRUE) %>% 
  dplyr::select(id, date, weight_diff, wet_season, incidentmalaria_season) #selecting the relevant variables

analysis_data <- readRDS((paste0(data_path, "outcome_zscores_data.RDS")))  

cov_data <- merge(analysis_data, comb_data, by = c("id", "date"), all.x = TRUE)

saveRDS(cov_data, paste0(data_path, "outcome_zscores_covariates.RDS"))


