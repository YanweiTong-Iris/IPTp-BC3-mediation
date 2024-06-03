library(anthro)
library(dplyr)
library(assertthat)
library(ggplot2)
library(mgcv)
library(tidyr)
library(mediation)
library(metafor)
library(msm)
library(tibble)
library(haven)
library(purrr)
library(readxl)
library(writexl)
library(arsenal)
library(expss)
library(geepack)
library(binom)
library(Hmisc)
library(OlinkAnalyze)
library(factoextra)
library(pls)
library(RColorBrewer)
library(forcats)
library(stringr)
library(ComplexHeatmap)
library(gridExtra)
library(forestplot)
library(knitr)
library(kableExtra)
library(magick)
library(webshot)
library(multcomp)
library(doParallel)
library(doRNG)
library(lubridate)

# FILL THE DATA DIRECTORY
if(Sys.getenv("USER") == "yanweitong"){ 
  box_path = "/Users/yanweitong/Library/CloudStorage/Box-Box/BC3 databases_Jade/"
}

data_path = paste0(box_path,"processed_data/")

#Load utility functions
util_functions = list.files(paste0(here::here(), "/0-base-functions/"), pattern = "*.R")
for (util in util_functions) {
  source(paste0(here::here(), "/0-base-functions/", util))
} 


figure_path = paste0(here::here(), "/5-figures/")
tables_path = paste0(here::here(), "/6-tables/")
results_path = paste0(here::here(), "/7-results/")

