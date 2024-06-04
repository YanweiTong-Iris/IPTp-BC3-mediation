# !/bin/bash

R CMD BATCH 1_z_score.R &
R CMD BATCH 2-outcome-zscores.R &
R CMD BATCH 3-incidence-outcomes.R &
R CMD BATCH 4-linear-growth-velocity.R &
R CMD BATCH 5-covariates.R &
R CMD BATCH 6-analysis-data.R 
