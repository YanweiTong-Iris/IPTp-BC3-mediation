# !/bin/bash

R CMD BATCH 1_z_score.R &
R CMD BATCH 2-outcome-zscores.R &
R CMD BATCH 3-incidence-outcomes.R &
R CMD BATCH 4-linear-growth-velocity.R &
R CMD BATCH 5-non-malarial-outcomes.R &
R CMD BATCH 6-covariates.R &
R CMD BATCH 7-medication.R &
R CMD BATCH 8-analysis-data.R 
