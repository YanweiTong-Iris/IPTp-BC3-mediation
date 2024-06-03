# !/bin/bash

R CMD BATCH 1-aim1-stratified.R &
R CMD BATCH 2a-intervention-mediator-analysis.R &
R CMD BATCH 2b-mediator-outcome-analysis.R &
R CMD BATCH 3-single-mediator-analysis.R &
R CMD BATCH 4-sensitivity-analysis.R &
R CMD BATCH 5-incidence-within-age.R
