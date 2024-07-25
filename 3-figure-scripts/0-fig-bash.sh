# !/bin/bash

R CMD BATCH 1a-fig-zscore-distributions.R &
R CMD BATCH 1b-fig-incidence-quarterly.R &
R CMD BATCH 2-fig-z-score-by-age.R &
R CMD BATCH 3-fig-intervention-mediator.R &
R CMD BATCH 4-fig-volcano-tx.R &
R CMD BATCH 5a-fig-mediator-outcome.R
R CMD BATCH 5b-fig-mediator-outcome-Olink.R &
R CMD BATCH 6a-fig-mediation.R &
R CMD BATCH 6b-fig-mediation-Olink-heatmap.R &
R CMD BATCH 6c-fig-mediation-Olink.R &
R CMD BATCH 7-fig-Olink-PCA.R
