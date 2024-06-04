# Pathways through which intermittent preventive treatment for malaria in pregnancy influences child growth faltering: a mediation analysis 

## Overview
Intermittent Preventive Treatment in Pregnancy (IPTp) with sulfadoxine-pyrimethamine (SP) is recommended by the WHO for regions with moderate-to-high malaria transmission. While SP is effective in reducing neonatal mortality and low birthweight, its efficacy has diminished in some areas of sub-Saharan Africa due to rising parasite resistance. Although IPTp with dihydroartemisinin-piperaquine (IPTp-DP) has demonstrated superior efficacy in reducing maternal malaria, its impact on birth outcomes has not significantly surpassed that of SP. The ultimate goal of IPTp extends beyond enhancing birth outcomes to include benefits during infancy and later stages. Yet, the effects of SP vs. DP in relation to infant growth post-birth and the underlying mechanisms remain unknown. Prior studies also found that different IPTp regimens worked through different pathways, with DP influencing birth outcomes by reducing placental malaria and SP influencing them through non-malarial pathways such as maternal weight gain. Here, we re-analyzed data from of a randomized trial in Uganda to explore the impacts of the two different IPTp regimens on infant growth and to understand potential mechanisms of impacts on infant growth.

## Additional Information
A pre-specified analysis plan is available on [Open Science Framework](https://osf.io/f8wy4/).

To run the scripts, the data directory for the user must be changed in **`0-config.R`**. This will allow for replication of study findings using scripts in this repository. Similar directory statement changes will be needed wherever output files are saved down (e.g., raw estimates, figures). To run all scripts required to reproduce all analysis datasets and results in the manuscript, run the bash script **`0-run-project`**.

## Directory structure
**`0-config.R` :** configuration file that sets data directories, sources base functions, and loads required libraries

**`0-base-functions` :** folder containing R scripts for general functions used across the analysis

**`1-dm` :** folder containing data management scripts. To rerun all scripts in this subdirectory, run the bash script `0-dm-bash.sh`.  
* `1-z-score.R` : calculate infant growth Z-scores using `anthro_zscores` package  
* `2-zscore-outcomes.R` : generate a cleaned Z-score outcome data set  
* `3-incidence-outcomes.R` : generate cleaned incidence outcome data sets   
* `4-linear-growth-velocity.R` : calculate linear growth velocity for both weight and length   
* `5-non-malarial-infections.R` : code maternal non-malarial infections   
* `6-covariates.R` : generate relevant covariates that were not in original trial data sets  
* `7-medication.R` : code maternal medications  
* `8-analysis-data.R` : with all potential mediators, outcomes, and covariates to create full analysis data sets for different age groups (monthly, bimonthly, quarterly, biannually, and annually)


**`2-analysis` :** folder containing analysis scripts. To rerun all scripts in this subdirectory, run the bash script `0-analysis-bash.sh`.    
* `1-aim1-stratified.R` : obtain total effects of IPTp-SP vs. DP on quarterly growth Z-scores and incidences of stunting and wasting (aim 1)   
* `2a-intervention-mediator-analysis.R` : obtain intervention-mediator effects   
* `2b-mediator-outcome-analysis.R` : obtain mediator-outcome effects   
* `3-single-mediator-analysis.R` : run mediation analysis to obtain mediated and direct effects    
* `4-sensitivity-analysis.R` : test unmeasured mediator-outcome confounding   


**`3-figure-scripts` :** folder containing figure scripts. First, run `2-analysis/0-analysis-bash.sh` to recreate the result data sets in `7-results`. Then, to rerun all scripts in this subdirectory, run the bash script `0-fig-bash.sh`.     
* `1a-fig-zscore-distributions.R` :  quarterly distributions of LAZ and WLZ   
* `1b-fig-incidence-quarterly.R` : quarterly distributions of stunting and wasting incidences    
* `2-fig-z-score-by-age.R` : monthly distributions of LAZ and WLZ by age and gestational age-corrected age (CA)    
* `3-fig-intervention-mediator.R` : forest plot of intervention-mediator associations  
* `4-fig-volcano-tx.R` : volcano plot for maternal inflammation-related proteins at delivery by IPTp-DP vs. SP among all gravidae     
* `5a-fig-mediator-outcome.R` : associations between potential non-inflammation-related mediators and mean LAZ and WLZ     
* `5b-fig-mediator-outcome-Olink.R` : heatmap for associations between maternal inflammation-related proteins at delivery and mean LAZ and WLZ    
* `6a-fig-mediation.R` : total effects and mediated effects on both Z-score and incidence outcomes    
* `6b-fig-mediation-Olink-heatmap.R` : heatmap for mediated effects through main Olink mediators on LAZ and WLZ    
* `6c-fig-mediation-Olink.R` : total effects and mediated effects on incidences of child stunting and wasting    


**`4-table-scripts` :** folder containing table scripts. Before running scripts in this subdirectory, first run `2-analysis/0-analysis-bash.sh` to recreate the result data sets in `7-results`.    
* `1-table-total-effect-binary.R` : unadjusted total effects of IPTp on incidence outcomes   
* `2-table-Olink.R` : table for Olink inflammation biomarker inclusion for intervention-mediator, mediator-outcome, and mediation analyses   
* `3-table-mediator-outcome.R` : mediator-outcome associations between non-Olink mediators and incidence outcomes   
* `4-table-mediated-effects.R` : IPTp's gravidity-stratified mediated (indirect) effects on Z-score outcomes through main mediators   


**`5-figures` :** folder containing figure files.

**`6-tables` :** folder containing table files.

**`7-results` :** folder containing analysis results objects.


Contributors: Jade Benjamin-Chung, Yanwei (Iris) Tong, Suhi Hanif, Anna T. Nguyen
