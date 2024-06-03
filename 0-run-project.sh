#!/bin/bash

cd 1-dm
bash 0-dm-bash.sh 

cd ..
cd 2-analysis
bash 0-analysis-bash.sh  

cd ..
cd 3-figure-scripts
bash 0-fig-bash.sh  
