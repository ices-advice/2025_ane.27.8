# Script information ------------------------------------------------------

# Title: Stock assessment and short term forecast for ane.27.8 (WGHANSA 2024)
#        3) Run stock assessment, short term forecast and retros

# Before running the script in folder data we have: 
#         input.RData 
# After running the script in folder model we have: 
#         folders with the assessment run and the retros

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2022/06/15, modified: 2024/11/21, modified: 2025/11/20

# Load libraries ----------------------------------------------------------

library(icesTAF)

# Working directory and folders -------------------------------------------

# check working directory

getwd()

# create model folder and subfolders using the function mkdir from icesTAF

mkdir("model")

# Run script for stock assessment -----------------------------------------

#source("model_01_run.R")

# Run script for retro ----------------------------------------------------

#source("model_02_retro.R")

# Run scripts for STF -----------------------------------------------------

#source("model_03_stf.R")

# Session info ------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------


