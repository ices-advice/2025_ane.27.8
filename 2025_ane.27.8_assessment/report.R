# Script information ------------------------------------------------------

# Title: Stock assessment and short term forecast for ane.27.8 (WGHANSA 2024)
#        4) Prepare tables and figures for the report

# Before running the script in folder output we have: 
#         XXXX
# After running the script in folder report we have: 
#         XXXX

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/22

# Load libraries ----------------------------------------------------------

library(icesTAF)
library(flextable)

# Working directory and folders -------------------------------------------

#check working directory

getwd()

# create model folder and subfolders using the function mkdir from icesTAF

mkdir("report")

# Create report figs and tabs ---------------------------------------------

# source("report_01_run.R")
# source("report_02_stf.R")
# source("report_03_retro.R")

# Session info ------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------

