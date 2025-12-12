# Script information ------------------------------------------------------

# Title: Stock assessment and short term forecast for ane.27.8 (WGHANSA 2024)
#        4) Extract results of interest, write TAF output tables

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/20

# Load libraries ----------------------------------------------------------

library(icesTAF)

# Working directory and folders -------------------------------------------

# check working directory

getwd()

# create model folder and subfolders using the function mkdir from icesTAF

mkdir("output")
mkdir("output/run")
mkdir("output/retro")
mkdir("output/stf")
mkdir("output/toFLR")

# # Outputs from stock assessment -------------------------------------------
# 
# source("output_01_run.R")
# 
# # Outputs from retro ------------------------------------------------------
# 
# source("output_02_retro.R")
# 
# # Outputs from stf --------------------------------------------------------
# 
# source("output_03_stf.R")
# 
# # Outputs: create FLR objects ---------------------------------------------
# 
# source("output_04_toFLR.R")

# Session info ------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------



