# Script information ------------------------------------------------------

# run stock assessment retros for ane.27.8 

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/20

# Load packages -----------------------------------------------------------

library(r4ss)
library(ss3diags)

# Run retros -----------------------------------------------------------

dir_retro<-"model/ane.27.8_ss3_2025assessment_model"

r4ss::retro(dir = dir_retro, exe = "ss", years = 0:-5, verbose = FALSE)


