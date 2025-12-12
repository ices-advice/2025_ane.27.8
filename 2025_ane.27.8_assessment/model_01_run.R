# Script information ------------------------------------------------------

# run stock assessment for ane.27.8 

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/23

# Load packages -----------------------------------------------------------

library(r4ss)

# Copy updated assessment folder from data ---------------------------------------------------------------

file.copy(from="data/ane.27.8_ss3_2025assessment_model",
          to="model", 
          overwrite=T,recursive=T) # SS model

# !!!!!!!!! IN 2025 changed manually AgeSel_ max and init values from 20 to 30
# DUE TO CONVERGENCE ISSUES

#second option, copy manually the .par file from 2024 assessment
# add values for yearly params, 0 value
#change starter file to use par file as starter

#the second option was used and presented

# Copy ss executable from boot/software ---------------------------------------------------------------

file.copy(from="boot/software/ss.exe",
          to="model/ane.27.8_ss3_2025assessment_model", 
          overwrite=T,recursive=T) # SS model

# Run SS model---------------------------------------------------------------

r4ss::run(dir="model/ane.27.8_ss3_2025assessment_model",exe="ss",skipfinished=F)


