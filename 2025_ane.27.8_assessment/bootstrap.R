# Script information ------------------------------------------------------

# Title: Stock assessment and short term forecast for ane.27.8 (WGHANSA 2025)
#        1) TAF boot procedure to set up required data and software

# See also: https://github.com/ices-taf/doc/wiki/Creating-a-TAF-analysis 

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/06/15, modified: 2025/11/20

# Load libraries ----------------------------------------------------------

library(icesTAF)
library(help=icesTAF)

# Check the working directory ---------------------------------------------

getwd()

# Create repository and R project -----------------------------------------

# Clone the empty github repository created by ICES secretariat
# The cloned repository is in C:/use/GitHub/ices-taf/2022_ane.27.8_assessment
# Create an R project in that folder

# Make the skeleton -------------------------------------------------------

# create initial directories and R scripts for a new TAF analysis
# 2022_ane.27.8_assessment    
# ¦--boot   
# ¦   °--initial 
# ¦       °--data
# ¦--data.R      
# ¦--model.R     
# ¦--output.R    
# °--report.R    

taf.skeleton()

# Upload initial data -----------------------------------------------------

# Include data in the folder: 2024_ane.27.8_assessment\boot\initial\data 
# It can be done by-hand or copying it using the following command: file.copy(from,to)
# 2024: uploaded directly by Leire Ibaibarriaga

# file.copy(from="C:/use/GitHub/WKBANSP2024/data/cage_bysem.csv",
#           to="./boot/initial/data/cage_bysem.csv", overwrite=T)
# 
# file.copy(from="C:/use/GitHub/WKBANSP2024/data/wage_bysem.csv",
#           to="./boot/initial/data/wage_bysem.csv", overwrite=T)
# 
# file.copy(from="C:/use/GitHub/WKBANSP2024/data/BIOMAN.xlsx",
#           to="./boot/initial/data/BIOMAN.xlsx", overwrite=T)
# 
# file.copy(from="C:/use/GitHub/WKBANSP2024/data/PELGAS.xlsx",
#           to="./boot/initial/data/PELGAS.xlsx", overwrite=T)
# 
# file.copy(from="C:/use/GitHub/WKBANSP2024/data/JUVENA.xlsx",
#           to="./boot/initial/data/JUVENA.xlsx", overwrite=T)

# Document data and create the data.bib file ------------------------------

# use function draft.data() to create the data.bib file
?draft.data

# NOTE that the columns of the csv files should be separated by comma (",")

# # data for stock assessment (unique file)
# 
# draft.data(originator="WGHANSA", 
#            year=2023, 
#            title="CBBM data for ane.27.8", 
#            period="1987-2023",
#            source="file",
#            data.files=c("realdata_1987-2023_nov2023.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=F)
# 
# # data with reference points
# 
# draft.data(originator="WKPELA", 
#            year=2013, 
#            title="Reference points for ane.27.8", 
#            period=" ",
#            source="file",
#            data.files=c("refpts_ane.27.8.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)
# data with reference points
# Note this is the first file, so: append=F

# # data with weights-at-age at spawning (needed to create the FLBiol and FLStock objects)
# 
# draft.data(originator="WGHANSA", 
#            year=2023, 
#            title="Weights at age at spawning", 
#            period="1987-2023",
#            source="file",
#            data.files=c("wage.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)
draft.data(originator="WKBANSP",
           year=2024,
           title="Reference points for ane.27.8",
           period=" ",
           source="file",
           data.files=c("refpts_ane.27.8.csv"),
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=F) 

# data with catch-at age and weights-at-age in the catch (needed to create the FLStock objects)
# data with catch-at age and weights-at-age in the catch

# draft.data(originator="WGHANSA", 
#            year=2025, 
#            title="Catch-at-age by semester", 
#            period="1987-2025",
#            source="file",
#            data.files=c("cage_bysem_2025.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)
# 
# draft.data(originator="WGHANSA", 
#            year=2025, 
#            title="Weight-at-age in the catch by semester", 
#            period="1987-2025",
#            source="file",
#            data.files=c("wage_bysem_2025.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)

draft.data(originator="WGHANSA", 
           year=2025, 
           title="Catch and Weight-at-age in the catch by semester", 
           period="1987-2025",
           source="file",
           data.files=c("cage_wage_2025.xlsx"), 
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)

draft.data(originator="WGHANSA", 
           year=2025, 
           title="Total catch by semester", 
           period="1987-2025",
           source="file",
           data.files=c("ctot_2025.xlsx"), 
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)

# data from BIOMAN, PELGAS and JUVENA surveys

draft.data(originator="WGACEGG", 
           year=2025, 
           title="Bioman survey data", 
           period="1987-2025",
           source="file",
           data.files=c("BIOMAN_1987-2025.xlsx"), 
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)

draft.data(originator="WGACEGG", 
           year=2025, 
           title="PELGAS survey data", 
           period="1987-2025",
           source="file",
           data.files=c("PELGAS_1987-2025.xlsx"), 
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)

draft.data(originator="WGACEGG", 
           year=2025, 
           title="JUVENA survey data", 
           period="2003-2025",
           source="file",
           data.files=c("JUVENA_2003-2025.xlsx"), 
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)
# 
# # data with stock assessment results from previous year
# 
# draft.data(originator="WGHANSA", 
#            year=2022, 
#            title="ICES summary table for ane.27.8 from WGHANSA 2022", 
#            period="1987-2022",
#            source="file",
#            data.files=c("table_summary_ices_wghansa2023nov.csv"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)

# example of data automatically downloaded from internet
# Data automatically downloaded from internet

# Official and historical catches from https://www.ices.dk/data/dataset-collections/Pages/Fish-catch-and-stock-assessment.aspx
# note that in this case it is a zip file that will need to be unzipped at some point

# draft.data(originator="FAO", 
#            year=2019, 
#            title="Historical landings", 
#            period="1950-2010",
#            source="https://www.ices.dk/data/Documents/CatchStats/HistoricalLandings1950-2010.zip", 
#            data.files=c("HistoricalLandings1950_2010.zip"), 
#            data.scripts=NULL, 
#            file="boot/data.bib", 
#            append=T)
draft.data(originator="FAO",
           year=2025,
           title="Official nominal catches",
           period="1950-2010",
           source="https://www.ices.dk/data/Documents/CatchStats/OfficialNominalCatches.zip",
           data.files=c("OfficialNominalCatches.zip"),
           data.scripts=NULL,
           file="boot/data.bib",
           append=T)
draft.data(originator="FAO",
           year=2019,
           title="Historical landings",
           period="1950-2010",
           source="https://www.ices.dk/data/Documents/CatchStats/HistoricalLandings1950-2010.zip",
           data.files=c("HistoricalLandings1950_2010.zip"),
           data.scripts=NULL,
           file="boot/data.bib",
           append=T)

# # ICES areas downloaded from internet using an script
# # the script should be located in the folder boot
# 
# draft.data(originator = "ICES",
#            title = "ICES Areas ESRI Shapefile",
#            period = FALSE,
#            source = "script",
#            data.files=NULL,
#            data.script = "icesareas",
#            file="boot/data.bib",
#            append=T)

# ICES report template downloaded from internet

draft.data(originator = "ICES",
           title = "ICES TAF Word template for report automation",
           period = FALSE,
           source = "https://github.com/ices-taf/doc/raw/master/reportTemplate.docx",
           data.files = "reportTemplate.docx",
           data.scripts = NULL,
           file = "boot/data.bib", 
           append=T)

# # ICES section 3 for anchovy ()
# 
# draft.data(originator="WGHANSA", 
# draft.data(originator = "ICES",
#            title = "ICES TAF Word template for report automation",
#            source="file",
#            data.files=c("sec_03_template.docx"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            period = FALSE,
#            source = "https://github.com/ices-taf/doc/raw/master/reportTemplate.docx",
#            data.files = "reportTemplate.docx",
#            data.scripts = NULL,
#            file = "boot/data.bib", 
#            append=T)

# ICES section 3 for anchovy (specific for section 3 for automatic numbering etc)

draft.data(originator="WGHANSA",
           title="ICES TAF Word template for report automation",
           source="file",
           data.files=c("sec_02_template.docx"),
           data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
           file="boot/data.bib",
           append=T)

# ICES csl file to format references 

# draft.data(originator="WGHANSA", 
#            title="ICES cls file", 
#            # source="https://github.com/citation-style-language/styles-distribution/blob/master/ices-journal-of-marine-science.csl",
#            # title="ICES csl file", 
#            # source="https://raw.githubusercontent.com/citation-style-language/styles-distribution/refs/heads/master/ices-journal-of-marine-science.csl",
#            data.files=c("ices-journal-of-marine-science.csl"), 
#            data.scripts=NULL, # we need to set this to NULL or it will automatically take all the files with extension .R
#            file="boot/data.bib",
#            append=T)

# Anchovy (Engraulis encrasicolus) downloaded from internet ---------------

draft.data(originator = "Gervais et Boulart",
           title = "Image anchovy (Engraulis encrasicolus)",
           period = FALSE,
           source = "https://upload.wikimedia.org/wikipedia/commons/2/2e/Engraulis_encrasicolus_Gervais.jpg",
           data.files = "anchovy.jpg",
           data.scripts = NULL,
           file = "boot/data.bib", 
           append=T)

# Upload software ---------------------------------------------------------

# create the folder "XX_ane.27.8_assessment\boot\initial\software"

if (!dir.exists("boot/initial/software")){
  dir.create("boot/initial/software")
}

if (!dir.exists("boot/initial/software/previous_year_model_files")){
  dir.create("boot/initial/software/previous_year_model_files")
}

# copy the SS3 model folder from 2024 benchmark to the folder boot/initial/software

file.copy(from="C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/ane.27.8_ss3_2024assessment_model/ane.ctl",
          to="./boot/initial/software/previous_year_model_files", 
          overwrite=T,recursive=T)
file.copy(from="C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/ane.27.8_ss3_2024assessment_model/ane.dat",
          to="./boot/initial/software/previous_year_model_files", 
          overwrite=T,recursive=T)
file.copy(from="C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/ane.27.8_ss3_2024assessment_model/forecast.ss",
          to="./boot/initial/software/previous_year_model_files", 
          overwrite=T,recursive=T)
file.copy(from="C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/ane.27.8_ss3_2024assessment_model/starter.ss",
          to="./boot/initial/software/previous_year_model_files", 
          overwrite=T,recursive=T)
file.copy(from="C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/ane.27.8_ss3_2024assessment_model/wtatage.ss",
          to="./boot/initial/software/previous_year_model_files", 
          overwrite=T,recursive=T)
#executable from benchmark
file.copy(from="C:/USE/OneDrive - AZTI/GitHub/WKBANSP2024/ane.27.8_ss3_2024benchmark_final_model/ss.exe",
          to="./boot/initial/software", 
          overwrite=T,recursive=T) # SS model


# Document software and create the software.bib file ----------------------

# use function draft.software to create the software.bib file

?draft.software

# document the JAGS model
# document the Stock Synthesis model

draft.software(package=c("boot/initial/software/ss.exe"), 
               author = "Rick Methot et al.", 
               year = 2024, 
               title = "SS executable",
               version = "V3.30.22.beta",
               source = NULL, 
               file = "boot/software.bib", 
               append = FALSE)


# document the r4ss and the ss3diags packages

library(r4ss)
draft.software("r4ss") # print in console
draft.software("r4ss",
               file = "boot/software.bib",
               append=T)


library(ss3diags)
draft.software("ss3diags") # print in console
draft.software("ss3diags",
               file = "boot/software.bib",
               append=T)

# Process the data.bib and software.bib metafiles -------------------------

# the function taf.boot() processes:
#   the data.bib file and creates/copies all the data files into the folder "boot/data"
#   and 
#   the software.bib file and creates/copies the software into the folder boot/software
?taf.boot

# apply taf.boot

taf.boot(taf=TRUE)

# Session info ------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------



