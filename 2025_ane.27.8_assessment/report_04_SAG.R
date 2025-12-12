################################################################################
#  WGHANSA ane.27.8 - Standard Graphs table                                    #
#------------------------------------------------------------------------------#
#   created by: Sonia Sanchez (AZTI-Tecnalia) for nothern hake                 #
#             09/05/2023                                                       #
#   modified: Leire Citores for anchovy in the Bay of Biscay                   #
################################################################################

# filling Standard Graphs table and saving it in ICES repository


# Copyright: AZTI, 2023
# Author: Sonia Sanchez (AZTI) (<ssanchez@azti.es>)
#
# Distributed under the terms of the GNU GPLv3

# Based on code from Hans Gerritsen modified: 2025/11/20

rm(list = ls())

#==============================================================================
# LIBRARIES                                                                ----
#==============================================================================

library(icesSAG)
library(openxlsx)


#==============================================================================
# LOAD data and assessment output                                          ----
#==============================================================================

load("data/inputData.RData") 

load("output/run/output_run.RData")


#==============================================================================
# TOKEN                                                                    ----
#==============================================================================

# You can generate a token like this:
# first log in on  https://standardgraphs.ices.dk/manage
# then go to sg.ices.dk/manage/CreateToken.aspx
# paste the token below after SG_PAT=

cat("# Standard Graphs personal access token",
    "SG_PAT=e9351534-20ac-4ad4-9752-98923e011213", # replace with your own token
    sep = "\n",
    file = "~/.Renviron_SG")
options(icesSAG.use_token = TRUE)



#==============================================================================
# DATA  needed                                                                   ----
#==============================================================================

stockcode <- 'ane.27.8'

#sumTab <- read.xlsx(file.path(getwd(),"report","table","Tab10.6.xlsx"))

summ_out<-SSsummarize(list(out,out))

SSB<-filter(summ_out$SpawnBio,Yr %in% 1987:ass.yr)[,1]
SSBupper<-filter(summ_out$SpawnBioUpper,Yr %in% 1987:ass.yr)[,1]
SSBlower<-filter(summ_out$SpawnBioLower,Yr %in% 1987:ass.yr)[,1]

Rec<-filter(summ_out$recruits,Yr %in% 1987:ass.yr)[,1]
Recupper<-filter(summ_out$recruitsUpper,Yr %in% 1987:ass.yr)[,1]
Reclower<-filter(summ_out$recruitsLower,Yr %in% 1987:ass.yr)[,1]

Fvalue<-filter(summ_out$Fvalue,Yr %in% 1987:ass.yr)[,c(1,4)]
Fvalueupper<-filter(summ_out$FvalueUpper,Yr %in% 1987:ass.yr)[,c(1,4)]
Fvaluelower<-filter(summ_out$FvalueLower,Yr %in% 1987:ass.yr)[,c(1,4)]

Catches<-round(ctot %>% group_by(Year) %>% summarize(ctot=sum(ctot)) %>% filter(Year %in% 1987:ass.yr),2)
Catches$Yr<-Catches$Year

Catch_F <- Catches %>%
  left_join(Fvaluelower, by='Yr')%>%
  left_join(Fvalue, by='Yr') %>%
  left_join(Fvalueupper, by='Yr') 

names(Catch_F)[4:6]<-c("Fvaluelower","Fvalue","Fvalueupper")
# 
# tab_sag<-cbind(Catch_F,Reclower,Rec,Recupper,SSBlower,SSB,SSBupper)
# 
# tab_for_sag<-write.csv(tab_sag,"report/tab_sag.csv")

# STOCK INFORMATION

refpts <- read.csv("boot/data/refpts_ane.27.8.csv", header=T)
Blim=refpts$Blim
Bpa=refpts$Bpa	


# icesSAG:::validNames("stockInfo")


info <- stockInfo( StockCode = 'ane.27.8',           # grep("hke", icesVocab::getCodeList("ICES_StockCode")$Key, value = TRUE)
                   AssessmentYear = ass.yr, 
                   ContactPerson = 'lcitores@azti.es', 
                   StockCategory = 1.0, 
                   ModelType = 'AL',                         # icesVocab::getCodeList("AssessmentModelType") 
                   ModelName = 'SS3',                        # icesVocab::getCodeList("AssessmentModelName") 
                   Blim = Blim, Bpa = Bpa,
                   Fage='apical', RecruitmentAge = 0, 
                   CatchesLandingsUnits = 't', 
                   ConfidenceIntervalDefinition = '95%', 
                   RecruitmentUnits = 'NE3', 
                   StockSizeUnits = 't', 
                   StockSizeDescription = 'SSB', # icesVocab::getCodeList("StockSizeIndicator") 
                   Purpose = 'Advice'
)

info$StockCategory <- 1 

# icesSAG:::validNames("stockFishdata")
Catch_F[is.na(Catch_F)]<-0
fishdata <- stockFishdata( Year = Catches$Year, 
                           Low_Recruitment = Reclower,Recruitment = Rec,  High_Recruitment = Recupper, 
                           Low_StockSize = SSBlower,StockSize = SSB,  High_StockSize = SSBupper,
                           Catches=Catches$ctot,Landings = Catches$ctot,Low_FishingPressure = Catch_F$Fvaluelower, FishingPressure = Catch_F$Fvalue, High_FishingPressure = Catch_F$Fvalueupper
)



write.csv(fishdata,"report/tab_sag.csv")


#==============================================================================
# UPLOAD                                                                   ----
#==============================================================================

key <- icesSAG::uploadStock(info, fishdata)

xmlfile <-icesSAG::createSAGxml(info, fishdata)

cat(xmlfile,  file= "output/myxmlfile.xml")

out <- readSAGxml(xmlfile)

sag_use_token(TRUE)
findAssessmentKey('ane.27.8', 2024, full = TRUE)
