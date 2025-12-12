# Script information ------------------------------------------------------

# Title: Stock assessment and short term forecast for ane.27.8 (WGHANSA 2024)
#        2) Prepare assessment data (after the bootstrap procedure) and write TAF data tables

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/20

# Load libraries ----------------------------------------------------------

library(icesTAF)
library(r4ss)
library(dplyr)
library(flextable)
library(officer)

# Working directory and folders -------------------------------------------

# check working directory

getwd()

# create data folder using the function mkdir in icesTAF
# the advantage of this function is that you can create more than one folder at the same time
# e.g. mkdir(c("data1","data2"))

mkdir("data")

# Read data ---------------------------------------------------------------

# list files in the bootstrap/file folder

lf<-list.files(taf.boot.path("data"))

# # read reference points from bootstrap/data

refpts <- read.csv(file.path(taf.boot.path("data"), lf[grep("refpts", lf)]))
refpts

# indices data

bioman<-readxl::read_xlsx("boot/data/BIOMAN_1987-2025.xlsx",na = "NA")

pelgas<-readxl::read_xlsx("boot/data/PELGAS_1987-2025.xlsx",na = "NA")

juvena<-readxl::read_xlsx("boot/data/JUVENA_2003-2025.xlsx",na = "NA")

#catch data

ctot<-readxl::read_xlsx("boot/data/ctot_2025.xlsx",na = "NA")

cage_wage<-readxl::read_xlsx("boot/data/cage_wage_2025.xlsx",na = "NA")


# Read and fill SS files ---------------------------------------------------------------

#set assesment year
ass.yr<-2025

# Read all SS inputs files from latest assessment folder

model_dir<-"boot/initial/software/previous_year_model_files"
inputs<-SS_read(dir=model_dir)
inputs_new<-inputs

## Update data file

inputs_new$dat$endyr<-ass.yr

#bioman and pelgas fleetnames
inputs_new$dat$fleetinfo$fleetname[3:4]<-c("Pelgas_survey","Bioman_survey")

#catch new year
ss_ctot1<-subset(inputs$dat$catch)
ss_ctot1_newy<-subset(ss_ctot1,year==ass.yr-1)
ss_ctot1_newy$year<-ass.yr
ss_ctot1_newy$catch<-round(subset(ctot,Year==ass.yr)$ctot,0)
ss_ctot1_new<-rbind(ss_ctot1,ss_ctot1_newy)

#catch previous year
ss_ctot1_new[ss_ctot1_new$year==ass.yr-1,"catch"]<-round(subset(ctot,Year==ass.yr-1)$ctot,0)

#ss_ctot1_new$catch[-c(1:2)]-round(ctot$ctot,0) check differences

inputs_new$dat$catch<-as.data.frame(ss_ctot1_new)

#pelgas
ss_index3<-subset(inputs$dat$CPUE,index==3)
ss_index3_newy<-subset(ss_index3,year==ass.yr-1)
ss_index3_newy$year<-ass.yr
ss_index3_newy$obs<-subset(pelgas,Year==ass.yr)$Biomass
ss_index3_newy$se_log<- sqrt(log(subset(pelgas,Year==ass.yr)$CV^2+1))
ss_index3_new<-rbind(ss_index3,ss_index3_newy)

#round((ss_index3_new$obs-na.omit(pelgas$Biomass))/na.omit(pelgas$Biomass),6) check differences

#bioman
ss_index4<-subset(inputs$dat$CPUE,index==4)
ss_index4_newy<-subset(ss_index4,year==ass.yr-1)
ss_index4_newy$year<-ass.yr
ss_index4_newy$obs<-subset(bioman,Year==ass.yr)$Biomass
ss_index4_newy$se_log<- sqrt(log(subset(bioman,Year==ass.yr)$CV^2+1))
ss_index4_new<-rbind(ss_index4,ss_index4_newy)

#round((ss_index4_new$obs-na.omit(bioman$Biomass))/na.omit(bioman$Biomass),6) check differences


#juvena
ss_index5<-subset(inputs$dat$CPUE,index==5)
ss_index5_newy<-subset(ss_index5,year==ass.yr-1)
ss_index5_newy$year<-ass.yr
ss_index5_newy$obs<-subset(juvena,Year==ass.yr)$Num_juveniles/1000
ss_index5_newy$se_log<- sqrt(log(0.25^2+1))
ss_index5_new<-rbind(ss_index5,ss_index5_newy)

#round((1000*ss_index5_new$obs-na.omit(juvena$Num_juveniles))/na.omit(juvena$Num_juveniles),6) check differences

inputs_new$dat$CPUE<-as.data.frame(rbind(ss_index3_new,ss_index4_new,ss_index5_new))

#age component

#catch at age for only sem 1 ass.yr
ss_age_fleet1<-subset(inputs$dat$agecomp,fleet==1)
ss_age_fleet1_newy<-subset(ss_age_fleet1,year==ass.yr-1)
ss_age_fleet1_newy$year<-ass.yr
ss_age_fleet1_newy[,c("a1","a2","a3")]<-subset(cage_wage,Year==ass.yr&Semester==1)[,c("Nage_Age1","Nage_Age2","Nage_Age3plus")]
ss_age_fleet1_new<-rbind(ss_age_fleet1,ss_age_fleet1_newy)

#catch at age for only sem 2 ass.yr-1
ss_age_fleet2<-subset(inputs$dat$agecomp,fleet==2)
ss_age_fleet2_newy<-subset(ss_age_fleet2,year==ass.yr-2)
ss_age_fleet2_newy$year<-ass.yr-1
ss_age_fleet2_newy[,c("a1","a2","a3")]<-subset(cage_wage,Year==ass.yr-1&Semester==2)[,c("Nage_Age1","Nage_Age2","Nage_Age3plus")]
ss_age_fleet2_new<-rbind(ss_age_fleet2,ss_age_fleet2_newy)

ss_age_fleet12_new<-rbind(ss_age_fleet1,ss_age_fleet1_newy,ss_age_fleet2,ss_age_fleet2_newy)
ss_age_fleet12_new<-ss_age_fleet12_new[order(ss_age_fleet12_new$year),]

#update ass.yr-1, catch has been corrected
ss_age_fleet12_new[ss_age_fleet12_new$year==ass.yr-1,c("a1","a2","a3")]<-subset(cage_wage,Year==ass.yr-1)[,c("Nage_Age1","Nage_Age2","Nage_Age3plus")]

# diff_cage<-which(round(ss_age_fleet12_new[,c("year","a1","a2","a3")] - subset(cage_wage,Nage_Age1>0)[,c("Year","Nage_Age1","Nage_Age2","Nage_Age3plus")],1)[,2]!=0)
# ss_age_fleet12_new[diff_cage,]

#pelgas age
ss_age_fleet3<-subset(inputs$dat$agecomp,fleet==3)
ss_age_fleet3_newy<-subset(ss_age_fleet3,year==ass.yr-1)
ss_age_fleet3_newy$year<-ass.yr
ss_age_fleet3_newy[,c("a1","a2")]<-subset(pelgas,Year==ass.yr)[,c("N1","N2")]
ss_age_fleet3_newy[,c("a3")]<-sum(subset(pelgas,Year==ass.yr)[,c("N3","N4plus")],na.rm=T)
ss_age_fleet3_new<-rbind(ss_age_fleet3,ss_age_fleet3_newy)

#bioman age
ss_age_fleet4<-subset(inputs$dat$agecomp,fleet==4)
ss_age_fleet4_newy<-subset(ss_age_fleet4,year==ass.yr-1)
ss_age_fleet4_newy$year<-ass.yr
ss_age_fleet4_newy[,c("a1","a2","a3")]<-subset(bioman,Year==ass.yr)[,c("N1","N2","N3plus")]
ss_age_fleet4_new<-rbind(ss_age_fleet4,ss_age_fleet4_newy)


inputs_new$dat$agecomp<-as.data.frame(rbind(ss_age_fleet12_new,ss_age_fleet3_new,ss_age_fleet4_new))

## Update wtatage file

# fleet 0 contains begin season pop WT
# fleet -1 contains mid season pop WT
# fleet -2 contains w*maturity*fecundity

# stock wt from BIOMAN fleet4=fleet3=fleet5=fleet0=fleet0=fleet-1=fleet=-2 (same both sem)
# catch wt by sem fleet1=fleet2 (different each sem)

wt_new<-subset(inputs$wtatage,year==ass.yr-1)
wt_new$year<-ass.yr
wt0<-subset(wt_new,fleet==0)
wt0[,c("1","2","3")]<-subset(bioman,Year==ass.yr)[,c("W1","W2","W3plus")]/1000

wt1<-subset(wt_new,fleet==1)
wt1[,c("1","2","3")]<-subset(cage_wage,Year==ass.yr)[,c("Wage_Age1","Wage_Age2","Wage_3plus")]/1000

wt2<-subset(wt_new,fleet==2)
wt2[,c("1","2","3")]<-wt1[,c("1","2","3")]

wt3<-subset(wt_new,fleet==3)
wt3[,c("1","2","3")]<-wt0[,c("1","2","3")]

wt4<-subset(wt_new,fleet==4)
wt4[,c("1","2","3")]<-wt0[,c("1","2","3")]

wt5<-subset(wt_new,fleet==5)
wt5[,c("1","2","3")]<-wt0[,c("1","2","3")]

wt6<-subset(wt_new,fleet==-1)
wt6[,c("1","2","3")]<-wt0[,c("1","2","3")]

wt7<-subset(wt_new,fleet==-2)
wt7[,c("1","2","3")]<-wt0[,c("1","2","3")]

wt_new<-rbind(wt0,wt1[1,],wt2[1,],wt3,wt4,wt5,wt6,wt7)

inputs_new$wtatage<-rbind(inputs$wtatage[inputs$wtatage$year>0,],wt_new)

# catch ass.yr-1
new_lines_seas2<-subset(inputs_new$wtatage,year==ass.yr-1&fleet %in% 1:2)
new_lines_seas2$seas<-2

inputs_new$wtatage<-rbind(inputs_new$wtatage,new_lines_seas2)

inputs_new$wtatage[inputs_new$wtatage$year==ass.yr-1&inputs_new$wtatage$fleet==1,c("1","2","3")]<-
subset(cage_wage,Year==ass.yr-1)[,c("Wage_Age1","Wage_Age2","Wage_3plus")]/1000

inputs_new$wtatage[inputs_new$wtatage$year==ass.yr-1&inputs_new$wtatage$fleet==2,c("1","2","3")]<-
  subset(cage_wage,Year==ass.yr-1)[,c("Wage_Age1","Wage_Age2","Wage_3plus")]/1000


#compute average weigths at age for the previous 3 year to the assessment years
wt_new_fore<-inputs_new$wtatage  %>% filter(year %in% c(ass.yr-(1:3)))%>% group_by (seas,sex,bio_pattern,birthseas,fleet)%>%
  summarize(year=-(ass.yr+1),`0`=mean(`0`,na.rm=T),`1`=mean(`1`,na.rm=T),`2`=mean(`2`,na.rm=T),`3`=mean(`3`,na.rm=T))

wt_new_fore[wt_new_fore$fleet %in%1:2 & wt_new_fore$seas==2,"year"]<- -ass.yr

inputs_new$wtatage<-rbind(inputs_new$wtatage,wt_new_fore)


#missing value age 3 sem 2, no cage, input with mean values for forecast

inputs_new$wtatage[which(is.na(inputs_new$wtatage[,"3"])),"3"]<-filter(wt_new_fore,seas==2&fleet %in% 1:2)[,"3"]


## Update control file

#blocks last year
inputs_new$ctl$Block_Design[[1]][4]<-ass.yr+2
inputs_new$ctl$Block_Design[[2]][4]<-ass.yr+2

#last main recdev
inputs_new$ctl$MainRdevYrLast<-ass.yr

#selectivity random walk last year
inputs_new$ctl$age_selex_parms$dev_maxyr[inputs_new$ctl$age_selex_parms$dev_maxyr==ass.yr-1]<-ass.yr


## Update starter file (nothing to update)
#inputs_new$start


## Update forecast file

inputs_new$fore$Ydecl<-ass.yr
inputs_new$fore$Yinit<-ass.yr

inputs_new$FirstYear_for_caps_and_allocations<-ass.yr+1


# write all new updated SS files to the new assessment folder

SS_write(inputs_new,dir="data/ane.27.8_ss3_2025assessment_model",overwrite=TRUE)

# read and write forecast file alone, it needs the readAll=T option
fore_file<-SS_readforecast("boot/initial/software/previous_year_model_files/forecast.ss",readAll=T)

fore_file$Ydecl<-ass.yr
fore_file$Yinit<-ass.yr

fore_file$FirstYear_for_caps_and_allocations<-ass.yr+1
fore_file$Fcast_years[1:4]<-c(-3,-1,-3,-1)
fore_file$Bmark_years[1:4]<-c(-3,-1,-3,-1)

SS_writeforecast(fore_file,dir="data/ane.27.8_ss3_2025assessment_model",writeAll=T,overwrite=TRUE)


# save updates on data

oldcage<-subset(inputs$dat$agecomp,fleet %in% 1:2 & year >ass.yr-3)[,c("year","fleet","a1","a2","a3")]
newcage<-subset(inputs_new$dat$agecomp,fleet %in% 1:2 & year >ass.yr-3)[,c("year","fleet","a1","a2","a3")]
oldcage$ass.yr<-ass.yr-1
newcage$ass.yr<-ass.yr
oldc<-subset(inputs$dat$catch,year > ass.yr-3)[,1:4]
newc<-subset(inputs_new$dat$catch,year > ass.yr-3)[,1:4]
oldc$ass.yr<-ass.yr-1
newc$ass.yr<-ass.yr

ft1 <- colformat_double(
  x = flextable(rbind(oldcage,newcage)),
  big.mark = "", digits = 0, na_str = "N/A"
)

ft2 <- colformat_double(
  x = flextable(rbind(oldc,newc)),
  big.mark = "", digits = 0, na_str = "N/A"
)

read_docx() |> 
  body_add_flextable(ft1) |> 
  body_add_break() |> 
  body_add_flextable(ft2) |> 
  body_add_break() |> 
  print(target="data/old_new_catch.docx")  

# Prepare TAF tables ------------------------------------------------------

# natural maturity table
natmort <- as.data.frame(inputs_new$ctl$MG_parms[c(1,2,4,6),"INIT"])

# historical numbers at age in the catch table
catage <- cage_wage[,c("Year","Semester","Nage_Age1","Nage_Age2","Nage_Age3plus")]

# catch in tonnes
catch <- ctot

# numbers at age in the acoustic survey 
natage_acoustic <- subset(inputs_new$dat$agecomp,fleet==3 & year%in% inputs_new$dat$styr:inputs_new$dat$endyr,c("year","a1","a2","a3"))

# total biomass in the acoustic survey 
survey_acoustic <- subset(inputs_new$dat$CPUE,index==3,c("year","obs"))

# numbers at age in the DEPM survey 
natage_DEPM <- subset(inputs_new$dat$agecomp,fleet==4 & year%in% inputs_new$dat$styr:inputs_new$dat$endyr,c("year","a1","a2","a3"))

# ssb estimated in the DEPM survey (biomass)
survey_DEPM <- subset(inputs_new$dat$CPUE,index==4,c("year","obs"))

# total numbers in the recruit survey 
survey_recruits <- subset(inputs_new$dat$CPUE,index==5,c("year","obs"))

# weight at age in the catch
waca <- cage_wage[,c("Year","Semester","Wage_Age1","Wage_Age2","Wage_3plus")]

# weight at age in the stock
west <- subset(inputs_new$wtatage, fleet=="3" & year %in% inputs_new$dat$styr:inputs_new$dat$endyr,c("year","0","1","2","3"))

# fecundity
fecundity <- NA

# maturity
maturity <- as.data.frame(c(0,1,1,1))



# Write TAF tables in data folder -----------------------------------------

write.taf(list(natmort=natmort, catage = catage, catch = catch, 
               natage_acoustic = natage_acoustic, survey_acoustic = survey_acoustic,
               survey_DEPM = survey_DEPM, survey_recruits = survey_recruits,
               waca = waca, west = west, maturity = maturity),dir="./data")

# Save data in RData file  -----------------------------------------

save.image("./data/inputData.RData")

# Script info -------------------------------------------------------------

sessionInfo()

# End of script -----------------------------------------------------------