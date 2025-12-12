# Script information ------------------------------------------------------

# run stochastic short-term forecast in SS for ane.27.8 

# based on code for the two hake stocks and southern horsr mackerel

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/25 modified: 2025/11/20

# Load libraries --------------------------------------------------------------

library(r4ss)
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
library(ggplot2)

library(icesAdvice) # for ICES rounding

# Set directories ---------------------------------------------------------

getwd()

# directory of the SS run

run.dir <- "model/ane.27.8_ss3_2025assessment_model"

# directories for STF

stf.dir   <- file.path("model","stf_SS3","stf")
stf_files.dir <- file.path(stf.dir, "stf_files")

dir.create(stf.dir, showWarnings=T, recursive=T)
dir.create(stf_files.dir, showWarnings=T, recursive=T)

# Load basic objects ------------------------------------------------------

# # load reference points
 refpts <- read.csv("boot/data/refpts_ane.27.8.csv", header=T)
 refpts
 Blim <- refpts$Blim
#Blim<-26600

# copy forecast file
file.copy(file.path(run.dir, "forecast.ss"), file.path(stf_files.dir, "forecast.ss"))

# read in the forecast file
fore <- r4ss::SS_readforecast(file = file.path(stf_files.dir, "forecast.ss"),verbose = FALSE)
fore

# Number of forecast years
Nfor <- fore$Nforecastyrs 
Nfor # 2

# read in assessment ouput
replist <- SS_output(dir = run.dir, verbose=TRUE, printstats=TRUE) ## read

# Intermediate year -------------------------------------------------------

# number of years (-1) over which to average F's for the interim year
Naver <- 0  # note that in this case there is no interim year

# exploitation by fleet. note this includes the forecast years
dat <- replist$exploitation
# keep year, seas and fleets
dat <- dat %>% select(-Seas_dur, -F_std, -annual_F, -annual_M)
head(dat) 
tail(dat)

# define start and end years
startyear <- max(dat$Yr)-Nfor-Naver # 2025: last assessment year
endyear   <- max(dat$Yr)-Nfor # 2025
# interyear <- endyear # In our case there is no intermediate year. it is the last assessment year
  
# average of the last (Naver+1) years across seasons and fleets

data_int <- dat %>% 
  filter(Yr>=startyear & Yr<=endyear) %>% 
  select(-Yr) %>% 
  group_by(Seas) %>% 
  summarise_all(mean)
data_int

# plot aver
ggplot(dat, aes(Yr, Commercial_vessels1))+
  geom_point()+
  geom_line()+
  geom_hline(data=data_int, aes(yintercept=Commercial_vessels1), col=2)+
  geom_vline(xintercept=c(startyear, endyear), lty=2)+
  facet_grid(~Seas)+
  ggtitle("Exploitation (F)")
ggplot(dat, aes(Yr, Commercial_vessels2))+
  geom_point()+
  geom_line()+
  geom_hline(data=data_int, aes(yintercept=Commercial_vessels2), col=2)+
  geom_vline(xintercept=c(startyear, endyear), lty=2)+
  facet_grid(~Seas)+
  ggtitle("Exploitation (F)")

## input for intermediate year data. In this case, it is the last year

dimen <- dim(data_int)
nfleet <- ncol(data_int)-1
Year  <- rep(endyear+1, dimen[1]*(dimen[2]-1))
fore_dat_int       <- data.frame(Year)
fore_dat_int$Seas  <- data_int$Seas
fore_dat_int$Fleet <- rep(1:nfleet, each = length(data_int$Seas))
fore_dat_int$F     <- as.vector(as.matrix(subset(data_int, select=-Seas)))
fore_dat_int

# F status quo ------------------------------------------------------------

# # Fsq is the last year F

datmul <- replist$exploitation %>%
  filter(Yr == endyear) %>%
  filter(!is.na(F_std)) %>% # whole year F stored in 1st season
  select(Yr, F_std)

Fsq <- datmul$F_std
Fsq

# F multipliers -----------------------------------------------------------

# Fishing mortality multipliers for F reference points
# Not available for ane.27.8 as it is a short-lived species

# fmult_msyl <- FmsyLower/Fsq
# fmult_msyu <- FmsyUpper/Fsq
# fmult_lim  <- Flim/Fsq
# fmult_msy  <- Fmsy/Fsq
# s1 <- s2 <- s3 <- 20
# Fmult <- c( seq(0, fmult_msyl, length.out=s1), 
#             seq(fmult_msyl+0.01, fmult_msyu, length.out=s2), 
#             seq(fmult_msyu+0.01, fmult_lim+0.05, length.out=s3), 
#             c(3, 4, 5, 6, 7, 1))


Fmult <- c(seq(0, 2, by=0.1),seq(5,15,by=0.5),c(16:20))
if(exists("interpolated_Fmults")){ #for when this script needs to run the second time for interpolated values
  Fmult<-c(Fmult,interpolated_Fmults)
}
#Fmult <- sort(unique(Fmult))
Fmult_names <- paste0("Fmult_", Fmult)

save(Fmult, Fmult_names,nfleet, file = file.path(stf.dir, "Fmult.RData"))

# Recruitment assumption --------------------------------------------------

# Get assessment outputs for geomean recruitment
ass.sum <- SSsummarize(SSgetoutput(dirvec=run.dir))
hist.rec <- as.data.frame(ass.sum$recruits) %>% 
  filter(Yr %in% (ass.sum$startyrs):(ass.sum$endyrs)) %>% .[,1]  # all series
virg.rec <- as.data.frame(ass.sum$recruits) %>% 
  filter(Label == "Recr_Virgin") %>% .[,1]

# GM recruitment
gmrec <- exp(mean(log(hist.rec)))
gmrec

# Create forecast file for default SR model -------------------------------

#Fix the proportion of F in sem 1 to apply in the forecast years
# as the average of the proportion of F in sem 1 in the previous 
# three years to the assessment year (i.e. in 2024, average 2021-2023)
 
Fs_toavg<-filter(dat,Yr %in% 2022:2024) %>% group_by(Yr) %>% 
  summarize(Commercial_vessels1=sum(Commercial_vessels1),Commercial_vessels2=sum(Commercial_vessels2))%>% 
  mutate(propF1=Commercial_vessels1/(Commercial_vessels1+Commercial_vessels2))

F1_prop<-mean(Fs_toavg$propF1)

## create data for following forecast years using Fmult
# note that there is no interim year

for (i in 1:length(Fmult)) {
  fore_dat <- NULL
  aux_fore <- fore_dat_int
  # Forecast year (alternative Fs based on multipliers for all years)
  for(j in 1:Nfor){
    aux_fore$Year <- endyear+j  # CHECK THESE   
    aux_fore$F    <- Fmult[i]*c(F1_prop,0,0,1-F1_prop)*(2*Fsq)
    fore_dat <- rbind(fore_dat, aux_fore)
  }

  # input 
  fore$InputBasis <- 99 # options: 99 for F, 2 for Catch
  fore$ForeCatch  <- fore_dat # input ForeCatch(orF) data

  # IF WE WANT TO CHANGE THE R SCENARIO, UNCOMMENT 
  #   # Mean of historical recruitments 
  #   # fore$fcast_rec_option <- 3 # Mean of last x years (set in Fcast_years)
  #   # fore$Fcast_years[c(5:6)] <- c(-999,0) # whole series (i.e. not excluding latest year) - as it is currently
  #   # as we want the geometric mean, we calculate value relative to virgin rec
  #   fore$fcast_rec_option <- 2 #= value*(virgin recruitment)
  #   fore$fcast_rec_val <- gmrec/virg.rec # geomean / virgin rec
  
  # write all forecast files/scenarios
  r4ss::SS_writeforecast(fore, dir = stf_files.dir, file = paste0("forecast_",Fmult_names[i], ".ss"), 
                         overwrite = TRUE, verbose = FALSE)
    
}



# Prepare folders for STF -------------------------------------------------

# create forecast folder and subfolders

for (i in 1:length(Fmult)){
  # create subfolders for each Fmult scenario  
  dir.Fmult <- file.path(stf.dir, Fmult_names[i])
  dir.create(path = dir.Fmult, showWarnings = F, recursive = T)
  
  # copy the SS base files in every TAC subfolder 
  file.copy(file.path(run.dir, "starter.ss"), file.path(dir.Fmult, "starter.ss"),overwrite = T)
  file.copy(paste(run.dir, "ane.ctl", sep="/"), paste(dir.Fmult, "ane.ctl", sep="/"),overwrite = T)
  file.copy(paste(run.dir, "ane.dat", sep="/"), paste(dir.Fmult, "ane.dat", sep="/"),overwrite = T)	
  file.copy(paste(run.dir, "ss3.par", sep="/"), paste(dir.Fmult, "ss3.par", sep="/"),overwrite = T) # is this NEEDED???
  file.copy(paste(run.dir, "ss.exe", sep="/"), paste(dir.Fmult, "ss.exe", sep="/"),overwrite = T)
  file.copy(paste(run.dir, "wtatage.ss", sep="/"), paste(dir.Fmult, "wtatage.ss", sep="/"),overwrite = T)
  
  #copy the forecast file
  file.copy(paste(stf_files.dir, paste0("forecast_",Fmult_names[i],".ss") , sep="/"),
            paste(dir.Fmult, "forecast.ss", sep="/"),overwrite = T)

  # Edit "starter.ss" 
  starter <- SS_readstarter(file.path(dir.Fmult, "starter.ss"), verbose = FALSE)
  # - use ss.par estimates (0=use init values in control file; 1=use ss.par)
  starter$init_values_src <- 1 # to use the parameters estimates 
  # - turn off estimation for parameters entering after this phase
  # starter$last_estimation_phase <- 0
  # # Range of years for SD report (set to projection period)
  # starter$minyr_sdreport <- dtyr+1
  # starter$maxyr_sdreport <- dtyr+3
  # Save file
  SS_writestarter(starter, dir = dir.Fmult, verbose = FALSE, overwrite = TRUE)
  
  #fmax is 4, if the mutplier leads to F>4 change this setting in the ctl file
  if (max(Fmult[i]*F1_prop*Fsq*2)>4){
    ctl_file<-SS_readctl(file=paste(dir.Fmult, "ane.ctl", sep="/"),datlist = paste(dir.Fmult, "ane.dat", sep="/"))
    ctl_file$maxF<-round(Fmult[i]*F1_prop*Fsq*2+1,0)
    SS_writectl(ctl_file,outfile=paste(dir.Fmult, "ane.ctl", sep="/"),overwrite = TRUE)
    
  }
}


# CHECK THIS LINE THAT WAS IN HUGO'S CODE
###### ------- change starter.ss! line 31 -2 # max yr for sdreport outputs (-1 for endyr+1; -2 for endyr+Nforecastyrs) 


# Check one run -----------------------------------------------------------

# i <- 1
# r4ss::run(dir = file.path(stf.dir, paste0("Fmult_",Fmult[i])),
#             exe = "ss.exe",
#             #extras="-nohess",
#             skipfinished = FALSE, verbose=T)



# Run STF -----------------------------------------------------------------

# run in parallel
set.seed(135)
mc.cores <- 1 # number of cores to use
parallel::mclapply( file.path(stf.dir, Fmult_names),
                    r4ss::run, 
                    # extras = "-nohess", 
                    extras = "-maxfn 0 -phase 99", 
                    exe = "ss.exe", skipfinished = TRUE,
                    mc.cores=mc.cores)
                    
# alternatively, run in a loop using a single core

# for (i in 1:length(Fmult)) {
#   print(paste("Running scenario", Fmult[i]))
#   r4ss::run(dir = file.path(stf.dir, paste0("Fmult_",Fmult[i])),
#             exe = "ss.exe",
#             #extras="-nohess",
#             skipfinished = FALSE)
# }


