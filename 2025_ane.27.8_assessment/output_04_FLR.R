# Script info -------------------------------------------------------------

# code to read the SS3 output for anchovy in the Bay of Biscay and 
# transform it into FLR

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2025/06/03 

# Install and load libraries ----------------------------------------------

# library(remotes)
# remotes::install_github("r4ss/r4ss")
# remotes::install_github("flr/ss3om")
# remotes::install_github("flr/ss3om", dependencies=F)

# alternatively: 
# install.packages("ss3om", repos=structure(
#   c(CRAN="https://cloud.r-project.org/", FLR="https://flr-project.org/R")))

# USE THIS REPOSITORY: most up-to-date version
# install.packages('ss3om', repos = c('https://flr.r-universe.dev', 'https://cloud.r-project.org'))

library(r4ss)
library(ss3om)

library(dplyr)
library(tidyr)
library(ggplot2)

theme_set(theme_bw())

# Main directories --------------------------------------------------------

getwd()

# directory of the SS run

run.name <- "ane.27.8_ss3_2025assessment_model"
run.dir <- file.path("model", run.name)

# Read SS input and output ------------------------------------------------

# create a list for the SS input

dat  <- SS_readdat_3.30(file.path(run.dir, "ane.dat"))

# create a list object for the SS output

out <- SS_output(run.dir) 

# Years -------------------------------------------------------------------

# assessment year
ass.yr <- 2025

# start & end year
ini.yr <- out$startyr
end.yr <- out$endyr

# range of years in assessment 
yrs <- c(ini.yr:end.yr)

# Create FLStock by semester using ss3om ----------------------------------

# create an FLStock object (by semesters) from a model run using the function readFLSss3()

stk.sem <- readFLSss3(run.dir, name = "ane.27.8", desc = paste0("WGHANSA", ass.yr))

# complete/correct the stock object (by semesters)

mat(stk.sem)[ac(0),] <- 0
mat(stk.sem)[,,,2] <- 0 # maturity in the second semester is 0, no SSB in semester2
harvest.spwn(stk.sem)[,,,1] <- 4.5/6 # spawning at mid-May. Divided by 6 because harvest is by semester in stk.sem
harvest.spwn(stk.sem)[,,,2] <- 0
m.spwn(stk.sem)[,,,1] <- 4.5/6 # spawning at mid-May. Divided by 6 because m is by semester in stk.sem
m.spwn(stk.sem)[,,,2] <- 0
range(stk.sem)["minfbar"] <- 2 # this is apical F (maximum at age 2)
range(stk.sem)["maxfbar"] <- 2 # this is apical F (maximum at age 2) 

# natural mortality in SS3
#   subset(out$Natural_Mortality, Seas==1 & Era=="TIME") # natural mortality per year in 1st semester
#   subset(out$Natural_Mortality, Seas==2 & Era=="TIME") # natural mortality per year in 2nd semester
# note natural mortality for age 0 in semester 1 is set at 0.8 (although this is not used because R enters at the beginning of season 2)
# while natural mortality for age 0 in semester 2 is set at 2.17 (per year)
# we set manually the natural mortality for age 0 in semester 1 equal to zero

m(stk.sem)[1,,,1] <- 0 # natural mortality for age 0 in semester 1
m(stk.sem)

# plot FLStock by semester

plot(stk.sem[,,,1,,]) # semester 1
plot(stk.sem[,,,2,,]) # semester 2
plot(FLStocks(Sem1=stk.sem[,,,1,,], Sem2=stk.sem[,,,2,,])) # both semesters

# Check stk.sem -----------------------------------------------------------

# here we'll check that the R, SSB and F's in the seasonal FLstock are the same as in the SS3 output

# R. they are the same.

plot(rec(stk.sem)[,,,2,,])
rec(stk.sem)[,,,2]
aux <- data.frame(V0=as.numeric(rec(stk.sem)[,,,2]),
                  V1=subset(out$recruit, Yr %in% yrs, select="pred_recr"),
                  V2=subset(out$timeseries, Yr %in% yrs & Seas==2, select="Recruit_0"),
                  V3=subset(out$derived_quants, Label %in% paste0("Recr_", yrs), select="Value")
)
all(aux$V0==aux$V1)
all(aux$V0==aux$V2)
all(aux$V0==aux$V3)
aux

# SSB. they are the same.

aux <- data.frame(V0=as.numeric(ssb(stk.sem)[,,,1]),
                  V1=subset(out$timeseries, Yr %in% yrs & Seas==1, select="SpawnBio"),
                  V2=subset(out$derived_quants, Label %in% paste0("SSB_", yrs), select="Value"))
all(aux$V0==aux$V1)
all(aux$V0==aux$V2)
aux

# F semester1. they are the same.

aux <- data.frame(V1=as.numeric(harvest(stk.sem)[ac(2),,,1,,]),
                  V2=subset(out$fatage, Yr %in% yrs & Fleet==1 & Seas==1, select="2", drop=T)/2)
all(aux$V1==aux$V2) 
round(aux$V1-aux$V2, 5) # harvest is fatage/2 in the 1st semester
aux

# F semester2. they are the same EXCEPT FOR the last year

aux <- data.frame(V1=as.numeric(harvest(stk.sem)[ac(2),,,2,,]),
                  V2=subset(out$fatage, Yr %in% yrs & Fleet==2 & Seas==2, select="2", drop=T)/2)
all(aux$V1==aux$V2) 
round(aux$V1-aux$V2, 5) # harvest is fatage/2 in the 2nd semester, except LAST YEAR
aux

# F apical (per year). they are the same EXCEPT FOR the last year
# note that in derived quantities F_2008 and F_2009 are not given (catch==0 in both semesters)

fapical <- subset(out$derived_quants, Label %in% paste0("F_", yrs)) %>% 
  separate(Label, into = c("Variable","Yr"),sep = "_",remove = TRUE,extra="drop") %>% 
  subset(select=c("Yr", "Value"))
fapical <- left_join(data.frame(Yr=as.character(yrs)), fapical) # because F for 2008 and 2009 are not given in derived_quants

aux <-  data.frame(V1=as.numeric(seasonSums(harvest(stk.sem)[ac(2)])), # sum because harvest is already multiplied by the season duration
                   V2=fapical$Value)
all(aux$V1==aux$V2) # approx TRUE, except LAST YEAR
round(aux$V1-aux$V2,5)
aux
tail(aux)

# note that F in the last year is not equal. In the FLStock is computed as:
# harvest(stock) <- harvest(stock.n(stock), catch = catch.n(stock), m = m(stock))
# can be related to that?

# numbers-at-age 1st semester. they are the same

aux1 <- as.data.frame(stock.n(stk.sem)[,,,1,,])
aux2 <- subset(out$natage, Yr %in% yrs & Seas==1 & `Beg/Mid`=="B", select=c("Yr",ac(0:3))) %>% 
  pivot_longer(cols=-c("Yr"), values_to="values", names_to="age") %>% 
  rename(year=Yr)
names(aux2)

aux <- merge(aux1, aux2)
summary(aux$data - aux$values)

# numbers-at-age 2nd semester. they are the same

aux1 <- as.data.frame(stock.n(stk.sem)[,,,2,,])
aux2 <- subset(out$natage, Yr %in% yrs & Seas==2 & `Beg/Mid`=="B", select=c("Yr",ac(0:3))) %>% 
  pivot_longer(cols=-c("Yr"), values_to="values", names_to="age") %>% 
  rename(year=Yr)
names(aux2)

aux <- merge(aux1, aux2)
summary(aux$data - aux$values)

# catch-at-age (modelled). they are the same

aux1 <- as.data.frame(catch.n(stk.sem))
aux2 <- subset(out$catage, Yr %in% yrs, select=c("Yr","Seas","Fleet", ac(0:3))) %>% 
  pivot_longer(cols=-c("Yr","Fleet","Seas"), values_to="values", names_to="age") %>% 
  rename(year=Yr,
         season=Seas) %>% 
  group_by(year, season, age) %>% 
  summarise(values=sum(values))
names(aux2)

aux <- merge(aux1, aux2, by=c("age","year","season"))
aux[which(abs(aux$data - aux$values)>0), ] 
summary(aux$data - aux$values)

subset(aux, year==2024)

# Change manually the harvest of 2nd semester in last year ----------------

harvest(stk.sem)[,ac(end.yr),,2,,]

tmp <- subset(out$fatage, Yr %in% ac(end.yr) & Fleet==2 & Seas==2, select=ac(0:3))/2
tmp

harvest(stk.sem)[,ac(end.yr),,2,,] <- as.numeric(tmp)

# check again

aux <- data.frame(V1=as.numeric(harvest(stk.sem)[ac(2),,,2,,]),
                  V2=subset(out$fatage, Yr %in% yrs & Fleet==2 & Seas==2, select="2", drop=T)/2)
all(aux$V1==aux$V2) 
round(aux$V1-aux$V2, 5) # harvest is fatage/2 in the 2nd semester has been corrected by hand
aux

# FLStock annual ----------------------------------------------------------

# create FLStock (no seasons)
# with weighted=T, the weights are weighted means (weighted by numbers)
# otherwise, the weights will be simple averages

stk <- noseason(stk.sem, spwn.season=1, rec.season=2, weighted=TRUE)   

# complete/correct the stock object (annual)

harvest.spwn(stk) <- 4.5/12 # spawning at mid-May
m.spwn(stk) <- 4.5/12 # spawning at mid-May

# note that the annual rate of natural mortality for age 0 corresponds only to the 2nd semester (2.17/2)
m(stk) 
seasonSums(m(stk.sem))

# note that harvest is equal to the harvest in the first semester

harvest(stk) == harvest(stk.sem)[,,,1,,] 

# to reproduce the calculation of ssb, as ssb occurs in 1st season and harvest.spwn is 4.5/12
# we multiply the harvest rate for the season
# but, note than the resulting harvest is NOT real (not the sum of the seasonal harvests)

harvest(stk)[] <- harvest(stk.sem)[,,,1,,] *2

# plot yearly object

plot(stk)

# alternatively, the object with the true harvest rate would be

stk2 <- stk
harvest(stk2)[] <- seasonSums(harvest(stk.sem))

plot(FLStocks(ssbok=stk, fok=stk2))

# Check stk ---------------------------------------------------------------

# here we'll check that the R, SSB and F's in the annual FLstock are the same as in the SS3 output

# R

aux <- data.frame(V0=as.numeric(rec(stk)),
                  V1=subset(out$recruit, Yr %in% yrs, select="pred_recr", drop=T),
                  V2=subset(out$timeseries, Yr %in% yrs & Seas==2, select="Recruit_0", drop=T),
                  V3=subset(out$derived_quants, Label %in% paste0("Recr_", yrs), select="Value", drop=T)
                  )
all(aux$V0==aux$V1)
all(aux$V0==aux$V2)
all(aux$V0==aux$V3)
aux

# SSB

aux <- data.frame(V0=as.numeric(ssb(stk)),
                  V1=as.numeric(ssb(stk2)),
                  V2=subset(out$timeseries, Yr %in% yrs & Seas==1, select="SpawnBio", drop=T),
                  V3=subset(out$derived_quants, Label %in% paste0("SSB_", yrs), select="Value", drop=T))
all(aux$V0==aux$V2) # stk has correct SSBs
all(aux$V0==aux$V3)
all(aux$V1==aux$V2) # stk2 has wrong SSBs
all(aux$V1==aux$V3) 
aux

# F

aux <- data.frame(V0=as.numeric(harvest(stk)[ac(2)]), 
                  V1=as.numeric(harvest(stk2)[ac(2)]), 
                  V2=fapical$Value)
all(aux$V0==aux$V2) 
all(aux$V1==aux$V2) 
round(aux$V0-aux$V2,5) # stk has wrong harvest (it was twice the harvest of semester1)
round(aux$V1-aux$V2,5) # stk2 has correct harvest
aux
tail(aux)

# To convert from Fapical to Fbar -----------------------------------------

# SS3 returns Fapical (which is F at age 2)

range(stk2)
fbar(stk2)
harvest(stk2)[ac(2),]
fapical

aux <- data.frame(V1=as.numeric(fbar(stk2)),
           V2=as.numeric(harvest(stk2)[ac(2),]),
           V3=as.numeric(fapical$Value)) # they are the same, except LAST YEAR
aux

# if we want to change to fbar(1-2), change minfbar and maxfbar in the ranges 

stk3 <- stk2
range(stk3)["minfbar"] <- 1
range(stk3)["maxfbar"] <- 2
fbar(stk3)

data.frame(V1=c(quantMeans(harvest(stk2)[ac(1:2),])),
           V2=c(quantMeans(harvest(stk2)[ac(2),])),
           V3=c(quantMeans(harvest(stk2)[ac(1:3),])),
           V4=c(fbar(stk3)),
           V5=c(fbar(stk2)))

plot(FLStocks(stk2=stk2, stk3=stk3))

# Save the FLStock object -------------------------------------------------

# FLStock by semester
save(stk.sem, file="output/toFLR/stksem.RData")

# End of script -----------------------------------------------------------


