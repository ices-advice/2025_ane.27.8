# Script information ------------------------------------------------------

# outputs of the short-term forecast in SS for ane.27.8 

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
library(gridExtra)
library(flextable)
library(icesAdvice) # for ICES rounding


# Get STF outputs ---------------------------------------------------------

stf.dir   <- file.path("model","stf_SS3","stf")

load(file.path(stf.dir,"Fmult.Rdata"))

forecastModels <- SSgetoutput(dirvec=file.path(stf.dir, Fmult_names), getcovar=F)
names(forecastModels) <- Fmult_names

save(forecastModels, file=file.path(stf.dir, "forecast.RData"))

# Summarise outputs -------------------------------------------------------

forecastSummary <- SSsummarize(forecastModels)

#plot selectivities
ggplot(subset(forecastSummary$agesel[1:600,],Fleet %in% 1:2& Yr>1986),aes(Yr,`1`,col=name))+geom_point()+geom_line()+facet_wrap(~Fleet,scale="free")

#look for catch at semester proportions
c1<-subset(forecastModels[[11]]$timeseries,Seas==1)[,c("Yr","obs_cat:_1","dead(B):_1","F:_1")]
c2<-subset(forecastModels[[11]]$timeseries,Seas==2)[,c("Yr","obs_cat:_2","dead(B):_2","F:_2")]
cs<-cbind(c1,c2)
cs$ctot<-cs$`dead(B):_1`+cs$`dead(B):_2`
cs$Ftot<-(cs$`F:_1`+cs$`F:_2`)

cs$crel1<-cs$`dead(B):_1`/cs$ctot
cs$Frel1<-cs$`F:_1`/cs$Ftot
subset(cs,Yr %in% 2022:2027)
p1<-ggplot(subset(cs[,-1],Yr<2026),aes(Yr,crel1))+geom_point()+geom_line()
p2<-ggplot(subset(cs[,-1],Yr<2026),aes(Yr,Frel1))+geom_point()+geom_line()
gridExtra::grid.arrange(p1,p2)

# Initialise tables for results -------------------------------------------
endyear<-2025
Table_fmult <- setNames(data.frame(matrix(0,ncol=15,nrow=length(Fmult_names))), 
                        c("Fmult", 
                          paste0("Rec",endyear),
                          paste0("Rec",endyear+1),
                          paste0("SSB",endyear+1),
                          paste0("SSB",endyear+2),
                          paste0("SSBsd",endyear+1),
                          paste0("SSBsd",endyear+2),
                          paste0("pBlim", endyear+1),
                          paste0("pBlim", endyear+2),
                          paste0("pBpa", endyear+1),
                          paste0("pBpa", endyear+2),
                          paste0("F",endyear+1),
                          paste0("Catches",endyear+1),
                          paste0("F",endyear+2),
                          paste0("Catches",endyear+2)))

Table_fmult$Fmult <- round(Fmult,4)

# Fill-in tables with SSB -------------------------------------------------

SSB   <- as.data.frame(forecastSummary[["SpawnBio"]])

Table_fmult[,paste0("SSB",endyear+1)] <- SSB %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("SSB",endyear+2)] <- SSB %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()

#plot SSB
SSB_long<-filter(SSB,Yr<2027 & Yr>1986) %>%
  pivot_longer(cols = starts_with("Fmult_"), names_to = "Fmult", values_to = "SSB")
ggplot(SSB_long,aes(Yr,SSB,col=Fmult))+geom_line()+geom_point()+theme(legend.position = "bottom")


# Fill-in tables with SSBsd -----------------------------------------------

# EMPTY BY NOW. I DON'T KNOW WHY

SSBsd <- as.data.frame(forecastSummary[["SpawnBioSD"]])

Table_fmult[,paste0("SSBsd",endyear+1)] <- SSBsd %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("SSBsd",endyear+2)] <- SSBsd %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()

#plot SSBsd
SSBsd_long<-filter(SSBsd,Yr<2027 & Yr>1986) %>%
  pivot_longer(cols = starts_with("Fmult_"), names_to = "Fmult", values_to = "SSBsd")
ggplot(SSBsd_long,aes(Yr,SSBsd,col=Fmult))+geom_line()+geom_point()+theme(legend.position = "bottom")

SSBsd_long$CV<-SSBsd_long$SSBsd/SSB_long$SSB
ggplot(SSBsd_long,aes(Yr,CV,col=Fmult))+geom_line()+geom_point()+theme(legend.position = "bottom")


# Fill-in tables with SpawnBioLower ---------------------------------------

# check what assumption lognormal or normal is closer

SSBLower <- as.data.frame(forecastSummary[["SpawnBioLower"]])

Table_fmult[,paste0("SSBLower",endyear+1)] <- SSBLower %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("SSBLower",endyear+2)] <- SSBLower %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()

# Fill-in tables with SpawnBioUpper ---------------------------------------

# check what assumption lognormal or normal is closer

SSBUpper <- as.data.frame(forecastSummary[["SpawnBioUpper"]])

Table_fmult[,paste0("SSBUpper",endyear+1)] <- SSBUpper %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("SSBUpper",endyear+2)] <- SSBUpper %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()

# Fill in tables with P(SSB<Blim) -----------------------------------------

# # load reference points
refpts <- read.csv("boot/data/refpts_ane.27.8.csv", header=T)
refpts
Blim <- refpts$Blim
Bpa <- refpts$Bpa

# # - normal distribution
# 

SSBly<- SSB %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
SSBly.sd<- SSBsd %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("pBlim",endyear+1)]  <- pnorm(Blim, mean = SSBly, sd = SSBly.sd)
Table_fmult[,paste0("pBpa",endyear+1)]  <- pnorm(Bpa, mean = SSBly, sd = SSBly.sd)

# - lognormal distribution

logSSBly.mu    <- log(SSBly^2 / sqrt(SSBly^2 + SSBly.sd^2))
logSSBly.sigma <- sqrt(log(1 + SSBly.sd^2/SSBly^2))
Table_fmult[,paste0("pBlim_LN_",endyear+1)] <- plnorm(Blim, mean = logSSBly.mu, sd = logSSBly.sigma)

#endyear+2
SSBly<- SSB %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()
SSBly.sd<- SSBsd %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("pBlim",endyear+2)]  <- pnorm(Blim, mean = SSBly, sd = SSBly.sd)
Table_fmult[,paste0("pBpa",endyear+2)]  <- pnorm(Bpa, mean = SSBly, sd = SSBly.sd)


logSSBly.mu    <- log(SSBly^2 / sqrt(SSBly^2 + SSBly.sd^2))
logSSBly.sigma <- sqrt(log(1 + SSBly.sd^2/SSBly^2))
Table_fmult[,paste0("pBlim_LN_",endyear+2)] <- plnorm(Blim, mean = logSSBly.mu, sd = logSSBly.sigma)


# Fill in tables with F's -------------------------------------------------

Fvalue <- as.data.frame(forecastSummary[["Fvalue"]])

Table_fmult[,paste0("F",endyear+1)] <- Fvalue %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("F",endyear+2)] <- Fvalue %>% filter(Yr == endyear+2) %>% select(-Label, -Yr) %>% unlist()

FvalueSD <- as.data.frame(forecastSummary[["FvalueSD"]])
FvalueSD

# Fill in tables with R's -------------------------------------------------

Recr <- as.data.frame(forecastSummary[["recruits"]])

Table_fmult[,paste0("Rec",endyear)]   <- Recr %>% filter(Yr == endyear) %>% select(-Label, -Yr) %>% unlist()
Table_fmult[,paste0("Rec",endyear+1)] <- Recr %>% filter(Yr == endyear+1) %>% select(-Label, -Yr) %>% unlist()

Recr_long<-filter(Recr, Yr<2027 & Yr >1986) %>%
  pivot_longer(cols = starts_with("Fmult_"), names_to = "Fmult", values_to = "Recr")
ggplot(Recr_long,aes(Yr,Recr,col=Fmult))+geom_line()+geom_point()+theme(legend.position = "bottom")

# Fill in table with Rsd's ------------------------------------------------

recruitsSD <- as.data.frame(forecastSummary[["recruitsSD"]])
recruitsSD

# Fill in tables with catches ---------------------------------------------

for (i in 1:length(Fmult_names)){
  output = forecastModels[[i]]
  fltnms <- setNames(output$definitions$Fleet_name, output$fleet_ID)
  
  catch <- output$timeseries %>% 
    filter(Era == "FORE") %>% 
    select("Yr", "Seas", starts_with("dead(B)")) 
  names(catch) <- c('year', 'season', fltnms[1:nfleet])
  
  catch <- catch %>%  
    tidyr::pivot_longer(cols = -c("year","season"), names_to = 'fleet', values_to = 'value') %>% 
    select('year', 'season', 'fleet', 'value') %>% 
    mutate(year = as.factor(year))
  
  Cat <- catch %>% 
    group_by(year) %>% 
    summarise(cat=sum(value))
  
  Table_fmult[i, paste0("Catches",Cat$year)] <- Cat$cat
  Table_fmult[i, paste0("CatchS1",Cat$year)] <- subset(catch,season==1)$value[c(1,3)]
  
}


#additional columns related to the HCR
Table_fmult$ruleC<- 0.4*Table_fmult$SSB2026 -2600
Table_fmult$difCs<-Table_fmult$Catches2026-Table_fmult$ruleC
Table_fmult$NEWcatch2026<-Table_fmult$CatchS12026/0.8
Table_fmult$NEWdifCs<-Table_fmult$NEWcatch2026-Table_fmult$ruleC

# Save output table -------------------------------------------------------


write.csv(Table_fmult, file.path("output/stf", "table_Fmult.csv"), row.names = FALSE)


# Interpolation -------------------------------------------------------------

#plot values to be interpolated
ggplot(Table_fmult,aes(SSB2026,Catches2026))+geom_point()+geom_line()+
  geom_line(aes(SSB2026,ruleC))+theme_bw()+geom_vline(xintercept=89000,linetype=2)+
  geom_hline(yintercept=33000,linetype=2)
ggplot(Table_fmult,aes(Fmult,Catches2026))+geom_point()+geom_line()+
geom_vline(xintercept=1,linetype=2)+
  geom_hline(yintercept=33000,linetype=2)
ggplot(Table_fmult,aes(Fmult,pBlim2026))+geom_point()+geom_line()+
  geom_hline(yintercept=0.05,linetype=2)
ggplot(Table_fmult,aes(Fmult,pBpa2026))+geom_point()+geom_line()+
  geom_hline(yintercept=0.05,linetype=2)
ggplot(Table_fmult,aes(Fmult,SSB2026))+geom_point()+geom_line()+
  geom_hline(yintercept=c(Blim,refpts$Bpa),linetype=2)

#compute interpolated multipiers to run again stf with this values
FmultRule_interp<-approx(Table_fmult$difCs,Table_fmult$Fmult,0)
FmulBlim_interp<-approx(Table_fmult$SSB2026,Table_fmult$Fmult,refpts$Blim)
FmulBpa_interp<-approx(Table_fmult$SSB2026,Table_fmult$Fmult,refpts$Bpa)
FmulpBlim005_interp<-approx(Table_fmult$pBlim2026,Table_fmult$Fmult,0.05)
FmulpBpa005_interp<-approx(Table_fmult$pBpa2026,Table_fmult$Fmult,0.05)
Fmult_c33<-approx(Table_fmult$Catches2026,Table_fmult$Fmult,33000)

interpolated_Fmults<-c(FmultRule_interp$y,FmulBlim_interp$y,FmulBpa_interp$y,FmulpBlim005_interp$y,FmulpBpa005_interp$y,Fmult_c33$y)
interpolated_Fmults[is.na(interpolated_Fmults)]<-0

# run again model_03_stf-R with interpolated multipliers -------------------------------------------------------------

source("model_03_stf.R")

#run this script again (only until write.csv(...))

x <- scan("output_03_stf.R", what = character(),
         sep = '\n', skip = 0,
         encoding = 'UTF-8', quiet = TRUE) 
x<-x[1:grep("write.csv",x)]

eval.parent(parse(text=x))

# Final table ----------------------------------------------------------------

#give names to computed options
Table_fmult$name<- c(paste("Fsq x",Table_fmult$Fmult[1:47]),"MP_nolim","SSB2026=Blim","SSB2026=Bpa","p(SSB2026<Blim)=0.05","p(SSB2026<Bpa)=0.05","Catch2026=33000")
#Table_fmult[,c("name","Fmult","F2025","Catches2025","SSB2025","pBlim2025","pBpa2025")]


#compute change %

Table_fmult$TAC_prev<-30663
Table_fmult[,paste0("SSB",endyear)] <- SSB %>% filter(Yr == endyear) %>% select(-Label, -Yr) %>% unlist()
Table_fmult$SSBchange<-100*(Table_fmult$SSB2026-Table_fmult$SSB2025)/Table_fmult$SSB2025
Table_fmult$TACchange<-100*(Table_fmult$Catches2026-Table_fmult$TAC_prev)/Table_fmult$TAC_prev
Table_fmult$ADVchange<-100*(Table_fmult$Catches2026-Table_fmult$TAC_prev)/Table_fmult$TAC_prev

## MP rule: max catch 33000
if( subset(Table_fmult,name=="MP_nolim")$Catches2026 > 33000 ){
  Table_fmult<-rbind(subset(Table_fmult,name=="Catch2026=33000"),Table_fmult)
  Table_fmult$name[1]<-"MP"
}

Final_table<-subset(Table_fmult,
  name %in% c("Fsq x 0","Fsq x 0.5","Fsq x 1","MP","Catch2026=33000","SSB2026=Blim",
              "SSB2026=Bpa","p(SSB2026<Blim)=0.05","p(SSB2026<Bpa)=0.05"))[,
          c("name","Fmult","F2026","Catches2026","SSB2026","pBlim2026","SSBchange","TACchange","ADVchange")]

# ices rounding

Final_table[,-1] <- Final_table[,-1] %>% 
  mutate(across(starts_with("F"), ~ as.numeric(icesRound(.x))), 
         #across(starts_with("p"), ~ as.numeric(icesRound(.x))), 
         across(starts_with("p"), ~ round(.x,4)),
         across(!starts_with("F") & !starts_with("p"), ~ round(.x)))


ft<-colformat_double(
  x = flextable(Final_table[c(1:4,9,8,7,6,5),c(1,4,5,3,7:9,6)]),
  big.mark = "",digits=4, na_str = "N/A"
)

flextable(Final_table[c(1:4,9,8,7,6,5),c(1,2,4,5,3,7:9,6)]) %>% save_as_docx( path = "output/stf/stf_final_table.docx")


save.image(file="output/stf/out_stf.RData")

# STF plot ----------------------------------------------------------------


# #huge catch in last scenarios, due to differences between catch w and poplation w
# outf<-forecastModels$Fmult_16.8871686556604
# subset(outf$natage,Yr==2026)[c(1,3),c("0","1","2","3")]*subset(outf$wtatage,fleet %in% c(4)&year==2026)[,c("0","1","2","3")]-
# subset(outf$catage,Yr==2026)[c(1,4),c("0","1","2","3")]*subset(outf$wtatage,fleet %in% c(1)&year==2026)[,c("0","1","2","3")]
# 
# subset(outf$natage,Yr==2026)[c(1,3),c("0","1","2","3")]-
# subset(outf$catage,Yr==2026)[c(1,4),c("0","1","2","3")]
# 
# sum((subset(outf$natage,Yr==2026)[c(1,3),c("0","1","2","3")]*subset(outf$wtatage,fleet %in% c(4)&year==2026)[,c("0","1","2","3")])[1,])
#   
# subset(outf$batage,Yr==2026)
# subset(outf$fatage,Yr==2026)

