# Script information ------------------------------------------------------

# change in advice
# script for ane.27.8 usgin ss3 output directly

# Authors: Leire Citores (lcitores@azti.es)

# Date: 2025/11/26

# Load packages -----------------------------------------------------------

library(ggplot2)
library(r4ss)

theme_set(theme_bw())

#select the forecast folders used to give advice in ass.yr and ass.yr-1
outf<-SS_output("C:/USE/OneDrive - AZTI/GitHub/2025_ane.27.8_assessment/model/stf_SS3/stf/Fmult_0.941013079087865")
outf_prev<-SS_output("C:/USE/OneDrive - AZTI/GitHub/2024_ane.27.8_assessment/model/stf_SS3/stf/Fmult_1.30155431515672")
ass.yr<-2025
outs<-SSsummarize(list(outf, outf_prev))
pdf("output/stf/change_in_stf.pdf")
#ssb, rec, F
par(mfrow=c(3,1))
SSplotComparisons(outs, subplots=2, endyrvec = c(ass.yr, ass.yr-1)+1, legendlabels = c("ass ass.yr","ass ass.yr-1"),new=F,xlim=c(2020,ass.yr+1))
SSplotComparisons(outs, subplots=10,endyrvec = c(ass.yr+0.5, ass.yr-1), legendlabels = c("ass ass.yr","ass ass.yr-1"),new=F,xlim=c(2020,ass.yr+1))
SSplotComparisons(outs, subplots=8,endyrvec = c(ass.yr+0.5, ass.yr-1)+1, legendlabels = c("ass ass.yr","ass ass.yr-1"),new=F,xlim=c(2020,ass.yr+1))
par(mfrow=c(1,1))

#catch
catch <- outf$timeseries %>% 
  filter(Yr>ass.yr-5,Yr<ass.yr+2) %>% 
  select("Yr", "Seas", starts_with("dead(B)"),"Bio_all") 
catch$catch<-catch[,3]+catch[,4]
catch$ass.yr<-ass.yr

catch_prev <- outf_prev$timeseries %>% 
  filter(Yr>ass.yr-5,Yr<ass.yr+1) %>% 
  select("Yr", "Seas", starts_with("dead(B)"),"Bio_all") 
catch_prev$catch<-catch_prev[,3]+catch_prev[,4]
catch_prev$ass.yr<-ass.yr-1

catches<-rbind(catch,catch_prev)
p1<-ggplot(catches,aes(Yr,catch,col=factor(ass.yr)))+geom_line(lwd=1)+facet_wrap(~Seas,labeller = label_both)
p2<-ggplot(catches,aes(Yr,catch,col=factor(ass.yr)))+stat_summary(fun=sum,geom="line",lwd=1.2)
gridExtra::grid.arrange(p1,p2)

#Biomass at the beginning of the year
pp<-ggplot(subset(catches,Seas==1),aes(Yr,Bio_all,col=factor(ass.yr)))+
geom_line(lwd=1)+geom_point()+facet_wrap(~Seas,labeller = label_both)+ylab("Biomass1+ 1st January")
print(pp)

# fishery selectivities
p1<-ggplot(subset(outs$agesel,Factor=="Asel"&Yr>ass.yr-5&Yr<ass.yr+2&Fleet<3),aes(Yr,`1`,col=name))+geom_line()+geom_point()+
  facet_wrap(~Fleet,labeller = label_both)+ylab("sel age 1")+ylim(0,1)+geom_vline(xintercept=ass.yr,lty=2)

p2<-ggplot(subset(outs$agesel,Factor=="Asel"&Yr>ass.yr-5&Yr<ass.yr+2&Fleet<3),aes(Yr,`3`,col=name))+geom_line()+geom_point()+
  facet_wrap(~Fleet,labeller = label_both)+ylab("sel age 3")+ylim(0,1)+geom_vline(xintercept=ass.yr,lty=2)

gridExtra::grid.arrange(p1,p2)

# weight at age
outw<-melt(subset(outf$wtatage,year<ass.yr+2),id.vars = names(outf$wtatage)[1:6])
names(outw)[7:8]<-c("age","w")
outw$ass.yr<-ass.yr
outw_prev<-melt(subset(outf_prev$wtatage,year<ass.yr+1),id.vars = names(outf_prev$wtatage)[1:6])
names(outw_prev)[7:8]<-c("age","w")
outw_prev$ass.yr<-ass.yr-1

weights<-rbind(outw,outw_prev)
pp<-ggplot(subset(weights,year>ass.yr-2&fleet %in% c(c(1,2,-1))),aes(as.numeric(age)-1,w,col=factor(ass.yr)))+
  geom_line()+geom_point()+facet_grid(fleet~year+seas,labeller = label_both)+
  theme(legend.position = "bottom")
print(pp)
# Numbres at age
outn<-melt(subset(outf$natage,Yr<ass.yr+2&`Beg/Mid`=="B"),id.vars = names(outf$natage)[1:12])
names(outn)[13:14]<-c("age","N")
outn$ass.yr<-ass.yr
outn_prev<-melt(subset(outf_prev$natage,Yr<ass.yr+1&`Beg/Mid`=="B"),id.vars = names(outf_prev$natage)[1:12])
names(outn_prev)[13:14]<-c("age","N")
outn_prev$ass.yr<-ass.yr-1

natages<-rbind(outn,outn_prev)
pp<-ggplot(subset(natages,Yr>ass.yr-2),aes(as.numeric(age)-1,N,col=factor(ass.yr)))+
  geom_line()+geom_point()+facet_grid(Yr~Seas,labeller = label_both)
print(pp)

#Biomas at age
outb<-melt(subset(outf$batage,Yr<ass.yr+2&`Beg/Mid`=="B"),id.vars = names(outf$batage)[1:12])
names(outb)[13:14]<-c("age","B")
outb$ass.yr<-ass.yr
outb_prev<-melt(subset(outf_prev$batage,Yr<ass.yr+1&`Beg/Mid`=="B"),id.vars = names(outf_prev$batage)[1:12])
names(outb_prev)[13:14]<-c("age","B")
outb_prev$ass.yr<-ass.yr-1

batages<-rbind(outb,outb_prev)
pp<-ggplot(subset(batages,Yr>ass.yr-2),aes(as.numeric(age)-1,B,col=factor(ass.yr)))+
  geom_line()+geom_point()+facet_grid(Yr~Seas,labeller = label_both)
print(pp)

dev.off()