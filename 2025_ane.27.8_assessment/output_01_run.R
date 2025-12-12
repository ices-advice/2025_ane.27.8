# Script information ------------------------------------------------------

# generate the outputs from stock assessment for ane.27.8 

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/20

# Load packages -----------------------------------------------------------

library(icesTAF)
library(r4ss)
library(reshape2)
library(ss3diags)


# basic settings for ggplot

# theme_set(theme_bw(base_size = 25))

# Load ss run ---------------------------------------------------------------

out<-SS_output("model/ane.27.8_ss3_2025assessment_model")

# # reference points
# 
# refpts <- read.csv("data/refpts.csv")
# refpts

# Tables and plots ----------------------------------------------------------

SS_plots(out,html=TRUE,maxyr=2026)

# # residuals test
# sspar(mfrow = c(3, 3))
# pd1<-SSplotRunstest(out,subplot=c("age") ,add = TRUE,verbose = F,print = F)
# pd2<-SSplotRunstest(out,subplot=c("cpue") ,add = TRUE,verbose = F,print = F)
# sspar(mfrow = c(1, 1))
# sspar(mfrow = c(1, 2))
# rmse1<-SSplotJABBAres(out, subplots = "cpue", add = TRUE)
# rmse2<-SSplotJABBAres(out, subplots = "age", add = TRUE)
# sspar(mfrow = c(1, 1))
# SSplotComps(out,subplots = 24,print = T,kind="AGE")
# 
# selec_alter<-melt(subset(out$ageselex,Factor=="Asel")[,c("Fleet","Yr","0","1","2","3")],id.var=c("Yr","Fleet"),variable.name = "age")
# selec_alter<-subset(selec_alter,Fleet!=5)
# selec_alter$Fleet<-factor(selec_alter$Fleet,labels=c("Comm. fl 1","Comm. fl 2","AC. survey","DEPM survey"))
# 
# ps1<-ggplot(selec_alter,aes(as.numeric(age)-1,value,col=factor(Yr)))+geom_line()+geom_point()+facet_wrap(~Fleet)+xlab("age")+ylab("")+ggtitle(paste(" selectivities"))#+theme(legend.position = "none")
# 
# ps2<-ggplot(selec_alter,aes(Yr,value))+geom_line()+geom_point()+facet_grid(Fleet~age)+
#   geom_vline(xintercept=2010,col=2)+ylab("selec")+ggtitle("selectivities in time")+xlim(1987,2025)
# 
# print(ps1)
# 
# print(ps2)
# 
# PLOT PARAMS
# params <- out$parameters
# head(params)
# 
# # Example: plot estimated vs. prior
# plot(params$Value, params$Pr_type, main = "Parameter Estimates", xlab = "Value", ylab = "Prior Type")
# 
# 
# # Example: Plot all parameter estimates with confidence intervals
# ggplot(params, aes(x = Label, y = Value)) +
#   geom_point(color = "blue", size = 3) +
#   geom_errorbar(aes(ymin = Value - Parm_StDev, ymax = Value + Parm_StDev),
#                 width = 0.2, color = "gray50") +
#   theme_minimal() + facet_wrap(~(Phase),scale="free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels
# 
# 
# library(ggplot2)
# library(dplyr)
# 
# 
# # Split into 4 roughly equal groups
# n <- nrow(params)
# params$group <- cut(seq_len(n), breaks = 4, labels = paste0("Group ", 1:4))
# 
# # Loop through groups and plot
# plots <- lapply(unique(params$group), function(g) {
#   ggplot(filter(params, group == g), aes(x = Label, y = Value)) +
#     geom_point(color = "blue", size = 3) +
#     geom_errorbar(aes(ymin = Value - Parm_StDev, ymax = Value + Parm_StDev),
#                   width = 0.2, color = "gray50") +
#     theme_minimal() +
#     labs(title = paste("SS3 Parameter Estimates -", g),
#          x = "Parameter",
#          y = "Estimated Value") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# })
# 
# # Display plots
# plots[[1]]
# plots[[2]]
# plots[[3]]
# plots[[4]]
# 
# 
#        
# 
# numjitter <- 5
# jit.likes <- jitter(
#   dir = "model/ane.27.8_ss3_2025assessment_model", Njitter = numjitter,
#   jitter_fraction = 0.1, init_values_src = 1,exe="ss.exe"
# )
# #### Read in results using other r4ss functions
# # (note that un-jittered model can be read using keyvec=0:numjitter)
# profilemodels <- SSgetoutput(dirvec = "model/ane.27.8_ss3_2025assessment_model", keyvec = 1:numjitter, getcovar = TRUE)
# # summarize output
# profilesummary <- SSsummarize(profilemodels)
# # Likelihoods
# profilesummary[["likelihoods"]][1, ]
# # Parameters
# profilesummary[["pars"]]
# 
# SSplotComparisons(profilesummary)

# Save output -------------------------------------------------------------

save(out, file="output/run/output_run.RData")
