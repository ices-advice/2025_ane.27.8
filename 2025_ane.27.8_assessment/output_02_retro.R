# Script information ------------------------------------------------------

# generate the retro outputs from stock assessment for ane.27.8 

# Authors: Leire Citores (lcitores@azti.es), Leire Ibaibarriaga (libaibarriaga@azti.es)

# Date: 2024/11/21 modified: 2025/11/20

# Load packages -----------------------------------------------------------

library(icesTAF)
library(r4ss)
library(ss3diags)

# basic settings for ggplot

# theme_set(theme_bw(base_size = 25))

# Load ss run ---------------------------------------------------------------

dir_retro<-"model/ane.27.8_ss3_2025assessment_model"

retro_mods <- r4ss::SSgetoutput(dirvec = file.path(dir_retro, "retrospectives", paste0("retro", seq(0, -5, by = -1))), verbose = F)
retroSummary <- r4ss::SSsummarize(retro_mods, verbose = F,SpawnOutputUnits="biomass")


# Tables and plots ----------------------------------------------------------

SSmohnsrho(retroSummary, verbose = TRUE)

SSplotRetro(retroSummary, subplots = "SSB",forecast = F,add = T)
SSplotRetro(retroSummary, subplots = "SSB",forecast = T,add = T)

SSplotComparisons(retroSummary,subplots=9,endyrvec = c(2023,2022,2021,2020,2019,2018)+2,legendlabels = paste0("retro", seq(0, -5, by = -1)),new=F)
SSplotComparisons(retroSummary,subplots=8,endyrvec = c(2023,2022,2021,2020,2019,2018)+3,legendlabels = paste0("retro", seq(0, -5, by = -1)),new=F)
SSplotComparisons(retroSummary,subplots=11,endyrvec = c(2023,2022,2021,2020,2019,2018)+2,legendlabels = paste0("retro", seq(0, -5, by = -1)),new=F)

llks<-data.frame(t(retroSummary$likelihoods[,1:6]))
colnames(llks)<-retroSummary$likelihoods[,7]
llks$total_pos<-as.numeric(apply(llks,1,function(x){sum(abs(x[-c(1,4)]))}))
llks<-llks/llks$total_pos
llks$name<-paste0(0:5)
llks$Survey<-llks$Survey* -1

llks<-reshape2::melt(llks)
ggplot(llks,aes(as.numeric(name),value))+geom_point()+facet_wrap(~variable,scale="free")

p1<-ggplot(subset(retroSummary$agesel,Factor=="Asel"&Yr>1984&Fleet<4),aes(Yr,`1`,col=name))+geom_line()+
  facet_wrap(Seas~Fleet)+ylab("sel age 1")

p2<-ggplot(subset(retroSummary$agesel,Factor=="Asel"&Yr>1984&Fleet<4),aes(Yr,`3`,col=name))+geom_line()+
  facet_wrap(Seas~Fleet)+ylab("sel age 3")

gridExtra::grid.arrange(p1,p2)
pars<-data.frame(t(retroSummary$pars[c(20,69:76),1:6]))
colnames(pars)<-retroSummary$pars[c(20,69:76),7]
pars$name<-paste0(0:5)
pars<-reshape2::melt(pars)
ggplot(pars,aes(as.numeric(name),value))+geom_point()+facet_wrap(~variable,scale="free")


r4ss::sspar(mfrow = c(2, 2))
hind1<-SSplotHCxval(retroSummary, subplots = "cpue", add = TRUE,verbose = F)
retroSummary_comps <- SSretroComps(retro_mods)
print(hind1)
r4ss::sspar(mfrow = c(2, 2))
hind2<-SSplotHCxval(retroSummary_comps, subplots = "age", add = TRUE,verbose=F)
print(hind2)
r4ss::sspar(mfrow = c(1, 1))

# Save output -------------------------------------------------------------

save(retro_mods,retroSummary, file="output/retro/output_retro.RData")



### DONT RUN - FAST TRIALS
### PLOT PARAMS FOR ALL RETROS
# 
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # Assume summary_output <- SSsummarize(list(model1, model2, model3))
# params <- retroSummary$pars  # Wide format: Label, replist1, replist2, ...
# 
# # Reshape from wide to long
# params_long <- params %>%
#   pivot_longer(cols = starts_with("replist"),
#                names_to = "Model",
#                values_to = "Value")
# 
# # Check structure
# head(params_long)
# # Columns: Label, Model, Value, Parm_StDev (if available)
# 
# # Split into 4 equal groups based on unique Labels
# n <- length(unique(params_long$Label))
# param_labels <- unique(params_long$Label)
# groups <- cut(seq_along(param_labels), breaks = 4, labels = paste0("Group ", 1:4))
# 
# # Map Label -> group
# label_group <- data.frame(Label = param_labels, group = groups)
# params_long <- left_join(params_long, label_group, by = "Label")
# 
# # Generate plots for each group
# plots <- lapply(unique(params_long$group), function(g) {
#   ggplot(subset(params_long, group == g & Model %in% unique(params_long$Model)), aes(x = Label, y = Value, color = Model)) +
#     geom_point(size = 3, position = position_dodge(width = 0.5)) +
#     theme_minimal() +
#     labs(title = paste("SS3 Parameter Estimates -", g),
#          x = "Parameter",
#          y = "Estimated Value",
#          color = "Model") +
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
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # Assume summary_output <- SSsummarize(list(model1, model2, model3))
# params <- retroSummary$pars
# 
# params[,1:6]<-params[,1:6]-params[,1]
# 
# # Reshape from wide to long
# params_long <- params %>%
#   pivot_longer(cols = starts_with("replist"),
#                names_to = "Model",
#                values_to = "Value")
# 
# # Compute min and max for each parameter across models
# range_info <- params_long %>%
#   group_by(Label) %>%
#   summarise(min_val = min(Value, na.rm = TRUE),
#             max_val = max(Value, na.rm = TRUE),
#             .groups = "drop")
# 
# # Assign groups:
# # Group 1: parameters with all values between -1 and 1
# range_info$group <- ifelse(range_info$min_val >= -1 & range_info$max_val <= 1,
#                            "Group: [-1, 1]",
#                            NA)
# 
# # For remaining parameters, split by max_val into 3 bins
# range_info$group[is.na(range_info$group)] <- cut(range_info$max_val[is.na(range_info$group)],
#                                                  breaks = 3,
#                                                  labels = c("Low", "Medium", "High"))
# 
# # Merge back
# params_long <- left_join(params_long, range_info[, c("Label", "group")], by = "Label")
# 
# # Generate plots for each group
# plots <- lapply(unique(params_long$group), function(g) {
#   ggplot(filter(params_long, group == g&Model=="replist3"), aes(x = Label, y = round(Value,4), color = Model)) +
#     geom_point(size = 3, position = position_dodge(width = 0.5)) +
#     theme_minimal() +
#     labs(title = paste("SS3 Parameter Estimates -", g),
#          x = "Parameter",
#          y = "Estimated Value",
#          color = "Model") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# })
# 
# 
# 
# # Display plots
# plots[[1]]
# plots[[2]]
# plots[[3]]
# plots[[4]]
# 
