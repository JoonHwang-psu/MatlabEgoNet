#####################################
#                                   #
#   Matlab Ego network analysis     #
#                                   #
#####################################

##### Sensitivity check with 20NN

library(readr)
library(reshape2)
library(tidyr)
library(netcom)
library(ineq)
library(ggplot2)
library(interplot)
library(stargazer)
library(gridExtra)
library(lattice)
library(dplyr)
library(ggeffects)
library(sjPlot)
library(forcats)
library(convey)
library(arm)
library(vcd)
library(MASS)
library(lmtest)
library(pscl)
library(patchwork)
library(broom)

setwd("C:/Users/junhw/OneDrive - The Pennsylvania State University/Graduate Study/0. DissertationProject/2. Chapter 2. EgoNet + Spatial/Data for dissertation/")

##### 1. creating contextual variables #####

Dist_20NN<-read_csv("MatlabEgoNet_HH_Dist_40NN_NoDuplicateRank.csv")
Dist_20NN_HHID<-Dist_20NN[,3:7]
Dist_20NN_wealth<-Dist_20NN[,c(3,5,6,10)]
Dist_20NN_distance<-Dist_20NN[,c(3,6,11)]

Dist_20NN_wide_HHID<-pivot_wider(Dist_20NN_HHID, names_from=NEAR_RANK, values_from=NEAR_FID)
Dist_20NN_wide_wealth<-pivot_wider(Dist_20NN_wealth, names_from=NEAR_RANK, values_from=HH_to_wealth)
Dist_20NN_wide_distance<-pivot_wider(Dist_20NN_distance, names_from=NEAR_RANK, values_from=NEAR_DIST)

n_HH<-nrow(Dist_20NN_wide_HHID)

Dist_20NN_wide_wealth$wealth_gini<-0
Dist_20NN_wide_distance$distance_mean<-0

# Here I only calculate neighborhood inequality among 20 nearest neighbors
for (i in 1:n_HH){
  Dist_20NN_wide_wealth$wealth_gini[i]<-ineq(as.numeric(Dist_20NN_wide_wealth[i,2:22]), type= "Gini", na.rm=TRUE)
  Dist_20NN_wide_distance$distance_mean[i]<-mean(as.numeric(Dist_20NN_wide_distance[i,2:22]), na.rm=TRUE)
}

write.csv(Dist_20NN_wide_HHID, file="MatlabEgoNet_20NN_HHID.csv")
write.csv(Dist_20NN_wide_wealth, file="MatlabEgoNet_20NN_wealth.csv")
write.csv(Dist_20NN_wide_distance, file="MatlabEgoNet_20NN_distance.csv")


summary(Dist_20NN_wide_wealth$wealth_gini)
hist(Dist_20NN_wide_wealth$wealth_gini)
summary(Dist_20NN_wide_distance$distance_mean)
hist(Dist_20NN_wide_distance$distance_mean)

context_20NN<-merge(Dist_20NN_wide_HHID, Dist_20NN_wide_distance[,c(1,42)], by="HH_from")
context_20NN<-merge(context_20NN, Dist_20NN_wide_wealth[,c(1,43)], by="HH_from")

write.csv(context_20NN, file = "MatlabEgoNet_20NN_context.csv", row.names = TRUE)

##### 2. connecting ego networks to contextual variables

EgoNet_total <- read_csv("MatlabEgoNet_250320.csv")
colnames(context_20NN)[1]<-c("HHID")
colnames(EgoNet_total)[1]<-c("HHID")

EgoNet_total_context <- merge(EgoNet_total, context_20NN[,c(1,44,45)], by="HHID")
EgoNet_total_context$wealth1000<-EgoNet_total_context$wealth_total/1000
EgoNet_total_context$age_2018_sq <- EgoNet_total_context$age_2018^2
EgoNet_total_context$sex <- as.factor(EgoNet_total_context$sex)

EgoNet_total_context <- EgoNet_total_context[EgoNet_total_context$distance_mean<3000,]
summary(EgoNet_total_context$distance_mean)
hist(EgoNet_total_context$distance_mean)

summary(EgoNet_total_context$wealth_gini)
quantile(EgoNet_total_context$wealth_gini, probs = c(0.10, 0.90), na.rm = TRUE) # 10% = 0.2296791, 90% = 0.3425539 


EgoNet_total_context$RecentEconRisk <- EgoNet_total_context$BusLossRecent+EgoNet_total_context$UnemployRecent+EgoNet_total_context$LgDebtRecent
EgoNet_total_context$RecentEnvRisk <- EgoNet_total_context$CropFailRecent+EgoNet_total_context$RiverErodeRecent+EgoNet_total_context$NatDisRecent
EgoNet_total_context$RecentHealthRisk <- EgoNet_total_context$SickSelfRecent+EgoNet_total_context$SickFamRecent
EgoNet_total_context$RecentSocRisk <- EgoNet_total_context$KinDispRecent+EgoNet_total_context$NeighDispRecent

EgoNet_total_context$wealth_gini0.1 <- 10*EgoNet_total_context$wealth_gini

##### 3. Data analysis

### please change the date of the folder for workind directory so that the results from different analysis can be saved separately

setwd("C:/Users/junhw/OneDrive - The Pennsylvania State University/Graduate Study/0. DissertationProject/2. Chapter 2. EgoNet + Spatial/Data for dissertation/risk experience")

##### 3.1. network size
prtn.main <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(prtn.main)

prtn.main.pois <- glm(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(prtn.main.pois)

dispersion <- sum(residuals(prtn.main.pois, type="pearson")^2) / df.residual(prtn.main.pois)
print(dispersion) ## overdispersion detected (=1.21)

lrtest(prtn.main.pois, prtn.main) ### this indicates that negative binomial model fits better

prtn.econrisk <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, data=EgoNet_total_context)
summary(prtn.econrisk)
prtn.healthrisk <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, data=EgoNet_total_context)
summary(prtn.healthrisk)
prtn.envrisk <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, data=EgoNet_total_context)
summary(prtn.envrisk)
prtn.socrisk <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, data=EgoNet_total_context)
summary(prtn.socrisk)

(prtn.econrisk.plot<-plot_model(prtn.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.econrisk$y)),ceiling(max(prtn.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.healthrisk.plot<-plot_model(prtn.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.healthrisk$y)),ceiling(max(prtn.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.envrisk.plot<-plot_model(prtn.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.envrisk$y)),ceiling(max(prtn.envrisk$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.socrisk.plot<-plot_model(prtn.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.socrisk$y)),ceiling(max(prtn.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(prtn.econrisk.plot$data)
print(prtn.healthrisk.plot$data)
print(prtn.envrisk.plot$data)
print(prtn.socrisk.plot$data)

summarize_model <- function(model, model_name) {
  coef <- coef(model)  
  exp_coef <- exp(coef)  
  ci <- confint(model) 
  exp_ci <- exp(ci) 
  
  result <- data.frame(
    Model = model_name,
    Variable = names(coef),
    Coefficient = coef,
    IRR_OR = exp_coef,
    CI_Lower = exp_ci[, 1],
    CI_Upper = exp_ci[, 2]
  )
  return(result)
}

models <- list(prtn.econrisk=prtn.econrisk, prtn.healthrisk=prtn.healthrisk, prtn.envrisk=prtn.envrisk, prtn.socrisk=prtn.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "Netsize_20NN.csv", row.names = FALSE)

# interaction with sex

prtn.sex <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(prtn.sex)
prtn.sex.pois <- glm(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(prtn.sex.pois)

dispersion <- sum(residuals(prtn.sex.pois, type="pearson")^2) / df.residual(prtn.sex.pois)
print(dispersion) ## overdispersion detected (=1.21)

lrtest(prtn.sex.pois, prtn.sex) ### this indicates that negative binomial model fits better

prtn.econ.sex <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, data=EgoNet_total_context)
summary(prtn.econ.sex)
prtn.health.sex <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, data=EgoNet_total_context)
summary(prtn.health.sex)
prtn.env.sex <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, data=EgoNet_total_context)
summary(prtn.env.sex)
prtn.soc.sex <- glm.nb(netsize~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, data=EgoNet_total_context)
summary(prtn.soc.sex)

(prtn.econ.sex.plot <-plot_model(prtn.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.econ.sex$y)),ceiling(max(prtn.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.health.sex.plot <-plot_model(prtn.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.health.sex$y)),ceiling(max(prtn.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.env.sex.plot <-plot_model(prtn.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.env.sex$y)),ceiling(max(prtn.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(prtn.soc.sex.plot <-plot_model(prtn.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(prtn.soc.sex$y)),ceiling(max(prtn.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

models <- list("Economic risk * Gini * Male"=prtn.econ.sex, "Health risk * Gini * Male"=prtn.health.sex, "Env. risk * Gini * Male"=prtn.env.sex, "Social risk * Gini * Male"=prtn.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3wayNetsize_20NN.csv", row.names = FALSE)

##### 3.2. reciprocal partners
mutual.main <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(mutual.main)

mutual.main.pois <- glm(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(mutual.main.pois)

dispersion <- sum(residuals(mutual.main.pois, type="pearson")^2) / df.residual(mutual.main.pois)
print(dispersion) ## overdispersion detected (=1.52)

lrtest(mutual.main.pois, mutual.main) # NB works better 

mutual.econrisk <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, data=EgoNet_total_context)
summary(mutual.econrisk)
mutual.healthrisk <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, data=EgoNet_total_context)
summary(mutual.healthrisk)
mutual.envrisk <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, data=EgoNet_total_context)
summary(mutual.envrisk)
mutual.socrisk <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, data=EgoNet_total_context)
summary(mutual.socrisk)

(mutual.econrisk.plot<-plot_model(mutual.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.econrisk$y)),ceiling(max(mutual.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.healthrisk.plot<-plot_model(mutual.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.healthrisk$y)),ceiling(max(mutual.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.envrisk.plot<-plot_model(mutual.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.envrisk$y)),ceiling(max(mutual.envrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.socrisk.plot<-plot_model(mutual.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.socrisk$y)),ceiling(max(mutual.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(mutual.econrisk.plot$data)
print(mutual.healthrisk.plot$data)
print(mutual.envrisk.plot$data)
print(mutual.socrisk.plot$data)

models <- list(mutual.econrisk=mutual.econrisk, mutual.healthrisk=mutual.healthrisk, mutual.envrisk=mutual.envrisk, mutual.socrisk=mutual.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Model,  # Each model becomes a column
    values_from = IRR_CI  # Fill cells with IRR and CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "Reciprocity_20NN.csv", row.names = FALSE)

# interaction with sex

mutual.sex <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(mutual.sex)
mutual.sex.pois <- glm(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(mutual.sex.pois)

dispersion <- sum(residuals(mutual.sex.pois, type="pearson")^2) / df.residual(mutual.sex.pois)
print(dispersion) ## overdispersion detected (=1.52)

lrtest(mutual.sex.pois, mutual.sex) ### NB fits better

mutual.econ.sex <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, data=EgoNet_total_context)
summary(mutual.econ.sex)
mutual.health.sex <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, data=EgoNet_total_context)
summary(mutual.health.sex)
mutual.env.sex <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, data=EgoNet_total_context)
summary(mutual.env.sex)
mutual.soc.sex <- glm.nb(mutual~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, data=EgoNet_total_context)
summary(mutual.soc.sex)

(mutual.econ.sex.plot <-plot_model(mutual.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.econ.sex$y)),ceiling(max(mutual.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.health.sex.plot <-plot_model(mutual.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.health.sex$y)),ceiling(max(mutual.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.env.sex.plot <-plot_model(mutual.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.env.sex$y)),ceiling(max(mutual.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(mutual.soc.sex.plot <-plot_model(mutual.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(mutual.soc.sex$y)),ceiling(max(mutual.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(mutual.econ.sex.plot$data)
print(mutual.health.sex.plot$data)
print(mutual.env.sex.plot$data)
print(mutual.soc.sex.plot$data)

models <- list("Economic risk * Gini * Male"=mutual.econ.sex, "Health risk * Gini * Male"=mutual.health.sex, "Env. risk * Gini * Male"=mutual.env.sex, "Social risk * Gini * Male"=mutual.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3wayReciprocity_20NN.csv", row.names = FALSE)


##### 3.3. long-distance partners
long.main <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(long.main)

long.main.pois <- glm(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(long.main.pois)

dispersion <- sum(residuals(long.main.pois, type="pearson")^2) / df.residual(long.main.pois)
print(dispersion) ## overdispersion detected (=2.04)

lrtest(long.main.pois, long.main) # NB fits better 

long.econrisk <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, data=EgoNet_total_context)
summary(long.econrisk)
long.healthrisk <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, data=EgoNet_total_context)
summary(long.healthrisk)
long.envrisk <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, data=EgoNet_total_context)
summary(long.envrisk)
long.socrisk <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, data=EgoNet_total_context)
summary(long.socrisk)

(long.econrisk.plot<-plot_model(long.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.econrisk$y)),ceiling(max(long.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.healthrisk.plot<-plot_model(long.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(long.healthrisk$y)),ceiling(max(long.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.envrisk.plot<-plot_model(long.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.envrisk$y)),ceiling(max(long.envrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.socrisk.plot<-plot_model(long.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(long.socrisk$y)),ceiling(max(long.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(long.econrisk.plot$data)
print(long.healthrisk.plot$data)
print(long.envrisk.plot$data)
print(long.socrisk.plot$data)

models <- list(long.econrisk=long.econrisk, long.healthrisk=long.healthrisk, long.envrisk=long.envrisk, long.socrisk=long.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model, 
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "LongDist_20NN.csv", row.names = FALSE)

# interaction with sex

long.sex <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(long.sex)
long.sex.pois <- glm(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(long.sex.pois)

dispersion <- sum(residuals(long.sex.pois, type="pearson")^2) / df.residual(long.sex.pois)
print(dispersion) ## overdispersion detected (=2.03)

lrtest(long.sex.pois, long.sex) ### NB fits better

long.econ.sex <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, data=EgoNet_total_context)
summary(long.econ.sex)
long.health.sex <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, data=EgoNet_total_context)
summary(long.health.sex)
long.env.sex <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, data=EgoNet_total_context)
summary(long.env.sex)
long.soc.sex <- glm.nb(long~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, data=EgoNet_total_context)
summary(long.soc.sex)

(long.econ.sex.plot <-plot_model(long.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.econ.sex$y)),ceiling(max(long.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.health.sex.plot <-plot_model(long.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.health.sex$y)),ceiling(max(long.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.env.sex.plot <-plot_model(long.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.env.sex$y)),ceiling(max(long.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(long.soc.sex.plot <-plot_model(long.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(long.soc.sex$y)),ceiling(max(long.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

models <- list("Economic risk * Gini * Male"=long.econ.sex, "Health risk * Gini * Male"=long.health.sex, "Env. risk * Gini * Male"=long.env.sex, "Social risk * Gini * Male"=long.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3wayLongDist_20NN.csv", row.names = FALSE)

##### 3.4. local partners
local.main <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(local.main)

local.main.pois <- glm(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(local.main.pois)

dispersion <- sum(residuals(local.main.pois, type="pearson")^2) / df.residual(local.main.pois)
print(dispersion) ## overdispersion detected (=1.13)

lrtest(local.main.pois, local.main) # NB fits better 

local.econrisk <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, data=EgoNet_total_context)
summary(local.econrisk)
local.healthrisk <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, data=EgoNet_total_context)
summary(local.healthrisk)
local.envrisk <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, data=EgoNet_total_context)
summary(local.envrisk)
local.socrisk <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, data=EgoNet_total_context)
summary(local.socrisk)

(local.econrisk.plot<-plot_model(local.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.econrisk$y)),ceiling(max(local.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.healthrisk.plot<-plot_model(local.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(local.healthrisk$y)),ceiling(max(local.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.envrisk.plot<-plot_model(local.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.envrisk$y)),ceiling(max(local.envrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.socrisk.plot<-plot_model(local.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(local.socrisk$y)),ceiling(max(local.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(local.econrisk.plot$data)
print(local.healthrisk.plot$data)
print(local.envrisk.plot$data)
print(local.socrisk.plot$data)

models <- list(local.econrisk=local.econrisk, local.healthrisk=local.healthrisk, local.envrisk=local.envrisk, local.socrisk=local.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),  # Round IRR to 2 decimal places
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Model,  # Each model becomes a column
    values_from = IRR_CI  # Fill cells with IRR and CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "Local_20NN.csv", row.names = FALSE)

# interaction with sex

local.sex <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(local.sex)
local.sex.pois <- glm(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(local.sex.pois)

dispersion <- sum(residuals(local.sex.pois, type="pearson")^2) / df.residual(local.sex.pois)
print(dispersion) ## overdispersion detected (=1.13)

lrtest(local.sex.pois, local.sex) ### NB fits better

local.econ.sex <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, data=EgoNet_total_context)
summary(local.econ.sex)
local.health.sex <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, data=EgoNet_total_context)
summary(local.health.sex)
local.env.sex <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, data=EgoNet_total_context)
summary(local.env.sex)
local.soc.sex <- glm.nb(local~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, data=EgoNet_total_context)
summary(local.soc.sex)

(local.econ.sex.plot <-plot_model(local.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.econ.sex$y)),ceiling(max(local.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.health.sex.plot <-plot_model(local.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.health.sex$y)),ceiling(max(local.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.env.sex.plot <-plot_model(local.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.env.sex$y)),ceiling(max(local.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(local.soc.sex.plot <-plot_model(local.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(local.soc.sex$y)),ceiling(max(local.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

models <- list("Economic risk * Gini * Male"=local.econ.sex, "Health risk * Gini * Male"=local.health.sex, "Env. risk * Gini * Male"=local.env.sex, "Social risk * Gini * Male"=local.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3wayLocal_20NN.csv", row.names = FALSE)

##### 3.5. support-receiving ties
supp.in.main <- glm.nb(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(supp.in.main)

supp.in.main.pois <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(supp.in.main.pois)

dispersion <- sum(residuals(supp.in.main.pois, type="pearson")^2) / df.residual(supp.in.main.pois)
print(dispersion) # no evidence of overdispersion (=1.02)

lrtest(supp.in.main.pois, supp.in.main) # no support that NB fits better

supp.in.econrisk <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, family="poisson", data=EgoNet_total_context)
summary(supp.in.econrisk)
supp.in.healthrisk <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, family="poisson", data=EgoNet_total_context)
summary(supp.in.healthrisk)
supp.in.envrisk <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, family="poisson", data=EgoNet_total_context)
summary(supp.in.envrisk)
supp.in.socrisk <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, family="poisson", data=EgoNet_total_context)
summary(supp.in.socrisk)

(supp.in.econrisk.plot<-plot_model(supp.in.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.econrisk$y)),ceiling(max(supp.in.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.healthrisk.plot<-plot_model(supp.in.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.healthrisk$y)),ceiling(max(supp.in.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.envrisk.plot<-plot_model(supp.in.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.envrisk$y)),ceiling(max(supp.in.envrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.socrisk.plot<-plot_model(supp.in.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.socrisk$y)),ceiling(max(supp.in.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(supp.in.econrisk.plot$data)
print(supp.in.healthrisk.plot$data)
print(supp.in.envrisk.plot$data)
print(supp.in.socrisk.plot$data)

models <- list(supp.in.econrisk=supp.in.econrisk, supp.in.healthrisk=supp.in.healthrisk, supp.in.envrisk=supp.in.envrisk, supp.in.socrisk=supp.in.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),  # Round IRR to 2 decimal places
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Model,  # Each model becomes a column
    values_from = IRR_CI  # Fill cells with IRR and CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "SuppIn_20NN.csv", row.names = FALSE)

# interaction with sex

supp.in.sex <- glm.nb(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(supp.in.sex)
supp.in.sex.pois <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(supp.in.sex.pois)

dispersion <- sum(residuals(supp.in.sex.pois, type="pearson")^2) / df.residual(supp.in.sex.pois)
print(dispersion) # no evidence of overdispersion (=1.02)

lrtest(supp.in.sex.pois, supp.in.sex) # no support that NB works better

supp.in.econ.sex <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, family="poisson", data=EgoNet_total_context)
summary(supp.in.econ.sex)
supp.in.health.sex <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, family="poisson", data=EgoNet_total_context)
summary(supp.in.health.sex)
supp.in.env.sex <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, family="poisson", data=EgoNet_total_context)
summary(supp.in.env.sex)
supp.in.soc.sex <- glm(support_in~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, family="poisson", data=EgoNet_total_context)
summary(supp.in.soc.sex)

(supp.in.econ.sex.plot <-plot_model(supp.in.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.econ.sex$y)),ceiling(max(supp.in.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.health.sex.plot <-plot_model(supp.in.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.health.sex$y)),ceiling(max(supp.in.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.env.sex.plot <-plot_model(supp.in.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.env.sex$y)),ceiling(max(supp.in.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.in.soc.sex.plot <-plot_model(supp.in.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.in.soc.sex$y)),ceiling(max(supp.in.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(supp.in.econ.sex.plot$data)
print(supp.in.health.sex.plot$data)
print(supp.in.env.sex.plot$data)
print(supp.in.soc.sex.plot$data)


models <- list("Economic risk * Gini * Male"=supp.in.econ.sex, "Health risk * Gini * Male"=supp.in.health.sex, "Env. risk * Gini * Male"=supp.in.env.sex, "Social risk * Gini * Male"=supp.in.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3waySuppIn_20NN.csv", row.names = FALSE)

##### 3.6. support-providing ties
supp.out.main <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1, data=EgoNet_total_context)
summary(supp.out.main)

supp.out.main.pois <- glm(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1, family="poisson",data=EgoNet_total_context)
summary(supp.out.main.pois)

dispersion <- sum(residuals(supp.out.main.pois, type="pearson")^2) / df.residual(supp.out.main.pois)
print(dispersion) # overdispersion detected (=1.42)

lrtest(supp.out.main.pois, supp.out.main) # NB fits better

supp.out.econrisk <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk, data=EgoNet_total_context)
summary(supp.out.econrisk)
supp.out.healthrisk <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk, data=EgoNet_total_context)
summary(supp.out.healthrisk)
supp.out.envrisk <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk, data=EgoNet_total_context)
summary(supp.out.envrisk)
supp.out.socrisk <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk, data=EgoNet_total_context)
summary(supp.out.socrisk)

(supp.out.econrisk.plot<-plot_model(supp.out.econrisk,type="pred",terms=c("RecentEconRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.econrisk$y)),ceiling(max(supp.out.econrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.healthrisk.plot<-plot_model(supp.out.healthrisk,type="pred",terms=c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.healthrisk$y)),ceiling(max(supp.out.healthrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.envrisk.plot<-plot_model(supp.out.envrisk,type="pred",terms=c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.envrisk$y)),ceiling(max(supp.out.envrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.socrisk.plot<-plot_model(supp.out.socrisk,type="pred",terms=c("RecentSocRisk","wealth_gini0.1[2.5,3.5]"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,2,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.socrisk$y)),ceiling(max(supp.out.socrisk$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

print(supp.out.econrisk.plot$data)
print(supp.out.healthrisk.plot$data)
print(supp.out.envrisk.plot$data)
print(supp.out.socrisk.plot$data)

models <- list(supp.out.econrisk=supp.out.econrisk, supp.out.healthrisk=supp.out.healthrisk, supp.out.envrisk=supp.out.envrisk, supp.out.socrisk=supp.out.socrisk)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),  # Round IRR to 2 decimal places
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Model,  # Each model becomes a column
    values_from = IRR_CI  # Fill cells with IRR and CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "SuppOut_251014.csv", row.names = FALSE)

# interaction with sex

supp.out.sex <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*sex, data=EgoNet_total_context)
summary(supp.out.sex)
supp.out.sex.pois <- glm(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentEnvRisk+RecentHealthRisk+RecentSocRisk+wealth_gini0.1*sex, family="poisson",data=EgoNet_total_context)
summary(supp.out.sex.pois)

dispersion <- sum(residuals(supp.out.sex.pois, type="pearson")^2) / df.residual(supp.out.sex.pois)
print(dispersion) # overdispersion detected (=1.42)

lrtest(supp.out.sex.pois, supp.out.sex) ### this indicates that negative binomial model fits better

supp.out.econ.sex <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEconRisk*sex, data=EgoNet_total_context)
summary(supp.out.econ.sex)
supp.out.health.sex <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentHealthRisk*sex, data=EgoNet_total_context)
summary(supp.out.health.sex)
supp.out.env.sex <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentEnvRisk*sex, data=EgoNet_total_context)
summary(supp.out.env.sex)
supp.out.soc.sex <- glm.nb(support_out~sex+muslim+age_2018+age_2018_sq+New_edu+wealth1000+RecentEconRisk+RecentHealthRisk+RecentEnvRisk+RecentSocRisk+wealth_gini0.1*RecentSocRisk*sex, data=EgoNet_total_context)
summary(supp.out.soc.sex)

(supp.out.econ.sex.plot <-plot_model(supp.out.econ.sex, type="pred", terms= c("RecentEconRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.econ.sex$y)),ceiling(max(supp.out.econ.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.health.sex.plot <-plot_model(supp.out.health.sex, type="pred", terms= c("RecentHealthRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.health.sex$y)),ceiling(max(supp.out.health.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.env.sex.plot <-plot_model(supp.out.env.sex, type="pred", terms= c("RecentEnvRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.env.sex$y)),ceiling(max(supp.out.env.sex$y)),by=2))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))
(supp.out.soc.sex.plot <-plot_model(supp.out.soc.sex, type="pred", terms= c("RecentSocRisk","wealth_gini0.1[2.5,3.5]","sex"),ci.lvl=0.89)+geom_line(aes(color=group))+scale_color_manual(values=c("3.5"="red","2.5"="blue"))+scale_fill_manual(values=c("3.5"="red","2.5"="blue"))+scale_x_continuous(breaks=seq(0,3,by=1))+scale_y_continuous(breaks=seq(floor(min(supp.out.soc.sex$y)),ceiling(max(supp.out.soc.sex$y)),by=1))+theme(axis.text=element_text(size=14),legend.position="none")+ggtitle(NULL)+labs(x=NULL,y=NULL))

models <- list("Economic risk * Gini * Male"=supp.out.econ.sex, "Health risk * Gini * Male"=supp.out.health.sex, "Env. risk * Gini * Male"=supp.out.env.sex, "Social risk * Gini * Male"=supp.out.soc.sex)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>% 
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "3waySuppOut_251014.csv", row.names = FALSE)

##### 3.7. Summarizing main effect models without interactions

models <- list("Network size"=prtn.main, "Reciprocity"=mutual.main, "Long-distance"=long.main, "Local"=local.main, "Support-receiving"=supp.in.main.pois, "Support-providing"=supp.out.main)

model_summaries <- do.call(rbind, lapply(names(models), function(name) {
  summarize_model(models[[name]], name)
}))

model_summaries <- model_summaries %>%
  mutate(
    IRR_CI = paste0(
      sprintf("%.2f", IRR_OR),
      " (", sprintf("%.2f", CI_Lower), ", ", sprintf("%.2f", CI_Upper), ")"
    )
  )

formatted_table <- model_summaries %>%
  dplyr::select(Model, Variable, IRR_CI) %>%
  pivot_wider(
    names_from = Model,
    values_from = IRR_CI
  )

formatted_table[is.na(formatted_table)] <- ""
print(formatted_table)
write.csv(formatted_table, "MainEffectModels_20NN.csv", row.names = FALSE)

# Prepare the data for plotting
plot_data <- model_summaries %>%
  filter(Variable %in% c("RecentEconRisk", "RecentHealthRisk", "RecentEnvRisk", "RecentSocRisk", "wealth_gini0.1")) %>%  # Include only specified variables
  mutate(
    Model = factor(Model, levels = unique(Model)),  # Order models
    Variable = factor(Variable, levels = rev(c("RecentEconRisk", "RecentHealthRisk", "RecentEnvRisk", "RecentSocRisk", "wealth_gini0.1"))),  # Specify variable order
    IRR_OR = as.numeric(IRR_OR),  # Ensure numeric IRR
    CI_Lower = as.numeric(CI_Lower),  # Ensure numeric CI_Lower
    CI_Upper = as.numeric(CI_Upper)  # Ensure numeric CI_Upper
  )

risk_data <- plot_data %>% filter(Variable != "wealth_gini0.1")
wealth_data <- plot_data %>% filter(Variable == "wealth_gini0.1")

# Risk-related variables plot
risk_plot <- ggplot(risk_data, aes(x = IRR_OR, y = Variable)) +
  geom_point(position = position_dodge(width = 0.6), size = 3, color = "black") +  # Black points
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                 position = position_dodge(width = 0.6), height = 0.2, color = "black") +  # Black error bars
  geom_vline(xintercept = 1.0, linetype = "solid", color = "black", size = 0.5) +
  facet_grid(~ Model) +  # Models aligned horizontally with shared x-axis
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Add outer lines for each panel
  )

# wealth_gini0.1 plot
wealth_plot <- ggplot(wealth_data, aes(x = IRR_OR, y = Variable)) +
  geom_point(position = position_dodge(width = 0.6), size = 3, color = "black") +  # Black points
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                 position = position_dodge(width = 0.6), height = 0.2, color = "black") +  # Black error bars
  geom_vline(xintercept = 1.0, linetype = "solid", color = "black", size = 0.5) +
  facet_grid(~ Model) +  # Separate x-axes for each model
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Add outer lines for each panel
  )

combined_plot <- risk_plot / wealth_plot +  # Stacks risk on top of gini
  plot_layout(heights = c(4, 1))  # Adjust height ratio (risk plot is larger)
print(combined_plot)
