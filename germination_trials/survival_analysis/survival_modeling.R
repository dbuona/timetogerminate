####survival analysis
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials")
library(tidyverse)
library(lubridate)
library("Hmisc")
library(drc)
library(survival)
library("survminer")

d<-read.csv("survival_analysis/surival_dat_nointerval.csv",header= TRUE)
###needa binary outcum
d$DAY<-ifelse(d$DAY==0,0.001,d$DAY)

realdat<-d
realdat$chill_time<-NA
realdat<- within(realdat, chill_time[COLD=="0" ]<-0)
realdat<- within(realdat, chill_time[COLD=="A" ]<-14)
realdat<- within(realdat, chill_time[COLD=="B" ]<-28)
realdat<- within(realdat, chill_time[COLD=="C" ]<-35)
realdat<- within(realdat, chill_time[COLD=="D" ]<-42)
realdat<- within(realdat, chill_time[COLD=="E" ]<-49)
realdat<- within(realdat, chill_time[COLD=="f" ]<-56)
realdat<- within(realdat, chill_time[COLD=="G" ]<-63)
realdat<- within(realdat, chill_time[COLD=="H" ]<-77)
realdat<- within(realdat, chill_time[COLD=="i" ]<-91)
realdat$chillweeks<-realdat$chill_time/7 # make chilling weeks instead of days

realdat$force<-NA # make forcing numeric
realdat<- within(realdat, force[INC=="L"]<-0)
realdat<- within(realdat, force[INC=="H"]<-1)

CC<-filter(realdat,Taxa=="Cryptotaenia canadensis")



Surv.Obj <- Surv(CC$DAY, CC$germinated,type = "right")

survregLogNormal <-survreg(Surv.Obj ~ chillweeks+force, CC, dist = "lognormal")
summary(survregLogNormal)

fit1 <- survfit(Surv.Obj ~ chillweeks+force, data = CC)
summary(fit1)
ggsurvplot(fit1, data = CC, pval = TRUE)
