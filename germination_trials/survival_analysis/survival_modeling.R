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
library(smcure)
d<-read.csv("survival_analysis/surival_dat_nointerval.csv",header= TRUE)
###needa binary outcum
d$DAY<-ifelse(d$DAY==0,0.001,d$DAY)

realdat<-d
realdat$chillweeks<-realdat$chill_time/7 # make chilling weeks instead of days

realdat$force<-NA # make forcing numeric
realdat<- within(realdat, force[INC=="L"]<-0)#
realdat<- within(realdat, force[INC=="H"]<-5)

CC<-filter(realdat,Taxa=="Cryptotaenia canadensis")



Surv.Obj <- Surv(CC$DAY, CC$germinated,type = "right")


goo<-smcure(Surv.Obj~chillweeks+force,~chillweeks+force,offset=NULL,CC,na.action=na.omit,model=c("aft"),link="probit",Var=TRUE,emmax=50,eps=1e-7,nboot=100)
printsmcure(goo)
?smcure()


 predgoo<-predictsmcure(goo,newX=cbind(c(0,4,6,8,12),c(0,0,0,0,0)),
                      newZ=cbind(c(0,4,6,8,12),c(0,0,0,0,0)),model="aft")


 plotpredictsmcure(predgoo,model="aft")

?plotpredictsmcure

survregLogNormal <-survreg(Surv.Obj ~ chillweeks+force, CC, dist = "lognormal")
summary(survregLogNormal)

fit1 <- survfit(Surv.Obj ~ chillweeks+force, data = CC)
summary(fit1)
ggsurvplot(fit1, data = CC, pval = TRUE)
