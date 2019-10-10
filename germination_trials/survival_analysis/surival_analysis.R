###a script to amke my data survival formated

rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(survival)
library(ggplot2)
library(dplyr)
library("ggfortify")
library("coxme")
library("icenReg")
library("emmeans")
library(drc)
library(lme4)
library(brms)
setwd("~/Documents/git/timetogerminate/germination_trials/input")

d<-read.csv("..//survival_analysis/surival_dat_nointerval.csv")

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

#d$DAY<-ifelse(d$DAY==0,0.001,d$DAY)

as<-filter(realdat,Taxa=="Asclepias syriaca")



priorz<-get_prior(DAY | cens(germinated)~force+chillweeks,data=as,family = lognormal(link = "identity", link_sigma = "log"))

fit2 <- brm(DAY | cens(germinated)~force+chillweeks, 
            data=as, family = lognormal(link = "identity", link_sigma = "log") ,prior=priorz,iter=6000,warmup = 5000) ##
?brmsfamily()
launch_shinystan(fit2)
summary(fit2)
###nonlinear mixed effect model: 
#1ignore those that didn't germinate


ggplot(germy,aes(DAY))+geom_histogram(bins=25)+facet_wrap(~Taxa)



