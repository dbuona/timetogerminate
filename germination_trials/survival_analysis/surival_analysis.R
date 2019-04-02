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
#d$DAY<-ifelse(d$DAY==0,0.001,d$DAY)

as<-filter(d,Taxa=="Asclepias syriaca")

as.warm.rest<-filter(as, COLD %in% c("D"))
priorz<-get_prior(DAY | cens(germinated)~INC,data=as.warm.rest,family = weibull())

fit2 <- brm(DAY | cens(germinated)~INC, 
            data=as.warm.rest, family = weibull(), ,prior=priorz) ##
launch_shinystan(fit2)
summary(fit2)
###nonlinear mixed effect model: 
#1ignore those that didn't germinate


ggplot(germy,aes(DAY))+geom_histogram(bins=25)+facet_wrap(~Taxa)



