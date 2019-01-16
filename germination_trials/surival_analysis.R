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
setwd("~/Documents/git/timetogerminate/germination_trials/input")

d<-read.csv("surival_dat_nointerval.csv")


###nonlinear mixed effect model: 
#1ignore those that didn't germinate
?survreg.distributions
germy<-filter(d,germinated==1)

mod <- survreg(Surv(DAY,germinated) ~ INC+COLD+frailty(Taxa), dist = "gaussian", data = d) 
summary(mod)

emmeans(mod,~COLD+INC,transform="response")

mod2 <- survreg(Surv(DAY,germinated) ~ INC+COLD, dist = "gaussian", data = germy) 
summary(mod2)
emmeans(mod2,~COLD+INC,transform="response")

ggplot(germy,aes(DAY))+geom_histogram(bins=25)+facet_wrap(~Taxa)



