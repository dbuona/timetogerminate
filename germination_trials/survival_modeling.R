####survival analysis
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library(lubridate)
library("Hmisc")
library(drc)
library(survival)

d<-read.csv("germ_data_forDRC.csv",header= TRUE)
###needa binary outcum

####use goop from drc_modeling
goop<- within(goop, DAY[DAY==-Inf]<-NA)
goop<- within(goop, END[END==Inf]<-NA)


Surv.Obj <- Surv(goop$DAY, goop$END,goop$germination)

goon<-survfit(formula = Surv.Obj ~ goop$plate_num, conf.type = "log-log")
?survfit.coxph.fit()
summary(goon)
goon
?Surv()
