rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

library(rstan)
library(tidyr)
library(drc)
library(dplyr)
library(shinystan)
library(extraDistr)
library(bayesplot)
library(ggplot2)




realdat<-read.csv("input/daily_dat_nointerval.csv")

##clean data
realdat$germ_perc<-NA
realdat$germ_perc<-realdat$germ_num/realdat$tot_seed
realdat$germ_perc<-ifelse(realdat$germ_perc>1,1,realdat$germ_perc)

##make chilling numeric
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

#write.csv(realdat,"germ_perc_data.csv", row.names = FALSE)


realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY) #elimiate 0 values for log logistic dist

realdatshorty<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))

##break yup datasheet
specieslist<-sort(unique(realdat$Taxa))
X<-split(realdat, with(realdat, realdat$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:20]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)




##Anemone
aa<-filter(`Anemone virginana`,INC=="L")
a<-drm(germ_perc~DAY,factor(chillweeks), data=aa,fct = LL.3(fixed = c(-60, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")

anem<-as.data.frame(coef(a))
anem$chillweeks<-rep(c(0,2,4,5,6,7,8,9,11,13),2)
anem$param<-rep(c("gmax","e50"),each=10)
colnames(anem)[1]<-"estimate"
anem.gmax<-filter(anem,param=="gmax")
anem.e50<-filter(anem,param!="gmax")
ggplot(anem.gmax,aes(chillweeks,estimate))+geom_point()
ggplot(anem.e50,aes(chillweeks,estimate))+geom_point()


b<-drm(germ_perc~DAY,factor(chillweeks):factor(INC), data=`Asclepias syriaca`,fct = LL.3(fixed = c(-50, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
plot(b)
ascp<-as.data.frame(coef(b))
ascp$chillweeks<-rep(c(0,2,4,5,6,7,8,9,11,13),each=2,2)
ascp$param<-rep(c("gmax","e50"),each=20)
ascp$INC<-rep(c("H","L"),20)

colnames(ascp)[1]<-"estimate"
ascp.gmax<-filter(ascp,param=="gmax")
ascp.e50<-filter(ascp,param!="gmax")
plot.asc.gmax<-ggplot(ascp.gmax,aes(chillweeks,estimate))+geom_point(aes(color=INC))+ylim(0,1)
plot.asc.e50<-ggplot(ascp.e50,aes(chillweeks,estimate))+geom_point(aes(color=INC))


f<-drm(germ_perc~DAY,factor(chillweeks):factor(INC), data=`Eurbia diviricata`,fct = LL.3(fixed = c(-60, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
plot(f)

eurb<-as.data.frame(coef(f))
eurb$chillweeks<-rep(c(0,2,4,5,6,7,8,9,11,13),each=2,2)
eurb$param<-rep(c("gmax","e50"),each=20)
eurb$INC<-rep(c("H","L"),20)

colnames(eurb)[1]<-"estimate"
eurb.gmax<-filter(eurb,param=="gmax")
eurb.e50<-filter(eurb,param!="gmax")


plot.eurb.gmax<-ggplot(eurb.gmax,aes(chillweeks,estimate))+geom_point(aes(color=INC))+ylim(0,1)
plot.eurb.e50<-ggplot(eurb.e50,aes(chillweeks,estimate))+geom_point(aes(color=INC))

ggpubr::ggarrange(plot.eurb.gmax,plot.asc.gmax)



c<-drm(germ_perc~DAY,factor(chillweeks):factor(INC), data=`Silene vulgaris`,fct = LL.3(fixed = c(-50, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")


sil<-as.data.frame(coef(c))
sil$chillweeks<-rep(c(0,2,4,5,6,7,8,9,11,13),each=2,2)
sil$param<-rep(c("gmax","e50"),each=20)
sil$INC<-rep(c("H","L"),20)

colnames(sil)[1]<-"estimate"
sil.gmax<-filter(sil,param=="gmax")
sil.e50<-filter(sil,param!="gmax")




coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(ascp.gmax,ascp.gmax$INC=="H"), start = list(g0 = 0, gmax =1,gamma_g=.01)))
coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(ascp.gmax,ascp.gmax$INC=="L"), start = list(g0 = 0, gmax =1,gamma_g=.01)))
coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(eurb.gmax,eurb.gmax$INC=="H"), start = list(g0 = 0, gmax =1,gamma_g=.025)))
coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(eurb.gmax,eurb.gmax$INC=="L"), start = list(g0 = 0, gmax =1,gamma_g=.025)))
coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(sil.gmax,sil.gmax$INC=="H"), start = list(g0 = 0, gmax =1,gamma_g=.025)))
coef(nls(estimate~g0+(gmax-g0)*(1-exp(1)^(-gamma_g*chillweeks)),data=subset(sil.gmax,sil.gmax$INC=="L"), start = list(g0 = 0, gmax =1,gamma_g=.025)))

####notes I don't feel like I exlored the space well
#but gamma_g ran from 0.2 to about 1.22 comepared to chill weeks.

