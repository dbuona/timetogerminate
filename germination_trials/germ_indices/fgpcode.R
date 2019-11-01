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
library(lme4)
library(tibble)

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

realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY) #elimiate 0 values for log logistic dist

realdatshorty<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis"))

fgp.dat<-filter(realdatshorty,DAY==25)

#full.fgp.mod<-brms::brm(germ_perc~force*chillweeks+(force*chillweeks|Taxa),data=fgp.dat,iter=4000,warmup=3000)

ggplot(fgp.dat,aes(chillweeks,germ_perc))+
  geom_smooth(method="loess",level=0.9,aes(color=as.factor(force),fill = as.factor(force)))+
  facet_wrap(~Taxa)+scale_fill_manual(values=c("royalblue","firebrick1"))+
  scale_color_manual(values=c("royalblue","firebrick1"))+geom_hline(aes(yintercept=0),color="black")+geom_hline(aes(yintercept=1),color="black")+
  ylim(-.3,1.3)+theme_linedraw()+geom_point(aes(color=as.factor(force)),size=0.5)+geom_vline(aes(xintercept=8),color="gray",size=2)


Cc.fgp<-filter(fgp.dat,Taxa=="Cryptotaenia canadensis")
ggplot(Cc.fgp,aes(chillweeks,germ_perc))+geom_point(aes(color=as.factor(force)),size=0.2,shape=1)+geom_smooth(method="loess", aes(color=as.factor(force)),size=1.2)+theme_minimal(base_size = 6)

unique(Cc.fgp$chillweeks)
Cc.0<-filter(Cc.fgp,chillweeks==c(0))
Cc.2<-filter(Cc.fgp,chillweeks==c(0,2))
Cc.4<-filter(Cc.fgp,chillweeks==c(0,2,4))
Cc.5<-filter(Cc.fgp,chillweeks==c(0,2,4,5))
Cc.6<-filter(Cc.fgp,chillweeks==c(0,2,3,4,5,6))
Cc.6

zero<-brms::brm(germ_perc~force*chillweeks,data=Cc.0,iter=4000,warmup=3000)





