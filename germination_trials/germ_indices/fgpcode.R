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

full.fgp.mod<-brms::brm(germ_perc~force*chillweeks+(force*chillweeks|Taxa),data=fgp.dat,iter=4000,warmup=3000)

ranef(full.fgp.mod)

ggplot(fgp.dat,aes(chillweeks,germ_perc))+
  geom_smooth(method="loess",level=0.9,aes(color=as.factor(force),fill = as.factor(force)))+
  facet_wrap(~Taxa)+scale_fill_manual(values=c("royalblue","firebrick1"))+
  scale_color_manual(values=c("royalblue","firebrick1"))+geom_hline(aes(yintercept=0),color="black")+geom_hline(aes(yintercept=1),color="black")+
  ylim(-.3,1.3)+theme_linedraw()+geom_point(aes(color=as.factor(force)),size=0.5)+geom_vline(aes(xintercept=8),color="gray",size=2)


unique(fgp.dat$Taxa)
Cc<-filter(fgp.dat,Taxa=="Cryptotaenia canadensis")

CC.C<-filter(Cc, force==0)
CC.C.loes<-ggplot(CC.C,aes(chillweeks,germ_perc))+geom_point(color="royalblue")+geom_smooth(method="loess",level=0.9,fill="royalblue")+theme_linedraw()
CC.C.lm<-ggplot(CC.C,aes(chillweeks,germ_perc))+geom_point(color="royalblue")+geom_smooth(method="lm",level=0.9,fill="royalblue")+theme_linedraw()

CC.C.mod.glm<-glm(germ_perc~chillweeks,CC.C,family=binomial(link = "logit"))


CC.C.break.glm<-segmented(CC.C.mod.glm,seg.Z=~chillweeks,npsi=2)
seg.glm<-plot.segmented(CC.C.break.glm,main="GLM")

CC.C.mod.lm<-lm(germ_perc~chillweeks,CC.C)



CC.C.break.lm<-segmented(CC.C.mod.lm,seg.Z=~chillweeks,npsi=2)
seg.lm<-plot.segmented(CC.C.break.lm,main="LM")

par(mfrow=c(1,2))


