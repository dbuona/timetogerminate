##### calculate the effect of soil warming on "stratification" condtions

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(ggplot2)
library(tidyr)

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

soiltemp<-read.csv("input/hf018-03-soil-temp.csv"  )
soiltemp.wint<-filter(soiltemp,doy>288) #oct15
soiltemp.spr<-filter(soiltemp,doy<106) #april 121
soiltemp<-rbind(soiltemp.spr,soiltemp.wint)


###first look at daily means
soiltempmeans<- soiltemp %>%group_by(year,doy) %>% summarise(meanctr=mean(ctrl.av),meanp1=mean(p1.av),
meanp2=mean(p2.av),meanp3=mean(p3.av))
soiltempmeans$ctrlcount<-ifelse(soiltempmeans$meanctr>=1 & soiltempmeans$meanctr<=6,1,0)
soiltempmeans$p1count<-ifelse(soiltempmeans$meanp1>=1 & soiltempmeans$meanp1<=6,1,0)
soiltempmeans$p2count<-ifelse(soiltempmeans$meanp2>=1 & soiltempmeans$meanp2<=6,1,0)
soiltempmeans$p3count<-ifelse(soiltempmeans$meanp3>=1 & soiltempmeans$meanp3<=6,1,0)

counts<- soiltempmeans %>% group_by(year) %>% summarise(ctrl=sum(ctrlcount,na.rm=TRUE),p3=sum(p3count,na.rm=TRUE))
counts$diff<-counts$ctrl-counts$p3
mean(counts$ctrl)
sd(counts$ctrl)

new.data.2<-data.frame(chill_time=c(rnorm(1000,79.2,26.63),rnorm(1000,47.73,25.7)),force=rep(c(0,1),each=1000))

mean(counts$p3)
sd(counts$p3)

##now min max
soiltempminmax<- soiltemp %>%group_by(year,doy) %>% summarise(minctr=min(ctrl.av),maxctr=max(ctrl.av),minp3=min(p3.av),maxp3=max(p3.av))
soiltempminmax$ctrlcount<-ifelse(soiltempminmax$minctr>=1 & soiltempminmax$maxctr<=6,1,0)
soiltempminmax$p3count<-ifelse(soiltempminmax$minp3>=1 & soiltempminmax$maxp3<=6,1,0)
counts2<- soiltempminmax %>% group_by(year) %>% summarise(ctrl=sum(ctrlcount,na.rm=TRUE),p3=sum(p3count,na.rm=TRUE))
counts2$diff<-counts2$ctrl-counts2$p3
mean(counts2$diff)

max(counts$ctrl)
min(counts$ctrl)
105-27
max(counts$p3)
min(counts$p3)

meanplot<-gather(counts,plot,estimate,2:3)
ggplot(meanplot,aes(plot,estimate))+stat_summary()
