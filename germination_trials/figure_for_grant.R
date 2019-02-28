### for figure
###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(tidyverse)
library(ggthemes)
library(drc)
setwd("~/Documents/git/timetogerminate/germination_trials")
real.data<-read.csv("input/daily_dat_nointerval.csv")
real.data<- real.data %>% group_by(plate_num)%>% mutate(germ.cum=cumsum(germ.daily))

unique(real.data$Taxa)
comp<-real.data %>% filter(Taxa %in% c("Cryptotaenia canadensis","Eurbia diviricata"))

comp<-comp %>% filter(COLD %in% c("B","f"))
#comp<-filter(comp,INC=="H")

comp$climate<-NA
comp$climate[comp$COLD=="f" & comp$INC=="H"]<-"warm"
comp$climate[comp$COLD=="B" & comp$INC=="L"]<-"cool"
comp<-filter(comp,!is.na(climate))


ggplot(comp, aes(x = DAY, y = germ.cum, color=Taxa)) + stat_summary(alpha=0.7)+theme_bw()+geom_line(stat = "summary", fun.y = mean)+facet_wrap(~climate)

comp1<-filter(comp, climate=="cool")
comp1.asc<-filter(comp1,Taxa=="Cryptotaenia canadensis")
com1.sv<-filter(comp1,Taxa=="Eurbia diviricata")
mod.1<-drm(germ.cum~DAY,,fct=LL.3(c(NA,NA,NA)),data=comp1.asc,type="continuous")
mod.2<-drm(germ.cum~DAY,,fct=LL.3(c(NA,NA,NA)),data=com1.sv,type="continuous")
summary(mod.1)


comp2<-filter(comp, climate=="warm")
comp2.asc<-filter(comp2,Taxa=="Cryptotaenia canadensis")
com2.sv<-filter(comp2,Taxa=="Eurbia diviricata")

mod.3<-drm(germ.cum~DAY,,fct=LL.3(c(NA,NA,NA)),data=comp2.asc,type="continuous")
mod.4<-drm(germ.cum~DAY,,fct=LL.3(c(NA,NA,NA)),data=com2.sv,type="continuous")

par(mfrow=c(1,2))
par(mar=c(4,4,1.8,0.5))
plot(mod.1,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="darkgreen",pch=20,cex=0.5,xlab="time (days)", ylab="cumulative germination",main="cool conditions")
plot(mod.1, broken = TRUE, type="confidence", add=TRUE, col="darkgreen")
plot(mod.2,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="navyblue",pch=20,add=TRUE,cex=0.5)
plot(mod.2, broken = FALSE, type="confidence", add=TRUE,col="navyblue")


plot(mod.3,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="darkgreen",pch=20,cex=0.5,ylab="",main="warm conditions",xlab="time (days)")
plot(mod.3, broken = FALSE, type="confidence", add=TRUE, col="darkgreen")
plot(mod.4,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="navyblue",pch=20,add=TRUE,cex=0.5)
plot(mod.4, broken = FALSE, type="confidence", add=TRUE,col="navyblue")


C.canadensis<-as.factor(c(1,0,4,3,2,1,0,2,1,0,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14))
E.divaricata<-as.factor(c(0,1,0,1,2,3,4,0,1,2,8,7,6,5,4,3,2,1,0,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0))

ggplot(,aes(C.canadensis,E.divaricata))+geom_point(shape=16,size=3)+theme_bw()+ylab("# seeds/pot of C. canadensis")+xlab("# seeds/pot of E. divaricata")

