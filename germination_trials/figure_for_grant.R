### for figure
###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(tidyverse)
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
par(mar=c(2,3,1.8,0.5))
plot(mod.1,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="darkgreen",pch=20,cex=0.5, ylab="cumulative germination",main="cool conditions")
plot(mod.1, broken = TRUE, type="confidence", add=TRUE, col="darkgreen")
plot(mod.2,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="navyblue",pch=20,add=TRUE,cex=0.5)
plot(mod.2, broken = FALSE, type="confidence", add=TRUE,col="navyblue")

plot(mod.3,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="darkgreen",pch=20,cex=0.5,ylab="",main="warm conditions")
plot(mod.3, broken = FALSE, type="confidence", add=TRUE, col="darkgreen")
plot(mod.4,ylim=c(0,20),xlim=c(0,25),log="",type="all",legendPos = c(15,5),col="navyblue",pch=20,add=TRUE,cex=0.5)
plot(mod.4, broken = FALSE, type="confidence", add=TRUE,col="navyblue")


