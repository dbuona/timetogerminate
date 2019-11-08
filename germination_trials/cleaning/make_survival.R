###a script to amke my data survival formated

rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library("ggfortify")
library("icenReg")

library("emmeans")

###make a small data set
#day<-c(0,2,4,6,8,10,0,2,4,6,8,10)
#end<-c(2,4,6,8,10,Inf,2,4,6,8,10,Inf)
#germination<-c(0,0,4,1,2,3,0,2,3,3,1,1)
#name<-c("A","A","A","A","A","A","B","B","B","B","B","B"


#dater<-uncount(dater,germination,.remove=TRUE)
#dater$germinationed<-ifelse(dater$day==10,0,1)
###This seems to work
setwd("~/Documents/git/timetogerminate/germination_trials/input")
daily.dat<-read.csv("daily_dat_nointerval.csv",header=TRUE)

###Note: not worrying about interval censoring right now so skip
#d<-read.csv("germ_data_forDRC.csv",header= TRUE)
#goober<-uncount(d,germination,.remove=FALSE)


#goober$END<-ifelse(goober$END==Inf,NA,goober$END)
#goober$censor<-ifelse(goober$END==0,2,1)
#goober<- within(goober, censor[DAY==25 ]<-0)
#goober$chill_time<-NA
#goober<- within(goober, chill_time[COLD=="O" ]<-0)
#goober<- within(goober, chill_time[COLD=="A" ]<-14)
#goober<- within(goober, chill_time[COLD=="B" ]<-28)
#goober<- within(goober, chill_time[COLD=="C" ]<-35)
#goober<- within(goober, chill_time[COLD=="D" ]<-42)
#goober<- within(goober, chill_time[COLD=="E" ]<-49)
#goober<- within(goober, chill_time[COLD=="f" ]<-56)
#goober<- within(goober, chill_time[COLD=="G" ]<-63)
#goober<- within(goober, chill_time[COLD=="H" ]<-77)
#goober<- within(goober, chill_time[COLD=="i" ]<-91)#

#goober$warmT<-NA
#goober<- within(goober, warmT[INC=="H" ]<-25)
#goober<- within(goober, warmT[INC=="L" ]<-20)
#write.csv(goober,"interval_survival_data.csv")

####converts toruvival anaysis
goo<-uncount(daily.dat,germ.daily,.remove=FALSE)
goo$germinated<-1

goo2<-filter(daily.dat,DAY==25 &germ.daily==0)
goo2<-uncount(goo2,viable,.remove=FALSE)
goo2$germinated<-0


data<-rbind(goo,goo2)
colnames(data)
surv.data<-dplyr::select(data,Taxa,INC,COLD,plate_num,date,DAY,germinated)
goober<-surv.data
goober$chill_time<-NA
goober<- within(goober,chill_time[COLD=="0" ]<-0)
goober<- within(goober, chill_time[COLD=="A" ]<-14)
goober<- within(goober, chill_time[COLD=="B" ]<-28)
goober<- within(goober, chill_time[COLD=="C" ]<-35)
goober<- within(goober, chill_time[COLD=="D" ]<-42)
goober<- within(goober, chill_time[COLD=="E" ]<-49)
goober<- within(goober, chill_time[COLD=="f" ]<-56)
goober<- within(goober, chill_time[COLD=="G" ]<-63)
goober<- within(goober, chill_time[COLD=="H" ]<-77)
goober<- within(goober, chill_time[COLD=="i" ]<-91)

goober$warmT<-NA
goober<- within(goober, warmT[INC=="H" ]<-25)
goober<- within(goober, warmT[INC=="L" ]<-20)

write.csv(goober,"surival_dat_nointerval.csv") ###not this data sheet throws out any seeds that were moldy
