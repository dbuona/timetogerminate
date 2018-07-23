setwd("~/Documents/git/timetogerminate")
temp<-read.csv("hf000-01-daily-m.csv",header=TRUE)
library(tidyverse)
temp<-separate(temp,date,c("year","month","day"),remove=FALSE)

temp<-filter(temp,month=="05")
temp<-filter(temp, year>1979)
a<-mean(temp$airt,na.rm = TRUE)
b<-mean(temp$airtmax,na.rm = TRUE)
c<-mean(temp$airtmin,na.rm=TRUE)

d<-1.8
a*d+32
b*d+32
c*d+32
22*d+32
a
b
c
max(temp$year)
