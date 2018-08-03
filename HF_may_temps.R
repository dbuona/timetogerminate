setwd("~/Documents/git/timetogerminate")
temp<-read.csv("hf000-01-daily-m.csv",header=TRUE)
library(tidyverse)
temp<-separate(temp,date,c("year","month","day"),remove=FALSE)

temp<-filter(temp,month=="05")
temp<-filter(temp, year>1979)
mean(temp$airt,na.rm = TRUE)
mean(temp$airtmax,na.rm = TRUE)
mean(temp$airtmin,na.rm=TRUE)

temp$fluct<-temp$airtmax-temp$airtmin
mean(temp$fluct,na.rm=TRUE)

quantile(temp$fluct, na.rm=TRUE)

quantile(temp$airtmax,na.rm=TRUE)

19*1.8+32

25*1.8+32
