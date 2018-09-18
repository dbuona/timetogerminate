
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")

data<-read.csv("nochill_practice_data.csv", header=TRUE)
colnames(data)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18"))
ncol(data)
data<-gather(data,"date","germination",7:18)
data$germination<-as.numeric(data$germination)

data$date<-as.Date(data$date,format =  "%m/%d/%y")
class(data$Date)
data$day<-yday(data$date)
unique(data$day)
start<-yday("2018/08/27")
data$DAY<-data$day-start

data$germination<-ifelse(is.na(data$germination),0,data$germination)





ggplot(data, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)


