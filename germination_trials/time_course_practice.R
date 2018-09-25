
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")

data<-read.csv("time_course_practice.csv", header=TRUE)
data2<-read.csv("time_course_practice2.csv", header=TRUE)

colnames(data)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18"))
ncol(data)
data2$"9/10/18"<-0
colnames(data2)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num", "9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/10/18"))

data<-gather(data,"date","germination",7:20)
data$germination<-as.numeric(data$germination)

data2<-gather(data2,"date","germination",7:13)
data2$germination<-as.numeric(data2$germination)

data$date<-as.Date(data$date,format =  "%m/%d/%y")
class(data$Date)
data$day<-yday(data$date)
unique(data$day)
start<-yday("2018/08/27")
data$DAY<-data$day-start

data2$date<-as.Date(data2$date,format =  "%m/%d/%y")
class(data2$date)
data2$day<-yday(data2$date)
unique(data2$day)
start2<-yday("2018/09/10")
data2$DAY<-data2$day-start2

data$germination<-ifelse(is.na(data$germination),0,data$germination)
data2$germination<-ifelse(is.na(data2$germination),0,data2$germination)




ggplot(data, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)
ggplot(data2, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)

full<-rbind(data,data2)

ggplot(full, aes(x = DAY, y = germination, color=INC, shape=COLD)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)


library(germinationmetrics)
d<-filter(data,Taxa=="Asclepias syriaca")
t50(germ.counts = d$DAY, intervals = d$germination, partial = FALSE, method = "farooq")
FourPHFfit(germ.counts = d$DAY, intervals = d$germination, total.seeds = 20, tmax = 20)
