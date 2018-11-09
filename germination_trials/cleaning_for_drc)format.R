###Reformatdata

rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")
data6<-read.csv("time_course_practice6.csv",header=TRUE)
data6<-dplyr::select(data6,1:20)
colnames(data6)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/15/18","10/17/18","10/19/18","10/21/18","10/22/18","10/24/18","10/26/18","10/28/18","10/29/18","10/31/18","11/2/18","11/4/18","11/5/18","11/7/18"))
data6<-data6 %>% replace(is.na(.), 0)


data6<-gather(data6,"date","germination",7:20)
data6$germination<-as.numeric(data6$germination)

data6$date<-as.Date(data6$date,format =  "%m/%d/%y")
class(data6$date)
data6$day<-yday(data6$date)
unique(data6$day)
start6<-yday("2018/10/15")
data6$start<-data6$day-start6

daty<- data6 %>% group_by(plate_num) %>% mutate(germ.daily = germination - lag(germination))
daty<- daty%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))

daty<- daty %>% group_by(plate_num) %>% mutate(end = lead(start))
daty$germ.daily[is.na(daty$germ.daily)]<-0 
colnames(daty)
dater<-dplyr::select(daty,tot_seed,Taxa,INC, COLD, plate_num,date,start,end,tru.daily)





                      