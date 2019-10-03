rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/Input")
library(ggplot2)
library(tidyr)
library(dplyr)
library(chillR)
library(brms)
library(tibble)

days<-42 #number of days of climate data
yr<-2018#year of climate data
JDay<-seq(1:days)#DOYs- in January in this case
Year<-rep(yr,length(JDay))
lat<-45.5#latitude

#2. Estimate chilling using Tmin and Tmax
#use min and max daily temp to create fake climate data frame for January 2018 
Tmin<- 4
Tmax<-4
minmaxdaily<-data.frame(Year,JDay,Tmin,Tmax)

#Convert daily tmin and tmax data to hourly data, using packaged functions in chillR
hrly<-stack_hourly_temps(latitude=lat,make_hourly_temps(lat, minmaxdaily, keep_sunrise_sunset = FALSE))$hourtemps

#Chilling calculations for lat and time period of interest (January 2018 in my example)
chillcalc <- chilling(hrly, hrly$JDay[1], hrly$JDay[nrow(hrly)]) 
#1008 Utah 30.35 Cp
###now when historically did we reach that?
weather<-read.csv("..//FLOBUDS/data/hf000-01-daily-m.csv",header = TRUE)
weather<-select(weather,c("date","airtmax","airtmin"))
weather<-separate(weather,date,c("Year","Month","Day"),sep="-",remove=TRUE)
colnames(weather)<-c("Year","Month","Day","Tmax","Tmin")
sapply(weather,mode) #mode(weather)
weather$Year<-as.numeric(weather$Year)
weather$Month<-as.numeric(weather$Month)
weather$Day<-as.numeric(weather$Day)

weather<-filter(weather,Year<1980)

weatherCC<-weather


weather<-make_all_day_table(weather)

hourtemps<-stack_hourly_temps(weather, latitude=42.5)$hourtemps ## make hourly
hourtemps$DATE<-ISOdate(hourtemps$Year,hourtemps$Month,hourtemps$Day,hourtemps$Hour)


##how much cold and GDD did they get before bud burst int aht year
ChillHF<-as.data.frame(chilling(hourtemps,Start_JDay=275,95))
mean(ChillHF$Utah_Model)

weatherCC<-weather
weatherCC$Tmin<-weather$Tmin+4
hourtempsCC<-stack_hourly_temps(weatherCC, latitude=42.5)$hourtemps ## make hourly
hourtempsCC$DATE<-ISOdate(hourtempsCC$Year,hourtempsCC$Month,hourtempsCC$Day,hourtempsCC$Hour)
ChillHFCC<-as.data.frame(chilling(hourtempsCC,Start_JDay=275,))
mean(ChillHFCC$Utah_Model)
