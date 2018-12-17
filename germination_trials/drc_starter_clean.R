###This formats cleaned data for DCR
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library(lubridate)
library("Hmisc")

source("germination_cleaning.R")

plateid<-unique(cum_data_full$plate_num)


daily.dat<- cum_data_full %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
daily.dat<- daily.dat%>% group_by(plate_num) %>% mutate(germination = lead(germ.daily))
daily.dat$germination<-ifelse(daily.dat$DAY==-Inf,daily.dat$germ_num,daily.dat$germination)
daily.dat<- daily.dat %>% group_by(plate_num) %>% mutate(END = lead(DAY))
daily.dat<-dplyr::select(daily.dat,Taxa,INC, COLD, plate_num,date,viable,DAY,END,germination,mold,adj_total)
daily.dat$END<-ifelse(is.na(daily.dat$END),Inf,daily.dat$END)
daily.dat$germination<-ifelse(is.na(daily.dat$germination),daily.dat$viable,daily.dat$germination)




unique(daily.dat$plate_num)
daily.dat$force<-ifelse(daily.dat$INC=="H",25,20)
daily.dat$germ_perc<-daily.dat$germination/daily.dat$adj_total

write.csv(daily.dat,"germ_data_forDRC.csv")



##this combines cold. but I cant get them to run together


specieslist<-unique(daily.dat$Taxa)
listE50<-list()
for (sp in seq_along(specieslist)){
  dataonesp <- subset(daily.dat, Taxa==specieslist[sp])
  mod <- drm(germination~start+end, factor(INC):factor(COLD), data = dataonesp, fct = LL.3(), type = "event")
  listE50[[paste(sp, specieslist[sp])]]<-list(ED(mod,c(50)))
}

###this combines species
listgoob<-list()
for (tr in seq_along(treatlist)){
  dataonesp <- subset(master.dat, COLD==treatlist[tr])
  mod <- drm(tru.daily~start+end, factor(INC), data = dataonesp, fct = LL.3(), type = "event")
  listgoob[[paste(tr, treatlist[tr])]]<-list(ED(mod,c(50,75)))
}

listgoob
listE50

