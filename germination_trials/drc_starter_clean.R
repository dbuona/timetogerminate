##This is the main cleaning file for germination work. Begin on november 12, 2018 by Dan B
###it will combine all cold trials into one data sheet format all data for dose response curve in the drc package.

rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")

chill.no<-read.csv("dat_CHILL_NO.CSV",header=TRUE)
colnames(chill.no)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/27/18","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/21/18","viable","chill_germ","mold","adj_total"))

chill.A<-read.csv("data_CHILL_A.csv",header=TRUE)
colnames(chill.A)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","9/10/18","9/12/18","9/14/18",	"9/16/18","9/17/18","9/18/18","9/20/18","9/21/18","9/23/18","9/26/18","9/28/18","9/30/18","10/3/18","10/5/18","viable","chill_germ","mold","adj_total"))

chill.B<-read.csv("dat_chill_B.csv",header=TRUE)
colnames(chill.B)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","9/24/18","9/26/18","9/28/18","9/30/18","10/3/18","10/5/18","10/7/18"  , "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18", "viable"     ,"chill_germ", "mold", "adj_total"))

chill.C<-read.csv("dat_chill_C.csv",header=TRUE)
colnames(chill.C)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/1/18","10/3/18","10/5/18","10/7/18"  , "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","viable"     ,"chill_germ", "mold", "adj_total"))

chill.D<-read.csv("dat_chill_D.csv",header=TRUE)
colnames(chill.D)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num", "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","10/28/18" ,"10/29/18",  "10/31/18" , "11/2/18","viable","chill_germ", "mold", "adj_total"))

chill.E<-read.csv("dat_chill_E.csv",header=TRUE) ####update when finish viability checking
colnames(chill.E)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","10/28/18" ,"10/29/18",  "10/31/18" , "11/2/18","11/4/18",  "11/5/18",   "11/7/18" ,  "11/9/18","viable","chill_germ", "mold", "adj_total"))

##make all nas 0
chill.no<-chill.no %>% replace(is.na(.), 0)
chill.A<-chill.A %>% replace(is.na(.), 0)
chill.B<-chill.B %>% replace(is.na(.), 0)
chill.C<-chill.C %>% replace(is.na(.), 0)
chill.D<-chill.D %>% replace(is.na(.), 0)
chill.E<-chill.E %>% replace(is.na(.), 0)

##3adjust germination percentages

####1] total germinated/ totalviable+total germ
###could also think of mold differently
chill.no$adj_total<-chill.no$"9/21/18"+chill.no$viable
chill.A$adj_total<-chill.A$"10/5/18"+chill.A$viable
chill.B$adj_total<-chill.B$"10/19/18"+chill.B$viable
chill.C$adj_total<-chill.C$"10/26/18"+chill.C$viable
chill.D$adj_total<-chill.D$"11/2/18"+chill.D$viable
chill.E$adj_total<-chill.E$"11/9/18"+chill.E$viable ###This might have to be updated after viability checking

##reformat data
chill.no<-gather(chill.no,"date","germ_num",7:22)
chill.no$germ_num<-as.numeric(chill.no$germ_num)
chill.A<-gather(chill.A,"date","germ_num",7:20)
chill.A$germ_num<-as.numeric(chill.A$germ_num)
chill.B<-gather(chill.B,"date","germ_num",7:20)
chill.B$germ_num<-as.numeric(chill.B$germ_num)
chill.C<-gather(chill.C,"date","germ_num",7:21)
chill.C$germ_num<-as.numeric(chill.C$germ_num)
chill.D<-gather(chill.D,"date","germ_num",7:21)
chill.D$germ_num<-as.numeric(chill.D$germ_num)
chill.E<-gather(chill.E,"date","germ_num",7:21)
chill.E$germ_num<-as.numeric(chill.E$germ_num)


chill.no$date<-as.Date(chill.no$date,format =  "%m/%d/%y")
chill.no$day<-yday(chill.no$date)
start<-yday("2018/08/27")
chill.no$start<-chill.no$day-start

chill.A$date<-as.Date(chill.A$date,format =  "%m/%d/%y")
chill.A$day<-yday(chill.A$date)
startA<-yday("2018/09/10")
chill.A$start<-chill.A$day-startA

chill.B$date<-as.Date(chill.B$date,format =  "%m/%d/%y")
chill.B$day<-yday(chill.B$date)
startB<-yday("2018/09/24")
chill.B$start<-chill.B$day-startB

chill.C$date<-as.Date(chill.C$date,format =  "%m/%d/%y")
chill.C$day<-yday(chill.C$date)
startC<-yday("2018/10/01")
chill.C$start<-chill.C$day-startC

chill.D$date<-as.Date(chill.D$date,format =  "%m/%d/%y")
chill.D$day<-yday(chill.D$date)
startD<-yday("2018/10/08")
chill.D$start<-chill.D$day-startD

chill.E$date<-as.Date(chill.E$date,format =  "%m/%d/%y")
chill.E$day<-yday(chill.E$date)
startE<-yday("2018/10/15")
chill.E$start<-chill.E$day-startE


####format for drc function
daty<- chill.no %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
daty<- daty%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
daty<- daty %>% group_by(plate_num) %>% mutate(end = lead(start))
daty$germ.daily[is.na(daty$germ.daily)]<-0 

datyA<- chill.A %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
datyA<- datyA%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
datyA<- datyA %>% group_by(plate_num) %>% mutate(end = lead(start))
datyA$germ.daily[is.na(datyA$germ.daily)]<-0 

datyB<- chill.B %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
datyB<- datyB%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
datyB<- datyB %>% group_by(plate_num) %>% mutate(end = lead(start))
datyB$germ.daily[is.na(datyB$germ.daily)]<-0 

datyC<- chill.C %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
datyC<- datyC%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
datyC<- datyC %>% group_by(plate_num) %>% mutate(end = lead(start))
datyC$germ.daily[is.na(datyC$germ.daily)]<-0 

datyD<- chill.D %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
datyD<- datyD%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
datyD<- datyD %>% group_by(plate_num) %>% mutate(end = lead(start))
datyD$germ.daily[is.na(datyD$germ.daily)]<-0 

datyE<- chill.E %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
datyE<- datyE%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))
datyE<- datyE %>% group_by(plate_num) %>% mutate(end = lead(start))
datyE$germ.daily[is.na(datyE$germ.daily)]<-0 

#dater.no<-dplyr::select(daty,tot_seed,Taxa,INC, COLD, plate_num,date,start,end,tru.daily)

daty$COLD<-"O"

master.dat<-rbind(daty,datyA,datyB,datyC,datyD,datyE)

###finish formatting for drc 
master.dat$end<-ifelse(is.na(master.dat$end),Inf,master.dat$end)
master.dat$tru.daily<-ifelse(is.na(master.dat$tru.daily),master.dat$viable,master.dat$tru.daily)

###select colums
master.dat<-dplyr::select(master.dat,Taxa,INC, COLD, plate_num,date,viable,start,end,tru.daily)

master.dat<-filter(master.dat,plate_num!=0)
master.dat$force<-ifelse(master.dat$INC=="H",25,20)

###Omit rows where start= 25 and germination is 0, iew complete germination
#master.dat<-filter(master.dat, start!=25 & tru.daily!=0)
master.dat<-filter(master.dat, Taxa!="Impatiens capensis")
master.dat<-filter(master.dat, Taxa!="Phlox cuspidata")
master.dat<-filter(master.dat, Taxa!="Carex grisea")
master.dat<-filter(master.dat, Taxa!="Oenethera biennis")
master.datL<-filter(master.dat,INC=="L")
master.datH<-filter(master.dat,INC=="H")
library(drc)
master.dat$indentifyer <- paste(master.dat$Taxa,master.dat$INC,master.dat$COLD)

specieslist <- unique(master.dat$Taxa)
treatlist <- unique(master.dat$COLD)

##this combines cold. but I cant get them to run together
listE50<-list()
for (sp in seq_along(specieslist)){
  dataonesp <- subset(master.dat, Taxa==specieslist[sp])
  mod <- drm(tru.daily~start+end, factor(INC), data = dataonesp, fct = LL.2(), type = "event")
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



df<-as.data.frame(listE50)

df2<-as.data.frame(df2<-t(df))
class(df2)
df2<-df2 %>% rownames_to_column(var="goo")
df2<-df2[grep("Estimate", df2$goo),] 

df2$scratch <- sapply(strsplit(as.character(df2$goo),"."), "[", 1)
df2$Genus <- sapply(strsplit(as.character(df2$goo),'.'), "[", 2)
df2$species <- sapply(strsplit(as.character(df2$goo),'.'), "[", 3)

df2<-df2 %>% separate(goo, c("scratch","Genus","species","treatment"))
unique(df2$Genus)

