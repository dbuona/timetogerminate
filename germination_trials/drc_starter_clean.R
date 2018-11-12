##This is the main cleaning file for germination work. Begin on november 12, 2018 by Dan B

rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")

chill.no<-read.csv("dat_CHILL_NO.CSV",header=TRUE)
colnames(chill.no)
colnames(chill.no)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/27/18","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/21/18","viable","chill_germ","mold","adj_total"))
colnames(chill.no)

##make all nas 0
chill.no<-chill.no %>% replace(is.na(.), 0)

##3adjust germination percentages

####1] total germinated/ totalviable+total germ
###could also think of mold differently
chill.no$adj_total<-chill.no$"9/21/18"+chill.no$viable

##reformat data
chill.no<-gather(chill.no,"date","germ_num",7:22)
chill.no$germ_num<-as.numeric(chill.no$germ_num)

chill.no$date<-as.Date(chill.no$date,format =  "%m/%d/%y")
class(chill.no$Date)
chill.no$day<-yday(chill.no$date)
unique(chill.no$day)
start<-yday("2018/08/27")
chill.no$start<-chill.no$day-start

###make NA's zero


###make germination percent
#chill.no$germ_perc<-chill.no$germ_num/chill.no$adj_total

####format for drc function
daty<- chill.no %>% group_by(plate_num) %>% mutate(germ.daily = germ_num - lag(germ_num))
daty<- daty%>% group_by(plate_num) %>% mutate(tru.daily = lead(germ.daily))

daty<- daty %>% group_by(plate_num) %>% mutate(end = lead(start))
daty$germ.daily[is.na(daty$germ.daily)]<-0 
colnames(daty)
#dater.no<-dplyr::select(daty,tot_seed,Taxa,INC, COLD, plate_num,date,start,end,tru.daily)
dater.no<-daty
###make 
dater.no$end<-ifelse(is.na(dater.no$end),Inf,dater.no$end)
dater.no$tru.daily<-ifelse(is.na(dater.no$tru.daily),dater.no$viable,dater.no$tru.daily)

###select colums
dater.no<-dplyr::select(dater.no,Taxa,INC, COLD, plate_num,date,start,end,tru.daily,viable)

dater.no<-filter(dater.no,plate_num!=0)
dater.no$force<-ifelse(dater.no$INC=="H",25,20)
#####this code will ultimately be moved to other scripts when modeling begins in earnest

goob<-filter( dater.no, Taxa %in% c("Asclepias syriaca","Hesperis matronalis","Silene vulgaris"))
germLL.2 <- drm(tru.daily ~ start + end, factor(Taxa):factor(INC), 
                data =goob, fct = LL.2(), type = "event")

plot(germLL.2, ylim=c(0, 1.5))
coef(germLL.2)

## first installing drc and drcData
devtools::install_github("DoseResponse/drcData")
#2devtools::install_github("DoseResponse/drc")
## then installing the development version of medrc
devtools::install_github("DoseResponse/medrc")

library(medrc)
help("broccoli")
bm <- medrm(tru.daily ~ start, curveid=b + c + e + d ~ factor(force),
            random=c + e + d ~ 1|Taxa,
            data=goob, fc=LL.()) 
drcData(broccoli)
