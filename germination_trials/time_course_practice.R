
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/")
library(tidyverse)
library("growthcurver")
library(lubridate)
library("Hmisc")

data<-read.csv("time_course_practice.csv", header=TRUE)
data2<-read.csv("time_course_practice2.csv", header=TRUE)
data3<-read.csv("time_course_practice3.csv", header=TRUE)
data4<-read.csv("time_course_practice4.csv", header=TRUE)
data5<-read.csv("time_course_practice5.csv", header=TRUE)
data6<-read.csv("time_course_practice6.csv",header=TRUE)

data<-dplyr::select(data,-c(X))
data3<-dplyr::select(data3,-c(X,X.1))
data5<-dplyr::select (data5,1:12)
data6<-dplyr::select(data6,1:9)

colnames(data)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/21/18"))
ncol(data)
data2$"9/10/18"<-0
ncol(data2)
colnames(data2)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num", "9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/21/18","9/23/18","9/26/18","9/28/18","9/30/18","10/03/18","10/05/18","9/10/18"))

data3$"9/24/18"<-0
colnames(data3)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","9/26/18","9/28/18","9/30/18","10/03/18","10/05/18","10/07/18","10/08/18","10/10/18","10/12/18","10/14/18","10/15/18","10/17/18","10/19/18","9/24/18"))

colnames(data4)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/01/18","10/03/18","10/05/18","10/07/18","10/08/18","10/10/18","10/12/18","10/14/18","10/15/18","10/17/18","10/19/18"))

data5$"10/08/18"<-0
colnames(data5)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/10/18","10/12/18","10/14/18","10/15/18","10/17/18","10/19/18","10/08/18"))

colnames(data6)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/15/18","10/17/18","10/19/18"))

data<-gather(data,"date","germination",7:21)
data$germination<-as.numeric(data$germination)

data2<-gather(data2,"date","germination",7:20)
data2$germination<-as.numeric(data2$germination)

data3<-gather(data3,"date","germination",7:20)
data3$germination<-as.numeric(data3$germination)

data4<-gather(data4,"date","germination",7:17)
data4$germination<-as.numeric(data4$germination)

data5<-gather(data5,"date","germination",7:13)
data5$germination<-as.numeric(data5$germination)

data6<-gather(data6,"date","germination",7:9)
data6$germination<-as.numeric(data6$germination)


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


data3$date<-as.Date(data3$date,format =  "%m/%d/%y")
class(data3$date)
data3$day<-yday(data3$date)
unique(data3$day)
start3<-yday("2018/09/24")
data3$DAY<-data3$day-start3

data4$date<-as.Date(data4$date,format =  "%m/%d/%y")
class(data4$date)
data4$day<-yday(data4$date)
unique(data4$day)
start4<-yday("2018/10/01")
data4$DAY<-data4$day-start4

data5$date<-as.Date(data5$date,format =  "%m/%d/%y")
class(data5$date)
data5$day<-yday(data5$date)
unique(data5$day)
start5<-yday("2018/10/08")
data5$DAY<-data5$day-start5

data6$date<-as.Date(data6$date,format =  "%m/%d/%y")
class(data6$date)
data6$day<-yday(data6$date)
unique(data6$day)
start6<-yday("2018/10/15")
data6$DAY<-data6$day-start6




data$germination<-ifelse(is.na(data$germination),0,data$germination)
data2$germination<-ifelse(is.na(data2$germination),0,data2$germination)
data3$germination<-ifelse(is.na(data3$germination),0,data3$germination)
data4$germination<-ifelse(is.na(data4$germination),0,data4$germination)
data5$germination<-ifelse(is.na(data5$germination),0,data5$germination)
data6$germination<-ifelse(is.na(data6$germination),0,data6$germination)

#ggplot(data, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)
#ggplot(data2, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)
#ggplot(data3, aes(x = DAY, y = germination, color=INC)) + stat_summary(alpha=0.7)+facet_wrap(~Taxa)+theme_bw()+geom_line(stat = "summary", fun.y = mean)

full<-rbind(data,data2,data3,data4,data5,data6)

full$INC<-ifelse(full$INC=="H", "High","Low")
table(full$COLD)
unique(full$Taxa)

full2<-dplyr::filter(full,Taxa %in% c("Oenethera biennis","Cryptotaenia canadensis", "Hesperis matronalis","Polygonum virginiatum","Asclepias syriaca","Silene stellata","Silene vulgaris","Eurbia diviricata","Thalictrum dioicum","Anemone virginana"))
ggplot(full, aes(x = DAY, y = germination, color=COLD, shape=INC)) + stat_summary(alpha=0.7)+facet_grid(Taxa~INC)+theme_bw()+geom_line(stat = "summary", fun.y = mean)+scale_color_manual(values=c("orange", "dodgerblue", "purple","darkgreen", "red","deeppink"))
unique(full$Taxa)

goodsp<-dplyr::filter(full, Taxa %in% c("Asclepias syriaca","Eurbia diviricata", "Cryptotaenia canadensis","Hesperis matronalis","Silene stellata","Polygonum virginiatum"))
#goodsp<-dplyr::filter(goodsp, COLD!="E") 
ggplot(goodsp, aes(x = DAY, y = germination, color=COLD, shape=INC))+facet_wrap(~Taxa)+stat_summary(alpha=0.7)+theme_linedraw()+geom_line(stat = "summary", fun.y = mean, aes(linetype=INC))+ylab("Germination (number of seeds)")+scale_linetype_manual(name="Incubation level",values=c("solid","dashed"))+scale_shape_manual(name="Incubation level",values=c(19,17))+scale_color_manual(values=c("black","sienna4","orangered1","orchid1","purple3","blue"),labels=c("0 Days", "14 days", "28 days", "35 days","42 days","49 days"),name="Stratification Level")

goodsp$germination_percent<-goodsp$germination/20

####plot for community rank
ggplot(goodsp, aes(x = DAY, y = germination, color=Taxa )) + stat_summary(alpha=0.7)+facet_grid(COLD~INC)+theme_bw()+geom_line(stat = "summary", fun.y = mean)
##This is where it gets wonky
AS<-filter(full, Taxa=="Asclepias syriaca")
ggplot(AS, aes(x = DAY, y = germination, color=COLD, shape=INC)) + stat_summary(alpha=0.7)+theme_bw()+geom_line(stat = "summary", fun.y = mean)+scale_color_manual(values=c("orange", "dodgerblue", "purple","darkgreen"), labels=c("0 Days", "14 days", "28 days", "35 days"),name="Stratification Level")

###

full2<-data%>% group_by(Taxa,plate_num,INC,COLD,DAY,tot_seed) %>% summarise(germ=mean(germination))
full2$perc<-full2$germ/full2$tot_seed
full2$uni<-paste(full2$Taxa,full2$INC,spe=".")



for(i in (unique(full2$plate_num))){
  XX<-data.frame(full2$plate_num,start = c(0,full2$DAY), end = c(full2$DAY,Inf))
  }
  
  nam <- paste( "data", sep = ".")
  (assign(nam, df[df$uni==i,]))
 

head$framz


list <- list(paste("data", unique(df$uni) , sep = "."))


lapply(list, function(x) data.frame(start = c(0,x$DAY), end = c(x$DAY,Inf)))






